/// +--------------------------------------------------------------------------------------+
/// | VectorFlow unit test. Generation + Propagation. --RUNNING IN VECTOR MODE--           |
/// | This test is intended to show the usage of the PipelineFlow                          |
/// | interface of VectorFlow. The example propagates charged particles in a               |
/// | constant field in several steps, until they reach a sphere with a                    |
/// | given radius centered in the event vertex. The pipeline                              |
/// | flow has only one task (TaskPropagator) and the data type to be used                 |
/// | are the tracks that are generated in each even of the event-loop. After              |
/// | each execution, the pipeline is cleared and the next tracks are processed.           |
/// +--------------------------------------------------------------------------------------+

#include <iostream>
#include <string>
#include <vector>
#include <vectorFlow/vectorFlow>
#include "CocktailGenerator.h"
#include "Event.h"
#include "Track.h"
#include "Particle.h"
#include "base/Vector3D.h"
#include "HelixPropagator.h"
#include "SimpleStepper.h"
#include "timer.h"

#ifdef USE_GPERFTOOLS
  #include <gperftools/profiler.h>
#endif

// for benchmarking purposes
#ifdef VECCORE_TIMER_CYCLES
  using time_unit = cycles;
  static const char *time_unit_name = "cycles";
#else
  using time_unit = milliseconds;
  static const char *time_unit_name = "ms";
#endif

using namespace vectorflow;

//=============================================================================
struct TaskPropagator : public Work<Track, std::vector<Track *>> {
  // Propagator task needs a stepper
  trackml::HelixPropagator *fPropagator  = nullptr;
  trackml::SimpleStepper   *fStepper     = nullptr;
  vecgeom::Vector3D<double> fPosChecksum = 0;
  size_t fNstepsSum = 0;
  bool fVerbose = false;

  // Task struct constructor
  TaskPropagator(vecgeom::Vector3D<double> const &bfield) {
    fPropagator = new trackml::HelixPropagator(bfield);
    fStepper = new trackml::SimpleStepper(fPropagator);
  }

  ~TaskPropagator() {
    delete fPropagator;
    delete fStepper;
  }

  // Scalar mode executor
  void Execute(Track *track) {
    // Propagate track to the boundary of a sphere of radius 20cm
    const double radius = 20. * geant::units::cm;
    fStepper->PropagateToR(radius, *track);
    fPosChecksum += track->Position();
    fNstepsSum += track->GetNsteps();
    if (fVerbose) {
      std::cout << "Track " << track->Index() << " made " << track->GetNsteps()
                << " steps. ";
      std::cout << "Exit position: " << track->Position() << "\n";
    }
  }

  // Vector mode executor
  void Execute(std::vector<Track *> const &tracks) {
    // Propagate track to the boundary of a sphere of radius 20cm
    const double radius = 20. * geant::units::cm;
    fStepper->PropagateToR(radius, tracks);
    for (const auto &track : tracks) {
      fPosChecksum += track->Position();
      fNstepsSum += track->GetNsteps();
      if (fVerbose) {
        std::cout << "Track " << track->Index() << " made "
                  << track->GetNsteps() << " steps. ";
        std::cout << "Exit position: " << track->Position() << "\n";
      }
    }
  }
};

//=============================================================================
void Propagate(std::vector<Track> &tracks,
               PipelineFlow<Track, std::vector<Track *>, 1> &flow) {
  for (auto &track : tracks) {
    // Only charged tracks will be added
    if (track.Charge() != 0)
      flow.AddData(&track);
  }
  
  // Process the flow
  flow.Execute();
}

//=============================================================================
unsigned long long RunTest(PipelineFlow<Track, std::vector<Track *>, 1> &flow,
                           CocktailGenerator::Event_t const *event,
                           std::vector<Track> &tracks) {
  auto nPrimaries = event->GetNprimaries();
  for (auto in = 0; in < nPrimaries; in++)
    tracks[in].Reset(*event->GetPrimary(in));

  Timer<time_unit> timer; // For benchmarking purposes
  timer.Start();

  Propagate(tracks, flow);

  unsigned long long t = timer.Elapsed();

  flow.Clear();
  return t;
}

//=============================================================================
int main(int argc, char *argv[]) {
  // Read number of events to be generated
  int nTracks      = 100000;
  int nTries       = 10;
  bool singletry   = false;
  bool vector_mode = false;

  if (argc == 1) {
    std::cout << "===  Usage:   testVectorTypesTrackPropagationPipeline "
                 "[ntracks] [nrepetitions] [s/v]\n";
  };

  if (argc > 1)
    nTracks = atoi(argv[1]);
  if (argc > 2)
    nTries = atoi(argv[2]);
  if (argc > 3) {
    singletry = true;
    if (!strcmp(argv[3], "v"))
      vector_mode = true;
  }

  // Add CocktailGenerator
  vecgeom::Vector3D<double> vertex(0., 0., 10.);
  CocktailGenerator cocktailGen(vertex);

  // Add particle species for the generator
  cocktailGen.AddPrimary("pi+", 0.3);
  cocktailGen.AddPrimary("pi-", 0.3);
  cocktailGen.AddPrimary("gamma", 0.4);

  // Set parameters of the generator
  cocktailGen.SetPrimaryEnergyRange(0.1 * geant::units::GeV, 10 * geant::units::GeV);
  cocktailGen.SetMaxPrimaryPerEvt(nTracks);
  cocktailGen.SetAvgPrimaryPerEvt(nTracks);
  cocktailGen.SetVertex(vertex);
  cocktailGen.SetMaxDepth(2);

  // Check if generator parameters are correct
  if (!cocktailGen.InitPrimaryGenerator()) {
    std::runtime_error("Failed to pass init check.");
  }

  // Constant magnetic field
  vecgeom::Vector3D<double> bfield(0., 0., 20. * geant::units::kilogauss);

  // Create and initialize the propagator task
  TaskPropagator tPropagate(bfield);

  // Create a pipeline flow with one stage
  using Pipeline_t = PipelineFlow<Track, std::vector<Track *>, 1>;
  Pipeline_t plFlow;

  // Add the task to the flow
  static constexpr size_t kPropagatorStage = 0;
  plFlow.AddWork(&tPropagate, kPropagatorStage);

  // Run one event for benchmarking purposes
  CocktailGenerator::Event_t *event = cocktailGen.NextEvent();

  // Make a copy of the tracks, so that we can run the test multiple times
  std::vector<Track> tracks(event->GetNprimaries());

  // SINGLETRY FOR PROFILING
  // =============================================================================
  if (singletry) {
    if (vector_mode) std::cout << "--EXECUTING SINGLE TRY IN VECTOR MODE--\n";
    else             std::cout << "--EXECUTING SINGLE TRY IN SCALAR MODE--\n";
    
    plFlow.SetVectorMode(kPropagatorStage, vector_mode);
    
    #ifdef USE_GPERFTOOLS
      // CPU profiling using gperftools
      std::cout << "=== Profiling using gperftools. Output goes to gperfprof.out\n";
      ProfilerStart("gperfprof.out");
    #endif
    
    auto time = RunTest(plFlow, event, tracks);
    
    #ifdef USE_GPERFTOOLS
      ProfilerStop();
      std::cout << "=== Profiling stopped\n";
    #endif
    
    std::cout << "Execution time:   " << time << " [" << time_unit_name << "]\n";
    return 0;
  }

  int nwarm_up = 20;
  double sum   = 0., sumv   = 0.;
  double sumsq = 0., sumsqv = 0.;
  
  // SCALAR FLOW
  // =============================================================================
  plFlow.SetVectorMode(kPropagatorStage, false);
  // Warm up
  for (auto i = 0; i < nwarm_up; ++i) {
    RunTest(plFlow, event, tracks);
  } // End of warm up
  tPropagate.fPosChecksum.Set(0., 0., 0.);
  tPropagate.fNstepsSum = 0;
  std::cout << "--EXECUTING IN SCALAR MODE--\n";
  for (auto i = 0; i < nTries; ++i) {
    auto t = RunTest(plFlow, event, tracks);
    sum   += (double)t;
    sumsq += (double)t * t;
  }
  vecgeom::Vector3D<double> posChecksum = tPropagate.fPosChecksum;
  size_t nstepssum = tPropagate.fNstepsSum;

  double average = sum / nTries;
  double sigma   = std::sqrt(nTries * sumsq - sum * sum) / nTries;
  double percent = sigma / average;
  std::cout << "Execution time:   " << average << " +/- " << sigma << " ["
            << time_unit_name << "] (" << 100 * percent << " %)\n";

  tPropagate.fPosChecksum.Set(0., 0., 0.);
  tPropagate.fNstepsSum = 0;

  // VECTOR FLOW
  // =============================================================================
  plFlow.SetVectorMode(kPropagatorStage, true);
  // Warm up
  for (auto i = 0; i < nwarm_up; ++i) {
    RunTest(plFlow, event, tracks);
  } // End of warm up
  tPropagate.fPosChecksum.Set(0., 0., 0.);
  tPropagate.fNstepsSum = 0;
  std::cout << "--EXECUTING IN VECTOR MODE--\n";
  for (auto i = 0; i < nTries; ++i) {
    auto t  = RunTest(plFlow, event, tracks);
    sumv   += (double)t;
    sumsqv += (double)t * t;
  }
  if ((posChecksum - tPropagate.fPosChecksum).Mag() > 1e-10 ||
      (nstepssum != tPropagate.fNstepsSum))
    std::runtime_error("Checksums not match between scalar and vector modes");

  double averagev = sumv / nTries;
  double sigmav   = std::sqrt(nTries * sumsqv - sumv * sumv) / nTries;
  double percentv = sigmav / averagev;
  std::cout << "Execution time:   " << averagev << " +/- " << sigmav << " ["
            << time_unit_name << "] (" << 100 * percentv << " %)\n";

  double speedup = average / averagev;
  double eps     = std::sqrt(percent * percent + percentv * percentv);
  double sigma_speedup = eps * speedup;

  std::cout << "=== Speedup:   " << speedup << " +/- " << sigma_speedup << " ("
            << 100 * eps << " %)\n";

  // Clear event
  event->Clear();

  // End of test
  std::cout << "\n--END OF TEST--\n";
  return 0;
}
