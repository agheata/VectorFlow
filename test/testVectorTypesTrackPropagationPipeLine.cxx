/// +--------------------------------------------------------------------------------------+
/// | VectorFlow unit test. Generation + Propagation. --RUNNING IN VECTOR MODE--
/// | | | | This test is intended to show the usage of the PipelineFlow
/// interface of VectorFlow. | | The example propagates charged particles in a
/// constant field in several steps,       | | until they reach a sphere with a
/// given radius centered in the event vertex.          | | | | The pipeline
/// flow has only one task (TaskPropagator) and the data type to be used    | |
/// are the tracks that are generated in each even of the event-loop. | | After
/// each execution, the pipeline is cleared and the next tracks are processed. |
/// +--------------------------------------------------------------------------------------+

#include "CocktailGenerator.h"
#include "Event.h"
#include "Geant/ConstFieldHelixStepper.h"
#include "HelixPropagator.h"
#include "Particle.h"
#include "SimpleStepper.h"
#include "Track.h"
#include "base/Vector3D.h"
#include "timer.h"
#include "vectorFlow/SystemOfUnits.h"
#include "vectorFlow/VectorTypes.h"
#include <VecCore/VecCore>
#include <benchmark/benchmark.h>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <vectorFlow/PipelineFlow.h>
#ifdef USE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

using namespace vectorflow;
using vecCore::Get;
using vecCore::Load;
using vecCore::Scalar;
using vecCore::Set;
using vecCore::Store;

// for benchmarking purposes
#ifdef VECCORE_TIMER_CYCLES
using time_unit = cycles;
static const char *time_unit_name = "cycles";
#else
using time_unit = milliseconds;
static const char *time_unit_name = "ms";
#endif

template <typename Data> struct Tracks_v {
  // Data in vecCore vector types
  vecgeom::Vector3D<Double_v> fPos_v;
  vecgeom::Vector3D<Double_v> fDir_v;
  Double_v fCharge_v, fMomentum_v, fNSteps_v;

  // Auxiliar vectorized methods
  Double_v Pt_v() const { return fMomentum_v * fDir_v.Perp(); }

  Double_v Curvature_v(double Bz) {
    using constants::kB2C;
    using constants::kTiny;
    Double_v qB_v = fCharge_v * Bz;
    Double_v curvature_v = vecCore::math::Abs(kB2C * qB_v / (Pt_v() + kTiny));
    return vecCore::Blend<Double_v>(vecCore::math::Abs(qB_v) < kTiny, kTiny,
                                    curvature_v);
  }

  // Check if all tracks are normalized
  VECGEOM_FORCE_INLINE
  bool IsNormalized() const {
    Double_v norm = fDir_v.Mag2();
    MaskD_v is_norm =
        1. - vecgeom::kTolerance < norm && norm < 1 + vecgeom::kTolerance;
    return vecCore::MaskFull(is_norm);
  }

  // Gather info to fill one SIMD lane
  void Gather(Data *track, const std::size_t lane) {
    Set(fPos_v.x(), lane, track->Position().x());
    Set(fPos_v.y(), lane, track->Position().y());
    Set(fPos_v.z(), lane, track->Position().z());
    Set(fDir_v.x(), lane, track->Direction().x());
    Set(fDir_v.y(), lane, track->Direction().y());
    Set(fDir_v.z(), lane, track->Direction().z());
    Set(fCharge_v, lane, track->Charge());
    Set(fMomentum_v, lane, track->P());
    Set(fNSteps_v, lane, track->GetNsteps());
  }

  // Scatter one lane info back to original structure
  void Scatter(const std::size_t lane, Data *track) {
    track->SetPosition(Get(fPos_v.x(), lane), Get(fPos_v.y(), lane),
                       Get(fPos_v.z(), lane));
    track->SetDirection(Get(fDir_v.x(), lane), Get(fDir_v.y(), lane),
                        Get(fDir_v.z(), lane));
    track->SetNsteps(Get(fNSteps_v, lane));
  }
};

//=============================================================================
struct TaskPropagator : public Work<Track, std::vector<Track *>> {
  // Propagator task needs a stepper
  trackml::HelixPropagator *fPropagator = nullptr;
  trackml::SimpleStepper *fStepper = nullptr;
  ConstFieldHelixStepper *fHelixStepper = nullptr;
  Tracks_v<Track> fTracks_v;
  vecgeom::Vector3D<double> fPosChecksum = 0;
  size_t fNstepsSum = 0;
  bool fVerbose = false;

  // Task struct constructor
  TaskPropagator(vecgeom::Vector3D<double> const &bfield) {
    fPropagator = new trackml::HelixPropagator(bfield);
    fStepper = new trackml::SimpleStepper(fPropagator);
    fHelixStepper = new ConstFieldHelixStepper(bfield);
  }

  ~TaskPropagator() {
    delete fPropagator;
    delete fStepper;
    delete fHelixStepper;
  }

  // Scalar mode executor
  void Execute(Track *track) {
    // Propagate track to the boundary of a sphere of radius 20cm
    const double radius = 20. * geant::units::cm;
    /*fStepper->*/ PropagateToR(radius, *track);
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
    PropagateToR(radius, tracks);
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

  // Function to check whether all lanes can still be propagated or not
  bool AllLanesCanDoStep(const Double_v &c_v) {
    return vecCore::MaskFull(c_v < 0);
  }

  // PropagateToR function for a track
  void PropagateToR(double radius, Track &track) const {
    // Compute safe distance for the track to getting outside the sphere of
    // given radius
    using namespace vecgeom;
    constexpr double epsilon = 1.E-4 * geant::units::mm;
    constexpr double toKiloGauss =
        1.0 / geant::units::kilogauss; // Converts to kilogauss

    Vector3D<double> pos = track.Position();
    Vector3D<double> dir = track.Direction();
    Vector3D<double> const bfield = fPropagator->GetBfield();
    const double radius2 = radius * radius;
    double rad2 = pos.Mag2();
    double c = rad2 - radius2;

    // Enter the stepping loop while the track is inside the sphere of given
    // radius
    while (c < 0) {
      double safety = radius - Sqrt(rad2);
      // Maximum allowed distance so that the sagitta of the track trajectory
      // along this distance is less than a fraction of epsilon of the distance
      double dmax = 8. * epsilon / track.Curvature(bfield.z() * toKiloGauss);

      // Compute distance along straight line to exit the sphere
      double pDotV = pos.Dot(dir);
      double d2 = pDotV * pDotV - c;
      double snext = -pDotV + Sqrt(vecCore::math::Abs(d2));

      track.SetSafety(safety);
      track.SetSnext(snext);

      double step_geom = Max(epsilon, track.GetSnext());
      double step_field = Max(dmax, safety);
      double step = Min(step_geom, step_field);

      // Propagate in field with step distance (update position and direction)
      fPropagator->Propagate(track, step);

      // Update track state other than position/direction
      track.DecreasePstep(step);
      track.DecreaseSnext(step);
      track.DecreaseSafety(step);
      track.IncreaseStep(step);
      track.IncrementNsteps();

      // Update local variables
      pos = track.Position();
      dir = track.Direction();
      rad2 = pos.Mag2();
      c = rad2 - radius2;
    }
  }

  // PropagateToR function for a vector of tracks
  void PropagateToR(const double radius, std::vector<Track *> const &tracks) {
    const std::size_t kTracksSize = tracks.size();
    const std::size_t kVectorSize = vecCore::VectorSize<Double_v>();

    // Execute scalar mode if there are less tracks than vector lanes
    if (kTracksSize < kVectorSize) {
      for (auto track : tracks)
        PropagateToR(radius, *track);
      return;
    }

    std::size_t lane[kVectorSize]; // Array to check which track is which lane
    std::size_t nextTrack = 0;     // Counter to dispatch tracks

    // Declare constants and auxiliar variables
    const double epsilon = 1.E-4 * geant::units::mm;
    const double toKiloGauss = 1.0 / geant::units::kilogauss;

    vecgeom::Vector3D<double> const bfield = fPropagator->GetBfield();
    const double radius2 = radius * radius;
    Double_v rad2_v, c_v, safety_v, dmax_v, pDotV_v, d2_v, snext_v, step_geom_v,
        step_field_v, step_v;

    // Set up first lanes
    for (auto i = 0; i < kVectorSize; i++) {
      lane[i] = i;
      fTracks_v.Gather(tracks[i], i);
    }
    nextTrack += kVectorSize;
    rad2_v = fTracks_v.fPos_v.Mag2();
    c_v = rad2_v - radius2;
    bool ongoing = vecCore::MaskFull(c_v < 0);

    // Execute until all lanes can at least do one more step
    while (ongoing) {
      safety_v = radius - vecCore::math::Sqrt(rad2_v);
      dmax_v = 8. * epsilon / fTracks_v.Curvature_v(bfield.z() * toKiloGauss);
      pDotV_v = fTracks_v.fPos_v.Dot(fTracks_v.fDir_v);
      d2_v = pDotV_v * pDotV_v - c_v;
      snext_v = -pDotV_v + vecCore::math::Sqrt(vecCore::math::Abs(d2_v));

      step_geom_v = vecCore::math::Max(static_cast<Double_v>(epsilon), snext_v);
      step_field_v = vecCore::math::Max(dmax_v, safety_v);
      step_v = vecCore::math::Min(step_geom_v, step_field_v);

      fHelixStepper->DoStep<Double_v>(
          fTracks_v.fPos_v, fTracks_v.fDir_v, fTracks_v.fCharge_v,
          fTracks_v.fMomentum_v, step_v, fTracks_v.fPos_v, fTracks_v.fDir_v);

      assert(fTracks_v.IsNormalized() &&
             "ERROR: Direction not normalized after field propagation");

      // Calculate c_v values again for the propagated tracks
      rad2_v = fTracks_v.fPos_v.Mag2();
      c_v = rad2_v - radius2;

      fTracks_v.fNSteps_v += 1;

      // Assign new track to a lane if previous was fully propagated
      if (!vecCore::MaskFull(c_v < 0)) {
        // There are finished lanes, check if we have enough tracks to refill
        bool refill = (kTracksSize - nextTrack) > kVectorSize;
        if (!refill)
          refill = (kTracksSize - nextTrack) > MaskCount(c_v >= 0);
        if (refill) {
          // We have enough tracks to refill the done lanes
          for (auto i = 0; i < kVectorSize; i++) {
            if (Get(c_v, i) >= 0) {
              // Scatter back the propagated track
              fTracks_v.Scatter(i, tracks[lane[i]]);

              // Set values for new track in its corresponding lane
              fTracks_v.Gather(tracks[nextTrack], i);
              lane[i] = nextTrack++;
            }
          }
          // Update rad2_v and c_v (even if possibly just one lane changed)
          rad2_v = fTracks_v.fPos_v.Mag2();
          c_v = rad2_v - radius2;
        } else {
          // Not enough tracks to refill the vector mode, scatter existing ones
          for (auto i = 0; i < kVectorSize; i++) {
            // Scatter back the propagated track
            fTracks_v.Scatter(i, tracks[lane[i]]);
            if (Get(c_v, i) < 0) {
              // Track not fully propagated yet, propagate in scalar mode
              PropagateToR(radius, *tracks[lane[i]]);
            }
          }
          ongoing = false;
        }
      }
    }

    // Execute remaining tracks in scalar mode
    while (nextTrack < kTracksSize)
      PropagateToR(radius, *tracks[nextTrack++]);
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

  Timer<time_unit> timer; // timer for benchmarking
  timer.Start();

  Propagate(tracks, flow);

  unsigned long long t = timer.Elapsed();

  flow.Clear();
  return t;
}

//=============================================================================
int main(int argc, char *argv[]) {
  // Read number of events to be generated
  int nTracks = 100000;
  int nTries = 10;
  bool singletry = false;
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
  cocktailGen.SetPrimaryEnergyRange(0.1 * geant::units::GeV,
                                    10 * geant::units::GeV);
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

  // Event loop
  //  int nEvents = 1; // benchmarking purpose
  //  for (auto i = 0; i < nEvents; i++) {
  CocktailGenerator::Event_t *event = cocktailGen.NextEvent();
  //  event->SetEvent(i);

  // Make a copy of the tracks, so that we can run multiple times
  std::vector<Track> tracks(event->GetNprimaries());

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
  double sum = 0., sumv = 0.;
  double sumsq = 0., sumsqv = 0.;

  // Set the task to be executed in scalar mode
  plFlow.SetVectorMode(kPropagatorStage, false);
  // std::cout << "--warming up ...\n";
  for (auto i = 0; i < nwarm_up; ++i)
    RunTest(plFlow, event, tracks);
  tPropagate.fPosChecksum.Set(0., 0., 0.);
  tPropagate.fNstepsSum = 0;

  std::cout << "--EXECUTING IN SCALAR MODE--\n";
  for (auto i = 0; i < nTries; ++i) {
    auto t = RunTest(plFlow, event, tracks);
    // std::cout << t << std::endl;
    sum += (double)t;
    sumsq += (double)t * t;
  }
  vecgeom::Vector3D<double> posChecksum = tPropagate.fPosChecksum;
  size_t nstepssum = tPropagate.fNstepsSum;

  double average = sum / nTries;
  double sigma = std::sqrt(nTries * sumsq - sum * sum) / nTries;
  double percent = sigma / average;
  std::cout << "Execution time:   " << average << " +/- " << sigma << " ["
            << time_unit_name << "] (" << 100 * percent << " %)\n";

  tPropagate.fPosChecksum.Set(0., 0., 0.);
  tPropagate.fNstepsSum = 0;

  // Set the task to be executed in vector mode
  plFlow.SetVectorMode(kPropagatorStage, true);
  // std::cout << "--warming up ...\n";
  for (auto i = 0; i < nwarm_up; ++i)
    RunTest(plFlow, event, tracks);
  tPropagate.fPosChecksum.Set(0., 0., 0.);
  tPropagate.fNstepsSum = 0;

  std::cout << "--EXECUTING IN VECTOR MODE--\n";
  for (auto i = 0; i < nTries; ++i) {
    auto t = RunTest(plFlow, event, tracks);
    // std::cout << t << std::endl;
    sumv += (double)t;
    sumsqv += (double)t * t;
  }
  if ((posChecksum - tPropagate.fPosChecksum).Mag() > 1e-10 ||
      (nstepssum != tPropagate.fNstepsSum))
    std::runtime_error("Checksums not match between scalar and vector modes");

  double averagev = sumv / nTries;
  double sigmav = std::sqrt(nTries * sumsqv - sumv * sumv) / nTries;
  double percentv = sigmav / averagev;
  std::cout << "Execution time:   " << averagev << " +/- " << sigmav << " ["
            << time_unit_name << "] (" << 100 * percentv << " %)\n";

  double speedup = average / averagev;
  double eps = std::sqrt(percent * percent + percentv * percentv);
  double sigma_speedup = eps * speedup;

  std::cout << "=== Speedup:   " << speedup << " +/- " << sigma_speedup << " ("
            << 100 * eps << " %)\n";

  // Clear event
  event->Clear();
  //  } // end of event-loop

  // End of test
  std::cout << "\n--END OF TEST--\n";
  return 0;
}
