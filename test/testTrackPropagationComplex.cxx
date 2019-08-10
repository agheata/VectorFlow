/// +--------------------------------------------------------------------------------------+
/// | VectorFlow unit test. Generation + Propagation.                                      |
/// |                                                                                      |
/// | This test is intended to show the usage of the ComplexFlow interface of VectorFlow.  |
/// | The example propagates charged particles in a constant field in several steps,       |
/// | in a geometry of adjacent concentric cylindrical layers                              |
/// |                                                                                      |
/// | The pipeline flow has only one task (TaskLayerPropagator) with multiple instances    |
/// | (one per layer). Each task has up to 2 client tasks (the inner and/or outer layer    |
/// | task). The data type to be used are the tracks that are generated in each even of    |
/// | the event-loop.                                                                      |
/// +--------------------------------------------------------------------------------------+

#include <iostream>
#include <string>
#include <vector>
#include <vectorFlow/ComplexFlow.h>
#include "CocktailGenerator.h"
#include "Event.h"
#include "Track.h"
#include "Particle.h"
#include "vectorFlow/SystemOfUnits.h"
#include "base/Vector3D.h"
#include "HelixPropagator.h"
#include "SimpleStepper.h"
#include "timer.h"

using namespace vectorflow;

// For benchmarking purposes
#ifdef VECCORE_TIMER_CYCLES
  using time_unit = cycles;
  static const char *time_unit_name = "cycles";
#else
  using time_unit = milliseconds;
  static const char *time_unit_name = "ms";
#endif

//=============================================================================
struct TaskLayerPropagator : public Work<Track, std::vector<Track*>> {
  // Propagator task needs a stepper
  int fLayer = -1;
  trackml::HelixPropagator *fPropagator = nullptr;
  trackml::SimpleStepper   *fStepper    = nullptr;
  bool fVerbose = true;

  // Struct constructor
  TaskLayerPropagator(int layer, vecgeom::Vector3D<double> const &bfield) {
    fLayer = layer;
    fPropagator = new trackml::HelixPropagator(bfield);
    fStepper    = new trackml::SimpleStepper(fPropagator);
  }

  ~TaskLayerPropagator() {
    delete fPropagator;
    delete fStepper;
  }

  // Scalar mode executor
  void Execute(Track* track) {
    // Propagate track to the boundary of the current cylindrical layer
    fStepper->PropagateInTube(fLayer, *track);
    if (fVerbose) {
      if (track->Status() == vectorflow::kKilled) {
        std::cout << "Track " << track->PrimaryParticleIndex() << " from event "
                  << track->Event() << " made " << track->GetNsteps() << " steps. ";
        std::cout << "Exit position: " << track->Position() << "\n";
        return;
      }
    }

    // Dispatch to next layer. Client 0 corresponds to the inner neighbour
    // layer, client 1 to the outer one
    bool move_outer = track->Position().Dot(track->Direction()) > 0;
    size_t client = 0;
    if (move_outer && fLayer > 0 && fLayer < 3) client = 1;

    Dispatch(track, client);
  }

  // Vector mode executor
  void Execute(std::vector<Track*> const &tracks) {
    std::cout << "TaskPropagator::Execute not implemented for vector mode.\n";
  }
};

//=============================================================================
void Propagate(std::vector<Track> &tracks, 
               ComplexFlow<Track, std::vector<Track*>, 4> &flow) {
  for (auto &track : tracks) {
    // Only charged tracks will be added
    if (track.Charge() != 0)
      flow.AddData(&track);
  }
  // Process the flow
  flow.Execute();
}

//=============================================================================
unsigned long long RunTest(ComplexFlow<Track, std::vector<Track *>, 4> &flow,
                           CocktailGenerator::Event_t const *event,
                           std::vector<Track> &tracks) {
  auto nPrimaries = event->GetNprimaries();
  for (auto i = 0; i < nPrimaries; i++)
    tracks[i].Reset(*event->GetPrimary(i));

  Timer<time_unit> timer; // for benchmarking purposes
  timer.Start();

  Propagate(tracks, flow);

  unsigned long long t = timer.Elapsed();

  // Clearing all stages of the flow
  for (auto i = 0; i < 4; ++i) flow.Clear(i);
  return t;
}

//=============================================================================
int main(int argc, char* argv[]) {
  using namespace geant::units;
  // Read number of events to be generated
  int nTracks = 50;
  int nTries  = 1;

  if (argc > 1)
    nTracks = atoi(argv[1]);
  if (argc > 2)
    nTries  = atoi(argv[2]);

  // Add CocktailGenerator
  vecgeom::Vector3D<double> vertex(0., 0., 10.);
  CocktailGenerator cocktailGen(vertex);

  // Add particle species for the generator
  cocktailGen.AddPrimary("pi+", 0.3);
  cocktailGen.AddPrimary("pi-", 0.3);
  cocktailGen.AddPrimary("gamma", 0.4);

  // Set parameters of the generator
  cocktailGen.SetPrimaryEnergyRange(0.1 * GeV, 10 * GeV);
  cocktailGen.SetMaxPrimaryPerEvt(nTracks);
  cocktailGen.SetAvgPrimaryPerEvt(nTracks);
  cocktailGen.SetVertex(vertex);
  cocktailGen.SetMaxDepth(2);

  // Check if generator parameters are correct
  if (!cocktailGen.InitPrimaryGenerator()) {
    std::cout << "Failed to pass init check.\n";
    exit(-1);
  }

  // Initialize an helix propagator in field and a simple stepper
  vecgeom::Vector3D<double> bfield(0., 0., 20. * kilogauss);

  // Create a complex flow where each geometry layer is a stage
  using LayerFlow_t = ComplexFlow<Track, std::vector<Track *>, 4>;

  LayerFlow_t flow;

  // Create and initialize the layer propagator tasks
  TaskLayerPropagator *propagators[4];
  for (auto i = 0; i < 4; ++i) {
    propagators[i] = new TaskLayerPropagator(i, bfield);
    flow.AddWork(propagators[i], i);
    flow.SetVectorMode(i, false); // scalar flow
  }

  // Connect clients to tasks
  for (auto i = 0; i < 4; ++i) {
    if (i > 0) propagators[i]->AddClient(&flow.InputData(i - 1));
    if (i < 3) propagators[i]->AddClient(&flow.InputData(i + 1));
  }

  // Run one event for benchmarking purposes
  CocktailGenerator::Event_t *event = cocktailGen.NextEvent();

  // Make a copy of the tracks, so that we can run the test multiple times
  std::vector<Track> tracks(event->GetNprimaries());

  std::cout << "\n--EXECUTING IN SCALAR MODE--\n";
  double sum = 0.;
  for (auto i = 0; i < nTries; i++) {
    auto t = RunTest(flow, event, tracks);
    sum += (double)t;
  }
  double average = sum / nTries;
  std::cout << "\nExecution time:   " << average << " [" << time_unit_name << "]\n";

  // Clearing created pointers
  event->Clear();
  for (auto i = 0; i < 4; ++i) delete propagators[i];

  // End of test
  std::cout << "\n--END OF TEST--\n";
  return 0;
}
