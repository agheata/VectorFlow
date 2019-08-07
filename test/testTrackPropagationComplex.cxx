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

using namespace vectorflow;

//struct TaskPropagator : public Work<Event_t, std::vector<Event_t*>> {
struct TaskLayerPropagator : public Work<Track, std::vector<Track*>> {
  // Propagator task needs a stepper
  int fLayer = -1;
  trackml::SimpleStepper* fStepper = nullptr;

  // Scalar mode executor
  void Execute(Track* track) {
    // Propagate track to the boundary of the current cylindrical layer
    fStepper->PropagateInTube(fLayer, *track);
    if (track->Status() == vectorflow::kKilled) {
      std::cout << "Track " << track->PrimaryParticleIndex() << " from event " << track->Event()
                << " made " << track->GetNsteps() << " steps. ";
      std::cout << "Exit position: " << track->Position() << "\n";
      return;
    }
    // Dispatch to next layer. Client 0 corresponds to the inner neighbour layer, client 1 to the outer one
    bool move_outer = track->Position().Dot(track->Direction()) > 0;
    size_t client = 0;
    if (move_outer && fLayer > 0 && fLayer < 3) client = 1;

    Dispatch(track, client);
  }

  // Vector mode executor
  void Execute(std::vector<Track*> const &tracks) {
    std::cout << "TaskPropagator::Execute not implemented for vector mode.\n";
  }

  // Struct constructor
  TaskLayerPropagator(int layer, trackml::SimpleStepper* stepper) : fLayer(layer), fStepper(stepper) {}
};

int main(int argc, char* argv[]) {
  using namespace geant::units;
  // Read number of events to be generated
  int nEvents = argv[1] == NULL ? 10 : atoi(argv[1]);

  // Add CocktailGenerator
  vecgeom::Vector3D<double> vertex(0., 0., 10.);
  CocktailGenerator cocktailGen(vertex);

  // Add particle species for the generator
  cocktailGen.AddPrimary("pi+", 0.3);
  cocktailGen.AddPrimary("pi-", 0.3);
  cocktailGen.AddPrimary("gamma", 0.4);

  // Set parameters of the generator
  cocktailGen.SetPrimaryEnergyRange(0.1 * GeV, 10 * GeV);
  cocktailGen.SetMaxPrimaryPerEvt(100);
  cocktailGen.SetAvgPrimaryPerEvt(70);
  cocktailGen.SetVertex(vertex);
  cocktailGen.SetMaxDepth(2);

  // Check if generator parameters are correct
  if (!cocktailGen.InitPrimaryGenerator()) {
    std::cout << "Failed to pass init check.\n";
    exit(-1);
  }

  // Initialize an helix propagator in field and a simple stepper
  vecgeom::Vector3D<double> bfield(0., 0., 20. * kilogauss);
  trackml::HelixPropagator* propagator = new trackml::HelixPropagator(bfield);
  trackml::SimpleStepper*   stepper    = new trackml::SimpleStepper(propagator);

  // Create a complex flow where each geometry layer is a stage
  using LayerFlow_t = ComplexFlow<Track, std::vector<Track *>, 4>;

  LayerFlow_t flow;

  // Create and initialize the layer propagator tasks
  TaskLayerPropagator *propagators[4];
  for (auto i = 0; i < 4; ++i) {
    propagators[i] = new TaskLayerPropagator(i, stepper);
    flow.AddWork(propagators[i], i);
    flow.SetVectorMode(i, false); // scalar flow
  }

  // Connect clients to tasks
  for (auto i = 0; i < 4; ++i) {
    if (i > 0) propagators[i]->AddClient(&flow.InputData(i - 1));
    if (i < 3) propagators[i]->AddClient(&flow.InputData(i + 1));
  }

  // Event loop
  for (auto i = 0; i < nEvents; i++) {
    CocktailGenerator::Event_t* event = cocktailGen.NextEvent();
    event->SetEvent(i);

    // Add data to the flow
    std::cout << "\n=== Propagating event " << i << "\n";
    for (auto j = 0; j < event->GetNprimaries(); j++) {
      Track* track = event->GetPrimary(j);
      // Only charged tracks will be added
      if (track->Charge() != 0) flow.AddData(track);
    }

    /// Process the complex flow as long as there is still data in the buffers
    while (flow.GetNstates())
      flow.Execute();

    // Clear the event
    event->Clear();
  }

  // Clearing created pointers
  delete stepper;
  delete propagator;
  for (auto i = 0; i < 4; ++i) delete propagators[i];

  // End of test
  std::cout << "\n--END OF TEST--\n";
  return 0;
}
