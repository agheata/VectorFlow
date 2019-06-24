/// -------------------------------------------------------------------------------------+
/// VectorFlow unit test. Generation + Propagation.                                      |
///                                                                                      |
/// This test is intended to show the usage of the PipelineFlow interface of VectorFlow. |
/// The example propagates charged particles in a constant field in several steps,       |
/// until they reach a sphere with a given radius centered in the event vertex.          |
///                                                                                      |
/// The pipeline flow has only one task (TaskPropagator) and the data type to be used    |
/// are the events that are generated in the event-loop.                                 |
/// After each execution, the pipeline is cleared and the next event is processed.       |
/// -------------------------------------------------------------------------------------+

#include <iostream>
#include <string>
#include <vector>
#include <vectorFlow/PipelineFlow.h>
#include "CocktailGenerator.h"
#include "Event.h"
#include "Track.h"
#include "Particle.h"
#include "vectorFlow/SystemOfUnits.h"
#include "base/Vector3D.h"
#include "HelixPropagator.h"
#include "SimpleStepper.h"

using namespace vectorflow;
using Event_t = CocktailGenerator::Event_t;

struct TaskPropagator : public Work<Event_t, std::vector<Event_t*>> {
  // Propagator task needs a stepper
  trackml::SimpleStepper* stepper;

  // Propagate function
  void PropagateTracks(Event_t* event) {
    std::cout << "\n=== Propagating event " << event->GetEvent() << ":\n";

    // Propagate all charged tracks to the boundary of a sphere of radius 20cm
    constexpr double radius = 20. * geant::units::cm;

    for (auto i = 0; i < event->GetNprimaries(); i++) {
      vectorflow::Track* track = event->GetPrimary(i);
      if (track->Charge() != 0) {
        stepper->PropagateToR(radius, *track);
        std::cout << "Track " << i << " made " << track->GetNsteps() << " steps. ";
        std::cout << "Exit position: " << track->Position() << "\n";
      }
    }
  } 

  // Scalar mode executor
  void Execute(Event_t* event) {
    PropagateTracks(event);
  }

  // Vector mode executor
  void Execute(std::vector<Event_t*> const &events) {
    std::cout << "TaskPropagator::Execute not implemented for vector mode\n";
  }

  // Struct constructor
  TaskPropagator(trackml::SimpleStepper* stpr) : stepper(stpr) {}
};

int main(int argc, char* argv[]) {
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
  cocktailGen.SetPrimaryEnergyRange(0.1 * geant::units::GeV, 10 * geant::units::GeV);
  cocktailGen.SetMaxPrimaryPerEvt(100);
  cocktailGen.SetAvgPrimaryPerEvt(70);
  cocktailGen.SetVertex(vertex);
  cocktailGen.SetMaxDepth(2);

  // Check if generator parameters are correct
  if (!cocktailGen.InitPrimaryGenerator()) {
    std::cout << "Failed to pass init check\n";
    exit(-1);
  }

  // Initialize an helix propagator in field and a simple stepper
  vecgeom::Vector3D<double> bfield(0., 0., 20. * geant::units::kilogauss);
  trackml::HelixPropagator* propagator = new trackml::HelixPropagator(bfield);
  trackml::SimpleStepper*   stepper    = new trackml::SimpleStepper(propagator);

  // Create a pipeline flow with one stage
  PipelineFlow<Event_t, std::vector<Event_t*>, 1> plFlow;

  // Create and initialize the propagator task
  TaskPropagator tPropagate(stepper);

  // Add the task to the flow
  static constexpr size_t kPropagatorStage = 0;
  plFlow.AddWork(&tPropagate, kPropagatorStage);

  // Event loop
  for (auto i = 0; i < nEvents; i++) {
    Event_t* event = cocktailGen.NextEvent();
    event->SetEvent(i);

    // Add data to the flow
    plFlow.AddData(event);

    // Process the flow
    plFlow.Execute();

    // Clear pipeline and event
    plFlow.Clear();
    event->Clear();
  }

  // Clearing created pointers
  delete stepper;
  delete propagator;

  // End of test
  std::cout << "\n--END OF TEST--\n";
  return 0;
}
