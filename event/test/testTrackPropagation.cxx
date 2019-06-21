#include <iostream>
#include "CocktailGenerator.h"
#include "Event.h"
#include "Track.h"
#include "Particle.h"
#include "vectorFlow/SystemOfUnits.h"
#include "base/Vector3D.h"
#include "HelixPropagator.h"
#include "SimpleStepper.h"
#include <string>

int main(int argc, char* argv[]) {
  // Read number of events to be generated
  int nEvents = argv[1] == NULL ? 10 : atoi(argv[1]);

  // Add CocktailGenerator
  vecgeom::Vector3D<double> vertex(0., 0., 10.);
  CocktailGenerator cocktailGen(vertex);

  // Add particle species
  cocktailGen.AddPrimary("pi+", 0.3);
  cocktailGen.AddPrimary("pi-", 0.3);
  cocktailGen.AddPrimary("gamma", 0.4);

  // Set parameters
  cocktailGen.SetPrimaryEnergyRange(0.1 * geant::units::GeV, 10. * geant::units::GeV);
  cocktailGen.SetMaxPrimaryPerEvt(100);
  cocktailGen.SetAvgPrimaryPerEvt(70);
  cocktailGen.SetVertex(vertex);
  cocktailGen.SetMaxDepth(2);

  // Check if parameters are correct
  if (!cocktailGen.InitPrimaryGenerator()) { printf("Failed to pass init check.\n"); exit(-1); }

  // Initialize a helix propagator in field and a simple stepper
  vecgeom::Vector3D<double> bfield(0., 0., 20. * geant::units::kilogauss);
  trackml::HelixPropagator *propagator = new trackml::HelixPropagator(bfield);
  trackml::SimpleStepper *stepper = new trackml::SimpleStepper(propagator);

  // Generate events
  for (auto i = 0; i < nEvents; i++) {
    CocktailGenerator::Event_t* event = cocktailGen.NextEvent();
    event->SetEvent(i);
    
    // Propagate all charged tracks to the boundary of a sphere of radius 20 cm
    constexpr double radius = 20. * geant::units::cm;

    std::cout << "=== Propagating event " << i << std::endl;
    for (auto j = 0; j < event->GetNprimaries(); j++) {
      vectorflow::Track *track = event->GetPrimary(j);
      if (track->Charge() != 0) {
        stepper->PropagateToR(radius, *track);
        std::cout << "track " << j << " made " << track->GetNsteps()
                  << " steps, exit position: " << track->Position() << std::endl;
      }
    }

  } // end of loop over events

  delete stepper;
  delete propagator;

  // End of test
  printf("\n--END OF TEST--\n");
  return 0;
}
