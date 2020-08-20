#include <iostream>
#include "CocktailGenerator.h"
#include "Event.h"
#include "Track.h"
#include "Particle.h"
#include "vectorFlow/SystemOfUnits.h"
#include "VecGeom/base/Vector3D.h"
#include <map>
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

  // Generate events
  for (auto i = 0; i < nEvents; i++) {
    CocktailGenerator::Event_t* event = cocktailGen.NextEvent();
    event->SetEvent(i);

    // Map that will count the number of different particles in the event
    std::map<std::string, int> observedParticles;
    for (auto i = 0; i < event->GetNprimaries(); i++) {
      const auto particleId    = (event->GetPrimary(i))->Particle();
      const geantphysics::Particle* particle = geantphysics::Particle::GetParticleByInternalCode(particleId);
      observedParticles[particle->GetName()]++;
    }
    
    // Print event and particles
    event->Print("");
    for (auto it = observedParticles.begin(); it != observedParticles.end(); it++)
      std::cout << "Particle: " << it->first << "\t\tCount: " << it->second << "\n";
    std::cout << "\n";
    if (i == nEvents - 1) { std::cout << "Last event was:\n"; event->Print("ALL"); }

    // Clear event
    event->Clear();
  }

  // End of test
  printf("\n--END OF TEST--\n");
  return 0;
}
