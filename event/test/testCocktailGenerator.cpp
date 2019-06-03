#include <iostream>
#include "CocktailGenerator.h"
#include "Event.h"
#include "Track.h"
#include "Particle.h"
#include "vectorFlow/PhysicalConstants.h"
#include "vectorFlow/SystemOfUnits.h"
#include "base/Vector3D.h"

#define N_EVENTS 10

int main(void) {
    // Add CocktailGenerator
    vecgeom::Vector3D<double> vertex(0., 0., 10.);
    CocktailGenerator* CG = new CocktailGenerator(vertex);

    // Add particle species
    CG->AddPrimary("pi+", 0.3);
    CG->AddPrimary("pi-", 0.3);
    CG->AddPrimary("gamma", 0.4);

    // Set parameters
    CG->SetPrimaryEnergyRange(0.1 * geant::units::GeV, 10. * geant::units::GeV);
    CG->SetMaxPrimaryPerEvt(100);
    CG->SetAvgPrimaryPerEvt(70);
    CG->SetVertex(vertex);
    CG->SetMaxDepth(2);

    // Check if parameters are correct
    if (CG->InitPrimaryGenerator()) printf("\nCocktail generator init checked.\n\n");

    // Generate events
    for (size_t i = 0; i < N_EVENTS; i++) {
        CocktailGenerator::Event_t* event = CG->NextEvent();
        event->SetEvent(i);
        event->Print("");
        if (i == N_EVENTS - 1) { printf("\nLast event was:\n"); event->Print("ALL"); }
        event->Clear();
    }

    // End of test
    printf("\n--END OF TEST--\n");
    return 0;
}
