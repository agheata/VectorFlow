/// +--------------------------------------------------------------------------------------+
/// | VectorFlow unit test. Generation + Propagation. --RUNNING IN VECTOR MODE--           |
/// |                                                                                      |
/// | This test is intended to show the usage of the PipelineFlow interface of VectorFlow. |
/// | The example propagates charged particles in a constant field in several steps,       |
/// | until they reach a sphere with a given radius centered in the event vertex.          |
/// |                                                                                      |
/// | The pipeline flow has only one task (TaskPropagator) and the data type to be used    |
/// | are the tracks that are generated in each even of the event-loop.                    |
/// | After each execution, the pipeline is cleared and the next tracks are processed.     |
/// +--------------------------------------------------------------------------------------+

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
#include "vectorFlow/AlignedArray.h"
#include "vectorFlow/VectorTypes.h"
#include "Geant/ConstFieldHelixStepper.h"

using namespace vectorflow;

template <typename Data, typename DataContainer>
struct PropagationSOA {
  std::size_t fSize;
  AlignedArray<double> fPosX   , fPosY    , fPosZ   ;
  AlignedArray<double> fDirX   , fDirY    , fDirZ   ;
  AlignedArray<double> fCharge , fMomentum, fStep   ;
  AlignedArray<double> fNewPosX, fNewPosY , fNewPosZ;
  AlignedArray<double> fNewDirX, fNewDirY , fNewDirZ;

  // Gather info to fill one SOA lane
  void Gather(Data* track, std::size_t lane) {
    fPosX[lane]     = track->Position().x();
    fPosY[lane]     = track->Position().y();
    fPosZ[lane]     = track->Position().z();
    fDirX[lane]     = track->Direction().x();
    fDirY[lane]     = track->Direction().y();
    fDirZ[lane]     = track->Direction().z();
    fCharge[lane]   = track->Charge();
    fMomentum[lane] = track->P();
    fStep[lane]     = track->GetStep();
  }

  // Scatter one lane info back to original structure
  void Scatter(Data* track, std::size_t lane) {
    track->SetPosition(fNewPosX[lane], fNewPosY[lane], fNewPosZ[lane]);
    track->SetDirection(fNewDirX[lane], fNewDirY[lane], fNewDirZ[lane]);
  }

  // Struct constructor
  PropagationSOA(std::size_t size) : fSize(size), fPosX(size), fPosY(size), fPosZ(size), fDirX(size), fDirY(size), fDirZ(size), fCharge(size), fMomentum(size), fStep(size), fNewPosX(size), fNewPosY(size), fNewPosZ(size), fNewDirX(size), fNewDirY(size), fNewDirZ(size) {}
};

struct TaskPropagator : public Work<Track, std::vector<Track*>> {
  // Propagator task needs a stepper
  trackml::HelixPropagator* fPropagator;
  trackml::SimpleStepper*   fStepper;
  ConstFieldHelixStepper*   fHelixStepper;

  // Scalar mode executor
  void Execute(Track* track) {
    // Propagate track to the boundary of a sphere of radius 20cm
    const double radius = 20. * geant::units::cm;
    fStepper->PropagateToR(radius, *track); 
  }

  // Check whether all lanes can still be propagated or not
  bool AllLanesCanDoStep(double* mask, const std::size_t& kVectorSize) {
    for (auto i = 0; i < kVectorSize; i++) 
      if (mask[i] >= 0.) return false;
    return true;
  }

  // PropagateToR function for a vector of tracks
  template <typename Real_v>
  void PropagateToR(std::vector<Track*> const& tracks) {
    using vecCore::Load;
    using vecCore::Store;
    std::size_t tracksSize = tracks.size();
    constexpr std::size_t kVectorSize = vecCore::VectorSize<Real_v>();
    
    // Execute scalar mode if there are less tracks than vector lanes
    if (tracksSize < kVectorSize) {
      for (auto i = 0; i < tracksSize; i++)
        Execute(tracks[i]);
      return;
    }

    PropagationSOA<Track, std::vector<Track*>> mySOA(kVectorSize);

    double mask[kVectorSize];      // Array to check if the track in the lane is fully propagated
    std::size_t lane[kVectorSize]; // Array to check which track is executing in which lane
    std::size_t nextTrack = 0;

    // Declare constants and auxiliar variables
    constexpr double radius      = 20.   * geant::units::cm;
    constexpr double epsilon     = 1.E-4 * geant::units::mm;
    constexpr double toKiloGauss = 1.0   / geant::units::kilogauss;

    vecgeom::Vector3D<double> pos;
    vecgeom::Vector3D<double> dir;
    vecgeom::Vector3D<double> const bfield = fPropagator->GetBfield();
    constexpr double radius2 = radius * radius;
    double rad2, c, safety, dmax, pDotV, d2, snext, step_geom, step_field, step;

    // Set up the first lanes 
    for (auto i = 0; i < kVectorSize; i++) {
      lane[i] = i;
      pos     = tracks[i]->Position();
      dir     = tracks[i]->Direction();
      rad2    = pos.Mag2();
      c       = rad2 - radius2;
      mask[i] = c;
      nextTrack++;
    }

    // Execute vector mode until there are no more tracks to dispatch or 
    // all lanes can still be propagated
    while (nextTrack < tracksSize || AllLanesCanDoStep(mask, kVectorSize)) {
      // Update tracks state and auxiliar variables
      for (auto i = 0; i < kVectorSize; i++) {
        // mask[i] = c >= 0 means fully propagation
        if (mask[i] >= 0. && nextTrack < tracksSize) { 
          lane[i] = nextTrack;
          pos  = tracks[nextTrack]->Position();
          dir  = tracks[nextTrack]->Direction();
          rad2 = pos.Mag2();
          c    = rad2 - radius2;
          mask[i] = c;
          nextTrack++;
        }
        safety = radius - vectorflow::Math::Sqrt(rad2);
        dmax   = 8. * epsilon / tracks[lane[i]]->Curvature(bfield.z() * toKiloGauss);
        pDotV  = pos.Dot(dir);
        d2     = pDotV * pDotV - mask[i];
        snext  = -pDotV + vectorflow::Math::Sqrt(vectorflow::Math::Abs(d2));
        tracks[lane[i]]->SetSafety(safety);
        tracks[lane[i]]->SetSnext(snext);
        step_geom  = vectorflow::Math::Max(epsilon, tracks[lane[i]]->GetSnext());
        step_field = vectorflow::Math::Max(dmax, safety);
        step       = vectorflow::Math::Min(step_geom, step_field);
        tracks[lane[i]]->SetStep(step);

        // Gather info into SOA type
        mySOA.Gather(tracks[lane[i]], i);
      }
    
      Real_v oldPosX_v, oldPosY_v, oldPosZ_v, oldDirX_v, oldDirY_v, oldDirZ_v;
      Real_v newPosX_v, newPosY_v, newPosZ_v, newDirX_v, newDirY_v, newDirZ_v;
      Real_v charge_v, momentum_v, stepSz_v;

      Load(oldPosX_v, mySOA.fPosX.Head());
      Load(oldPosY_v, mySOA.fPosY.Head());
      Load(oldPosZ_v, mySOA.fPosZ.Head());
      Load(oldDirX_v, mySOA.fDirX.Head());
      Load(oldDirY_v, mySOA.fDirY.Head());
      Load(oldDirZ_v, mySOA.fDirZ.Head());
      Load(charge_v, mySOA.fCharge.Head());
      Load(momentum_v, mySOA.fMomentum.Head());
      Load(stepSz_v, mySOA.fStep.Head());

      fHelixStepper->DoStep<Real_v>(oldPosX_v, oldPosY_v, oldPosZ_v, oldDirX_v, oldDirY_v, oldDirZ_v, charge_v, momentum_v, stepSz_v, newPosX_v, newPosY_v, newPosZ_v, newDirX_v, newDirY_v, newDirZ_v);

      Store(newPosX_v, mySOA.fNewPosX.Head());
      Store(newPosY_v, mySOA.fNewPosY.Head());
      Store(newPosZ_v, mySOA.fNewPosZ.Head());
      Store(newDirX_v, mySOA.fNewDirX.Head());
      Store(newDirY_v, mySOA.fNewDirY.Head());
      Store(newDirZ_v, mySOA.fNewDirZ.Head());

      // Update tracks state and auxiliar variables
      for (auto i = 0; i < kVectorSize; i++) {
        // Scatter back info into original structure
        mySOA.Scatter(tracks[lane[i]], i);

        tracks[lane[i]]->DecreasePstep(tracks[lane[i]]->GetStep());
        tracks[lane[i]]->DecreaseSnext(tracks[lane[i]]->GetStep());
        tracks[lane[i]]->DecreaseSafety(tracks[lane[i]]->GetStep());
        tracks[lane[i]]->IncreaseStep(tracks[lane[i]]->GetStep());
        tracks[lane[i]]->IncrementNsteps();

        pos = tracks[lane[i]]->Position();
        dir = tracks[lane[i]]->Direction();
        rad2 = pos.Mag2();
        mask[i] = rad2 - radius2;
      }
    }

    // Execute remaining tracks in scalar mode
    for (auto i = 0; i < kVectorSize; i++) Execute(tracks[lane[i]]);

    // The step was determined following: examples/trackML/src/SimpleStepper.cxx 
  }

  // Vector mode executor
  void Execute(std::vector<Track*> const& tracks) {
    std::cout << "Execution in vector mode.\n";
    PropagateToR<Double_v>(tracks);
    for (const auto& track : tracks) {
      std::cout << "Track " << track->Event() << " made " << track->GetNsteps() << " steps. ";
      std::cout << "Exit position: " << track->Position() << "\n";
    }
  }

  // Struct constructor
  TaskPropagator(trackml::HelixPropagator* propagator, trackml::SimpleStepper* stepper, ConstFieldHelixStepper* helixstepper) : fPropagator(propagator), fStepper(stepper), fHelixStepper(helixstepper) {}
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
    std::cout << "Failed to pass init check.\n";
    exit(-1);
  }

  // Initialize an helix propagator in field and a simple stepper
  vecgeom::Vector3D<double> bfield(0., 0., 20. * geant::units::kilogauss);
  trackml::HelixPropagator* propagator   = new trackml::HelixPropagator(bfield);
  trackml::SimpleStepper*   stepper      = new trackml::SimpleStepper(propagator);
  ConstFieldHelixStepper*   helixstepper = new ConstFieldHelixStepper(bfield);

  // Create a pipeline flow with one stage
  PipelineFlow<Track, std::vector<Track*>, 1> plFlow;

  // Create and initialize the propagator task
  TaskPropagator tPropagate(propagator, stepper, helixstepper);

  // Add the task to the flow
  static constexpr size_t kPropagatorStage = 0;
  plFlow.AddWork(&tPropagate, kPropagatorStage);

  // Set the task to be executed in vector mode
  plFlow.SetVectorMode(kPropagatorStage, true);

  // Event loop
  for (auto i = 0; i < nEvents; i++) {
    CocktailGenerator::Event_t* event = cocktailGen.NextEvent();
    event->SetEvent(i);

    // Add data to the flow
    std::cout << "\n=== Propagating event " << i << "\n";
    for (auto j = 0; j < event->GetNprimaries(); j++) {
      Track* track = event->GetPrimary(j);
      // Only charged tracks will be added
      if (track->Charge() != 0) plFlow.AddData(track);
    }

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
