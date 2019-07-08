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
#include <vector>
#include <vectorFlow/PipelineFlow.h>
#include <vecCore/vecCore>
#include "CocktailGenerator.h"
#include "Event.h"
#include "Track.h"
#include "Particle.h"
#include "vectorFlow/SystemOfUnits.h"
#include "base/Vector3D.h"
#include "HelixPropagator.h"
#include "SimpleStepper.h"
#include "vectorFlow/VectorTypes.h"
#include "Geant/ConstFieldHelixStepper.h"

using namespace vectorflow;
using vecCore::Get;
using vecCore::Set;
using vecCore::Load;
using vecCore::Store;
using vecCore::Scalar;

template <typename Data>
struct SIMDTracks {
  // Data in vecCore vector types
  Double_v fPosX_v   , fPosY_v    , fPosZ_v   ;
  Double_v fDirX_v   , fDirY_v    , fDirZ_v   ;
  Double_v fCharge_v , fMomentum_v, fStep_v   ;
  Double_v fNewPosX_v, fNewPosY_v , fNewPosZ_v;
  Double_v fNewDirX_v, fNewDirY_v , fNewDirZ_v;

  // Gather info to fill one SIMD lane
  void Gather(Data* track, std::size_t lane) {
    Set(fPosX_v    , lane, track->Position().x());
    Set(fPosY_v    , lane, track->Position().y());
    Set(fPosZ_v    , lane, track->Position().z());
    Set(fDirX_v    , lane, track->Direction().x());
    Set(fDirY_v    , lane, track->Direction().y());
    Set(fDirZ_v    , lane, track->Direction().z());
    Set(fCharge_v  , lane, track->Charge());
    Set(fMomentum_v, lane, track->P());
    Set(fStep_v    , lane, track->GetStep());
  }

  // Scatter one lane info back to original structure
  void Scatter(Data* track, std::size_t lane) {
    track->SetPosition(static_cast<Scalar<Double_v>>(Get(fNewPosX_v, lane)), static_cast<Scalar<Double_v>>(Get(fNewPosY_v, lane)), static_cast<Scalar<Double_v>>(Get(fNewPosZ_v, lane)));
    track->SetDirection(static_cast<Scalar<Double_v>>(Get(fNewDirX_v, lane)), static_cast<Scalar<Double_v>>(Get(fNewDirY_v, lane)), static_cast<Scalar<Double_v>>(Get(fNewDirZ_v, lane)));
  }
};

struct TaskPropagator : public Work<Track, std::vector<Track*>> {
  // Propagator task needs a stepper
  trackml::HelixPropagator* fPropagator;
  trackml::SimpleStepper*   fStepper;
  ConstFieldHelixStepper*   fHelixStepper;
  SIMDTracks<Track> fSIMDStruct;

  // Scalar mode executor
  void Execute(Track* track) {
    // Propagate track to the boundary of a sphere of radius 20cm
    const double radius = 20. * geant::units::cm;
    fStepper->PropagateToR(radius, *track); 
  }

  // function to check whether all lanes can still be propagated or not
  bool AllLanesCanDoStep(double* cvalues, const std::size_t& kVectorSize) {
    for (auto i = 0; i < kVectorSize; i++) 
      if (cvalues[i] >= 0.) return false;
    return true;
  }

  // PropagateToR function for a vector of tracks
  void PropagateToR(const double& radius, std::vector<Track*> const& tracks) {
    const std::size_t kTracksSize = tracks.size();
    const std::size_t kVectorSize = vecCore::VectorSize<Double_v>();
    
    // Execute scalar mode if there are less tracks than vector lanes
    if (kTracksSize < kVectorSize) {
      for (auto i = 0; i < kTracksSize; i++)
        Execute(tracks[i]);
      return;
    }

    double cvalues[kVectorSize];   // Auxiliar array to check if the track in the lane is fully propagated
    std::size_t lane[kVectorSize]; // Auxiliar array to check which track is executing in which lane
    std::size_t nextTrack = 0;     // Auxiliar counter to keep record of which track will be dispatched 

    // Declare constants and auxiliar variables
    const double epsilon     = 1.E-4 * geant::units::mm;
    const double toKiloGauss = 1.0   / geant::units::kilogauss;

    vecgeom::Vector3D<double> pos;
    vecgeom::Vector3D<double> dir;
    vecgeom::Vector3D<double> const bfield = fPropagator->GetBfield();
    const double radius2 = radius * radius;
    double rad2, c, safety, dmax, pDotV, d2, snext, step_geom, step_field, step;

    // Set up the first lanes 
    for (auto i = 0; i < kVectorSize; i++) {
      lane[i]    = i;
      pos        = tracks[i]->Position();
      dir        = tracks[i]->Direction();
      rad2       = pos.Mag2();
      c          = rad2 - radius2;
      cvalues[i] = c;
      nextTrack++;
    }

    // Execute in vector mode until there are no more tracks to dispatch
    // or all lanes in SIMD struct can still be propagated
    while (nextTrack < kTracksSize || AllLanesCanDoStep(cvalues, kVectorSize)) {
      // Calculate track step
      for (auto i = 0; i < kVectorSize; i++) {
        // cvalues[i] = c >= 0 means fully propagation
        if (cvalues[i] >= 0. && nextTrack < kTracksSize) { 
          // Assign new track to a lane 
          lane[i]    = nextTrack;
          pos        = tracks[nextTrack]->Position();
          dir        = tracks[nextTrack]->Direction();
          rad2       = pos.Mag2();
          c          = rad2 - radius2;
          cvalues[i] = c;
          nextTrack++;
        }
        safety = radius - vecCore::math::Sqrt(rad2);
        dmax   = 8. * epsilon / tracks[lane[i]]->Curvature(bfield.z() * toKiloGauss);
        pDotV  = pos.Dot(dir);
        d2     = pDotV * pDotV - cvalues[i];
        snext  = -pDotV + vecCore::math::Sqrt(vecCore::math::Abs(d2));
        
        tracks[lane[i]]->SetSafety(safety);
        tracks[lane[i]]->SetSnext(snext);
        
        step_geom  = vecCore::math::Max(epsilon, snext);
        step_field = vecCore::math::Max(dmax, safety);
        step       = vecCore::math::Min(step_geom, step_field);
        
        tracks[lane[i]]->SetStep(step);

        // Gather info into SIMD struct
        fSIMDStruct.Gather(tracks[lane[i]], i);
      }

      fHelixStepper->DoStep<Double_v>(fSIMDStruct.fPosX_v, fSIMDStruct.fPosY_v, fSIMDStruct.fPosZ_v, fSIMDStruct.fDirX_v, fSIMDStruct.fDirY_v, fSIMDStruct.fDirZ_v, fSIMDStruct.fCharge_v, fSIMDStruct.fMomentum_v, fSIMDStruct.fStep_v, fSIMDStruct.fNewPosX_v, fSIMDStruct.fNewPosY_v, fSIMDStruct.fNewPosZ_v, fSIMDStruct.fNewDirX_v, fSIMDStruct.fNewDirY_v, fSIMDStruct.fNewDirZ_v);

      // Update tracks state and auxiliar variables
      for (auto i = 0; i < kVectorSize; i++) {
        // Scatter back info into original structure
        fSIMDStruct.Scatter(tracks[lane[i]], i);

        tracks[lane[i]]->DecreasePstep(tracks[lane[i]]->GetStep());
        tracks[lane[i]]->DecreaseSnext(tracks[lane[i]]->GetStep());
        tracks[lane[i]]->DecreaseSafety(tracks[lane[i]]->GetStep());
        tracks[lane[i]]->IncreaseStep(tracks[lane[i]]->GetStep());
        tracks[lane[i]]->IncrementNsteps();
        
        pos = tracks[lane[i]]->Position();
        dir = tracks[lane[i]]->Direction();

        // Check direction is normalized
        assert(dir.IsNormalized() && "ERROR: Direction not normalized after field propagation");

        rad2       = pos.Mag2();
        c          = rad2 - radius2;
        cvalues[i] = rad2 - radius2;
      }
    }

    // Execute remaining tracks in scalar mode
    for (auto i = 0; i < kVectorSize; i++) Execute(tracks[lane[i]]);

    // The step was determined following: examples/trackML/src/SimpleStepper.cxx 
  }

  // Vector mode executor
  void Execute(std::vector<Track*> const& tracks) {
    std::cout << "--EXECUTION IN VECTOR MODE--\n";
    // Propagate track to the boundary of a sphere of radius 20cm
    const double radius = 20. * geant::units::cm;
    PropagateToR(radius, tracks);
    for (const auto& track : tracks) {
      std::cout << "Track " << track->Event() << " made " << track->GetNsteps() << " steps. ";
      std::cout << "Exit position: " << track->Position() << "\n";
    }
  }

  // Task struct constructor
  TaskPropagator(trackml::HelixPropagator* propagator, trackml::SimpleStepper* stepper, ConstFieldHelixStepper* helixstepper, SIMDTracks<Track> simd) : fPropagator(propagator), fStepper(stepper), fHelixStepper(helixstepper), fSIMDStruct(simd) {}
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

  // Create SIMD struct to handle vectorization
  SIMDTracks<Track> mySIMDStruct;

  // Create a pipeline flow with one stage
  PipelineFlow<Track, std::vector<Track*>, 1> plFlow;

  // Create and initialize the propagator task
  TaskPropagator tPropagate(propagator, stepper, helixstepper, mySIMDStruct);

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

  // Clear created pointers
  delete stepper;
  delete propagator;
  delete helixstepper;

  // End of test
  std::cout << "\n--END OF TEST--\n";
  return 0;
}
