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
  Double_v fCharge_v , fMomentum_v, fNSteps_v ;
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
    Set(fNSteps_v  , lane, 0);
  }

  // Scatter one lane info back to original structure
  void Scatter(Data* track, std::size_t lane) {
    track->SetPosition(static_cast<Scalar<Double_v>>(Get(fNewPosX_v, lane)),
        static_cast<Scalar<Double_v>>(Get(fNewPosY_v, lane)),
        static_cast<Scalar<Double_v>>(Get(fNewPosZ_v, lane)));
    track->SetDirection(static_cast<Scalar<Double_v>>(Get(fNewDirX_v, lane)),
        static_cast<Scalar<Double_v>>(Get(fNewDirY_v, lane)),
        static_cast<Scalar<Double_v>>(Get(fNewDirZ_v, lane)));
    track->SetNsteps(static_cast<Scalar<Double_v>>(Get(fNSteps_v, lane)));
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
  bool AllLanesCanDoStep(Double_v c_v, const std::size_t& kVectorSize) {
    for (auto i = 0; i < kVectorSize; i++)
      if (static_cast<Scalar<Double_v>>(Get(c_v, i)) >= 0.) return false;
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
    
    std::size_t lane[kVectorSize]; // Array to check which track is executing in which lane
    std::size_t nextTrack = 0;     // Counter to keep record of which track will be dispatched 

    // Declare constants and auxiliar variables
    const double epsilon     = 1.E-4 * geant::units::mm;
    const double toKiloGauss = 1.0   / geant::units::kilogauss;

    vecgeom::Vector3D<Double_v> pos_v;
    vecgeom::Vector3D<Double_v> dir_v;
    vecgeom::Vector3D<double> const bfield = fPropagator->GetBfield();
    const double radius2 = radius * radius;
    Double_v rad2_v, c_v, safety_v, dmax_v, pDotV_v, d2_v, snext_v,
             step_geom_v, step_field_v, step_v;

    // Set up first lanes
    for (auto i = 0; i < kVectorSize; i++) {
      lane[i] = i;
      fSIMDStruct.Gather(tracks[i], i);
    }
    nextTrack += kVectorSize;

    pos_v.Set(fSIMDStruct.fPosX_v, fSIMDStruct.fPosY_v, fSIMDStruct.fPosZ_v);
    dir_v.Set(fSIMDStruct.fDirX_v, fSIMDStruct.fDirY_v, fSIMDStruct.fDirZ_v);

    rad2_v = pos_v.Mag2();
    c_v    = rad2_v - radius2;

    // Execute in vector mode until there are no more tracks to dispatch
    // or all lanes in SIMD struct can still do at least one more step
    while (nextTrack < kTracksSize || AllLanesCanDoStep(c_v, kVectorSize)) {
      // Assign new track to a lane
      for (auto i = 0; i < kVectorSize; i++) {
        if (static_cast<Scalar<Double_v>>(Get(c_v, i)) >= 0. && nextTrack < kTracksSize) {
          // Scatter back the propagated track
          fSIMDStruct.Scatter(tracks[lane[i]], i);

          // Set values for new track in its corresponding lane
          lane[i] = nextTrack;
          Set(pos_v.x(), i, tracks[nextTrack]->Position().x());
          Set(pos_v.y(), i, tracks[nextTrack]->Position().y());
          Set(pos_v.z(), i, tracks[nextTrack]->Position().z());
          Set(dir_v.x(), i, tracks[nextTrack]->Direction().x());
          Set(dir_v.y(), i, tracks[nextTrack]->Direction().y());
          Set(dir_v.z(), i, tracks[nextTrack]->Direction().z());
          Set(c_v, i, tracks[nextTrack]->Position().Mag2() - radius2);

          // Gather rest of track info
          // However, it also copies unusued data like fPos's and fDir's
          // At this point, there is room for optimization since the only
          // data required here will be the charge, momentum, and steps
          fSIMDStruct.Gather(tracks[nextTrack], i);
          nextTrack++;
        }
      }
      safety_v = radius - vecCore::math::Sqrt(rad2_v);
      // Get curvature for each lane
      for (auto i = 0; i < kVectorSize; i++) Set(dmax_v, i, 8. * epsilon /
          tracks[lane[i]]->Curvature(bfield.z() * toKiloGauss));
      pDotV_v = pos_v.Dot(dir_v);
      d2_v    = pDotV_v * pDotV_v - c_v;
      snext_v = -pDotV_v + vecCore::math::Sqrt(vecCore::math::Abs(d2_v));

      Double_v epsilon_v = epsilon;
      step_geom_v  = vecCore::math::Max(epsilon_v, snext_v);
      step_field_v = vecCore::math::Max(dmax_v, safety_v);
      step_v       = vecCore::math::Min(step_geom_v, step_field_v);

      fHelixStepper->DoStep<Double_v>(pos_v.x(), pos_v.y(), pos_v.z(),
          dir_v.x(), dir_v.y(), dir_v.z(), fSIMDStruct.fCharge_v,
          fSIMDStruct.fMomentum_v, step_v, fSIMDStruct.fNewPosX_v,
          fSIMDStruct.fNewPosY_v, fSIMDStruct.fNewPosZ_v,
          fSIMDStruct.fNewDirX_v, fSIMDStruct.fNewDirY_v,
          fSIMDStruct.fNewDirZ_v);

      pos_v.Set(fSIMDStruct.fNewPosX_v, fSIMDStruct.fNewPosY_v, fSIMDStruct.fNewPosZ_v);
      dir_v.Set(fSIMDStruct.fNewDirX_v, fSIMDStruct.fNewDirY_v, fSIMDStruct.fNewDirZ_v);
      
      assert(dir_v.IsNormalized() && "ERROR: Direction not normalized after field propagation");

      rad2_v = pos_v.Mag2();
      c_v    = rad2_v - radius2;

      fSIMDStruct.fNSteps_v += 1;
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
  TaskPropagator(trackml::HelixPropagator* propagator, trackml::SimpleStepper*
      stepper, ConstFieldHelixStepper* helixstepper, SIMDTracks<Track> simd) :
    fPropagator(propagator), fStepper(stepper), fHelixStepper(helixstepper),
    fSIMDStruct(simd) {}
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
