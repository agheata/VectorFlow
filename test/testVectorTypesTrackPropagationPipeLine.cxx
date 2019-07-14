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
#include <cassert>
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
struct Tracks_v {
  // Data in vecCore vector types
  vecgeom::Vector3D<Double_v> fPos_v;
  vecgeom::Vector3D<Double_v> fDir_v;
  Double_v fCharge_v , fMomentum_v, fNSteps_v ;

  // Auxiliar vectorized methods
  Double_v Pt_v() const { return fMomentum_v * fDir_v.Perp(); }

  Double_v Curvature_v(double Bz) { 
    Double_v qB_v = fCharge_v * Bz; 
    Double_v curvature_v = vecCore::math::Abs(Track::kB2C * qB_v / (Pt_v() +
          Track::kTiny)); 
    return vecCore::Blend<Double_v>(vecCore::math::Abs(qB_v) < Track::kTiny,
        Track::kTiny, curvature_v);
  }

  // Gather info to fill one SIMD lane
  void Gather(Data* track, const std::size_t lane) {
    Set(fPos_v.x() , lane, track->Position().x());
    Set(fPos_v.y() , lane, track->Position().y());
    Set(fPos_v.z() , lane, track->Position().z());
    Set(fDir_v.x() , lane, track->Direction().x());
    Set(fDir_v.y() , lane, track->Direction().y());
    Set(fDir_v.z() , lane, track->Direction().z());
    Set(fCharge_v  , lane, track->Charge());
    Set(fMomentum_v, lane, track->P());
    Set(fNSteps_v  , lane, 0);
  }

  // Scatter one lane info back to original structure
  void Scatter(const std::size_t lane, Data* track) {
    track->SetPosition(Get(fPos_v.x(), lane),
        Get(fPos_v.y(), lane),
        Get(fPos_v.z(), lane));
    track->SetDirection(Get(fDir_v.x(), lane),
        Get(fDir_v.y(), lane),
        Get(fDir_v.z(), lane));
    track->SetNsteps(Get(fNSteps_v, lane));
  }
};

struct TaskPropagator : public Work<Track, std::vector<Track*>> {
  // Propagator task needs a stepper
  trackml::HelixPropagator* fPropagator;
  trackml::SimpleStepper*   fStepper;
  ConstFieldHelixStepper*   fHelixStepper;
  Tracks_v<Track> fTracks_v;

  // Scalar mode executor
  void Execute(Track* track) {
    // Propagate track to the boundary of a sphere of radius 20cm
    const double radius = 20. * geant::units::cm;
    fStepper->PropagateToR(radius, *track); 
  }

  // Function to check whether all lanes can still be propagated or not
  bool AllLanesCanDoStep(const Double_v& c_v) {
    return vecCore::MaskFull(c_v < 0);
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
    
    std::size_t lane[kVectorSize]; // Array to check which track is which lane
    std::size_t nextTrack = 0;     // Counter to dispatch tracks 

    // Declare constants and auxiliar variables
    const double epsilon     = 1.E-4 * geant::units::mm;
    const double toKiloGauss = 1.0   / geant::units::kilogauss;

    vecgeom::Vector3D<double> const bfield = fPropagator->GetBfield();
    const double radius2 = radius * radius;
    Double_v rad2_v, c_v, safety_v, dmax_v, pDotV_v, d2_v, snext_v,
             step_geom_v, step_field_v, step_v;

    // Set up first lanes
    for (auto i = 0; i < kVectorSize; i++) {
      lane[i] = i;
      fTracks_v.Gather(tracks[i], i);
    }
    nextTrack += kVectorSize;
    rad2_v     = fTracks_v.fPos_v.Mag2();
    c_v        = rad2_v - radius2;

    // Execute until all lanes can at least do one more step
    while (AllLanesCanDoStep(c_v)) {      
      safety_v = radius - vecCore::math::Sqrt(rad2_v);
      dmax_v   = 8. * epsilon / fTracks_v.Curvature_v(bfield.z() * toKiloGauss);
      pDotV_v  = fTracks_v.fPos_v.Dot(fTracks_v.fDir_v);
      d2_v     = pDotV_v * pDotV_v - c_v;
      snext_v  = -pDotV_v + vecCore::math::Sqrt(vecCore::math::Abs(d2_v));

      step_geom_v  = vecCore::math::Max(static_cast<Double_v>(epsilon), snext_v);
      step_field_v = vecCore::math::Max(dmax_v, safety_v);
      step_v       = vecCore::math::Min(step_geom_v, step_field_v);

      fHelixStepper->DoStep<Double_v>(fTracks_v.fPos_v, fTracks_v.fDir_v,
          fTracks_v.fCharge_v, fTracks_v.fMomentum_v, step_v, fTracks_v.fPos_v,
          fTracks_v.fDir_v);
      
      assert(fTracks_v.fDir_v.IsNormalized() && 
          "ERROR: Direction not normalized after field propagation");

      rad2_v = fTracks_v.fPos_v.Mag2();
      c_v    = rad2_v - radius2;

      fTracks_v.fNSteps_v += 1;

      // Assign new track to a lane if previous was fully propagated
      if (!AllLanesCanDoStep(c_v)) {
        for (auto i = 0; i < kVectorSize; i++) {
          if (Get(c_v, i) >= 0 && nextTrack < kTracksSize) {
            // Scatter back the propagated track
            fTracks_v.Scatter(i, tracks[lane[i]]);

            // Set values for new track in its corresponding lane
            fTracks_v.Gather(tracks[nextTrack], i);
            lane[i] = nextTrack++;
          }
        }
        // Calculate c_v values again with new tracks
        rad2_v = fTracks_v.fPos_v.Mag2();
        c_v    = rad2_v - radius2;
      }
    }

    // Execute remaining tracks in scalar mode
    for (auto i = 0; i < kVectorSize; i++) Execute(tracks[lane[i]]);
  }

  // Vector mode executor
  void Execute(std::vector<Track*> const& tracks) {
    std::cout << "--EXECUTION IN VECTOR MODE--\n";
    // Propagate track to the boundary of a sphere of radius 20cm
    const double radius = 20. * geant::units::cm;
    PropagateToR(radius, tracks);
    for (const auto& track : tracks) {
      std::cout << "Track " << track->Event() << " made " << track->GetNsteps()
        << " steps. ";
      std::cout << "Exit position: " << track->Position() << "\n";
    }
  }

  // Task struct constructor
  TaskPropagator(trackml::HelixPropagator* propagator, trackml::SimpleStepper*
      stepper, ConstFieldHelixStepper* helixstepper, Tracks_v<Track> simd) :
    fPropagator(propagator), fStepper(stepper), fHelixStepper(helixstepper),
    fTracks_v(simd) {}
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
  Tracks_v<Track> myTracks_v;

  // Create a pipeline flow with one stage
  PipelineFlow<Track, std::vector<Track*>, 1> plFlow;

  // Create and initialize the propagator task
  TaskPropagator tPropagate(propagator, stepper, helixstepper, myTracks_v);

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
