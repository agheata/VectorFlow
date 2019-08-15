#include "SimpleStepper.h"
#include <management/GeoManager.h>
#include <volumes/Tube.h>
#include "Track.h"
#include "HelixPropagator.h"
#include "vectorFlow/SystemOfUnits.h"
#include "Geant/ConstFieldHelixStepper.h"
#include "Tracks_v.h"

namespace trackml {

SimpleStepper::SimpleStepper(HelixPropagator *prop) : fPropagator(prop)
{
  using namespace vecgeom;
  fLayers.reserve(4);
  fLayers.push_back(GeoManager::MakeInstance<UnplacedTube>(0., kPixVolRmin, kWorldDz, 0., kTwoPi));
  fLayers.push_back(GeoManager::MakeInstance<UnplacedTube>(kPixVolRmin, kPixVolRmax, kWorldDz, 0., kTwoPi));
  fLayers.push_back(GeoManager::MakeInstance<UnplacedTube>(kPixVolRmax, kStrips1VolRmax, kWorldDz, 0., kTwoPi));
  fLayers.push_back(GeoManager::MakeInstance<UnplacedTube>(kStrips1VolRmax, kStrips2VolRmax, kWorldDz, 0., kTwoPi));
}

void SimpleStepper::PropagateToR(const double radius, vectorflow::Track &track) const
{
  // Compute safe distance for the track to getting outside the sphere of given radius
  using namespace vecgeom;
  constexpr double tolerance = 1.E-9 * geant::units::mm;
  constexpr double epsilon = 1.E-4 * geant::units::mm;
  constexpr double toKiloGauss = 1.0 / geant::units::kilogauss; // Converts to kilogauss
  
  Vector3D<double> pos = track.Position();
  Vector3D<double> dir = track.Direction();
  Vector3D<double> const bfield = fPropagator->GetBfield();
  const double radius2 = radius * radius;
  double rad2 = pos.Mag2();
  double c = rad2 - radius2;

  // Enter the stepping loop while the track is inside the sphere of given radius
  while (c < 0) {
    double safety = radius - Sqrt(rad2);
    // Maximum allowed distance so that the sagitta of the track trajectory along this distance
    // is less than a fraction of epsilon of the distance
    double dmax = 8. * epsilon / track.Curvature(bfield.z() * toKiloGauss);

    // Compute distance along straight line to exit the sphere
    double pDotV   = pos.Dot(dir);
    double d2      = pDotV * pDotV - c;
    double snext   = -pDotV + Sqrt(vecCore::math::Abs(d2));

    track.SetSafety(safety);
    track.SetSnext(snext);

    double step_geom = Max(epsilon, track.GetSnext());
    double step_field = Max(dmax, safety);
    double step = Min(step_geom, step_field);

    // Propagate in field with step distance (update position and direction)
    fPropagator->Propagate(track, step);

    // Update track state other than position/direction
    track.DecreasePstep(step);
    track.DecreaseSnext(step);
    track.DecreaseSafety(step);
    track.IncreaseStep(step);
    track.IncrementNsteps();

    // Update local variables
    pos = track.Position();
    dir = track.Direction();
    rad2 = pos.Mag2();
    c = rad2 - radius2;
  }
}

void SimpleStepper::PropagateInTube(int layer, vectorflow::Track &track) const
{
  // Propagate along a helix inside a tube 
  using namespace vecgeom;
  using namespace vecCore::math;
  constexpr double tolerance = 1.E-9 * geant::units::mm;
  constexpr double epsilon = 1.E-4 * geant::units::mm;
  constexpr double toKiloGauss = 1.0 / geant::units::kilogauss; // Converts to kilogauss
  const auto tube = fLayers[layer];

  Vector3D<double> const bfield = fPropagator->GetBfield();

  bool inside = true;
  // Keep propagating while the particle is still inside the tube
  while (inside) {
    Vector3D<double> const &pos = track.Position();
    Vector3D<double> const &dir = track.Direction();
    double safety = tube->SafetyToOut(pos);
    // Maximum allowed distance so that the sagitta of the track trajectory along this distance
    // is less than a fraction of epsilon of the distance
    double dmax = 8. * epsilon / track.Curvature(bfield.z() * toKiloGauss);
    // Compute distance along straight line to exit the tube
    double snext = tube->DistanceToOut(pos, dir);

    track.SetSafety(safety);
    track.SetSnext(snext);

    double step_geom = Max(epsilon, snext);
    double step_field = Max(dmax, safety);
    double step = Min(step_geom, step_field);

    // Propagate in field with step distance (update position and direction)
    fPropagator->Propagate(track, step);

    // Update track state other than position/direction
    track.DecreasePstep(step);
    track.DecreaseSnext(step);
    track.DecreaseSafety(step);
    track.IncreaseStep(step);
    track.IncrementNsteps();

    inside = (step <= safety) ? true : tube->Contains(track.Position());
  }

  // Track is now on a boundary
  track.SetStatus(vectorflow::kBoundary);
  
  // Exiting setup?
  if (Abs(track.Position().z()) > tube->z() ||
      track.Position().Perp() > kStrips2VolRmax) {
    track.SetStatus(vectorflow::kKilled);
  }
}

void SimpleStepper::PropagateInTube(int layer, std::vector<vectorflow::Track*> const &tracks) const
{
  const std::size_t kTracksSize = tracks.size();
  const std::size_t kVectorSize = vecCore::VectorSize<Double_v>();

  // Execute in scalar mode if there are less tracks than vector lanes
  if (kTracksSize < kVectorSize) {
    for (auto track : tracks)
      PropagateInTube(layer, *track);
    return;
  }

  // VECTOR TYPES IMPLEMENTATION
  Tracks_v<Track> tracks_v;  

  std::size_t lane[kVectorSize]; // Array to check which track is in which lane
  std::size_t nextTrack = 0;     // Counter to keep track and dispatch tracks

  using namespace vecgeom;
  using namespace vecCore::math;
  using vecCore::Set;
  using vecCore::Get;

  const double tolerance   = 1.E-9 * geant::units::mm;
  const double epsilon     = 1.E-4 * geant::units::mm;
  const double toKiloGauss = 1.0   / geant::units::kilogauss; // Converts to kilogauss
  const auto   tube = fLayers[layer];

  Vector3D<double> const bfield = fPropagator->GetBfield();
  ConstFieldHelixStepper *helixStepper = new ConstFieldHelixStepper(bfield);

  // Setup first lanes
  for (auto i = 0; i < kVectorSize; i++) {
    lane[i] = i;
    tracks_v.Gather(tracks[i], i);
  }
  nextTrack += kVectorSize;

  Double_v safety_v, dmax_v, snext_v, step_geom_v, step_field_v, step_v;

  bool ongoing = true;
  while (ongoing) {
    safety_v = tube->SafetyToOutVec(tracks_v.fPos_v);
    dmax_v   = 8. * epsilon / tracks_v.Curvature_v(bfield.z() * toKiloGauss);
    snext_v  = tube->DistanceToOutVec(tracks_v.fPos_v, tracks_v.fDir_v, vecgeom::kInfLength);
    
    step_geom_v  = Max(static_cast<Double_v>(epsilon), snext_v);
    step_field_v = Max(dmax_v, safety_v);
    step_v       = Min(step_geom_v, step_field_v);
    
    helixStepper->DoStep<Double_v>(
        tracks_v.fPos_v, tracks_v.fDir_v, tracks_v.fCharge_v,
        tracks_v.fMomentum_v, step_v, tracks_v.fPos_v, tracks_v.fDir_v);

    tracks_v.fNSteps_v += 1;
    
    // Assign new track to a lane if previous was fully propagated
    int propagatedCount = 0;
    bool propagatedTrack[kVectorSize] = {false};
    if (vecCore::MaskFull(step_v <= safety_v)) continue;
    for (auto i = 0; i < kVectorSize; i++) {
      vecgeom::Vector3D<double> pos(Get(tracks_v.fPos_v.x(), i),
                                    Get(tracks_v.fPos_v.y(), i), 
                                    Get(tracks_v.fPos_v.z(), i));
      bool inside = tube->Contains(pos);
      propagatedTrack[i] = !inside;
      if (propagatedTrack[i]) propagatedCount++;
    }
    bool refill = (kTracksSize - nextTrack) > kVectorSize;
    if (!refill) refill = (kTracksSize - nextTrack) > propagatedCount;
    if (refill) {
      // There are enough tracks to refill the done lanes
      for (auto i = 0; i < kVectorSize; i++) {
        if (propagatedTrack[i]) {
          // Scatter back the propagated track
          tracks_v.Scatter(i, tracks[lane[i]]);
          
          // Track is now on boundary
          tracks[lane[i]]->SetStatus(vectorflow::kBoundary);

          // Exiting setup?
          if (Abs(tracks[lane[i]]->Position().z()) > tube->z() || 
              tracks[lane[i]]->Position().Perp() > kStrips2VolRmax) {
            tracks[lane[i]]->SetStatus(vectorflow::kKilled);
          }

          // Set values for new track in its corresponding lane
          tracks_v.Gather(tracks[nextTrack], i);
          lane[i] = nextTrack++;
        }
      }
    } else {
      // Not enough tracks to refill the vector mode, scatter existing ones
      for (auto i = 0; i < kVectorSize; i++) {
        tracks_v.Scatter(i, tracks[lane[i]]);
        if (!propagatedTrack[i]) {
          // Track is not fully propagated yet, propagate in scalar mode
          PropagateInTube(layer, *tracks[lane[i]]);
        } else {
          // Remaining track got propagated, check the boundary
          // Track is now on boundary
          tracks[lane[i]]->SetStatus(vectorflow::kBoundary);

          // Exiting setup?
          if (Abs(tracks[lane[i]]->Position().z()) > tube->z() || 
              tracks[lane[i]]->Position().Perp() > kStrips2VolRmax) {
            tracks[lane[i]]->SetStatus(vectorflow::kKilled);
          }
        }
      }
      ongoing = false;
    }
  }

  // Execute remaining tracks in scalar mode
  while (nextTrack < kTracksSize) 
    PropagateInTube(layer, *tracks[nextTrack++]);

  delete helixStepper;
}

} // namespace trackml
