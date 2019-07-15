#include "SimpleStepper.h"

#include "Track.h"
#include "HelixPropagator.h"
#include "vectorFlow/SystemOfUnits.h"

namespace trackml {

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

} // namespace trackml
