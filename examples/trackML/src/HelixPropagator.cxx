#include "HelixPropagator.h"
#include "Track.h"
#include "Geant/ConstFieldHelixStepper.h"

namespace trackml {

HelixPropagator::HelixPropagator(vecgeom::Vector3D<double> const &bfield) : fBfield(bfield)
{
  fHelixStepper = new ConstFieldHelixStepper(bfield);
}

HelixPropagator::~HelixPropagator()
{
  delete fHelixStepper;
}

void HelixPropagator::Propagate(vectorflow::Track &track, double step) const
{
  using Vector3 = vecgeom::Vector3D<double>;
  Vector3 newPosition, newDirection;

  fHelixStepper->DoStep<double>(track.Position(), track.Direction(), track.Charge(), track.P(), step,
                        newPosition, newDirection);

  track.SetPosition(newPosition);
  track.SetDirection(newDirection);
}

} // namespace trackml
