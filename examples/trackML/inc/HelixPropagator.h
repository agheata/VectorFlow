#ifndef TRACKML_HELIX_PROPAGATOR_H
#define TRACKML_HELIX_PROPAGATOR_H

#include <VecGeom/base/Vector3D.h>

class ConstFieldHelixStepper;

namespace vectorflow {
  class Track;
}

namespace trackml {
class HelixPropagator {
private:
  vecgeom::Vector3D<double>  fBfield;                 ///< Uniform magnetic field
  ConstFieldHelixStepper    *fHelixStepper = nullptr; ///< The helix stepper

public:
  HelixPropagator(vecgeom::Vector3D<double> const &bfield);
  ~HelixPropagator();
  
  /// Getter for the field
  vecgeom::Vector3D<double> GetBfield() const { return fBfield; }

  /// Method that propagates a track in field along a given step distance
  /// Updates only track position and direction
  void Propagate(vectorflow::Track &track, double step) const;
};
} // namespace trackml

#endif
