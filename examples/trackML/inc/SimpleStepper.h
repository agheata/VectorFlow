#ifndef TRACKML_SIMPLE_STEPPER_H
#define TRACKML_SIMPLE_STEPPER_H

#include <vector>
#include <vectorFlow/Typedefs.h>
#include "GeomData.h"

namespace vectorflow {
  class Track;
}

namespace vecgeom {
  namespace VECGEOM_IMPL_NAMESPACE {
    class UnplacedTube;
  }
}

namespace trackml {
/// Very simple stepper that takes a track from vertex and propagates it untill a given distance to the vertex.
/// The propagator makes sure that the error for the last step is lesser than an epsilon value.
class HelixPropagator;

class SimpleStepper {
private:
  HelixPropagator *fPropagator    = nullptr;
  std::vector<vecgeom::UnplacedTube *> fLayers;

public:
  SimpleStepper(HelixPropagator *prop);
  ~SimpleStepper() {}

  void PropagateToR(double radius, vectorflow::Track &track) const;

  void PropagateInTube(int layer, vectorflow::Track &track) const;

  void PropagateInTube(int layer, std::vector<vectorflow::Track*> const &tracks) const;
};

} // namespace trackml

#endif
