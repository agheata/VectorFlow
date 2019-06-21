#ifndef TRACKML_SIMPLE_STEPPER_H
#define TRACKML_SIMPLE_STEPPER_H

namespace vectorflow {
  class Track;
}

namespace trackml {
/// Very simple stepper that takes a track from vertex and propagates it untill a given distance to the vertex.
/// The propagator makes sure that the error for the last step is lesser than an epsilon value.
class HelixPropagator;

class SimpleStepper {
private:
  HelixPropagator *fPropagator = nullptr;

public:
  SimpleStepper(HelixPropagator *prop) : fPropagator(prop) {}
  ~SimpleStepper() {}

  void PropagateToR(double radius, vectorflow::Track &track) const;
};

} // namespace trackml

#endif
