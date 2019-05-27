#ifndef VECTORFLOW_TYPEDEFS_H
#define VECTORFLOW_TYPEDEFS_H

//#include "Geant/VectorTypes.h"

#include <vector>
template <class T>
using vector_t = std::vector<T>;

namespace geantphysics {
class Particle;
} // namespace geantphysics

typedef geantphysics::Particle Particle_t;

#include "navigation/NavigationState.h"
typedef VECGEOM_NAMESPACE::NavigationState VolumePath_t;
#include "volumes/LogicalVolume.h"
typedef VECGEOM_NAMESPACE::LogicalVolume Volume_t;
#include "volumes/PlacedVolume.h"
typedef VECGEOM_NAMESPACE::VPlacedVolume Node_t;
#endif
