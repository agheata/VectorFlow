//===--- PrimaryGenerator.h - Geant-V ---------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file PrimaryGenerator.h
 * @brief Implementation of primary generators for Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef VECTORFLOW_PrimaryGenerator_h
#define VECTORFLOW_PrimaryGenerator_h

#include <VecGeom/base/Vector3D.h>

template <typename DerivedT> struct Generator_traits;

namespace vectorflow {

using vecgeom::kRadToDeg;
using vecgeom::kDegToRad;

/**
 * @brief Class of primary generators
 */
template <typename DerivedT>
class PrimaryGenerator {
protected:
  using Event_t = typename Generator_traits<DerivedT>::Event_t;
  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

  Vector3D<double> fVertex;                      ///< Primary vertex position
public:
  PrimaryGenerator(Vector3D<double> const &vertex) : fVertex(vertex) {}
  
  ~PrimaryGenerator() {}
  
  static double EtaToTheta(double eta) { return (2. * atan(exp(-eta)) * kRadToDeg); }
  static double ThetaToEta(double theta) { return (-log(tan(0.5 * theta * kDegToRad))); }

  /// Interface for initialization of primary generator
  void InitPrimaryGenerator() { static_cast<DerivedT *>(this)->template InitPrimaryGenerator(); }

  ///  Method to produce roduce the next event
  Event_t *NextEvent() { return (static_cast<DerivedT *>(this)->template NextEvent()); }
}; // PrimaryGenerator

} // namespace vectorflow

#endif
