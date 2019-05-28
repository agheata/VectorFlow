#ifndef GEANT_VECTOR_TYPES_H
#define GEANT_VECTOR_TYPES_H

#include <base/Global.h>
#include <base/Vector3D.h>

namespace vectorflow {

using VectorBackend = vecgeom::VectorBackend;
typedef VectorBackend::Float_v Float_v;
typedef VectorBackend::Double_v Double_v;
typedef VectorBackend::Int_v Int_v;

typedef vecCore::Mask_v<Float_v> MaskF_v;
typedef vecCore::Mask_v<Double_v> MaskD_v;
typedef vecCore::Mask_v<Int_v> MaskI_v;

typedef vecCore::Index<Double_v> IndexD_v;
typedef vecCore::Index<Float_v> IndexF_v;

const int kVecLenF = vecCore::VectorSize<Float_v>();
const int kVecLenD = vecCore::VectorSize<Double_v>();

const auto kVecAlignD = sizeof(Double_v);
const auto kVecAlignF = sizeof(Float_v);

inline
void CopyFltToDbl(Float_v const &flt_v, Double_v &dbl1_v, Double_v &dbl2_v)
{
  // Copy the float SIMD lanes into 2 Double_v variables
  for (size_t lane = 0; lane < kVecLenD; ++lane) {
    vecCore::Set(dbl1_v, lane, (double)vecCore::Get(flt_v, lane));
    vecCore::Set(dbl2_v, lane, (double)vecCore::Get(flt_v, lane + kVecLenD));
  }
}

inline
void CopyFltToDbl(vecgeom::Vector3D<Float_v> const &flt_v, vecgeom::Vector3D<Double_v> &dbl1_v,
                  vecgeom::Vector3D<Double_v> &dbl2_v)
{
  // Copy the float SIMD lanes into 2 Double_v variables
  for (size_t lane = 0; lane < kVecLenD; ++lane) {
    for (size_t i = 0; i < 3; ++i) {
      vecCore::Set(dbl1_v[i], lane, (double)vecCore::Get(flt_v[i], lane));
      vecCore::Set(dbl2_v[i], lane, (double)vecCore::Get(flt_v[i], lane + kVecLenD));
    }
  }
}

inline
void CopyDblToFlt(Double_v const &dbl1_v, Double_v const &dbl2_v, Float_v &flt_v)
{
  // Copy the 2 Double_v SIMD lanes into one Float_v variable
  for (size_t lane = 0; lane < kVecLenD; ++lane) {
    vecCore::Set(flt_v, lane, (float)vecCore::Get(dbl1_v, lane));
    vecCore::Set(flt_v, lane + kVecLenD, (double)vecCore::Get(dbl2_v, lane));
  }
}

inline
void CopyDblToFlt(vecgeom::Vector3D<Double_v> const &dbl1_v, vecgeom::Vector3D<Double_v> const &dbl2_v,
                  vecgeom::Vector3D<Float_v> &flt_v)
{
  // Copy the 2 Double_v SIMD lanes into one Float_v variable
  for (size_t lane = 0; lane < kVecLenD; ++lane) {
    for (size_t i = 0; i < 3; ++i) {
      vecCore::Set(flt_v[i], lane, (float)vecCore::Get(dbl1_v[i], lane));
      vecCore::Set(flt_v[i], lane + kVecLenD, (double)vecCore::Get(dbl2_v[i], lane));
    }
  }
}

} // namespace vectorflow

#endif
