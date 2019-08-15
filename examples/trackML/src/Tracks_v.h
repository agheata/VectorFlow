#ifndef VECTORFLOW_TRACKS_V_H
#define VECTORFLOW_TRACKS_V_H

#include <VecCore/VecCore>
#include "base/Vector3D.h"
#include "vectorFlow/VectorTypes.h"

using namespace vectorflow;
using vecCore::Set;
using vecCore::Get;

template <typename Data>
struct Tracks_v {  
  // Data in vecCore vector types
  vecgeom::Vector3D<Double_v> fPos_v;
  vecgeom::Vector3D<Double_v> fDir_v;
  Double_v fCharge_v, fMomentum_v, fNSteps_v;
  
  // Gather info to fill one SIMD lane
  void Gather(Data *track, const std::size_t lane) {
    Set(fPos_v.x(),  lane, track->PosX());
    Set(fPos_v.y(),  lane, track->PosY());
    Set(fPos_v.z(),  lane, track->PosZ());
    Set(fDir_v.x(),  lane, track->DirX());
    Set(fDir_v.y(),  lane, track->DirY());
    Set(fDir_v.z(),  lane, track->DirZ());
    Set(fCharge_v,   lane, track->Charge());
    Set(fMomentum_v, lane, track->P());
    Set(fNSteps_v,   lane, track->GetNsteps());
  }

  // Scatter one lane info back to original structure 
  void Scatter(const std::size_t lane, Data *track) {
    track->SetPosition(Get(fPos_v.x(), lane),  Get(fPos_v.y(), lane), Get(fPos_v.z(), lane));
    track->SetDirection(Get(fDir_v.x(), lane), Get(fDir_v.y(), lane), Get(fDir_v.z(), lane));
    track->SetNsteps(Get(fNSteps_v, lane));
  }

  // Auxiliar vectorized methods
  Double_v Pt_v() const { return fMomentum_v * fDir_v.Perp(); }

  Double_v Curvature_v(double Bz) {
    using constants::kB2C;
    using constants::kTiny;
    using namespace vecCore::math;

    Double_v qB_v = fCharge_v * Bz;
    Double_v curvature_v = Abs(kB2C * qB_v / (Pt_v() + kTiny));
    
    return vecCore::Blend<Double_v>(Abs(qB_v) < kTiny, kTiny, curvature_v);
  }
};

#endif
