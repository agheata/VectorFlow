//===--- Track.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Track.h
 * @brief Implementation of track for GeantV prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef VECTORFLOW_TRACK
#define VECTORFLOW_TRACK

#include "base/Vector3D.h"
#include "vectorFlow/math_wrappers.h"

#include <functional>
#include <climits>
#include <float.h>
#include <atomic>
#include <mutex>
#include "vectorFlow/Typedefs.h"

#ifndef GEANT_ALIGN_PADDING
#define GEANT_ALIGN_PADDING 64
#endif

namespace vectorflow {

namespace constants {
  constexpr double kB2C  = -0.299792458e-3;
  constexpr double kTiny = 1.E-50;
}

/// The definition of track status flags.
enum TrackStatus_t { kAlive, kKilled, kInFlight, kBoundary, kExitingSetup, kPhysics, kPostponed, kNew };

/// The definition of different transport actions
enum TransportAction_t {
  kDone     = 0, ///< Return immediately - no tracks left
  kPostpone = 1, ///< Return imediately and postpone whatever tracks left
  kSingle   = 2, ///< Perform remaining loop in single track mode
  kVector   = 3  ///< Perform remaining loop in vectorized mode
};

class Track;

using InplaceConstructFunc_t = std::function<void(void *)>;
using PrintUserDataFunc_t    = std::function<void(void *)>;

typedef std::vector<Track *> TrackVec_t;

/// Track class.
class Track {

  template <typename T>
  using Vector3D = vecgeom::Vector3D<T>;

private:
  int fEvent               = -1;      ///< Event number
  int fParticle            = -1;      ///< Particle type
  int fIndex               = -1;      ///< Index of corresponding particle
  int fPrimaryIndx         = -1;      ///< Index of the primary particle in the current event
  int fMother              = -1;      ///< Index of mother particle
  int fCharge              = 0;       ///< Particle charge
  int fNsteps              = 0;       ///< Number of steps made
  int fMaxDepth            = 0;       ///< Maximum geometry depth
  int fGeneration          = 0;       ///< Track generation: 0=primary
  TrackStatus_t fStatus    = kAlive;  ///< Track status
  double fMass             = 0;       ///< Particle mass
  Vector3D<double> fPos;              ///< Particle position
  Vector3D<double> fDir;              ///< Particle direction
  double fP                = 0;       ///< Momentum
  double fE                = 0;       ///< Energy
  double fTime             = 0;       ///< Time
  double fEdep             = 0;       ///< Energy deposition in the step
  double fPstep            = 1.e+20;  ///< Selected physical step
  double fStep             = 0;       ///< Current step
  double fSnext            = 0;       ///< Straight distance to next boundary
  double fSafety           = 0;       ///< Safe distance to any boundary
  bool fBoundary           = false;   ///< True if starting from boundary
  Volume_t const *fVolume  = nullptr; ///< Current volume the particle is in


private:
  VolumePath_t *fPath     = nullptr; ///< Paths for the particle in the geometry */
  VolumePath_t *fNextpath = nullptr; ///< Path for next volume */

public:

  /// Track constructor
  Track(size_t maxdepth = 0);

public:
  /// Track destructor
  ~Track();

  /// Track copy constructor
  Track(const Track &other) = default; // = delete;

  /// @brief Assignment
  Track &operator=(const Track &other);

  /// Getter for the event number
  VECCORE_FORCE_INLINE
  int Event() const { return fEvent; }

  /// Getter for the index of corresponding particle
  VECCORE_FORCE_INLINE
  int Particle() const { return fParticle; }

  /// Getter for the index of the primary particle in the current event
  VECCORE_FORCE_INLINE
  int PrimaryParticleIndex() const { return fPrimaryIndx; }

  /// Getter for the index of the track
  VECCORE_FORCE_INLINE
  int Index() const { return fIndex; }

  /// Getter for thes index of mother particle
  VECCORE_FORCE_INLINE
  int Mother() const { return fMother; }

  /// Getter for the charge value
  VECCORE_FORCE_INLINE
  int Charge() const { return fCharge; }

  /** @brief Fills scalar charge from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param e_v SIMD charges to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_FORCE_INLINE
  static size_t GetCharge_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks,
                            Real_v &charge_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem              = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem);
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(charge_v, i, tracks[offset + i]->Charge());
    return nelem;
  }

  /// Getter for the number of physical step made
  VECCORE_FORCE_INLINE
  int GetNsteps() const { return fNsteps; }

  /// Maximum geometry depth
  VECCORE_FORCE_INLINE
  int GetMaxDepth() const { return fMaxDepth; }

  /// Getter for track generation (0 = primary)
  VECCORE_FORCE_INLINE
  int GetGeneration() const { return fGeneration; }

  /// Getter for the track status
  TrackStatus_t Status() const { return fStatus; }

  /// Getter for the rest mass value
  VECCORE_FORCE_INLINE
  double Mass() const { return fMass; }

  /// Position
  VECCORE_FORCE_INLINE
  Vector3D<double> const &Position() const { return fPos; }

  /// Position x component
  VECCORE_FORCE_INLINE
  double PosX() const { return fPos[0]; }
  
  /// Position y component
  VECCORE_FORCE_INLINE
  double PosY() const { return fPos[1]; }

  /// Position z component
  VECCORE_FORCE_INLINE
  double PosZ() const { return fPos[2]; }
  
  /** Fills scalar position components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param pos_v SIMD position to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_FORCE_INLINE
  static size_t GetPos_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks,
                         Vector3D<Real_v> &pos_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem              = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem);
    for (size_t i = 0; i < nelem; ++i) {
      vecCore::Set(pos_v[0], i, tracks[offset + i]->Position()[0]);
      vecCore::Set(pos_v[1], i, tracks[offset + i]->Position()[1]);
      vecCore::Set(pos_v[2], i, tracks[offset + i]->Position()[2]);
    }
    return nelem;
  }

  /// Direction
  VECCORE_FORCE_INLINE
  Vector3D<double> const &Direction() const { return fDir; }
 
  /// Direction x component
  VECCORE_FORCE_INLINE
  double DirX() const { return fDir[0]; }
  
  /// Direction y component
  VECCORE_FORCE_INLINE
  double DirY() const { return fDir[1]; }

  /// Direction z component
  VECCORE_FORCE_INLINE
  double DirZ() const { return fDir[2]; }
  
  /** Fills scalar direction components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param dir_v SIMD direction to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_FORCE_INLINE
  static size_t GetDir_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks,
                         Vector3D<Real_v> &dir_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem              = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem);
    for (size_t i = 0; i < nelem; ++i) {
      vecCore::Set(dir_v[0], i, tracks[offset + i]->Direction()[0]);
      vecCore::Set(dir_v[1], i, tracks[offset + i]->Direction()[1]);
      vecCore::Set(dir_v[2], i, tracks[offset + i]->Direction()[2]);
    }
    return nelem;
  }

  /// Getter for the momentum value */
  
  
  VECCORE_FORCE_INLINE
  double P() const { return fP; }

  /// Getter for the momentum X component */
  
  
  VECCORE_FORCE_INLINE
  double Px() const { return fP * fDir[0]; }

  /// Getter for the momentum Y component */
  
  
  VECCORE_FORCE_INLINE
  double Py() const { return fP * fDir[1]; }

  /// Getter for the momentum Z component */
  
  
  VECCORE_FORCE_INLINE
  double Pz() const { return fP * fDir[2]; }

  /** Fills scalar momentum components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param pos_v SIMD momentum to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   */
  template <typename Real_v, bool Tail = false>
  VECCORE_FORCE_INLINE
  static size_t GetP_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks,
                       Vector3D<Real_v> &mom_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem              = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem);
    for (size_t i = 0; i < nelem; ++i) {
      vecCore::Set(mom_v[0], i, tracks[offset + i]->Px());
      vecCore::Set(mom_v[1], i, tracks[offset + i]->Py());
      vecCore::Set(mom_v[2], i, tracks[offset + i]->Pz());
    }
    return nelem;
  }

  /// Getter for the module momentum's value
  VECCORE_FORCE_INLINE
  double Pt() const { return fP * fDir.Perp(); }

  /// Getter for the energy value
  VECCORE_FORCE_INLINE
  double E() const { return fE; }

  /** Fills scalar energy components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param e_v SIMD energies to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_FORCE_INLINE
  static size_t GetE_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks,
                       Real_v &e_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem              = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem);
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(e_v, i, tracks[offset + i]->E());
    return nelem;
  }

  /// Getter for the kinetic energy value
  VECCORE_FORCE_INLINE
  double Ekin() const { return (fE - fMass); }

  /** Fills scalar kinetic energy components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param t_v SIMD kinetic energies to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_FORCE_INLINE
  static size_t GetEkin_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks,
                          Real_v &t_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem              = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem);
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(t_v, i, tracks[offset + i]->Ekin());
    return nelem;
  }

  ///< Getter for the time
  VECCORE_FORCE_INLINE
  double Time() const { return fTime; }

  /// Getter for the energy deposition value
  VECCORE_FORCE_INLINE
  double Edep() const { return fEdep; }

  /// Getter for the selected physical step
  VECCORE_FORCE_INLINE
  double GetPstep() const { return fPstep; }

  /// Getter for the physical step
  VECCORE_FORCE_INLINE
  double GetStep() const { return fStep; }

  /// Getter for the time traveled in the step
  VECCORE_FORCE_INLINE
  double TimeStep(double step) const { return fE * step / fP; }

  /// Getter for the straight distance to next boundary
  VECCORE_FORCE_INLINE
  double GetSnext() const { return fSnext; }

  /// Getter for the safe distance to any boundary
  VECCORE_FORCE_INLINE
  double GetSafety() const { return fSafety; }

  /// Getter for the true if starting from boundary
  VECCORE_FORCE_INLINE
  bool Boundary() const { return fBoundary; }

  /// Getter for the volume
  VECCORE_FORCE_INLINE
  Volume_t const *GetVolume() const { return fVolume; }

  /// Getter for thes current path
  VECCORE_FORCE_INLINE
  VolumePath_t *Path() const { return fPath; }

  /// Getter for thes next path
  VECCORE_FORCE_INLINE
  VolumePath_t *NextPath() const { return fNextpath; }

  /// Getter for the next volume (no backup)
  VECCORE_FORCE_INLINE
  Volume_t const *GetNextVolume() const { return (fNextpath->Top()->GetLogicalVolume()); }

  /// Getter for the beta value
  VECCORE_FORCE_INLINE
  double Beta() const { return fP / fE; }

  /** Fills scalar beta from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param beta_v SIMD betas to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_FORCE_INLINE
  static size_t GetBeta_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks,
                          Real_v &beta_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem              = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem);
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(beta_v, i, tracks[offset + i]->Beta());
    return nelem;
  }

  /// Getter for the curvature. To be changed when handling properly field
  VECCORE_FORCE_INLINE
  double Curvature(double Bz) const
  {
    using constants::kB2C;
    using constants::kTiny;
    // Curvature
    double qB = fCharge * Bz;
    if (fabs(qB) < kTiny) return kTiny;
    return fabs(kB2C * qB / (Pt() + kTiny));
  }

  /// Getter for the gamma value
  VECCORE_FORCE_INLINE
  double Gamma() const { return fMass ? fE / fMass : DBL_MAX; }

  /** Fills scalar gamma components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param e_v SIMD gammas to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_FORCE_INLINE
  static size_t GetGamma_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks,
                           Real_v &gamma_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem              = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem);
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(gamma_v, i, tracks[offset + i]->Gamma());
    return nelem;
  }

  /// Function that check if track is alive
  VECCORE_FORCE_INLINE
  bool IsAlive() const { return (fStatus != kKilled); }

  /// Function that check if track is on boundary
  VECCORE_FORCE_INLINE
  bool IsOnBoundary() const { return (fStatus == kBoundary); }

  ///  Check direction normalization within tolerance
  VECCORE_FORCE_INLINE
  bool IsNormalized() const { return fDir.IsNormalized(); }

  /// Setter for event number
  VECCORE_FORCE_INLINE
  void SetEvent(int event) { fEvent = event; }

  /// Setter for particle id
  VECCORE_FORCE_INLINE
  void SetParticle(int particle) { fParticle = particle; }

  /// Setter for primary particle index in the current event
  VECCORE_FORCE_INLINE
  void SetPrimaryParticleIndex(int primaryindx) { fPrimaryIndx = primaryindx; }

  /// Setter for the track index in the current event
  VECCORE_FORCE_INLINE
  void SetIndex(int indx) { fIndex = indx; }

  /// Setter for mother index
  VECCORE_FORCE_INLINE
  void SetMother(int mother) { fMother = mother; }

  /// Setter for charge
  VECCORE_FORCE_INLINE
  void SetCharge(int charge) { fCharge = charge; }

  /// Setter for number of steps
  VECCORE_FORCE_INLINE
  void SetNsteps(int nsteps) { fNsteps = nsteps; }

  /// Increment number of steps
  VECCORE_FORCE_INLINE
  void IncrementNsteps(int nsteps = 1) { fNsteps += nsteps; }

  /// Setter for particle generation
  VECCORE_FORCE_INLINE
  void SetGeneration(int generation) { fGeneration = generation; }

  /// Setter for track status
  VECCORE_FORCE_INLINE
  void SetStatus(TrackStatus_t status) { fStatus = status; }

  /// Setter for particle mass
  VECCORE_FORCE_INLINE
  void SetMass(double mass) { fMass = mass; }

  /// Setter for the maximum geometry depth
  VECCORE_FORCE_INLINE
  void SetMaxDepth(int depth);

  /// Setter for position from components
  VECCORE_FORCE_INLINE
  void SetPosition(double x, double y, double z)  { fPos.Set(x, y, z); }

  /// Setter for position from vector
  VECCORE_FORCE_INLINE
  void SetPosition(Vector3D<double> const &pos) { fPos = pos; }

  /// Setter for direction from components
  VECCORE_FORCE_INLINE
  void SetDirection(double dx, double dy, double dz) { fDir.Set(dx, dy, dz); }

  /// Setter for direction from components
  VECCORE_FORCE_INLINE
  void SetDirection(Vector3D<double> const &dir) { fDir = dir; }

  /// Setter for momentum
  VECCORE_FORCE_INLINE
  void SetP(double p) { fP = p; }

  /// Setter for energy
  VECCORE_FORCE_INLINE
  void SetE(double e) { fE = e; }

  /// Setter for the kinetic energy
  VECCORE_FORCE_INLINE
  void SetEkin(double ekin) { fE = ekin + fMass; }

  /// Setter for energy
  VECCORE_FORCE_INLINE
  void DecreaseE(double e) { fE -= e; }

  /// Setter for time
  VECCORE_FORCE_INLINE
  void SetTime(double time) { fTime = time; }

  /// Increase time
  VECCORE_FORCE_INLINE
  void IncreaseTime(double time) { fTime += time; }

  /// Setter for energy deposition
  VECCORE_FORCE_INLINE
  void SetEdep(double edep) { fEdep = edep; }

  /// Setter for energy deposition
  VECCORE_FORCE_INLINE
  void IncreaseEdep(double edep) { fEdep += edep; }

  /// Setter for proposed physics step
  VECCORE_FORCE_INLINE
  void SetPstep(double pstep) { fPstep = pstep; }

  /// Decrease proposed physics step after some propagation
  VECCORE_FORCE_INLINE
  void DecreasePstep(double len) { fPstep -= len; }

  /// Setter for the physics step
  VECCORE_FORCE_INLINE
  void SetStep(double step) { fStep = step; }

  /// Increase progression step after some propagation
  VECCORE_FORCE_INLINE
  void IncreaseStep(double len) { fStep += len; }

  /// Setter for straight distance to next boundary
  VECCORE_FORCE_INLINE
  void SetSnext(double snext) { fSnext = snext; }

  /// Decrease snext after some propagation
  VECCORE_FORCE_INLINE
  void DecreaseSnext(double len) { fSnext -= len; }

  /// Setter for the isotropic safe distance to volumes
  VECCORE_FORCE_INLINE
  void SetSafety(double safety) { fSafety = safety; }

  /// Decrease safety after some propagation
  VECCORE_FORCE_INLINE
  void DecreaseSafety(double len) { fSafety -= len; }

  /// Setter for the starting from boundary flag
  VECCORE_FORCE_INLINE
  void SetBoundary(bool flag) { fBoundary = flag; }

  /// Setter for the current geometry path
  VECCORE_FORCE_INLINE
  void SetPath(VolumePath_t const *const path);

  /// Setter for the next geometry path
  VECCORE_FORCE_INLINE
  void SetNextPath(VolumePath_t const *const path);

  /// Function that stops the track depositing its kinetic energy
  VECCORE_FORCE_INLINE
  void Stop()
  {
    fEdep   = Ekin();
    fE      = fMass;
    fP      = 0;
    fStatus = kKilled;
  }

  /// Setter for the status killed to track
  VECCORE_FORCE_INLINE
  void Kill() { fStatus = kKilled; }

  ///< Clear function
  VECCORE_FORCE_INLINE
  void Clear(const char *option = "");

  ///< Fast reset function
  void Reset(Track const &blueprint);

  /// Print function
  void Print(const char *msg = "") const;

  /// Print function for a container of tracks
  static void PrintTracks(TrackVec_t &tracks);

  /// Function that swaps path and next path
  VECCORE_FORCE_INLINE
  void UpdateSwapPath()
  {
    VolumePath_t *tmp = fNextpath;
    fNextpath         = fPath;
    fPath             = tmp;
    UpdateVolume();
  }

  /// Function that makes next path have the same content as the current.
  VECCORE_FORCE_INLINE
  void UpdateSameNextPath() { *fNextpath = *fPath; }

  /// Function that updates the current volume the particle is in
  VECCORE_FORCE_INLINE
  void UpdateVolume() { fVolume = fPath->Top()->GetLogicalVolume(); }

  /// Function to normalize direction
  VECCORE_FORCE_INLINE
  void Normalize() { fDir.Normalize(); }
  
  /// Function to fast normalize direction
  VECCORE_FORCE_INLINE
  void NormalizeFast()
  {
    // Liniarized normalization function, working in the assumption that the direction is quasi-normalized
    fDir *= 1.5 - 0.5 * fDir.Mag2();
  }

  /// Function to make a step along the current direction
  VECCORE_FORCE_INLINE
  void MakeStep(double step)
  {
    fPstep -= step;
    fStep += step;
    fSafety -= step;
    fSnext -= step;
    fPos += step * fDir;
  }
}; // Track
} // namespace vectorflow

#endif
