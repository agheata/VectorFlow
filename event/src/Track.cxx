#include "Track.h"

namespace vectorflow {

//______________________________________________________________________________
Track::Track(size_t maxdepth)
{
  if (maxdepth) {
    fPath = VolumePath_t::MakeInstance(maxdepth);
    fNextpath = VolumePath_t::MakeInstance(maxdepth);
  }
}

//______________________________________________________________________________
Track::~Track()
{
  VolumePath_t::ReleaseInstance(fPath);
  VolumePath_t::ReleaseInstance(fNextpath);
}

//______________________________________________________________________________
Track &Track::operator=(const Track &other)
{
  // Assignment
  if (&other != this) {
    fEvent              = other.fEvent;
    fParticle           = other.fParticle;
    fPrimaryIndx        = other.fPrimaryIndx;
    fMother             = other.fMother;
    fCharge             = other.fCharge;
    fNsteps             = other.fNsteps;
    fMaxDepth           = other.fMaxDepth;
    fGeneration         = other.fGeneration;
    fSpecies            = other.fSpecies;
    fStatus             = other.fStatus;
    fMass               = other.fMass;
    fPos                = other.fPos;
    fDir                = other.fDir;
    fP                  = other.fP;
    fE                  = other.fE;
    fTime               = other.fTime;
    fEdep               = other.fEdep;
    fPstep              = other.fPstep;
    fStep               = other.fStep;
    fSnext              = other.fSnext;
    fSafety             = other.fSafety;
    fBoundary           = other.fBoundary;
    fVolume             = other.fVolume;

    *fPath     = *other.fPath;
    *fNextpath = *other.fNextpath;
  }
  return *this;
}

//______________________________________________________________________________
void Track::Clear(const char *)
{
  // Resets track content.
  fEvent       = -1;
  fParticle    = -1;
  fPrimaryIndx = -1;
  fMother      = 0;
  fCharge      = 0;
  fNsteps      = 0;
  fSpecies     = kHadron;
  fStatus      = kAlive;
  fMass        = 0.;
  fPos         = 0.;
  fDir         = 0.;
  fP           = 0.;
  fE           = 0.;
  fTime        = 0.;
  fEdep        = 0;
  fPstep       = 1.E20;
  fStep        = 0.;
  fSnext       = 0.;
  fSafety      = 0.;
  fBoundary    = false;
  fMaxDepth    = 0;
  fGeneration  = 0;
  fVolume      = nullptr;
  fPath->Clear();
  fNextpath->Clear();
}

//______________________________________________________________________________
void Track::Reset(Track const &blueprint)
{
  // Fast reset of the content of an existing track, avoiding in-place construction

  // Copy main content from blueprint, except the pointers to geometry states and user data
  memcpy(&fEvent, &blueprint.fEvent, sizeof(Track) - 3 * sizeof(void *));

  // Clear Geometry path
  fPath->Clear();
  fNextpath->Clear();
}

//______________________________________________________________________________
void Track::SetMaxDepth(int depth)
{
  if (fPath) {
    VolumePath_t::ReleaseInstance(fPath);
    VolumePath_t::ReleaseInstance(fNextpath);
  }
  fPath = VolumePath_t::MakeInstance(depth);
  fNextpath = VolumePath_t::MakeInstance(depth);
}

//______________________________________________________________________________
void Track::SetPath(VolumePath_t const *const path)
{
  // Set path.
  *fPath = *path;
  UpdateVolume();
}

//______________________________________________________________________________
void Track::SetNextPath(VolumePath_t const *const path)
{
  // Set next path.
  *fNextpath = *path;
}

//______________________________________________________________________________
void Track::Print(const char *msg) const
{
  const char *status[8] = {"alive", "killed", "inflight", "boundary", "exitSetup", "physics", "postponed", "new"};

  printf("%s: evt=%d part=%d prim=%d mth=%d chg=%d nstp=%d spc=%d status=%s mass=%g "
         "xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g time=%g pstp=%g stp=%g snxt=%g saf=%g bdr=%d\n",
         msg, fEvent, fParticle, fPrimaryIndx, fMother, fCharge, fNsteps,
         (int)fSpecies, status[int(fStatus)], fMass, fPos[0], fPos[1], fPos[2], fDir[0], fDir[1], fDir[2], fP, fE, fTime, fPstep,
         fStep, fSnext, fSafety, fBoundary);
}

//______________________________________________________________________________
void Track::PrintTracks(TrackVec_t &tracks)
{
  for (auto track : tracks)
    track->Print("xxx");
}

} // namespace vectorflow
