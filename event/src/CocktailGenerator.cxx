#include "CocktailGenerator.h"
#include <cassert>
#include "Event.h"
#include "Track.h"
#include "Particle.h"
#include "vectorFlow/PhysicalConstants.h"
#include <map>
#include <string>

CocktailGenerator::CocktailGenerator(vecgeom::Vector3D<double> const &vertex)
  : vectorflow::PrimaryGenerator<CocktailGenerator>(vertex)
{
  geantphysics::Particle::DefineParticles();
  fCurrentEvent.store(0);
  fRndm = new vectorflow::RngWrapper;
}

CocktailGenerator::~CocktailGenerator()
{
  delete fRndm;
}

bool CocktailGenerator::InitPrimaryGenerator()
{
  // Check validity.
  if (fInitialized) return true;
  using namespace vecCore::math;
  double sumw = 0;
  for (auto part : fSpecies) sumw += part.second;
  if (Abs(sumw - 1.) > 1.E-6) {
    printf("CocktailGenerator: species weights not normalized\n");
    return false;
  }
  if (fMaxBeamEnergy <= fMinBeamEnergy) {
    printf("CocktailGenerator: min/max beam energy not set or not consistent\n");
    return false;
  }
  if (fMaxPrimaries <= 0 || fAveragePrimaries <= 0 || fAveragePrimaries > fMaxPrimaries) {
    printf("CocktailGenerator: max/average number of primary tracks not set or not consistent\n");
    return false;
  }
  if (fMaxDepth == 0) {
    printf("CocktailGenerator: maximum geometry depth not set\n");
    return false;
  }
  return true;
}

void CocktailGenerator::AddPrimary(const char *name, float weight)
{
  const geantphysics::Particle *particle = geantphysics::Particle::GetParticleByName(std::string(name));
  if (!particle) {
    printf("AddPrimary: Error: particle %s was not defined\n", name);
    return;
  }
  auto p = std::make_pair(particle->GetInternalCode(), weight);
  if (std::find(fSpecies.begin(), fSpecies.end(), p) == fSpecies.end())
    fSpecies.push_back(p);
}

CocktailGenerator::Event_t *CocktailGenerator::NextEvent()
{
  using namespace vecCore::math;
  using Vector3 = vecgeom::Vector3D<double>;
  const int kGamma = 0;
  const int kElectron = 1;
  const int kPositron = 2;
  Event_t *event = new Event_t();
  // Sample the number of tracks
  int mintracks = Max(0, 2*fAveragePrimaries - fMaxPrimaries);
  int maxtracks = fMaxPrimaries;
  int ntracks = fRndm->uniform(mintracks, maxtracks);

  // Generate random tracks having 
  for (int i = 0; i < ntracks; ++i) {
    // Generate particle species
    double prob = fRndm->uniform();
    double probi = 0.;
    unsigned int pid = 0;
    for (auto part : fSpecies) {
      probi += part.second;
      if (prob < probi) {
        pid = part.first;
        break;
      }
    }
    const geantphysics::Particle *particle = geantphysics::Particle::GetParticleByInternalCode(pid);
    assert(particle && "Fatal: wrong particle pid");
    
    // Generate particle direction
    Vector3 dir;
    double phi = 2 * M_PI * fRndm->uniform();
    double theta = ACos(1. - 2. * fRndm->uniform());
    dir.Set(Sin(theta) * Cos(phi), Sin(theta) * Sin(phi), Cos(theta));
    dir.Normalize();

    // Generate particle energy
    double energy = fRndm->uniform(fMinBeamEnergy, fMaxBeamEnergy);
    vectorflow::Track *track = new vectorflow::Track(fMaxDepth);
    track->SetEvent(fCurrentEvent++);
    track->SetParticle(pid);
    track->SetPrimaryParticleIndex(i);
    track->SetCharge((int)particle->GetPDGCharge());
    track->SetMass(particle->GetPDGMass());
    track->SetPosition(fVertex);
    track->SetDirection(dir);
    track->SetEkin(energy);
    double e = track->E();
    double m = track->Mass();
    track->SetP(std::sqrt((e - m) * (e + m)));
    // geometry state not yet initialized
    // Add track to the current event
    event->AddPrimary(track);
  }

  return event;
}
