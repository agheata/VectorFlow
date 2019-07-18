
#include "Particle.h"

#include "vectorFlow/PhysicalConstants.h"
#include "vectorFlow/SystemOfUnits.h"

namespace geantphysics {

std::vector<Particle *> Particle::gTheParticleTable;
std::vector<Particle *> Particle::gInternalParticleCodes;
std::map<unsigned int, unsigned int> Particle::gPDGtoInternalCode;
std::map<const std::string, Particle *> Particle::gNametoParticle;

Particle::Particle(const std::string &name, int pdgcode, int intcode, double mass, double charge)
    : fName(name), fIndex(-1), fInternalCode(intcode), fPDGCode(pdgcode), fPDGMass(mass), fPDGCharge(charge)
{
  fIndex = gTheParticleTable.size();
  gTheParticleTable.push_back(this);
  //
  unsigned long icode = intcode;
  if (gInternalParticleCodes.size() < icode + 1) {
    gInternalParticleCodes.resize(icode + 1, nullptr);
  }
  gInternalParticleCodes[icode] = this;
  //
  gPDGtoInternalCode[pdgcode] = icode;
  gNametoParticle[name]       = this;
}

int Particle::DefineParticles()
{
  new Particle("gamma",     22,  0, 0.0, 0.0);
  new Particle("e-",        11,  1, geant::units::kElectronMassC2, -1.0 * geant::units::eplus);
  new Particle("e+",       -11,  2, geant::units::kElectronMassC2,  1.0 * geant::units::eplus);
  new Particle("pi-",     -211,  3, 0.13957 * geant::units::GeV,   -1.0 * geant::units::eplus);
  new Particle("pi+",      211,  4, 0.13957 * geant::units::GeV,    1.0 * geant::units::eplus);
  new Particle("pi0",      111,  5, 0.134977 * geant::units::GeV,   0);
  new Particle("proton",  2212,  6, geant::units::kProtonMassC2,    1.0 * geant::units::eplus);
  new Particle("neutron", 2112,  7, geant::units::kNeutronMassC2,   0);
  new Particle("K_L0",     130,  8, 0.497614 * geant::units::GeV,   0);
  new Particle("K-",      -321,  9, 0.493677 * geant::units::GeV,  -1.0 * geant::units::eplus);
  new Particle("K+",       321, 10, 0.493677 * geant::units::GeV,   1.0 * geant::units::eplus);
  new Particle("K_S0",     310, 11, 0.497614 * geant::units::GeV,   0);
  new Particle("K0",       311, 12, 0.497614 * geant::units::GeV,   0);
  // printf("=== DefineParticles: %zu particles registered to the particle table\n", Particle::GetTheParticleTable().size());
}

} // namespace geantphysics
