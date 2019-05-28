#ifndef VECTORFLOW_PARTICLE
#define VECTORFLOW_PARTICLE

#include <string>
#include <vector>
#include <map>

namespace geantphysics {
/// Simple class holding few particle properties
class Particle {
public:
  Particle(const std::string &name, int pdgcode, int intcode, double mass, double charge);
  ~Particle() {}

  const std::string &GetName() const { return fName; }
  int GetIndex() const { return fIndex; }
  int GetInternalCode() const { return fInternalCode; }
  int GetPDGCode() const { return fPDGCode; }
  double GetPDGMass() const { return fPDGMass; }
  double GetPDGCharge() const { return fPDGCharge; }

  static int DefineParticles();

  static const Particle *GetParticleByInternalCode(unsigned int intercode)
  {
    if (intercode < gInternalParticleCodes.size()) return gInternalParticleCodes[intercode];
    return nullptr;
  }

  static const Particle *GetParticleByPDGCode(unsigned int pdgcode)
  {
    auto search = gPDGtoInternalCode.find(pdgcode);
    if (search != gPDGtoInternalCode.end()) return gInternalParticleCodes[search->second];
    return nullptr;
  }

  static const Particle *GetParticleByName(const std::string pname)
  {
    auto search = gNametoParticle.find(pname);
    if (search != gNametoParticle.end()) return search->second;
    return nullptr;
  }

  static const std::vector<Particle *> &GetTheParticleTable() { return gTheParticleTable; }

  static const std::vector<Particle *> &GetInternalParticleTable() { return gInternalParticleCodes; }

private:
  std::string fName;
  int fIndex; // in the global particle table
  int fInternalCode;
  int fPDGCode;
  double fPDGMass;
  double fPDGCharge;

  // the particle table
  static std::vector<Particle *> gTheParticleTable;
  static std::vector<Particle *> gInternalParticleCodes;

  // map of PDG codes to internal codes
  static std::map<unsigned int, unsigned int> gPDGtoInternalCode;
  // map of name to particle ptr
  static std::map<const std::string, Particle *> gNametoParticle;
};

} // namespace geantphysics

#endif // PARTICLE_H
