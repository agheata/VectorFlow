
#ifndef VECTORFLOW_COCKTAIL_GENERATOR_H
#define VECTORFLOW_COCKTAIL_GENERATOR_H

#include "PrimaryGenerator.h"
#include "Typedefs.h"
#include "VecMath/Rng.h"

namespace vectorflow {
  class Event;
}

// geantphysics classdef
namespace geantphysics {
  class Particle;
}

#include <string>
#include <utility>

class CocktailGenerator;

template <> struct Generator_traits<CocktailGenerator> {
  using Event_t = vectorflow::Event;
};

class CocktailGenerator : public vectorflow::PrimaryGenerator<CocktailGenerator> {
public:
  using Base_t  = typename vectorflow::PrimaryGenerator<CocktailGenerator>;
  using Event_t = vectorflow::Event;
  using Backend_t = vecCore::backend::Scalar;
  // CTR DTR
  CocktailGenerator(vecgeom::Vector3D<double> const &vertex) : Base_t(vertex) {}
  ~CocktailGenerator() {}

  // interface methods
  bool InitPrimaryGenerator();

  Event_t *NextEvent();

  // public setters/getters
  void SetMaxPrimaryPerEvt(const int pperevt) { fMaxPrimaries = pperevt; }
  int GetMaxPrimaryPerEvt() const { return fMaxPrimaries; }
  void SetAvgPrimaryPerEvt(const int average) { fAveragePrimaries = average; }
  int GetAvgPrimaryPerEvt() const { return fAveragePrimaries; }
  void SetPrimaryEnergyRange(const double emin, const double emax) { fMinBeamEnergy = emin; fMaxBeamEnergy = emax; }
  void AddPrimary(const char *name, float weight);
  /// Get primary particle species by index, retrieving its weight in the generator
  int GetPrimary(int index, double &weight);
  void SetVertex(vecgeom::Vector3D<double> const &vertex) { fVertex = vertex; }

  void Print();
  int GetNumberOfPrimaryTypes() { return fSpecies.size(); }

private:
  // CocktailGenerator() = delete;
  CocktailGenerator(const CocktailGenerator &) = delete;
  CocktailGenerator &operator=(const CocktailGenerator &) = delete;

private:
  bool fInitialized;                               ///< Generator initialized
  std::vector<std::pair<int, double>> fSpecies;    ///< Species of particles and their weights in the generator
  int fMaxPrimaries            = 0;                ///< Maximum number of primaries per event
  int fAveragePrimaries        = 0;                ///< average number of primaries per event
  
//  std::string gNameParticlesVector[];
//  std::map<std::string, int> gPrimaryNameToIndexMap;
  
  double fMinBeamEnergy        = 0.;               ///< Minimum particle energy
  double fMaxBeamEnergy        = 0.;               ///< Maximum particle energy

  vecgeom::Vector3D<double> fVertex;               ///< Vertex position
  vecRng::MRG32k3a<Backend_t> fRng; ///< Random number generator
};

#endif // CocktailGenerator
