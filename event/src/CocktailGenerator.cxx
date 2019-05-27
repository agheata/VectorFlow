#include "CocktailGenerator.h"
#include <cassert>
#include "Event.h"

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
  return true;
}

CocktailGenerator::Event_t *CocktailGenerator::NextEvent()
{
  using namespace vecCore::math;
  Event_t *event = new Event_t();
  // Sample the number of tracks
  int mintracks = Max(0, 2*fAveragePrimaries - fMaxPrimaries);
  int maxtracks = fMaxPrimaries;
  int ntracks = mintracks - (maxtracks - mintracks) * fRng.Uniform<Backend_t>();
  return event;
}
