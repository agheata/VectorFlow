#include "Event.h"
#include "Track.h"
#include <iostream>
#include <string>

namespace vectorflow {

//______________________________________________________________________________
int Event::AddTrack()
{
  // Thread safe track addition
  int ntracks   = ++fNtracks;
  return (ntracks - 1);
}

//______________________________________________________________________________
void Event::Clear()
{
  // Clear the event.
  fEvent       = 0;
  fNtracks.store(0);
  // Release primary tracks
  for (auto track : fPrimaries) delete track;
  fPrimaries.clear();
}

//______________________________________________________________________________
void Event::Print(const char * str) const
{
  // Print events content
  std::cout << "Event " << fEvent << ": " << GetNprimaries() << " primaries, " << GetNtracks() << " tracks." << std::endl;
  if (strcmp(str, "ALL") == 0) vectorflow::Track::PrintTracks(fPrimaries);
}

} // namespace geant
