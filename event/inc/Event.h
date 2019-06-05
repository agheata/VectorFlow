//===--- Event.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file Event.h
 * @brief Implementation of event for GeantV prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef VECTORFLOW_EVENT
#define VECTORFLOW_EVENT

#include <atomic>
#include <vector>
#include "base/Vector3D.h"

namespace vectorflow {

class Track;

/// Class Event that decribes interaction events with tracks emerging from the same vertex
class Event {

private:
  vecgeom::Vector3D<double> fVertex;         ///< Vertex position
  int fEvent         = 0;                    ///< Event number
  std::atomic_int fNtracks;                  ///< Total number of tracks
  mutable std::vector<Track *> fPrimaries;   ///< Vector containing all primary tracks
public:
  /// Event default constructor
  Event() : fNtracks(0) {}

  /// Event destructor
  ~Event() { Clear(); }

  /// Function for accounting adding a new track.
  int AddTrack();

  /// Function for accounting adding a new primary track.
  int AddPrimary(Track *track)
  {
    fPrimaries.push_back(track);
    return AddTrack();
  }

  /// Clear the event and release all primaries.
  void Clear();

  /// Function for retrieving a primary. No range check.
  Track *GetPrimary(int i) { return fPrimaries[i]; }

  /// Get the number of primary tracks.
  int GetNprimaries() const { return fPrimaries.size(); }

  /// Function for reserving a number of primaries.
  void ReservePrimaries(int nprim) { fPrimaries.reserve(nprim); }

  /// Function that returns the event vertex.
  vecgeom::Vector3D<double> GetVertex() const { return fVertex; }

  /// Function that returns the event number.
  int GetEvent() const { return fEvent; }

  /// Function that returns the number of tracks.
  int GetNtracks() const { return fNtracks.load(); }

  /// Setter for the event number
  void SetEvent(int event) { fEvent = event; }

  /// Function to set the vertex.
  void SetVertex(double x, double y, double z) { fVertex.Set(x, y, z); }

  /// Print function.
  void Print(const char *option = "") const;
}; // Event

} // namespace vectorflow

#endif
