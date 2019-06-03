#ifndef TRACKML_TRACKER_GEOM_DATA
#define TRACKML_TRACKER_GEOM_DATA

#include <vectorFlow/SystemOfUnits.h>
#include <vectorFlow/Typedefs.h>

namespace trackml {

/// Geometry data for the TrackML tracker detector.
// World volume size
static constexpr double kWorldDx = 1100. * geant::units::mm;
static constexpr double kWorldDy = 1100. * geant::units::mm;
static constexpr double kWorldDz = 3100. * geant::units::mm;

// Mother volume for pixel detector
static constexpr double kPixVolDz = 1830. * geant::units::mm;
static constexpr double kPixVolRmin = 25. * geant::units::mm;
static constexpr double kPixVolRmax = 200. * geant::units::mm;

// Pixel volume 7
static constexpr double kPixVol7Z = -1130. * geant::units::mm;
static constexpr double kPixVol7Dz = 1830. * geant::units::mm;
static constexpr double kPixVol7Rmin = 25. * geant::units::mm;
static constexpr double kPixVol7Rmax = 200. * geant::units::mm;

// Pixel volume 8
static constexpr double kPixVol8Z = 0. * geant::units::mm;
static constexpr double kPixVol8Dz = 450. * geant::units::mm;
static constexpr double kPixVol8Rmin = 25. * geant::units::mm;
static constexpr double kPixVol8Rmax = 200. * geant::units::mm;

// Pixel volume 9
static constexpr double kPixVol9Z = 1130. * geant::units::mm;
static constexpr double kPixVol9Dz = 1830. * geant::units::mm;
static constexpr double kPixVol9Rmin = 25. * geant::units::mm;
static constexpr double kPixVol9Rmax = 200. * geant::units::mm;

// Mother volume for Strips1 detector
static constexpr double kStrips1VolDz = 3000. * geant::units::mm;
static constexpr double kStrips1VolRmin = 225. * geant::units::mm;
static constexpr double kStrips1VolRmax = 730. * geant::units::mm;

// Strips1 volume 12
static constexpr double kStrips1Vol12Z = -2075. * geant::units::mm;
static constexpr double kStrips1Vol12Dz = 925. * geant::units::mm;
static constexpr double kStrips1Vol12Rmin = 250. * geant::units::mm;
static constexpr double kStrips1Vol12Rmax = 730. * geant::units::mm;

// Strips1 volume 13
static constexpr double kStrips1Vol13Z = 0. * geant::units::mm;
static constexpr double kStrips1Vol13Dz = 1140. * geant::units::mm;
static constexpr double kStrips1Vol13Rmin = 250. * geant::units::mm;
static constexpr double kStrips1Vol13Rmax = 730. * geant::units::mm;

// Strips1 volume 14
static constexpr double kStrips1Vol14Z = 2075. * geant::units::mm;
static constexpr double kStrips1Vol14Dz = 925. * geant::units::mm;
static constexpr double kStrips1Vol14Rmin = 250. * geant::units::mm;
static constexpr double kStrips1Vol14Rmax = 730. * geant::units::mm;

// Mother volume for Strips2 detector
static constexpr double kStrips2VolDz = 3000. * geant::units::mm;
static constexpr double kStrips2VolRmin = 750. * geant::units::mm;
static constexpr double kStrips2VolRmax = 1050. * geant::units::mm;

// Strips2 volume 16
static constexpr double kStrips2Vol16Z = -2075. * geant::units::mm;
static constexpr double kStrips2Vol16Dz = 925. * geant::units::mm;
static constexpr double kStrips2Vol16Rmin = 750. * geant::units::mm;
static constexpr double kStrips2Vol16Rmax = 1050. * geant::units::mm;

// Strips2 volume 17
static constexpr double kStrips2Vol17Z = 0. * geant::units::mm;
static constexpr double kStrips2Vol17Dz = 1140. * geant::units::mm;
static constexpr double kStrips2Vol17Rmin = 750. * geant::units::mm;
static constexpr double kStrips2Vol17Rmax = 1050. * geant::units::mm;

// Strips2 volume 18
static constexpr double kStrips2Vol18Z = 2075. * geant::units::mm;
static constexpr double kStrips2Vol18Dz = 925. * geant::units::mm;
static constexpr double kStrips2Vol18Rmin = 750. * geant::units::mm;
static constexpr double kStrips2Vol18Rmax = 1050. * geant::units::mm;
} // namespace trackml
#endif
