#ifndef TRACKML_TRACKER_GEOM_DATA
#define TRACKML_TRACKER_GEOM_DATA

#include <vectorFlow/SystemOfUnits.h>
#include <vectorFlow/Typedefs.h>

namespace trackml {

using geant::units::mm;
/// Geometry data for the TrackML tracker detector.
// World volume size
static constexpr double kWorldDx = 1100 * mm;
static constexpr double kWorldDy = 1100 * mm;
static constexpr double kWorldDz = 3100 * mm;

// Mother volume for pixel detector
static constexpr double kPixVolDz = 1830 * mm;
static constexpr double kPixVolRmin = 25 * mm;
static constexpr double kPixVolRmax = 200 * mm;

// Pixel volume 7
static constexpr double kPixVol7Z = -1130 * mm;
static constexpr double kPixVol7Dz = 670 * mm;
static constexpr double kPixVol7Rmin = 25 * mm;
static constexpr double kPixVol7Rmax = 200 * mm;
static constexpr int kPixVol7Nlayers = 7;
static constexpr double kPixVol7LayerDz = 3 * mm;
static constexpr double kPixVol7LayersZ[7] = {
    -1500 * mm, -1300 * mm, -1100 * mm, -960 * mm,
    -820 * mm,  -700 * mm,  -600 * mm};

// Pixel volume 8
static constexpr double kPixVol8Z = 0 * mm;
static constexpr double kPixVol8Dz = 500 * mm;
static constexpr double kPixVol8Rmin = 25 * mm;
static constexpr double kPixVol8Rmax = 200 * mm;
static constexpr int kPixVol8Nlayers = 4;
static constexpr double kPixVol8LayerDr = 5 * mm;
static constexpr double kPixVol8LayersR[4] = {32 * mm, 72 * mm, 116 * mm,
                                              172 * mm};

// Pixel volume 9
static constexpr double kPixVol9Z = 1130 * mm;
static constexpr double kPixVol9Dz = 670 * mm;
static constexpr double kPixVol9Rmin = 25 * mm;
static constexpr double kPixVol9Rmax = 200 * mm;
static constexpr int kPixVol9Nlayers = 7;
static constexpr double kPixVol9LayerDz = 3 * mm;
static constexpr double kPixVol9LayersZ[7] = {
    600 * mm, 700 * mm, 820 * mm, 960 * mm, 1100 * mm, 1300 * mm, 1500 * mm};

// Mother volume for Strips1 detector
static constexpr double kStrips1VolDz = 3000 * mm;
static constexpr double kStrips1VolRmin = 225 * mm;
static constexpr double kStrips1VolRmax = 730 * mm;

// Strips1 volume 12
static constexpr double kStrips1Vol12Z = -2075 * mm;
static constexpr double kStrips1Vol12Dz = 925 * mm;
static constexpr double kStrips1Vol12Rmin = 235 * mm;
static constexpr double kStrips1Vol12Rmax = 730 * mm;
static constexpr int kStrips1Vol12Nlayers = 6;
static constexpr double kStrips1Vol12LayerDz = 5 * mm;
static constexpr double kStrips1Vol12LayersZ[6] = {
    -2950 * mm, -2550 * mm, -2150 * mm, -1800 * mm, -1500 * mm, -1220 * mm};

// Strips1 volume 13
static constexpr double kStrips1Vol13Z = 0 * mm;
static constexpr double kStrips1Vol13Dz = 1140 * mm;
static constexpr double kStrips1Vol13Rmin = 235 * mm;
static constexpr double kStrips1Vol13Rmax = 730 * mm;
static constexpr int kStrips1Vol13Nlayers = 4;
static constexpr double kStrips1Vol13LayerDr = 7 * mm;
static constexpr double kStrips1Vol13LayersR[4] = {260 * mm, 360 * mm, 500 * mm,
                                                   660 * mm};

// Strips1 volume 14
static constexpr double kStrips1Vol14Z = 2075 * mm;
static constexpr double kStrips1Vol14Dz = 925 * mm;
static constexpr double kStrips1Vol14Rmin = 235 * mm;
static constexpr double kStrips1Vol14Rmax = 730 * mm;
static constexpr int kStrips1Vol14Nlayers = 6;
static constexpr double kStrips1Vol14LayerDz = 5 * mm;
static constexpr double kStrips1Vol14LayersZ[6] = {
    1220 * mm, 1500 * mm, 1800 * mm, 2150 * mm, 2550 * mm, 2950 * mm};

// Mother volume for Strips2 detector
static constexpr double kStrips2VolDz = 3000 * mm;
static constexpr double kStrips2VolRmin = 750 * mm;
static constexpr double kStrips2VolRmax = 1050 * mm;

// Strips2 volume 16
static constexpr double kStrips2Vol16Z = -2075 * mm;
static constexpr double kStrips2Vol16Dz = 925 * mm;
static constexpr double kStrips2Vol16Rmin = 750 * mm;
static constexpr double kStrips2Vol16Rmax = 1050 * mm;
static constexpr int kStrips2Vol16Nlayers = 6;
static constexpr double kStrips2Vol16LayerDz = 6 * mm;
static constexpr double kStrips2Vol16LayersZ[6] = {
    -2950 * mm, -2550 * mm, -2150 * mm, -1800 * mm, -1500 * mm, -1220 * mm};

// Strips2 volume 17
static constexpr double kStrips2Vol17Z = 0 * mm;
static constexpr double kStrips2Vol17Dz = 1140 * mm;
static constexpr double kStrips2Vol17Rmin = 750 * mm;
static constexpr double kStrips2Vol17Rmax = 1050 * mm;
static constexpr int kStrips2Vol17Nlayers = 2;
static constexpr double kStrips2Vol17LayerDr = 5 * mm;
static constexpr double kStrips2Vol17LayersR[2] = {820 * mm, 1020 * mm};

// Strips2 volume 18
static constexpr double kStrips2Vol18Z = 2075 * mm;
static constexpr double kStrips2Vol18Dz = 925 * mm;
static constexpr double kStrips2Vol18Rmin = 750 * mm;
static constexpr double kStrips2Vol18Rmax = 1050 * mm;
static constexpr int kStrips2Vol18Nlayers = 6;
static constexpr double kStrips2Vol18LayerDz = 6 * mm;
static constexpr double kStrips2Vol18LayersZ[6] = {
    1220 * mm, 1500 * mm, 1800 * mm, 2150 * mm, 2550 * mm, 2950 * mm};
} // namespace trackml
#endif
