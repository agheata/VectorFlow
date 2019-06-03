#include "TrackMLGeom.h"
#include "GeomData.h"

#include <iostream>
#include <sstream>
#include <string>

#include <base/Global.h>
#include <base/Vector3D.h>
#include <management/GeoManager.h>
#include <management/RootGeoManager.h>
#include <volumes/Box.h>
#include <volumes/LogicalVolume.h>
#include <volumes/Trd.h>
#include <volumes/Tube.h>

namespace trackml {

TrackMLGeom::TrackMLGeom() {}

TrackMLGeom::~TrackMLGeom() {}

vecgeom::VPlacedVolume *TrackMLGeom::CreateGeometry(const char *filename) {
  using namespace vecgeom;
  std::ifstream infile(filename);
  try {
    if (!infile.good()) {
      std::string err =
          "=== ERROR === CreateGeometry: Cannot open input file: ";
      err += filename;
      throw std::runtime_error(err);
    }
  } catch (const std::exception &e) {
    std::cout << e.what() << "\n";
    return nullptr;
  }
  // Create world volume.
  auto uWorld =
      GeoManager::MakeInstance<UnplacedBox>(kWorldDx, kWorldDy, kWorldDz);
  Transformation3D *ident = new Transformation3D();
  LogicalVolume *lWorld = new LogicalVolume("world", uWorld);

  // Pixels
  auto uPixVol = GeoManager::MakeInstance<UnplacedTube>(
      kPixVolRmin, kPixVolRmax, kPixVolDz, 0., kTwoPi);
  auto lPixVol = new LogicalVolume("Pix", uPixVol);
  auto *pPixVol = lWorld->PlaceDaughter("Pix", lPixVol, ident);

  auto uPixVol7 = GeoManager::MakeInstance<UnplacedTube>(
      kPixVol7Rmin, kPixVol7Rmax, kPixVol7Dz, 0., kTwoPi);
  auto lPixVol7 = new LogicalVolume("Pix7", uPixVol7);
  fVolData[7].fVol = lPixVol7;
  Transformation3D *trPixVol7 = new Transformation3D(0., 0., kPixVol7Z);
  auto *pPixVol7 = lPixVol->PlaceDaughter("Pix7", lPixVol7, trPixVol7);

  auto uPixVol8 = GeoManager::MakeInstance<UnplacedTube>(
      kPixVol8Rmin, kPixVol8Rmax, kPixVol8Dz, 0., kTwoPi);
  auto lPixVol8 = new LogicalVolume("Pix8", uPixVol8);
  fVolData[8].fVol = lPixVol8;
  auto *pPixVol8 = lPixVol->PlaceDaughter("Pix8", lPixVol8, ident);

  auto uPixVol9 = GeoManager::MakeInstance<UnplacedTube>(
      kPixVol9Rmin, kPixVol9Rmax, kPixVol9Dz, 0., kTwoPi);
  auto lPixVol9 = new LogicalVolume("Pix9", uPixVol9);
  fVolData[9].fVol = lPixVol9;
  Transformation3D *trPixVol9 = new Transformation3D(0., 0., kPixVol9Z);
  auto *pPixVol9 = lPixVol->PlaceDaughter("Pix9", lPixVol9, trPixVol9);

  // Strips1
  auto uStrips1Vol = GeoManager::MakeInstance<UnplacedTube>(
      kStrips1VolRmin, kStrips1VolRmax, kStrips1VolDz, 0., kTwoPi);
  auto lStrips1Vol = new LogicalVolume("Strips1", uStrips1Vol);
  auto *pStrips1Vol = lWorld->PlaceDaughter("Strips1", lStrips1Vol, ident);

  auto uStrips1Vol12 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips1Vol12Rmin, kStrips1Vol12Rmax, kStrips1Vol12Dz, 0., kTwoPi);
  auto lStrips1Vol12 = new LogicalVolume("Strips1_12", uStrips1Vol12);
  fVolData[12].fVol = lStrips1Vol12;
  Transformation3D *trStrips1Vol12 =
      new Transformation3D(0., 0., kStrips1Vol12Z);
  auto *pStrips1Vol12 =
      lStrips1Vol->PlaceDaughter("Strips1_12", lStrips1Vol12, trStrips1Vol12);

  auto uStrips1Vol13 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips1Vol13Rmin, kStrips1Vol13Rmax, kStrips1Vol13Dz, 0., kTwoPi);
  auto lStrips1Vol13 = new LogicalVolume("Strips1_13", uStrips1Vol13);
  fVolData[13].fVol = lStrips1Vol13;
  auto *pStrips1Vol13 =
      lStrips1Vol->PlaceDaughter("Strips1_13", lStrips1Vol13, ident);

  auto uStrips1Vol14 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips1Vol14Rmin, kStrips1Vol14Rmax, kStrips1Vol14Dz, 0., kTwoPi);
  auto lStrips1Vol14 = new LogicalVolume("Strips1_14", uStrips1Vol14);
  fVolData[14].fVol = lStrips1Vol14;
  Transformation3D *trStrips1Vol14 =
      new Transformation3D(0., 0., kStrips1Vol14Z);
  auto *pStrips1Vol14 =
      lStrips1Vol->PlaceDaughter("Strips1_14", lStrips1Vol14, trStrips1Vol14);

  // Strips2
  auto uStrips2Vol = GeoManager::MakeInstance<UnplacedTube>(
      kStrips2VolRmin, kStrips2VolRmax, kStrips2VolDz, 0., kTwoPi);
  auto lStrips2Vol = new LogicalVolume("Strips2", uStrips2Vol);
  auto *pStrips2Vol = lWorld->PlaceDaughter("Strips2", lStrips2Vol, ident);

  auto uStrips2Vol16 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips2Vol16Rmin, kStrips2Vol16Rmax, kStrips2Vol16Dz, 0., kTwoPi);
  auto lStrips2Vol16 = new LogicalVolume("Strips2_16", uStrips2Vol16);
  fVolData[16].fVol = lStrips2Vol16;
  Transformation3D *trStrips2Vol16 =
      new Transformation3D(0., 0., kStrips2Vol16Z);
  auto *pStrips2Vol16 =
      lStrips2Vol->PlaceDaughter("Strips2_16", lStrips2Vol16, trStrips2Vol16);

  auto uStrips2Vol17 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips2Vol17Rmin, kStrips2Vol17Rmax, kStrips2Vol17Dz, 0., kTwoPi);
  auto lStrips2Vol17 = new LogicalVolume("Strips2_17", uStrips2Vol17);
  fVolData[17].fVol = lStrips2Vol17;
  auto *pStrips2Vol17 =
      lStrips2Vol->PlaceDaughter("Strips2_17", lStrips2Vol17, ident);

  auto uStrips2Vol18 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips2Vol18Rmin, kStrips2Vol18Rmax, kStrips2Vol18Dz, 0., kTwoPi);
  auto lStrips2Vol18 = new LogicalVolume("Strips2_18", uStrips2Vol18);
  fVolData[18].fVol = lStrips2Vol18;
  Transformation3D *trStrips2Vol18 =
      new Transformation3D(0., 0., kStrips2Vol18Z);
  auto *pStrips2Vol18 =
      lStrips2Vol->PlaceDaughter("Strips2_18", lStrips2Vol18, trStrips2Vol18);

  // Now read in the detailed module geometry
  std::string line;
  Module_t module;
  std::getline(infile, line); // skip header
  int volume_ref = 0;
  int n_modules = 0;
  Record_t rec, rec_ref;
  while (std::getline(infile, line)) {
    GetRecord(line, rec);
    module.t = rec.module_t;
    module.minhu = rec.module_minhu;
    module.maxhu = rec.module_maxhu;
    module.hv = rec.module_hv;
    if (rec.volume_id != rec_ref.volume_id) {
      // This is the first record for a given volume: print some statistics for
      // the previous volume
      if (rec_ref.volume_id > 0) {
        std::cout << "Volume " << rec_ref.volume_id << ":  #modules types = "
                  << fVolData[rec_ref.volume_id].fModules.size()
                  << "  #modules: " << n_modules << std::endl;
      }
      n_modules = 0;
      rec_ref = rec;
      fVolPixErr[rec.volume_id].error_u = rec.pitch_u;
      fVolPixErr[rec.volume_id].error_v = rec.pitch_v;
    }

    // Check if the module/pixels have same geometry as the reference
    if (!fVolData[rec.volume_id].HasModule(module)) {
      std::stringstream modName;
      modName << "Module" << rec.volume_id << "_"
              << fVolData[rec.volume_id].fModules.size();

      LogicalVolume *lMod = nullptr;
      if (module.minhu == module.maxhu) {
        // Rectangular module
        auto uMod = GeoManager::MakeInstance<UnplacedBox>(module.minhu,
                                                          module.hv, module.t);
        lMod = new LogicalVolume(modName.str().c_str(), uMod);
      } else {
        // Trapezoidal module
        auto uMod = GeoManager::MakeInstance<UnplacedTrd>(
            module.minhu, module.maxhu, module.hv, module.t);
        lMod = new LogicalVolume(modName.str().c_str(), uMod);
      }
      module.vol = lMod;
      fVolData[rec.volume_id].fModules.push_back(module);
    }
    Transformation3D *tr = new Transformation3D(
        rec.module_trans[0], rec.module_trans[1], rec.module_trans[2],
        rec.module_rot[0], rec.module_rot[1], rec.module_rot[2],
        rec.module_rot[3], rec.module_rot[4], rec.module_rot[5],
        rec.module_rot[6], rec.module_rot[7], rec.module_rot[8]);
    // Place the module in its mothe container
    fVolData[rec.volume_id].fVol->PlaceDaughter(module.vol, tr);
    n_modules++;
  }
  std::cout << "Volume " << rec_ref.volume_id << ":  #module types = "
            << fVolData[rec_ref.volume_id].fModules.size()
            << "  #modules: " << n_modules << std::endl;

  VPlacedVolume *pWorld = lWorld->Place();
  GeoManager::Instance().SetWorldAndClose(pWorld);
  fWorld = pWorld;
  return pWorld;
}

bool TrackMLGeom::ExportToROOT(const char *filename) {
  return vecgeom::RootGeoManager::Instance().ExportToROOTGeometry(
      fWorld, std::string(filename));
}

void TrackMLGeom::GetRecord(std::string const &line, Record_t &rec) const {
  std::stringstream linestream(line);
  std::string value;
  // volume_id,layer_id,module_id,cx,cy,cz,rot_xu,rot_xv,rot_xw,rot_yu,rot_yv,rot_yw,rot_zu,rot_zv,rot_zw,module_t,module_minhu,module_maxhu,module_hv,pitch_u,pitch_v
  std::getline(linestream, value, ',');
  rec.volume_id = std::stoi(value);
  std::getline(linestream, value, ',');
  rec.layer_id = std::stoi(value);
  std::getline(linestream, value, ',');
  rec.module_id = std::stoi(value);
  std::getline(linestream, value, ',');
  rec.module_trans[0] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_trans[1] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_trans[2] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[0] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[1] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[2] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[3] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[4] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[5] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[6] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[7] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_rot[8] = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_t = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_minhu = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_maxhu = std::stod(value);
  std::getline(linestream, value, ',');
  rec.module_hv = std::stod(value);
  std::getline(linestream, value, ',');
  rec.pitch_u = std::stod(value);
  std::getline(linestream, value, ',');
  rec.pitch_v = std::stod(value);
}

} // namespace trackml