#include "TrackMLGeom.h"
#include "GeomData.h"

#include <iostream>
#include <sstream>
#include <string>

#include <VecGeom/base/Global.h>
#include <VecGeom/base/Vector3D.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/management/RootGeoManager.h>
#include <VecGeom/volumes/Box.h>
#include <VecGeom/volumes/LogicalVolume.h>
#include <VecGeom/volumes/Trd.h>
#include <VecGeom/volumes/Tube.h>
#include <VecGeom/navigation/VNavigator.h>
#include <VecGeom/navigation/NewSimpleNavigator.h>
#include <VecGeom/navigation/SimpleABBoxNavigator.h>
#include <VecGeom/navigation/SimpleABBoxLevelLocator.h>
#include <VecGeom/navigation/HybridNavigator2.h>

namespace trackml {

TrackMLGeom::TrackMLGeom() {}

TrackMLGeom::~TrackMLGeom() {}

//_____________________________________________________________________________
vecgeom::VPlacedVolume *TrackMLGeom::CreateGeometry(const char *filename) {
  using namespace vecgeom;
  using geant::units::mm;
  using vecCore::math::Sqrt;

  std::ifstream infile(filename);
  // Make sure the file is good
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

  // Create logical volume containers for the world volume, inner silicon pixels and two outer layers of silicon strips
  LogicalVolume *lWorld = CreateContainers();

  // Now read in the detailed module geometry
  std::string line;
  Module_t module;
  std::getline(infile, line); // skip header
  int n_modules = 0;
  Record_t rec, rec_ref;
  int layer_ref = -1;
  double zmin_layer, zmax_layer, rmin_layer, rmax_layer;
  while (std::getline(infile, line)) {
    GetRecord(line, rec);
    module.t = rec.module_t;
    module.minhu = rec.module_minhu;
    module.maxhu = rec.module_maxhu;
    module.hv = rec.module_hv;
    if (rec.layer_id != layer_ref) {
      if (layer_ref >= 0) {
        // Print info about the previous layer
        std::cout << "      volume " << rec_ref.volume_id << "  layer "
                  << layer_ref
                  << ": rad = " << 0.5 * (rmin_layer + rmax_layer) / mm
                  << " drad = " << 0.5 * (rmax_layer - rmin_layer) / mm
                  << " z = " << 0.5 * (zmin_layer + zmax_layer) / mm
                  << " dz = " << 0.5 * (zmax_layer - zmin_layer) / mm
                  << std::endl;
      }
      layer_ref = rec.layer_id;
      zmin_layer = 99999.;
      rmin_layer = 99999.;
      zmax_layer = -99999.;
      rmax_layer = -99999.;
    }
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
    double ztr = 0.;
    int layer_id = rec.layer_id / 2 - 1;
    switch (rec.volume_id) {
    case 7:
      ztr = kPixVol7LayersZ[layer_id];
      break;
    case 9:
      ztr = kPixVol9LayersZ[layer_id];
      break;
    case 12:
      ztr = kStrips1Vol12LayersZ[layer_id];
      break;
    case 14:
      ztr = kStrips1Vol14LayersZ[layer_id];
      break;
    case 16:
      ztr = kStrips2Vol16LayersZ[layer_id];
      break;
    case 18:
      ztr = kStrips2Vol18LayersZ[layer_id];
      break;
    }
    Transformation3D *tr = new Transformation3D(
        rec.module_trans[0], rec.module_trans[1], rec.module_trans[2] - ztr,
        rec.module_rot[0], rec.module_rot[1], rec.module_rot[2],
        rec.module_rot[3], rec.module_rot[4], rec.module_rot[5],
        rec.module_rot[6], rec.module_rot[7], rec.module_rot[8]);

    // Find position of layers
    double radius = Sqrt(rec.module_trans[0] * rec.module_trans[0] +
                         rec.module_trans[1] * rec.module_trans[1]);
    if (radius < rmin_layer)
      rmin_layer = radius;
    if (radius > rmax_layer)
      rmax_layer = radius;
    if (rec.module_trans[2] < zmin_layer)
      zmin_layer = rec.module_trans[2];
    if (rec.module_trans[2] > zmax_layer)
      zmax_layer = rec.module_trans[2];
    // Place the module in its mothe container
    fVolData[rec.volume_id].fLayers[layer_id]->PlaceDaughter(module.vol, tr);
    n_modules++;
  }
  // Print info about the last layer
  std::cout << "      volume " << rec_ref.volume_id << "  layer " << layer_ref
            << ": rad = " << 0.5 * (rmin_layer + rmax_layer) / mm
            << " drad = " << 0.5 * (rmax_layer - rmin_layer) / mm
            << " z = " << 0.5 * (zmin_layer + zmax_layer) / mm
            << " dz = " << 0.5 * (zmax_layer - zmin_layer) / mm << std::endl;
  std::cout << "Volume " << rec_ref.volume_id << ":  #module types = "
            << fVolData[rec_ref.volume_id].fModules.size()
            << "  #modules: " << n_modules << std::endl;

  VPlacedVolume *pWorld = lWorld->Place();
  GeoManager::Instance().SetWorldAndClose(pWorld);
  CreateNavigators();
  fWorld = pWorld;
  return pWorld;
}

//_____________________________________________________________________________
void TrackMLGeom::CreateNavigators()
{
// Create all navigators.
  for (auto &lvol : vecgeom::GeoManager::Instance().GetLogicalVolumesMap()) {
    if (lvol.second->GetDaughtersp()->size() < 4) {
      lvol.second->SetNavigator(vecgeom::NewSimpleNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 5) {
      lvol.second->SetNavigator(vecgeom::SimpleABBoxNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 10) {
      lvol.second->SetNavigator(vecgeom::HybridNavigator<>::Instance());
      vecgeom::HybridManager2::Instance().InitStructure((lvol.second));
    }
    lvol.second->SetLevelLocator(vecgeom::SimpleABBoxLevelLocator::GetInstance());
  } 
}

//_____________________________________________________________________________
bool TrackMLGeom::ExportToROOT(const char *filename) {
  return vecgeom::RootGeoManager::Instance().ExportToROOTGeometry(
      fWorld, std::string(filename));
}

//_____________________________________________________________________________
vecgeom::LogicalVolume *TrackMLGeom::CreateContainers() {
  using namespace vecgeom;
  // Create world volume.
  auto uWorld =
      GeoManager::MakeInstance<UnplacedBox>(kWorldDx, kWorldDy, kWorldDz);
  Transformation3D *ident = new Transformation3D();
  LogicalVolume *lWorld = new LogicalVolume("world", uWorld);

  // Pixels container volume (volumes 7, 8 and 9)
  auto uPixVol = GeoManager::MakeInstance<UnplacedTube>(
      kPixVolRmin, kPixVolRmax, kPixVolDz, 0., kTwoPi);
  auto lPixVol = new LogicalVolume("Pix", uPixVol);
  auto *pPixVol = lWorld->PlaceDaughter("Pix", lPixVol, ident);

  // Pixels volume 7 container
  auto uPixVol7 = GeoManager::MakeInstance<UnplacedTube>(
      kPixVol7Rmin, kPixVol7Rmax, kPixVol7Dz, 0., kTwoPi);
  auto lPixVol7 = new LogicalVolume("Pix7", uPixVol7);
  fVolData[7].fVol = lPixVol7;
  Transformation3D *trPixVol7 = new Transformation3D(0., 0., kPixVol7Z);
  auto *pPixVol7 = lPixVol->PlaceDaughter("Pix7", lPixVol7, trPixVol7);
  // Layer containers in Pixel volume 7 (disks)
  for (auto i = 0; i < kPixVol7Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer7 = GeoManager::MakeInstance<UnplacedTube>(
        kPixVol7Rmin, kPixVol7Rmax, kPixVol7LayerDz, 0., kTwoPi);
    auto lLayerVol7 = new LogicalVolume(name.str().c_str(), uLayer7);
    fVolData[7].fLayers.push_back(lLayerVol7);
    Transformation3D *trLayer =
        new Transformation3D(0., 0., kPixVol7LayersZ[i] - kPixVol7Z);
    auto *pPixLayer7 =
        lPixVol7->PlaceDaughter(name.str().c_str(), lLayerVol7, trLayer);
  }

  // Pixels volume 8 container
  auto uPixVol8 = GeoManager::MakeInstance<UnplacedTube>(
      kPixVol8Rmin, kPixVol8Rmax, kPixVol8Dz, 0., kTwoPi);
  auto lPixVol8 = new LogicalVolume("Pix8", uPixVol8);
  fVolData[8].fVol = lPixVol8;
  auto *pPixVol8 = lPixVol->PlaceDaughter("Pix8", lPixVol8, ident);
  // Layers in Pixel volume 8 (tube layers)
  for (auto i = 0; i < kPixVol8Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer8 = GeoManager::MakeInstance<UnplacedTube>(
        kPixVol8LayersR[i] - kPixVol8LayerDr,
        kPixVol8LayersR[i] + kPixVol8LayerDr, kPixVol8Dz, 0., kTwoPi);
    auto lLayerVol8 = new LogicalVolume(name.str().c_str(), uLayer8);
    fVolData[8].fLayers.push_back(lLayerVol8);
    auto *pPixLayer8 =
        lPixVol8->PlaceDaughter(name.str().c_str(), lLayerVol8, ident);
  }

  // Pixels volume 9 container
  auto uPixVol9 = GeoManager::MakeInstance<UnplacedTube>(
      kPixVol9Rmin, kPixVol9Rmax, kPixVol9Dz, 0., kTwoPi);
  auto lPixVol9 = new LogicalVolume("Pix9", uPixVol9);
  fVolData[9].fVol = lPixVol9;
  Transformation3D *trPixVol9 = new Transformation3D(0., 0., kPixVol9Z);
  auto *pPixVol9 = lPixVol->PlaceDaughter("Pix9", lPixVol9, trPixVol9);
  // Layers in Pixel volume 9 (disks)
  for (auto i = 0; i < kPixVol9Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer9 = GeoManager::MakeInstance<UnplacedTube>(
        kPixVol9Rmin, kPixVol9Rmax, kPixVol9LayerDz, 0., kTwoPi);
    auto lLayerVol9 = new LogicalVolume(name.str().c_str(), uLayer9);
    fVolData[9].fLayers.push_back(lLayerVol9);
    Transformation3D *trLayer =
        new Transformation3D(0., 0., kPixVol9LayersZ[i] - kPixVol9Z);
    auto *pPixLayer9 =
        lPixVol9->PlaceDaughter(name.str().c_str(), lLayerVol9, trLayer);
  }

  // Strips1 container (volumes 12, 13 and 14)
  auto uStrips1Vol = GeoManager::MakeInstance<UnplacedTube>(
      kStrips1VolRmin, kStrips1VolRmax, kStrips1VolDz, 0., kTwoPi);
  auto lStrips1Vol = new LogicalVolume("Strips1", uStrips1Vol);
  auto *pStrips1Vol = lWorld->PlaceDaughter("Strips1", lStrips1Vol, ident);

  // Volume 12 container
  auto uStrips1Vol12 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips1Vol12Rmin, kStrips1Vol12Rmax, kStrips1Vol12Dz, 0., kTwoPi);
  auto lStrips1Vol12 = new LogicalVolume("Strips1_12", uStrips1Vol12);
  fVolData[12].fVol = lStrips1Vol12;
  Transformation3D *trStrips1Vol12 =
      new Transformation3D(0., 0., kStrips1Vol12Z);
  auto *pStrips1Vol12 =
      lStrips1Vol->PlaceDaughter("Strips1_12", lStrips1Vol12, trStrips1Vol12);
  // Layers in Strips1 volume 12 (disks)
  for (auto i = 0; i < kStrips1Vol12Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer12 = GeoManager::MakeInstance<UnplacedTube>(
        kStrips1Vol12Rmin, kStrips1Vol12Rmax, kStrips1Vol12LayerDz, 0., kTwoPi);
    auto lLayerVol12 = new LogicalVolume(name.str().c_str(), uLayer12);
    fVolData[12].fLayers.push_back(lLayerVol12);
    Transformation3D *trLayer =
        new Transformation3D(0., 0., kStrips1Vol12LayersZ[i] - kStrips1Vol12Z);
    auto *pStrips1Layer12 =
        lStrips1Vol12->PlaceDaughter(name.str().c_str(), lLayerVol12, trLayer);
  }

  // Volume 13 container
  auto uStrips1Vol13 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips1Vol13Rmin, kStrips1Vol13Rmax, kStrips1Vol13Dz, 0., kTwoPi);
  auto lStrips1Vol13 = new LogicalVolume("Strips1_13", uStrips1Vol13);
  fVolData[13].fVol = lStrips1Vol13;
  auto *pStrips1Vol13 =
      lStrips1Vol->PlaceDaughter("Strips1_13", lStrips1Vol13, ident);
  // Layers in Strips1 volume 13 (tube layers)
  for (auto i = 0; i < kStrips1Vol13Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer13 = GeoManager::MakeInstance<UnplacedTube>(
        kStrips1Vol13LayersR[i] - kStrips1Vol13LayerDr,
        kStrips1Vol13LayersR[i] + kStrips1Vol13LayerDr, kStrips1Vol13Dz, 0.,
        kTwoPi);
    auto lLayerVol13 = new LogicalVolume(name.str().c_str(), uLayer13);
    fVolData[13].fLayers.push_back(lLayerVol13);
    auto *pStrips1Layer13 =
        lStrips1Vol13->PlaceDaughter(name.str().c_str(), lLayerVol13, ident);
  }

  // Volume 14 container
  auto uStrips1Vol14 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips1Vol14Rmin, kStrips1Vol14Rmax, kStrips1Vol14Dz, 0., kTwoPi);
  auto lStrips1Vol14 = new LogicalVolume("Strips1_14", uStrips1Vol14);
  fVolData[14].fVol = lStrips1Vol14;
  Transformation3D *trStrips1Vol14 =
      new Transformation3D(0., 0., kStrips1Vol14Z);
  auto *pStrips1Vol14 =
      lStrips1Vol->PlaceDaughter("Strips1_14", lStrips1Vol14, trStrips1Vol14);
  // Layers in Strips1 volume 14 (disks)
  for (auto i = 0; i < kStrips1Vol14Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer14 = GeoManager::MakeInstance<UnplacedTube>(
        kStrips1Vol14Rmin, kStrips1Vol14Rmax, kStrips1Vol14LayerDz, 0., kTwoPi);
    auto lLayerVol14 = new LogicalVolume(name.str().c_str(), uLayer14);
    fVolData[14].fLayers.push_back(lLayerVol14);
    Transformation3D *trLayer =
        new Transformation3D(0., 0., kStrips1Vol14LayersZ[i] - kStrips1Vol14Z);
    auto *pStrips1Layer14 =
        lStrips1Vol14->PlaceDaughter(name.str().c_str(), lLayerVol14, trLayer);
  }

  // Strips2 container (volumes 16, 17 and 18)
  auto uStrips2Vol = GeoManager::MakeInstance<UnplacedTube>(
      kStrips2VolRmin, kStrips2VolRmax, kStrips2VolDz, 0., kTwoPi);
  auto lStrips2Vol = new LogicalVolume("Strips2", uStrips2Vol);
  auto *pStrips2Vol = lWorld->PlaceDaughter("Strips2", lStrips2Vol, ident);

  // Volume 16 container
  auto uStrips2Vol16 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips2Vol16Rmin, kStrips2Vol16Rmax, kStrips2Vol16Dz, 0., kTwoPi);
  auto lStrips2Vol16 = new LogicalVolume("Strips2_16", uStrips2Vol16);
  fVolData[16].fVol = lStrips2Vol16;
  Transformation3D *trStrips2Vol16 =
      new Transformation3D(0., 0., kStrips2Vol16Z);
  auto *pStrips2Vol16 =
      lStrips2Vol->PlaceDaughter("Strips2_16", lStrips2Vol16, trStrips2Vol16);
  // Layers in Strips2 volume 16 (disks)
  for (auto i = 0; i < kStrips2Vol16Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer16 = GeoManager::MakeInstance<UnplacedTube>(
        kStrips2Vol16Rmin, kStrips2Vol16Rmax, kStrips2Vol16LayerDz, 0., kTwoPi);
    auto lLayerVol16 = new LogicalVolume(name.str().c_str(), uLayer16);
    fVolData[16].fLayers.push_back(lLayerVol16);
    Transformation3D *trLayer =
        new Transformation3D(0., 0., kStrips2Vol16LayersZ[i] - kStrips2Vol16Z);
    auto *pStrips2Layer16 =
        lStrips2Vol16->PlaceDaughter(name.str().c_str(), lLayerVol16, trLayer);
  }

  // Volume 17 container
  auto uStrips2Vol17 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips2Vol17Rmin, kStrips2Vol17Rmax, kStrips2Vol17Dz, 0., kTwoPi);
  auto lStrips2Vol17 = new LogicalVolume("Strips2_17", uStrips2Vol17);
  fVolData[17].fVol = lStrips2Vol17;
  auto *pStrips2Vol17 =
      lStrips2Vol->PlaceDaughter("Strips2_17", lStrips2Vol17, ident);
  // Layers in Strips2 volume 17 (tube layers)
  for (auto i = 0; i < kStrips2Vol17Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer17 = GeoManager::MakeInstance<UnplacedTube>(
        kStrips2Vol17LayersR[i] - kStrips2Vol17LayerDr,
        kStrips2Vol17LayersR[i] + kStrips2Vol17LayerDr, kStrips2Vol17Dz, 0.,
        kTwoPi);
    auto lLayerVol17 = new LogicalVolume(name.str().c_str(), uLayer17);
    fVolData[17].fLayers.push_back(lLayerVol17);
    auto *pStrips2Layer17 =
        lStrips2Vol17->PlaceDaughter(name.str().c_str(), lLayerVol17, ident);
  }

  // Volume 18 container
  auto uStrips2Vol18 = GeoManager::MakeInstance<UnplacedTube>(
      kStrips2Vol18Rmin, kStrips2Vol18Rmax, kStrips2Vol18Dz, 0., kTwoPi);
  auto lStrips2Vol18 = new LogicalVolume("Strips2_18", uStrips2Vol18);
  fVolData[18].fVol = lStrips2Vol18;
  Transformation3D *trStrips2Vol18 =
      new Transformation3D(0., 0., kStrips2Vol18Z);
  auto *pStrips2Vol18 =
      lStrips2Vol->PlaceDaughter("Strips2_18", lStrips2Vol18, trStrips2Vol18);
  // Layers in Strips2 volume 18 (disks)
  for (auto i = 0; i < kStrips2Vol18Nlayers; ++i) {
    std::stringstream name;
    name << "Layer" << 2 * i;
    auto uLayer18 = GeoManager::MakeInstance<UnplacedTube>(
        kStrips2Vol18Rmin, kStrips2Vol18Rmax, kStrips2Vol18LayerDz, 0., kTwoPi);
    auto lLayerVol18 = new LogicalVolume(name.str().c_str(), uLayer18);
    fVolData[18].fLayers.push_back(lLayerVol18);
    Transformation3D *trLayer =
        new Transformation3D(0., 0., kStrips2Vol18LayersZ[i] - kStrips2Vol18Z);
    auto *pStrips2Layer18 =
        lStrips2Vol18->PlaceDaughter(name.str().c_str(), lLayerVol18, trLayer);
  }
  return lWorld;
}

//_____________________________________________________________________________
void TrackMLGeom::GetRecord(std::string const &line, Record_t &rec) const {
  using geant::units::mm;
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
  rec.module_trans[0] = std::stod(value) * mm;
  std::getline(linestream, value, ',');
  rec.module_trans[1] = std::stod(value) * mm;
  std::getline(linestream, value, ',');
  rec.module_trans[2] = std::stod(value) * mm;
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
  rec.module_t = std::stod(value) * mm;
  std::getline(linestream, value, ',');
  rec.module_minhu = std::stod(value) * mm;
  std::getline(linestream, value, ',');
  rec.module_maxhu = std::stod(value) * mm;
  std::getline(linestream, value, ',');
  rec.module_hv = std::stod(value) * mm;
  std::getline(linestream, value, ',');
  rec.pitch_u = std::stod(value) * mm;
  std::getline(linestream, value, ',');
  rec.pitch_v = std::stod(value) * mm;
}

} // namespace trackml