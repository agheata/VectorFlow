#ifndef TRACKML_TRACKER_GEOM
#define TRACKML_TRACKER_GEOM

#include <vectorFlow/SystemOfUnits.h>
#include <vectorFlow/Typedefs.h>

namespace trackml {

/// Geometry representation of the TrackML tracker detector.
/**
  <a href="https://asalzbur.web.cern.ch/asalzbur/work/tml/Detector.png"></a>

  The detector is built from silicon slabs (or modules, rectangular or
  trapezoïdal), arranged in cylinders and disks, which measure the position (or
  hits) of the particles that cross them. The detector modules are organized
  into detector groups or volumes identified by a volume id. Inside a volume
  they are further grouped into layers identified by a layer id. Each layer can
  contain an arbitrary number of detector modules, the smallest geometrically
  distinct detector object, each identified by a module_id. Within each group,
  detector modules are of the same type have e.g. the same granularity. All
  simulated detector modules are so-called semiconductor sensors that are build
  from thin silicon sensor chips. Each module can be represented by a
  two-dimensional, planar, bounded sensitive surface. These sensitive surfaces
  are subdivided into regular grids that define the detectors cells, the
  smallest granularity within the detector.

  <a
  href="https://asalzbur.web.cern.ch/asalzbur/work/tml/localToGlobal.png"></a>

  Each module has a different position and orientation described in the
  detectors file. A local, right-handed coordinate system is defined on each
  sensitive surface such that the first two coordinates u and v are on the
  sensitive surface and the third coordinate w is normal to the surface. The
  orientation and position are defined by the following transformation:

    pos_xyz = rotation_matrix * pos_uvw + translation

  that transform a position described in local coordinates u,v,w into the
  equivalent position x,y,z in global coordinates using a rotation matrix and an
  translation vector (cx,cy,cz).

  - volume_id: numerical identifier of the detector group.
  - layer_id: numerical identifier of the detector layer inside the group.
  - module_id: numerical identifier of the detector module inside the layer.
  - cx, cy, cz: position of the local origin in the described in the global
  coordinate system (in millimeter).
  - rot_xu, rot_xv, rot_xw, rot_yu, ...: components of the rotation matrix to
  rotate from local u,v,w to global x,y,z coordinates.
  - module_t: half thickness of the detector module (in millimeter).
  - module_minhu, module_maxhu: the minimum/maximum half-length of the module
  boundary along the local u direction (in millimeter).
  - module_hv: the half-length of the module boundary along the local v
  direction (in millimeter).
  - pitch_u, pitch_v: the size of detector cells along the local u and v
  direction (in millimeter).

  There are two different module shapes in the detector, rectangular and
  trapezoidal. The pixel detector ( with volume_id = 7,8,9) is fully built from
  rectangular modules, and so are the cylindrical barrels in volume_id=13,17.
  The remaining layers are made out disks that need trapezoidal shapes to cover
  the full disk.

  <a href="https://asalzbur.web.cern.ch/asalzbur/work/tml/ModuleTypes.png"></a>
  */

/// Structure representing a record in the geometry csv file
struct Record_t {
  int volume_id = 0; ///< Numerical identifier of the detector group.
  int layer_id =
      0; ///< Numerical identifier of the detector layer inside the group.
  int module_id =
      0; ///< Numerical identifier of the detector module inside the layer.
  double module_trans[3] = {
      0.}; ///< Position of the local origin in the described in the global
           ///< coordinate system (in millimeter)
  double module_rot[9] = {
      0.}; ///< Components of the rotation matrix to rotate from local u,v,w to
           ///< global x,y,z coordinates.
  double module_t =
      0.; ///< Half thickness of the detector module (in millimeter).
  double module_minhu = 0.; ///< The minimum half-length of the module boundary
                            ///< along the local u direction (in millimeter).
  double module_maxhu = 0.; ///< The maximum half-length of the module boundary
                            ///< along the local u direction (in millimeter).
  double module_hv = 0.; ///< The half-length of the module boundary along the
                         ///< local v direction (in millimeter).
  double pitch_u = 0.;   ///< The size of detector cells along the local u
                         ///< direction (in millimeter).
  double pitch_v = 0.;   ///< The size of detector cells along the local v
                         ///< direction (in millimeter).
};

/// Structure storing the dimensions of one module, used to create a single
/// logical volume for one set of dimensions.
struct Module_t {
  double t = 0.;     ///< Half thickness of the detector module (in millimeter).
  double minhu = 0.; ///< The minimum half-length of the module boundary along
                     ///< the local u direction (in millimeter).
  double maxhu = 0.; ///< The maximum half-length of the module boundary along
                     ///< the local u direction (in millimeter).
  double hv = 0.; ///< The half-length of the module boundary along the local v
                  ///< direction (in millimeter).
  vecgeom::LogicalVolume *vol =
      nullptr; ///< Pointer to corresponding logical volume

  bool operator==(const Module_t &other) const {
    return !(t != other.t || minhu != other.minhu || maxhu != other.maxhu ||
             hv != other.hv);
  }
};

/// Structure storing the hit u/v errors for hits in one volume, matching the
/// size of the pixels for the modules in that volume
struct ErrorUV_t {
  double error_u = 0.;
  double error_v = 0.;
};

struct VolumeData_t {
  std::vector<Module_t> fModules; ///< List of different modules in the volume
  vecgeom::LogicalVolume *fVol =
      nullptr;      ///< Pointer to corresponding logical volume
  ErrorUV_t fErrUV; ///< UV components of pixel errors

  bool HasModule(Module_t &mod) const {
    auto it = std::find(fModules.begin(), fModules.end(), mod);
    if (it == fModules.end())
      return false;
    mod = *it;
    return true;
  }
};

/// Class creating the VecGeom geometry representation of TrackML tracker
/// detector
class TrackMLGeom {
public:
  // using namespace geant::units;
  using LogicalVolume = vecgeom::LogicalVolume;
  using VPlacedVolume = vecgeom::VPlacedVolume;
  using GeoManager = vecgeom::GeoManager;

  TrackMLGeom();
  ~TrackMLGeom();

  TrackMLGeom(TrackMLGeom const &) = delete;
  TrackMLGeom const &operator=(TrackMLGeom const &) = delete;

  VPlacedVolume *CreateGeometry(const char *filename = "detectors.csv");
  bool ExportToROOT(const char *filename = "trackML.root");

private:
  VPlacedVolume *fWorld;     ///< The vecgeom world volume
  VolumeData_t fVolData[20]; ///< Array of logical volumes for each detector
  ErrorUV_t fVolPixErr[20];  ///< Array of pixel errors for each detector

  /// Fill one record from one line read from the geometry file.
  void GetRecord(std::string const &line, Record_t &rec) const;
};
} // namespace trackml
#endif
