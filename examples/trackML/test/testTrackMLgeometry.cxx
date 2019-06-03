#include "TrackMLGeom.h"
#include "base/Vector3D.h"
#include "vectorFlow/PhysicalConstants.h"
#include "vectorFlow/SystemOfUnits.h"
#include <iostream>
/// This executable creates the TrackMK geometry representation in VecGeom, then
/// exports it to a ROOT file.
/** To visualize the geometry in ROOT (ROOT must be compiled with OpenGL
  support: cd $VECTORFLOW_INSRALL/bin/tests
  ./testTrackMLgeometry
  root
  (inside the root seddion)
  root [0] auto geom = TGeoManager::Import("trackML.root");
  root [1] geom->DefaultColors();
  root [2] geom->SetVisLevel(4);
  root [3] geom->GetTopVolume()->Draw("ogl");
*/

int main() {
  trackml::TrackMLGeom *geom = new trackml::TrackMLGeom();
  auto world = geom->CreateGeometry("data/detectors.csv");

  if (!world) {
    std::cout << "example testTrackMLgeometry FAILED\n";
    delete geom;
    return 1;
  }
  geom->ExportToROOT();
  delete geom;
  std::cout << "example testTrackMLgeometry OK\n";
  return 0;
}
