Assuming $VFLOW_ROOT is the base installation directory. The CMAKE_BUILD_TYPE can be changed to Debug or Release.
The following dependencies have to be installed first:

#veccore
cd $VFLOW_ROOT
git clone https://github.com/root-project/veccore.git
mkdir $VFLOW_ROOT/veccore_build && cd $VFLOW_ROOT/veccore_build
cmake $GEANT_ROOT/veccore -DCMAKE_INSTALL_PREFIX=$VFLOW_ROOT/veccore_install -DBUILD_VC=ON -DBUILD_UMESIMD=ON
make -j install

#vecmath
cd $VFLOW_ROOT
git clone https://github.com/root-project/vecmath.git
mkdir $VFLOW_ROOT/vecmath_build && cd $VFLOW_ROOT/vecmath_build
cmake $VFLOW_ROOT/vecmath -DCMAKE_INSTALL_PREFIX=$VFLOW_ROOT/vecmath_install -DCMAKE_BUILD_TYPE=RelWithDebInfo -DVecCore_DIR=$VFLOW_ROOT/veccore_install/share/VecCore/cmake

#vecgeom
cd $VFLOW_ROOT
git clone https://gitlab.cern.ch/VecGeom/VecGeom.git
mkdir $VFLOW_ROOT/vecgeom_build && cd $VFLOW_ROOT/vecgeom_build
cmake  $VFLOW_ROOT/VecGeom -DCMAKE_INSTALL_PREFIX=$VFLOW_ROOT/vecgeom_install -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_PREFIX_PATH=$VFLOW_ROOT/veccore_install -DBACKEND=Vc

The package can be installed using:

cd $VFLOW_ROOT
mkdir $VFLOW_ROOT/vectorflow_build && cd $VFLOW_ROOT/vectorflow_build
cmake $VFLOW_ROOT/vectorflow -DCMAKE_INSTALL_PREFIX=$VFLOW_ROOT/vectorflow_install -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_PREFIX_PATH=$VFLOW_ROOT/veccore_install -DVecCore_DIR=$VFLOW_ROOT/veccore_install/share/VecCore/cmake -DVecMath_DIR=$VFLOW_ROOT/vecmath_install/lib/cmake/VecMath -DVecGeom_DIR=$VFLOW_ROOT/vecgeom_install/lib/cmake/VecGeom
