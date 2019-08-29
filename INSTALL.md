**VectorFlow** has package dependencies that need to be installed first, those are VECCORE, VECMATH, and VECGEOM. Also, these packages has other dependencies that need to be addressed too. Here is a list:

## Tools needed for multiple purposes

* sudo apt-get install git dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev libgd-dev doxygen curl graphviz

## ROOT (VecGeom dependency)

For Ubuntu OS users, the binaries can be found:

* wget https://root.cern.ch/download/root_v6.16.00.Linux-ubuntu18-x86_64-gcc7.3.tar.gz
* tar -zxf root_v6.16.00.Linux-ubuntu18-x86_64-gcc7.3.tar.gz

**NOTE 1:** To use ROOT, add the libraries to your LD_LIBRARY_PATH:

* export LD_LIBRARY_PATH=$LIB_PATH:$VFLOW_ROOT/root/lib

**NOTE 2:** This is a version of ROOT, for more information and/or different versions, please visit:

* https://root.cern.ch/content/release-61600
* https://root.cern.ch/build-prerequisites

# BUILDING THE PROJECT

Assuming $VFLOW_ROOT is the base installation directory. The CMAKE_BUILD_TYPE can be changed to Debug or Release. The following dependencies have to be installed first:

## VECCORE
* cd $VFLOW_ROOT
* git clone https://github.com/root-project/veccore.git
* mkdir $VFLOW_ROOT/veccore_build
* cd $VFLOW_ROOT/veccore_build
* cmake $VFLOW_ROOT/veccore -DCMAKE_INSTALL_PREFIX=$VFLOW_ROOT/veccore_install -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-O2 -DBUILD_VC=ON -DBUILD_UMESIMD=ON
* make -j install

**NOTE:** To build VecCore benchmarks, add to cmake the option -DBUILD_BENCHMARKS=ON

## VECMATH
* cd $VFLOW_ROOT
* git clone https://github.com/root-project/vecmath.git
* mkdir $VFLOW_ROOT/vecmath_build
* cd $VFLOW_ROOT/vecmath_build
* cmake $VFLOW_ROOT/vecmath -DCMAKE_INSTALL_PREFIX=$VFLOW_ROOT/vecmath_install -DCMAKE_BUILD_TYPE=RelWithDebInfo -DVecCore_DIR=$VFLOW_ROOT/veccore_install/lib/cmake/VecCore
* make -j install

## VECGEOM 
* cd $VFLOW_ROOT
* git clone https://gitlab.cern.ch/VecGeom/VecGeom.git
* mkdir $VFLOW_ROOT/vecgeom_build
* cd $VFLOW_ROOT/vecgeom_build
* cmake $VFLOW_ROOT/VecGeom -DCMAKE_INSTALL_PREFIX=$VFLOW_ROOT/vecgeom_install -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_PREFIX_PATH=$VFLOW_ROOT/veccore_install -DBACKEND=Vc -DROOT=ON -DROOT_DIR=$VFLOW_ROOT/root/cmake -DVECGEOM_VECTOR=avx2
* make -j install

## VECTORFLOW
* cd $VFLOW_ROOT
* mkdir $VFLOW_ROOT/vectorflow_build
* cd $VFLOW_ROOT/vectorflow_build
* cmake $VFLOW_ROOT/vectorflow -DCMAKE_INSTALL_PREFIX=$VFLOW_ROOT/vectorflow_install -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_PREFIX_PATH=/usr/local -DVecCore_DIR=$VFLOW_ROOT/veccore_install/lib/cmake/VecCore -DVc_DIR=$VFLOW_ROOT/veccore_install/lib/cmake/Vc -DVecMath_DIR=$VFLOW_ROOT/vecmath_install/lib/cmake/VecMath -DVecGeom_DIR=$VFLOW_ROOT/vecgeom_install/lib/cmake/VecGeom -DROOT_DIR=$VFLOW_ROOT/root/cmake -DROOT=ON
* make -j install

**NOTE 1:** To build VectorFlow examples and benchmarks, add to cmake the options -DBUILD_TESTING=ON &&-DBUILD_EXAMPLES=ON

**NOTE 2:** To enable GPERFTOOLS for profiling, add to cmake the option -DGPERFTOOLS=ON (see howto_gperftools.md file for instructions in how to install it)
