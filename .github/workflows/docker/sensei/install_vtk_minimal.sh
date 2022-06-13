#!/bin/bash

set -x
set -e

# Include common build tools
. $(dirname ${BASH_SOURCE[0]})/tools.sh

DOWNLOAD_ROOT=/tmp/downloads
EXTRACT_ROOT=/tmp/downloads/extracted
mkdir -p ${DOWNLOAD_ROOT}
mkdir -p ${EXTRACT_ROOT}

# Download and install VTK (Common + Python bindings)
export CMAKE_GENERATOR=Ninja
export CC=clang
export CXX=clang++

if [[ -z ${VTK_VERSION} ]]; then
  VTK_VERSION=9.1.0
fi

parse_version_env VTK_VERSION
curl https://vtk.org/files/release/${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}/VTK-${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}.${VTK_VERSION_PATCH}.tar.gz \
  --output ${DOWNLOAD_ROOT}/vtk.tar.gz \
  --silent
extract ${DOWNLOAD_ROOT}/vtk.tar.gz -C ${EXTRACT_ROOT}/vtk

cd ${EXTRACT_ROOT}/vtk/VTK-${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}.${VTK_VERSION_PATCH}
patch -i /tmp/vtk_use_mpi.patch
mkdir -p build
cd build
cmake .. \
 -DCMAKE_BUILD_TYPE=Release \
 -DVTK_BUILD_EXAMPLES=OFF \
 -DVTK_BUILD_TESTING=OFF \
 -DVTK_BUILD_DOCUMENTATION=OFF \
 -DVTK_BUILD_SHARED_LIBS=ON \
 -DVTK_GROUP_ENABLE_StandAlone=DONT_WANT \
 -DVTK_GROUP_ENABLE_Imaging=NO \
 -DVTK_GROUP_ENABLE_Qt=NO \
 -DVTK_GROUP_ENABLE_Rendering=NO \
 -DVTK_GROUP_ENABLE_Views=NO \
 -DVTK_GROUP_ENABLE_Web=NO \
 -DVTK_GROUP_ENABLE_MPI=DONT_WANT \
 -DVTK_USE_MPI=ON \
 -DVTK_MODULE_ENABLE_VTK_ParallelMPI=YES \
 -DVTK_MODULE_ENABLE_VTK_IOParallelXML=YES \
 -DVTK_MODULE_ENABLE_VTK_IOXML=YES \
 -DVTK_MODULE_ENABLE_VTK_IOLegacy=YES \
 -DVTK_ENABLE_KITS=OFF \
 -DVTK_WRAP_PYTHON=YES \
 -DVTK_PYTHON_VERSION=3
cmake --build .
cmake --install . --prefix ${VTK_DIR}

rm -rf ${DOWNLOAD_ROOT}
rm -rf ${EXTRACT_ROOT}
