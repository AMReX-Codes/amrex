#!/bin/bash

# Include common build tools
. $(dirname ${BASH_SOURCE[0]})/tools.sh

DOWNLOAD_ROOT=/root/downloads
EXTRACT_ROOT=/root/downloads/extracted
mkdir -p ${DOWNLOAD_ROOT}
mkdir -p ${EXTRACT_ROOT}

# Download and install VTK (Common + Python bindings)
export CMAKE_GENERATOR=Ninja
export CC=clang
export CXX=clang++

parse_version_env VTK_VERSION
curl https://vtk.org/files/release/${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}/VTK-${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}.${VTK_VERSION_PATCH}.tar.gz \
  --output ${DOWNLOAD_ROOT}/vtk.tar.gz \
  --silent
extract ${DOWNLOAD_ROOT}/vtk.tar.gz -C ${EXTRACT_ROOT}/vtk

cd ${EXTRACT_ROOT}/vtk/VTK-${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}.${VTK_VERSION_PATCH} &&
mkdir -p build &&
cd build &&
cmake .. \
 -DCMAKE_BUILD_TYPE=Release \
 -DVTK_BUILD_EXAMPLES=OFF \
 -DVTK_BUILD_TESTING=OFF \
 -DVTK_BUILD_DOCUMENTATION=OFF \
 -DVTK_BUILD_SHARED_LIBS=ON \
 -DVTK_USE_MPI=ON \
 -DVTK_GROUP_ENABLE_StandAlone=NO \
 -DVTK_GROUP_ENABLE_Imaging=NO \
 -DVTK_GROUP_ENABLE_Qt=NO \
 -DVTK_GROUP_ENABLE_Rendering=NO \
 -DVTK_GROUP_ENABLE_Views=NO \
 -DVTK_GROUP_ENABLE_Web=NO \
 -DVTK_WRAP_PYTHON=YES \
 -DVTK_PYTHON_VERSION=3 \
 -DVTK_MODULE_ENABLE_VTK_PythonInterpreter=YES \
 -DVTK_MODULE_ENABLE_VTK_WrappingPythonCore=YES \
 -DVTK_MODULE_ENABLE_VTK_CommonDataModel=YES \
 -DVTK_MODULE_ENABLE_VTK_CommonMisc=YES \
 -DVTK_MODULE_ENABLE_VTK_CommonExecutionModel=YES \
 -DVTK_MODULE_ENABLE_VTK_IOParallelNetCDF=NO \
 -DVTK_MODULE_ENABLE_VTK_IOMPIParallel=NO \
 -DVTK_MODULE_ENABLE_VTK_IOMPIImage=NO \
 -DVTK_MODULE_ENABLE_VTK_FiltersParallelVerdict=NO \
 -DVTK_MODULE_ENABLE_VTK_FiltersParallelGeometry=NO \
 -DVTK_MODULE_ENABLE_VTK_DomainsParallelChemistry=NO \
 -DVTK_MODULE_ENABLE_VTK_FiltersParallelMPI=NO \
 -DVTK_MODULE_ENABLE_VTK_ParallelCore=YES \
 -DVTK_MODULE_ENABLE_VTK_CommonSystem=YES \
 -DVTK_MODULE_ENABLE_VTK_IOLegacy=YES \
 -DVTK_MODULE_ENABLE_VTK_IOCore=YES \
 -DVTK_MODULE_ENABLE_VTK_IOXMLParser=YES \
 -DVTK_MODULE_ENABLE_VTK_IOXML=YES \
 -DVTK_MODULE_ENABLE_VTK_CommonMath=YES \
 -DVTK_MODULE_ENABLE_VTK_CommonTransforms=YES &&
cmake --build . &&
cmake --install . --prefix ${VTK_DIR}

rm -rf ${DOWNLOAD_ROOT}
rm -rf ${EXTRACT_ROOT}
