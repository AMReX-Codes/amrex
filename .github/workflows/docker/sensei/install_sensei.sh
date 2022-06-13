#!/bin/bash

set -x
set -e

DOWNLOAD_ROOT=/tmp/downloads
EXTRACT_ROOT=/tmp/downloads/extracted
mkdir -p ${DOWNLOAD_ROOT}
mkdir -p ${EXTRACT_ROOT}

# Download and install Sensei (VTK + MPI + Python)
export CMAKE_GENERATOR=Ninja

cd ${DOWNLOAD_ROOT}
mkdir sensei
cd sensei
git clone https://github.com/SENSEI-insitu/sensei.git
cd sensei
git checkout develop
mkdir build
cd build
cmake .. \
  -D SENSEI_VERSION=${SENSEI_VERSION} \
  -D BUILD_TESTING=OFF \
  -D ENABLE_PYTHON=ON \
  -D ENABLE_MANDELBROT=OFF \
  -D ENABLE_OSCILLATORS=OFF \
  -D ENABLE_VTK_IO=ON \
  -D ENABLE_VTK_MPI=ON \
  -D MPI4PY_INCLUDE_DIR=/usr/lib64/python3.10/site-packages/openmpi/mpi4py/include
cmake --build .
cmake --install . --prefix ${SENSEI_DIR}

rm -rf ${DOWNLOAD_ROOT}
rm -rf ${EXTRACT_ROOT}
