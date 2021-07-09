#!/bin/bash

DOWNLOAD_ROOT=/root/downloads
EXTRACT_ROOT=/root/downloads/extracted
mkdir -p ${DOWNLOAD_ROOT}
mkdir -p ${EXTRACT_ROOT}

# Download and install Sensei (VTK + MPI + Python)
export CMAKE_GENERATOR=Ninja

cd ${DOWNLOAD_ROOT}
mkdir sensei
cd sensei
git clone https://gitlab.kitware.com/sensei/sensei.git &&
cd sensei &&
git checkout ${SENSEI_VERSION} &&
mkdir build &&
cd build &&
cmake .. \
  -D BUILD_TESTING=OFF \
  -D ENABLE_PYTHON=ON \
  -D ENABLE_MANDELBROT=OFF \
  -D ENABLE_OSCILLATORS=OFF \
  -D ENABLE_VTK_MPI=ON \
  -D MPI4PY_INCLUDE_DIR=/usr/lib64/python3.9/site-packages/openmpi/mpi4py/include &&
cmake --build . &&
cmake --install . --prefix ${SENSEI_DIR}

rm -rf ${DOWNLOAD_ROOT}
rm -rf ${EXTRACT_ROOT}
