#!/bin/bash

BoxLib_CONFIG=0
BoxLib_Root=${PWD}
BoxLib_INSTALL_PREFIX=${BoxLib_Root}/install
BoxLib_MAKE=0
BoxLib_CLOBBER=0
BoxLib_TEST=0
BoxLib_INSTALL=0
BoxLib_VERBOSE=0

ENABLE_MPI=1
ENABLE_OpenMP=1
MPI_PREFIX=/usr

# 1, 2, 3
SPACEDIM=2

# FLOAT, DOUBLE
PRECISION=DOUBLE

#debug, release
CMAKE_BUILD_TYPE=debug

# cmake, ctest executables
CMAKE=cmake
CTEST=ctest

# Relative to BoxLib_Root
Build_Dir=build

# cores for build
BoxLib_NPROCS=1

while getopts "ab:cfin:mtv" flag
do
  case $flag in
    a) BoxLib_CONFIG=1; BoxLib_MAKE=1; BoxLib_TEST=1; BoxLib_INSTALL=1;;
    b) Build_Dir=${OPTARG};;
    c) BoxLib_CONFIG=1;;
    f) BoxLib_CLOBBER=1;;
    i) BoxLib_INSTALL=1;;
    n) BoxLib_NPROCS=${OPTARG};;
    m) BoxLib_MAKE=1;;
    t) BoxLib_TEST=1;;
    v) BoxLib_VERBOSE=1;;
  esac
done

BoxLib_Build=${BoxLib_Root}/${Build_Dir}

echo BoxLib_Root=$BoxLib_Root
echo BoxLib_Build=$BoxLib_Build
echo BoxLib_CONFIG=$BoxLib_CONFIG
echo BoxLib_MAKE=$BoxLib_MAKE
echo BoxLib_TEST=$BoxLib_TEST
echo BoxLib_INSTALL=$BoxLib_INSTALL
echo BoxLib_NPROCS=$BoxLib_NPROCS

if [ $BoxLib_CLOBBER -eq 1 ]; then
    rm -rf ${BoxLib_Build}
fi

if [ $BoxLib_CONFIG -eq 1 ]; then
    mkdir -p ${BoxLib_Build}
    cd ${BoxLib_Build}

    CMAKE_ARGS="-D CMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE} \
                -D ENABLE_Config_Report:BOOL=ON \
                -D BL_SPACEDIM:INT=${SPACEDIM} \
                -D BL_PRECISION:STRING=${PRECISION} \
                -D ENABLE_TESTS:BOOL=ON \
                -D CMAKE_INSTALL_PREFIX:FILEPATH=${BoxLib_INSTALL_PREFIX} "

    if [ ${BoxLib_VERBOSE} -eq 1 ]; then
        CMAKE_ARGS=${CMAKE_ARGS} " --debug-output"
    fi

    if [ ${ENABLE_MPI} -eq 1 ]; then
        CMAKE_ARGS="${CMAKE_ARGS} -D ENABLE_MPI:BOOL=${ENABLE_MPI} \
                                   -D MPI_PREFIX:FILEPATH=${MPI_PREFIX}"

        if [ ${BoxLib_TEST} -eq 1 ]; then
            CMAKE_ARGS="${CMAKE_ARGS} -D MPI_EXEC:FILEPATH=${MPI_PREFIX}/bin/mpiexec \
                                      -D MPI_EXEC_NUMPROCS_FLAG:STRING=-np \
                                      -D MPI_EXEC_ARGS:STRING='-mca mpi_yield_when_idle 1'"
        fi
    fi

    if [ ${ENABLE_OpenMP} -eq 1 ]; then
        CMAKE_ARGS="${CMAKE_ARGS} -D ENABLE_OpenMP:BOOL=${ENABLE_OpenMP}"
    fi

    ${CMAKE} ${CMAKE_ARGS} ${BoxLib_Root}

    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

if [ $BoxLib_MAKE -eq 1 ]; then

    cd ${BoxLib_Build}

    MAKE_ARGS=

    if [ ${BoxLib_NPROCS} -ne 1 ]; then
        MAKE_ARGS="${MAKE_ARGS} -j ${BoxLib_NPROCS}"
    fi

    if [ ${BoxLib_VERBOSE} -eq 1 ]; then
        MAKE_ARGS="${MAKE_ARGS} VERBOSE=ON"
    fi

    make ${MAKE_ARGS}

    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

if [ $BoxLib_INSTALL -eq 1 ]; then
    cd ${BoxLib_Build}

    make install

    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

if [ $BoxLib_TEST -eq 1 ]; then
    cd ${BoxLib_Build}

    ${CTEST} --timeout 60 --output-on-failure

    if [ $? -ne 0 ]; then
        exit 1
    fi
fi

