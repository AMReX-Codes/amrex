# A build environment with AMReX compiled and installed, for reproducing user
# reports and for working through the tutorials without setting up a toolchain.
#
#   docker build -t amrex .
#   docker run --rm -it -v "$PWD":/work amrex
#
# The build is configured through --build-arg, e.g. a 2D build without MPI:
#
#   docker build -t amrex --build-arg AMREX_SPACEDIM=2 --build-arg AMREX_MPI=OFF .
#
# AMReX is installed into /usr/local, so a downstream project can pick it up
# with find_package(AMReX REQUIRED). The compilers and CMake are kept in the
# image so applications can be built against it.

# Ubuntu 24.04 provides CMake 3.28 and GCC 13, meeting AMReX's CMake 3.25 and
# C++20 requirements. Pinned by digest so the image is reproducible.
FROM ubuntu:24.04@sha256:008173c23f95b170204355c12626cb5a965d779a7e1283b09e9cffbb1bf33ca3

ARG AMREX_SPACEDIM=3
ARG AMREX_MPI=ON
ARG AMREX_OMP=OFF
ARG AMREX_PRECISION=DOUBLE
ARG BUILD_JOBS=4

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        git \
        libopenmpi-dev \
        openmpi-bin \
        python3 \
    && rm -rf /var/lib/apt/lists/*

COPY . /opt/amrex-src

RUN cmake -S /opt/amrex-src -B /tmp/amrex-build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DAMReX_SPACEDIM="${AMREX_SPACEDIM}" \
        -DAMReX_MPI="${AMREX_MPI}" \
        -DAMReX_OMP="${AMREX_OMP}" \
        -DAMReX_PRECISION="${AMREX_PRECISION}" \
    && cmake --build /tmp/amrex-build -j "${BUILD_JOBS}" --target install \
    && rm -rf /tmp/amrex-build

# Running as root inside the container is avoided; /work is the mount point for
# the caller's own sources.
RUN useradd --create-home --shell /bin/bash amrex \
    && mkdir -p /work \
    && chown amrex:amrex /work
USER amrex
WORKDIR /work

CMD ["/bin/bash"]
