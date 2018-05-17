# This Dockerfile is used for automated testing of WarpX

FROM ubuntu:14.04

# Install a few packages, as root
RUN apt-get update \
    && apt-get install -y \
    wget \
    make \
    git \
    gcc \
    gfortran \
    g++ \
    python

WORKDIR /home/

# Install MPI
RUN apt-get install -y \
    openmpi-bin \
libopenmpi-dev

# Install FFTW
RUN apt-get install -y \
    libfftw3-dev \
    libfftw3-mpi-dev

# Clone amrex
RUN git clone https://github.com/AMReX-Codes/amrex.git \
    && cd amrex/ \
    && git checkout development \
    && cd ..

# Clone their regression test utility
RUN git clone https://github.com/AMReX-Codes/regression_testing.git

# Clone picsar
RUN git clone https://bitbucket.org/berkeleylab/picsar.git

# Copy warpx
RUN mkdir -p /home/warpx/
COPY ./ /home/warpx/

# Prepare regression tests
RUN mkdir -p rt-WarpX \
    && cp warpx/Regression/WarpX-tests.ini rt-WarpX \
    && sed -i 's\regtester/\\g' rt-WarpX/WarpX-tests.ini \
    && sed -i 's/sendEmailWhenFail = 1/sendEmailWhenFail = 0/g' rt-WarpX/WarpX-tests.ini \
    && sed -i 's\AMReX_RegTesting/\\g' rt-WarpX/WarpX-tests.ini

