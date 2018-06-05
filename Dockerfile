# This Dockerfile is used for automated testing of WarpX on Shippable
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
    ccache \
    openmpi-bin \
    libopenmpi-dev \
    libfftw3-dev \
    libfftw3-mpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Create new user
RUN useradd --create-home regtester
RUN mkdir -p /home/regtester/AMReX_RegTesting/warpx/
COPY ./ /home/regtester/AMReX_RegTesting/warpx/
RUN chown -R regtester /home/regtester/
# Grant sudo access without password
RUN echo 'regtester ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers

USER regtester
WORKDIR /home/regtester/AMReX_RegTesting/

# Install miniconda
RUN cd /home/regtester \
    && wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh \
    && bash Miniconda-latest-Linux-x86_64.sh -b
ENV PATH /home/regtester/miniconda2/bin:$PATH

# Install required python packages
RUN pip install --upgrade pip \
    && pip install numpy scipy matplotlib yt

# Clone amrex
RUN git clone https://github.com/AMReX-Codes/amrex.git \
    && cd amrex/ \
    && git checkout development \
    && cd ..

# Clone their regression test utility
RUN git clone https://github.com/AMReX-Codes/regression_testing.git

# Clone picsar
RUN git clone https://bitbucket.org/berkeleylab/picsar.git

# Prepare regression tests
RUN mkdir -p rt-WarpX/WarpX-benchmarks \
    && cd warpx/Regression \
    && python prepare_file_shippable.py \
    && cp shippable-tests.ini ../../rt-WarpX
