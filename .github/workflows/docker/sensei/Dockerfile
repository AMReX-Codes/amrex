FROM fedora:35

# Set install root
ENV PACKAGE_ROOT=/root/install

COPY tools.sh /root/bin/tools.sh

# Copy and run the install script
COPY install_deps.sh /root/bin/install_deps.sh
RUN /root/bin/install_deps.sh

# Configure MPI environment
ENV MPI_HOME=/usr/lib64/openmpi/

# Configure Python environment
ENV PYTHONPATH=/usr/lib64/python3.10/site-packages/openmpi

# Configure VTK environment
ENV VTK_VERSION=9.1.0
ENV VTK_DIR=${PACKAGE_ROOT}/vtk/${VTK_VERSION}
COPY install_vtk_minimal.sh /root/bin/install_vtk.sh
COPY vtk_use_mpi.patch /tmp/vtk_use_mpi.patch
RUN /root/bin/install_vtk.sh

# Configure Sensei Environment
ENV SENSEI_VERSION=v4.0.0
ENV SENSEI_DIR=${PACKAGE_ROOT}/sensei/${SENSEI_VERSION}
COPY install_sensei.sh /root/bin/install_sensei.sh
RUN /root/bin/install_sensei.sh
