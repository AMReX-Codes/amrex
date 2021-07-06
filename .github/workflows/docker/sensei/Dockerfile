FROM fedora:33

# Set install root
ENV PACKAGE_ROOT=/root/install

# Configure MPI environment
ENV MPI_HOME=/usr/lib64/openmpi/

# Configure Python environment
ENV PYTHONPATH=/usr/lib64/python3.9/site-packages/openmpi

# Configure VTK environtment
ENV VTK_VERSION=9.0.1
ENV VTK_DIR=${PACKAGE_ROOT}/vtk/${VTK_VERSION}

# Configure Sensei Environment
ENV SENSEI_VERSION=develop
ENV SENSEI_DIR=${PACKAGE_ROOT}/sensei/${SENSEI_VERSION}

# Copy and run the install script
COPY install_deps.sh /root/bin/install_deps.sh
RUN /root/bin/install_deps.sh

COPY tools.sh /root/bin/tools.sh

COPY install_vtk_minimal.sh /root/bin/install_vtk.sh
RUN /root/bin/install_vtk.sh

COPY install_sensei.sh /root/bin/install_sensei.sh
RUN /root/bin/install_sensei.sh
