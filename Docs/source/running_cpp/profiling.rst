Profiling the code
==================

Profiling with AMREX's built-in profiling tools
-----------------------------------------------
See `this page <https://amrex-codes.github.io/amrex/docs_html/Chapter12.html>`__ in the AMReX documentation.


Profiling the code with Intel Advisor on NERSC
----------------------------------------------

Follow these steps:

- Instrument the code during compilation

  ::

     module swap craype-haswell craype-mic-knl
     make -j 16 COMP=intel USE_VTUNE=TRUE

  (where the first line is only needed for KNL)

- In your SLURM submission script, use the following
  lines in order to run the executable. (In addition
  to setting the usual ``OMP`` environment variables.)

  ::

     module load advisor
     export ADVIXE_EXPERIMENTAL=roofline
     srun -n <n_mpi> -c <n_logical_cores_per_mpi> --cpu_bind=cores advixe-cl -collect survey -project-dir advisor -trace-mpi -- <warpx_executable> inputs
     srun -n <n_mpi> -c <n_logical_cores_per_mpi> --cpu_bind=cores advixe-cl -collect tripcounts -flop -project-dir advisor -trace-mpi -- <warpx_executable> inputs

  where ``<n_mpi>`` and ``<n_logical_cores_per_mpi>`` should be replaced by
  the proper values, and ``<warpx_executable>`` should be replaced by the
  name of the WarpX executable.

- Launch the Intel Advisor GUI

  ::

     module load advisor
     advixe-gui

  (Note: this requires to use ``ssh -XY`` when connecting to Cori.)

