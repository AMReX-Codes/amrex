.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran

.. _sec:IO:

Plotfile
========

AMReX has its own native plotfile format. Many visualization tools are
available for AMReX plotfiles (see the chapter on :ref:`Chap:Visualization`).
AMReX provides the following two functions for writing a generic AMReX plotfile.
Many AMReX application codes may have their own plotfile routines that store
additional information such as compiler options, git hashes of the
source codes and :cpp:`ParmParse` runtime parameters.

.. highlight:: c++

::

      void WriteSingleLevelPlotfile (const std::string &plotfilename,
                                     const MultiFab &mf,
                                     const Vector<std::string> &varnames,
                                     const Geometry &geom,
                                     Real time,
                                     int level_step);

      void WriteMultiLevelPlotfile (const std::string &plotfilename,
                                    int nlevels,
                                    const Vector<const MultiFab*> &mf,
                                    const Vector<std::string> &varnames,
                                    const Vector<Geometry> &geom,
                                    Real time,
                                    const Vector<int> &level_steps,
                                    const Vector<IntVect> &ref_ratio);

:cpp:`WriteSingleLevelPlotfile` is for single level runs and
:cpp:`WriteMultiLevelPlotfile` is for multiple levels. The name of the
plotfile is specified by the plotfilename argument. This is the
top level directory name for the plotfile. In AMReX convention, the
plotfile name consist of letters followed by numbers (e.g.,
plt00258). :cpp:`amrex::Concatenate` is a useful helper function for
making such strings.

.. highlight:: c++

::

      int istep = 258;
      const std::string& pfname = amrex::Concatenate("plt",istep); // plt00258

      // By default there are 5 digits, but we can change it to say 4.
      const std::string& pfname2 = amrex::Concatenate("plt",istep,4); // plt0258

      istep =1234567;  // Having more than 5 digits is OK.
      const std::string& pfname3 = amrex::Concatenate("plt",istep); // plt1234567

The argument :cpp:`mf` above (:cpp:`MultiFab` for single level and
:cpp:`Vector<const MultiFab*>` for multi-level) is the data to be written
to the disk. Note that many visualization tools expect this to be
cell-centered data. So for nodal data, we need to convert them to
cell-centered data through some kind of averaging. Also note that if
you have data at each AMR level in several MultiFabs, you need
to build a new MultiFab at each level to hold all the data on
that level. This involves local data copy in memory and is not
expected to significantly increase the total wall time for writing
plotfiles. For the multi-level version, the function expects
:cpp:`Vector<const MultiFab*>`, whereas the multi-level data are often
stored as :cpp:`Vector<std::unique_ptr<MultiFab>>`. AMReX has a
helper function for this and one can use it as follows,

.. highlight:: c++

::

       WriteMultiLevelPlotfile(......, amrex::GetVecOfConstPtrs(mf), ......);

The argument :cpp:`varnames` has the names for each component of the
MultiFab data. The size of the Vector should be equal to the
number of components. The argument :cpp:`geom` is for passing
:cpp:`Geometry` objects that contain the physical domain
information. The argument :cpp:`time` is for the time associated with the
data. The argument :cpp:`level_step` is for the current time step
associated with the data. For multi-level plotfiles, the argument
:cpp:`nlevels` is the total number of levels, and we also need to provide
the refinement ratio via an :cpp:`Vector` of size nlevels-1.

We note that AMReX does not overwrite old plotfiles if the new
plotfile has the same name. The old plotfiles will be renamed to
new directories named like plt00350.old.46576787980.

Async Output
============

AMReX provides the ability to print MultiFabs, plotfiles and
particle data asynchronously.  Asynchronous output works by creating
a copy of the data at the time of the call, which is written to disk
by a persistent thread created during AMReX's initialization.  This allows
the calculation to continue immediately, which can drastically reduce
walltime spent writing to disk.

If the number of output files is less than the number of MPI ranks,
AMReX's async output requires MPI to be initialized with THREAD_MULTIPLE
support. THREAD_MULTIPLE support allows multiple unique threads to run unique
MPI calls simultaneously.  This support is required to allow AMReX applications
to perform MPI work while the Async Output concurrently pings ranks to signal
that they can safely begin writing to their assigned files.  However,
THREAD_MULTIPLE can introduce additional overhead as each threads' MPI operations
must be scheduled safely around each other. Therefore, AMReX uses a lower level
of support, SERIALIZED, by default and applications have to turn on THREAD_MULTIPLE
support.

To turn on Async Output, use the input flag ``amrex.async_out=1``.  The number
of output files can also be set, using ``amrex.async_out_nfiles``.  The default
number of files is ``64``. If the number of ranks is larger than the number of
files, THREAD_MULTIPLE must be turned on by adding
``MPI_THREAD_MULTIPLE=TRUE`` to the GNUMakefile. Otherwise, AMReX
will throw an error.

Async Output works for a wide range of AMReX calls, including:

* ``amrex::WriteSingleLevelPlotfile()``
* ``amrex::WriteMultiLevelPlotfile()``
* ``amrex::WriteMLMF()``
* ``VisMF::AsyncWrite()``
* ``ParticleContainer::Checkpoint()``
* ``ParticleContainer::WritePlotFile()``
* ``Amr::writePlotFile()``
* ``Amr::writeSmallPlotFile()``
* ``Amr::checkpoint()``
* ``AmrLevel::writePlotFile()``
* ``StateData::checkPoint()``
* ``FabSet::write()``

Be aware: when using Async Output, a thread is spawned and exclusively used
to perform output throughout the runtime.  As such, you may oversubscribe
resources if you launch an AMReX application that assigns all available
hardware threads in another way, such as OpenMP.  If you see any degradation
when using Async Output and OpenMP, try using one less thread in
``OMP_NUM_THREADS`` to prevent oversubscription and get more consistent
results.

Checkpoint File
===============

Checkpoint files are used for restarting simulations from where the
checkpoints are written. Each application code has its own set of
data needed for restart. AMReX provides I/O functions for basic
data structures like :cpp:`MultiFab` and :cpp:`BoxArray`. These
functions can be used to build codes for reading and writing
checkpoint files. Since each application code has its own
requirement, there is no standard AMReX checkpoint format.
However we have provided an example restart capability in the tutorial
`Advection AmrCore`_.
Refer to the functions :cpp:`ReadCheckpointFile()` and
:cpp:`WriteCheckpointFile()` in this tutorial.

.. _`Advection AmrCore`: https://amrex-codes.github.io/amrex/tutorials_html/AMR_Tutorial.html#advection-amrcore

A checkpoint file is actually a directory with name, e.g.,
``chk00010`` containing a ``Header`` (text) file, along with
subdirectories ``Level_0``, ``Level_1``, etc. containing the
:cpp:`MultiFab` data at each level of refinement.
The ``Header`` file contains problem-specific data (such as the
finest level, simulation time, time step, etc.), along with a printout
of the :cpp:`BoxArray` at each level of refinement.

When starting a simulation from a checkpoint file, a typical sequence in the code
could be:

- Read in the ``Header`` file data (except for the :cpp:`BoxArray` data).

- For each level of refinement, do the following in order:

  -- Read in the :cpp:`BoxArray`

  -- Build a :cpp:`DistributionMapping`

  -- Define any :cpp:`MultiFab`, :cpp:`FluxRegister`, etc. objects that are built upon the
  :cpp:`BoxArray` and the :cpp:`DistributionMapping`

  -- Read in the :cpp:`MultiFab` data

We do this one level at a time because when you create a distribution map,
it checks how much allocated :cpp:`MultiFab` data already exists before assigning
grids to processors.

Typically a checkpoint file is a directory containing some text files
and sub-directories (e.g., ``Level_0`` and ``Level_1``)
containing various data. It is a good idea that we fist make these
directories ready for subsequently writing to the disk. For example,
to build directories ``chk00010``, ``chk00010/Level_0``, and
``chk00010/Level_1``, you could write:

.. highlight:: c++

::

   const std::string& checkpointname = amrex::Concatenate("chk",10);

   amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

   const int nlevels = 2;

   bool callBarrier = true;

   // ---- prebuild a hierarchy of directories
   // ---- dirName is built first.  if dirName exists, it is renamed.  then build
   // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
   // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
   // ---- after all directories are built
   // ---- ParallelDescriptor::IOProcessor() creates the directories
   amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, callBarrier);

A checkpoint file of AMReX application codes often has a clear text
Header file that only the I/O process writes to it using
:cpp:`std::ofstream`. The Header file contains problem-dependent
information such as
the time, the physical domain size, grids, etc. that are necessary for
restarting the simulation. To guarantee that precision is not lost
for storing floating point number like time in clear text file, the
file stream's precision needs to be set properly. And a stream buffer
can also be used. For example,

.. highlight:: c++

::

   // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       std::string HeaderFileName(checkpointname + "/Header");
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                               std::ofstream::trunc |
                                               std::ofstream::binary);

       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for AmrCoreAdv\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }


:cpp:`amrex::VisMF` is a class that can be used to perform
:cpp:`MultiFab` I/O in parallel. How many processes are allowed to
perform I/O simultaneously can be set via

::

      VisMF::SetNOutFiles(64);  // up to 64 processes, which is also the default.

The optimal number is of course system dependent. The following code
shows how to write a :cpp:`MultiFab`.

.. highlight:: c++

::

   // write the MultiFab data to, e.g., chk00010/Level_0/
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
   }

It should also be noted that all the
data including those in ghost cells are written/read by
:cpp:`VisMF::Write/Read`.

For reading the Header file, AMReX can have the I/O process
read the file from the disk and broadcast it to others as
:cpp:`Vector<char>`. Then all processes can read the information with
:cpp:`std::istringstream`. For example,

.. highlight:: c++

::

    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

The following code how to read in a :cpp:`BoxArray`, create a
:cpp:`DistributionMapping`, build :cpp:`MultiFab` and :cpp:`FluxRegister` data,
and read in a :cpp:`MultiFab` from a checkpoint file, on a level-by-level basis:

.. highlight:: c++

::

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp = 1;
        int nghost = 0;
        phi_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        phi_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        if (lev > 0 && do_reflux) {
            flux_reg[lev] = std::make_unique<FluxRegister>(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp);
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

It should be emphasized that calling :cpp:`VisMF::Read` with an empty
:cpp:`MultiFab` (i.e., no memory allocated for floating point data)
will result in a :cpp:`MultiFab` with a new :cpp:`DistributionMapping`
that could be different from any other existing
:cpp:`DistributionMapping` objects and is not recommended.
