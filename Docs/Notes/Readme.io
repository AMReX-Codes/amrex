
plotfiles and checkpoints are directories.
each mpi process writes its own data.
each file is written to by only one process at a time.
c++ binary i/o.
data format independent of byte order, precision, nprocs.
nfiles [1,nprocs], user settable.

resiliency--stream retry
nooverwrite of existing data, existing directories are renamed.
directories are named *.temp until complete.

demand driven reads.
headers contain min/max and seek for each grid (VisMF Header Versions 1 and 3).
data is addressable to a single component of a single grid.
no restriction on the relationship between nprocs and nfiles for reading.
stream throttling for reading to prevent thrashing.



// ---- i/o parameters
vismf.headerversion           (def:  Version_v1  (1) )
vismf.groupsets               (def:  false)
vismf.setbuf                  (def:  true)
vismf.usesingleread           (def:  false)
vismf.usesinglewrite          (def:  false)
vismf.checkfilepositions      (def:  false)
vismf.usepersistentifstreams  (def:  true)
vismf.usesynchronousreads     (def:  false)
vismf.usedynamicsetselection  (def:  true)
vismf.iobuffersize            (def:  VisMF::IO_Buffer_Size)
amr.plot_nfiles               (def:  64)
amr.checkpoint_nfiles         (def:  64)
amr.mffile_nstreams           (def:  1)
amr.plot_headerversion        (def:  Version_v1  (1) )
amr.checkpoint_headerversion  (def:  Version_v1  (1) )
amr.prereadFAHeaders          (def:  true)
amr.precreateDirectories      (def:  true)

particles.particles_nfiles = 1024

you can also call these to set fabconv buffer sizes:
RealDescriptor::SetReadBufferSize(rbs);
RealDescriptor::SetWriteBufferSize(wbs);


// ---- recommended settings
amr.checkpoint_headerversion = 2
vismf.usesingleread          = true
vismf.usesinglewrite         = true
amr.mffile_nstreams          = 4


