
plotfiles and checkpoints are directories.
each mpi process writes its own data.
each file is written to by only one process at a time.
c++ binary i/o.
data format independent of byte order, precision, nprocs.
nfiles [1,nprocs], user settable.

resiliency--stream retry
nooverwrite of existing data
.temp files

demand driven reads.
headers contain min/max and seek for each grid
data is addressable to a single component of a single grid.
no restriction on the relationship between nprocs and nfiles for reading.
stream throttling for reading to prevent thrashing.


compressed metadata
preallocation
dynamic set selection
striping tests
more combining writes
tests for meta operations:  seek open close mkdir rename
//make set selection a function
//recursive mkdir
dot in path issue
remove how
fix VisMF::Check for new formats.
remove duplicate writeplotfile code.
add more retry code to vismf
test set_ghost
check nfilesiter with long paths
special profiling
support new shared boxarrays (test restart)
check vismf reads for copies and non-contiguous
check vismf colors
fix types on stream functions (streamsize in read, streampos, etc.)
add NFiles files to CMakeLists.txt.
tests for copy multifab speed.
new fabarray read algorithm.
tests for buffer sizes (copy buffer [shouldread], pubsetbuf, no buffer).
check if pubsetbuf works
buffer size adjustments
check visit and yt with new formats.
test performance of rank order vs. file order reads.
byte order swap for integers.
persistent streams
nstreams in new read
pre make plot and check directories.  ParticleContainer::Checkpoint
  also makes directories.
tutorial for nfiles.
