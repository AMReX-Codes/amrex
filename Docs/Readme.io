
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
persistent streams
dynamic set selection
buffer size adjustments
striping.
more combining writes.
tests for meta operations:  seek open close mkdir rename
make set selection a function.
recursive mkdir
dot in path issue
remove how
check if pubsetbuf works
add more retry code to vismf
