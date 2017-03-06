#!/usr/bin/env bash

BOXLIB_HEADERS="\
ABec2_F\.H \
ABec2\.H \
ABec4_F\.H \
ABec4\.H \
ABec_F\.H \
ABecLaplacian\.H \
AmrCore\.H \
AmrData\.H \
Amr\.H \
AmrLevel\.H \
AmrParGDB\.H \
AmrParticles\.H \
AmrvisConstants\.H \
Arena\.H \
Array\.H \
ArrayLim\.H \
AuxBoundaryData\.H \
BArena\.H \
BaseFab_f\.H \
BaseFab\.H \
BCRec\.H \
BC_TYPES\.H \
BLassert\.H \
BLBackTrace\.H \
BLFort\.H \
BLPgas\.H \
BLProfiler\.H \
BndryData\.H \
BndryRegister\.H \
BoundCond\.H \
BoxArray\.H \
BoxDomain\.H \
Box\.H \
BoxLib\.H \
BoxList\.H \
CArena\.H \
ccse-mpi\.H \
CGSolver\.H \
Cluster\.H \
CONSTANTS\.H \
COORDSYS_F\.H \
CoordSys\.H \
DataServices\.H \
Derive\.H \
DistributionMapping\.H \
DivVis_F\.H \
DivVis\.H \
ErrorList\.H \
Extrapolater\.H \
FabArray\.H \
FabConv\.H \
FabSet\.H \
FArrayBox\.H \
FillPatchUtil\.H \
FLUSH_F\.H \
FLUXREG_F\.H \
FluxRegister\.H \
FMultiGrid\.H \
FPC\.H \
Geometry\.H \
IArrayBox\.H \
iMultiFab\.H \
IndexType\.H \
INTERPBNDRYDATA_F\.H \
InterpBndryData\.H \
INTERP_F\.H \
Interpolater\.H \
IntVect\.H \
Laplacian\.H \
Lazy\.H \
LevelBld\.H \
LinOp\.H \
LO_BCTYPES\.H \
LO_F\.H \
Looping\.H \
LP_F\.H \
MacBndry\.H \
MAKESLICE_F\.H \
Mask\.H \
MCCGSolver\.H \
MCINTERPBNDRYDATA_F\.H \
MCInterpBndryData\.H \
MCLinOp\.H \
MCLO_F\.H \
MCMultiGrid\.H \
MemPool\.H \
MemProfiler\.H \
MG_F\.H \
MGT_Solver\.H \
MultiFab\.H \
MultiFabUtil_F\.H \
MultiFabUtil\.H \
MultiGrid\.H \
MultiMask\.H \
NFiles\.H \
Orientation\.H \
ParallelDescriptor\.H \
ParGDB\.H \
ParmParse\.H \
ParticleInit\.H \
Particles_F\.H \
Particles\.H \
Periodicity\.H \
PhysBCFunct\.H \
PList\.H \
PlotFileUtil\.H \
Pointers\.H \
PROB_AMR_F\.H \
RealBox\.H \
REAL\.H \
SLABSTAT_F\.H \
SlabStat\.H \
SPACE_F\.H \
SPACE\.H \
StateData\.H \
StateDescriptor\.H \
StationData\.H \
stencil_types\.H \
TagBox\.H \
TinyProfiler\.H \
TracerParticles\.H \
Tuple\.H \
UseCount\.H \
Utility\.H \
VisMF\.H \
winstd\.H"

NUM_HEADERS=126
I_HEADER=0

# include 'bc_types.fi' in fortran
old="bc_types\.fi"
find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.f" -o -name "*.f90" -o -name "*.F" -o -name "*.F90" \) -exec grep -Iq . {} \; -exec sed -i 's/\(include\s*\(\x27\|\"\)\s*\)'"${old}"'\(\s*\(\x27\|\"\)\)/\1AMReX_'"${old}"'\3/g' {} +

for header in ${BOXLIB_HEADERS}; do
  find . -type d \( -name .git -o -path ./Tools/Migration -o -path ./Docs/Migration \) -prune -o -type f \( -name "*.f" -o -name "*.f90" -o -name "*.F" -o -name "*.F90" -o -name "*.H" -o -name "*.h" -o -name "*.cpp" -o -name "*.c" -o -name "*.cc" \) -exec grep -Iq . {} \; -exec sed -i 's/\(#include\s*\(<\|\"\)\s*\)'"${header}"'\(\s*\(>\|\"\)\)/\1AMReX_'"${header}"'\3/g' {} +
  I_HEADER=$((I_HEADER+1))
  echo -ne "Progress: ${I_HEADER}/${NUM_HEADERS}"\\r
done
