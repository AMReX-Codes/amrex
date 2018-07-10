#define XDIR 0
#define YDIR 1
#define ZDIR 2
#define DIM  3

using namespace amrex;

int  nProcs;
int  myProc;
int  IOProc;

int  verbose;

int  nVars;

int  div_free;
int  transpose_dp;
int  density_weighting;

int  use_cutoff_density;
Real cutoff_density;

int  BLix, BLjx, BLkx;
int  FTix, FTjx, FTkx;
int  FTis, FTjs, FTks;
int  FTmx, FThkxpo;
int  wavenumbers;

int  Qix, Qjx, Qkx;
Vector<Real*> spectrum;
Vector<Real*> spectrumS;
Vector<Real*> spectrumC;
Vector<Real*> Qx, Qy, Qz;
Real *sum, *sum2;

// Box probDomain;
Vector<Real> probLo;
Vector<Real> probHi;
Real Lx, Ly, Lz;
Real dx, dy, dz;

std::string infile;
Vector<std::string> whichVar;
Real Time;
int  timeSteps;

int  local_ix, local_i_start;
int  local_jx_after_transpose, local_j_start_after_transpose;
int  total_local_size;
Vector<Real> local_xlo;
Vector<Real> local_xhi;

Vector<fftw_real*>    local_data;
Vector<fftw_complex*> local_data_c;

rfftwnd_mpi_plan plan_real2cplx, plan_cplx2real; 

void plan_ffts(Box &probDomain, BoxArray &domainBoxArray, Vector<int> &pmap);

void Spectra(MultiFab &mf, Box &probDomain);

int        do_filter;
int        nFilters;
Vector<int> filterWN;
void       filter(MultiFab &mf, Box &probDomain);
