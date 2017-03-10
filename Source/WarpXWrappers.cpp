
#include <AMReX.H>
#include <AMReX_BLProfiler.H>

#include <WarpXWrappers.h>
#include <WarpX.H>

extern "C"
{
    void amrex_init (int argc, char* argv[])
    {
	amrex::Initialize(argc,argv);
    }

#ifdef BL_USE_MPI
    void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm mpicomm)
    {
	amrex::Initialize(argc,argv,true,mpicomm);	
    }
#endif

    void amrex_finalize (int finalize_mpi)
    {
	amrex::Finalize(finalize_mpi);
    }

    void warpx_init ()
    {
	WarpX& warpx = WarpX::GetInstance();
	warpx.InitData();
    }

    void warpx_finalize ()
    {
	WarpX::ResetInstance();
    }

    void warpx_evolve (int numsteps)
    {
	WarpX& warpx = WarpX::GetInstance();
	warpx.Evolve(numsteps);
    }

    void addNParticles(int speciesnumber, int lenx, double* x, double* y, double* z, double* vx, double* vy, double* vz, int nattr, double* attr, int uniqueparticles)
    {
	auto & mypc = WarpX::GetInstance().GetPartContainer();
	auto & myspc = mypc.GetParticleContainer(speciesnumber);
	myspc.AddNParticles(lenx, x, y, z, vx, vy, vz, nattr, attr, uniqueparticles);
    }

    double warpx_getProbLo(int dir)
    {
      WarpX& warpx = WarpX::GetInstance();
      const amrex::Geometry& geom = warpx.Geom(0);
      return geom.ProbLo(dir);
    }

    double warpx_getProbHi(int dir)
    {
      WarpX& warpx = WarpX::GetInstance();
      const amrex::Geometry& geom = warpx.Geom(0);
      return geom.ProbHi(dir);
    }

  long warpx_getNumParticles(int speciesnumber) {
    auto & mypc = WarpX::GetInstance().GetPartContainer();
    auto & myspc = mypc.GetParticleContainer(speciesnumber);
    return myspc.TotalNumberOfParticles();
  }

  double* warpx_getParticlePositions(int speciesnumber) {
    amrex::Array<amrex::Real> part_data;
    auto & mypc = WarpX::GetInstance().GetPartContainer();
    auto & myspc = mypc.GetParticleContainer(speciesnumber);
    myspc.GetParticleLocations(part_data);
    double* buffer = (double*) malloc(part_data.size()*sizeof(double));
    for (long i = 0; i < part_data.size(); ++i) {
      buffer[i] = part_data[i];
    }
    return buffer;
  }
}

