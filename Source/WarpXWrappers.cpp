
#include <AMReX.H>
#include <AMReX_BLProfiler.H>

#include <WarpXWrappers.h>
#include <WarpX.H>

namespace 
{
    double** getMultiFabPointers(const amrex::MultiFab& mf, int *num_boxes, int *ngrow, int **shapes)
    {
        *ngrow = mf.nGrow();
        *num_boxes = mf.local_size();
        *shapes = (int*) malloc(BL_SPACEDIM * (*num_boxes) * sizeof(int));
        double** data = (double**) malloc((*num_boxes) * sizeof(double*));
        
        int i = 0;
        for ( amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi, ++i ) {
            data[i] = (double*) mf[mfi].dataPtr();
            for (int j = 0; j < BL_SPACEDIM; ++j) {
                (*shapes)[BL_SPACEDIM*i+j] = mf[mfi].box().length(j); 
            }
        }
        return data;
    }
}

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

    void warpx_addNParticles(int speciesnumber, int lenx,
                             double* x, double* y, double* z,
                             double* vx, double* vy, double* vz,
                             int nattr, double* attr, int uniqueparticles)
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

    double** warpx_getEfield(int lev, int direction,
                             int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getEfield(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    double** warpx_getBfield(int lev, int direction,
                             int *return_size, int *ngrow, int **shapes) {

        auto & mf = WarpX::GetInstance().getBfield(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    double** warpx_getCurrentDensity(int lev, int direction,
                                     int *return_size, int *ngrow, int **shapes) {

        auto & mf = WarpX::GetInstance().getcurrent(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }
    
    double** warpx_getParticleStructs(int speciesnumber,
                                      int* num_tiles, int** particles_per_tile) {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);

        const int level = 0;
        *num_tiles = myspc.numLocalTilesAtLevel(level);
        *particles_per_tile = (int*) malloc(*num_tiles*sizeof(int));
        
        double** data = (double**) malloc(*num_tiles*sizeof(typename WarpXParticleContainer::ParticleType*));
        int i = 0;
        for (WarpXParIter pti(myspc, level); pti.isValid(); ++pti, ++i) {
            auto& aos = pti.GetArrayOfStructs();
            data[i] = (double*) aos.data();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    double** warpx_getParticleArrays(int speciesnumber, int comp,
                                     int* num_tiles, int** particles_per_tile) {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);

        const int level = 0;
        *num_tiles = myspc.numLocalTilesAtLevel(level);
        *particles_per_tile = (int*) malloc(*num_tiles*sizeof(int));
        
        double** data = (double**) malloc(*num_tiles*sizeof(double*));
        int i = 0;
        for (WarpXParIter pti(myspc, level); pti.isValid(); ++pti, ++i) {
            auto& soa = pti.GetStructOfArrays();
            data[i] = (double*) soa[comp].dataPtr();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

}

