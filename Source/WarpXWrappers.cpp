
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

    void addNParticles(int speciesnumber, int lenx,
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
                             int *return_size, int **shapes) {

        auto & mf = WarpX::GetInstance().getEfield(lev, direction);

        int num_boxes = mf.local_size();
        *return_size = num_boxes;
        *shapes = (int*) malloc(3*num_boxes*sizeof(int));
        
        double** data = (double**) malloc(num_boxes*sizeof(double*));

        int i = 0;
        for ( amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi, ++i ) {
            data[i] = (double*) mf[mfi].dataPtr();
            for (int j = 0; j < 3; ++j) {
                (*shapes)[3*i+j] = mf[mfi].box().length(j); 
            }
        }
        return data;
    }
    
    double* warpx_getParticlePositions(int speciesnumber) {
        amrex::Array<amrex::Real> *part_data_ptr = new amrex::Array<amrex::Real>;
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);
        myspc.GetParticleLocations(*part_data_ptr);
        return (double*) part_data_ptr->dataPtr();
    }

    double* warpx_getParticleData(int speciesnumber, int start_comp, int num_comp) {
        amrex::Array<amrex::Real> *part_data_ptr = new amrex::Array<amrex::Real>;
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);
        myspc.GetArrayData(*part_data_ptr, start_comp, num_comp);
        return (double*) part_data_ptr->dataPtr();        
    }

    int* warpx_getParticleIDs(int speciesnumber) {
        amrex::Array<int> *part_data_ptr = new amrex::Array<int>;
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);
        myspc.GetParticleIDs(*part_data_ptr);
        return (int*) part_data_ptr->dataPtr();        
    }

    int* warpx_getParticleCPU(int speciesnumber) {
        amrex::Array<int> *part_data_ptr = new amrex::Array<int>;
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);
        myspc.GetParticleCPU(*part_data_ptr);
        return (int*) part_data_ptr->dataPtr();        
    }

}

