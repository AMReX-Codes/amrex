
#include <AMReX.H>
#include <AMReX_BLProfiler.H>

#include <WarpXWrappers.h>
#include <WarpXParticleContainer.H>
#include <WarpX.H>
#include <WarpX_py.H>

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

    int warpx_nSpecies()
    {
	auto & mypc = WarpX::GetInstance().GetPartContainer();
        return mypc.nSpecies();
    }

    int warpx_nComps() 
    {
        return PIdx::nattribs;        
    }

    int warpx_SpaceDim() 
    {
        return BL_SPACEDIM;
    }

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

    void warpx_set_callback_py_funcs (WARPX_CALLBACK_PY_FUNC_1 print_step)
    {
        warpx_py_print_step = print_step;
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
        const int lev = 0;
        amrex::Abort("warpx_addNParticles needs to be updated to include lev argument");
	myspc.AddNParticles(lev, lenx, x, y, z, vx, vy, vz, nattr, attr, uniqueparticles);
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
            data[i] = (double*) soa.GetRealData(comp).dataPtr();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    void warpx_ComputeDt () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.ComputeDt ();
    }
    void warpx_MoveWindow () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.MoveWindow (true);
    }

    void warpx_EvolveE (int lev, double dt) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveE (lev, dt);
    }
    void warpx_EvolveB (int lev, double dt) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveB (lev, dt);
    }
    void warpx_FillBoundaryE (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.FillBoundaryE (lev);
    }
    void warpx_FillBoundaryB (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.FillBoundaryB (lev);
    }
    void warpx_PushParticlesandDepose (int lev, double cur_time) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.PushParticlesandDepose (lev, cur_time);
    }

    int warpx_getistep (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.getistep (lev);
    }
    void warpx_setistep (int lev, int ii) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.setistep (lev, ii);
    }
    double warpx_gett_new (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.gett_new (lev);
    }
    void warpx_sett_new (int lev, double time) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.sett_new (lev, time);
    }
    double warpx_getdt (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.getdt (lev);
    }

    int warpx_maxStep () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.maxStep ();
    }
    double warpx_stopTime () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.stopTime ();
    }

    int warpx_checkInt () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.checkInt ();
    }
    int warpx_plotInt () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.plotInt ();
    }

    void warpx_WriteCheckPointFile () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.WriteCheckPointFile ();
    }
    void warpx_WritePlotFile () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.WritePlotFile ();
    }

    int warpx_finestLevel () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.finestLevel ();
    }

    void mypc_Redistribute () {
	    auto & mypc = WarpX::GetInstance().GetPartContainer();
        mypc.Redistribute();
    }

}

