
#include <AMReX.H>
#include <AMReX_BLProfiler.H>

#include <WarpXWrappers.h>
#include <WarpXParticleContainer.H>
#include <WarpX.H>
#include <WarpXUtil.H>
#include <WarpX_py.H>

namespace 
{
    double** getMultiFabPointers(const amrex::MultiFab& mf, int *num_boxes, int *ngrow, int **shapes)
    {
        *ngrow = mf.nGrow();
        *num_boxes = mf.local_size();
        *shapes = (int*) malloc(AMREX_SPACEDIM * (*num_boxes) * sizeof(int));
        double** data = (double**) malloc((*num_boxes) * sizeof(double*));
        
        int i = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi, ++i ) {
            data[i] = (double*) mf[mfi].dataPtr();
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                (*shapes)[AMREX_SPACEDIM*i+j] = mf[mfi].box().length(j); 
            }
        }
        return data;
    }
    int* getMultiFabLoVects(const amrex::MultiFab& mf, int *num_boxes, int *ngrow)
    {
        *ngrow = mf.nGrow();
        *num_boxes = mf.local_size();
        int *loVects = (int*) malloc((*num_boxes)*AMREX_SPACEDIM * sizeof(int));

        int i = 0;
        for ( amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi, ++i ) {
            const int* loVect = mf[mfi].loVect();
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                loVects[AMREX_SPACEDIM*i+j] = loVect[j];
            }
        }
        return loVects;
    }
}

extern "C"
{

    int warpx_nSpecies()
    {
	auto & mypc = WarpX::GetInstance().GetPartContainer();
        return mypc.nSpecies();
    }

    bool warpx_use_fdtd_nci_corr()
    {
	return WarpX::use_fdtd_nci_corr;
    }

    int warpx_l_lower_order_in_v()
    {
	return WarpX::l_lower_order_in_v;
    }

    int warpx_nComps() 
    {
        return PIdx::nattribs;        
    }

    int warpx_SpaceDim() 
    {
        return AMREX_SPACEDIM;
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
        if (warpx_py_afterinit) warpx_py_afterinit();
        if (warpx_py_particleloader) warpx_py_particleloader();
    }

    void warpx_finalize ()
    {
	WarpX::ResetInstance();
    }

    void warpx_set_callback_py_afterinit (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_afterinit = callback;
    }
    void warpx_set_callback_py_beforeEsolve (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_beforeEsolve = callback;
    }
    void warpx_set_callback_py_afterEsolve (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_afterEsolve = callback;
    }
    void warpx_set_callback_py_beforedeposition (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_beforedeposition = callback;
    }
    void warpx_set_callback_py_afterdeposition (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_afterdeposition = callback;
    }
    void warpx_set_callback_py_particlescraper (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_particlescraper = callback;
    }
    void warpx_set_callback_py_particleloader (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_particleloader = callback;
    }
    void warpx_set_callback_py_beforestep (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_beforestep = callback;
    }
    void warpx_set_callback_py_afterstep (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_afterstep = callback;
    }
    void warpx_set_callback_py_afterrestart (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_afterrestart = callback;
    }
    void warpx_set_callback_py_particleinjection (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_particleinjection = callback;
    }
    void warpx_set_callback_py_appliedfields (WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        warpx_py_appliedfields = callback;
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
	myspc.AddNParticles(lev, lenx, x, y, z, vx, vy, vz, nattr, attr, uniqueparticles);
    }

    void warpx_ConvertLabParamsToBoost()
    {
      ConvertLabParamsToBoost();
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

    int* warpx_getEfieldLoVects(int lev, int direction,
                                int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getEfield(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getEfieldCP(int lev, int direction,
                               int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getEfield_cp(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    int* warpx_getEfieldCPLoVects(int lev, int direction,
                                  int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getEfield_cp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getEfieldFP(int lev, int direction,
                               int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getEfield_fp(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    int* warpx_getEfieldFPLoVects(int lev, int direction,
                                  int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getEfield_fp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getBfield(int lev, int direction,
                             int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getBfield(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    int* warpx_getBfieldLoVects(int lev, int direction,
                                int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getBfield(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getBfieldCP(int lev, int direction,
                               int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getBfield_cp(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    int* warpx_getBfieldCPLoVects(int lev, int direction,
                                  int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getBfield_cp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getBfieldFP(int lev, int direction,
                               int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getBfield_fp(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    int* warpx_getBfieldFPLoVects(int lev, int direction,
                                  int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getBfield_fp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getCurrentDensity(int lev, int direction,
                                     int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getcurrent(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    int* warpx_getCurrentDensityLoVects(int lev, int direction,
                                        int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getcurrent(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getCurrentDensityCP(int lev, int direction,
                                       int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getcurrent_cp(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    int* warpx_getCurrentDensityCPLoVects(int lev, int direction,
                                          int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getcurrent_cp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getCurrentDensityFP(int lev, int direction,
                                       int *return_size, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getcurrent_fp(lev, direction);
        return getMultiFabPointers(mf, return_size, ngrow, shapes);
    }

    int* warpx_getCurrentDensityFPLoVects(int lev, int direction,
                                          int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getcurrent_fp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    double** warpx_getParticleStructs(int speciesnumber,
                                      int* num_tiles, int** particles_per_tile) {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);

        const int level = 0;

        int i = 0;
        for (WarpXParIter pti(myspc, level); pti.isValid(); ++pti, ++i) {}

        // *num_tiles = myspc.numLocalTilesAtLevel(level);
        *num_tiles = i;
        *particles_per_tile = (int*) malloc(*num_tiles*sizeof(int));

        double** data = (double**) malloc(*num_tiles*sizeof(typename WarpXParticleContainer::ParticleType*));
        i = 0;
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

        int i = 0;
        for (WarpXParIter pti(myspc, level); pti.isValid(); ++pti, ++i) {}

        // *num_tiles = myspc.numLocalTilesAtLevel(level);
        *num_tiles = i;
        *particles_per_tile = (int*) malloc(*num_tiles*sizeof(int));

        double** data = (double**) malloc(*num_tiles*sizeof(double*));
        i = 0;
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

    void warpx_EvolveE (double dt) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveE (dt);
    }
    void warpx_EvolveB (double dt) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveB (dt);
    }
    void warpx_FillBoundaryE () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.FillBoundaryE ();
    }
    void warpx_FillBoundaryB () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.FillBoundaryB ();
    }
    void warpx_SyncCurrent () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.SyncCurrent ();
    }
    void warpx_UpdateAuxilaryData () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.UpdateAuxilaryData ();
    }
    void warpx_PushParticlesandDepose (double cur_time) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.PushParticlesandDepose (cur_time);
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

