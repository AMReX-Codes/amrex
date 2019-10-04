
#include <WarpXWrappers.h>
#include <WarpXParticleContainer.H>
#include <WarpX.H>
#include <WarpXUtil.H>
#include <WarpX_py.H>

namespace
{
    amrex::Real** getMultiFabPointers(const amrex::MultiFab& mf, int *num_boxes, int *ncomps, int *ngrow, int **shapes)
    {
        *ncomps = mf.nComp();
        *ngrow = mf.nGrow();
        *num_boxes = mf.local_size();
        int shapesize = AMREX_SPACEDIM;
        if (mf.nComp() > 1) shapesize += 1;
        *shapes = (int*) malloc(shapesize * (*num_boxes) * sizeof(int));
        amrex::Real** data = (amrex::Real**) malloc((*num_boxes) * sizeof(amrex::Real*));

#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi ) {
            int i = mfi.LocalIndex();
            data[i] = (amrex::Real*) mf[mfi].dataPtr();
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                (*shapes)[shapesize*i+j] = mf[mfi].box().length(j);
            }
            if (mf.nComp() > 1) (*shapes)[shapesize*i+AMREX_SPACEDIM] = mf.nComp();
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

    int warpx_Real_size()
    {
        return (int)sizeof(amrex::Real);
    }

    int warpx_ParticleReal_size()
    {
        return (int)sizeof(amrex::ParticleReal);
    }

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
        amrex::Finalize();
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
                             amrex::ParticleReal* x, amrex::ParticleReal* y, amrex::ParticleReal* z,
                             amrex::ParticleReal* vx, amrex::ParticleReal* vy, amrex::ParticleReal* vz,
                             int nattr, amrex::ParticleReal* attr, int uniqueparticles)
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

    amrex::Real warpx_getProbLo(int dir)
    {
      WarpX& warpx = WarpX::GetInstance();
      const amrex::Geometry& geom = warpx.Geom(0);
      return geom.ProbLo(dir);
    }

    amrex::Real warpx_getProbHi(int dir)
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

    amrex::Real** warpx_getEfield(int lev, int direction,
                             int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getEfield(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getEfieldLoVects(int lev, int direction,
                                int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getEfield(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::Real** warpx_getEfieldCP(int lev, int direction,
                               int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getEfield_cp(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getEfieldCPLoVects(int lev, int direction,
                                  int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getEfield_cp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::Real** warpx_getEfieldFP(int lev, int direction,
                               int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getEfield_fp(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getEfieldFPLoVects(int lev, int direction,
                                  int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getEfield_fp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::Real** warpx_getBfield(int lev, int direction,
                             int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getBfield(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getBfieldLoVects(int lev, int direction,
                                int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getBfield(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::Real** warpx_getBfieldCP(int lev, int direction,
                               int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getBfield_cp(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getBfieldCPLoVects(int lev, int direction,
                                  int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getBfield_cp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::Real** warpx_getBfieldFP(int lev, int direction,
                               int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getBfield_fp(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getBfieldFPLoVects(int lev, int direction,
                                  int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getBfield_fp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::Real** warpx_getCurrentDensity(int lev, int direction,
                                     int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getcurrent(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getCurrentDensityLoVects(int lev, int direction,
                                        int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getcurrent(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::Real** warpx_getCurrentDensityCP(int lev, int direction,
                                       int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getcurrent_cp(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getCurrentDensityCPLoVects(int lev, int direction,
                                          int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getcurrent_cp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::Real** warpx_getCurrentDensityFP(int lev, int direction,
                                       int *return_size, int *ncomps, int *ngrow, int **shapes) {
        auto & mf = WarpX::GetInstance().getcurrent_fp(lev, direction);
        return getMultiFabPointers(mf, return_size, ncomps, ngrow, shapes);
    }

    int* warpx_getCurrentDensityFPLoVects(int lev, int direction,
                                          int *return_size, int *ngrow) {
        auto & mf = WarpX::GetInstance().getcurrent_fp(lev, direction);
        return getMultiFabLoVects(mf, return_size, ngrow);
    }

    amrex::ParticleReal** warpx_getParticleStructs(int speciesnumber, int lev,
                                      int* num_tiles, int** particles_per_tile) {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);

        int i = 0;
        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti, ++i) {}

        // *num_tiles = myspc.numLocalTilesAtLevel(lev);
        *num_tiles = i;
        *particles_per_tile = (int*) malloc(*num_tiles*sizeof(int));

        amrex::ParticleReal** data = (amrex::ParticleReal**) malloc(*num_tiles*sizeof(typename WarpXParticleContainer::ParticleType*));
        i = 0;
        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti, ++i) {
            auto& aos = pti.GetArrayOfStructs();
            data[i] = (amrex::ParticleReal*) aos.data();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    amrex::ParticleReal** warpx_getParticleArrays(int speciesnumber, int comp, int lev,
                                     int* num_tiles, int** particles_per_tile) {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        auto & myspc = mypc.GetParticleContainer(speciesnumber);

        int i = 0;
        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti, ++i) {}

        // *num_tiles = myspc.numLocalTilesAtLevel(lev);
        *num_tiles = i;
        *particles_per_tile = (int*) malloc(*num_tiles*sizeof(int));

        amrex::ParticleReal** data = (amrex::ParticleReal**) malloc(*num_tiles*sizeof(amrex::ParticleReal*));
        i = 0;
        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti, ++i) {
            auto& soa = pti.GetStructOfArrays();
            data[i] = (amrex::ParticleReal*) soa.GetRealData(comp).dataPtr();
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

    void warpx_EvolveE (amrex::Real dt) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveE (dt);
    }
    void warpx_EvolveB (amrex::Real dt) {
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
    void warpx_PushParticlesandDepose (amrex::Real cur_time) {
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
    amrex::Real warpx_gett_new (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.gett_new (lev);
    }
    void warpx_sett_new (int lev, amrex::Real time) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.sett_new (lev, time);
    }
    amrex::Real warpx_getdt (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.getdt (lev);
    }

    int warpx_maxStep () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.maxStep ();
    }
    amrex::Real warpx_stopTime () {
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

