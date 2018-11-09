
// #include <AMReX_Vector.H>
// #include <AMReX_AmrCore.H>
// #include <AMReX_ParallelDescriptor.H>
// #include <AMReX_ParmParse.H>
// #include <AMReX_MultiFabUtil.H>
// #include <AMReX_FillPatchUtil.H>
// #include <AMReX_PlotFileUtil.H>
// #include <AMReX_VisMF.H>
// #include <AMReX_PhysBCFunct.H>
//
// #ifdef BL_USE_SENSEI_INSITU
// #include <AMReX_AmrMeshInSituBridge.H>
// #endif
//
// #ifdef BL_MEM_PROFILING
// #include <AMReX_MemProfiler.H>
// #endif


//using namespace amrex;
//namespace amrex {

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
// NOTE: Since super class (LSCoreBase) constructor has no arguments, it is
// called by default.
template<class F>
LSCore<F>::LSCore (EB2::GeometryShop<F> a_gshop)
    : m_gshop(a_gshop)
{

}


template<class F>
LSCore<F>::~LSCore () {}


// Delete level data overrides the pure virtual function in AmrCore
template<class F>
void LSCore<F>::ClearLevel (int lev) {
    level_set[lev].clear();
    // NOTE: maybe I shouldn't do this (this method is being called by regrid...)
    //ls_factory[lev].clear();
}


// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. overrides the pure
// virtual function in AmrCore
template<class F>
void LSCore<F>::MakeNewLevelFromScratch (int lev, Real time, const BoxArray & ba,
                                      const DistributionMapping & dm) {
    const int ncomp  = 1;
    const int nghost = 1;

    level_set[lev].define(ba, dm, ncomp, nghost);
    // NOTE: use move semantics since LSFactory contains unique_ptrs.
    int levelset_pad    = 1;
    int levelset_eb_pad = 1;
    ls_factory[lev] = std::unique_ptr<LSFactory>(
        new LSFactory(lev, 1, 1, levelset_pad, levelset_eb_pad, ba, geom[lev], dm)
        );

    const Real * dx = geom[lev].CellSize();

    // Number of particles
    //int np = particle_data.size();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(level_set[lev], true); mfi.isValid(); ++mfi) {
        const Box &  tile_box = mfi.tilebox();
              auto & phi_tile = level_set[lev][mfi];

        // fill_levelset_ib ( tile_box.loVect(), tile_box.hiVect(),
        //                    particle_data.dataPtr(), & np,
        //                    BL_TO_FORTRAN_3D(phi_tile),
        //                    dx                                     );

    }

    level_set[lev].FillBoundary(geom[lev].periodicity());

}



//}
