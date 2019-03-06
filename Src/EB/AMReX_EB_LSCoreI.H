
//! NOTE: Since super class (LSCoreBase) constructor has no arguments, it is
//! called by default.
template<class F>
LSCore<F>::LSCore (EB2::GeometryShop<F> a_gshop)
: m_gshop(a_gshop)
{}


//! NOTE: here we need to explicitly call the super call (LSCoreBase)
//! constructor due to the extra arguments
template<class F>
LSCore<F>::LSCore (EB2::GeometryShop<F> a_gshop,
                   const RealBox * rb, int max_level_in, const Vector<int> & n_cell_in, int coord)
: LSCoreBase(rb, max_level_in, n_cell_in, coord), m_gshop(a_gshop)
{}



template<class F>
LSCore<F>::~LSCore () {}


//! Delete level data overrides the pure virtual function in AmrCore
template<class F>
void LSCore<F>::ClearLevel (int lev) {
    level_set[lev].clear();
    level_set_valid[lev].clear();

    // NOTE: maybe I shouldn't do this (this method is being called by regrid...)
    //ls_factory[lev].clear();
}


template<class F>
const EB2::Level * LSCore<F>::BuildEBLevel ( const EB2::GeometryShop<F> gshop,
                                             int lev, const Vector<Geometry> & geom,
                                             int max_crse_level, int required_crse_lev) {
    EB2::Build(gshop, geom[lev], required_crse_lev, max_crse_level);
    const EB2::IndexSpace & ebis = EB2::IndexSpace::top();

    return & ebis.getLevel(geom[lev]);
}



//! Build a new EB level and return a pointer to this level
template<class F>
const EB2::Level * LSCore<F>::BuildEBLevel (const EB2::GeometryShop<F> & gshop,
                                            int lev, const Vector<Geometry> & geom ) {
    int max_crse_level = lev; //!< Important! this is needed to allow EBFactory to
                              //!< be built on coarser levels (which is needed by
                              //!< FillCoarsePatch)
    int required_crse_lev = std::max(0, lev - 1);
    //int required_crse_lev = lev;

    return LSCore<F>::BuildEBLevel (gshop, lev, geom, max_crse_level, required_crse_lev);
}


//! Make a new level from scratch using provided BoxArray and
//! DistributionMapping. Only used during initialization. Overrides the pure
//! virtual function in AmrCore
template<class F>
void LSCore<F>::MakeNewLevelFromScratch (int lev, Real time, const BoxArray & ba,
                                         const DistributionMapping & dm) {

    const int ncomp  = 1;
    const int nghost = 1;

    Real strt_time = ParallelDescriptor::second();

    /****************************************************************************
     * Make sure EB geometry exists (build if necessary)                        *
     ***************************************************************************/

    if (rebuild_eb[lev] == 1) {
        eb_levels[lev] = LSCore<F>::BuildEBLevel(m_gshop, lev, geom);

        rebuild_eb[lev] = 0;
    }

    /****************************************************************************
     * Define Basic data MultiFabs                                              *
     ***************************************************************************/

    BoxArray ba_nd = amrex::convert(ba, IntVect::TheUnitVector());

    // Define MultiFab data (will store final level set for this level)
    level_set[lev].define(ba_nd, dm, ncomp, nghost);

    // At this point, we consider _everywhere_ as valid. This is maintained for
    // legacy reasons. TODO: There might be a better way of doing things.
    level_set_valid[lev].define(ba_nd, dm, ncomp, nghost);
    level_set_valid[lev].setVal(1);

    /****************************************************************************
     * Construct LSFactory (does the actual level-set building on lev 0)        *
     ***************************************************************************/

    // TODO: more light-weight structure here? (since we don't need all the
    // LSFactory features for LSCore)
    int levelset_pad    = nghost;
    int levelset_eb_pad = eb_pad;

    // Basic search radius
    for (int i = lev - 1; i >= 0; i--){
        IntVect rr = refRatio(i);
        levelset_eb_pad *= std::max({AMREX_D_DECL(rr[0], rr[1], rr[2])});
    }
    if (levelset_eb_pad > max_eb_pad) levelset_eb_pad = max_eb_pad;

    ls_factory[lev] = std::unique_ptr<LSFactory>(
        new LSFactory(lev, 1, 1, levelset_pad, levelset_eb_pad, ba, geom[lev], dm)
        );

    GShopLSFactory<F> ls_gshop(m_gshop, * ls_factory[lev]);
    std::unique_ptr<MultiFab> mf_impfunc = ls_gshop.fill_impfunc();

    /****************************************************************************
     * Fill the actual level-set MultiFabs (one for each level)                 *
     ***************************************************************************/

    if (lev == 0) {

        //_______________________________________________________________________
        // Compute level 0 using LSFactory

        int eb_pad = ls_factory[lev]->get_eb_pad();
        EBFArrayBoxFactory eb_factory(* eb_levels[lev], geom[lev], ba, dm,
                                      {eb_pad, eb_pad, eb_pad}, EBSupport::full);

        int ng = ls_factory[lev]->get_ls_pad();

        ls_factory[lev]->Fill(eb_factory, * mf_impfunc, levelset_eb_pad);
        level_set[lev].copy(* ls_factory[lev]->get_data(), 0, 0, 1, ng, ng);

    } else {

        //_______________________________________________________________________
        // Build finer levels locally, using coarse level set

        // Fills level-set[lev] with coarse data
        MakeNewLevelFromCoarse(lev, time, ba, dm);

        EBFArrayBoxFactory eb_factory(* eb_levels[lev], geom[lev], ba, dm,
                                      {levelset_eb_pad + 1, levelset_eb_pad + 1, levelset_eb_pad + 1},
                                      EBSupport::full);

        IntVect ebt_size{AMREX_D_DECL(levelset_eb_pad, levelset_eb_pad, levelset_eb_pad)};
        LSCoreBase::FillLevelSet(level_set[lev], level_set[lev], eb_factory, (* mf_impfunc),
                                 ebt_size, levelset_eb_pad, geom[lev]);
    }

    Real end_time = ParallelDescriptor::second() - strt_time;
    amrex::Print() << "** Time to build the level set at level " << lev
                   << " " << end_time << std::endl;
}
