
// NOTE: Since super class (LSCoreBase) constructor has no arguments, it is
// called by default.
template<class F>
LSCore<F>::LSCore (EB2::GeometryShop<F> a_gshop)
    : m_gshop(a_gshop)
{}


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

    // Make sure EB geometry exists (build if necessary)
    if (rebuild_eb[lev] == 1) {
        int max_level = 0;
        EB2::Build(m_gshop, geom[lev], max_level, max_level);
        const EB2::IndexSpace & cylinder_ebis = EB2::IndexSpace::top();
        eb_levels[lev] = & cylinder_ebis.getLevel(geom[lev]);

        rebuild_eb[lev] = 0;
    }


    BoxArray ba_nd = amrex::convert(ba, IntVect{1, 1, 1});


    // Define MultiFab data (will store final level set for this level)
    level_set[lev].define(ba_nd, dm, ncomp, nghost);
    level_set_cc[lev].define(ba, dm, ncomp, nghost);


    // Construct LSFactory (does the actual level-set building)
    // TODO: more light-weight structure here? (since we don't need all the
    // LSFactory features for LSCore)
    int levelset_pad    = nghost;
    int levelset_eb_pad = eb_pad * 2 * lev;
    ls_factory[lev] = std::unique_ptr<LSFactory>(
        new LSFactory(lev, 1, 1, levelset_pad, levelset_eb_pad, ba, geom[lev], dm)
        );


    // Compute LSFactory
    int eb_pad = ls_factory[lev]->get_eb_pad();
    EBFArrayBoxFactory eb_factory(* eb_levels[lev], geom[lev], ba, dm,
                                  {eb_pad, eb_pad, eb_pad}, EBSupport::full);


    GShopLSFactory<F> ls_gshop(m_gshop, * ls_factory[lev]);
    std::unique_ptr<MultiFab> cylinder_mf_impfunc = ls_gshop.fill_impfunc();
    ls_factory[lev]->intersection_ebf(eb_factory, * cylinder_mf_impfunc);

    int ng = ls_factory[lev]->get_ls_pad();
    level_set[lev].copy(* ls_factory[lev]->get_data(), 0, 0, 1, ng, ng);

    // VisMF::Write(level_set[lev], "LevelSet_L");

    Print() << "made level: " << lev << std::endl;

    level_set[lev].FillBoundary(geom[lev].periodicity());

}
