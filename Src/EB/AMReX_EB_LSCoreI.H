
// NOTE: Since super class (LSCoreBase) constructor has no arguments, it is
// called by default.
template<class F>
LSCore<F>::LSCore (EB2::GeometryShop<F> a_gshop)
: m_gshop(a_gshop)
{}


// NOTE: here we need to explicitly call the super call (LSCoreBase)
// constructor due to the extra arguments
template<class F>
LSCore<F>::LSCore (EB2::GeometryShop<F> a_gshop,
                   const RealBox * rb, int max_level_in, const Vector<int> & n_cell_in, int coord)
: LSCoreBase(rb, max_level_in, n_cell_in, coord), m_gshop(a_gshop)
{}



template<class F>
LSCore<F>::~LSCore () {}


// Delete level data overrides the pure virtual function in AmrCore
template<class F>
void LSCore<F>::ClearLevel (int lev) {
    level_set[lev].clear();
    level_set_valid[lev].clear();

    // NOTE: maybe I shouldn't do this (this method is being called by regrid...)
    //ls_factory[lev].clear();
}

// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. Overrides the pure
// virtual function in AmrCore
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
        int max_crse_level = lev; // Important! this is needed to allow
                                  // EBFactory to be built on coarser levels
                                  // (which is needed by FillCoarsePatch)
        int required_crse_lev = std::max(0, lev - 1);
        //int required_crse_lev = lev;
        EB2::Build(m_gshop, geom[lev], required_crse_lev, max_crse_level);
        const EB2::IndexSpace & cylinder_ebis = EB2::IndexSpace::top();
        eb_levels[lev] = & cylinder_ebis.getLevel(geom[lev]);

        rebuild_eb[lev] = 0;
    }

    /****************************************************************************
     * Define Basic data MultiFabs                                              *
     ***************************************************************************/

    BoxArray ba_nd = amrex::convert(ba, IntVect{1, 1, 1});

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

        ls_factory[lev]->fill(eb_factory, * mf_impfunc, levelset_eb_pad);
        level_set[lev].copy(* ls_factory[lev]->get_data(), 0, 0, 1, ng, ng);

    } else {

        //_______________________________________________________________________
        // Build finer levels locally, using coarse level set

        MakeNewLevelFromCoarse(lev, time, ba, dm);

        EBFArrayBoxFactory eb_factory(* eb_levels[lev], geom[lev], ba, dm,
                                      {max_eb_pad + 1, max_eb_pad + 1, max_eb_pad + 1},
                                      EBSupport::full);

        // EB boundary-centre data
        const MultiCutFab & bndrycent = eb_factory.getBndryCent();
        const auto & flags = eb_factory.getMultiEBCellFlagFab();

        MultiFab normal(ba, dm, 3, max_eb_pad + 1);
        FillEBNormals(normal, eb_factory, geom[lev]);

        iMultiFab eb_valid(ba, dm, 1, max_eb_pad + 1);
        eb_valid.setVal(0);

        // TODO: some of this might be shifted to the base class?

#if (AMREX_SPACEDIM == 1)
        Real min_dx = geom[lev].CellSize(0);
#elif (AMREX_SPACEDIM == 2)
        Real min_dx = std::min(geom[lev].CellSize(0), geom[lev].CellSize(1));
#elif (AMREX_SPACEDIM == 3)
        Real min_dx = std::min(geom[lev].CellSize(0),
                               std::min(geom[lev].CellSize(1), geom[lev].CellSize(2)));
#endif
        Real ls_threshold = min_dx * levelset_eb_pad;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(level_set[lev], IntVect{AMREX_D_DECL(levelset_eb_pad,
                                                             levelset_eb_pad,
                                                             levelset_eb_pad)});
             mfi.isValid();
             ++mfi)
        {
            const auto & ls_tile = level_set[lev][mfi];
                  bool bail      = false;
                  Box tile_box   = mfi.tilebox();
                  Box eb_search  = EBSearchBox(tile_box, ls_tile, geom[lev], bail);

            if (bail) continue;

            int n_facets = 0;
            const auto & flag = flags[mfi];
            // Need to count number of eb-facets (in order to allocate facet_list)
            amrex_eb_count_facets(BL_TO_FORTRAN_BOX(eb_search),
                                  BL_TO_FORTRAN_3D(flag),
                                  & n_facets);

            int facet_list_size = 6 * n_facets;
            Vector<Real> facet_list(facet_list_size);


                  auto & ls_tile_w = level_set[lev][mfi];
                  auto & v_tile    = eb_valid[mfi];
            const auto & if_tile   = (* mf_impfunc)[mfi];

            if (n_facets > 0) {
                const auto & norm_tile = normal[mfi];
                const auto & bcent_tile = bndrycent[mfi];

                int c_facets = 0;
                amrex_eb_as_list(BL_TO_FORTRAN_BOX(eb_search), & c_facets,
                                 BL_TO_FORTRAN_3D(flag),
                                 BL_TO_FORTRAN_3D(norm_tile),
                                 BL_TO_FORTRAN_3D(bcent_tile),
                                 facet_list.dataPtr(), & facet_list_size,
                                 geom[lev].CellSize()                      );

                amrex_eb_fill_levelset_loc(BL_TO_FORTRAN_BOX(tile_box),
                                           facet_list.dataPtr(), & facet_list_size,
                                           BL_TO_FORTRAN_3D(v_tile),
                                           BL_TO_FORTRAN_3D(ls_tile_w),
                                           BL_TO_FORTRAN_3D(ls_tile), & ls_threshold,
                                           geom[lev].CellSize(), geom[lev].CellSize() );

            }


            //_______________________________________________________________________
            // Enforce threshold of local level-set
            amrex_eb_threshold_levelset(BL_TO_FORTRAN_BOX(tile_box), & ls_threshold,
                                        BL_TO_FORTRAN_3D(ls_tile_w));


            //_______________________________________________________________________
            // Validate level-set
            const int ls_grid_ref = 1;
            amrex_eb_validate_levelset(BL_TO_FORTRAN_BOX(tile_box), & ls_grid_ref,
                                       BL_TO_FORTRAN_3D(if_tile),
                                       BL_TO_FORTRAN_3D(v_tile),
                                       BL_TO_FORTRAN_3D(ls_tile_w)   );
        }
    }

    level_set[lev].FillBoundary(geom[lev].periodicity());

    Real end_time = ParallelDescriptor::second() - strt_time;
    amrex::Print() << "** Time to build the level set at level " << lev
                   << " " << end_time << std::endl;
}
