#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_REAL.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_Redistribution.H>
#include <AMReX_EBMultiFabUtil.H>

namespace amrex {

#if (AMREX_SPACEDIM > 1)
    //
    // Do small cell redistribution on one FAB
    //
    void apply_eb_redistribution ( const Box& bx,
                                   MultiFab& div_mf,
                                   MultiFab& divc_mf,
                                   const MultiFab& weights,
                                   MFIter* mfi,
                                   int icomp,
                                   int ncomp,
                                   const EBCellFlagFab& flags_fab,
                                   const MultiFab* volfrac,
                                   Box& /*domain*/,
                                   const Geometry & geom,
                                   bool use_wts_in_divnc)
    {
        //
        // Check that grid is uniform
        //
        const Real* dx = geom.CellSize();

#if (AMREX_SPACEDIM == 2)
        if (! amrex::almostEqual(dx[0], dx[1])) {
            amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
        }
#elif (AMREX_SPACEDIM == 3)
        if( ! amrex::almostEqual(dx[0],dx[1]) ||
            ! amrex::almostEqual(dx[1],dx[2]) ) {
            amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
        }
#endif

        //
        // Get array4 from arguments
        //
        Array4<Real> const& div  = div_mf.array(*mfi);
        Array4<Real> const& divc = divc_mf.array(*mfi);
        auto const&         wt   = weights.array(*mfi);
        auto const&        flags = flags_fab.array();
        auto const&        vfrac = volfrac->array(*mfi);

        apply_flux_redistribution ( bx, div, divc, wt, icomp, ncomp, flags, vfrac, geom, use_wts_in_divnc);
    }

    //
    // Do small cell redistribution on a MultiFab -- with a weighting function
    //
    void single_level_weighted_redistribute (MultiFab& div_tmp_in, MultiFab& div_out, const MultiFab& weights,
                                             int div_comp, int ncomp, const Geometry& geom, bool use_wts_in_divnc)
    {
        Box domain(geom.Domain());

        int nghost = 2;
        AMREX_ASSERT(div_tmp_in.nGrowVect().allGE(nghost));

        EB_set_covered(div_tmp_in, 0, ncomp, div_tmp_in.nGrow(), eb_covered_val);

        div_tmp_in.FillBoundary(geom.periodicity());

        // Here we take care of both the regular and covered cases ... all we do below is the cut cell cases
        MultiFab::Copy(div_out, div_tmp_in, 0, div_comp, ncomp, 0);

        // Get EB geometric info
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(div_out.Factory());
        const MultiFab* volfrac   = &(ebfactory. getVolFrac());

        for (MFIter mfi(div_out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            const Box& bx = mfi.tilebox ();

            // this is to check efficiently if this tile contains any eb stuff
            const auto&  div_fab = static_cast<EBFArrayBox const&>(div_out[mfi]);
            const EBCellFlagFab&  flags = div_fab.getEBCellFlagFab();

            if ( !(flags.getType(amrex::grow(bx,     0)) == FabType::covered) &&
                 !(flags.getType(amrex::grow(bx,nghost)) == FabType::regular) )
            {
                // Compute div(tau) with EB algorithm
                apply_eb_redistribution(bx, div_out, div_tmp_in, weights, &mfi,
                                        div_comp, ncomp, flags, volfrac, domain, geom, use_wts_in_divnc);

            }
        }
    }

    //
    // Do small cell redistribution on a MultiFab -- without a weighting function
    //
    void single_level_redistribute (MultiFab& div_tmp_in, MultiFab& div_out,
                                    int div_comp, int ncomp, const Geometry& geom)
    {
        // We create a weighting array to use inside the redistribution array
        MultiFab weights(div_out.boxArray(), div_out.DistributionMap(), 1, div_tmp_in.nGrow());
        weights.setVal(1.0);

        bool use_wts_in_divnc = false;
        single_level_weighted_redistribute (div_tmp_in, div_out, weights, div_comp, ncomp, geom, use_wts_in_divnc);
    }
#endif

} // end namespace
