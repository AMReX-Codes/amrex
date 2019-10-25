#include <AMReX_EB_F.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB_utils.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>

namespace amrex {

    // 
    // Fill EB normals
    // 
    void FillEBNormals(MultiFab & normals, const EBFArrayBoxFactory & eb_factory,
                       const Geometry & geom) {

        BL_PROFILE("amrex::FillEBNormals()")

        BoxArray ba = normals.boxArray();
        DistributionMapping dm = normals.DistributionMap();
        int n_grow = normals.nGrow();

        // Dummy array for MFIter
        MultiFab dummy(ba, dm, 1, n_grow, MFInfo(), eb_factory);
        // Area fraction data
        std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac = eb_factory.getAreaFrac();

        const auto & flags = eb_factory.getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for(MFIter mfi(dummy, true); mfi.isValid(); ++mfi) {
            Box tile_box = mfi.growntilebox();
            const int * lo = tile_box.loVect();
            const int * hi = tile_box.hiVect();

            const auto & flag = flags[mfi];

            if (flag.getType(tile_box) == FabType::singlevalued) {
                // Target for compute_normals(...)
                auto & norm_tile = normals[mfi];
                // Area fractions in x, y, and z directions
                const auto & af_x_tile = (* areafrac[0])[mfi];
                const auto & af_y_tile = (* areafrac[1])[mfi];
                const auto & af_z_tile = (* areafrac[2])[mfi];

                amrex_eb_compute_normals(lo, hi,
                                         BL_TO_FORTRAN_3D(flag),
                                         BL_TO_FORTRAN_3D(norm_tile),
                                         BL_TO_FORTRAN_3D(af_x_tile),
                                         BL_TO_FORTRAN_3D(af_y_tile),
                                         BL_TO_FORTRAN_3D(af_z_tile)  );
            }
        }

        normals.FillBoundary(geom.periodicity());
    }

#if (AMREX_SPACEDIM > 1)
    // 
    // Do small cell redistribution on one FAB
    // 
    void apply_eb_redistribution ( Box& bx,
                                   MultiFab& div_mf,
                                   MultiFab& divc_mf,
                                   const MultiFab& weights,
                                   MFIter* mfi,
                                   const int icomp,
                                   const int ncomp,
                                   const EBCellFlagFab& flags_fab,
                                   const MultiFab* volfrac,
                                   Box& domain,
                                   const int cyclic_x,
                                   const int cyclic_y,
                                   const int cyclic_z,
                                   const Real* dx)
    {
        //
        // Check that grid is uniform
        //
        const Real tolerance = std::numeric_limits<Real>::epsilon();

#if (AMREX_SPACEDIM == 2)
        if (std::abs(dx[0] - dx[1]) > tolerance) 
            amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
#elif (AMREX_SPACEDIM == 3)
        if( (std::abs(dx[0] - dx[1]) > tolerance) or
            (std::abs(dx[0] - dx[2]) > tolerance) or
            (std::abs(dx[1] - dx[2]) > tolerance) )
            amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
#endif

        const amrex::Dim3 dom_low  = amrex::lbound(domain);
        const amrex::Dim3 dom_high = amrex::ubound(domain);

        //
        // Get array from arguments
        //
        Array4<Real> const& div  = div_mf.array(*mfi);
        Array4<Real> const& divc = divc_mf.array(*mfi);
        auto const&         wt   = weights.array(*mfi);
        auto const&        flags = flags_fab.array();
        auto const&        vfrac = volfrac->array(*mfi);

        const Box& grown1_bx = amrex::grow(bx,1);
        const Box& grown2_bx = amrex::grow(bx,2);

        //
        // Working arrays
        //
        FArrayBox  delm_fab(grown1_bx,ncomp);
        FArrayBox  optmp_fab(grown2_bx,ncomp);
        FArrayBox  mask_fab(grown2_bx);

        Array4<Real> const& optmp = optmp_fab.array();
        Array4<Real> const& mask  = mask_fab.array();
        Array4<Real> const& delm  = delm_fab.array();

        //
        // Array "mask" is used to sever the link to ghost cells when the BCs 
        // are not periodic
        // It is set to 1 when a cell can be used in computations, 0 otherwise
        //
        AMREX_FOR_3D(grown2_bx, i, j, k,
        {
            if(((not cyclic_x) and (i < dom_low.x or i > dom_high.x)) or
               ((not cyclic_y) and (j < dom_low.y or j > dom_high.y)) or
               ((not cyclic_z) and (k < dom_low.z or k > dom_high.z)))
                mask(i,j,k) = 0;
            else
                mask(i,j,k) = 1;
        });

        //
        // Init to zero tmp array
        //
        AMREX_FOR_4D(grown2_bx, ncomp, i, j, k, n,
        {
            optmp(i,j,k,n) = 0;
        });

        //
        // Step 2: compute delta M (mass gain or loss) on (lo-1,lo+1)
        //
        AMREX_FOR_4D(grown1_bx, ncomp, i, j, k, n,
        {
        if(flags(i,j,k).isSingleValued())
        {
            Real divnc(0.0);
            Real vtot(0.0);
            Real wted_frac(0.0);

            for(int ii(-1); ii <= 1; ii++)
                for(int jj(-1); jj <= 1; jj++)
                    for(int kk(-1); kk <= 1; kk++)
                        if( (ii != 0 or jj != 0 or kk != 0) and
                            (flags(i,j,k).isConnected({AMREX_D_DECL(ii,jj,kk)}) == 1))
                        {
                            wted_frac = vfrac(i+ii,j+jj,k+kk) * wt(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                            vtot   += wted_frac;
                            divnc  += wted_frac * divc(i+ii,j+jj,k+kk,n);
                        }

            divnc /= vtot;

            // We need to multiply divc by mask to make sure optmp is zero for cells
            // outside the domain for non-cyclic BCs
            optmp(i,j,k,n) =  (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n) * mask(i,j,k));
            delm(i,j,k,n)  = -(    vfrac(i,j,k)) * optmp(i,j,k,n);
        }
        else
        {
            delm(i,j,k,n) = 0;
        }
        });


        //
        // Step 3: redistribute excess/loss of mass
        //
        AMREX_FOR_4D(grown1_bx, ncomp, i, j, k, n,
        {
        if(flags(i,j,k).isSingleValued())
        {
            Real wtot(0.0);

            for(int ii(-1); ii <= 1; ii++)
                for(int jj(-1); jj <= 1; jj++)
                    for(int kk(-1); kk <= 1; kk++)
                        if((ii != 0 or jj != 0 or kk != 0) and
                            (flags(i,j,k).isConnected({AMREX_D_DECL(ii,jj,kk)}) == 1))
                        {
                            wtot += wt(i+ii,j+jj,k+kk) * vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                        }


            wtot = 1.0/wtot;

            for(int ii(-1); ii <= 1; ii++)
                for(int jj(-1); jj <= 1; jj++)
                    for(int kk(-1); kk <= 1; kk++)
                        if((ii != 0 or jj != 0 or kk != 0) and
                            (flags(i,j,k).isConnected({AMREX_D_DECL(ii,jj,kk)}) == 1))
                        {
#ifdef AMREX_USE_CUDA
                            Gpu::Atomic::Add(&optmp(i+ii,j+jj,k+kk,n),
                                             delm(i,j,k,n) * wtot * mask(i+ii,j+jj,k+kk) * wt(i+ii,j+jj,k+kk));
#else
                            optmp(i+ii,j+jj,k+kk,n) += delm(i,j,k,n) * wtot * mask(i+ii,j+jj,k+kk) * wt(i+ii,j+jj,k+kk);
#endif
                        }
        }
        });

        //
        // Resume the correct sign, AKA return the negative
        //
        AMREX_FOR_4D(bx, ncomp, i, j, k, n,
        {
            div(i,j,k,icomp+n) = divc(i,j,k,n) + optmp(i,j,k,n);
        });

        Gpu::synchronize();
    }

    // 
    // Do small cell redistribution on a MultiFab -- with a weighting function
    // 
    void single_level_weighted_redistribute ( int lev, MultiFab& div_tmp_in, MultiFab& div_out, const MultiFab& weights,
                                              int div_comp, int ncomp, const Vector<Geometry>& geom)
    {
        Box domain(geom[lev].Domain());

        const int cyclic_x = geom[0].isPeriodic(0) ? 1 : 0;
        const int cyclic_y = geom[0].isPeriodic(1) ? 1 : 0;
#if (AMREX_SPACEDIM == 2)
        const int cyclic_z = 0;
#elif (AMREX_SPACEDIM == 3)
        const int cyclic_z = geom[0].isPeriodic(2) ? 1 : 0;
#endif

        Real covered_val = 1.e40;

        int nghost = div_tmp_in.nGrow();

        EB_set_covered(div_tmp_in, div_comp, ncomp, nghost, covered_val);

        div_tmp_in.FillBoundary(geom[lev].periodicity());

        // Here we take care of both the regular and covered cases ... all we do below is the cut cell cases
        MultiFab::Copy(div_out, div_tmp_in, 0, div_comp, ncomp, div_out.nGrow());

        // Get EB geometric info
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(div_out.Factory());
        const MultiFab* volfrac   = &(ebfactory. getVolFrac());

        for (MFIter mfi(div_out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tilebox
            Box bx = mfi.tilebox ();
    
            // this is to check efficiently if this tile contains any eb stuff
            const EBFArrayBox&  div_fab = static_cast<EBFArrayBox const&>(div_out[mfi]);
            const EBCellFlagFab&  flags = div_fab.getEBCellFlagFab();
    
            if ( !(flags.getType(amrex::grow(bx,     0)) == FabType::covered) &&
                 !(flags.getType(amrex::grow(bx,nghost)) == FabType::regular) )
            {
                // Compute div(tau) with EB algorithm
                apply_eb_redistribution(bx, div_out, div_tmp_in, weights, &mfi,
                                               div_comp, ncomp, flags, volfrac, domain,
                                               cyclic_x, cyclic_y, cyclic_z,
                                               geom[lev].CellSize());

            }
        }
    }

    // 
    // Do small cell redistribution on a MultiFab -- without a weighting function
    // 
    void single_level_redistribute ( int lev, MultiFab& div_tmp_in, MultiFab& div_out,
                                     int div_comp, int ncomp, const Vector<Geometry>& geom)
    {
        // We create a weighting array to use inside the redistribution array
        MultiFab weights(div_out.boxArray(), div_out.DistributionMap(), 1, div_tmp_in.nGrow());
        weights.setVal(1.0);

        single_level_weighted_redistribute ( lev, div_tmp_in, div_out, weights, div_comp, ncomp, geom);
    }
#endif

} // end namespace
