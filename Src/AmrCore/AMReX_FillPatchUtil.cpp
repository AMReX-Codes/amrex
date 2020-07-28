#include <AMReX_FillPatchUtil.H>

#ifndef BL_NO_FORT
#include <AMReX_FillPatchUtil_F.H>
#endif

namespace amrex
{
#ifndef BL_NO_FORT
    // B fields are assumed to be on staggered grids.
    void InterpCrseFineBndryEMfield (InterpEM_t interp_type,
                                     const Array<MultiFab,AMREX_SPACEDIM>& crse,
                                     Array<MultiFab,AMREX_SPACEDIM>& fine,
                                     const Geometry& cgeom, const Geometry& fgeom,
                                     int ref_ratio)
    {
        InterpCrseFineBndryEMfield(interp_type,
                                   {{AMREX_D_DECL(&crse[0],&crse[1],&crse[2])}},
                                   {{AMREX_D_DECL(&fine[0],&fine[1],&fine[2])}},
                                   cgeom, fgeom, ref_ratio);
    }

    void InterpCrseFineBndryEMfield (InterpEM_t interp_type,
                                     const Array<MultiFab const*,AMREX_SPACEDIM>& crse,
                                     const Array<MultiFab*,AMREX_SPACEDIM>& fine,
                                     const Geometry& cgeom, const Geometry& fgeom,
                                     int ref_ratio)
    {
        // TODO gpu

        BL_ASSERT(ref_ratio == 2);

        const IntVect& ngrow = fine[0]->nGrowVect();
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
            BL_ASSERT(ngrow == fine[idim]->nGrowVect());
        }

        if (ngrow.max() == 0) return;

        bool include_periodic = true;
        bool include_physbndry = false;
        const auto& cfinfo = FabArrayBase::TheCFinfo(*fine[0], fgeom, ngrow,
                                                     include_periodic, include_physbndry);

        if (! cfinfo.ba_cfb.empty())
        {
            std::array<MultiFab, AMREX_SPACEDIM> cmf;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& fine_ba = fine[idim]->boxArray();
                BoxArray fba = cfinfo.ba_cfb;
                fba.convert(fine_ba.ixType());
                BoxArray cba = fba;
                cba.coarsen(ref_ratio);
                const DistributionMapping& dm = cfinfo.dm_cfb;

#ifdef AMREX_USE_EB
                amrex::Abort("InterpCrseFineBndryEMfield: EB is allowed");
#endif

                cmf[idim].define(cba, dm, 1, 1, MFInfo(), crse[0]->Factory());

                cmf[idim].copy(*crse[idim], 0, 0, 1, 0, 1, cgeom.periodicity());
            }

            const Real* dx = cgeom.CellSize();

            const int use_limiter = 0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                std::array<FArrayBox,AMREX_SPACEDIM> bfab;
                for (MFIter mfi(cmf[0]); mfi.isValid(); ++mfi)
                {
                    const int fi = cfinfo.fine_grid_idx[mfi.LocalIndex()];

                    Box ccbx = amrex::grow(fine[0]->boxArray()[fi], ngrow);
                    ccbx.enclosedCells();
                    ccbx.coarsen(ref_ratio).refine(ref_ratio);  // so that ccbx is coarsenable

                    const FArrayBox& cxfab = cmf[0][mfi];
                    bfab[0].resize(amrex::convert(ccbx,fine[0]->ixType()));
#if (AMREX_SPACEDIM > 1)
                    const FArrayBox& cyfab = cmf[1][mfi];
                    bfab[1].resize(amrex::convert(ccbx,fine[1]->ixType()));
#endif
#if (AMREX_SPACEDIM > 2)
                    const FArrayBox& czfab = cmf[2][mfi];
                    bfab[2].resize(amrex::convert(ccbx,fine[2]->ixType()));
#endif
                    // interpolate from cmf to fmf
                    if (interp_type == InterpB)
                    {
                        amrex_interp_div_free_bfield(BL_TO_FORTRAN_BOX(ccbx),
                                                     AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bfab[0]),
                                                                  BL_TO_FORTRAN_ANYD(bfab[1]),
                                                                  BL_TO_FORTRAN_ANYD(bfab[2])),
                                                     AMREX_D_DECL(BL_TO_FORTRAN_ANYD(cxfab),
                                                                  BL_TO_FORTRAN_ANYD(cyfab),
                                                                  BL_TO_FORTRAN_ANYD(czfab)),
                                                     dx, &ref_ratio, &use_limiter);
                    }
                    else if (interp_type == InterpE)
                    {
                        amrex_interp_efield(BL_TO_FORTRAN_BOX(ccbx),
                                            AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bfab[0]),
                                                         BL_TO_FORTRAN_ANYD(bfab[1]),
                                                         BL_TO_FORTRAN_ANYD(bfab[2])),
                                            AMREX_D_DECL(BL_TO_FORTRAN_ANYD(cxfab),
                                                         BL_TO_FORTRAN_ANYD(cyfab),
                                                         BL_TO_FORTRAN_ANYD(czfab)),
                                            &ref_ratio, &use_limiter);
                    }
                    else
                    {
                        amrex::Abort("InterpCrseFineBndryEMfield: unknown interp_type");
                    }

                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        const BoxArray& fine_ba = fine[idim]->boxArray();
                        const Box& fine_valid_box = fine_ba[fi];
                        Box b = bfab[idim].box();
                        b &= fine_valid_box;
                        const BoxList& diff = amrex::boxDiff(b, fine_valid_box); // skip valid cells
                        FArrayBox& fine_fab = (*fine[idim])[fi];
                        for (const auto& x : diff)
                        {
                            fine_fab.copy<RunOn::Host>(bfab[idim], x, 0, x, 0, 1);
                        }
                    }
                }
            }
        }
    }
#endif
}
