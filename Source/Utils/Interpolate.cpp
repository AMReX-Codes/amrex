#include "Interpolate.H"
#include <AMReX_FillPatchUtil_F.H>

namespace Interpolate
{
    using namespace amrex;

    std::unique_ptr<MultiFab>
    getInterpolatedScalar(
        const MultiFab& F_cp, const MultiFab& F_fp,
        const DistributionMapping& dm, const int r_ratio,
        const Real* /*dx*/, const int ngrow )
    {
        // Prepare the structure that will contain the returned fields
        std::unique_ptr<MultiFab> interpolated_F;
        interpolated_F.reset( new MultiFab(F_fp.boxArray(), dm, 1, ngrow) );
        interpolated_F->setVal(0.);

        // Loop through the boxes and interpolate the values from the _cp data
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox ffab; // Temporary array ; contains interpolated fields
            for (MFIter mfi(*interpolated_F); mfi.isValid(); ++mfi)
            {
                Box finebx = mfi.fabbox();
                finebx.coarsen(r_ratio).refine(r_ratio); // so that finebx is coarsenable

                const FArrayBox& cfab = (F_cp)[mfi];
                ffab.resize(finebx);

                // - Fully nodal
                if ( F_fp.is_nodal() ){
                    IntVect refinement_vector{AMREX_D_DECL(r_ratio, r_ratio, r_ratio)};
                    node_bilinear_interp.interp(cfab, 0, ffab, 0, 1,
                                                finebx, refinement_vector, {}, {}, {}, 0, 0, RunOn::Cpu);
                } else {
                    amrex::Abort("Unknown field staggering.");
                }

                // Add temporary array to the returned structure
                const Box& bx = (*interpolated_F)[mfi].box();
                (*interpolated_F)[mfi].plus<RunOn::Host>(ffab, bx, bx, 0, 0, 1);
            }
        }
        return interpolated_F;
    }

    std::array<std::unique_ptr<MultiFab>, 3>
    getInterpolatedVector(
        const MultiFab* Fx_cp,
        const MultiFab* Fy_cp,
        const MultiFab* Fz_cp,
        const MultiFab* Fx_fp,
        const MultiFab* Fy_fp,
        const MultiFab* Fz_fp,
        const DistributionMapping& dm, const int r_ratio,
        const Real* dx, const int ngrow )
    {

        // Prepare the structure that will contain the returned fields
        std::array<std::unique_ptr<MultiFab>, 3> interpolated_F;
        interpolated_F[0].reset( new MultiFab(Fx_fp->boxArray(), dm, 1, ngrow) );
        interpolated_F[1].reset( new MultiFab(Fy_fp->boxArray(), dm, 1, ngrow) );
        interpolated_F[2].reset( new MultiFab(Fz_fp->boxArray(), dm, 1, ngrow) );
        for (int i=0; i<3; i++) interpolated_F[i]->setVal(0.);

        // Loop through the boxes and interpolate the values from the _cp data
        const int use_limiter = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::array<FArrayBox,3> ffab; // Temporary array ; contains interpolated fields
            for (MFIter mfi(*interpolated_F[0]); mfi.isValid(); ++mfi)
            {
                Box ccbx = mfi.fabbox();
                ccbx.enclosedCells();
                ccbx.coarsen(r_ratio).refine(r_ratio); // so that ccbx is coarsenable

                const FArrayBox& cxfab = (*Fx_cp)[mfi];
                const FArrayBox& cyfab = (*Fy_cp)[mfi];
                const FArrayBox& czfab = (*Fz_cp)[mfi];
                ffab[0].resize(amrex::convert(ccbx,(*Fx_fp)[mfi].box().type()));
                ffab[1].resize(amrex::convert(ccbx,(*Fy_fp)[mfi].box().type()));
                ffab[2].resize(amrex::convert(ccbx,(*Fz_fp)[mfi].box().type()));

                // - Face centered, in the same way as B on a Yee grid
                if ( (*Fx_fp)[mfi].box().type() == IntVect{AMREX_D_DECL(1,0,0)} ){
#if (AMREX_SPACEDIM == 3)
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 BL_TO_FORTRAN_ANYD(ffab[0]),
                                                 BL_TO_FORTRAN_ANYD(ffab[1]),
                                                 BL_TO_FORTRAN_ANYD(ffab[2]),
                                                 BL_TO_FORTRAN_ANYD(cxfab),
                                                 BL_TO_FORTRAN_ANYD(cyfab),
                                                 BL_TO_FORTRAN_ANYD(czfab),
                                                 dx, &r_ratio, &use_limiter);
#else
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 BL_TO_FORTRAN_ANYD(ffab[0]),
                                                 BL_TO_FORTRAN_ANYD(ffab[2]),
                                                 BL_TO_FORTRAN_ANYD(cxfab),
                                                 BL_TO_FORTRAN_ANYD(czfab),
                                                 dx, &r_ratio, &use_limiter);
                    amrex_interp_cc_bfield(ccbx.loVect(), ccbx.hiVect(),
                                           BL_TO_FORTRAN_ANYD(ffab[1]),
                                           BL_TO_FORTRAN_ANYD(cyfab),
                                           &r_ratio, &use_limiter);
#endif
                    // - Edge centered, in the same way as E on a Yee grid
                } else if ( (*Fx_fp)[mfi].box().type() == IntVect{AMREX_D_DECL(0,1,1)} ){
#if (AMREX_SPACEDIM == 3)
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        BL_TO_FORTRAN_ANYD(ffab[0]),
                                        BL_TO_FORTRAN_ANYD(ffab[1]),
                                        BL_TO_FORTRAN_ANYD(ffab[2]),
                                        BL_TO_FORTRAN_ANYD(cxfab),
                                        BL_TO_FORTRAN_ANYD(cyfab),
                                        BL_TO_FORTRAN_ANYD(czfab),
                                        &r_ratio, &use_limiter);
#else
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        BL_TO_FORTRAN_ANYD(ffab[0]),
                                        BL_TO_FORTRAN_ANYD(ffab[2]),
                                        BL_TO_FORTRAN_ANYD(cxfab),
                                        BL_TO_FORTRAN_ANYD(czfab),
                                        &r_ratio,&use_limiter);
                    amrex_interp_nd_efield(ccbx.loVect(), ccbx.hiVect(),
                                           BL_TO_FORTRAN_ANYD(ffab[1]),
                                           BL_TO_FORTRAN_ANYD(cyfab),
                                           &r_ratio);
#endif
                } else {
                    amrex::Abort("Unknown field staggering.");
                }

                // Add temporary array to the returned structure
                for (int i = 0; i < 3; ++i) {
                    const Box& bx = (*interpolated_F[i])[mfi].box();
                    (*interpolated_F[i])[mfi].plus<RunOn::Host>(ffab[i], bx, bx, 0, 0, 1);
                }
            }
        }
        return interpolated_F;
    }
}
