#include <AMReX_BCRec.H>
#include <AMReX_EBFluxRegister.H>
#include <AMReX_YAFluxRegister.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_YAFluxRegister_K.H>
#include <AMReX_EBMultiFabUtil_3D_C.H>
#include <AMReX_EB_Redistribution.H>

namespace amrex {

void
amrex_flux_redistribute (
    const Box& bx,
    Array4<Real            > const& dqdt,
    Array4<Real       const> const& divc,
    Array4<Real       const> const& wt,
    Array4<Real       const> const& vfrac,
    Array4<EBCellFlag const> const& flag,
    int as_crse,
    Array4<Real            > const& rr_drho_crse,
    Array4<int        const> const& rr_flag_crse,
    int as_fine,
    Array4<Real            > const& dm_as_fine,
    Array4<int        const> const& levmsk,
    const Geometry& geom,
    bool use_wts_in_divnc,
    int level_mask_not_covered,
    int icomp, int ncomp, Real dt)
{
    //
    // Check that grid is uniform
    //
    const Real* dx = geom.CellSize();

#if (AMREX_SPACEDIM == 2)
    if (! amrex::almostEqual(dx[0], dx[1]))
#elif (AMREX_SPACEDIM == 3)
    if( ! amrex::almostEqual(dx[0],dx[1]) ||
        ! amrex::almostEqual(dx[1],dx[2]) )
#endif
    {
        amrex::Abort("apply_eb_redistribution(): grid spacing must be uniform");
    }

    const Box dbox1 = geom.growPeriodicDomain(1);
    const Box dbox2 = geom.growPeriodicDomain(2);

    const Box& grown1_bx = amrex::grow(bx,1);
    const Box& grown2_bx = amrex::grow(bx,2);

    Real reredistribution_threshold = amrex_eb_get_reredistribution_threshold();

    int bx_ilo = bx.smallEnd()[0];
    int bx_ihi = bx.bigEnd()[0];
    int bx_jlo = bx.smallEnd()[1];
    int bx_jhi = bx.bigEnd()[1];
#if (AMREX_SPACEDIM == 3)
    int bx_klo = bx.smallEnd()[2];
    int bx_khi = bx.bigEnd()[2];
#endif

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
        mask(i,j,k) = (dbox2.contains(IntVect(AMREX_D_DECL(i,j,k)))) ? 1.0 : 0.0;
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
    //Gives an optoion to use weights to compute divnc or not

    if (use_wts_in_divnc) {
        amrex::ParallelFor(grown1_bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isSingleValued())
            {
                Real vtot(0.);
                Real divnc(0.);
#if (AMREX_SPACEDIM == 2)
                int kk(0);
#else
                for (int kk = -1; kk <= 1; kk++) {
#endif
                for (int jj = -1; jj <= 1; jj++) {
                for (int ii = -1; ii <= 1; ii++) {
                    if ( (ii != 0 || jj != 0 || kk != 0) && flag(i,j,k).isConnected(ii,jj,kk) &&
                         dbox2.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))))
                    {
                        Real wted_frac = vfrac(i+ii,j+jj,k+kk) * wt(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                        vtot   += wted_frac;
                        divnc  += wted_frac * divc(i+ii,j+jj,k+kk,n);
                    }
                AMREX_D_TERM(},},})

                divnc /= vtot;

                // We need to multiply by mask to make sure optmp is zero for cells
                // outside the domain for non-cyclic BCs
                optmp(i,j,k,n) =  (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n)) * mask(i,j,k);
                delm(i,j,k,n)  = -(    vfrac(i,j,k)) * optmp(i,j,k,n);

            } else {
                delm(i,j,k,n) = 0.;
            }
        });
    } else {
        amrex::ParallelFor(grown1_bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (flag(i,j,k).isSingleValued())
            {
                Real vtot(0.);
                Real divnc(0.);
#if (AMREX_SPACEDIM == 2)
                int kk(0);
#else
                for (int kk = -1; kk <= 1; kk++) {
#endif
                for (int jj = -1; jj <= 1; jj++) {
                for (int ii = -1; ii <= 1; ii++) {
                    if ( (ii != 0 || jj != 0 || kk != 0) && flag(i,j,k).isConnected(ii,jj,kk) &&
                         dbox2.contains(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk))))
                    {
                        Real unwted_frac = vfrac(i+ii,j+jj,k+kk) * mask(i+ii,j+jj,k+kk);
                        vtot  += unwted_frac;
                        divnc += unwted_frac*divc(i+ii,j+jj,k+kk,n);
                    }
                AMREX_D_TERM(},},})

                divnc /= vtot;

                // We need to multiply by mask to make sure optmp is zero for cells
                // outside the domain for non-cyclic BCs
                optmp(i,j,k,n) =  (1 - vfrac(i,j,k)) * (divnc - divc(i,j,k,n)) * mask(i,j,k);
                delm(i,j,k,n)  = -(    vfrac(i,j,k)) * optmp(i,j,k,n);
            } else {
                delm(i,j,k,n) = 0.;
            }
        });
    }

    //
    // Step 3: redistribute excess/loss of mass
    //

    amrex::ParallelFor(grown1_bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        bool valid_dst_cell;
        if (flag(i,j,k).isSingleValued())
        {
            Real wtot = Real(0.);
#if (AMREX_SPACEDIM == 2)
            int kk(0);
#else
            for (int kk = -1; kk <= 1; kk++) {
#endif
            for (int jj = -1; jj <= 1; jj++) {
            for (int ii = -1; ii <= 1; ii++) {
                if ( (ii != 0 || jj != 0 || kk != 0) && flag(i,j,k).isConnected(ii,jj,kk) )
                {
                    wtot += vfrac(i+ii,j+jj,k+kk)*wt(i+ii,j+jj,k+kk)* mask(i+ii,j+jj,k+kk);
                }
            AMREX_D_TERM(},},})

#ifdef AMREX_USE_FLOAT
            wtot = Real(1.0)/(wtot + Real(1.e-30));
#else
            wtot = Real(1.0)/(wtot + Real(1.e-80));
#endif

            bool as_crse_crse_cell    = false;
            bool as_crse_covered_cell = false;

            if (as_crse)
            {
                bool inside =
#if (AMREX_SPACEDIM == 2)
                    ( (i >= bx_ilo) && (i <= bx_ihi) && (j >= bx_jlo) && (j <= bx_jhi) );
#else
                    ( (i >= bx_ilo) && (i <= bx_ihi) && (j >= bx_jlo) && (j <= bx_jhi) && (k >= bx_klo) && (k <= bx_khi) );
#endif
                as_crse_crse_cell    = inside && (rr_flag_crse(i,j,k) == amrex_yafluxreg_crse_fine_boundary_cell);
                as_crse_covered_cell = (rr_flag_crse(i,j,k) == amrex_yafluxreg_fine_cell);
            }

            bool as_fine_valid_cell = false;  // valid cells near box boundary
            bool as_fine_ghost_cell = false;  // ghost cells just outside valid region

            if (as_fine)
            {
                bool inside =
#if (AMREX_SPACEDIM == 2)
                    ( (i >= bx_ilo) && (i <= bx_ihi) && (j >= bx_jlo) && (j <= bx_jhi) );
#else
                    ( (i >= bx_ilo) && (i <= bx_ihi) && (j >= bx_jlo) && (j <= bx_jhi) && (k >= bx_klo) && (k <= bx_khi) );
#endif
                if (inside) { as_fine_valid_cell = true; }
                as_fine_ghost_cell = (levmsk(i,j,k) == level_mask_not_covered); // not covered by other grids at this level
            }

#if (AMREX_SPACEDIM == 2)
            kk = 0;
#else
            for (int kk = -1; kk <= 1; kk++) {
#endif
            for (int jj = -1; jj <= 1; jj++) {
            for (int ii = -1; ii <= 1; ii++) {
                if ( (ii != 0 || jj != 0 || kk != 0) && flag(i,j,k).isConnected(ii,jj,kk) )
                {
                    int iii = i + ii;
                    int jjj = j + jj;
                    int kkk = k + kk;

                    Real drho = delm(i,j,k,n)*wtot*wt(iii,jjj,kkk)* mask(iii,jjj,kkk) ;
                    Gpu::Atomic::Add(&optmp(iii,jjj,kkk,n), drho);

                    valid_dst_cell = ( (iii >= bx_ilo) && (iii <= bx_ihi) &&
                                       (jjj >= bx_jlo) && (jjj <= bx_jhi) );
#if (AMREX_SPACEDIM == 3)
                    valid_dst_cell &= ( (kkk >= bx_klo) && (kkk <= bx_khi) );
#endif

                    if (as_crse_crse_cell)
                    {
                        if ( (rr_flag_crse(iii,jjj,kkk) == amrex_yafluxreg_fine_cell) &&
                             (vfrac(i,j,k) > reredistribution_threshold) )
                        {
                            Gpu::Atomic::Add(&rr_drho_crse(i,j,k,n),
                                             dt*drho*(vfrac(iii,jjj,kkk)/vfrac(i,j,k)));
                        }
                    }

                    if (as_crse_covered_cell && valid_dst_cell)
                    {
                        if ( (rr_flag_crse(iii,jjj,kkk) == amrex_yafluxreg_crse_fine_boundary_cell) &&
                             (vfrac(iii,jjj,kkk) > reredistribution_threshold) )
                        {
                            // recipient is a crse/fine boundary cell
                            Gpu::Atomic::Add(&rr_drho_crse(iii,jjj,kkk,n), -dt*drho);
                        }
                    }

                    if (as_fine_valid_cell && !valid_dst_cell && dbox1.contains(IntVect(AMREX_D_DECL(iii,jjj,kkk))))
                    {
                        Gpu::Atomic::Add(&dm_as_fine(iii,jjj,kkk,n), dt*drho*vfrac(iii,jjj,kkk));
                    }

                    if (as_fine_ghost_cell && valid_dst_cell && dbox1.contains(IntVect(AMREX_D_DECL(i,j,k))))
                    {
                        Gpu::Atomic::Add(&dm_as_fine(i,j,k,n), -dt*drho*vfrac(iii,jjj,kkk));
                    }
                } // isConnected
            AMREX_D_TERM(},},})
        } // isSingleValued
    });

    amrex::ParallelFor(bx, ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (!flag(i,j,k).isCovered()) {
            dqdt(i,j,k,icomp+n) = divc(i,j,k,n) + optmp(i,j,k,n);
        }
    });


    Gpu::streamSynchronize(); // because of FArrayBoxes defined in this function
} // end amrex_flux_redistribute

//
// Do small cell redistribution on one FAB with the Array4's already passed in
//
void
apply_flux_redistribution ( const Box& bx,
                            Array4<Real      > const& div,
                            Array4<Real const> const& divc,
                            Array4<Real const> const& wt,
                            int icomp, int ncomp,
                            Array4<EBCellFlag const> const& flag_arr,
                            Array4<Real       const> const& vfrac,
                            const Geometry & geom,
                            bool use_wts_in_divnc)
{
    int as_crse = 0;
    int as_fine = 0;
    Real dummy_dt  = 0.0;
    Array4<Real> dummy_drho_crse  = Array4<Real>();
    Array4<Real> dummy_dm_as_fine = Array4<Real>();
    Array4<int > dummy_levmsk     = Array4<int>();
    Array4<int > dummy_flag_crse  = Array4<int>();
    int dummy_not_covered = 0;
    amrex_flux_redistribute (bx, div, divc, wt, vfrac, flag_arr,
                             as_crse, dummy_drho_crse , dummy_flag_crse,
                             as_fine, dummy_dm_as_fine, dummy_levmsk,
                             geom, use_wts_in_divnc, dummy_not_covered,
                             icomp, ncomp, dummy_dt);
}
} // end namespace
