/**
 * \file AMReX_EB_StateRedistribute.cpp
 */

#include <AMReX_EB_Redistribution.H>
#include <AMReX_EB_StateRedistSlopeLimiter_K.H>
#include <AMReX_EB_Slopes_K.H>
#include <AMReX_YAFluxRegister_K.H>

// Neighbor index lookup tables in constant memory for efficient GPU access
// These map neighbor indices to (i,j,k) offsets
#if (AMREX_SPACEDIM == 2)
namespace {
    AMREX_GPU_CONSTANT int imap[9] = {0,-1, 0, 1,-1, 1,-1, 0, 1};
    AMREX_GPU_CONSTANT int jmap[9] = {0,-1,-1,-1, 0, 0, 1, 1, 1};
    AMREX_GPU_CONSTANT int kmap[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
}
#elif (AMREX_SPACEDIM == 3)
namespace {
    AMREX_GPU_CONSTANT int imap[27] = {0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    AMREX_GPU_CONSTANT int jmap[27] = {0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    AMREX_GPU_CONSTANT int kmap[27] = {0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
}
#endif

namespace amrex {

void
MLStateRedistribute ( Box const& bx, int ncomp,
                      Array4<Real> const& U_out,
                      Array4<Real> const& U_in,
                      Array4<EBCellFlag const> const& flag,
                      Array4<Real const> const& vfrac,
                      AMREX_D_DECL(Array4<Real const> const& fcx,
                                   Array4<Real const> const& fcy,
                                   Array4<Real const> const& fcz),
                      Array4<Real const> const& ccent,
                      amrex::BCRec  const* d_bcrec_ptr,
                      Array4< int const> const& itracker,
                      Array4<Real const> const& nrs,
                      Array4<Real const> const& alpha,
                      Array4<Real const> const& nbhd_vol,
                      Array4<Real const> const& cent_hat,
                      Geometry const& lev_geom,
                      // Starting ML stuff
                      int as_crse,
                      Array4<Real            > const& drho_as_crse,
                      Array4<int        const> const& flag_as_crse,
                      int as_fine,
                      Array4<Real            > const& dm_as_fine,
                      Array4<int        const> const& levmsk,
                      int is_ghost_cell,
                      Real fac_for_deltaR,
                      int max_order)
{
    BL_PROFILE_VAR("MLStateRedist", MLStateRedist);
    // Note that itracker has {4 in 2D, 8 in 3D} components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    //
    // In 2D, we identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->
    //
    // In 3D, We identify the cells in the remaining three components with the following ordering
    //
    //    at k-1   |   at k  |   at k+1
    //
    // ^  15 16 17 |  6 7 8  |  24 25 26
    // |  12 13 14 |  4   5  |  21 22 23
    // j  9  10 11 |  1 2 3  |  18 19 20
    //   i --->
    //
    // Note the first component of each of these arrays should never be used
    // imap, jmap, kmap are now in constant memory (see file scope declarations above)
    //

    const Box domain = lev_geom.Domain();
    const int domain_ilo = domain.smallEnd(0);
    const int domain_ihi = domain.bigEnd(0);
    const int domain_jlo = domain.smallEnd(1);
    const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);
#endif

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2););

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);
    Box const& bxg3 = amrex::grow(bx,3);

    // Qhat is formed for every cell in bxg2, and the neighborhood of a cell in bxg2 can
    //    reach one cell further out, so we must grow by 3 -- not 2 -- in the periodic
    //    directions in order for the neighbors counted in nbhd_vol (see MakeStateRedistUtils,
    //    which grows by 5) to also contribute their state to Qhat.  U_in and vfrac are
    //    available on bxg3, so the data we need is there.
    Box domain_per_grown = domain;
    if (is_periodic_x) { domain_per_grown.grow(0,3); }
    if (is_periodic_y) { domain_per_grown.grow(1,3); }
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) { domain_per_grown.grow(2,3); }
#endif

    // Solution at the centroid of my nbhd
    FArrayBox    Qhat_fab (bxg3,ncomp,The_Async_Arena());
    Array4<Real> Qhat = Qhat_fab.array();

#if (AMREX_SPACEDIM == 2)
    FArrayBox qtracker_fab (bxg3,4,The_Async_Arena());
#elif (AMREX_SPACEDIM == 3)
    FArrayBox qtracker_fab (bxg3,8,The_Async_Arena());
#endif
    Array4<Real> qt = qtracker_fab.array();

    // Raw slopes, limiting, and redistribution updates use separate kernels so
    // their live values do not contribute to one another's register pressure.
    // The limiter overwrites the raw slopes in place before the scatter kernel.
    // The final component records which cells satisfy the limiter predicate so
    // the limiter kernel does not need to capture vfrac and itracker.
    FArrayBox slopes_eb_fab (bxg1,AMREX_SPACEDIM+1,The_Async_Arena());
    Array4<Real> slopes_eb_arr = slopes_eb_fab.array();

    // We do not initialize dm_as_fine here because we also add to cells outside
    // bx, which for a tiled MFIter belong to another tile. Zeroing it is the
    // caller's responsibility.

    for (int n = 0; n < ncomp; n++)
    {
        // Define Qhat (from Berger and Guliani)
        // Here we initialize Qhat to equal U_in on all cells in bxg3 so that
        //      in the event we need to use Qhat 3 cells out from the bx limits
        //      in a modified slope computation, we have a value of Qhat to use.
        //      But we only modify Qhat inside bxg2
        amrex::ParallelFor(bxg3,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (vfrac(i,j,k) > 0.0 && bxg2.contains(IntVect(AMREX_D_DECL(i,j,k)))
                                   && domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))))
            {
                if (itracker(i,j,k,0) == 0)
                {
                    qt(i,j,k,0  ) = alpha(i,j,k,0) * U_in(i,j,k,n) * vfrac(i,j,k) / nbhd_vol(i,j,k);
                    Qhat(i,j,k,n) = qt(i,j,k,0);
                }
                else
                {
                    Qhat(i,j,k,n) = 0.;

                    // This loops over (i,j,k) and the neighbors of (i,j,k)
                    for (int i_nbor = 0; i_nbor <= itracker(i,j,k,0); i_nbor++)
                    {
                        int r = i; int s = j; int t = k;
                        Real fac = alpha(i,j,k,0);
                        if (i_nbor > 0) {
                            r = i+imap[itracker(i,j,k,i_nbor)];
                            s = j+jmap[itracker(i,j,k,i_nbor)];
                            t = k+kmap[itracker(i,j,k,i_nbor)];
                            fac = alpha(i,j,k,1) / nrs(r,s,t);
                        }

                        if (domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))))
                        {
                            qt(i,j,k,i_nbor) = fac * U_in(r,s,t,n) * vfrac(r,s,t) / nbhd_vol(i,j,k);
                            Qhat(i,j,k,n) += qt(i,j,k,i_nbor);
                        }
                    }
                }
            } else {
                Qhat(i,j,k,n) = U_in(i,j,k,n);
            }
        });

       //
       // ****************************************************************************************
       //
        amrex::ParallelFor(bxg1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (vfrac(i,j,k) > 0.0 && itracker(i,j,k,0) == 0)
            {
                if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    U_out(i,j,k,n) += alpha(i,j,k,0) * qt(i,j,k,0);
                }
            }
        });

        // Step 1: Combined stencil and boundary classification kernel
        // Classify cells as inactive (0), interior (1-8), or boundary (11-18)
        // The mask encodes the grown-axis requirements in bits 0-2
        // Encoding: 0=inactive, 1-8=interior (mask+1), 11-18=boundary (mask+11)
        {
        amrex::Gpu::streamSynchronize();
        BL_PROFILE("SK:Kernel2-classify");
        amrex::ParallelFor(bxg1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool const has_redist_nbors = vfrac(i,j,k) > 0.0 && itracker(i,j,k,0) > 0;

            if (!has_redist_nbors) {
                slopes_eb_arr(i,j,k,AMREX_SPACEDIM) = 0.0_rt;
                return;
            }

            // Do we have enough extent in each coordinate direction to use the 3x3x3 stencil
            //    or do we need to enlarge it?
            AMREX_D_TERM(Real x_max = -Real(1.e30); Real x_min = Real(1.e30);,
                         Real y_max = -Real(1.e30); Real y_min = Real(1.e30);,
                         Real z_max = -Real(1.e30); Real z_min = Real(1.e30););

            constexpr Real slope_stencil_min_width = Real(0.5);
#if (AMREX_SPACEDIM == 2)
            int kkk = 0;
#elif (AMREX_SPACEDIM == 3)
            // Prevent loop unrolling to reduce register pressure
#ifdef AMREX_USE_CUDA
#pragma unroll 1
#endif
            for(int kkk(-1); kkk<=1; kkk++) {
#endif
#ifdef AMREX_USE_CUDA
#pragma unroll 1
#endif
            for(int jjj(-1); jjj<=1; jjj++) {
#ifdef AMREX_USE_CUDA
#pragma unroll 1
#endif
            for(int iii(-1); iii<=1; iii++) {
                if (flag(i,j,k).isConnected(iii,jjj,kkk))
                {
                    int rr = i+iii; int ss = j+jjj; int tt = k+kkk;

                    x_max = amrex::max(x_max, cent_hat(rr,ss,tt,0)+static_cast<Real>(iii));
                    x_min = amrex::min(x_min, cent_hat(rr,ss,tt,0)+static_cast<Real>(iii));
                    y_max = amrex::max(y_max, cent_hat(rr,ss,tt,1)+static_cast<Real>(jjj));
                    y_min = amrex::min(y_min, cent_hat(rr,ss,tt,1)+static_cast<Real>(jjj));
#if (AMREX_SPACEDIM == 3)
                    z_max = amrex::max(z_max, cent_hat(rr,ss,tt,2)+static_cast<Real>(kkk));
                    z_min = amrex::min(z_min, cent_hat(rr,ss,tt,2)+static_cast<Real>(kkk));
#endif
                }
            AMREX_D_TERM(},},})

            // Build grown-axis mask as an int
            int mask = 0;
            if ( (x_max-x_min) < slope_stencil_min_width ) { mask |= 1; }
            if ( (y_max-y_min) < slope_stencil_min_width ) { mask |= 2; }
#if (AMREX_SPACEDIM == 3)
            if ( (z_max-z_min) < slope_stencil_min_width ) { mask |= 4; }
#endif

            // Base code: mask+1 gives values 1-8
            int code = mask + 1;

            // Check if this cell needs boundary stencil (reproduces needs_bdry_stencil test)
            // Compute boundary condition flags on device
            auto const& bc = d_bcrec_ptr[n];
            bool const extdir_ilo = (bc.lo(0) == amrex::BCType::ext_dir ||
                                     bc.lo(0) == amrex::BCType::hoextrap);
            bool const extdir_ihi = (bc.hi(0) == amrex::BCType::ext_dir ||
                                     bc.hi(0) == amrex::BCType::hoextrap);
            bool const extdir_jlo = (bc.lo(1) == amrex::BCType::ext_dir ||
                                     bc.lo(1) == amrex::BCType::hoextrap);
            bool const extdir_jhi = (bc.hi(1) == amrex::BCType::ext_dir ||
                                     bc.hi(1) == amrex::BCType::hoextrap);
#if (AMREX_SPACEDIM == 3)
            bool const extdir_klo = (bc.lo(2) == amrex::BCType::ext_dir ||
                                     bc.lo(2) == amrex::BCType::hoextrap);
            bool const extdir_khi = (bc.hi(2) == amrex::BCType::ext_dir ||
                                     bc.hi(2) == amrex::BCType::hoextrap);
#endif

            // Note: includes cells one in from the boundary because 4th-order stencils reach domlo+1/domhi-1
            bool const needs_bdry = (extdir_ilo && i <= domain_ilo+1) || (extdir_ihi && i >= domain_ihi-1) ||
                                    (extdir_jlo && j <= domain_jlo+1) || (extdir_jhi && j >= domain_jhi-1)
#if (AMREX_SPACEDIM == 3)
                                 || (extdir_klo && k <= domain_klo+1) || (extdir_khi && k >= domain_khi-1)
#endif
                                    ;

            // Store classification: interior cells (1-8) or boundary cells (11-18)
            if (needs_bdry) {
                code += 10;  // Remap to boundary range: 1-8 -> 11-18
            }
            slopes_eb_arr(i,j,k,AMREX_SPACEDIM) = Real(code);
        });
        amrex::Gpu::streamSynchronize();
        } // Kernel2-classify

        // Step 2b: Interior slope computation kernel (no boundary logic)
        {
        amrex::Gpu::streamSynchronize();
        BL_PROFILE("SK:Kernel2-slope-interior");
        amrex::ParallelFor(bxg1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const code_real = slopes_eb_arr(i,j,k,AMREX_SPACEDIM);
            // Process only interior cells: code in range [1, 8]
            if (code_real < 1.0_rt || code_real > 8.0_rt) { return; }

            int const code = static_cast<int>(code_real);
            int const mask = code - 1;

            // Interior path: call the boundary-free slope functions directly
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes_eb;
            if (mask == 0) {
                // Standard 3x3x3 stencil (nx=ny=nz=1)
                slopes_eb = amrex_calc_slopes_eb(i,j,k,n,Qhat,cent_hat,vfrac,flag,max_order);
            } else {
                // Grown stencil path
                int const nx = (mask & 1) ? 2 : 1;
                int const ny = (mask & 2) ? 2 : 1;
#if (AMREX_SPACEDIM == 3)
                int const nz = (mask & 4) ? 2 : 1;
#else
                int const nz = 1;
#endif
                slopes_eb = amrex_calc_slopes_eb_grown(i,j,k,n,AMREX_D_DECL(nx,ny,nz),
                                                       Qhat,cent_hat,vfrac,flag,max_order);
            }

            AMREX_D_TERM(slopes_eb_arr(i,j,k,0) = slopes_eb[0];,
                         slopes_eb_arr(i,j,k,1) = slopes_eb[1];,
                         slopes_eb_arr(i,j,k,2) = slopes_eb[2];);
        });
        amrex::Gpu::streamSynchronize();
        } // Kernel2-slope-interior

        // Step 2c: Boundary slope computation kernel (with boundary logic)
        {
        amrex::Gpu::streamSynchronize();
        BL_PROFILE("SK:Kernel2-slope-boundary");

        amrex::ParallelFor(bxg1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const code_real = slopes_eb_arr(i,j,k,AMREX_SPACEDIM);
            // Process only boundary cells: code in range [11, 18]
            if (code_real < 11.0_rt || code_real > 18.0_rt) { return; }

            // Decode: boundary codes [11,18] -> masks [0,7]
            int const code = static_cast<int>(code_real) - 10;
            int const mask = code - 1;

            // Compute boundary condition flags on device
            auto const& bc = d_bcrec_ptr[n];
            bool const extdir_ilo = (bc.lo(0) == amrex::BCType::ext_dir ||
                                     bc.lo(0) == amrex::BCType::hoextrap);
            bool const extdir_ihi = (bc.hi(0) == amrex::BCType::ext_dir ||
                                     bc.hi(0) == amrex::BCType::hoextrap);
            bool const extdir_jlo = (bc.lo(1) == amrex::BCType::ext_dir ||
                                     bc.lo(1) == amrex::BCType::hoextrap);
            bool const extdir_jhi = (bc.hi(1) == amrex::BCType::ext_dir ||
                                     bc.hi(1) == amrex::BCType::hoextrap);
#if (AMREX_SPACEDIM == 3)
            bool const extdir_klo = (bc.lo(2) == amrex::BCType::ext_dir ||
                                     bc.lo(2) == amrex::BCType::hoextrap);
            bool const extdir_khi = (bc.hi(2) == amrex::BCType::ext_dir ||
                                     bc.hi(2) == amrex::BCType::hoextrap);
#endif

            // Boundary path: call functions with full boundary logic
            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes_eb;
            if (mask == 0) {
                // Standard 3x3x3 stencil (nx=ny=nz=1)
                slopes_eb = amrex_calc_slopes_extdir_eb(
                                            i,j,k,n,Qhat,cent_hat,vfrac,
                                            AMREX_D_DECL(fcx,fcy,fcz),flag,
                                            AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                                            AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                                            AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                            AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                            max_order);
            } else {
                // Grown stencil path
                int const nx = (mask & 1) ? 2 : 1;
                int const ny = (mask & 2) ? 2 : 1;
#if (AMREX_SPACEDIM == 3)
                int const nz = (mask & 4) ? 2 : 1;
#else
                int const nz = 1;
#endif
                slopes_eb = amrex_calc_slopes_extdir_eb_grown(
                                            i,j,k,n,AMREX_D_DECL(nx,ny,nz),
                                            Qhat,cent_hat,vfrac,
                                            AMREX_D_DECL(fcx,fcy,fcz),flag,
                                            AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                                            AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                                            AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                            AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                            max_order);
            }

            AMREX_D_TERM(slopes_eb_arr(i,j,k,0) = slopes_eb[0];,
                         slopes_eb_arr(i,j,k,1) = slopes_eb[1];,
                         slopes_eb_arr(i,j,k,2) = slopes_eb[2];);
        });
        amrex::Gpu::streamSynchronize();
        } // Kernel2-slope-boundary

        {
        amrex::Gpu::streamSynchronize();
        BL_PROFILE("SK:Kernel2-limiter");
        amrex::ParallelFor(bxg1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (slopes_eb_arr(i,j,k,AMREX_SPACEDIM) != 0.0_rt)
            {
                AMREX_D_TERM(Real slope_x = slopes_eb_arr(i,j,k,0);,
                             Real slope_y = slopes_eb_arr(i,j,k,1);,
                             Real slope_z = slopes_eb_arr(i,j,k,2););

                // We do the limiting separately because this limiter limits the slope based on the values
                //    extrapolated to the cell centroid (cent_hat) locations - unlike the limiter in amrex
                //    which bases the limiting on values extrapolated to the face centroids.
                amrex_apply_centroid_limiter(i,j,k,n,Qhat,flag,
                                             AMREX_D_DECL(slope_x,slope_y,slope_z),
                                             cent_hat);

                AMREX_D_TERM(slopes_eb_arr(i,j,k,0) = slope_x;,
                             slopes_eb_arr(i,j,k,1) = slope_y;,
                             slopes_eb_arr(i,j,k,2) = slope_z;);
            } // activity marker
        });
        amrex::Gpu::streamSynchronize();
        } // Kernel2-limiter


        {
        amrex::Gpu::streamSynchronize();
        BL_PROFILE("SK:Kernel2-scatter");        
        amrex::ParallelFor(bxg1,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (vfrac(i,j,k) > 0.0 && itracker(i,j,k,0) > 0)
            {
                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> lim_slope{
                    AMREX_D_DECL(slopes_eb_arr(i,j,k,0),
                                 slopes_eb_arr(i,j,k,1),
                                 slopes_eb_arr(i,j,k,2))};

                // Now loop over neighbors using the precomputed slopes
                for (int i_nbor = 0; i_nbor <= itracker(i,j,k,0); i_nbor++)
                {
                    int r = i; int s = j; int t = k;
                    Real fac = alpha(i,j,k,0) * nrs(i,j,k);
                    if (i_nbor > 0) {
                        r += imap[itracker(i,j,k,i_nbor)];
                        s += jmap[itracker(i,j,k,i_nbor)];
                        t += kmap[itracker(i,j,k,i_nbor)];
                        fac = alpha(i,j,k,1);
                    }

                    if (domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))))
                    {

                        for (int r_nbor = 0; r_nbor <= itracker(i,j,k,0); r_nbor++)
                        {
                            //
                            // Identify which U_in(ii,jj,kk) is contributing to this U_out(r,s,t) via Qhat(i,j,k)
                            //
                            int ii = i; int jj = j; int kk = k;
                            Real fac2 = alpha(i,j,k,0);
                            if (r_nbor > 0) {
                                ii += imap[itracker(i,j,k,r_nbor)];
                                jj += jmap[itracker(i,j,k,r_nbor)];
                                kk += kmap[itracker(i,j,k,r_nbor)];
                                fac2 = alpha(i,j,k,1) / nrs(ii,jj,kk);
                            }

                            //
                            // Scale the slopes for each qt
                            //
                            Real q_over_Q = fac2*vfrac(ii,jj,kk)/nbhd_vol(i,j,k);

                            // Skip if neighbor (ii,jj,kk) is outside domain_per_grown
                            // This matches the write condition in the first loop
                            if (!domain_per_grown.contains(IntVect(AMREX_D_DECL(ii,jj,kk)))) {
                                 continue;
                            }

                            Real update = qt(i,j,k,r_nbor);
                            AMREX_D_TERM(update += q_over_Q * lim_slope[0] *
                                                   (ccent(r,s,t,0)-cent_hat(i,j,k,0) + static_cast<Real>(r-i));,
                                         update += q_over_Q * lim_slope[1] *
                                                   (ccent(r,s,t,1)-cent_hat(i,j,k,1) + static_cast<Real>(s-j));,
                                         update += q_over_Q * lim_slope[2] *
                                                   (ccent(r,s,t,2)-cent_hat(i,j,k,2) + static_cast<Real>(t-k)););

                            //
                            // Update U_out(r,s,t) with qt(i,j,k,r_nbor)
                            //
                            if (bx.contains(IntVect(AMREX_D_DECL(r,s,t)))) {
                                amrex::HostDevice::Atomic::Add(&U_out(r,s,t,n),fac*update/nrs(r,s,t));
                            }

                            if (as_crse) {

                               // Covered (by fine) to uncovered (by fine)
                               if ( bx.contains(IntVect(AMREX_D_DECL(r,s,t))) &&
                                    flag_as_crse( r, s, t) == amrex_yafluxreg_crse_fine_boundary_cell &&
                                    flag_as_crse(ii,jj,kk) == amrex_yafluxreg_fine_cell )
                               {
                                   amrex::HostDevice::Atomic::Add(&drho_as_crse(r,s,t,n),
                                                                  -fac*update/nrs(r,s,t) * fac_for_deltaR);
                               }

                               // Uncovered (by fine) to covered (by fine)
                               if ( bx.contains(IntVect(AMREX_D_DECL(ii,jj,kk))) &&
                                    flag_as_crse( r, s, t) == amrex_yafluxreg_fine_cell  &&
                                    flag_as_crse(ii,jj,kk) == amrex_yafluxreg_crse_fine_boundary_cell )
                               {
                                   amrex::HostDevice::Atomic::Add(&drho_as_crse(ii,jj,kk,n),
                                                                  fac * update / nrs(r,s,t) *
                                                                  (vfrac(r,s,t) / vfrac(ii,jj,kk)) * fac_for_deltaR);
                               }
                            } // as_crse

                            if (as_fine) {

                               // Ghost (ii,jj,kk) to valid (r,s,t)
                               if (levmsk(ii,jj,kk) == is_ghost_cell && bx.contains(IntVect(AMREX_D_DECL(r,s,t)))) {

                                   amrex::HostDevice::Atomic::Add(&dm_as_fine(ii,jj,kk,n),
                                                                  -fac*update/nrs(r,s,t) * vfrac(r,s,t) * fac_for_deltaR);
                               }

                               // Valid (ii,jj,kk) to ghost (r,s,t)
                               if (bx.contains(IntVect(AMREX_D_DECL(ii,jj,kk))) && levmsk(r,s,t) == is_ghost_cell) {

                                   amrex::HostDevice::Atomic::Add(&dm_as_fine(r,s,t,n),
                                                                  fac*update/nrs(r,s,t) * vfrac(r,s,t) * fac_for_deltaR);
                               }
                            } // as_fine

                        } // r_nbor
                    } // bx contains
                } // i_nbor
            } // vfrac > 0 && itracker > 0

        });
        amrex::Gpu::streamSynchronize();
        } // Kernel2-scatter
    } // loop over n = 0, ncomp-1

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (flag(i,j,k).isCovered())
        {
            U_out(i,j,k,n) = Real(1.e30);
        }
    });
}

void
StateRedistribute ( Box const& bx, int ncomp,
                    Array4<Real> const& U_out,
                    Array4<Real> const& U_in,
                    Array4<EBCellFlag const> const& flag,
                    Array4<Real const> const& vfrac,
                    AMREX_D_DECL(Array4<Real const> const& fcx,
                                 Array4<Real const> const& fcy,
                                 Array4<Real const> const& fcz),
                    Array4<Real const> const& ccent,
                    amrex::BCRec  const* d_bcrec_ptr,
                    Array4< int const> const& itracker,
                    Array4<Real const> const& nrs,
                    Array4<Real const> const& alpha,
                    Array4<Real const> const& nbhd_vol,
                    Array4<Real const> const& cent_hat,
                    Geometry const& lev_geom,
                    int max_order)
{
    // Note that itracker has {4 in 2D, 8 in 3D} components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    //
    // In 2D, we identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->
    //
    // In 3D, We identify the cells in the remaining three components with the following ordering
    //
    //    at k-1   |   at k  |   at k+1
    //
    // ^  15 16 17 |  6 7 8  |  24 25 26
    // |  12 13 14 |  4   5  |  21 22 23
    // j  9  10 11 |  1 2 3  |  18 19 20
    //   i --->
    //
    // Note the first component of each of these arrays should never be used
    // imap, jmap, kmap are now in constant memory (see file scope declarations above)
    //

    const Box domain = lev_geom.Domain();
    const int domain_ilo = domain.smallEnd(0);
    const int domain_ihi = domain.bigEnd(0);
    const int domain_jlo = domain.smallEnd(1);
    const int domain_jhi = domain.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int domain_klo = domain.smallEnd(2);
    const int domain_khi = domain.bigEnd(2);
#endif

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2););

    Box const& bxg1 = amrex::grow(bx,1);
    Box const& bxg2 = amrex::grow(bx,2);
    Box const& bxg3 = amrex::grow(bx,3);

    // Qhat is formed for every cell in bxg2, and the neighborhood of a cell in bxg2 can
    //    reach one cell further out, so we must grow by 3 -- not 2 -- in the periodic
    //    directions in order for the neighbors counted in nbhd_vol (see MakeStateRedistUtils,
    //    which grows by 5) to also contribute their state to Qhat.  U_in and vfrac are
    //    available on bxg3, so the data we need is there.
    Box domain_per_grown = domain;
    if (is_periodic_x) { domain_per_grown.grow(0,3); }
    if (is_periodic_y) { domain_per_grown.grow(1,3); }
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) { domain_per_grown.grow(2,3); }
#endif

    // Solution at the centroid of my nbhd
    FArrayBox    Qhat_fab (bxg3,ncomp,The_Async_Arena());
    Array4<Real> Qhat = Qhat_fab.array();

    // Define Qhat (from Berger and Guliani)
    // Here we initialize Qhat to equal U_in on all cells in bxg3 so that
    //      in the event we need to use Qhat 3 cells out from the bx limits
    //      in a modified slope computation, we have a value of Qhat to use.
    //      But we only modify Qhat inside bxg2
    amrex::ParallelFor(bxg3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        for (int n = 0; n < ncomp; n++) {
            Qhat(i,j,k,n) = U_in(i,j,k,n);
        }

        if (vfrac(i,j,k) > 0.0 && bxg2.contains(IntVect(AMREX_D_DECL(i,j,k)))
                               && domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k)))) {

            // Start with U_in(i,j,k) itself
            for (int n = 0; n < ncomp; n++) {
                Qhat(i,j,k,n) = U_in(i,j,k,n) * alpha(i,j,k,0) * vfrac(i,j,k);
            }

            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int r = i+imap[itracker(i,j,k,i_nbor)];
                int s = j+jmap[itracker(i,j,k,i_nbor)];
                int t = k+kmap[itracker(i,j,k,i_nbor)];

                if (domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))))
                {
                    for (int n = 0; n < ncomp; n++) {
                        Qhat(i,j,k,n) += U_in(r,s,t,n) * alpha(i,j,k,1) * vfrac(r,s,t) / nrs(r,s,t);
                    }
                }
            }
            for (int n = 0; n < ncomp; n++) {
                Qhat(i,j,k,n) /= nbhd_vol(i,j,k);
            }
        }
    });

    amrex::ParallelFor(bxg1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac(i,j,k) > 0.0)
        {
            int num_nbors = itracker(i,j,k,0);

            if (itracker(i,j,k,0) == 0)
            {
                if (bx.contains(IntVect(AMREX_D_DECL(i,j,k))))
                {
                    for (int n = 0; n < ncomp; n++) {
                        amrex::Gpu::Atomic::Add(&U_out(i,j,k,n),alpha(i,j,k,0)*nrs(i,j,k)*Qhat(i,j,k,n));
                    }
                }

            } else {

                for (int n = 0; n < ncomp; n++)
                {
                    bool extdir_ilo = (d_bcrec_ptr[n].lo(0) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(0) == amrex::BCType::hoextrap);
                    bool extdir_ihi = (d_bcrec_ptr[n].hi(0) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(0) == amrex::BCType::hoextrap);
                    bool extdir_jlo = (d_bcrec_ptr[n].lo(1) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(1) == amrex::BCType::hoextrap);
                    bool extdir_jhi = (d_bcrec_ptr[n].hi(1) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(1) == amrex::BCType::hoextrap);
#if (AMREX_SPACEDIM == 3)
                    bool extdir_klo = (d_bcrec_ptr[n].lo(2) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].lo(2) == amrex::BCType::hoextrap);
                    bool extdir_khi = (d_bcrec_ptr[n].hi(2) == amrex::BCType::ext_dir ||
                                       d_bcrec_ptr[n].hi(2) == amrex::BCType::hoextrap);
#endif
                    // Initialize so that the slope stencil goes from -1:1 in each direction
                    int nx = 1; int ny = 1; int nz = 1;

                    // Do we have enough extent in each coordinate direction to use the 3x3x3 stencil
                    //    or do we need to enlarge it?
                    AMREX_D_TERM(Real x_max = -Real(1.e30); Real x_min = Real(1.e30);,
                                 Real y_max = -Real(1.e30); Real y_min = Real(1.e30);,
                                 Real z_max = -Real(1.e30); Real z_min = Real(1.e30););

                    Real slope_stencil_min_width = 0.5;
#if (AMREX_SPACEDIM == 2)
                    int kk = 0;
#elif (AMREX_SPACEDIM == 3)
                    for(int kk(-1); kk<=1; kk++) {
#endif
                    for(int jj(-1); jj<=1; jj++) {
                    for(int ii(-1); ii<=1; ii++) {
                        if (flag(i,j,k).isConnected(ii,jj,kk))
                        {
                            int r = i+ii; int s = j+jj; int t = k+kk;

                            x_max = amrex::max(x_max, cent_hat(r,s,t,0)+static_cast<Real>(ii));
                            x_min = amrex::min(x_min, cent_hat(r,s,t,0)+static_cast<Real>(ii));
                            y_max = amrex::max(y_max, cent_hat(r,s,t,1)+static_cast<Real>(jj));
                            y_min = amrex::min(y_min, cent_hat(r,s,t,1)+static_cast<Real>(jj));
#if (AMREX_SPACEDIM == 3)
                            z_max = amrex::max(z_max, cent_hat(r,s,t,2)+static_cast<Real>(kk));
                            z_min = amrex::min(z_min, cent_hat(r,s,t,2)+static_cast<Real>(kk));
#endif
                        } // isConnected
                    AMREX_D_TERM(},},})

                    // If we need to grow the stencil, we let it be -nx:nx in the x-direction,
                    //    for example.   Note that nx,ny,nz are either 1 or 2
                    if ( (x_max-x_min) < slope_stencil_min_width ) { nx = 2; }
                    if ( (y_max-y_min) < slope_stencil_min_width ) { ny = 2; }
#if (AMREX_SPACEDIM == 3)
                    if ( (z_max-z_min) < slope_stencil_min_width ) { nz = 2; }
#endif

                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes_eb;
                    if (nx*ny*nz == 1) {
                        // Compute slope using 3x3x3 stencil
                        slopes_eb = amrex_calc_slopes_extdir_eb(
                                                    i,j,k,n,Qhat,cent_hat,vfrac,
                                                    AMREX_D_DECL(fcx,fcy,fcz),flag,
                                                    AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                                                    AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                                                    AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                                    AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                                    max_order);

                    } else {
                        // Compute slope using grown stencil (no larger than 5x5x5)
                        slopes_eb = amrex_calc_slopes_extdir_eb_grown(
                                                    i,j,k,n,AMREX_D_DECL(nx,ny,nz),
                                                    Qhat,cent_hat,vfrac,
                                                    AMREX_D_DECL(fcx,fcy,fcz),flag,
                                                    AMREX_D_DECL(extdir_ilo, extdir_jlo, extdir_klo),
                                                    AMREX_D_DECL(extdir_ihi, extdir_jhi, extdir_khi),
                                                    AMREX_D_DECL(domain_ilo, domain_jlo, domain_klo),
                                                    AMREX_D_DECL(domain_ihi, domain_jhi, domain_khi),
                                                    max_order);
                    }

                    // We do the limiting separately because this limiter limits the slope based on the values
                    //    extrapolated to the cell centroid (cent_hat) locations - unlike the limiter in amrex
                    //    which bases the limiting on values extrapolated to the face centroids.
                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> lim_slope = slopes_eb;
                    amrex_apply_centroid_limiter(i,j,k,n,Qhat,flag,
                                                 AMREX_D_DECL(lim_slope[0],
                                                              lim_slope[1],
                                                              lim_slope[2]),
                                                 cent_hat);

                    // Add to the cell itself
                    if (bx.contains(IntVect(AMREX_D_DECL(i,j,k))))
                    {
                        Real update = Qhat(i,j,k,n);
                        AMREX_D_TERM(update += lim_slope[0] * (ccent(i,j,k,0)-cent_hat(i,j,k,0));,
                                     update += lim_slope[1] * (ccent(i,j,k,1)-cent_hat(i,j,k,1));,
                                     update += lim_slope[2] * (ccent(i,j,k,2)-cent_hat(i,j,k,2)););
                        amrex::Gpu::Atomic::Add(&U_out(i,j,k,n),alpha(i,j,k,0)*nrs(i,j,k)*update);
                    } // if bx contains

                    // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
                    for (int i_nbor = 1; i_nbor <= num_nbors; i_nbor++)
                    {
                        int r = i+imap[itracker(i,j,k,i_nbor)];
                        int s = j+jmap[itracker(i,j,k,i_nbor)];
                        int t = k+kmap[itracker(i,j,k,i_nbor)];

                        if (bx.contains(IntVect(AMREX_D_DECL(r,s,t))))
                        {
                            Real update = Qhat(i,j,k,n);
                            AMREX_D_TERM(update += lim_slope[0] * (ccent(r,s,t,0)-cent_hat(i,j,k,0) + static_cast<Real>(r-i));,
                                         update += lim_slope[1] * (ccent(r,s,t,1)-cent_hat(i,j,k,1) + static_cast<Real>(s-j));,
                                         update += lim_slope[2] * (ccent(r,s,t,2)-cent_hat(i,j,k,2) + static_cast<Real>(t-k)););
                            amrex::Gpu::Atomic::Add(&U_out(r,s,t,n),alpha(i,j,k,1)*update);
                        } // if bx contains
                    } // i_nbor
                } // n
            } // num_nbors
        } // vfrac
    });

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (!flag(i,j,k).isCovered())
        {
            // This seems to help with a compiler issue ...
            auto denom = Real(1.) / (nrs(i,j,k) + Real(1.e-30));
            U_out(i,j,k,n) *= denom;
        }
        else
        {
            U_out(i,j,k,n) = Real(1.e30);
        }
    });
}
}
