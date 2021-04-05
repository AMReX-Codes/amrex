#ifndef AMREX_EB_SLOPES_K_H_
#define AMREX_EB_SLOPES_K_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#else
#include <AMReX_FArrayBox.H>
#endif

#include <AMReX_Slopes_K.H>

namespace {

#if (AMREX_SPACEDIM > 1)
#ifdef AMREX_USE_EB

// amrex_calc_slopes_eb calculates the slope in each coordinate direction using a
// least squares linear fit to the 9 in 2D / 27 in 3D nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.  Note that this routine
// is called either by amrex_calc_slopes_eb or amrex_calc_slopes_extdir_eb; in the former
// A is defined with the cell centroids; in the latter, the A values corresponding to values
// defined on faces use the face centroid location.
//
// All the slopes are returned in one call.

#if (AMREX_SPACEDIM == 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb_given_A (int i, int j, int /*k*/, int n,
                              amrex::Real A[9][AMREX_SPACEDIM],
                              amrex::Array4<amrex::Real const> const& state,
                              amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    amrex::Real du[9];

    {
        int lc=0;
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){

            if( flag(i,j,0).isConnected(ii,jj,0) &&
                ! (ii==0 && jj==0)) {

              du[lc] = state(i+ii,j+jj,0,n) - state(i,j,0,n);

            } else {

              du[lc] = 0.0;
            }
            lc++;
          }
        }
    }

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj)
    {
      for(int ii(0); ii<AMREX_SPACEDIM; ++ii){
        AtA[ii][jj] = 0.0;
      }
      Atb[jj] = 0.0;
    }

    for(int lc(0); lc<9; ++lc)
    {
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
        AtA[1][1] += A[lc][1]* A[lc][1];

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];

    amrex::Real detAtA =
      (AtA[0][0]*AtA[1][1])-
      (AtA[0][1]*AtA[1][0]);

    amrex::Real detAtA_x =
        (Atb[0]   *AtA[1][1]) -
        (AtA[0][1]*Atb[1]);

    // Slope at centroid of (i,j,k)
    amrex::Real xs = detAtA_x / detAtA;

    amrex::Real detAtA_y =
        (AtA[0][0]*Atb[1]) -
        (Atb[0] * AtA[1][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real ys = detAtA_y / detAtA;

   return {xs,ys};
}
#elif (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb_given_A (int i, int j, int k, int n,
                              amrex::Real A[27][AMREX_SPACEDIM],
                              amrex::Array4<amrex::Real const> const& state,
                              amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    amrex::Real du[27];

    {
        int lc=0;
        for(int kk(-1); kk<=1; kk++)
        {
            for(int jj(-1); jj<=1; jj++){
              for(int ii(-1); ii<=1; ii++){

                if( flag(i,j,k).isConnected(ii,jj,kk) &&
                    ! (ii==0 && jj==0 && kk==0)) {

                  du[lc] = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);

                } else {

                  du[lc] = 0.0;
                }
                lc++;
              }
            }
        }
    }

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj){
      for(int ii(0); ii<AMREX_SPACEDIM; ++ii){
        AtA[ii][jj] = 0.0;
      }
      Atb[jj] = 0.0;
    }

    for(int lc(0); lc < 27; ++lc)
    {
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
        AtA[0][2] += A[lc][0]* A[lc][2];
        AtA[1][1] += A[lc][1]* A[lc][1];
        AtA[1][2] += A[lc][1]* A[lc][2];
        AtA[2][2] += A[lc][2]* A[lc][2];

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
        Atb[2] += A[lc][2]*du[lc];
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];
    AtA[2][0] = AtA[0][2];
    AtA[2][1] = AtA[1][2];

    amrex::Real detAtA =
      AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[2][1]) -
      AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
      AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    amrex::Real detAtA_x =
        Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
        AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );

    // Slope at centroid of (i,j,k)
    amrex::Real xs = detAtA_x / detAtA;

    amrex::Real detAtA_y =
        AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
        Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real ys = detAtA_y / detAtA;

    amrex::Real detAtA_z =
        AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
        Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

    // Slope at centroid of (i,j,k)
    amrex::Real zs = detAtA_z / detAtA;

    return {xs,ys,zs};
}
#endif

// amrex_calc_slopes_eb calculates the slope in each coordinate direction using a
// 1) standard 2nd order limited slope if all three cells in the stencil are regular cells
// OR
// 2) least squares linear fit to the at-most 9 in 2D / 27 in 3D nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.  Note that calc_slopes_eb
// defines the matrix A and passes this A to amrex_calc_slopes_eb_given_A.
//
// This routine assumes that there are no relevant hoextrap/extdir domain boundary conditions for this cell --
//     it does not test for them so this should not be called if such boundaries might be present
//
// All the slopes are returned in one call.

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb (int i, int j, int k, int n,
                      amrex::Array4<amrex::Real const> const& state,
                      amrex::Array4<amrex::Real const> const& ccent,
                      amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
#if (AMREX_SPACEDIM == 2)
    constexpr int dim_a = 9;
#elif (AMREX_SPACEDIM == 3)
    constexpr int dim_a = 27;
#endif
    amrex::Real A[dim_a][AMREX_SPACEDIM];

    int lc=0;
#if (AMREX_SPACEDIM == 3)
    for(int kk(-1); kk<=1; kk++)
#else
    int kk = 0;
#endif
    {
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){

            if( flag(i,j,k).isConnected(ii,jj,kk) &&
                ! (ii==0 && jj==0 && kk==0)) {

            // Not multplying by dx to be consistent with how the
            // slope is stored. Also not including the global shift
            // wrt plo or i,j,k. We only need relative distance.

              A[lc][0] = ii + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
              A[lc][1] = jj + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
              A[lc][2] = kk + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2);
#endif

            } else {

              AMREX_D_TERM(A[lc][0] = 0.0;,
                           A[lc][1] = 0.0;,
                           A[lc][2] = 0.0;);
            }
            lc++;
          }
        }
    }

    //
    // These slopes use the EB stencil without testing whether it is actually needed
    //
    const auto& slopes = amrex_calc_slopes_eb_given_A (i,j,k,n,A,state,flag);
    AMREX_D_TERM(amrex::Real xslope = slopes[0];,
                 amrex::Real yslope = slopes[1];,
                 amrex::Real zslope = slopes[2];);

    //
    // Here we over-write -- if possible -- with a stencil not using the EB stencil
    //

    // X direction
    if (flag(i,j,k).isRegular() && flag(i-1,j,k).isRegular() && flag(i+1,j,k).isRegular())
    {
        int order = 2;
        xslope = amrex_calc_xslope(i,j,k,n,order,state);
    }

    // Y direction
    if (flag(i,j,k).isRegular() && flag(i,j-1,k).isRegular() && flag(i,j+1,k).isRegular())
    {
        int order = 2;
        yslope = amrex_calc_yslope(i,j,k,n,order,state);
    }

#if (AMREX_SPACEDIM == 3)
    // Z direction
    if (flag(i,j,k).isRegular() && flag(i,j,k-1).isRegular() && flag(i,j,k+1).isRegular())
    {
        int order = 2;
        zslope = amrex_calc_zslope(i,j,k,n,order,state);
    }
#endif

    return {AMREX_D_DECL(xslope,yslope,zslope)};
}

// amrex_calc_slopes_extdir_eb calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
//    (this stencil sees the extdir/hoextrap boundary condition if there is one)
// OR
// 2) least squares linear fit to the at-most 9 in 2D / 27 in 3D nearest neighbors, with the function
// least squares linear fit to the 9 in 2D / 27 in 3D nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.
//
// All the slopes are returned in one call.
//

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_extdir_eb (int i, int j, int k, int n,
                             amrex::Array4<amrex::Real const> const& state,
                             amrex::Array4<amrex::Real const> const& ccent,
                             AMREX_D_DECL(
                             amrex::Array4<amrex::Real const> const& fcx,
                             amrex::Array4<amrex::Real const> const& fcy,
                             amrex::Array4<amrex::Real const> const& fcz),
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                             AMREX_D_DECL(bool edlo_x, bool edlo_y, bool edlo_z),
                             AMREX_D_DECL(bool edhi_x, bool edhi_y, bool edhi_z),
                             AMREX_D_DECL(int domlo_x, int domlo_y, int domlo_z),
                             AMREX_D_DECL(int domhi_x, int domhi_y, int domhi_z)) noexcept
{
#if (AMREX_SPACEDIM == 2)
    constexpr int dim_a = 9;
#elif (AMREX_SPACEDIM == 3)
    constexpr int dim_a = 27;
#endif

    AMREX_D_TERM(amrex::Real xslope = 0.0;,
                 amrex::Real yslope = 0.0;,
                 amrex::Real zslope = 0.0;);

    // First get EB-aware slope that doesn't know about extdir
    bool needs_bdry_stencil = (edlo_x && i <= domlo_x) || (edhi_x && i >= domhi_x) ||
                              (edlo_y && j <= domlo_y) || (edhi_y && j >= domhi_y);
#if (AMREX_SPACEDIM == 3)
         needs_bdry_stencil = needs_bdry_stencil ||
                              (edlo_z && k <= domlo_z) || (edhi_z && k >= domhi_z);
#endif

    //
    // This call does not have any knowledge of extdir / hoextrap boundary conditions
    //
    if (!needs_bdry_stencil)
    {
        // This returns slopes calculated with the regular 1-d approach if all cells in the stencil
        //      are regular.  If not, it uses the EB-aware least squares approach to fit a linear profile
        //      using the neighboring un-covered cells.
        const auto& slopes = amrex_calc_slopes_eb (i,j,k,n,state,ccent,flag);
        return slopes;

    } else {

    amrex::Real A[dim_a][AMREX_SPACEDIM];

    int lc=0;
#if (AMREX_SPACEDIM == 3)
    for(int kk(-1); kk<=1; kk++)
#else
    int kk = 0;
#endif
    {
        for(int jj(-1); jj<=1; jj++)
          for(int ii(-1); ii<=1; ii++)
          {
            if( flag(i,j,k).isConnected(ii,jj,kk) &&
                ! (ii==0 && jj==0 && kk==0))
                {
                    bool ilo_test = ( edlo_x && (i == domlo_x) && ii == -1);
                    bool ihi_test = ( edhi_x && (i == domhi_x) && ii ==  1);

                    bool jlo_test = ( edlo_y && (j == domlo_y) && jj == -1);
                    bool jhi_test = ( edhi_y && (j == domhi_y) && jj ==  1);

#if (AMREX_SPACEDIM == 3)
                    bool klo_test = ( edlo_z && (k == domlo_z) && kk == -1);
                    bool khi_test = ( edhi_z && (k == domhi_z) && kk ==  1);
#else
                    bool klo_test = false;
                    bool khi_test = false;
#endif

                    // These are the default values if no physical boundary
                    A[lc][0] = ii + ccent(i+ii,j+jj,k+kk,0);
                    A[lc][1] = jj + ccent(i+ii,j+jj,k+kk,1);
#if (AMREX_SPACEDIM == 3)
                    A[lc][2] = kk + ccent(i+ii,j+jj,k+kk,2);
#endif
                   // Do corrections for entire x-face
                   if (ilo_test)
                   {
                        if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                        {
                            A[lc][1] = jj + fcx(i   ,j+jj,k+kk,0);
#if (AMREX_SPACEDIM == 3)
                            A[lc][2] = kk + fcx(i   ,j+jj,k+kk,1);
#endif
                        }
                        A[lc][0] = -0.5                      ;
                    } else if (ihi_test) {

                        if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                        {
                            A[lc][1] = jj + fcx(i+ii,j+jj,k+kk,0);
#if (AMREX_SPACEDIM == 3)
                            A[lc][2] = kk + fcx(i+ii,j+jj,k+kk,1);
#endif
                        }
                        A[lc][0] = 0.5                       ;
                   }

                   // Do corrections for entire y-face
                   if (jlo_test)
                   {
                        if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                        {
                            A[lc][0] = ii + fcy(i+ii,j   ,k+kk,0);
#if (AMREX_SPACEDIM == 3)
                            A[lc][2] = kk + fcy(i+ii,j   ,k+kk,1);
#endif
                        }
                        A[lc][1] = -0.5                      ;

                   } else if (jhi_test) {

                        if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                        {
                            A[lc][0] = ii + fcy(i+ii,j+jj,k+kk,0);
#if (AMREX_SPACEDIM == 3)
                            A[lc][2] = kk + fcy(i+ii,j+jj,k+kk,1);
#endif
                        }
                        A[lc][1] = 0.5                       ;
                    }

#if (AMREX_SPACEDIM == 3)
                   // Do corrections for entire z-face
                   if (klo_test)
                   {
                        if (!ilo_test && !ihi_test && !jlo_test && !jhi_test)
                        {
                            A[lc][0] = ii + fcz(i+ii,j+jj,k   ,0);
                            A[lc][1] = jj + fcz(i+ii,j+jj,k   ,1);
                        }
                        A[lc][2] = -0.5                      ;

                   } else if (khi_test) {
                        if (!ilo_test && !ihi_test && !jlo_test && !jhi_test)
                        {
                            A[lc][0] = ii + fcz(i+ii,j+jj,k+kk,0);
                            A[lc][1] = jj + fcz(i+ii,j+jj,k+kk,1);
                        }
                        A[lc][2] = 0.5                       ;
                   }
#endif
                   AMREX_D_TERM(A[lc][0] -= ccent(i,j,k,0);,
                                A[lc][1] -= ccent(i,j,k,1);,
                                A[lc][2] -= ccent(i,j,k,2););

                } else {
                    AMREX_D_TERM(A[lc][0] = 0.0;,
                                 A[lc][1] = 0.0;,
                                 A[lc][2] = 0.0;);
                }
                lc++;
            } // i,j
    } // k

    const auto& slopes = amrex_calc_slopes_eb_given_A (i,j,k,n,A,state,flag);
    AMREX_D_TERM(xslope = slopes[0];,
                 yslope = slopes[1];,
                 zslope = slopes[2];);

    //
    // Correct only those cells which are affected by extdir but not by EB:
    //    2) If all the cells are regular we use the "regular slope" in the extdir direction
    //

    int order = 2;

    // Overwrite the tangential slope with the regular stencils if we can compute from non-EB cells
    if ( flag(i,j,k).isRegular() && flag(i-1,j,k).isRegular() && flag(i+1,j,k).isRegular() )
      xslope = amrex_calc_xslope_extdir(i,j,k,n,order,state,edlo_x,edhi_x,domlo_x,domhi_x);

    if ( flag(i,j,k).isRegular() && flag(i,j-1,k).isRegular() && flag(i,j+1,k).isRegular() )
      yslope = amrex_calc_yslope_extdir(i,j,k,n,order,state,edlo_y,edhi_y,domlo_y,domhi_y);

#if (AMREX_SPACEDIM == 3)
    if ( flag(i,j,k).isRegular() && flag(i,j,k-1).isRegular() && flag(i,j,k+1).isRegular() )
      zslope = amrex_calc_zslope_extdir(i,j,k,n,order,state,edlo_z,edhi_z,domlo_z,domhi_z);
#endif
    } // end of needs_bndry_stencil

    // Zero out slopes outside of an extdir (or hoextrap) boundary
    // TODO:  is this the right thing to do at a HOEXTRAP boundary??
#if (AMREX_SPACEDIM == 2)
    if ( (edlo_x && i < domlo_x) || (edhi_x && i > domhi_x)   ||
         (edlo_y && j < domlo_y) || (edhi_y && j > domhi_y) )
    {
          AMREX_D_TERM(xslope = 0.;,
                       xslope = 0.; yslope = 0.;,
                       xslope = 0.; yslope = 0.; zslope = 0.;);
    }
#elif (AMREX_SPACEDIM == 3)
    if ( (edlo_x && i < domlo_x) || (edhi_x && i > domhi_x) ||
         (edlo_y && j < domlo_y) || (edhi_y && j > domhi_y) ||
         (edlo_z && k < domlo_z) || (edhi_z && k > domhi_z) )
    {
          AMREX_D_TERM(xslope = 0.;,
                       xslope = 0.; yslope = 0.;,
                       xslope = 0.; yslope = 0.; zslope = 0.;);
    }
#endif

    return {AMREX_D_DECL(xslope,yslope,zslope)};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
amrex_calc_alpha_stencil(amrex::Real q_hat, amrex::Real q_max,
                         amrex::Real q_min, amrex::Real state, amrex::Real alpha) noexcept
{
    amrex::Real alpha_temp = 0.0;
    amrex::Real small  = 1.0e-13;

    if ((q_hat-state) > small) {
        alpha_temp = amrex::min(1.0,(q_max-state)/(q_hat-state));
    } else if ((q_hat-state) < -small) {
        alpha_temp = amrex::min(1.0,(q_min-state)/(q_hat-state));
    } else {
        alpha_temp = 1.0;
    }
    return amrex::min(alpha, alpha_temp);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_alpha_limiter(int i, int j, int k, int n,
                         amrex::Array4<amrex::Real const> const& state,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& slopes,
                         AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                      amrex::Array4<amrex::Real const> const& fcy,
                                      amrex::Array4<amrex::Real const> const& fcz),
                         amrex::Array4<amrex::Real const> const& ccent) noexcept
{
    amrex::Real alpha = 2.0;
    amrex::Real q_min = state(i,j,k,n);
    amrex::Real q_max = state(i,j,k,n);

    AMREX_D_TERM(int cuts_x = 0;,
                 int cuts_y = 0;,
                 int cuts_z = 0;);

    //Compute max and min values in a 3x3 stencil and how many cut or regular faces do we have
#if (AMREX_SPACEDIM == 3)
    for(int kk(-1); kk<=1; kk++)
#else
    int kk = 0;
#endif
    {
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){
            if( flag(i,j,k).isConnected(ii,jj,kk) and
                not (ii==0 and jj==0 and kk==0)) {
                if (state(i+ii,j+jj,k+kk,n) > q_max) q_max = state(i+ii,j+jj,k+kk,n);
                if (state(i+ii,j+jj,k+kk,n) < q_min) q_min = state(i+ii,j+jj,k+kk,n);

                if ((ii==-1 or ii==1) and jj==0 and kk==0) cuts_x++;
                if ((jj==-1 or jj==1) and ii==0 and kk==0) cuts_y++;
#if (AMREX_SPACEDIM == 3)
                if ((kk==-1 or kk==1) and ii==0 and jj==0) cuts_z++;
#endif
            }
          }
        }
    }

    AMREX_D_TERM(amrex::Real xc = ccent(i,j,k,0);, // centroid of cell (i,j,k)
                 amrex::Real yc = ccent(i,j,k,1);,
                 amrex::Real zc = ccent(i,j,k,2););

    //Reconstruct values at the face centroids and compute the limiter
    if(flag(i,j,k).isConnected(0,1,0)) {
           amrex::Real xf = fcy(i,j+1,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
           amrex::Real zf = fcy(i,j+1,k,1); // local (x,z) of centroid of y-face we are extrapolating to
#endif

           AMREX_D_TERM(amrex::Real delta_x = xf  - xc;,
                        amrex::Real delta_y = 0.5 - yc;,
                        amrex::Real delta_z = zf  - zc;);

#if (AMREX_SPACEDIM == 3)
           amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                              + delta_y * slopes[1]
                                              + delta_z * slopes[2];
#else
           amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                              + delta_y * slopes[1];
#endif

           alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if (flag(i,j,k).isConnected(0,-1,0)){
           amrex::Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
           amrex::Real zf = fcy(i,j,k,1); // local (x,z) of centroid of y-face we are extrapolating to
#endif
           AMREX_D_TERM(amrex::Real delta_x = xf  - xc;,
                        amrex::Real delta_y = 0.5 + yc;,
                        amrex::Real delta_z = zf  - zc;);

#if (AMREX_SPACEDIM == 3)
           amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                              - delta_y * slopes[1]
                                              + delta_z * slopes[2];
#else
           amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                              - delta_y * slopes[1];
#endif

           alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if(flag(i,j,k).isConnected(1,0,0)) {
           amrex::Real yf = fcx(i+1,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
           amrex::Real zf = fcx(i+1,j,k,1); // local (y,z) of centroid of x-face we are extrapolating to
#endif
           AMREX_D_TERM(amrex::Real delta_x = 0.5 - xc;,
                        amrex::Real delta_y = yf  - yc;,
                        amrex::Real delta_z = zf  - zc;);

#if (AMREX_SPACEDIM == 3)
           amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                              + delta_y * slopes[1]
                                              + delta_z * slopes[2];
#else
           amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                              + delta_y * slopes[1];
#endif

           alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if(flag(i,j,k).isConnected(-1,0,0)) {
           amrex::Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
#if (AMREX_SPACEDIM == 3)
           amrex::Real zf = fcx(i,j,k,1); // local (y,z) of centroid of x-face we are extrapolating to
#endif
           AMREX_D_TERM(amrex::Real delta_x = 0.5 + xc;,
                        amrex::Real delta_y = yf  - yc;,
                        amrex::Real delta_z = zf  - zc;);
#if (AMREX_SPACEDIM == 3)
           amrex::Real q_hat = state(i,j,k,n) - delta_x * slopes[0]
                                              + delta_y * slopes[1]
                                              + delta_z * slopes[2];
#else
           amrex::Real q_hat = state(i,j,k,n) - delta_x * slopes[0]
                                              + delta_y * slopes[1];
#endif
           alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
#if (AMREX_SPACEDIM == 3)
    if(flag(i,j,k).isConnected(0,0,1)) {
           amrex::Real xf = fcz(i,j,k+1,0); // local (x,y) of centroid of z-face we are extrapolating to
           amrex::Real yf = fcz(i,j,k+1,1); // local (x,y) of centroid of z-face we are extrapolating to

           AMREX_D_TERM(amrex::Real delta_x = xf  - xc;,
                        amrex::Real delta_y = yf  - yc;,
                        amrex::Real delta_z = 0.5 - zc;);

           amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                              + delta_y * slopes[1]
                                              + delta_z * slopes[2];

           alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if(flag(i,j,k).isConnected(0,0,-1)) {
           amrex::Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
           amrex::Real yf = fcz(i,j,k,1); // local (x,y) of centroid of z-face we are extrapolating to

           AMREX_D_TERM(amrex::Real delta_x = xf  - xc;,
                        amrex::Real delta_y = yf  - yc;,
                        amrex::Real delta_z = 0.5 + zc;);

           amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                              + delta_y * slopes[1]
                                              - delta_z * slopes[2];

           alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
#endif

    AMREX_D_TERM(amrex::Real xalpha = alpha;,
                 amrex::Real yalpha = alpha;,
                 amrex::Real zalpha = alpha;);

    //Zeroing out the slopes in the direction where a covered face exists.
    if (cuts_x<2) xalpha = 0;
    if (cuts_y<2) yalpha = 0;
#if (AMREX_SPACEDIM == 3)
    if (cuts_z<2) zalpha = 0;
#endif

    return {AMREX_D_DECL(xalpha,yalpha,zalpha)};
}

//amrex_lim_slopes_eb computes the slopes calling amrex_calc_slopes_eb, and then each slope component
//is multiplied by a limiter based on the work of Barth-Jespersen.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_lim_slopes_eb (int i, int j, int k, int n,
                     amrex::Array4<amrex::Real const> const& state,
                     amrex::Array4<amrex::Real const> const& ccent,
                     AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                  amrex::Array4<amrex::Real const> const& fcy,
                                  amrex::Array4<amrex::Real const> const& fcz),
                     amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> alpha_lim;

    slopes = amrex_calc_slopes_eb(i,j,k,n,state,ccent,flag);
    alpha_lim = amrex_calc_alpha_limiter(i,j,k,n,state,flag,slopes,AMREX_D_DECL(fcx,fcy,fcz),ccent);

    // Setting limiter to 1 for stencils that just consists of non-EB cells because
    // amrex_calc_slopes_eb routine will call the slope routine for non-EB stencils that has already a limiter
    if ( flag(i,j,k).isRegular() and flag(i-1,j,k).isRegular() and flag(i+1,j,k).isRegular() )
      alpha_lim[0] = 1.0;

    if ( flag(i,j,k).isRegular() and flag(i,j-1,k).isRegular() and flag(i,j+1,k).isRegular() )
      alpha_lim[1] = 1.0;

#if (AMREX_SPACEDIM == 3)
    if ( flag(i,j,k).isRegular() and flag(i,j,k-1).isRegular() and flag(i,j,k+1).isRegular() )
      alpha_lim[2] = 1.0;
#endif

    return {AMREX_D_DECL(alpha_lim[0]*slopes[0],alpha_lim[1]*slopes[1],alpha_lim[2]*slopes[2])};
}

//amrex_lim_slopes_extdir_eb computes the slopes calling amrex_calc_slopes_extdir_eb, and then each slope component
//is multiplied by a limiter based on the work of Barth-Jespersen.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_lim_slopes_extdir_eb (int i, int j, int k, int n,
                            amrex::Array4<amrex::Real const> const& state,
                            amrex::Array4<amrex::Real const> const& ccent,
                            AMREX_D_DECL(
                            amrex::Array4<amrex::Real const> const& fcx,
                            amrex::Array4<amrex::Real const> const& fcy,
                            amrex::Array4<amrex::Real const> const& fcz),
                            amrex::Array4<amrex::EBCellFlag const> const& flag,
                            AMREX_D_DECL(bool edlo_x, bool edlo_y, bool edlo_z),
                            AMREX_D_DECL(bool edhi_x, bool edhi_y, bool edhi_z),
                            AMREX_D_DECL(int domlo_x, int domlo_y, int domlo_z),
                            AMREX_D_DECL(int domhi_x, int domhi_y, int domhi_z)) noexcept
{

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> alpha_lim;

    slopes = amrex_calc_slopes_extdir_eb(i,j,k,n,state,ccent,AMREX_D_DECL(fcx,fcy,fcz),flag,
                                         AMREX_D_DECL(edlo_x,edlo_y,edlo_z),AMREX_D_DECL(edhi_x,edhi_y,edhi_z),
                                         AMREX_D_DECL(domlo_x,domlo_y,domlo_z),AMREX_D_DECL(domhi_x,domhi_y,domhi_z));
    alpha_lim = amrex_calc_alpha_limiter(i,j,k,n,state,flag,slopes,AMREX_D_DECL(fcx,fcy,fcz),ccent);

    // Setting limiter to 1 for stencils that just consists of non-EB cells because
    // amrex_calc_slopes_extdir_eb routine will call the slope routine for non-EB stencils that has already a limiter
    if ( flag(i,j,k).isRegular() and flag(i-1,j,k).isRegular() and flag(i+1,j,k).isRegular() )
      alpha_lim[0] = 1.0;

    if ( flag(i,j,k).isRegular() and flag(i,j-1,k).isRegular() and flag(i,j+1,k).isRegular() )
      alpha_lim[1] = 1.0;

#if (AMREX_SPACEDIM == 3)
    if ( flag(i,j,k).isRegular() and flag(i,j,k-1).isRegular() and flag(i,j,k+1).isRegular() )
      alpha_lim[2] = 1.0;
#endif

    return {AMREX_D_DECL(alpha_lim[0]*slopes[0],alpha_lim[1]*slopes[1],alpha_lim[2]*slopes[2])};
}

#endif
#endif

}
#endif
