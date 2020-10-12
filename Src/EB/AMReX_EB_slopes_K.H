#ifndef AMREX_EB_SLOPES_K_H_
#define AMREX_EB_SLOPES_K_H_

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
amrex_calc_slopes_eb_given_A (int i, int j, int k, int n,
                              amrex::Real A[9][AMREX_SPACEDIM],
                              amrex::Array4<amrex::Real const> const& state,
                              amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    amrex::Real du[9];

    int lc=0;
    {
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){

            if( flag(i,j,0).isConnected(ii,jj,0) and
                not (ii==0 and jj==0)) {

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

    int lc=0;
    for(int kk(-1); kk<=1; kk++)
    {
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){

            if( flag(i,j,k).isConnected(ii,jj,kk) and
                not (ii==0 and jj==0 and kk==0)) {

              du[lc] = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);

            } else {

              du[lc] = 0.0;
            }
            lc++;
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

            if( flag(i,j,k).isConnected(ii,jj,kk) and
                not (ii==0 and jj==0 and kk==0)) {

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
    if (flag(i,j,k).isRegular() and flag(i-1,j,k).isRegular() and flag(i+1,j,k).isRegular()) 
    {
        int order = 2;
        xslope = amrex_calc_xslope(i,j,k,n,order,state);
    }

    // Y direction
    if (flag(i,j,k).isRegular() and flag(i,j-1,k).isRegular() and flag(i,j+1,k).isRegular()) 
    {
        int order = 2;
        yslope = amrex_calc_yslope(i,j,k,n,order,state);
    }

#if (AMREX_SPACEDIM == 3)
    // Z direction
    if (flag(i,j,k).isRegular() and flag(i,j,k-1).isRegular() and flag(i,j,k+1).isRegular()) 
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
    bool needs_bdry_stencil = (edlo_x and i <= domlo_x) or (edhi_x and i >= domhi_x) or
                              (edlo_y and j <= domlo_y) or (edhi_y and j >= domhi_y);
#if (AMREX_SPACEDIM == 3)
         needs_bdry_stencil = needs_bdry_stencil or 
                              (edlo_z and k <= domlo_z) or (edhi_z and k >= domhi_z);
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

    //
    // Correct only those cells which are affected by extdir:
    //    1) Here if any of the cells are not regular we use the face and cell centroids
    //

    if ( (edlo_x and i == domlo_x) or (edhi_x and i == domhi_x) )
    {
        int lc=0;
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
    
                    // Not multplying by dx to be consistent with how the
                    // slope is stored. Also not including the global shift
                    // wrt plo or i,j,k. We only need relative distance.

                    if ( (edlo_x and i == domlo_x) and ii == -1) {
                        A[lc][0] = -0.5                       - ccent(i,j,k,0);
                        A[lc][1] = jj + fcx(i   ,j+jj,k+kk,0) - ccent(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
                        A[lc][2] = kk + fcx(i   ,j+jj,k+kk,1) - ccent(i,j,k,2);
#endif
                    } 
                    else if ( (edhi_x and i == domhi_x) and ii == 1) 
                    {
                        A[lc][0] = 0.5                        - ccent(i,j,k,0);
                        A[lc][1] = jj + fcx(i+ii,j+jj,k+kk,0) - ccent(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
                        A[lc][2] = kk + fcx(i+ii,j+jj,k+kk,1) - ccent(i,j,k,2);
#endif
                    } else {
                        A[lc][0] = ii + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
                        A[lc][1] = jj + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
                        A[lc][2] = kk + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2);
#endif
                    }

                 } else {
                    AMREX_D_TERM(A[lc][0] = 0.0;,
                                 A[lc][1] = 0.0;,
                                 A[lc][2] = 0.0;);
                }
                lc++;
              }
            }
        }
        const auto& slopes = amrex_calc_slopes_eb_given_A (i,j,k,n,A,state,flag);
        AMREX_D_TERM(xslope = slopes[0];,
                     yslope = slopes[1];,
                     zslope = slopes[2];);
    }

    if ( (edlo_y and j == domlo_y) or (edhi_y and j == domhi_y) )
    {
        int lc=0;
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
    
                    // Not multplying by dx to be consistent with how the
                    // slope is stored. Also not including the global shift
                    // wrt plo or i,j,k. We only need relative distance.

                    if (edlo_y and j == domlo_y and jj == -1) {
                        A[lc][0] = ii + fcy(i+ii,j   ,k+kk,0) - ccent(i,j,k,0);
                        A[lc][1] = -0.5                       - ccent(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
                        A[lc][2] = kk + fcy(i+ii,j   ,k+kk,1) - ccent(i,j,k,2);
#endif
                    } 
                    else if (edhi_y and j == domhi_y and jj == 1) 
                    {
                        A[lc][0] = ii + fcy(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
                        A[lc][1] = 0.5                        - ccent(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
                        A[lc][2] = kk + fcy(i+ii,j+jj,k+kk,1) - ccent(i,j,k,2);
#endif
                    } else {
                        A[lc][0] = ii + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
                        A[lc][1] = jj + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
                        A[lc][2] = kk + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2);
#endif
                    }

                 } else {
                    AMREX_D_TERM(A[lc][0] = 0.0;,
                                 A[lc][1] = 0.0;,
                                 A[lc][2] = 0.0;);
                }

                lc++;
              }
            }
        }

        const auto& slopes = amrex_calc_slopes_eb_given_A (i,j,k,n,A,state,flag);
        AMREX_D_TERM(xslope = slopes[0];,
                     yslope = slopes[1];,
                     zslope = slopes[2];);
    }

#if (AMREX_SPACEDIM == 3)
    if ( (edlo_z and k == domlo_z) or (edhi_z and k == domhi_z) )
    {
        int lc=0;
        for(int kk(-1); kk<=1; kk++) {
            for(int jj(-1); jj<=1; jj++){
              for(int ii(-1); ii<=1; ii++){

                if( flag(i,j,k).isConnected(ii,jj,kk) and
                    not (ii==0 and jj==0 and kk==0)) 
                {
                    // Not multplying by dx to be consistent with how the
                    // slope is stored. Also not including the global shift
                    // wrt plo or i,j,k. We only need relative distance.

                    if (edlo_z and k == domlo_z and kk == -1) {
                        A[lc][0] = ii + fcz(i+ii,j+jj,k   ,0) - ccent(i,j,k,0);
                        A[lc][1] = jj + fcz(i+ii,j+jj,k   ,1) - ccent(i,j,k,1);
                        A[lc][2] = -0.5                       - ccent(i,j,k,2);
                    } 
                    else if (edhi_z and k == domhi_z and kk == 1) 
                    {
                        A[lc][0] = ii + fcz(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
                        A[lc][1] = jj + fcz(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
                        A[lc][2] = 0.5                        - ccent(i,j,k,2);
                    } else {
                        A[lc][0] = ii + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
                        A[lc][1] = jj + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
                        A[lc][2] = kk + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2);
                    }

                 } else {
                    AMREX_D_TERM(A[lc][0] = 0.0;,
                                 A[lc][1] = 0.0;,
                                 A[lc][2] = 0.0;);
                }

                lc++;
              }
            }
        }
        const auto& slopes = amrex_calc_slopes_eb_given_A (i,j,k,n,A,state,flag);

        AMREX_D_TERM(xslope = slopes[0];,
                     yslope = slopes[1];,
                     zslope = slopes[2];);
    }
#endif

    //
    // Correct only those cells which are affected by extdir but not by EB:
    //    2) If all the cells are regular we use the "regular slope" in the extdir direction
    //

    int order = 2;

    // Overwrite the tangential slope with the regular stencils if we can compute from non-EB cells
    if ( flag(i,j,k).isRegular() and flag(i-1,j,k).isRegular() and flag(i+1,j,k).isRegular() )
      xslope = amrex_calc_xslope_extdir(i,j,k,n,order,state,edlo_x,edhi_x,domlo_x,domhi_x);

    if ( flag(i,j,k).isRegular() and flag(i,j-1,k).isRegular() and flag(i,j+1,k).isRegular() ) 
      yslope = amrex_calc_yslope_extdir(i,j,k,n,order,state,edlo_y,edhi_y,domlo_y,domhi_y);

#if (AMREX_SPACEDIM == 3)
    if ( flag(i,j,k).isRegular() and flag(i,j,k-1).isRegular() and flag(i,j,k+1).isRegular() ) 
      zslope = amrex_calc_zslope_extdir(i,j,k,n,order,state,edlo_z,edhi_z,domlo_z,domhi_z);
#endif
    }

    // Zero out slopes outside of an extdir (or hoextrap) boundary
    // TODO:  is this the right thing to do at a HOEXTRAP boundary??
#if (AMREX_SPACEDIM == 2)
    if ( (edlo_x and i < domlo_x) or (edhi_x and i > domhi_x)   or
         (edlo_y and j < domlo_y) or (edhi_y and j > domhi_y) )  
    { 
          AMREX_D_TERM(xslope = 0.;,  
                       xslope = 0.; yslope = 0.;, 
                       xslope = 0.; yslope = 0.; zslope = 0.;);
    } 
#elif (AMREX_SPACEDIM == 3)
    if ( (edlo_x and i < domlo_x) or (edhi_x and i > domhi_x) or
         (edlo_y and j < domlo_y) or (edhi_y and j > domhi_y) or
         (edlo_z and k < domlo_z) or (edhi_z and k > domhi_z) )
    {
          AMREX_D_TERM(xslope = 0.;,  
                       xslope = 0.; yslope = 0.;, 
                       xslope = 0.; yslope = 0.; zslope = 0.;);
    }
#endif

    return {AMREX_D_DECL(xslope,yslope,zslope)};
}

#endif
#endif

}
#endif
