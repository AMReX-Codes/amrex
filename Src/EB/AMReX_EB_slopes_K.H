#ifndef AMREX_EB_SLOPES_K_H_
#define AMREX_EB_SLOPES_K_H_

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#else
#include <AMReX_FArrayBox.H>
#endif

namespace {

#if (AMREX_SPACEDIM > 1)
#ifdef AMREX_USE_EB

// amrex_calc_slopes_eb calculates the slope in each coordinate direction using a 
// least squares linear fit to the 9 in 2D / 27 in 3D nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.
//
// No boundary conditions are used in this call.
//
// All the slopes are returned in one call.

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb (int i, int j, int k, int n,
                      amrex::Array4<amrex::Real const> const& state,
                      amrex::Array4<amrex::Real const> const& ccent,
                      amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{

#if (AMREX_SPACEDIM == 3)
    constexpr int m_size = 27;
#else
    constexpr int m_size = 9;
#endif

    amrex::Real A[m_size][AMREX_SPACEDIM];
    amrex::Real du[m_size];

    {
      int lc=0;
#if (AMREX_SPACEDIM == 3)
      for(int kk(-1); kk<=1; kk++){
#else
      int kk = 0;
#endif
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
              du[lc] = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);

            } else {

              AMREX_D_TERM(A[lc][0] = 0.0;,
                           A[lc][1] = 0.0;,
                           A[lc][2] = 0.0;);

              du[lc] = 0.0;
            }

            lc++;
          }
        }
#if (AMREX_SPACEDIM == 3)
      }
#endif
    }

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj){
      for(int ii(0); ii<AMREX_SPACEDIM; ++ii){
        AtA[ii][jj] = 0.0;
      }
      Atb[jj] = 0.0;
    }

    {

      for(int lc(0); lc<m_size; ++lc){
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
#if (AMREX_SPACEDIM == 3)
        AtA[0][2] += A[lc][0]* A[lc][2];
#endif
        AtA[1][1] += A[lc][1]* A[lc][1];
#if (AMREX_SPACEDIM == 3)
        AtA[1][2] += A[lc][1]* A[lc][2];
        AtA[2][2] += A[lc][2]* A[lc][2];
#endif

        Atb[0] += A[lc][0]*du[lc];
        Atb[1] += A[lc][1]*du[lc];
#if (AMREX_SPACEDIM == 3)
        Atb[2] += A[lc][2]*du[lc];
#endif
      }
    }

    // Fill in symmetric
    AtA[1][0] = AtA[0][1];
#if (AMREX_SPACEDIM == 3)
    AtA[2][0] = AtA[0][2];
    AtA[2][1] = AtA[1][2];
#endif

#if (AMREX_SPACEDIM == 3)
    amrex::Real detAtA =
      AtA[0][0]*(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[2][1]) -
      AtA[0][1]*(AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
      AtA[0][2]*(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);
#else
    amrex::Real detAtA =
      (AtA[0][0]*AtA[1][1])-
      (AtA[0][1]*AtA[1][0]);
#endif

    AMREX_D_TERM(amrex::Real xs = 0.0;,
                 amrex::Real ys = 0.0;,
                 amrex::Real zs = 0.0;);

    // X direction
    if(flag(i  ,j,k).isSingleValued() or
      (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
      (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) 
    {

#if (AMREX_SPACEDIM == 3)
      amrex::Real detAtA_x =
        Atb[0]   *(AtA[1][1]*AtA[2][2] - AtA[1][2]*AtA[1][2]) -
        AtA[0][1]*(Atb[1] *  AtA[2][2] - AtA[1][2]*Atb[2]   ) +
        AtA[0][2]*(Atb[1] *  AtA[2][1] - AtA[1][1]*Atb[2]   );
#else
      amrex::Real detAtA_x =
        (Atb[0]   *AtA[1][1]) -
        (AtA[0][1]*Atb[1]); 
#endif

      // Slope at centroid of (i,j,k)
      xs = detAtA_x / detAtA;

    } else {

      amrex::Real du_xl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
      amrex::Real du_xr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
      amrex::Real du_xc = 0.5*(state(i+1,j,k,n) - state(i-1,j,k,n));

      amrex::Real xslope = amrex::min(amrex::Math::abs(du_xl), amrex::Math::abs(du_xc), amrex::Math::abs(du_xr));
      xslope = (du_xr*du_xl > 0.0) ? xslope : 0.0;
      xs = (du_xc > 0.0) ? xslope : -xslope;

    }

    // Y direction
    if(flag(i,j  ,k).isSingleValued() or
      (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
      (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) 
    {

#if (AMREX_SPACEDIM == 3)
      amrex::Real detAtA_y =
        AtA[0][0]*(Atb[1]  * AtA[2][2] - AtA[1][2]*Atb[2]   ) -
        Atb[0] *  (AtA[1][0]*AtA[2][2] - AtA[1][2]*AtA[2][0]) +
        AtA[0][2]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]);
#else
      amrex::Real detAtA_y =
        (AtA[0][0]*Atb[1]) -
        (Atb[0] * AtA[1][0]);
#endif

      // Slope at centroid of (i,j,k)
      ys = detAtA_y / detAtA;

    } else {

      amrex::Real du_yl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
      amrex::Real du_yr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
      amrex::Real du_yc = 0.5*(state(i,j+1,k,n) - state(i,j-1,k,n));

      amrex::Real yslope = amrex::min(amrex::Math::abs(du_yl), amrex::Math::abs(du_yc), amrex::Math::abs(du_yr));
      yslope = (du_yr*du_yl > 0.0) ? yslope : 0.0;
      ys = (du_yc > 0.0) ? yslope : -yslope;
    }
#if (AMREX_SPACEDIM == 3)
    // Z direction
    if(flag(i,j,k  ).isSingleValued() or
      (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
      (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) 
    {

      amrex::Real detAtA_z =
        AtA[0][0]*(AtA[1][1]*Atb[2]    - Atb[1]   *AtA[1][2]) -
        AtA[0][1]*(AtA[1][0]*Atb[2]    - Atb[1]   *AtA[2][0]) +
        Atb[0]   *(AtA[1][0]*AtA[2][1] - AtA[1][1]*AtA[2][0]);

      zs = detAtA_z / detAtA;

    } else {

      amrex::Real du_zl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
      amrex::Real du_zr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
      amrex::Real du_zc = 0.5*(state(i,j,k+1,n) - state(i,j,k-1,n));

      amrex::Real zslope = amrex::min(amrex::Math::abs(du_zl), amrex::Math::abs(du_zc), amrex::Math::abs(du_zr));
      zslope = (du_zr*du_zl > 0.0) ? zslope : 0.0;
      zs = (du_zc > 0.0) ? zslope : -zslope;
    }
#endif

   return {AMREX_D_DECL(xs,ys,zs)};
}

// amrex_calc_slopes_extdir_eb calculates the slope in each coordinate direction using a 
// least squares linear fit to the 9 in 2D / 27 in 3D nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.
//
// No boundary conditions are used in this call.
//
// All the slopes are returned in one call.
//
// This differs from the previous call in that it over-rides the slope in the direction
// normal to a domain boundary when the boundary condition is imposed as external Dirichlet.
// We note that this is only correct in the case where the EB normals are parallel to the domain 
// boundary so that the transverse locations of the centroids in the cells being used are identical.

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_extdir_eb (int i, int j, int k, int n,
                             amrex::Array4<amrex::Real const> const& state,
                             amrex::Array4<amrex::Real const> const& ccent,
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                             AMREX_D_DECL(bool edlo_x, bool edlo_y, bool edlo_z),
                             AMREX_D_DECL(bool edhi_x, bool edhi_y, bool edhi_z), 
                             AMREX_D_DECL(int domlo_x, int domlo_y, int domlo_z),
                             AMREX_D_DECL(int domhi_x, int domhi_y, int domhi_z)) noexcept
{
    // First get EB-aware slope that doesn't know about extdir
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> eb_slope = amrex_calc_slopes_eb (i,j,k,n,state,ccent,flag);
    AMREX_D_TERM(amrex::Real xslope = eb_slope[0];,
                 amrex::Real yslope = eb_slope[1];,
                 amrex::Real zslope = eb_slope[2];);

    // Now correct only those cells which are not affected by EB but are by extdir
    if        (edlo_x and i == domlo_x) {
        if( !( flag(i  ,j,k).isSingleValued() or
              (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
              (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) ) 
        {
           amrex::Real dl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
           amrex::Real dr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
           amrex::Real dc = 
              (state(i+1,j,k,n)+3.0*state(i,j,k,n)-4.0*state(i-1,j,k,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           xslope = (dc > 0.0) ? slope : -slope;
        }

    } else if (edhi_x and i == domhi_x) {
        if( !( flag(i  ,j,k).isSingleValued() or
              (flag(i-1,j,k).isSingleValued() or not flag(i,j,k).isConnected(-1,0,0)) or
              (flag(i+1,j,k).isSingleValued() or not flag(i,j,k).isConnected( 1,0,0))) )
        {
           amrex::Real dl = 2.0*(state(i  ,j,k,n) - state(i-1,j,k,n));
           amrex::Real dr = 2.0*(state(i+1,j,k,n) - state(i  ,j,k,n));
           amrex::Real dc = 
              (4.0*state(i+1,j,k,n)-3.0*state(i,j,k,n)-state(i-1,j,k,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           xslope = (dc > 0.0) ? slope : -slope;
        }
    }

    if        (edlo_y and j == domlo_y) {
        if( !( flag(i,j  ,k).isSingleValued() or
              (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
              (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) )
        {
           amrex::Real dl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
           amrex::Real dr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
           amrex::Real dc = 
              (state(i,j+1,k,n)+3.0*state(i,j,k,n)-4.0*state(i,j-1,k,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           yslope = (dc > 0.0) ? slope : -slope;
        }
    } else if (edhi_y and j == domhi_y) {
        if( !(flag(i,j  ,k).isSingleValued() or
           (flag(i,j-1,k).isSingleValued() or not flag(i,j,k).isConnected(0,-1,0)) or
           (flag(i,j+1,k).isSingleValued() or not flag(i,j,k).isConnected(0, 1,0))) )
        {
           amrex::Real dl = 2.0*(state(i,j  ,k,n) - state(i,j-1,k,n));
           amrex::Real dr = 2.0*(state(i,j+1,k,n) - state(i,j  ,k,n));
           amrex::Real dc = 
              (4.0*state(i,j+1,k,n)-3.0*state(i,j,k,n)-state(i,j-1,k,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           yslope = (dc > 0.0) ? slope : -slope;
        }
    }
#if (AMREX_SPACEDIM == 3)
    if        (edlo_z and k == domlo_z) {
        if( !( flag(i,j,k  ).isSingleValued() or
              (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
              (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) )
        {
           amrex::Real dl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
           amrex::Real dr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
           amrex::Real dc = 
              (state(i,j,k+1,n)+3.0*state(i,j,k,n)-4.0*state(i,j,k-1,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           zslope = (dc > 0.0) ? slope : -slope;
        }
    } else if (edhi_z and k == domhi_z) {
        if( !( flag(i,j,k  ).isSingleValued() or
             (flag(i,j,k-1).isSingleValued() or not flag(i,j,k).isConnected(0,0,-1)) or
             (flag(i,j,k+1).isSingleValued() or not flag(i,j,k).isConnected(0,0, 1))) )
        {
           amrex::Real dl = 2.0*(state(i,j,k  ,n) - state(i,j,k-1,n));
           amrex::Real dr = 2.0*(state(i,j,k+1,n) - state(i,j,k  ,n));
           amrex::Real dc = 
              (4.0*state(i,j,k+1,n)-3.0*state(i,j,k,n)-state(i,j,k-1,n))/3.0;

           amrex::Real slope = amrex::min(amrex::Math::abs(dl),amrex::Math::abs(dc),amrex::Math::abs(dr));
           slope = (dr*dl > 0.0) ? slope : 0.0;
           zslope = (dc > 0.0) ? slope : -slope;
        }
    }
#endif

    return {AMREX_D_DECL(xslope,yslope,zslope)};
}

#endif
#endif

}
#endif
