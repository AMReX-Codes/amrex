#ifndef AMREX_HYDRO_EB_SLOPES_2D_K_H_
#define AMREX_HYDRO_EB_SLOPES_2D_K_H_

#include <AMReX_Config.H>
#include <AMReX_Slopes_K.H>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

// amrex_calc_slopes_eb calculates the slope in each coordinate direction using a
// least squares linear fit to the 8 nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.  Note that this routine
// is called either by amrex_calc_slopes_eb or amrex_calc_slopes_extdir_eb; in the former
// A is defined with the cell centroids; in the latter, the A values corresponding to values
// defined on faces use the face centroid location.
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb_given_A (int i, int j, int /*k*/, int n,
                              amrex::Real A[9][AMREX_SPACEDIM],
                              amrex::Array4<amrex::Real const> const& state,
                              amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    constexpr int dim_a = 9;
    amrex::Real du[dim_a];

    int ll=0;
    for(int jj(-1); jj<=1; jj++){
      for(int ii(-1); ii<=1; ii++){
        if( flag(i,j,0).isConnected(ii,jj,0) &&
            ! (ii==0 && jj==0)) {
          du[ll] = state(i+ii,j+jj,0,n) - state(i,j,0,n);
        } else {
          du[ll] = amrex::Real(0.0);
        }
        ll++;
      }
    }

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj)
    {
      for(int ii(0); ii<AMREX_SPACEDIM; ++ii){ // NOLINT(modernize-loop-convert)
        AtA[jj][ii] = amrex::Real(0.0);
      }
      Atb[jj] = amrex::Real(0.0);
    }

    for(int lc(0); lc<dim_a; ++lc)
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

// amrex_calc_slopes_eb_grown calculates the slope in each coordinate direction using a
// least squares linear fit to the (up to) 24 nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.  Note that this routine
// is called either by amrex_calc_slopes_eb_grown or amrex_calc_slopes_extdir_eb_grown; in the former
// A is defined with the cell centroids; in the latter, the A values corresponding to values
// defined on faces use the face centroid location.
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb_given_A_grown (int i, int j, int /*k*/, int n, int nx, int ny,
                                    amrex::Real A[25][AMREX_SPACEDIM],
                                    amrex::Array4<amrex::Real const> const& state,
                                    amrex::Array4<amrex::EBCellFlag const> const& flag) noexcept
{
    AMREX_ASSERT(nx == 1 || nx == 2);
    AMREX_ASSERT(ny == 1 || ny == 2);

    constexpr int dim_a = 25;

    amrex::Real AtA[AMREX_SPACEDIM][AMREX_SPACEDIM];
    amrex::Real Atb[AMREX_SPACEDIM];

    for(int jj(0); jj<AMREX_SPACEDIM; ++jj)
    {
        for(int ii(0); ii<AMREX_SPACEDIM; ++ii){ // NOLINT(modernize-loop-convert)
          AtA[jj][ii] = amrex::Real(0.0);
        }
        Atb[jj] = amrex::Real(0.0);
    }

    amrex::Real du;
    for(int jj(-ny); jj<=ny; jj++){
      for(int ii(-nx); ii<=nx; ii++){
          int lc = (jj+ny)*(2*nx+1) + ii+nx;
          if(!flag(i+ii,j+jj,0).isCovered() && !(ii==0 && jj==0)) {
              du = state(i+ii,j+jj,0,n) - state(i,j,0,n);
              Atb[0] += A[lc][0]*du;
              Atb[1] += A[lc][1]*du;
          }
          lc++;
      }
    }

    for(int lc(0); lc<dim_a; ++lc)
    {
        AtA[0][0] += A[lc][0]* A[lc][0];
        AtA[0][1] += A[lc][0]* A[lc][1];
        AtA[1][1] += A[lc][1]* A[lc][1];
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

//
// amrex_overwrite_with_regular_slopes calculates the slope in each coordinate direction
// with a standard non-EB slope calculation (that depends on max_order)
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
amrex_overwrite_with_regular_slopes(int i, int j, int k, int n,
                                    amrex::Real& xslope, amrex::Real& yslope,
                                    amrex::Array4<amrex::Real const> const& state,
                                    amrex::Array4<amrex::Real const> const& vfrac,
                                    int max_order) noexcept
{
    //
    // Here we over-write -- if possible -- with a stencil not using the EB stencil
    //
    AMREX_ALWAYS_ASSERT( max_order == 0 || max_order == 2 || max_order == 4);

    // We have enough cells in the x-direction to do 4th order slopes
    // centered on (i,j,k) with all values at cell centers
    if (max_order == 4 &&
            vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. &&
                                  vfrac(i+1,j,k) == 1. && vfrac(i+2,j,k) == 1.)
    {
        int order = 4;
        xslope = amrex_calc_xslope(i,j,k,n,order,state);
    }
    // We have enough cells in the x-direction to do 2nd order slopes
    // centered on (i,j,k) with all values at cell centers
    else if (max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i+1,j,k) == 1.)
    {
        int order = 2;
        xslope = amrex_calc_xslope(i,j,k,n,order,state);
    } else if (max_order == 0) {
        xslope = 0.;
    }

    // We have enough cells in the y-direction to do 4th order slopes
    // centered on (i,j,k) with all values at cell centers
    if (max_order == 4 &&
            vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. &&
                                  vfrac(i,j+1,k) == 1. && vfrac(i,j+2,k) == 1.)
    {
        int order = 4;
        yslope = amrex_calc_yslope(i,j,k,n,order,state);
    }
    // We have enough cells in the y-direction to do 2nd order slopes
    // centered on (i,j,k) with all values at cell centers
    else if (max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j+1,k) == 1.)
    {
        int order = 2;
        yslope = amrex_calc_yslope(i,j,k,n,order,state);
    } else if (max_order == 0) {
        yslope = 0.;
    }
}

// amrex_calc_slopes_eb calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
// OR
// 2) least squares linear fit to the at-most 8 nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.  Note that calc_slopes_eb
// defines the matrix A and passes this A to amrex_calc_slopes_eb_given_A.
//
// This routine assumes that there are no relevant hoextrap/extdir domain boundary conditions for this cell --
//     it does not test for them so this should not be called if such boundaries might be present
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb (int i, int j, int k, int n,
                      amrex::Array4<amrex::Real const> const& state,
                      amrex::Array4<amrex::Real const> const& ccent,
                      amrex::Array4<amrex::Real const> const& vfrac,
                      amrex::Array4<amrex::EBCellFlag const> const& flag,
                      int max_order) noexcept
{
    constexpr int dim_a = 9;
    amrex::Real A[dim_a][AMREX_SPACEDIM];

    int lc=0;
    int kk = 0;
    {
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){

            if( flag(i,j,k).isConnected(ii,jj,kk) &&
                ! (ii==0 && jj==0 && kk==0)) {

            // Not multiplying by dx to be consistent with how the
            // slope is stored. Also not including the global shift
            // wrt plo or i,j,k. We only need relative distance.

              A[lc][0] = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
              A[lc][1] = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);

            } else {

              A[lc][0] = amrex::Real(0.0);
              A[lc][1] = amrex::Real(0.0);
            }
            lc++;
          }
        }
    }

    const auto&  eb_slopes = amrex_calc_slopes_eb_given_A (i,j,k,n,A,state,flag);

    amrex::Real xslope = eb_slopes[0];
    amrex::Real yslope = eb_slopes[1];

    // This will over-write the values of xslope and yslope if appropriate
    amrex_overwrite_with_regular_slopes(i,j,k,n,xslope,yslope,state,vfrac,max_order);

    return {xslope,yslope};
}

// amrex_calc_slopes_eb_grown calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
// OR
// 2) least squares linear fit to the at-most 24 nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.  Note that calc_slopes_eb_grown
// defines the matrix A and passes this A to amrex_calc_slopes_eb_given_A_grown.
//
// This routine assumes that there are no relevant hoextrap/extdir domain boundary conditions for this cell --
//     it does not test for them so this should not be called if such boundaries might be present
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb_grown (int i, int j, int k, int n, int nx, int ny,
                            amrex::Array4<amrex::Real const> const& state,
                            amrex::Array4<amrex::Real const> const& ccent,
                            amrex::Array4<amrex::Real const> const& vfrac,
                            amrex::Array4<amrex::EBCellFlag const> const& flag,
                            int max_order) noexcept
{
    AMREX_ASSERT(nx == 1 || nx == 2);
    AMREX_ASSERT(ny == 1 || ny == 2);

    constexpr int dim_a = 25;
    amrex::Real A[dim_a][AMREX_SPACEDIM];

    int kk = 0;
    {
        // Make sure to zero all the entries in A (since the loop below may not cover all 25)
        int lc=0;
        for(int jj(-2); jj<=2; jj++){
          for(int ii(-2); ii<=2; ii++){
            A[lc][0] = amrex::Real(0.0);
            A[lc][1] = amrex::Real(0.0);
            lc++;
          }
        }

        lc=0;
        for(int jj(-ny); jj<=ny; jj++){
          for(int ii(-nx); ii<=nx; ii++){

            if(!flag(i+ii,j+jj,k+kk).isCovered() && !(ii==0 && jj==0 && kk==0))
            {
                A[lc][0] = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
                A[lc][1] = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
            }
            lc++;
          }
        }
    }

    const auto&  eb_slopes = amrex_calc_slopes_eb_given_A_grown (i,j,k,n,nx,ny,A,state,flag);

    amrex::Real xslope = eb_slopes[0];
    amrex::Real yslope = eb_slopes[1];

    // This will over-write the values of xslope and yslope if appropriate
    amrex_overwrite_with_regular_slopes(i,j,k,n,xslope,yslope,state,vfrac,max_order);

    return {xslope,yslope};
}

//
// amrex_overwrite_with_regular_slopes_extdir calculates the slope in each coordinate direction
// with a standard non-EB slope calculation (that depends on max_order)
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
amrex_overwrite_with_regular_slopes_extdir(int i, int j, int k, int n,
                                           amrex::Real& xslope, amrex::Real& yslope,
                                           amrex::Array4<amrex::Real const> const& state,
                                           amrex::Array4<amrex::Real const> const& vfrac,
                                           bool edlo_x, bool edlo_y, bool edhi_x, bool edhi_y,
                                           int domlo_x, int domlo_y, int domhi_x, int domhi_y,
                                           int max_order) noexcept
{
    //
    // Correct only those cells which are affected by extdir but not by EB:
    //    2) If all the cells are regular we use the "regular slope" in the extdir direction
    //

    // We have enough cells in the x-direction to do 4th order slopes
    // centered on (i,j,k) with all values at cell centers
    if (max_order == 4 &&
            vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. &&
                                  vfrac(i+1,j,k) == 1. && vfrac(i+2,j,k) == 1.)
    {
        int order = 4;
        xslope = amrex_calc_xslope_extdir(i,j,k,n,order,state,edlo_x,edhi_x,domlo_x,domhi_x);
    }
    // We have enough cells in the x-direction to do 2nd order slopes
    // centered on (i,j,k) with all values at cell centers
    else if (max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i+1,j,k) == 1.)
    {
        int order = 2;
        xslope = amrex_calc_xslope_extdir(i,j,k,n,order,state,edlo_x,edhi_x,domlo_x,domhi_x);
    } else if (max_order == 0) {
        xslope = 0.;
    }

    // We have enough cells in the y-direction to do 4th order slopes
    // centered on (i,j,k) with all values at cell centers
    if (max_order == 4 &&
            vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. &&
                                  vfrac(i,j+1,k) == 1. && vfrac(i,j+2,k) == 1.)
    {
        int order = 4;
        yslope = amrex_calc_yslope_extdir(i,j,k,n,order,state,edlo_y,edhi_y,domlo_y,domhi_y);
    }
    // We have enough cells in the y-direction to do 2nd order slopes
    // centered on (i,j,k) with all values at cell centers
    else if (max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j+1,k) == 1.)
    {
        int order = 2;
        yslope = amrex_calc_yslope_extdir(i,j,k,n,order,state,edlo_y,edhi_y,domlo_y,domhi_y);
    } else if (max_order == 0) {
        yslope = 0.;
    }
}

// amrex_calc_slopes_extdir_eb calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
//    (this stencil sees the extdir/hoextrap boundary condition if there is one)
// OR
// 2) least squares linear fit to the at-most 8 nearest neighbors, with the function
//    going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
//    where the data is assume to live, are the same as cell centers.
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_extdir_eb (int i, int j, int k, int n,
                             amrex::Array4<amrex::Real const> const& state,
                             amrex::Array4<amrex::Real const> const& ccent,
                             amrex::Array4<amrex::Real const> const& vfrac,
                             amrex::Array4<amrex::Real const> const& fcx,
                             amrex::Array4<amrex::Real const> const& fcy,
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                             bool edlo_x, bool edlo_y, bool edhi_x, bool edhi_y,
                             int domlo_x, int domlo_y, int domhi_x, int domhi_y,
                             int max_order) noexcept
{
    constexpr int dim_a = 9;

    auto xslope = amrex::Real(0.0);
    auto yslope = amrex::Real(0.0);

    // First get EB-aware slope that doesn't know about extdir
    bool needs_bdry_stencil = (edlo_x && i <= domlo_x) || (edhi_x && i >= domhi_x) ||
                              (edlo_y && j <= domlo_y) || (edhi_y && j >= domhi_y);

    //
    // This call does not have any knowledge of extdir / hoextrap boundary conditions
    //
    if (!needs_bdry_stencil)
    {
        // This returns slopes calculated with the regular 1-d approach if all cells in the stencil
        //      are regular.  If not, it uses the EB-aware least squares approach to fit a linear profile
        //      using the neighboring un-covered cells.
        const auto& slopes = amrex_calc_slopes_eb (i,j,k,n,state,ccent,vfrac,flag,max_order);
        return slopes;

    } else {

        amrex::Real A[dim_a][AMREX_SPACEDIM];

        int lc=0;
        int kk = 0;

        for(int jj(-1); jj<=1; jj++) {
        for(int ii(-1); ii<=1; ii++)
        {
            if( flag(i,j,k).isConnected(ii,jj,kk) &&
                ! (ii==0 && jj==0 && kk==0))
            {
                bool ilo_test = ( edlo_x && (i == domlo_x) && ii == -1);
                bool ihi_test = ( edhi_x && (i == domhi_x) && ii ==  1);

                bool jlo_test = ( edlo_y && (j == domlo_y) && jj == -1);
                bool jhi_test = ( edhi_y && (j == domhi_y) && jj ==  1);

                bool klo_test = false;
                bool khi_test = false;

                // These are the default values if no physical boundary
                A[lc][0] = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0);
                A[lc][1] = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1);
                // Do corrections for entire x-face
                if (ilo_test)
                {
                    if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                    {
                        A[lc][1] = amrex::Real(jj) + fcx(i   ,j+jj,k+kk,0);
                    }
                    A[lc][0] = amrex::Real(-0.5)                      ;
                } else if (ihi_test) {

                    if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                    {
                        A[lc][1] = amrex::Real(jj) + fcx(i+ii,j+jj,k+kk,0);
                    }
                    A[lc][0] = amrex::Real(0.5)                       ;
                }

                // Do corrections for entire y-face
                if (jlo_test)
                {
                    if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                    {
                        A[lc][0] = amrex::Real(ii) + fcy(i+ii,j   ,k+kk,0);
                    }
                    A[lc][1] = amrex::Real(-0.5)                      ;

                } else if (jhi_test) {

                    if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                    {
                        A[lc][0] = amrex::Real(ii) + fcy(i+ii,j+jj,k+kk,0);
                    }
                    A[lc][1] = amrex::Real(0.5)                       ;
                }

                A[lc][0] -= ccent(i,j,k,0);
                A[lc][1] -= ccent(i,j,k,1);

            } else {
                A[lc][0] = amrex::Real(0.0);
                A[lc][1] = amrex::Real(0.0);
            }
            lc++;
        }} // i,j

        const auto& slopes = amrex_calc_slopes_eb_given_A (i,j,k,n,A,state,flag);
        xslope = slopes[0];
        yslope = slopes[1];

        // This will over-write the values of xslope and yslope if appropriate
        amrex_overwrite_with_regular_slopes_extdir(i,j,k,n,xslope,yslope,state,vfrac,
                                                   edlo_x,edlo_y,edhi_x,edhi_y,
                                                   domlo_x,domlo_y,domhi_x,domhi_y,max_order);

        // Zero out slopes outside of an extdir (or hoextrap) boundary
        // TODO:  is this the right thing to do at a HOEXTRAP boundary??
        if ( (edlo_x && i < domlo_x) || (edhi_x && i > domhi_x)   ||
             (edlo_y && j < domlo_y) || (edhi_y && j > domhi_y) )
        {
            xslope = 0.;
            yslope = 0.;
        }
        return {xslope,yslope};
    }
}

// amrex_calc_slopes_extdir_eb_grown calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
//    (this stencil sees the extdir/hoextrap boundary condition if there is one)
// OR
// 2) least squares linear fit to the at-most 24  nearest neighbors, with the function
//    going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
//    where the data is assume to live, are the same as cell centers.
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_extdir_eb_grown (int i, int j, int k, int n, int nx, int ny,
                                   amrex::Array4<amrex::Real const> const& state,
                                   amrex::Array4<amrex::Real const> const& ccent,
                                   amrex::Array4<amrex::Real const> const& vfrac,
                                   amrex::Array4<amrex::Real const> const& fcx,
                                   amrex::Array4<amrex::Real const> const& fcy,
                                   amrex::Array4<amrex::EBCellFlag const> const& flag,
                                   bool edlo_x, bool edlo_y, bool edhi_x, bool edhi_y,
                                   int domlo_x, int domlo_y, int domhi_x, int domhi_y,
                                   int max_order) noexcept
{
    constexpr int dim_a = 25;

    auto xslope = amrex::Real(0.0);
    auto yslope = amrex::Real(0.0);

    // First get EB-aware slope that doesn't know about extdir
    bool needs_bdry_stencil = (edlo_x && i <= domlo_x) || (edhi_x && i >= domhi_x) ||
                              (edlo_y && j <= domlo_y) || (edhi_y && j >= domhi_y);

    //
    // This call does not have any knowledge of extdir / hoextrap boundary conditions
    //
    if (!needs_bdry_stencil)
    {
        // This returns slopes calculated with the regular 1-d approach if all cells in the stencil
        //      are regular.  If not, it uses the EB-aware least squares approach to fit a linear profile
        //      using the neighboring un-covered cells.
        const auto& slopes = amrex_calc_slopes_eb_grown (i,j,k,n,nx,ny,state,ccent,vfrac,flag,max_order);
        return slopes;

    } else {

        amrex::Real A[dim_a][AMREX_SPACEDIM];

        int lc=0;
        int kk = 0;

        for(int jj(-ny); jj<=ny; jj++) {
        for(int ii(-nx); ii<=nx; ii++)
        {
            if ( !flag(i+ii,j+jj,k).isCovered() && !(ii==0 && jj==0 && kk==0) )
            {
                bool ilo_test = ( edlo_x && (i == domlo_x) && ii == -1);
                bool ihi_test = ( edhi_x && (i == domhi_x) && ii ==  1);

                bool jlo_test = ( edlo_y && (j == domlo_y) && jj == -1);
                bool jhi_test = ( edhi_y && (j == domhi_y) && jj ==  1);

                bool klo_test = false;
                bool khi_test = false;

                // These are the default values if no physical boundary
                A[lc][0] = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0);
                A[lc][1] = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1);
                // Do corrections for entire x-face
                if (ilo_test)
                {
                    if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                    {
                        A[lc][1] = amrex::Real(jj) + fcx(i   ,j+jj,k+kk,0);
                    }
                    A[lc][0] = amrex::Real(-0.5)                      ;
                } else if (ihi_test) {

                    if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                    {
                        A[lc][1] = amrex::Real(jj) + fcx(i+ii,j+jj,k+kk,0);
                    }
                    A[lc][0] = amrex::Real(0.5)                       ;
                }

                // Do corrections for entire y-face
                if (jlo_test)
                {
                    if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                    {
                        A[lc][0] = amrex::Real(ii) + fcy(i+ii,j   ,k+kk,0);
                    }
                    A[lc][1] = amrex::Real(-0.5)                      ;

                } else if (jhi_test) {

                    if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                    {
                        A[lc][0] = amrex::Real(ii) + fcy(i+ii,j+jj,k+kk,0);
                    }
                    A[lc][1] = amrex::Real(0.5)                       ;
                }

                A[lc][0] -= ccent(i,j,k,0);
                A[lc][1] -= ccent(i,j,k,1);

            } else {
                A[lc][0] = amrex::Real(0.0);
                A[lc][1] = amrex::Real(0.0);
            }
            lc++;
        }} // i,j

        const auto& slopes = amrex_calc_slopes_eb_given_A_grown (i,j,k,n,nx,ny,A,state,flag);
        xslope = slopes[0];
        yslope = slopes[1];

        // This will over-write the values of xslope and yslope if appropriate
        amrex_overwrite_with_regular_slopes_extdir(i,j,k,n,xslope,yslope,state,vfrac,
                                                   edlo_x,edlo_y,edhi_x,edhi_y,
                                                   domlo_x,domlo_y,domhi_x,domhi_y,max_order);

        // Zero out slopes outside of an extdir (or hoextrap) boundary
        // TODO:  is this the right thing to do at a HOEXTRAP boundary??
        if ( (edlo_x && i < domlo_x) || (edhi_x && i > domhi_x)   ||
             (edlo_y && j < domlo_y) || (edhi_y && j > domhi_y) )
        {
            xslope = 0.;
            yslope = 0.;
        }
        return {xslope,yslope};
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
amrex_calc_alpha_stencil(amrex::Real q_hat, amrex::Real q_max,
                         amrex::Real q_min, amrex::Real state, amrex::Real alpha) noexcept
{
    using namespace amrex::literals;

    auto alpha_temp = 0.0_rt;
    auto small  = 1.0e-13_rt;

    if ((q_hat-state) > small) {
        alpha_temp = amrex::min(1.0_rt,(q_max-state)/(q_hat-state));
    } else if ((q_hat-state) < -small) {
        alpha_temp = amrex::min(1.0_rt,(q_min-state)/(q_hat-state));
    } else {
        alpha_temp = 1.0_rt;
    }
    return amrex::min(alpha, alpha_temp);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_alpha_limiter(int i, int j, int k, int n,
                         amrex::Array4<amrex::Real const> const& state,
                         amrex::Array4<amrex::EBCellFlag const> const& flag,
                         const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>& slopes,
                         amrex::Array4<amrex::Real const> const& fcx,
                         amrex::Array4<amrex::Real const> const& fcy,
                         amrex::Array4<amrex::Real const> const& ccent) noexcept
{
    amrex::Real alpha = 2.0;

    int cuts_x = 0;
    int cuts_y = 0;

    // Compute how many cut or regular faces do we have in 3x3 block
    int kk = 0;
    for(int jj(-1); jj<=1; jj++){
      for(int ii(-1); ii<=1; ii++){
        if( flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0))
        {
            if ((ii==-1 || ii==1) && jj==0) { cuts_x++; }
            if ((jj==-1 || jj==1) && ii==0) { cuts_y++; }
        }
      }
    }

    amrex::Real xc = ccent(i,j,k,0); // centroid of cell (i,j,k)
    amrex::Real yc = ccent(i,j,k,1);

    //Reconstruct values at the face centroids and compute the limiter
    if(flag(i,j,k).isConnected(0,1,0)) {
        amrex::Real xf = fcy(i,j+1,k,0); // local (x,z) of centroid of y-face we are extrapolating to

        amrex::Real delta_x = xf  - xc;
        amrex::Real delta_y = amrex::Real(0.5) - yc;

        amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0] + delta_y * slopes[1];

        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);

        // Compute max and min values in a 3x2 stencil
        for(int jj(0); jj<=1; jj++){
            for(int ii(-1); ii<=1; ii++){
                if ( flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0) ) {
                    if (state(i+ii,j+jj,k+kk,n) > q_max) { q_max = state(i+ii,j+jj,k+kk,n); }
                    if (state(i+ii,j+jj,k+kk,n) < q_min) { q_min = state(i+ii,j+jj,k+kk,n); }
                }
            }
        }

        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if (flag(i,j,k).isConnected(0,-1,0)){
        amrex::Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to

        amrex::Real delta_x = xf  - xc;
        amrex::Real delta_y = amrex::Real(0.5) + yc;

        amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0] - delta_y * slopes[1];

        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);

        // Compute max and min values in a 3x2 stencil
        for(int jj(-1); jj<=0; jj++){
            for(int ii(-1); ii<=1; ii++){
                if ( flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0) ) {
                    if (state(i+ii,j+jj,k+kk,n) > q_max) { q_max = state(i+ii,j+jj,k+kk,n); }
                    if (state(i+ii,j+jj,k+kk,n) < q_min) { q_min = state(i+ii,j+jj,k+kk,n); }
                }
            }
        }

        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if(flag(i,j,k).isConnected(1,0,0)) {
        amrex::Real yf = fcx(i+1,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to

        amrex::Real delta_x = amrex::Real(0.5) - xc;
        amrex::Real delta_y = yf  - yc;

        amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0] + delta_y * slopes[1];

        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);

        for(int jj(-1); jj<=1; jj++){
            for(int ii(0); ii<=1; ii++){
                if ( flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0) ) {
                    if (state(i+ii,j+jj,k+kk,n) > q_max) { q_max = state(i+ii,j+jj,k+kk,n); }
                    if (state(i+ii,j+jj,k+kk,n) < q_min) { q_min = state(i+ii,j+jj,k+kk,n); }
                }
            }
        }

        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if(flag(i,j,k).isConnected(-1,0,0)) {
        amrex::Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to

        amrex::Real delta_x = amrex::Real(0.5) + xc;
        amrex::Real delta_y = yf  - yc;

        amrex::Real q_hat = state(i,j,k,n) - delta_x * slopes[0] + delta_y * slopes[1];

        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);

        for(int jj(-1); jj<=1; jj++){
            for(int ii(-1); ii<=0; ii++){
                if( flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0)) {
                    if (state(i+ii,j+jj,k+kk,n) > q_max) { q_max = state(i+ii,j+jj,k+kk,n); }
                    if (state(i+ii,j+jj,k+kk,n) < q_min) { q_min = state(i+ii,j+jj,k+kk,n); }
                }
            }
        }
        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }

    amrex::Real xalpha = alpha;
    amrex::Real yalpha = alpha;

    //Zeroing out the slopes in the direction where a covered face exists.
    if (cuts_x<2) { xalpha = 0; }
    if (cuts_y<2) { yalpha = 0; }

    return {xalpha,yalpha};
}

//amrex_lim_slopes_eb computes the slopes calling amrex_calc_slopes_eb, and then each slope component
//is multiplied by a limiter based on the work of Barth-Jespersen.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_lim_slopes_eb (int i, int j, int k, int n,
                     amrex::Array4<amrex::Real const> const& state,
                     amrex::Array4<amrex::Real const> const& ccent,
                     amrex::Array4<amrex::Real const> const& vfrac,
                     amrex::Array4<amrex::Real const> const& fcx,
                     amrex::Array4<amrex::Real const> const& fcy,
                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                     int max_order) noexcept
{

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> alpha_lim;

    slopes = amrex_calc_slopes_eb(i,j,k,n,state,ccent,vfrac,flag,max_order);

    alpha_lim = amrex_calc_alpha_limiter(i,j,k,n,state,flag,slopes,fcx,fcy,ccent);

    // Setting limiter to 1 for stencils that just consists of non-EB cells because
    // amrex_calc_slopes_eb routine will call the slope routine for non-EB stencils that has already a limiter
    if ( max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i+1,j,k) == 1.) {
        alpha_lim[0] = 1.0;
    }

    if ( max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j+1,k) == 1.) {
        alpha_lim[1] = 1.0;
    }

    return {alpha_lim[0]*slopes[0],alpha_lim[1]*slopes[1]};
}

//amrex_lim_slopes_extdir_eb computes the slopes calling amrex_calc_slopes_extdir_eb, and then each slope component
//is multiplied by a limiter based on the work of Barth-Jespersen.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_lim_slopes_extdir_eb (int i, int j, int k, int n,
                            amrex::Array4<amrex::Real const> const& state,
                            amrex::Array4<amrex::Real const> const& ccent,
                            amrex::Array4<amrex::Real const> const& vfrac,
                            amrex::Array4<amrex::Real const> const& fcx,
                            amrex::Array4<amrex::Real const> const& fcy,
                            amrex::Array4<amrex::EBCellFlag const> const& flag,
                            bool edlo_x, bool edlo_y, bool edhi_x, bool edhi_y,
                            int domlo_x, int domlo_y, int domhi_x, int domhi_y,
                            int max_order) noexcept
{

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> alpha_lim;

    slopes = amrex_calc_slopes_extdir_eb(i,j,k,n,state,ccent,vfrac,fcx,fcy,flag,
                                         edlo_x,edlo_y,edhi_x,edhi_y,
                                         domlo_x,domlo_y,domhi_x,domhi_y,max_order);
    alpha_lim = amrex_calc_alpha_limiter(i,j,k,n,state,flag,slopes,fcx,fcy,ccent);

    // Setting limiter to 1 for stencils that just consists of non-EB cells because
    // amrex_calc_slopes_extdir_eb routine will call the slope routine for non-EB stencils that has already a limiter
    if ( max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i-1,j,k) == 1. && vfrac(i+1,j,k) == 1.) {
        alpha_lim[0] = 1.0;
    }

    if ( max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i,j-1,k) == 1. && vfrac(i,j+1,k) == 1.) {
        alpha_lim[1] = 1.0;
    }

    return {alpha_lim[0]*slopes[0],alpha_lim[1]*slopes[1]};
}

#endif
