#ifndef AMREX_HYDRO_EB_SLOPES_3D_K_H_
#define AMREX_HYDRO_EB_SLOPES_3D_K_H_

#include <AMReX_Config.H>
#include <AMReX_Slopes_K.H>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

struct EBSlopeNormalEquation
{
    Real a00;
    Real a01;
    Real a02;
    Real a11;
    Real a12;
    Real a22;
    Real b0;
    Real b1;
    Real b2;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    EBSlopeNormalEquation() noexcept
        : a00(0.0_rt), a01(0.0_rt), a02(0.0_rt),
          a11(0.0_rt), a12(0.0_rt), a22(0.0_rt),
          b0(0.0_rt), b1(0.0_rt), b2(0.0_rt)
    {}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void add (Real offset_x, Real offset_y, Real offset_z, Real du) noexcept
    {
        a00 += offset_x*offset_x;
        a01 += offset_x*offset_y;
        a02 += offset_x*offset_z;
        a11 += offset_y*offset_y;
        a12 += offset_y*offset_z;
        a22 += offset_z*offset_z;

        b0 += offset_x*du;
        b1 += offset_y*du;
        b2 += offset_z*du;
    }

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    GpuArray<Real,AMREX_SPACEDIM> solve () const noexcept
    {
        Real detAtA =
          a00*(a11*a22 - a12*a12) -
          a01*(a01*a22 - a12*a02) +
          a02*(a01*a12 - a11*a02);

        Real detAtA_x =
            b0 *(a11*a22 - a12*a12) -
            a01*(b1*a22 - a12*b2) +
            a02*(b1*a12 - a11*b2);

        // Slope at centroid of (i,j,k)
        Real xs = detAtA_x / detAtA;

        Real detAtA_y =
            a00*(b1*a22 - a12*b2) -
            b0 *(a01*a22 - a12*a02) +
            a02*(a01*b2 - b1*a02);

        // Slope at centroid of (i,j,k)
        Real ys = detAtA_y / detAtA;

        Real detAtA_z =
            a00*(a11*b2 - b1*a12) -
            a01*(a01*b2 - b1*a02) +
            b0 *(a01*a12 - a11*a02);

        // Slope at centroid of (i,j,k)
        Real zs = detAtA_z / detAtA;

        return {xs,ys,zs};
    }
};

}

//
// amrex_overwrite_with_regular_slopes calculates the slope in each coordinate direction
// with a standard non-EB slope calculation (that depends on max_order)
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
amrex_overwrite_with_regular_slopes(int i, int j, int k, int n,
                                    amrex::Real& xslope, amrex::Real& yslope,
                                    amrex::Real& zslope,
                                    amrex::Array4<amrex::Real const> const& state,
                                    amrex::Array4<amrex::Real const> const& vfrac,
                                    int max_order)
{
    //
    // Here we over-write -- if possible -- with a stencil not using the EB stencil
    //

    if (max_order == 0) {
        xslope = 0.;
        yslope = 0.;
        zslope = 0.;
        return;
    }

    if (vfrac(i,j,k) != 1.) { return; }

    // We have enough cells in the x-direction to do 4th order slopes
    // centered on (i,j,k) with all values at cell centers
    if (max_order == 4 &&
            vfrac(i-1,j,k) == 1. && vfrac(i-2,j,k) == 1. &&
            vfrac(i+1,j,k) == 1. && vfrac(i+2,j,k) == 1.)
    {
        int order = 4;
        xslope = amrex_calc_xslope(i,j,k,n,order,state);
    }
    // We have enough cells in the x-direction to do 2nd order slopes
    // centered on (i,j,k) with all values at cell centers
    else if (vfrac(i-1,j,k) == 1. && vfrac(i+1,j,k) == 1.)
    {
        int order = 2;
        xslope = amrex_calc_xslope(i,j,k,n,order,state);
    }

    // We have enough cells in the y-direction to do 4th order slopes
    // centered on (i,j,k) with all values at cell centers
    if (max_order == 4 &&
            vfrac(i,j-1,k) == 1. && vfrac(i,j-2,k) == 1. &&
            vfrac(i,j+1,k) == 1. && vfrac(i,j+2,k) == 1.)
    {
        int order = 4;
        yslope = amrex_calc_yslope(i,j,k,n,order,state);
    }
    // We have enough cells in the y-direction to do 2nd order slopes
    // centered on (i,j,k) with all values at cell centers
    else if (vfrac(i,j-1,k) == 1. && vfrac(i,j+1,k) == 1.)
    {
        int order = 2;
        yslope = amrex_calc_yslope(i,j,k,n,order,state);
    }

#if (AMREX_SPACEDIM == 3)
    // We have enough cells in the z-direction to do 4th order slopes
    // centered on (i,j,k) with all values at cell centers
    if (max_order == 4 &&
            vfrac(i,j,k-1) == 1. && vfrac(i,j,k-2) == 1. &&
            vfrac(i,j,k+1) == 1. && vfrac(i,j,k+2) == 1.)
    {
        int order = 4;
        zslope = amrex_calc_zslope(i,j,k,n,order,state);
    }
    // We have enough cells in the z-direction to do 2nd order slopes
    // centered on (i,j,k) with all values at cell centers
    else if (vfrac(i,j,k-1) == 1. && vfrac(i,j,k+1) == 1.)
    {
        int order = 2;
        zslope = amrex_calc_zslope(i,j,k,n,order,state);
    }
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool
amrex_has_all_regular_slope_stencils (
    int i, int j, int k,
    amrex::Array4<amrex::Real const> const& vfrac) noexcept
{
    return vfrac(i,j,k) == 1. &&
           vfrac(i-1,j,k) == 1. && vfrac(i+1,j,k) == 1. &&
           vfrac(i,j-1,k) == 1. && vfrac(i,j+1,k) == 1. &&
           vfrac(i,j,k-1) == 1. && vfrac(i,j,k+1) == 1.;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_regular_slopes (int i, int j, int k, int n,
                           amrex::Array4<amrex::Real const> const& state,
                           amrex::Array4<amrex::Real const> const& vfrac,
                           int max_order) noexcept
{
    int order = 2;
    if (max_order == 4 && vfrac(i-2,j,k) == 1. && vfrac(i+2,j,k) == 1.) {
        order = 4;
    }
    amrex::Real const xslope = amrex_calc_xslope(i,j,k,n,order,state);

    order = 2;
    if (max_order == 4 && vfrac(i,j-2,k) == 1. && vfrac(i,j+2,k) == 1.) {
        order = 4;
    }
    amrex::Real const yslope = amrex_calc_yslope(i,j,k,n,order,state);

    order = 2;
    if (max_order == 4 && vfrac(i,j,k-2) == 1. && vfrac(i,j,k+2) == 1.) {
        order = 4;
    }
    amrex::Real const zslope = amrex_calc_zslope(i,j,k,n,order,state);

    return {xslope,yslope,zslope};
}

// amrex_calc_slopes_eb calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
// OR
// 2) least squares linear fit to the at-most 26 nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.
//
// This routine assumes that there are no relevant hoextrap/extdir domain boundary conditions for this cell --
//     it does not test for them so this should not be called if such boundaries might be present
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb (int i, int j, int k, int n,
                      amrex::Array4<amrex::Real const> const& state,
                      amrex::Array4<amrex::Real const> const& ccent,
                      amrex::Array4<amrex::Real const> const& vfrac,
                      amrex::Array4<amrex::EBCellFlag const> const& flag,
                      int max_order) noexcept
{
    if (max_order == 0) { return {0.,0.,0.}; }

    if (amrex_has_all_regular_slope_stencils(i,j,k,vfrac)) {
        return amrex_calc_regular_slopes(i,j,k,n,state,vfrac,max_order);
    }

    amrex::EBSlopeNormalEquation normal_eq;

    // Prevent loop unrolling to reduce register pressure
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for(int kk(-1); kk<=1; kk++) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for(int jj(-1); jj<=1; jj++) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for(int ii(-1); ii<=1; ii++)
    {
        if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0))
        {
            amrex::Real offset_x = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
            amrex::Real offset_y = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
            amrex::Real offset_z = amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2);
            // Do not read state for an excluded zero row: zero times NaN is still NaN.
            if (offset_x != amrex::Real(0) || offset_y != amrex::Real(0) ||
                offset_z != amrex::Real(0)) {
                amrex::Real du = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);
                normal_eq.add(offset_x,offset_y,offset_z,du);
            }
        }
    }}}

    //
    // These slopes use the EB stencil without testing whether it is actually needed
    //
    const auto slopes = normal_eq.solve();
    amrex::Real xslope = slopes[0];
    amrex::Real yslope = slopes[1];
    amrex::Real zslope = slopes[2];

    // This will over-write the values of xslope, yslope and/or zslope if appropriate
    amrex_overwrite_with_regular_slopes(i,j,k,n,xslope,yslope,zslope,state,vfrac,max_order);

    return {xslope,yslope,zslope};
}

// amrex_calc_slopes_eb_grown calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
// OR
// 2) least squares linear fit to the at-most 124 nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.
//
// This routine assumes that there are no relevant hoextrap/extdir domain boundary conditions for this cell --
//     it does not test for them so this should not be called if such boundaries might be present
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_eb_grown (int i, int j, int k, int n, int nx, int ny, int nz,
                            amrex::Array4<amrex::Real const> const& state,
                            amrex::Array4<amrex::Real const> const& ccent,
                            amrex::Array4<amrex::Real const> const& vfrac,
                            amrex::Array4<amrex::EBCellFlag const> const& flag,
                            int max_order) noexcept
{
    AMREX_ASSERT(nx == 1 || nx == 2);
    AMREX_ASSERT(ny == 1 || ny == 2);
    AMREX_ASSERT(nz == 1 || nz == 2);

    if (max_order == 0) { return {0.,0.,0.}; }

    if (amrex_has_all_regular_slope_stencils(i,j,k,vfrac)) {
        return amrex_calc_regular_slopes(i,j,k,n,state,vfrac,max_order);
    }

    amrex::EBSlopeNormalEquation normal_eq;

    // Prevent loop unrolling to reduce register pressure
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for(int kk(-nz); kk<=nz; kk++) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for(int jj(-ny); jj<=ny; jj++) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
    for(int ii(-nx); ii<=nx; ii++)
    {
        if (!flag(i+ii,j+jj,k+kk).isCovered() && !(ii==0 && jj==0 && kk==0))
        {
            amrex::Real offset_x = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0) - ccent(i,j,k,0);
            amrex::Real offset_y = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1) - ccent(i,j,k,1);
            amrex::Real offset_z = amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2) - ccent(i,j,k,2);
            // Do not read state for an excluded zero row: zero times NaN is still NaN.
            if (offset_x != amrex::Real(0) || offset_y != amrex::Real(0) ||
                offset_z != amrex::Real(0)) {
                amrex::Real du = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);
                normal_eq.add(offset_x,offset_y,offset_z,du);
            }
        }
    }}}
    //
    // These slopes use the EB stencil without testing whether it is actually needed
    //
    const auto slopes = normal_eq.solve();
    amrex::Real xslope = slopes[0];
    amrex::Real yslope = slopes[1];
    amrex::Real zslope = slopes[2];

    // This will over-write the values of xslope, yslope and/or zslope if appropriate
    amrex_overwrite_with_regular_slopes(i,j,k,n,xslope,yslope,zslope,state,vfrac,max_order);

    return {xslope,yslope,zslope};
}

//
// amrex_overwrite_with_regular_slopes_extdir calculates the slope in each coordinate direction
// with a standard non-EB slope calculation (that depends on max_order)
//
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
amrex_overwrite_with_regular_slopes_extdir(int i, int j, int k, int n,
                                           amrex::Real& xslope, amrex::Real& yslope,
                                           amrex::Real& zslope,
                                           amrex::Array4<amrex::Real const> const& state,
                                           amrex::Array4<amrex::Real const> const& vfrac,
                                           bool edlo_x, bool edlo_y, bool edlo_z,
                                           bool edhi_x, bool edhi_y, bool edhi_z,
                                           int domlo_x, int domlo_y, int domlo_z,
                                           int domhi_x, int domhi_y, int domhi_z,
                                           int max_order) noexcept
{
    //
    // Correct only those cells which are affected by extdir but not by EB:
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

    // We have enough cells in the z-direction to do 4th order slopes
    // centered on (i,j,k) with all values at cell centers
    if (max_order == 4 &&
            vfrac(i,j,k) == 1. && vfrac(i,j,k-1) == 1. && vfrac(i,j,k-2) == 1. &&
                                  vfrac(i,j,k+1) == 1. && vfrac(i,j,k+2) == 1.)
    {
        int order = 4;
        zslope = amrex_calc_zslope_extdir(i,j,k,n,order,state,edlo_z,edhi_z,domlo_z,domhi_z);
    }
    // We have enough cells in the z-direction to do 2nd order slopes
    // centered on (i,j,k) with all values at cell centers
    else if (max_order > 0 && vfrac(i,j,k) == 1. && vfrac(i,j,k-1) == 1. && vfrac(i,j,k+1) == 1.)
    {
        int order = 2;
        zslope = amrex_calc_zslope_extdir(i,j,k,n,order,state,edlo_z,edhi_z,domlo_z,domhi_z);
    } else if (max_order == 0) {
        zslope = 0.;
    }
}

// amrex_calc_slopes_extdir_eb calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
//    (this stencil sees the extdir/hoextrap boundary condition if there is one)
// OR
// 2) least squares linear fit to the at-most 26 nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.
//
// All the slopes are returned in one call.
//

template <bool BoundaryOnly = false>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_extdir_eb (int i, int j, int k, int n,
                             amrex::Array4<amrex::Real const> const& state,
                             amrex::Array4<amrex::Real const> const& ccent,
                             amrex::Array4<amrex::Real const> const& vfrac,
                             amrex::Array4<amrex::Real const> const& fcx,
                             amrex::Array4<amrex::Real const> const& fcy,
                             amrex::Array4<amrex::Real const> const& fcz,
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                             bool edlo_x, bool edlo_y, bool edlo_z,
                             bool edhi_x, bool edhi_y, bool edhi_z,
                             int domlo_x, int domlo_y, int domlo_z,
                             int domhi_x, int domhi_y, int domhi_z,
                             int max_order) noexcept
{
    if constexpr (!BoundaryOnly) {
        // Include cells one in from the boundary because fourth-order one-sided
        // stencils reach domlo+1 / domhi-1.
        bool const needs_bdry_stencil =
            (edlo_x && i <= domlo_x+1) || (edhi_x && i >= domhi_x-1) ||
            (edlo_y && j <= domlo_y+1) || (edhi_y && j >= domhi_y-1) ||
            (edlo_z && k <= domlo_z+1) || (edhi_z && k >= domhi_z-1);

        if (!needs_bdry_stencil) {
            return amrex_calc_slopes_eb(i,j,k,n,state,ccent,vfrac,flag,max_order);
        }
    }

    amrex::EBSlopeNormalEquation normal_eq;

        // Prevent loop unrolling to reduce register pressure
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
        for(int kk(-1); kk<=1; kk++) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
        for(int jj(-1); jj<=1; jj++) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
        for(int ii(-1); ii<=1; ii++)
        {
            if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0))
            {
                bool ilo_test = ( edlo_x && (i == domlo_x) && ii == -1);
                bool ihi_test = ( edhi_x && (i == domhi_x) && ii ==  1);

                bool jlo_test = ( edlo_y && (j == domlo_y) && jj == -1);
                bool jhi_test = ( edhi_y && (j == domhi_y) && jj ==  1);

                bool klo_test = ( edlo_z && (k == domlo_z) && kk == -1);
                bool khi_test = ( edhi_z && (k == domhi_z) && kk ==  1);

                // These are the default values if no physical boundary
                amrex::Real offset_x = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0);
                amrex::Real offset_y = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1);
                amrex::Real offset_z = amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2);
                // Do corrections for entire x-face
                if (ilo_test)
                {
                    if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                    {
                        offset_y = amrex::Real(jj) + fcx(i   ,j+jj,k+kk,0);
                        offset_z = amrex::Real(kk) + fcx(i   ,j+jj,k+kk,1);
                    }
                    offset_x = -amrex::Real(0.5);
                } else if (ihi_test) {

                    if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                    {
                        offset_y = amrex::Real(jj) + fcx(i+ii,j+jj,k+kk,0);
                        offset_z = amrex::Real(kk) + fcx(i+ii,j+jj,k+kk,1);
                    }
                    offset_x = amrex::Real(0.5);
                }

                // Do corrections for entire y-face
                if (jlo_test)
                {
                    if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                    {
                        offset_x = amrex::Real(ii) + fcy(i+ii,j   ,k+kk,0);
                        offset_z = amrex::Real(kk) + fcy(i+ii,j   ,k+kk,1);
                    }
                    offset_y = -amrex::Real(0.5);

                } else if (jhi_test) {

                    if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                    {
                        offset_x = amrex::Real(ii) + fcy(i+ii,j+jj,k+kk,0);
                        offset_z = amrex::Real(kk) + fcy(i+ii,j+jj,k+kk,1);
                    }
                    offset_y = amrex::Real(0.5);
                }

                // Do corrections for entire z-face
                if (klo_test)
                {
                    if (!ilo_test && !ihi_test && !jlo_test && !jhi_test)
                    {
                        offset_x = amrex::Real(ii) + fcz(i+ii,j+jj,k   ,0);
                        offset_y = amrex::Real(jj) + fcz(i+ii,j+jj,k   ,1);
                    }
                    offset_z = -amrex::Real(0.5);

                } else if (khi_test) {
                    if (!ilo_test && !ihi_test && !jlo_test && !jhi_test)
                    {
                        offset_x = amrex::Real(ii) + fcz(i+ii,j+jj,k+kk,0);
                        offset_y = amrex::Real(jj) + fcz(i+ii,j+jj,k+kk,1);
                    }
                    offset_z = amrex::Real(0.5);
                }

                offset_x -= ccent(i,j,k,0);
                offset_y -= ccent(i,j,k,1);
                offset_z -= ccent(i,j,k,2);

                // Do not read state for an excluded zero row: zero times NaN is still NaN.
                if (offset_x != amrex::Real(0) || offset_y != amrex::Real(0) ||
                    offset_z != amrex::Real(0)) {
                    amrex::Real du = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);
                    normal_eq.add(offset_x,offset_y,offset_z,du);
                }
            }
        }}}

    const auto slopes = normal_eq.solve();
    auto xslope = slopes[0];
    auto yslope = slopes[1];
    auto zslope = slopes[2];

        // This will over-write the values of xslope and yslope if appropriate
    amrex_overwrite_with_regular_slopes_extdir(i,j,k,n,xslope,yslope,zslope,state,vfrac,
                                               edlo_x,edlo_y,edlo_z,edhi_x,edhi_y,edhi_z,
                                               domlo_x,domlo_y,domlo_z,domhi_x,domhi_y,domhi_z,
                                               max_order);

    // Zero out slopes outside of an extdir (or hoextrap) boundary
    // TODO:  is this the right thing to do at a HOEXTRAP boundary??
    if ( (edlo_x && i < domlo_x) || (edhi_x && i > domhi_x) ||
         (edlo_y && j < domlo_y) || (edhi_y && j > domhi_y) ||
         (edlo_z && k < domlo_z) || (edhi_z && k > domhi_z) )
    {
          xslope = 0.; yslope = 0.; zslope = 0.;
    }

    return {xslope,yslope,zslope};
}

// amrex_calc_slopes_extdir_eb_grown calculates the slope in each coordinate direction using a
// 1) standard limited slope if all three cells in the stencil are regular cells
//    (this stencil sees the extdir/hoextrap boundary condition if there is one)
// OR
// 2) least squares linear fit to the at-most 124 nearest neighbors, with the function
// going through the centroid of cell(i,j,k).  This does not assume that the cell centroids,
// where the data is assume to live, are the same as cell centers.
//
template <bool BoundaryOnly = false>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
amrex_calc_slopes_extdir_eb_grown (int i, int j, int k, int n,
                                   int nx, int ny, int nz,
                                   amrex::Array4<amrex::Real const> const& state,
                                   amrex::Array4<amrex::Real const> const& ccent,
                                   amrex::Array4<amrex::Real const> const& vfrac,
                                   amrex::Array4<amrex::Real const> const& fcx,
                                   amrex::Array4<amrex::Real const> const& fcy,
                                   amrex::Array4<amrex::Real const> const& fcz,
                                   amrex::Array4<amrex::EBCellFlag const> const& flag,
                                   bool edlo_x, bool edlo_y, bool edlo_z,
                                   bool edhi_x, bool edhi_y, bool edhi_z,
                                   int domlo_x, int domlo_y, int domlo_z,
                                   int domhi_x, int domhi_y, int domhi_z,
                                   int max_order) noexcept
{
    if constexpr (!BoundaryOnly) {
        // Include cells one in from the boundary because fourth-order one-sided
        // and grown least-squares stencils can reach that far.
        bool const needs_bdry_stencil =
            (edlo_x && i <= domlo_x+1) || (edhi_x && i >= domhi_x-1) ||
            (edlo_y && j <= domlo_y+1) || (edhi_y && j >= domhi_y-1) ||
            (edlo_z && k <= domlo_z+1) || (edhi_z && k >= domhi_z-1);

        if (!needs_bdry_stencil) {
            return amrex_calc_slopes_eb_grown(
                i,j,k,n,nx,ny,nz,state,ccent,vfrac,flag,max_order);
        }
    }

    amrex::EBSlopeNormalEquation normal_eq;

        // Prevent loop unrolling to reduce register pressure
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
        for(int kk(-nz); kk<=nz; kk++) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
        for(int jj(-ny); jj<=ny; jj++) {
#ifdef AMREX_USE_GPU
#pragma unroll 1
#endif
        for(int ii(-nx); ii<=nx; ii++)
        {
            // Cells more than one cell outside an extdir/hoextrap boundary hold no usable data,
            //    so we exclude them from the fit just as we exclude covered cells
            bool outside = (edlo_x && i+ii < domlo_x-1) || (edhi_x && i+ii > domhi_x+1) ||
                           (edlo_y && j+jj < domlo_y-1) || (edhi_y && j+jj > domhi_y+1) ||
                           (edlo_z && k+kk < domlo_z-1) || (edhi_z && k+kk > domhi_z+1);

            if (!outside && !flag(i+ii,j+jj,k+kk).isCovered() && !(ii==0 && jj==0 && kk==0))
            {
                // The value in the first ghost cell outside an extdir/hoextrap boundary lives on
                //    the domain face, not at the cell centroid.  With the grown (n == 2) stencil
                //    that ghost cell can be reached from a cell one in from the boundary, so we
                //    key these tests on the neighbor index rather than on the offset alone.
                bool ilo_test = ( edlo_x && (i+ii == domlo_x-1) );
                bool ihi_test = ( edhi_x && (i+ii == domhi_x+1) );

                bool jlo_test = ( edlo_y && (j+jj == domlo_y-1) );
                bool jhi_test = ( edhi_y && (j+jj == domhi_y+1) );

                bool klo_test = ( edlo_z && (k+kk == domlo_z-1) );
                bool khi_test = ( edhi_z && (k+kk == domhi_z+1) );

                // These are the default values if no physical boundary
                amrex::Real offset_x = amrex::Real(ii) + ccent(i+ii,j+jj,k+kk,0);
                amrex::Real offset_y = amrex::Real(jj) + ccent(i+ii,j+jj,k+kk,1);
                amrex::Real offset_z = amrex::Real(kk) + ccent(i+ii,j+jj,k+kk,2);
                // Do corrections for entire x-face
                if (ilo_test)
                {
                    if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                    {
                        offset_y = amrex::Real(jj) + fcx(domlo_x,j+jj,k+kk,0);
                        offset_z = amrex::Real(kk) + fcx(domlo_x,j+jj,k+kk,1);
                    }
                    offset_x = amrex::Real(domlo_x-i) - amrex::Real(0.5);
                } else if (ihi_test) {

                    if (!jlo_test && !jhi_test && !klo_test && !khi_test)
                    {
                        offset_y = amrex::Real(jj) + fcx(domhi_x+1,j+jj,k+kk,0);
                        offset_z = amrex::Real(kk) + fcx(domhi_x+1,j+jj,k+kk,1);
                    }
                    offset_x = amrex::Real(domhi_x-i) + amrex::Real(0.5);
                }

                // Do corrections for entire y-face
                if (jlo_test)
                {
                    if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                    {
                        offset_x = amrex::Real(ii) + fcy(i+ii,domlo_y,k+kk,0);
                        offset_z = amrex::Real(kk) + fcy(i+ii,domlo_y,k+kk,1);
                    }
                    offset_y = amrex::Real(domlo_y-j) - amrex::Real(0.5);

                } else if (jhi_test) {

                    if (!ilo_test && !ihi_test && !klo_test && !khi_test)
                    {
                        offset_x = amrex::Real(ii) + fcy(i+ii,domhi_y+1,k+kk,0);
                        offset_z = amrex::Real(kk) + fcy(i+ii,domhi_y+1,k+kk,1);
                    }
                    offset_y = amrex::Real(domhi_y-j) + amrex::Real(0.5);
                }

                // Do corrections for entire z-face
                if (klo_test)
                {
                    if (!ilo_test && !ihi_test && !jlo_test && !jhi_test)
                    {
                        offset_x = amrex::Real(ii) + fcz(i+ii,j+jj,domlo_z,0);
                        offset_y = amrex::Real(jj) + fcz(i+ii,j+jj,domlo_z,1);
                    }
                    offset_z = amrex::Real(domlo_z-k) - amrex::Real(0.5);

                } else if (khi_test) {
                    if (!ilo_test && !ihi_test && !jlo_test && !jhi_test)
                    {
                        offset_x = amrex::Real(ii) + fcz(i+ii,j+jj,domhi_z+1,0);
                        offset_y = amrex::Real(jj) + fcz(i+ii,j+jj,domhi_z+1,1);
                    }
                    offset_z = amrex::Real(domhi_z-k) + amrex::Real(0.5);
                }

                offset_x -= ccent(i,j,k,0);
                offset_y -= ccent(i,j,k,1);
                offset_z -= ccent(i,j,k,2);

                // Do not read state for an excluded zero row: zero times NaN is still NaN.
                if (offset_x != amrex::Real(0) || offset_y != amrex::Real(0) ||
                    offset_z != amrex::Real(0)) {
                    amrex::Real du = state(i+ii,j+jj,k+kk,n) - state(i,j,k,n);
                    normal_eq.add(offset_x,offset_y,offset_z,du);
                }
            }
        }}}

    const auto slopes = normal_eq.solve();
    auto xslope = slopes[0];
    auto yslope = slopes[1];
    auto zslope = slopes[2];

        // This will over-write the values of xslope and yslope if appropriate
    amrex_overwrite_with_regular_slopes_extdir(i,j,k,n,xslope,yslope,zslope,state,vfrac,
                                               edlo_x,edlo_y,edlo_z,edhi_x,edhi_y,edhi_z,
                                               domlo_x,domlo_y,domlo_z,domhi_x,domhi_y,domhi_z,
                                               max_order);

    // Zero out slopes outside of an extdir (or hoextrap) boundary
    // TODO:  is this the right thing to do at a HOEXTRAP boundary??
    if ( (edlo_x && i < domlo_x) || (edhi_x && i > domhi_x) ||
         (edlo_y && j < domlo_y) || (edhi_y && j > domhi_y) ||
         (edlo_z && k < domlo_z) || (edhi_z && k > domhi_z) )
    {
          xslope = 0.; yslope = 0.; zslope = 0.;
    }

    return {xslope,yslope,zslope};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real
amrex_calc_alpha_stencil(amrex::Real q_hat, amrex::Real q_max,
                         amrex::Real q_min, amrex::Real state, amrex::Real alpha) noexcept
{
    using namespace amrex::literals;

    auto alpha_temp = amrex::Real(0.0);
#ifdef AMREX_USE_FLOAT
    auto sml  = amrex::Real(1.0e-5);
#else
    auto sml  = amrex::Real(1.0e-13);
#endif

    if ((q_hat-state) > sml) {
        alpha_temp = amrex::min(1.0_rt,(q_max-state)/(q_hat-state));
    } else if ((q_hat-state) < -sml) {
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
                         amrex::Array4<amrex::Real const> const& fcz,
                         amrex::Array4<amrex::Real const> const& ccent) noexcept
{
    auto alpha = amrex::Real(2.0);

    int cuts_x = 0; int cuts_y = 0; int cuts_z = 0;

    // Compute how many cut or regular faces do we have in 3x3 block
    for(int kk(-1); kk<=1; kk++)
    {
        for(int jj(-1); jj<=1; jj++){
          for(int ii(-1); ii<=1; ii++){
            if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0))
            {
                if ((ii==-1 || ii==1) && jj==0 && kk==0) { cuts_x++; }
                if ((jj==-1 || jj==1) && ii==0 && kk==0) { cuts_y++; }
                if ((kk==-1 || kk==1) && ii==0 && jj==0) { cuts_z++; }
            }
          }
        }
    }

    amrex::Real xc = ccent(i,j,k,0); // centroid of cell (i,j,k)
    amrex::Real yc = ccent(i,j,k,1);
    amrex::Real zc = ccent(i,j,k,2);


    //Reconstruct values at the face centroids and compute the limiter
    if (flag(i,j,k).isConnected(0,1,0)) {
        amrex::Real xf = fcy(i,j+1,k,0); // local (x,z) of centroid of y-face we are extrapolating to
        amrex::Real zf = fcy(i,j+1,k,1); // local (x,z) of centroid of y-face we are extrapolating to

        amrex::Real delta_x = xf  - xc;
        amrex::Real delta_y = amrex::Real(0.5) - yc;
        amrex::Real delta_z = zf  - zc;

        amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                           + delta_y * slopes[1]
                                           + delta_z * slopes[2];

        // Compute max and min values in a 3x2x3 stencil
        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);
        for(int kk(-1); kk<=1; kk++) {
            for(int jj(0); jj<=1; jj++){
                for(int ii(-1); ii<=1; ii++){
                    if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0)) {
                        q_max = std::max(state(i+ii,j+jj,k+kk,n), q_max);
                        q_min = std::min(state(i+ii,j+jj,k+kk,n), q_min);
                    }
                }
            }
        }

        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if (flag(i,j,k).isConnected(0,-1,0)){
        amrex::Real xf = fcy(i,j,k,0); // local (x,z) of centroid of y-face we are extrapolating to
        amrex::Real zf = fcy(i,j,k,1); // local (x,z) of centroid of y-face we are extrapolating to

        amrex::Real delta_x = xf  - xc;
        amrex::Real delta_y = amrex::Real(0.5) + yc;
        amrex::Real delta_z = zf  - zc;

        amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                           - delta_y * slopes[1]
                                           + delta_z * slopes[2];

        // Compute max and min values in a 3x2x3 stencil
        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);
        for(int kk(-1); kk<=1; kk++) {
            for(int jj(-1); jj<=0; jj++){
                for(int ii(-1); ii<=1; ii++){
                    if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0)) {
                        q_max = std::max(state(i+ii,j+jj,k+kk,n), q_max);
                        q_min = std::min(state(i+ii,j+jj,k+kk,n), q_min);
                    }
                }
            }
        }

        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if (flag(i,j,k).isConnected(1,0,0)) {
        amrex::Real yf = fcx(i+1,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
        amrex::Real zf = fcx(i+1,j,k,1); // local (y,z) of centroid of x-face we are extrapolating to

        amrex::Real delta_x = amrex::Real(0.5) - xc;
        amrex::Real delta_y = yf  - yc;
        amrex::Real delta_z = zf  - zc;

        amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                           + delta_y * slopes[1]
                                           + delta_z * slopes[2];

        // Compute max and min values in a 2x3x3 stencil
        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);
        for(int kk(-1); kk<=1; kk++) {
            for(int jj(-1); jj<=1; jj++){
                for(int ii(0); ii<=1; ii++){
                    if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0)) {
                        q_max = std::max(state(i+ii,j+jj,k+kk,n), q_max);
                        q_min = std::min(state(i+ii,j+jj,k+kk,n), q_min);
                    }
                }
            }
        }

        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if (flag(i,j,k).isConnected(-1,0,0)) {
        amrex::Real yf = fcx(i,j,k,0); // local (y,z) of centroid of x-face we are extrapolating to
        amrex::Real zf = fcx(i,j,k,1); // local (y,z) of centroid of x-face we are extrapolating to

        amrex::Real delta_x = amrex::Real(0.5) + xc;
        amrex::Real delta_y = yf  - yc;
        amrex::Real delta_z = zf  - zc;

        amrex::Real q_hat = state(i,j,k,n) - delta_x * slopes[0]
                                           + delta_y * slopes[1]
                                           + delta_z * slopes[2];

        // Compute max and min values in a 3x2x3 stencil
        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);
        for(int kk(-1); kk<=1; kk++) {
            for(int jj(-1); jj<=1; jj++){
                for(int ii(-1); ii<=0; ii++){
                    if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0)) {
                        q_max = std::max(state(i+ii,j+jj,k+kk,n), q_max);
                        q_min = std::min(state(i+ii,j+jj,k+kk,n), q_min);
                    }
                }
            }
        }
        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if (flag(i,j,k).isConnected(0,0,1)) {
        amrex::Real xf = fcz(i,j,k+1,0); // local (x,y) of centroid of z-face we are extrapolating to
        amrex::Real yf = fcz(i,j,k+1,1); // local (x,y) of centroid of z-face we are extrapolating to

        amrex::Real delta_x = xf  - xc;
        amrex::Real delta_y = yf  - yc;
        amrex::Real delta_z = amrex::Real(0.5) - zc;

        amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                           + delta_y * slopes[1]
                                           + delta_z * slopes[2];

        // Compute max and min values in a 3x3x2 stencil
        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);
        for(int kk(0); kk<=1; kk++) {
            for(int jj(-1); jj<=1; jj++){
                for(int ii(-1); ii<=1; ii++){
                    if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0)) {
                        q_max = std::max(state(i+ii,j+jj,k+kk,n), q_max);
                        q_min = std::min(state(i+ii,j+jj,k+kk,n), q_min);
                    }
                }
            }
        }
        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }
    if (flag(i,j,k).isConnected(0,0,-1)) {
        amrex::Real xf = fcz(i,j,k,0); // local (x,y) of centroid of z-face we are extrapolating to
        amrex::Real yf = fcz(i,j,k,1); // local (x,y) of centroid of z-face we are extrapolating to

        amrex::Real delta_x = xf  - xc;
        amrex::Real delta_y = yf  - yc;
        amrex::Real delta_z = amrex::Real(0.5) + zc;

        amrex::Real q_hat = state(i,j,k,n) + delta_x * slopes[0]
                                           + delta_y * slopes[1]
                                           - delta_z * slopes[2];

        // Compute max and min values in a 3x2x3 stencil
        amrex::Real q_min = state(i,j,k,n);
        amrex::Real q_max = state(i,j,k,n);
        for(int kk(-1); kk<=0; kk++) {
            for(int jj(-1); jj<=1; jj++){
                for(int ii(-1); ii<=1; ii++){
                    if (flag(i,j,k).isConnected(ii,jj,kk) && !(ii==0 && jj==0 && kk==0)) {
                        q_max = std::max(state(i+ii,j+jj,k+kk,n), q_max);
                        q_min = std::min(state(i+ii,j+jj,k+kk,n), q_min);
                    }
                }
            }
        }
        alpha = amrex_calc_alpha_stencil(q_hat, q_max, q_min, state(i,j,k,n), alpha);
    }

    amrex::Real xalpha = alpha;
    amrex::Real yalpha = alpha;
    amrex::Real zalpha = alpha;

    //Zeroing out the slopes in the direction where a covered face exists.
    if (cuts_x<2) { xalpha = 0; }
    if (cuts_y<2) { yalpha = 0; }
    if (cuts_z<2) { zalpha = 0; }

    return {xalpha,yalpha,zalpha};
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
                     amrex::Array4<amrex::Real const> const& fcz,
                     amrex::Array4<amrex::EBCellFlag const> const& flag,
                     int max_order) noexcept
{
    using namespace amrex::literals;

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> alpha_lim;

    slopes = amrex_calc_slopes_eb(i,j,k,n,state,ccent,vfrac,flag,max_order);

    alpha_lim = amrex_calc_alpha_limiter(i,j,k,n,state,flag,slopes,fcx,fcy,fcz,ccent);

    // Setting limiter to 1 for stencils that just consists of non-EB cells because
    // amrex_calc_slopes_eb routine will call the slope routine for non-EB stencils that has already a limiter
    if ( max_order > 0 && vfrac(i,j,k) == 1.0_rt && vfrac(i-1,j,k) == 1.0_rt && vfrac(i+1,j,k) == 1.0_rt) {
        alpha_lim[0] = 1.0_rt;
    }

    if ( max_order > 0 && vfrac(i,j,k) == 1.0_rt && vfrac(i,j-1,k) == 1.0_rt && vfrac(i,j+1,k) == 1.0_rt) {
        alpha_lim[1] = 1.0_rt;
    }

    if ( max_order > 0 && vfrac(i,j,k) == 1.0_rt && vfrac(i,j,k-1) == 1.0_rt && vfrac(i,j,k+1) == 1.0_rt) {
        alpha_lim[2] = 1.0_rt;
    }

    return {alpha_lim[0]*slopes[0],alpha_lim[1]*slopes[1],alpha_lim[2]*slopes[2]};
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
                            amrex::Array4<amrex::Real const> const& fcz,
                            amrex::Array4<amrex::EBCellFlag const> const& flag,
                            bool edlo_x, bool edlo_y, bool edlo_z,
                            bool edhi_x, bool edhi_y, bool edhi_z,
                            int domlo_x, int domlo_y, int domlo_z,
                            int domhi_x, int domhi_y, int domhi_z,
                            int max_order) noexcept
{
    using namespace amrex::literals;

    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> slopes;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> alpha_lim;

    slopes = amrex_calc_slopes_extdir_eb(i,j,k,n,state,ccent,vfrac,fcx,fcy,fcz,flag,
                                         edlo_x,edlo_y,edlo_z,edhi_x,edhi_y,edhi_z,
                                         domlo_x,domlo_y,domlo_z,domhi_x,domhi_y,domhi_z,max_order);
    alpha_lim = amrex_calc_alpha_limiter(i,j,k,n,state,flag,slopes,fcx,fcy,fcz,ccent);

    // Setting limiter to 1 for stencils that just consists of non-EB cells because
    // amrex_calc_slopes_extdir_eb routine will call the slope routine for non-EB stencils that has already a limiter
    if ( max_order > 0 && vfrac(i,j,k) == 1.0_rt && vfrac(i-1,j,k) == 1.0_rt && vfrac(i+1,j,k) == 1.0_rt) {
        alpha_lim[0] = 1.0_rt;
    }

    if ( max_order > 0 && vfrac(i,j,k) == 1.0_rt && vfrac(i,j-1,k) == 1.0_rt && vfrac(i,j+1,k) == 1.0_rt) {
        alpha_lim[1] = 1.0_rt;
    }

    if ( max_order > 0 && vfrac(i,j,k) == 1.0_rt && vfrac(i,j,k-1) == 1.0_rt && vfrac(i,j,k+1) == 1.0_rt) {
        alpha_lim[2] = 1.0_rt;
    }

    return {alpha_lim[0]*slopes[0],alpha_lim[1]*slopes[1],alpha_lim[2]*slopes[2]};
}

#endif
