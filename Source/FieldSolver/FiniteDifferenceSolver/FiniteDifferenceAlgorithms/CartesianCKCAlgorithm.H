/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_CKC_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_CKC_H_

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>

#include <algorithm>
#include <array>

/**
 * This struct contains only static functions to initialize the stencil coefficients
 * and to compute finite-difference derivatives for the Cartesian CKC algorithm.
 */
struct CartesianCKCAlgorithm {

    static void InitializeStencilCoefficients (
        std::array<amrex::Real,3>& cell_size,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_x,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_y,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_z ) {

        using namespace amrex;

        // Compute Cole-Karkkainen-Cowan coefficients according
        // to Cowan - PRST-AB 16, 041303 (2013)
        Real const inv_dx = 1._rt/cell_size[0];
        Real const inv_dy = 1._rt/cell_size[1];
        Real const inv_dz = 1._rt/cell_size[2];
#if defined WARPX_DIM_3D
        Real const delta = std::max( { inv_dx,inv_dy,inv_dz } );
        Real const rx = (inv_dx/delta)*(inv_dx/delta);
        Real const ry = (inv_dy/delta)*(inv_dy/delta);
        Real const rz = (inv_dz/delta)*(inv_dz/delta);
        Real const beta = 0.125_rt*(1._rt - rx*ry*rz/(ry*rz + rz*rx + rx*ry));
        Real const betaxy = ry*beta*inv_dx;
        Real const betaxz = rz*beta*inv_dx;
        Real const betayx = rx*beta*inv_dy;
        Real const betayz = rz*beta*inv_dy;
        Real const betazx = rx*beta*inv_dz;
        Real const betazy = ry*beta*inv_dz;
        Real const inv_r_fac = (1._rt/(ry*rz + rz*rx + rx*ry));
        Real const gammax = ry*rz*(0.0625_rt - 0.125_rt*ry*rz*inv_r_fac);
        Real const gammay = rx*rz*(0.0625_rt - 0.125_rt*rx*rz*inv_r_fac);
        Real const gammaz = rx*ry*(0.0625_rt - 0.125_rt*rx*ry*inv_r_fac);
        Real const alphax = (1._rt - 2._rt*ry*beta - 2._rt*rz*beta - 4._rt*gammax)*inv_dx;
        Real const alphay = (1._rt - 2._rt*rx*beta - 2._rt*rz*beta - 4._rt*gammay)*inv_dy;
        Real const alphaz = (1._rt - 2._rt*rx*beta - 2._rt*ry*beta - 4._rt*gammaz)*inv_dz;
#elif defined WARPX_DIM_XZ
        Real const delta = std::max(inv_dx,inv_dz);
        Real const rx = (inv_dx/delta)*(inv_dx/delta);
        Real const rz = (inv_dz/delta)*(inv_dz/delta);
        constexpr Real beta = 0.125_rt;
        Real const betaxz = beta*rz*inv_dx;
        Real const betazx = beta*rx*inv_dz;
        Real const alphax = (1._rt - 2._rt*rz*beta)*inv_dx;
        Real const alphaz = (1._rt - 2._rt*rx*beta)*inv_dz;
        // Other coefficients are 0 in 2D Cartesian
        // (and will actually not be used in the stencil)
        constexpr Real gammax=0._rt, gammay=0._rt, gammaz=0._rt;
        constexpr Real betaxy=0._rt, betazy=0._rt, betayx=0._rt, betayz=0._rt;
        constexpr Real alphay=0._rt;
#endif

        // Store the coefficients in array `stencil_coefs`, in prescribed order
        stencil_coefs_x.resize(6);
        stencil_coefs_x[0] = inv_dx;
        stencil_coefs_x[1] = alphax;
        stencil_coefs_x[2] = betaxy;
        stencil_coefs_x[3] = betaxz;
        stencil_coefs_x[4] = gammax*inv_dx;
        stencil_coefs_y.resize(6);
        stencil_coefs_y[0] = inv_dy;
        stencil_coefs_y[1] = alphay;
        stencil_coefs_y[2] = betayz;
        stencil_coefs_y[3] = betayx;
        stencil_coefs_y[4] = gammay*inv_dy;
        stencil_coefs_z.resize(6);
        stencil_coefs_z[0] = inv_dz;
        stencil_coefs_z[1] = alphaz;
        stencil_coefs_z[2] = betazx;
        stencil_coefs_z[3] = betazy;
        stencil_coefs_z[4] = gammaz*inv_dz;
    }

    /**
     * Perform derivative along x on a cell-centered grid, from a nodal field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDx (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_x, int const n_coefs_x,
        int const i, int const j, int const k ) {

        amrex::Real const alphax = coefs_x[1];
        amrex::Real const betaxy = coefs_x[2];
        amrex::Real const betaxz = coefs_x[3];
        amrex::Real const gammax = coefs_x[4];
#if defined WARPX_DIM_3D
        return alphax * (F(i+1,j  ,k  ) - F(i,  j,  k  ))
             + betaxy * (F(i+1,j+1,k  ) - F(i  ,j+1,k  )
                      +  F(i+1,j-1,k  ) - F(i  ,j-1,k  ))
             + betaxz * (F(i+1,j  ,k+1) - F(i  ,j  ,k+1)
                      +  F(i+1,j  ,k-1) - F(i  ,j  ,k-1))
             + gammax * (F(i+1,j+1,k+1) - F(i  ,j+1,k+1)
                      +  F(i+1,j-1,k+1) - F(i  ,j-1,k+1)
                      +  F(i+1,j+1,k-1) - F(i  ,j+1,k-1)
                      +  F(i+1,j-1,k-1) - F(i  ,j-1,k-1));
#elif (defined WARPX_DIM_XZ)
        return alphax * (F(i+1,j  ,k  ) - F(i,  j,  k  ))
             + betaxz * (F(i+1,j+1,k  ) - F(i  ,j+1,k  )
                      +  F(i+1,j-1,k  ) - F(i  ,j-1,k  ));
#endif
    };

    /**
     * Perform derivative along x on a nodal grid, from a cell-centered field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDx (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_x, int const n_coefs_x,
        int const i, int const j, int const k ) {

        amrex::Real const inv_dx = coefs_x[0];
        return inv_dx*( F(i,j,k) - F(i-1,j,k) );
    };

    /**
     * Perform derivative along y on a cell-centered grid, from a nodal field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDy (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k ) {

        using namespace amrex;
#if defined WARPX_DIM_3D
        Real const alphay = coefs_y[1];
        Real const betayz = coefs_y[2];
        Real const betayx = coefs_y[3];
        Real const gammay = coefs_y[4];
        return alphay * (F(i  ,j+1,k  ) - F(i  ,j  ,k  ))
             + betayx * (F(i+1,j+1,k  ) - F(i+1,j  ,k  )
                      +  F(i-1,j+1,k  ) - F(i-1,j  ,k  ))
             + betayz * (F(i  ,j+1,k+1) - F(i  ,j  ,k+1)
                      +  F(i  ,j+1,k-1) - F(i  ,j  ,k-1))
             + gammay * (F(i+1,j+1,k+1) - F(i+1,j  ,k+1)
                      +  F(i-1,j+1,k+1) - F(i-1,j  ,k+1)
                      +  F(i+1,j+1,k-1) - F(i+1,j  ,k-1)
                       +  F(i-1,j+1,k-1) - F(i-1,j  ,k-1));
#elif (defined WARPX_DIM_XZ)
            return 0._rt; // 2D Cartesian: derivative along y is 0
#endif
    };

    /**
     * Perform derivative along y on a nodal grid, from a cell-centered field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDy (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k ) {

        using namespace amrex;
#if defined WARPX_DIM_3D
        Real const inv_dy = coefs_y[0];
        return inv_dy*( F(i,j,k) - F(i,j-1,k) );
#elif (defined WARPX_DIM_XZ)
        return 0._rt; // 2D Cartesian: derivative along y is 0
#endif
    };

    /**
     * Perform derivative along z on a cell-centered grid, from a nodal field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k ) {

        using namespace amrex;
        Real const alphaz = coefs_z[1];
        Real const betazx = coefs_z[2];
        Real const betazy = coefs_z[3];
        Real const gammaz = coefs_z[4];
#if defined WARPX_DIM_3D
        return alphaz * (F(i  ,j  ,k+1) - F(i  ,j  ,k  ))
             + betazx * (F(i+1,j  ,k+1) - F(i+1,j  ,k  )
                      +  F(i-1,j  ,k+1) - F(i-1,j  ,k  ))
             + betazy * (F(i  ,j+1,k+1) - F(i  ,j+1,k  )
                      +  F(i  ,j-1,k+1) - F(i  ,j-1,k  ))
             + gammaz * (F(i+1,j+1,k+1) - F(i+1,j+1,k  )
                      +  F(i-1,j+1,k+1) - F(i-1,j+1,k  )
                      +  F(i+1,j-1,k+1) - F(i+1,j-1,k  )
                      +  F(i-1,j-1,k+1) - F(i-1,j-1,k  ));
#elif (defined WARPX_DIM_XZ)
        return alphaz * (F(i  ,j+1,k  ) - F(i  ,j  ,k  ))
             + betazx * (F(i+1,j+1,k  ) - F(i+1,j  ,k  )
                      +  F(i-1,j+1,k  ) - F(i-1,j  ,k  ));
#endif
    };

    /**
     * Perform derivative along z on a nodal grid, from a cell-centered field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k ) {

        amrex::Real const inv_dz = coefs_z[0];
#if defined WARPX_DIM_3D
        return inv_dz*( F(i,j,k) - F(i,j,k-1) );
#elif (defined WARPX_DIM_XZ)
        return inv_dz*( F(i,j,k) - F(i,j-1,k) );
#endif
    };

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_CKC_H_
