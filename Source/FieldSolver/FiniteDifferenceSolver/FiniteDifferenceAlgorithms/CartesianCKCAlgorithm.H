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
        Real const inv_dx = 1./cell_size[0];
        Real const inv_dy = 1./cell_size[1];
        Real const inv_dz = 1./cell_size[2];
#if defined WARPX_DIM_3D
            Real const delta = std::max( { inv_dx,inv_dy,inv_dz } );
            Real const rx = (inv_dx/delta)*(inv_dx/delta);
            Real const ry = (inv_dy/delta)*(inv_dy/delta);
            Real const rz = (inv_dz/delta)*(inv_dz/delta);
            Real const beta = 0.125*(1. - rx*ry*rz/(ry*rz + rz*rx + rx*ry));
            Real const betaxy = ry*beta*inv_dx;
            Real const betaxz = rz*beta*inv_dx;
            Real const betayx = rx*beta*inv_dy;
            Real const betayz = rz*beta*inv_dy;
            Real const betazx = rx*beta*inv_dz;
            Real const betazy = ry*beta*inv_dz;
            Real const gamma = (0.0625 - 0.125*ry*rz/(ry*rz + rz*rx + rx*ry));
            Real const gammax = ry*rz*gamma;
            Real const gammay = rx*rz*gamma;
            Real const gammaz = rx*ry*gamma;
            Real const alphax = (1. - 2.*ry*beta - 2.*rz*beta - 4.*ry*rz*gamma)*inv_dx;
            Real const alphay = (1. - 2.*rx*beta - 2.*rz*beta - 4.*rx*rz*gamma)*inv_dy;
            Real const alphaz = (1. - 2.*rx*beta - 2.*ry*beta - 4.*rx*ry*gamma)*inv_dz;
#elif defined WARPX_DIM_XZ
            Real const delta = std::max(inv_dx,inv_dz);
            Real const rx = (inv_dx/delta)*(inv_dx/delta);
            Real const rz = (inv_dz/delta)*(inv_dz/delta);
            Real const beta = 0.125;
            Real const betaxz = beta*rz*inv_dx;
            Real const betazx = beta*rx*inv_dz;
            Real const alphax = (1. - 2.*rz*beta)*inv_dx;
            Real const alphaz = (1. - 2.*rx*beta)*inv_dz;
            // Other coefficients are 0 in 2D Cartesian
            // (and will actually not be used in the stencil)
            Real const gammax=0, gammay=0, gammaz=0;
            Real const betaxy=0, betazy=0, betayx=0, betayz=0;
            Real const alphay=0;
#endif

        // Store the coefficients in array `stencil_coefs`, in prescribed order
        stencil_coefs_x.resize(6);
        stencil_coefs_x[0] = inv_dx;
        stencil_coefs_x[1] = alphax;
        stencil_coefs_x[2] = betaxy;
        stencil_coefs_x[3] = betaxz;
        stencil_coefs_x[4] = gammax;
        stencil_coefs_y.resize(6);
        stencil_coefs_y[0] = inv_dy;
        stencil_coefs_y[1] = alphay;
        stencil_coefs_y[2] = betayz;
        stencil_coefs_y[3] = betayx;
        stencil_coefs_y[4] = gammay;
        stencil_coefs_z.resize(6);
        stencil_coefs_z[0] = inv_dz;
        stencil_coefs_z[1] = alphaz;
        stencil_coefs_z[2] = betazx;
        stencil_coefs_z[3] = betazy;
        stencil_coefs_z[4] = gammaz;
    }

    /**
    /* Perform derivative along x on a cell-centered grid, from a nodal field `F`*/
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
    /* Perform derivative along x on a nodal grid, from a cell-centered field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDx (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_x, int const n_coefs_x,
        int const i, int const j, int const k ) {

        amrex::Real const inv_dx = coefs_x[0];
        return inv_dx*( F(i,j,k) - F(i-1,j,k) );
    };

    /**
    /* Perform derivative along y on a cell-centered grid, from a nodal field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDy (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k ) {

#if defined WARPX_DIM_3D
            amrex::Real const alphay = coefs_y[1];
            amrex::Real const betayz = coefs_y[2];
            amrex::Real const betayx = coefs_y[3];
            amrex::Real const gammay = coefs_y[4];
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
            return 0; // 2D Cartesian: derivative along y is 0
#endif
    };

    /**
    /* Perform derivative along y on a nodal grid, from a cell-centered field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDy (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k ) {

#if defined WARPX_DIM_3D
            amrex::Real const inv_dy = coefs_y[0];
            return inv_dy*( F(i,j,k) - F(i,j-1,k) );
#elif (defined WARPX_DIM_XZ)
            return 0; // 2D Cartesian: derivative along y is 0
#endif
    };

    /**
    /* Perform derivative along z on a cell-centered grid, from a nodal field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k ) {

        amrex::Real const alphaz = coefs_z[1];
        amrex::Real const betazx = coefs_z[2];
        amrex::Real const betazy = coefs_z[3];
        amrex::Real const gammaz = coefs_z[4];
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
    /* Perform derivative along z on a nodal grid, from a cell-centered field `F`*/
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
