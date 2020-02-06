/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CKC_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CKC_H_

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>

struct CKCAlgorithm {

    static void InitializeStencilCoefficients (
        std::array<amrex::Real,3>& cell_size,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_x,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_y,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_z ) {

        using namespace amrex;

        // Compute Cole-Karkkainen-Cowan coefficients according
        // to Cowan - PRST-AB 16, 041303 (2013)
        Real alphax=0, alphay=0, alphaz=0;
        Real betaxy=0, betaxz=0, betayx=0, betayz=0, betazx=0, betazy=0;
        Real gammax=0, gammay=0, gammaz=0;
        Real const inv_dx = 1./cell_size[0];
        Real const inv_dy = 1./cell_size[1];
        Real const inv_dz = 1./cell_size[2];
#if defined WARPX_DIM_3D
            const Real delta = std::max( { inv_dx,inv_dy,inv_dz } );
            const Real rx = (inv_dx/delta)*(inv_dx/delta);
            const Real ry = (inv_dy/delta)*(inv_dy/delta);
            const Real rz = (inv_dz/delta)*(inv_dz/delta);
            const Real beta = 0.125*(1. - rx*ry*rz/(ry*rz + rz*rx + rx*ry));
            betaxy = ry*beta;
            betaxz = rz*beta;
            betayx = rx*beta;
            betayz = rz*beta;
            betazx = rx*beta;
            betazy = ry*beta;
            gammax = ry*rz*(0.0625 - 0.125*ry*rz/(ry*rz + rz*rx + rx*ry));
            gammay = rx*rz*(0.0625 - 0.125*rx*rz/(ry*rz + rz*rx + rx*ry));
            gammaz = rx*ry*(0.0625 - 0.125*rx*ry/(ry*rz + rz*rx + rx*ry));
            alphax = 1. - 2.*betaxy - 2.*betaxz - 4.*gammax;
            alphay = 1. - 2.*betayx - 2.*betayz - 4.*gammay;
            alphaz = 1. - 2.*betazx - 2.*betazy - 4.*gammaz;
            betaxy *= inv_dx;
            betaxz *= inv_dx;
            betayx *= inv_dy;
            betayz *= inv_dy;
            betazx *= inv_dz;
            betazy *= inv_dz;
            alphax *= inv_dx;
            alphay *= inv_dy;
            alphaz *= inv_dz;
            gammax *= inv_dx;
            gammay *= inv_dy;
            gammaz *= inv_dz;
#elif defined WARPX_DIM_XZ
            const Real delta = std::max(inv_dx,inv_dz);
            const Real rx = (inv_dx/delta)*(inv_dx/delta);
            const Real rz = (inv_dz/delta)*(inv_dz/delta);
            betaxz = 0.125*rz;
            betazx = 0.125*rx;
            alphax = 1. - 2.*betaxz;
            alphaz = 1. - 2.*betazx;
            betaxz *= inv_dx;
            betazx *= inv_dz;
            alphax *= inv_dx;
            alphaz *= inv_dz;
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
    static amrex::Real UpwardDx(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_x, int const n_coefs_x,
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
    static amrex::Real DownwardDx(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_x, int const n_coefs_x,
        int const i, int const j, int const k ) {

        amrex::Real const inv_dx = coefs_x[0];
        return inv_dx*( F(i,j,k) - F(i-1,j,k) );
    };

    /**
    /* Perform derivative along y on a cell-centered grid, from a nodal field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDy(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_y, int const n_coefs_y,
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
    static amrex::Real DownwardDy(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_y, int const n_coefs_y,
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
    static amrex::Real UpwardDz(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_z, int const n_coefs_z,
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
    static amrex::Real DownwardDz(
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const* coefs_z, int const n_coefs_z,
        int const i, int const j, int const k ) {

        amrex::Real const inv_dz = coefs_z[0];
#if defined WARPX_DIM_3D
            return inv_dz*( F(i,j,k) - F(i,j,k-1) );
#elif (defined WARPX_DIM_XZ)
            return inv_dz*( F(i,j,k) - F(i,j-1,k) );
#endif
    };

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CKC_H_
