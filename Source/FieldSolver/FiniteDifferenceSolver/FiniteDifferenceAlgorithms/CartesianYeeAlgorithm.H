/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_YEE_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_YEE_H_

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>

#include <array>

/**
 * This struct contains only static functions to initialize the stencil coefficients
 * and to compute finite-difference derivatives for the Cartesian Yee algorithm.
 */
struct CartesianYeeAlgorithm {

    static void InitializeStencilCoefficients (
        std::array<amrex::Real,3>& cell_size,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_x,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_y,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_z ) {

        using namespace amrex;
        // Store the inverse cell size along each direction in the coefficients
        stencil_coefs_x.resize(1);
        stencil_coefs_x[0] = 1._rt/cell_size[0];
        stencil_coefs_y.resize(1);
        stencil_coefs_y[0] = 1._rt/cell_size[1];
        stencil_coefs_z.resize(1);
        stencil_coefs_z[0] = 1._rt/cell_size[2];
    }

    /**
     * Perform derivative along x on a cell-centered grid, from a nodal field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDx (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_x, int const n_coefs_x,
        int const i, int const j, int const k ) {

        amrex::Real const inv_dx = coefs_x[0];
        return inv_dx*( F(i+1,j,k) - F(i,j,k) );
    };

    /**
     * Perform derivative along x on a nodal grid, from a cell-centered field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDx (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_x, int const n_coefs_x,
        int const i, int const j, int const k ) {

        amrex::Real const inv_dx = coefs_x[0];
        return inv_dx*( F(i,j,k) - F(i-1,j,k) );
    };

    /**
     * Perform derivative along y on a cell-centered grid, from a nodal field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDy (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k ) {

        using namespace amrex;
#if defined WARPX_DIM_3D
        Real const inv_dy = coefs_y[0];
        return inv_dy*( F(i,j+1,k) - F(i,j,k) );
#elif (defined WARPX_DIM_XZ)
        return 0._rt; // 2D Cartesian: derivative along y is 0
#endif
    };

    /**
     * Perform derivative along y on a nodal grid, from a cell-centered field `F`*/
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
     * Perform derivative along z on a cell-centered grid, from a nodal field `F`*/
   AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k ) {

        using namespace amrex;
        Real const inv_dz = coefs_z[0];
#if defined WARPX_DIM_3D
        return inv_dz*( F(i,j,k+1) - F(i,j,k) );
#elif (defined WARPX_DIM_XZ)
        return inv_dz*( F(i,j+1,k) - F(i,j,k) );
#endif
    };

    /**
     * Perform derivative along z on a nodal grid, from a cell-centered field `F`*/
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k ) {

        using namespace amrex;
        Real const inv_dz = coefs_z[0];
#if defined WARPX_DIM_3D
        return inv_dz*( F(i,j,k) - F(i,j,k-1) );
#elif (defined WARPX_DIM_XZ)
        return inv_dz*( F(i,j,k) - F(i,j-1,k) );
#endif
    };

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_YEE_H_
