/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_NODAL_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_NODAL_H_

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>

#include <array>

/**
 * This struct contains only static functions to initialize the stencil coefficients
 * and to compute finite-difference derivatives for the Cartesian nodal algorithm.
 */
struct CartesianNodalAlgorithm {

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
     * Perform derivative along x
     * (For a solver on a staggered grid, `UpwardDx` and `DownwardDx` take into
     * account the staggering; but for `CartesianNodalAlgorithm`, they are equivalent) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDx (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_x, int const n_coefs_x,
        int const i, int const j, int const k ) {

        using namespace amrex;
        Real const inv_dx = coefs_x[0];
        return 0.5_rt*inv_dx*( F(i+1,j,k) - F(i-1,j,k) );
    };

    /**
     * Perform derivative along x
     * (For a solver on a staggered grid, `UpwardDx` and `DownwardDx` take into
     * account the staggering; but for `CartesianNodalAlgorithm`, they are equivalent) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDx (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_x, int const n_coefs_x,
        int const i, int const j, int const k ) {

        return UpwardDx( F, coefs_x, n_coefs_x, i, j, k );
        // For CartesianNodalAlgorithm, UpwardDx and DownwardDx are equivalent
    };

    /**
     * Perform derivative along y
     * (For a solver on a staggered grid, `UpwardDy` and `DownwardDy` take into
     * account the staggering; but for `CartesianNodalAlgorithm`, they are equivalent) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDy (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k ) {

        using namespace amrex;
#if defined WARPX_DIM_3D
        Real const inv_dy = coefs_y[0];
        return 0.5_rt*inv_dy*( F(i,j+1,k) - F(i,j-1,k) );
#elif (defined WARPX_DIM_XZ)
        return 0._rt; // 2D Cartesian: derivative along y is 0
#endif
    };

    /**
     * Perform derivative along y
     * (For a solver on a staggered grid, `UpwardDy` and `DownwardDy` take into
     * account the staggering; but for `CartesianNodalAlgorithm`, they are equivalent) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDy (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_y, int const n_coefs_y,
        int const i, int const j, int const k ) {

        return UpwardDy( F, coefs_y, n_coefs_y, i, j, k );
        // For CartesianNodalAlgorithm, UpwardDy and DownwardDy are equivalent
    };

    /**
     * Perform derivative along z
     * (For a solver on a staggered grid, `UpwardDz` and `DownwardDz` take into
     * account the staggering; but for `CartesianNodalAlgorithm`, they are equivalent) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k ) {

        using namespace amrex;
        Real const inv_dz = coefs_z[0];
#if defined WARPX_DIM_3D
        return 0.5_rt*inv_dz*( F(i,j,k+1) - F(i,j,k-1) );
#elif (defined WARPX_DIM_XZ)
        return 0.5_rt*inv_dz*( F(i,j+1,k) - F(i,j-1,k) );
#endif
    };

    /**
     * Perform derivative along z
     * (For a solver on a staggered grid, `UpwardDz` and `DownwardDz` take into
     * account the staggering; but for `CartesianNodalAlgorithm`, they are equivalent) */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k ) {

        return UpwardDz( F, coefs_z, n_coefs_z, i, j, k );
        // For CartesianNodalAlgorithm, UpwardDz and DownwardDz are equivalent
    };

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CARTESIAN_NODAL_H_
