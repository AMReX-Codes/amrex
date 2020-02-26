/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>

#include <array>

/**
 * This struct contains only static functions to initialize the stencil coefficients
 * and to compute finite-difference derivatives for the Cartesian Yee algorithm.
 */
struct CylindricalYeeAlgorithm {

    static void InitializeStencilCoefficients (
        std::array<amrex::Real,3>& cell_size,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_r,
        amrex::Gpu::ManagedVector<amrex::Real>& stencil_coefs_z ) {

        using namespace amrex;
        // Store the inverse cell size along each direction in the coefficients
        stencil_coefs_r.resize(1);
        stencil_coefs_r[0] = 1._rt/cell_size[0];  // 1./dr
        stencil_coefs_z.resize(1);
        stencil_coefs_z[0] = 1._rt/cell_size[2];  // 1./dz
    }

    /** Applies the differential operator `1/r * d(rF)/dr`,
     * where `F` is on a *nodal* grid in `r`
     * and the differential operator is evaluated on a *cell-centered* grid.
     * The input parameter `r` is given at the cell-centered position */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDrr_over_r (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const r, amrex::Real const dr,
        amrex::Real const * const coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        using namespace amrex;
        Real const inv_dr = coefs_r[0];
        return 1._rt/r * inv_dr*( (r+0.5_rt*dr)*F(i+1,j,k,comp) - (r-0.5_rt*dr)*F(i,j,k,comp) );
    };

    /** Applies the differential operator `1/r * d(rF)/dr`,
     * where `F` is on a *cell-centered* grid in `r`
     * and the differential operator is evaluated on a *nodal* grid.
     * The input parameter `r` is given at the cell-centered position */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDrr_over_r (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const r, amrex::Real const dr,
        amrex::Real const * const coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        using namespace amrex;
        Real const inv_dr = coefs_r[0];
        return 1._rt/r * inv_dr*( (r+0.5_rt*dr)*F(i,j,k,comp) - (r-0.5_rt*dr)*F(i-1,j,k,comp) );
    };

    /**
     * Perform derivative along r on a cell-centered grid, from a nodal field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDr (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        using namespace amrex;
        Real const inv_dr = coefs_r[0];
        return inv_dr*( F(i+1,j,k,comp) - F(i,j,k,comp) );
    };

    /**
     * Perform derivative along r on a nodal grid, from a cell-centered field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDr (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_r, int const n_coefs_r,
        int const i, int const j, int const k, int const comp ) {

        using namespace amrex;
        Real const inv_dr = coefs_r[0];
        return inv_dr*( F(i,j,k,comp) - F(i-1,j,k,comp) );
    };

    /**
     * Perform derivative along z on a cell-centered grid, from a nodal field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real UpwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k, int const comp ) {

        amrex::Real const inv_dz = coefs_z[0];
        return inv_dz*( F(i,j+1,k,comp) - F(i,j,k,comp) );
    };

    /**
     * Perform derivative along z on a nodal grid, from a cell-centered field `F` */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real DownwardDz (
        amrex::Array4<amrex::Real> const& F,
        amrex::Real const * const coefs_z, int const n_coefs_z,
        int const i, int const j, int const k, int const comp ) {

        amrex::Real const inv_dz = coefs_z[0];
        return inv_dz*( F(i,j,k,comp) - F(i,j-1,k,comp) );
    };

};

#endif // WARPX_FINITE_DIFFERENCE_ALGORITHM_CYLINDRICAL_YEE_H_
