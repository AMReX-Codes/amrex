/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_SOLVER_H_
#define WARPX_FINITE_DIFFERENCE_SOLVER_H_

#include <AMReX_MultiFab.H>

/**
 * \brief Top-level class for the electromagnetic finite-difference solver
 *
 * Stores the coefficients of the finite-difference stencils,
 * and has member functions to update fields over one time step.
 */
class FiniteDifferenceSolver
{
    public:

        // Constructor
        /** \brief Initialize the finite-difference Maxwell solver (for a given refinement level)
         *
         * This function initializes the stencil coefficients for the chosen finite-difference algorithm
         *
         * \param fdtd_algo Identifies the chosen algorithm, as defined in WarpXAlgorithmSelection.H
         * \param cell_size Cell size along each dimension, for the chosen refinement level
         * \param do_nodal  Whether the solver is applied to a nodal or staggered grid
         */
        FiniteDifferenceSolver (
            int const fdtd_algo,
            std::array<amrex::Real,3> cell_size,
            bool const do_nodal );

        void EvolveB ( std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
                       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
                       amrex::Real const dt );

        void EvolveE ( std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
                       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
                       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
                       std::unique_ptr<amrex::MultiFab> const& Ffield,
                       amrex::Real const dt );

        void EvolveF ( std::unique_ptr<amrex::MultiFab>& Ffield,
                       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
                       std::unique_ptr<amrex::MultiFab> const& rhofield,
                       int const rhocomp,
                       amrex::Real const dt );

        void ComputeDivE ( const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
                           amrex::MultiFab& divE );


    private:

        int m_fdtd_algo;
        bool m_do_nodal;

#ifdef WARPX_DIM_RZ
        amrex::Real m_dr, m_rmin;
        amrex::Real m_nmodes;
        amrex::Gpu::ManagedVector<amrex::Real> m_stencil_coefs_r;
        amrex::Gpu::ManagedVector<amrex::Real> m_stencil_coefs_z;
#else
        amrex::Gpu::ManagedVector<amrex::Real> m_stencil_coefs_x;
        amrex::Gpu::ManagedVector<amrex::Real> m_stencil_coefs_y;
        amrex::Gpu::ManagedVector<amrex::Real> m_stencil_coefs_z;
#endif

    public:
        // The member functions below contain extended __device__ lambda.
        // In order to compile with nvcc, they need to be public.

#ifdef WARPX_DIM_RZ
        template< typename T_Algo >
        void EvolveBCylindrical (
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
            amrex::Real const dt );

        template< typename T_Algo >
        void EvolveECylindrical (
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
            std::unique_ptr<amrex::MultiFab> const& Ffield,
            amrex::Real const dt );

        template< typename T_Algo >
        void EvolveFCylindrical (
            std::unique_ptr<amrex::MultiFab>& Ffield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
            std::unique_ptr<amrex::MultiFab> const& rhofield,
            int const rhocomp,
            amrex::Real const dt );

        template< typename T_Algo >
        void ComputeDivECylindrical (
            const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
            amrex::MultiFab& divE );

#else
        template< typename T_Algo >
        void EvolveBCartesian (
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
            amrex::Real const dt );

        template< typename T_Algo >
        void EvolveECartesian (
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
            std::unique_ptr<amrex::MultiFab> const& Ffield,
            amrex::Real const dt );

        template< typename T_Algo >
        void EvolveFCartesian (
            std::unique_ptr<amrex::MultiFab>& Ffield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
            std::unique_ptr<amrex::MultiFab> const& rhofield,
            int const rhocomp,
            amrex::Real const dt );

        template< typename T_Algo >
        void ComputeDivECartesian (
            const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
            amrex::MultiFab& divE );

#endif

};

#endif // WARPX_FINITE_DIFFERENCE_SOLVER_H_
