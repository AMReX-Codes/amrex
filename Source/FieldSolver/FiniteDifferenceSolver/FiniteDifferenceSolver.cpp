/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Utils/WarpXAlgorithmSelection.H"
#ifdef WARPX_DIM_RZ
#    include "FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#else
#    include "FiniteDifferenceAlgorithms/CartesianYeeAlgorithm.H"
#    include "FiniteDifferenceAlgorithms/CartesianCKCAlgorithm.H"
#    include "FiniteDifferenceAlgorithms/CartesianNodalAlgorithm.H"
#endif
#include "FiniteDifferenceSolver.H"
#include "WarpX.H"

/* This function initializes the stencil coefficients for the chosen finite-difference algorithm */
FiniteDifferenceSolver::FiniteDifferenceSolver (
    int const fdtd_algo,
    std::array<amrex::Real,3> cell_size,
    bool do_nodal ) {

    // Register the type of finite-difference algorithm
    m_fdtd_algo = fdtd_algo;
    m_do_nodal = do_nodal;

    // Calculate coefficients of finite-difference stencil
#ifdef WARPX_DIM_RZ
    m_dr = cell_size[0];
    m_nmodes = WarpX::GetInstance().n_rz_azimuthal_modes;
    m_rmin = WarpX::GetInstance().Geom(0).ProbLo(0);
    if (fdtd_algo == MaxwellSolverAlgo::Yee) {

        CylindricalYeeAlgorithm::InitializeStencilCoefficients( cell_size,
            m_stencil_coefs_r, m_stencil_coefs_z );
#else
    if (do_nodal) {

        CartesianNodalAlgorithm::InitializeStencilCoefficients( cell_size,
            m_stencil_coefs_x, m_stencil_coefs_y, m_stencil_coefs_z );

    } else if (fdtd_algo == MaxwellSolverAlgo::Yee) {

        CartesianYeeAlgorithm::InitializeStencilCoefficients( cell_size,
            m_stencil_coefs_x, m_stencil_coefs_y, m_stencil_coefs_z );

    } else if (fdtd_algo == MaxwellSolverAlgo::CKC) {

        CartesianCKCAlgorithm::InitializeStencilCoefficients( cell_size,
            m_stencil_coefs_x, m_stencil_coefs_y, m_stencil_coefs_z );

#endif
    } else {
        amrex::Abort("Unknown algorithm");
    }
};
