#include "FiniteDifferenceSolver.H"
#include "WarpX.H"

// Constructor
FiniteDifferenceSolver::FiniteDifferenceSolver ( int const fdtd_algo,
    std::array<amrex::Real,3> cell_size ) {

    // Register the type of finite-difference algorithm
    m_fdtd_algo = fdtd_algo;

    // Calculate coefficients of finite-difference stencil
#ifdef WARPX_DIM_RZ
    m_dr = cell_size[0];
    m_nmodes = WarpX::GetInstance().n_rz_azimuthal_modes;
    m_rmin = WarpX::GetInstance().Geom(0).ProbLo(0);
    if (fdtd_algo == MaxwellSolverAlgo::Yee){
        CylindricalYeeAlgorithm::InitializeStencilCoefficients( cell_size,
            stencil_coefs_r, stencil_coefs_z );
#else
    if (fdtd_algo == MaxwellSolverAlgo::Yee){
        YeeAlgorithm::InitializeStencilCoefficients( cell_size,
            stencil_coefs_x, stencil_coefs_y, stencil_coefs_z );
    } else if (fdtd_algo == MaxwellSolverAlgo::CKC) {
        CKCAlgorithm::InitializeStencilCoefficients( cell_size,
            stencil_coefs_x, stencil_coefs_y, stencil_coefs_z );
#endif
    } else {
        amrex::Abort("Unknown algorithm");
    }
};
