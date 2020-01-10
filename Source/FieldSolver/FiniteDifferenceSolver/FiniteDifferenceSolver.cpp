

void FiniteDifferenceSolver::FiniteDifferenceSolver ( std::array<Real,3> cell_size ) {
    // Select algorithm (The choice of algorithm is a runtime option,
    // but we compile code for each algorithm, using templates)
    if (fdtd_algo == MaxwellSolverAlgo::Yee){
        YeeAlgorithm::InitializeStencilCoefficients( cell_size,
            stencil_coefs_x, stencil_coefs_y, stencil_coefs_z );
    } else if (fdtd_algo == MaxwellSolverAlgo::CKC) {
        CKCAlgorithm::InitializeStencilCoefficients( cell_size,
            stencil_coefs_x, stencil_coefs_y, stencil_coefs_z );
    } else {
        amrex::Abort("Unknown algorithm");
    }
};
