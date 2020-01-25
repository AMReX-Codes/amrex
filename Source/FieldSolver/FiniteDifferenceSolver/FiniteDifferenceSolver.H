#ifndef WARPX_FINITE_DIFFERENCE_SOLVER_H_
#define WARPX_FINITE_DIFFERENCE_SOLVER_H_

#include "WarpXAlgorithmSelection.H"
#include "FiniteDifferenceAlgorithms/YeeAlgorithm.H"
#include "FiniteDifferenceAlgorithms/CKCAlgorithm.H"
#include <AMReX_MultiFab.H>

/**
 * \brief Top-level class for the electromagnetic finite-difference solver
 *
 * TODO
 */
class FiniteDifferenceSolver
{
    public:

        using VectorField = std::array< std::unique_ptr<amrex::MultiFab>, 3 >;

        // Constructor
        FiniteDifferenceSolver ( int const fdtd_algo,
            std::array<amrex::Real,3> cell_size ) {

            // Register the type of finite-difference algorithm
            m_fdtd_algo = fdtd_algo;

            // Calculate coefficients of finite-difference stencil
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

        void EvolveB ( VectorField& Bfield,
                       VectorField const& Efield,
                       amrex::Real const dt );
    private:

        int m_fdtd_algo;

        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_x;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_y;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_z;

        template< typename T_Algo >
        void EvolveBwithAlgo ( VectorField& Bfield,
                               VectorField const& Efield,
                               amrex::Real const dt );

};

#endif // WARPX_FINITE_DIFFERENCE_SOLVER_H_
