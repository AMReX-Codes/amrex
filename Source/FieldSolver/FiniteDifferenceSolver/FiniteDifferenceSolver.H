#ifndef WARPX_FINITE_DIFFERENCE_SOLVER_H_
#define WARPX_FINITE_DIFFERENCE_SOLVER_H_

#include<FiniteDifferenceAlgorithms/YeeAlgorithm.H>

/**
 * \brief Top-level class for the electromagnetic finite-difference solver
 *
 * TODO
 */
class FiniteDifferenceSolver
{
    public:

        using VectorField = std::array< std::unique_ptr<amrex::MultiFab>, 3 >;
        using ConstVectorField = std::array< std::unique_ptr<amrex::MultiFab const>, 3 >;

        // Constructor
        void FiniteDifferenceSolver::FiniteDifferenceSolver ( std::array<Real,3> cell_size ) {

            // Register the type of finite-difference algorithm
            fdtd_algo = WarpX::maxwell_fdtd_solver_id;

            // Calculate coefficients of finite-difference stencil
            if (fdtd_algo == MaxwellSolverAlgo::Yee){
                YeeAlgorithm::InitializeStencilCoefficients( cell_size,
                    stencil_coefs_x, stencil_coefs_y, stencil_coefs_z );
//            } else if (fdtd_algo == MaxwellSolverAlgo::CKC) {
//                CKCAlgorithm::InitializeStencilCoefficients( cell_size,
//                    stencil_coefs_x, stencil_coefs_y, stencil_coefs_z );
            } else {
                amrex::Abort("Unknown algorithm");
            }
        };

        void EvolveB ( VectorField Bfield,
                       ConstVectorField Efield,
                       amrex::Real dt ) const;
    private:

        int fdtd_algo;

        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_x;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_y;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_z;

        template< typename fdtd_algo >
        void EvolveBwithAlgo ( VectorField Bfield,
                               ConstVectorField Efield,
                               amrex::Real dt ) const;

};

#endif // WARPX_FINITE_DIFFERENCE_SOLVER_H_
