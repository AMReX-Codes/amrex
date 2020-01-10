#ifndef WARPX_FINITE_DIFFERENCE_SOLVER_H_
#define WARPX_FINITE_DIFFERENCE_SOLVER_H_

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
        void FiniteDifferenceSolver( std::array<Real,3> dx );

        void EvolveB( VectorField Bfield,
                    ConstVectorField Efield,
                    amrex::Real dt ) const;

        };

    private:

        int fdtd_algo;

        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_x;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_y;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_z;

        template< typename fdtd_algo >
        void EvolveB

};

#endif // WARPX_FINITE_DIFFERENCE_SOLVER_H_
