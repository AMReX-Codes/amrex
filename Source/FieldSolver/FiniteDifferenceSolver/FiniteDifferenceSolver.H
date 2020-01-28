#ifndef WARPX_FINITE_DIFFERENCE_SOLVER_H_
#define WARPX_FINITE_DIFFERENCE_SOLVER_H_

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
        FiniteDifferenceSolver (
            int const fdtd_algo,
            std::array<amrex::Real,3> cell_size,
            bool const do_nodal );

        void EvolveB ( VectorField& Bfield,
                       VectorField const& Efield,
                       amrex::Real const dt );
    private:

        int m_fdtd_algo;
        bool m_do_nodal;

#ifdef WARPX_DIM_RZ
        amrex::Real m_dr, m_rmin;
        amrex::Real m_nmodes;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_r;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_z;
#else
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_x;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_y;
        amrex::Gpu::ManagedVector<amrex::Real> stencil_coefs_z;
#endif

#ifdef WARPX_DIM_RZ
        template< typename T_Algo >
        void EvolveBCylindrical ( VectorField& Bfield,
                               VectorField const& Efield,
                               amrex::Real const dt );
#else
        template< typename T_Algo >
        void EvolveBCartesian ( VectorField& Bfield,
                               VectorField const& Efield,
                               amrex::Real const dt );
#endif

};

#endif // WARPX_FINITE_DIFFERENCE_SOLVER_H_
