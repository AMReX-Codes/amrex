#include "WarpXAlgorithmSelection.H"
#include "FiniteDifferenceSolver.H"
#ifdef WARPX_DIM_RZ
    #include "FiniteDifferenceAlgorithms/CylindricalYeeAlgorithm.H"
#else
    #include "FiniteDifferenceAlgorithms/YeeAlgorithm.H"
    #include "FiniteDifferenceAlgorithms/CKCAlgorithm.H"
    #include "FiniteDifferenceAlgorithms/NodalAlgorithm.H"
#endif
#include <AMReX_Gpu.H>

using namespace amrex;

void FiniteDifferenceSolver::EvolveB ( VectorField& Bfield,
                                       VectorField const& Efield,
                                       amrex::Real const dt ) {

   // Select algorithm (The choice of algorithm is a runtime option,
   // but we compile code for each algorithm, using templates)
#ifdef WARPX_DIM_RZ
    if (m_fdtd_algo == MaxwellSolverAlgo::Yee){
        EvolveBCylindrical <CylindricalYeeAlgorithm> ( Bfield, Efield, dt );
#else
    if (m_do_nodal) {
        EvolveBCartesian <NodalAlgorithm> ( Bfield, Efield, dt );
    } else if (m_fdtd_algo == MaxwellSolverAlgo::Yee) {
        EvolveBCartesian <YeeAlgorithm> ( Bfield, Efield, dt );
    } else if (m_fdtd_algo == MaxwellSolverAlgo::CKC) {
        EvolveBCartesian <CKCAlgorithm> ( Bfield, Efield, dt );
#endif
    } else {
        amrex::Abort("Unknown algorithm");
    }

}

#ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveBCartesian ( VectorField& Bfield,
                                               VectorField const& Efield,
                                               amrex::Real const dt ) {

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bfield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        auto const& Bx = Bfield[0]->array(mfi);
        auto const& By = Bfield[1]->array(mfi);
        auto const& Bz = Bfield[2]->array(mfi);
        auto const& Ex = Efield[0]->array(mfi);
        auto const& Ey = Efield[1]->array(mfi);
        auto const& Ez = Efield[2]->array(mfi);

        // Extract stencil coefficients
        Real const* AMREX_RESTRICT coefs_x = stencil_coefs_x.dataPtr();
        int const n_coefs_x = stencil_coefs_x.size();
        Real const* AMREX_RESTRICT coefs_y = stencil_coefs_y.dataPtr();
        int const n_coefs_y = stencil_coefs_y.size();
        Real const* AMREX_RESTRICT coefs_z = stencil_coefs_z.dataPtr();
        int const n_coefs_z = stencil_coefs_z.size();

        // Extract tileboxes for which to loop
        const Box& tbx  = mfi.tilebox(Bfield[0]->ixType().ixType());
        const Box& tby  = mfi.tilebox(Bfield[1]->ixType().ixType());
        const Box& tbz  = mfi.tilebox(Bfield[2]->ixType().ixType());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tbx, tby, tbz,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bx(i, j, k) += dt * T_Algo::UpwardDz(Ey, coefs_z, n_coefs_z, i, j, k)
                             - dt * T_Algo::UpwardDy(Ez, coefs_y, n_coefs_y, i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                By(i, j, k) += dt * T_Algo::UpwardDx(Ez, coefs_x, n_coefs_x, i, j, k)
                             - dt * T_Algo::UpwardDz(Ex, coefs_z, n_coefs_z, i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bz(i, j, k) += dt * T_Algo::UpwardDy(Ex, coefs_y, n_coefs_y, i, j, k)
                             - dt * T_Algo::UpwardDx(Ey, coefs_x, n_coefs_x, i, j, k);
            }

        );

    }

}

#else // corresponds to ifndef WARPX_DIM_RZ

template<typename T_Algo>
void FiniteDifferenceSolver::EvolveBCylindrical ( VectorField& Bfield,
                                               VectorField const& Efield,
                                               amrex::Real const dt ) {

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bfield[0], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        auto const& Br = Bfield[0]->array(mfi);
        auto const& Bt = Bfield[1]->array(mfi);
        auto const& Bz = Bfield[2]->array(mfi);
        auto const& Er = Efield[0]->array(mfi);
        auto const& Et = Efield[1]->array(mfi);
        auto const& Ez = Efield[2]->array(mfi);

        // Extract stencil coefficients
        Real const* AMREX_RESTRICT coefs_r = stencil_coefs_r.dataPtr();
        int const n_coefs_r = stencil_coefs_r.size();
        Real const* AMREX_RESTRICT coefs_z = stencil_coefs_z.dataPtr();
        int const n_coefs_z = stencil_coefs_z.size();

        // Extract cylindrical specific parameters
        Real const dr = m_dr;
        int const nmodes = m_nmodes;
        Real const rmin = m_rmin;

        // Extract tileboxes for which to loop
        const Box& tbr  = mfi.tilebox(Bfield[0]->ixType().ixType());
        const Box& tbt  = mfi.tilebox(Bfield[1]->ixType().ixType());
        const Box& tbz  = mfi.tilebox(Bfield[2]->ixType().ixType());

        // Loop over the cells and update the fields
        amrex::ParallelFor(tbr, tbt, tbz,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Real const r = rmin + i*dr; // r on nodal point (Br is nodal in r)
                Br(i, j, 0, 0) += dt * T_Algo::UpwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 0); // Mode m=0
                for (int m=1; m<nmodes; m++) { // Higher-order modes
                    Br(i, j, 0, 2*m-1) += dt*(
                        T_Algo::UpwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 2*m-1)
                        - m * T_Algo::DivideByR(Ez, r, dr, m, i, j, 0, 2*m  ));  // Real part
                    Br(i, j, 0, 2*m  ) += dt*(
                        T_Algo::UpwardDz(Et, coefs_z, n_coefs_z, i, j, 0, 2*m  )
                        + m * T_Algo::DivideByR(Ez, r, dr, m, i, j, 0, 2*m-1)); // Imaginary part
                }
                // Ensure that Br remains 0 on axis (except for m=1)
                if (r==0) { // On axis
                    Br(i, j, 0, 0) = 0.; // Mode m=0
                    for (int m=2; m<nmodes; m++) { // Higher-order modes (but not m=1)
                        Br(i, j, 0, 2*m-1) = 0.;
                        Br(i, j, 0, 2*m  ) = 0.;
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bt(i, j, 0, 0) += dt*(
                    T_Algo::UpwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 0)
                    - T_Algo::UpwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 0)); // Mode m=0
                for (int m=1 ; m<nmodes ; m++) { // Higher-order modes
                    Bt(i, j, 0, 2*m-1) += dt*(
                        T_Algo::UpwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 2*m-1)
                        - T_Algo::UpwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 2*m-1)); // Real part
                    Bt(i, j, 0, 2*m  ) += dt*(
                        T_Algo::UpwardDr(Ez, coefs_r, n_coefs_r, i, j, 0, 2*m  )
                        - T_Algo::UpwardDz(Er, coefs_z, n_coefs_z, i, j, 0, 2*m  )); // Imaginary part
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Real const r = rmin + (i + 0.5)*dr; // r on a cell-centered grid (Bz is cell-centered in r)
                Bz(i, j, 0, 0) += dt*( - T_Algo::UpwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 0));
                for (int m=1 ; m<nmodes ; m++) { // Higher-order modes
                    Bz(i, j, 0, 2*m-1) += dt*( m * Er(i, j, 0, 2*m  )/r
                        - T_Algo::UpwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 2*m-1)); // Real part
                    Bz(i, j, 0, 2*m  ) += dt*(-m * Er(i, j, 0, 2*m-1)/r
                        - T_Algo::UpwardDrr_over_r(Et, r, dr, coefs_r, n_coefs_r, i, j, 0, 2*m  )); // Imaginary part
                }
            }

        );

    }

}

#endif // corresponds to ifndef WARPX_DIM_RZ
