#include <WarpXAlgorithmSelection.H>
#include<FiniteDifferenceAlgorithms/YeeAlgorithm.H>
#include<FiniteDifferenceSolver.H>
#include<AMReX_Gpu.H>

using namespace amrex;

void FiniteDifferenceSolver::EvolveB ( VectorField& Bfield,
                                       VectorField& Efield,
                                       amrex::Real dt ) {
    // Select algorithm (The choice of algorithm is a runtime option,
    // but we compile code for each algorithm, using templates)
    if (m_fdtd_algo == MaxwellSolverAlgo::Yee){
        EvolveBwithAlgo <YeeAlgorithm> ( Bfield, Efield, dt );
//    } else if (fdtd_algo == MaxwellSolverAlgo::CKC) {
//       EvolveBwithAlgo <CKCAlgorithm> ( Bfield, Efield, dt );
    } else {
        amrex::Abort("Unknown algorithm");
    }
}

template<typename algo>
void FiniteDifferenceSolver::EvolveBwithAlgo ( VectorField& Bfield,
                                               VectorField& Efield,
                                               amrex::Real dt ) {

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
                Bx(i, j, k) += dt * algo::UpwardDz(Ey, coefs_z, n_coefs_z, i, j, k)
                             - dt * algo::UpwardDy(Ez, coefs_y, n_coefs_y, i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                By(i, j, k) += dt * algo::UpwardDx(Ez, coefs_x, n_coefs_x, i, j, k)
                             - dt * algo::UpwardDz(Ex, coefs_z, n_coefs_z, i, j, k);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bz(i, j, k) += dt * algo::UpwardDy(Ex, coefs_y, n_coefs_y, i, j, k)
                             - dt * algo::UpwardDx(Ey, coefs_x, n_coefs_x, i, j, k);
            }

        );

    }

}
