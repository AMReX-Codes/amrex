// TODO include statements

FiniteDifferenceSolver::EvolveB ( VectorField Bfield,
                                  ConstVectorField Efield,
                                  amrex::Real dt ) {
    // Select algorithm (The choice of algorithm is a runtime option,
    // but we compile code for each algorithm, using templates)
    if (fdtd_algo == MaxwellSolverAlgo::Yee){
        EvolveBwithAlgo <YeeAlgorithm> ( Bfield, Efield, dt );
    } else if (fdtd_algo == MaxwellSolverAlgo::CKC) {
        EvolveBwithAlgo <CKCAlgorithm> ( Bfield, Efield, dt );
    } else {
        amrex::Abort("Unknown algorithm");
    }
)

template<typename algo>
FiniteDifferenceSolver::EvolveBwithAlgo ( VectorField Bfield,
                                          ConstVectorField Efield,
                                          amrex::Real dt ) {

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bx, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Extract field data for this grid/tile
        auto const& Bx = Bfield[0]->array(mfi);
        auto const& By = Bfield[1]->array(mfi);
        auto const& Bz = Bfield[2]->array(mfi);
        auto const& Ex = Efield[0]->array(mfi);
        auto const& Ey = Efield[1]->array(mfi);
        auto const& Ez = Efield[2]->array(mfi);

        // Extract stencil coefficients
        Real const* AMREX_RESTRICT coefs_x = stencil_coefs_x.dataPtr();
        Real const* AMREX_RESTRICT coefs_y = stencil_coefs_y.dataPtr();
        Real const* AMREX_RESTRICT coefs_z = stencil_coefs_z.dataPtr();

        // Extract tileboxes for which to loop
        const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
        const Box& tby  = mfi.tilebox(By_nodal_flag);
        const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

        // Loop over the cells and update the fields
        amrex::ParallelFor(tbx, tby, tbz,

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bx(i, j, k) += dt * algo::UpwardDz( Ey, i, j, k, coefs_z)
                             - dt * algo::UpwardDy( Ez, i, j, k, coefs_y);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                By(i, j, k) += dt * algo::UpwardDx( Ez, i, j, k, coefs_x)
                             - dt * algo::UpwardDz( Ex, i, j, k, coefs_z);
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k){
                Bz(i, j, k) += dt * algo::UpwardDy( Ex, i, j, k, coefs_y)
                             - dt * algo::UpwardDx( Ey, i, j, k, coefs_x);
            }

        );

    }

};
