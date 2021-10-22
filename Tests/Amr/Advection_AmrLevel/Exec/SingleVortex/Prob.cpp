
/**
 * Initalize problem
 *
 * \param level
 * \param time
 * \param phi
 * \param dx
 * \param prob_lo
 */


initdata(
    int level,
    Real time,
    const Array4 <Real>& phi,
    GpuArray dx,
    GpuArray prob_lo){

    int dm;

    // Detemine dimension
    // 2 if 
    // 3 if 
    if (  ) {
        dm = 2;
    } else {
        dm = 3;
    } 

    // Populate MultiFab data 
    Real x,y,z; 
    ParallelFor(box, [=] AMREX_GPU_DEVICE ( int i, int j, int k) noexcept 
    {
        z = prob_lo[2] + (k + 0.5) * dx[2];
        y = prob_lo[1] + (j + 0.5) * dx[1];
        x = prob_lo[0] + (i + 0.5) * dx[0];

        if ( dm == 2 ) {
            r2 = ((x-0.5)*(x-0.5) + (y-0.75)*(y-0.75)) / 0.01;
            phi(i,j,k) = 1.0 + std::exp(-r2);
        } else {
           r2 = ((x-0.5)*(x-0.5) + (y-0.75)*(y-0.75) + (z-0.5)*(z-0.5)) / 0.01;
           phi(i,j,k) = 1.0 + std::exp(-r2);
        }

    });

}

