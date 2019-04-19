#ifndef WARPX_SPECTRAL_SOLVER_H_
#define WARPX_SPECTRAL_SOLVER_H_

/* \brief
 * TODO
 */
class SpectralSolver
{
    private:
        SpectralKSpace k_space; // Contains size of each box in spectral space,
                                // and corresponding values of the k vectors
        SpectralFieldData field_data; // Store field in spectral space
                                      // and perform the Fourier transforms
        PsatdAlgorithm algorithm; // Contains Psatd coefficients
                                  // and field update equation
    public:
        SpectralSolver();
        void ForwardTransform( const MultiFab& mf, const int field_index ){
            field_data.ForwardTransform( mf, field_index );
        };
        void BackwardTransform( MultiFab& mf, const int field_index ){
            field_data.BackwardTransform( mf, field_index );
        };
        void pushSpectralFields(){
            algorithm.pushSpectralFields( field_data );
        };
};

SpectralSolver::SpectralSolver( const SpectralKSpace& realspace_ba,
                                const DistributionMapping& dm,
                                const int norder_x, const int norder_y,
                                const int norder_z, const Real* dx,
                                const Real dt )
{
    // Initialize all structures using the same distribution mapping dm
    // - Initialize k space values
    k_space = SpectralKSpace( realspace_ba, dm, dx );
    // - Initialize algorithm (coefficients) over this space
    algorithm = PsatdAlgorithm( k_space, dm, norder_x, norder_y, norder_z );
    // - Initialize arrays for fields in Fourier space + FFT plans
    field_data = SpectralFieldData( realspace_ba, k_space, dm );
};

#endif // WARPX_SPECTRAL_SOLVER_H_
