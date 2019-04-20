#ifndef WARPX_SPECTRAL_SOLVER_H_
#define WARPX_SPECTRAL_SOLVER_H_

#include <SpectralKSpace.H>
#include <PsatdAlgorithm.H>
#include <SpectralFieldData.H>

/* \brief
 * TODO
 */
class SpectralSolver
{
    private:
        SpectralFieldData field_data; // Store field in spectral space
                                      // and perform the Fourier transforms
        PsatdAlgorithm algorithm; // Contains Psatd coefficients
                                  // and field update equation
    public:
        SpectralSolver( const amrex::BoxArray& realspace_ba,
                        const amrex::DistributionMapping& dm,
                        const int norder_x, const int norder_y, const int norder_z,
                        const amrex::Array<amrex::Real,3> dx, const amrex::Real dt )
        {
            // Initialize all structures using the same distribution mapping dm
            // - Initialize k space (Contains size of each box in spectral space,
            // and corresponding values of the k vectors)
            const SpectralKSpace k_space = SpectralKSpace( realspace_ba, dm, dx );
            // - Initialize algorithm (coefficients) over this space
            algorithm = PsatdAlgorithm( k_space, dm, norder_x, norder_y, norder_z, dt );
            // - Initialize arrays for fields in Fourier space + FFT plans
            field_data = SpectralFieldData( realspace_ba, k_space, dm );
        };
        void ForwardTransform( const amrex::MultiFab& mf,
                               const int field_index, const int i_comp=0 ){
            field_data.ForwardTransform( mf, field_index, i_comp );
        };
        void BackwardTransform( amrex::MultiFab& mf,
                                const int field_index, const int i_comp=0 ){
            field_data.BackwardTransform( mf, field_index, i_comp );
        };
        void pushSpectralFields(){
            algorithm.pushSpectralFields( field_data );
        };
};

#endif // WARPX_SPECTRAL_SOLVER_H_
