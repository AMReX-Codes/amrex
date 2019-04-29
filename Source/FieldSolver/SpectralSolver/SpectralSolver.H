#ifndef WARPX_SPECTRAL_SOLVER_H_
#define WARPX_SPECTRAL_SOLVER_H_

#include <SpectralKSpace.H>
#include <PsatdAlgorithm.H>
#include <SpectralFieldData.H>

/* \brief Top-level class for the electromagnetic spectral solver
 *
 * Stores the field in spectral space, and has member functions
 * to Fourier-transform the fields between real space and spectral space
 * and to update fields in spectral space over one time step.
 */
class SpectralSolver
{
    public:
        // Inline definition of the member functions of `SpectralSolver`
        // The body of these functions is short, since the work is done in the
        // underlying classes `SpectralFieldData` and `PsatdAlgorithm`

        /* \brief Initialize the spectral solver */
        SpectralSolver( const amrex::BoxArray& realspace_ba,
                        const amrex::DistributionMapping& dm,
                        const int norder_x, const int norder_y,
                        const int norder_z, const bool nodal,
                        const amrex::RealVect dx, const amrex::Real dt ) {
            // Initialize all structures using the same distribution mapping dm

            // - Initialize k space object (Contains info about the size of
            // the spectral space corresponding to each box in `realspace_ba`,
            // as well as the value of the corresponding k coordinates)
            const SpectralKSpace k_space= SpectralKSpace(realspace_ba, dm, dx);
            // - Initialize the algorithm (coefficients) over this space
            algorithm = PsatdAlgorithm( k_space, dm, norder_x, norder_y,
                                        norder_z, nodal, dt );
            // - Initialize arrays for fields in Fourier space + FFT plans
            field_data = SpectralFieldData( realspace_ba, k_space, dm );
        };

        /* \brief Transform the component `i_comp` of MultiFab `mf`
         *  to spectral space, and store the corresponding result internally
         *  (in the spectral field specified by `field_index`) */
        void ForwardTransform( const amrex::MultiFab& mf,
                               const int field_index,
                               const int i_comp=0 ){
            BL_PROFILE("SpectralSolver::ForwardTransform");
            field_data.ForwardTransform( mf, field_index, i_comp );
        };

        /* \brief Transform spectral field specified by `field_index` back to
         * real space, and store it in the component `i_comp` of `mf` */
        void BackwardTransform( amrex::MultiFab& mf,
                                const int field_index,
                                const int i_comp=0 ){
            BL_PROFILE("SpectralSolver::BackwardTransform");
            field_data.BackwardTransform( mf, field_index, i_comp );
        };

        /* \brief Update the fields in spectral space, over one timestep */
        void pushSpectralFields(){
            BL_PROFILE("SpectralSolver::pushSpectralFields");
            algorithm.pushSpectralFields( field_data );
        };

    private:
        SpectralFieldData field_data; // Store field in spectral space
                                      // and perform the Fourier transforms
        PsatdAlgorithm algorithm; // Contains Psatd coefficients
                                  // and field update equation
};

#endif // WARPX_SPECTRAL_SOLVER_H_
