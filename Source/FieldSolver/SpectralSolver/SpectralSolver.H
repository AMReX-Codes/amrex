#ifndef WARPX_SPECTRAL_SOLVER_H_
#define WARPX_SPECTRAL_SOLVER_H_

#include <SpectralBaseAlgorithm.H>
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
        // Inline definition of the member functions of `SpectralSolver`,
        // except the constructor (see `SpectralSolver.cpp`)
        // The body of these functions is short, since the work is done in the
        // underlying classes `SpectralFieldData` and `PsatdAlgorithm`

        // Constructor
        SpectralSolver( const amrex::BoxArray& realspace_ba,
                        const amrex::DistributionMapping& dm,
                        const int norder_x, const int norder_y,
                        const int norder_z, const bool nodal,
                        const amrex::RealVect dx, const amrex::Real dt );

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
            // Virtual function: the actual function used here depends
            // on the sub-class of `SpectralBaseAlgorithm` that was
            // initialized in the constructor of `SpectralSolver`
            algorithm->pushSpectralFields( field_data );
        };

    private:
        SpectralFieldData field_data; // Store field in spectral space
                                      // and perform the Fourier transforms
        std::unique_ptr<SpectralBaseAlgorithm> algorithm;
        // Defines field update equation in spectral space,
        // and the associated coefficients.
        // SpectralBaseAlgorithm is a base class ; this pointer is meant
        // to point an instance of a *sub-class* defining a specific algorithm
};

#endif // WARPX_SPECTRAL_SOLVER_H_
