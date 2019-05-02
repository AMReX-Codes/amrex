#ifndef WARPX_SPECTRAL_BASE_ALGORITHM_H_
#define WARPX_SPECTRAL_BASE_ALGORITHM_H_

#include <SpectralKSpace.H>
#include <SpectralFieldData.H>

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 *
 * `SpectralBaseAlgorithm` is only a base class and cannot be used directly.
 * Instead use its subclasses, which implement the specific field update
 * equations for a given spectral algorithm.
 */
class SpectralBaseAlgorithm
{
    public:
        // Member function that updates the fields in spectral space ;
        // meant to be overridden in subclasses
        virtual void pushSpectralFields(SpectralFieldData& f) const = 0;
        // The destructor should also be a virtual function, so that
        // a pointer to subclass of `SpectraBaseAlgorithm` actually
        // calls the subclass's destructor.
        virtual ~SpectralBaseAlgorithm() {};

    protected: // Meant to be used in the subclasses

        using SpectralCoefficients = amrex::FabArray< amrex::BaseFab <amrex::Real> >;

        // Constructor
        SpectralBaseAlgorithm(const SpectralKSpace& spectral_kspace,
                              const amrex::DistributionMapping& dm,
                              const int norder_x, const int norder_y,
                              const int norder_z, const bool nodal)
          // Compute and assign the modified k vectors
          : modified_kx_vec(spectral_kspace.getModifiedKComponent(dm,0,norder_x,nodal)),
#if (AMREX_SPACEDIM==3)
            modified_ky_vec(spectral_kspace.getModifiedKComponent(dm,1,norder_y,nodal)),
            modified_kz_vec(spectral_kspace.getModifiedKComponent(dm,2,norder_z,nodal))
#else
            modified_kz_vec(spectral_kspace.getModifiedKComponent(dm,1,norder_z,nodal))
#endif
          {};

        // Modified finite-order vectors
        KVectorComponent modified_kx_vec, modified_kz_vec;
#if (AMREX_SPACEDIM==3)
        KVectorComponent modified_ky_vec;
#endif
};

#endif // WARPX_SPECTRAL_BASE_ALGORITHM_H_
