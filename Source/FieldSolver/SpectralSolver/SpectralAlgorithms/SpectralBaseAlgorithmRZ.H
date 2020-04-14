/* Copyright 2019 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_BASE_ALGORITHM_RZ_H_
#define WARPX_SPECTRAL_BASE_ALGORITHM_RZ_H_

#include "FieldSolver/SpectralSolver/SpectralKSpaceRZ.H"
#include "FieldSolver/SpectralSolver/SpectralFieldDataRZ.H"

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 *
 * `SpectralBaseAlgorithmRZ` is only a base class and cannot be used directly.
 * Instead use its subclasses, which implement the specific field update
 * equations for a given spectral algorithm.
 */
class SpectralBaseAlgorithmRZ
{
    public:
        // Virtual member function ; meant to be overridden in subclasses
        virtual void pushSpectralFields(SpectralFieldDataRZ & f) = 0;
        virtual int getRequiredNumberOfFields() const = 0;
        // The destructor should also be a virtual function, so that
        // a pointer to subclass of `SpectraBaseAlgorithm` actually
        // calls the subclass's destructor.
        virtual ~SpectralBaseAlgorithmRZ() {};

        /**
         * \brief Compute spectral divergence of E
         */
        void ComputeSpectralDivE ( SpectralFieldDataRZ& field_data,
                                   const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
                                   amrex::MultiFab& divE );

    protected: // Meant to be used in the subclasses

        using SpectralCoefficients = amrex::FabArray< amrex::BaseFab <amrex::Real> >;

        // Constructor
        SpectralBaseAlgorithmRZ(SpectralKSpaceRZ const & spectral_kspace,
                                amrex::DistributionMapping const & dm,
                                int const norder_z, bool const nodal)
          // Compute and assign the modified k vectors
          : modified_kz_vec(spectral_kspace.getModifiedKComponent(dm, 1, norder_z, nodal))
          {};

        // Modified finite-order vectors
        KVectorComponent modified_kz_vec;
};

#endif // WARPX_SPECTRAL_BASE_ALGORITHM_RZ_H_
