/* Copyright 2019-2020 Maxence Thevenet, Remi Lehe, Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_SOLVER_H_
#define WARPX_SPECTRAL_SOLVER_H_

#include "SpectralAlgorithms/SpectralBaseAlgorithm.H"
#include "SpectralFieldData.H"


#ifdef WARPX_USE_PSATD
/**
 * \brief Top-level class for the electromagnetic spectral solver
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
                        const amrex::Array<amrex::Real,3>& v_galilean,
                        const amrex::RealVect dx, const amrex::Real dt,
                        const bool pml=false );

        /**
         * \brief Transform the component `i_comp` of MultiFab `mf`
         *  to spectral space, and store the corresponding result internally
         *  (in the spectral field specified by `field_index`) */
        void ForwardTransform( const amrex::MultiFab& mf,
                               const int field_index,
                               const int i_comp=0 );

        /**
         * \brief Transform spectral field specified by `field_index` back to
         * real space, and store it in the component `i_comp` of `mf`
         */
        void BackwardTransform( amrex::MultiFab& mf,
                                const int field_index,
                                const int i_comp=0 );

        /**
         * \brief Update the fields in spectral space, over one timestep
         */
        void pushSpectralFields();

        /**
          * \brief Public interface to call the member function ComputeSpectralDivE
          * of the base class SpectralBaseAlgorithm from objects of class SpectralSolver
          */
        void ComputeSpectralDivE ( const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
                                   amrex::MultiFab& divE ) {
            algorithm->ComputeSpectralDivE( field_data, Efield, divE );
        };

    private:

        // Store field in spectral space and perform the Fourier transforms
        SpectralFieldData field_data;

        // Defines field update equation in spectral space and the associated coefficients.
        // SpectralBaseAlgorithm is a base class; this pointer is meant to point
        // to an instance of a sub-class defining a specific algorithm
        std::unique_ptr<SpectralBaseAlgorithm> algorithm;
};
#endif // WARPX_USE_PSATD
#endif // WARPX_SPECTRAL_SOLVER_H_
