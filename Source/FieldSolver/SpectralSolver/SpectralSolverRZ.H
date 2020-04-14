/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_SOLVER_RZ_H_
#define WARPX_SPECTRAL_SOLVER_RZ_H_

#include "SpectralAlgorithms/SpectralBaseAlgorithmRZ.H"
#include "SpectralFieldDataRZ.H"

/* \brief Top-level class for the electromagnetic spectral solver
 *
 * Stores the field in spectral space, and has member functions
 * to Fourier-transform the fields between real space and spectral space
 * and to update fields in spectral space over one time step.
 */
class SpectralSolverRZ
{
    public:
        // Inline definition of the member functions of `SpectralSolverRZ`,
        // except the constructor (see `SpectralSolverRZ.cpp`)
        // The body of these functions is short, since the work is done in the
        // underlying classes `SpectralFieldData` and `PsatdAlgorithm`

        // Constructor
        SpectralSolverRZ(amrex::BoxArray const & realspace_ba,
                         amrex::DistributionMapping const & dm,
                         int const n_rz_azimuthal_modes,
                         int const norder_z, bool const nodal,
                         amrex::RealVect const dx, amrex::Real const dt,
                         int const lev,
                         bool const pml=false);

        /* \brief Transform the component `i_comp` of MultiFab `field_mf`
         *  to spectral space, and store the corresponding result internally
         *  (in the spectral field specified by `field_index`) */
        void ForwardTransform (amrex::MultiFab const & field_mf,
                               int const field_index,
                               int const i_comp=0 ) {
            BL_PROFILE("SpectralSolverRZ::ForwardTransform");
            field_data.ForwardTransform(field_mf, field_index, i_comp);
        };

        /* \brief Transform the two MultiFabs `field_mf1` and `field_mf2`
         *  to spectral space, and store the corresponding results internally
         *  (in the spectral field specified by `field_index1` and `field_index2`) */
        void ForwardTransform (amrex::MultiFab const & field_mf1, int const field_index1,
                               amrex::MultiFab const & field_mf2, int const field_index2) {
            BL_PROFILE("SpectralSolverRZ::ForwardTransform");
            field_data.ForwardTransform(field_mf1, field_index1,
                                        field_mf2, field_index2);
        };

        /* \brief Transform spectral field specified by `field_index` back to
         * real space, and store it in the component `i_comp` of `field_mf` */
        void BackwardTransform (amrex::MultiFab& field_mf,
                                int const field_index,
                                int const i_comp=0) {
            BL_PROFILE("SpectralSolverRZ::BackwardTransform");
            field_data.BackwardTransform(field_mf, field_index, i_comp);
        };

        /* \brief Transform spectral fields specified by `field_index1` and `field_index2`
         * back to real space, and store it in `field_mf1` and `field_mf2`*/
        void BackwardTransform (amrex::MultiFab& field_mf1, int const field_index1,
                                amrex::MultiFab& field_mf2, int const field_index2) {
            BL_PROFILE("SpectralSolverRZ::BackwardTransform");
            field_data.BackwardTransform(field_mf1, field_index1,
                                         field_mf2, field_index2);
        };

        /* \brief Update the fields in spectral space, over one timestep */
        void pushSpectralFields () {
            BL_PROFILE("SpectralSolverRZ::pushSpectralFields");
            // Virtual function: the actual function used here depends
            // on the sub-class of `SpectralBaseAlgorithm` that was
            // initialized in the constructor of `SpectralSolverRZ`
            algorithm->pushSpectralFields(field_data);
        };

        /**
          * \brief Public interface to call the member function ComputeSpectralDivE
          * of the base class SpectralBaseAlgorithmRZ from objects of class SpectralSolverRZ
          */
        void ComputeSpectralDivE ( const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
                                   amrex::MultiFab& divE ) {
            algorithm->ComputeSpectralDivE( field_data, Efield, divE );
        };

    private:

        SpectralFieldDataRZ field_data; // Store field in spectral space
                                        // and perform the Fourier transforms
        std::unique_ptr<SpectralBaseAlgorithmRZ> algorithm;
        // Defines field update equation in spectral space,
        // and the associated coefficients.
        // SpectralBaseAlgorithmRZ is a base class ; this pointer is meant
        // to point an instance of a *sub-class* defining a specific algorithm

};

#endif // WARPX_SPECTRAL_SOLVER_RZ_H_
