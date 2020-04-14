/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRALHANKELTRANSFORMER_H_
#define WARPX_SPECTRALHANKELTRANSFORMER_H_

#include <AMReX_FArrayBox.H>
#include "HankelTransform.H"

/* \brief Object that allows to transform the fields back and forth between the
 *  spectral and interpolation grid.
 *
 *  Attributes :
 *  - dht0, dhtm, dhtp : the discrete Hankel transform objects for the modes,
 *     operating along r
*/

class SpectralHankelTransformer
{

    public:

        SpectralHankelTransformer() {};

        SpectralHankelTransformer(const int nr,
                                  const int n_rz_azimuthal_modes,
                                  const amrex::Real rmax);

        void
        ExtractKrArray();

        // Returns an array that holds the kr for all of the modes
        HankelTransform::RealVector const & getKrArray() const {return m_kr;}

        // Converts a scalar field from the physical to the spectral space
        void
        PhysicalToSpectral_Scalar(amrex::Box const & box,
                                  amrex::FArrayBox const & F_physical,
                                  amrex::FArrayBox       & G_spectral);

        // Converts a vector field from the physical to the spectral space
        void
        PhysicalToSpectral_Vector(amrex::Box const & box,
                                  amrex::FArrayBox & F_r_physical,
                                  amrex::FArrayBox & F_t_physical,
                                  amrex::FArrayBox & G_p_spectral,
                                  amrex::FArrayBox & G_m_spectral);

        // Converts a scalar field from the spectral to the physical space
        void
        SpectralToPhysical_Scalar(amrex::Box const & box,
                                  amrex::FArrayBox const & G_spectral,
                                  amrex::FArrayBox       & F_physical);

        // Converts a vector field from the spectral to the physical space
        void
        SpectralToPhysical_Vector(amrex::Box const & box,
                                  amrex::FArrayBox const & G_p_spectral,
                                  amrex::FArrayBox const & G_m_spectral,
                                  amrex::FArrayBox       & F_r_physical,
                                  amrex::FArrayBox       & F_t_physical);

    private:

        int m_nr;
        int m_n_rz_azimuthal_modes;
        HankelTransform::RealVector m_kr;

        amrex::Vector< std::unique_ptr<HankelTransform> > dht0;
        amrex::Vector< std::unique_ptr<HankelTransform> > dhtm;
        amrex::Vector< std::unique_ptr<HankelTransform> > dhtp;

};

#endif
