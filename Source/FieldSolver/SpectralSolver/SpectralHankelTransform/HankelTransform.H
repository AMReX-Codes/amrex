/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_HANKEL_TRANSFORM_H_
#define WARPX_HANKEL_TRANSFORM_H_

#include <AMReX_FArrayBox.H>

/* \brief This defines the class that performs the Hankel transform.
 * Original authors: Remi Lehe, Manuel Kirchen
 *
 *
 * Definition of the Hankel forward and backward transform of order p:
 * g(kr) = \int_0^\infty f(r) J_p(kr r) r dr
 * f(r ) = \int_0^\infty g(kr) J_p(kr r) kr dkr
*/
class HankelTransform
{
    public:

        using RealVector = amrex::Gpu::ManagedDeviceVector<amrex::Real>;

        // Constructor
        HankelTransform(const int hankel_order,
                        const int azimuthal_mode,
                        const int nr,
                        const amrex::Real rmax);

        const RealVector & getSpectralWavenumbers() {return m_kr;}

        void HankelForwardTransform(amrex::FArrayBox const& F, int const F_icomp,
                                    amrex::FArrayBox      & G, int const G_icomp);

        void HankelInverseTransform(amrex::FArrayBox const& G, int const G_icomp,
                                    amrex::FArrayBox      & F, int const F_icomp);

    private:
        // Even though nk == nr always, use a seperate variable for clarity.
        int m_nr, m_nk;

        RealVector m_kr;

        RealVector invM;
        RealVector M;
};

#endif
