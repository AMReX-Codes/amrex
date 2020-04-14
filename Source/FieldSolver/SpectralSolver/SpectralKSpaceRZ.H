/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_K_SPACE_RZ_H_
#define WARPX_SPECTRAL_K_SPACE_RZ_H_

#include "SpectralKSpace.H"

/* \brief Class that represents the spectral, Hankel/FFT, space.
 *
 * (Contains info about the size of the spectral space corresponding
 * to each box in `realspace_ba`, as well as the value of the
 * corresponding kz coordinates)
 */
class SpectralKSpaceRZ
    :
    public SpectralKSpace
{
    public:
        SpectralKSpaceRZ(const amrex::BoxArray& realspace_ba,
                         const amrex::DistributionMapping& dm,
                         const amrex::RealVect realspace_dx);

        KVectorComponent const & getKzArray () const {return k_vec[1];}
        amrex::RealVect const & getCellSize () const {return dx;}

};

#endif // WARPX_SPECTRAL_K_SPACE_RZ_H_
