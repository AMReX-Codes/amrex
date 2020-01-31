/* Copyright 2019 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_
#define WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_

#include <AMReX_Random.H>

/* \brief Shuffle array according to Fisher-Yates algorithm.
 *        Only shuffle the part between is <= i < ie, n = ie-is.
 *        T_index shall be
 *        amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type
*/

template <typename T_index>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ShuffleFisherYates (T_index *array, T_index const is, T_index const ie)
{
    int     j;
    T_index buf;
    for (int i = ie-1; i >= is+1; --i)
    {
        // get random number j: is <= j <= i
        j = amrex::Random_int(i-is+1) + is;
        // swop the ith array element with the jth
        buf      = array[i];
        array[i] = array[j];
        array[j] = buf;
    }
}

#endif // WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_
