#ifndef WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_
#define WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_

#include "WarpXParticleContainer.H"
//#include <WarpX.H>
#include <AMReX_DenseBins.H>
#include <AMReX_REAL.H>

/* \brief Shuffle array according to Fisher-Yates algorithm.
 *        Only shuffle the part between is <= i < ie, n = ie-is.*/

template <typename T_index>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ShuffleFisherYates(T_index *array, T_index const is, T_index const ie)
{
    int j;
    T_index buf;
    for (int i = ie-1; i >= is+1; --i)
    {
        // get random number is <= j <= i
        while ( true ) {
            j = (int) floor(amrex::Random()*(i-is+1)) + is;
            if ( j >= is && j <= i ) { break; }
        }
        buf = array[i];
        array[i] = array[j];
        array[j] = buf;
    }
}

//AMREX_GPU_HOST_DEVICE AMREX_INLINE
//void ShuffleFisherYates(
//    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type *array,
//    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const is,
//    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type const ie)
//{
//    int j;
//    amrex::DenseBins<WarpXParticleContainer::ParticleType>::index_type buf;
//    for (int i = (int) ie-1; i >= (int) is+1; --i)
//    {
//        // get random number is <= j <= i
//        while ( true ) {
//            j = (int) floor(amrex::Random()*(i-is+1)) + is;
//            if ( j >= (int) is && j <= (int) i ) { break; }
//        }
//        buf = array[i];
//        array[i] = array[j];
//        array[j] = buf;
//    }
//}

#endif // WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_
