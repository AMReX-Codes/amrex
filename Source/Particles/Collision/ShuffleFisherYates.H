#ifndef WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_
#define WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_

#include <AMReX_REAL.H>

/* \brief Shuffle array according to Fisher-Yates algorithm.
 *        Only shuffle the part between is <= i < ie, n = ie-is.*/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ShuffleFisherYates(int *array, const int is, const int ie)
{
    int j; int buf;
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

#endif // WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_
