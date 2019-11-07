#ifndef WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_
#define WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_

#include <AMReX_REAL.H>

/* \brief Shuffle array according to Fisher-Yates algorithm.*/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ShuffleFisherYates(int *array, const int n)
{
    int j; int buf;
    for (int i = n-1; i >= 1; --i)
    {
        // get random number 0 <= j <= i
        while ( true ) {
            j = (int) floor( amrex::Random() * (i+1) );
	    if ( j >= 0 && j <= i ) { break; }
        }
        buf = array[i];
        array[i] = array[j];
        array[j] = buf;
    }
}

#endif // WARPX_PARTICLES_COLLISION_SHUFFLE_FISHER_YATES_H_
