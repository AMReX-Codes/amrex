#ifdef AMREX_USE_HDF5_ASYNC
#include "hdf5.h"

namespace amrex {

/**
 * \file AMReX_ParticleUtilHDF5.H
 *
 * Async-HDF5 helpers shared by the particle checkpoint/write routines.
 */

//! Global async VOL event-set used to synchronize particle HDF5 I/O.
hid_t es_par_g = 0;

/**
 * \brief Block until all pending particle async VOL operations complete.
 */
void async_vol_es_wait_particle()
{
    size_t num_in_progress;
    hbool_t op_failed;
    if (es_par_g != 0) {
        H5ESwait(es_par_g, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
        if (num_in_progress != 0)
            std::cout << "After H5ESwait, still has async operations in progress!" << std::endl;
    }
    return;
}

/**
 * \brief Wait for pending async VOL work and destroy the shared event-set.
 */
void async_vol_es_wait_close_particle()
{
    size_t num_in_progress;
    hbool_t op_failed;
    if (es_par_g != 0) {
        H5ESwait(es_par_g, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
        if (num_in_progress != 0)
            std::cout << "After H5ESwait, still has async operations in progress!" << std::endl;
        H5ESclose(es_par_g);
        /* std::cout << "es_par_g closed!" << std::endl; */
        es_par_g = 0;
    }
    return;
}

}

#endif
