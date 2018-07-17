#ifndef BL_MEMPOOL_H
#define BL_MEMPOOL_H

extern "C" {
    void  amrex_mempool_init ();
    void  amrex_mempool_finalize ();
    void* amrex_mempool_alloc (size_t n);
    void  amrex_mempool_free (void* p);
    void  amrex_mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot);  // min, max & tot in MB
    // We should eventually use amre_real instead of double.
    // We cannot do it now because of F_BaseLib.
    void  amrex_real_array_init (double* p, size_t nelems); 
    void  amrex_array_init_snan (double* p, size_t nelems);
}

#endif
