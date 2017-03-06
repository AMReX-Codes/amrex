#ifndef BL_MEMPOOL_H
#define BL_MEMPOOL_H

extern "C" {
    void  mempool_init ();
    void* mempool_alloc (size_t n);
    void  mempool_free (void* p);
    void  mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot);  // min, max & tot in MB
    // We should eventually use amre_real instead of double.
    // We cannot do it now because of F_BaseLib.
    void  real_array_init (double* p, size_t nelems); 
    void  array_init_snan (double* p, size_t nelems);
}

#endif
