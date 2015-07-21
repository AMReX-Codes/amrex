#ifndef BL_MEMPOOL_H
#define BL_MEMPOOL_H

extern "C" void mempool_init ();
extern "C" void* mempool_alloc (size_t n);
extern "C" void mempool_free (void* p);
extern "C" void mempool_get_stats (int& mp_min, int& mp_max, int& mp_tot);  // min, max & tot in MB

#endif
