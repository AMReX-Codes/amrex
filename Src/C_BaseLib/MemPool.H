#ifndef BL_MEMPOOL_H
#define BL_MEMPOOL_H

extern void init_memmory_pool ();
extern "C" void* bl_allocate_c (size_t n);
extern "C" void bl_deallocate_c (void* p);

#endif
