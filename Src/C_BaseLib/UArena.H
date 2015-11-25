
#ifndef BL_UARENA_H
#define BL_UARENA_H

#include <Arena.H>
#include <upcxx.h>

//
// A Concrete Class for Dynamic Memory Management
//
// This uses UPC++'s allocator if the UPC++ team size is greater than 1,
// otherwise it just uses C++ new and delete.
//

class UArena
    :
    public Arena
{
public:
    //
    // Allocates a dynamic memory arena of size sz.
    // Returns a pointer to this memory.
    //
    virtual void* alloc (std::size_t sz) override;
    //
    // Deletes the arena pointed to by pt.
    //
    virtual void free (void* pt) override;
};

#endif /*BL_UARENA_H*/
