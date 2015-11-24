
#ifndef BL_BARENA_H
#define BL_BARENA_H

#include <Arena.H>
#include <BL_CXX11.H>

//
// A Concrete Class for Dynamic Memory Management
//
// This is the simplest dynamic memory management class derived from Arena.
//
// Makes calls to ::operator new() and ::operator delete().
//

class BArena
    :
    public Arena
{
public:
    //
    // Allocates a dynamic memory arena of size sz.
    // Returns a pointer to this memory.
    //
    virtual void* alloc (std::size_t sz) BL_OVERRIDE;
    //
    // Deletes the arena pointed to by pt.
    //
    virtual void free (void* pt) BL_OVERRIDE;
};

#endif /*BL_BARENA_H*/
