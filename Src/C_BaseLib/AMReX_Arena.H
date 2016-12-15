
#ifndef BL_ARENA_H
#define BL_ARENA_H

#include <winstd.H>
#include <cstddef>

class Arena;

namespace BoxLib
{
    Arena* The_Arena ();
}

//
// A Virtual Base Class for Dynamic Memory Management
//
// This is a virtual base class for objects that manage their own dynamic
// memory allocation.  Since it is a virtual base class, you have to derive
// something from it to use it.
//

class Arena
{
public:

    virtual ~Arena ();
    //
    // Allocate a dynamic memory arena of size sz.
    // A pointer to this memory should be returned.
    //
    virtual void* alloc (std::size_t sz) = 0;
    //
    // A pure virtual function for deleting the arena pointed to by pt.
    //
    virtual void free (void* pt) = 0;
    //
    // Given a minimum required arena size of sz bytes, this returns
    // the next largest arena size that will align to align_size bytes.
    //
    static std::size_t align (std::size_t sz);

protected:

#if 0
    union Word
    {
        void*  p;
        double d;
        long   l;
        void (*f) ();
    };
    static const unsigned int align_size = sizeof(Word);
#endif

    static const unsigned int align_size = 16;
};

#endif /*BL_ARENA_H*/
