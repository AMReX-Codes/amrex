
#ifndef BL_ARENA_H
#define BL_ARENA_H

#include <cstddef>

namespace amrex {

class Arena;
class BArena;
class CArena;

void Arena_Initialize();
void Arena_Finalize();

#ifdef BL_COALESCE_FABS
CArena* The_Arena ();
#else
BArena* The_Arena ();
#endif

/**
* \brief 
* A virtual base class for objects that manage their own dynamic
* memory allocation.  
*/

class Arena
{
public:

    ~Arena ();
    /**
    * \brief Given a minimum required arena size of sz bytes, this returns
    * the next largest arena size that will align to align_size bytes
    */
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

}

#endif /*BL_ARENA_H*/
