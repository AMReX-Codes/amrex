//BL_COPYRIGHT_NOTICE

#ifndef BL_ARENA
#define BL_ARENA

//
// $Id: Arena.H,v 1.5 1998-02-09 20:46:52 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstddef>
#else
#include <stddef.h>
#endif

//
//@Man:
//@Memo: A Virtual Base Class for Dynamic Memory Management
/*@Doc:

  This is a virtual base class for objects that manage their own dynamic
  memory allocation.  Since it is a virtual base class, you have to derive
  something from it to use it.
*/

class Arena
{
public:
    /*@ManDoc: Allocate a dynamic memory arena of size sz.
               A pointer to this memory should be returned.
    */
    virtual void* alloc (size_t sz) = 0;

    /*@ManDoc: A pure virtual function for deleting the arena pointed
               to by pt.
    */
    virtual void free (void* pt) = 0;

    /*@ManDoc: Given a minimum required arena size of sz bytes, this returns
               the next largest arena size that will hold an integral number
               of objects of the largest of the types void*, long,
               double and function pointer.
    */
    static size_t align (size_t sz);

protected:
    //
    // Types used by align().
    //
    typedef void (*FP) ();

    union Word
    {
        void*  p;
        double d;
        long   l;
        FP     f;
    };
};

//
// Inlines.
//

inline
size_t
Arena::align (size_t s)
{
    size_t x = s + sizeof(Word) - 1;
    x -= x%sizeof(Word);
    return x;
}

#endif /*BL_ARENA*/
