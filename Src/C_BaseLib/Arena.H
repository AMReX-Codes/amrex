//BL_COPYRIGHT_NOTICE

#ifndef BL_ARENA
#define BL_ARENA

//
// $Id: Arena.H,v 1.3 1997-09-18 20:12:42 lijewski Exp $
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

    union Word
    {
        void*  p;
        double d;
        long   l;
    };

    /*@ManDoc: A pure virtual function for allocating a dynamic memory arena
               of size sz.  A pointer to this memory should be returned.
               Also the *pt is set to the returned pointer if pt != 0;
    */
    virtual void *alloc (size_t sz, void** pt = 0) = 0;

    /*@ManDoc: A pure virtual function for deleting the arena pointed
               to by pt.
    */
    virtual void free (void* pt) = 0;

    /*@ManDoc: Given a minimum required arena size of sz bytes, this returns
               the next largest arena size that will hold an integral number
               of objects of the largest of the three types void*, long
               and double.
    */
    static size_t align (size_t sz);
};

//
// Inlines.
//

inline
size_t
Arena::align (size_t s)
{
    int x = s + sizeof(Word) - 1;
    x -= x%sizeof(Word);
    return x;
}

#endif /*BL_ARENA*/
