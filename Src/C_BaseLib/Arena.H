#ifndef BL_ARENA_H
#define BL_ARENA_H
//
// $Id: Arena.H,v 1.13 2001-10-16 19:59:51 lijewski Exp $
//
#include <cstddef>

class Arena;

namespace BoxLib
{
    Arena* The_Arena ();
}

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
    virtual ~Arena ();

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
    union Word
    {
        void*  p;
        double d;
        long   l;
        void (*f) ();
    };
};

#endif /*BL_ARENA_H*/
