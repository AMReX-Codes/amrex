#ifndef BL_ARENA_H
#define BL_ARENA_H
//
// $Id: Arena.H,v 1.16 2004-01-07 21:44:20 car Exp $
//
#include <winstd.H>

#include <cstddef>

class Arena;

namespace BoxLib
{
    Arena* The_Arena ();

    Arena* ResetArena (Arena*);
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
    virtual void* alloc (std::size_t sz) = 0;

    /*@ManDoc: A pure virtual function for deleting the arena pointed
               to by pt.
    */
    virtual void free (void* pt) = 0;

    /*@ManDoc: Given a minimum required arena size of sz bytes, this returns
               the next largest arena size that will hold an integral number
               of objects of the largest of the types void*, long,
               double and function pointer.
    */
    static std::size_t align (std::size_t sz);

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
