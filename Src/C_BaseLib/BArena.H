#ifndef BL_BARENA_H
#define BL_BARENA_H

//
// $Id: BArena.H,v 1.9 2001-07-19 16:57:30 lijewski Exp $
//

#include <Arena.H>

//
//@Man:
//@Memo: A Concrete Class for Dynamic Memory Management
/*@Doc:

  This is the simplest dynamic memory management class derived from Arena.

  Makes calls to ::operator new() and ::operator delete().
*/

class BArena
    :
    public Arena
{
public:

    /*@ManDoc: Allocates a dynamic memory arena of size sz.  Returns a
               pointer to this memory.
    */
    virtual void* alloc (size_t sz);
    //
    //@ManDoc: Deletes the arena pointed to by pt.
    //
    virtual void free (void* pt);
};

#endif /*BL_BARENA_H*/
