//BL_COPYRIGHT_NOTICE

#ifndef BL_BARENA_H
#define BL_BARENA_H

//
// $Id: BArena.H,v 1.3 1997-09-17 22:51:12 lijewski Exp $
//

#include <Arena.H>

//
//@Man:
//@Memo: A Concrete Class for Dynamic Memory Management
/*@Doc:

  This is the simplest dynamic memory management class derived from Arena.
*/

class BArena
    : public Arena
{
public:

    /*@ManDoc: Allocates a dynamic memory arena of size sz.  Returns a
               pointer to this memory.  Also stores in it in *pt, if pt != 0.
    */
    virtual void* alloc (size_t sz,
                         void** pt);
    //
    //@ManDoc: Deletes the arena pointed to by pt.
    //
    virtual void free (void* pt);
};

#endif /*BL_BARENA_H*/
