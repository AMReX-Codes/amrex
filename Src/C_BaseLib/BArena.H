//BL_COPYRIGHT_NOTICE

#ifndef BL_BARENA_H
#define BL_BARENA_H

//
// $Id: BArena.H,v 1.4 1998-02-09 20:46:52 lijewski Exp $
//

#include <Arena.H>

//
//@Man:
//@Memo: A Concrete Class for Dynamic Memory Management
/*@Doc:

  This is the simplest dynamic memory management class derived from Arena.
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
