//BL_COPYRIGHT_NOTICE

#ifndef _BARENA_H_
#define _BARENA_H_

//
// $Id: BArena.H,v 1.1 1997-09-12 18:00:01 lijewski Exp $
//

#include <Arena.H>

//
//@Man:
//@Memo: A Concrete Class for Dynamic Memory Management
/*@Doc:

  This is the simplest dynamic memory management class derived from
  Arena.  It is used to manage the memory for class SimpleDPtr.
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

#endif /*_BARENA_H_*/
