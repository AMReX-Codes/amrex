
#ifndef BL_BARENA_H
#define BL_BARENA_H

//
// $Id: BArena.H,v 1.7 2000-10-02 20:52:32 lijewski Exp $
//

#include <Arena.H>

//
//@Man:
//@Memo: A Concrete Class for Dynamic Memory Management
/*@Doc:

  This is the simplest dynamic memory management class derived from Arena.

  Makes calls to ::operator new() and ::operator delete().
*/

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

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

#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif

#endif /*BL_BARENA_H*/
