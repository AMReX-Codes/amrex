#ifndef BL_USECOUNT_H
#define BL_USECOUNT_H
//
// $Id: UseCount.H,v 1.13 2001-07-31 22:43:19 lijewski Exp $
//
#include <BLassert.H>
#include <BoxLib.H>

//
//@Man:
//@Memo: A Class Encapsulating Reference Counts
/*@Doc:

  This class encapsulates reference counts.
*/

class UseCount
{
public:
    //
    //@ManDoc: Construct a UseCount initialized to one.
    //
    UseCount ();
    //
    //@ManDoc: The copy constructor -- bumps reference count.
    //
    UseCount (const UseCount& rhs);
 
    /*@ManDoc: The assignment operator.  Increments the reference count on
               rhs, decrements the reference count on this, and then makes
               this and rhs point to the same reference count.
    */
    UseCount& operator= (const UseCount& rhs);
    //
    //@ManDoc: The destructor -- decrements the reference count.
    //
    ~UseCount ();
    //
    //@ManDoc: The UseCount is unique if the reference count == 1.
    //
    bool unique () const;
    //
    //@ManDoc: Returns the reference count.
    //
    int linkCount () const;

private:
    //
    // A pointer to the reference count.
    //
    unsigned int* cnt;
    //
    // Decrement the reference count and delete the reference
    // counter if there are no more references.
    //
    void decrement ();
};

inline
bool
UseCount::unique () const
{
    return *cnt == 1;
}

inline
int
UseCount::linkCount () const
{
    return *cnt;
}

#endif /*BL_USECOUNT_H*/
