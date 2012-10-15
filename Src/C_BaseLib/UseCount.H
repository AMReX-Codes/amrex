
#ifndef BL_USECOUNT_H
#define BL_USECOUNT_H

#include <BLassert.H>
#include <BoxLib.H>

//
// A Class Encapsulating Reference Counts
//
// This class encapsulates reference counts.
//

class UseCount
{
public:
    //
    // Construct a UseCount initialized to one.
    //
    UseCount ();
    //
    // The copy constructor -- bumps reference count.
    //
    UseCount (const UseCount& rhs);
    //
    // The assignment operator.  Increments the reference count on
    // rhs, decrements the reference count on this, and then makes
    // this and rhs point to the same reference count.
    //
    UseCount& operator= (const UseCount& rhs);
    //
    // The destructor -- decrements the reference count.
    //
    ~UseCount ();
    //
    // The UseCount is unique if the reference count == 1.
    //
    bool unique () const { return *cnt == 1; }
    //
    // Returns the reference count.
    //
    int linkCount () const { return *cnt; }

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

#endif /*BL_USECOUNT_H*/
