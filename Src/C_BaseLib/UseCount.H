//BL_COPYRIGHT_NOTICE

#ifndef _USECOUNT_H_
#define _USECOUNT_H_

//
// $Id: UseCount.H,v 1.1 1997-09-12 18:00:17 lijewski Exp $
//

#include <Assert.H>
#include <Boolean.H>
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

//
// Inlines.
//

inline
UseCount::UseCount ()
    : cnt(new unsigned int(1))
{
    if (cnt == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
}

inline
UseCount::UseCount (const UseCount& rhs)
    : cnt(rhs.cnt)
{
    ++*cnt;
}

inline
bool
UseCount::unique () const
{
    return *cnt == 1;
}

inline
void
UseCount::decrement ()
{
    if (unique())
    {
        delete cnt;
        cnt = 0;
    }
    else
        --*cnt;
}

inline
UseCount&
UseCount::operator= (const UseCount& rhs)
{
    ++*rhs.cnt;
    decrement();
    cnt = rhs.cnt;
    return *this;
}

inline
UseCount::~UseCount ()
{
    decrement();
}

inline
int
UseCount::linkCount () const
{
    return *cnt;
}

#endif /*_USECOUNT_H_*/
