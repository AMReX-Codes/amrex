//BL_COPYRIGHT_NOTICE

#ifndef _ALIASEDDPTR_H_
#define _ALIASEDDPTR_H_

//
// $Id: AliasedDPtr.H,v 1.1 1997-09-12 17:59:59 lijewski Exp $
//

#include <Assert.H>
#include <stddef.h>
#include <DPtr.H>
#include <SimpleDPtr.H>
#include <Pointers.H>

//
//@Man:
//@Memo: A Shared and Counted DPtr
/*@Doc:

  This class uses a `LnClassPtr<SimpleDPtr<T>>' to allow several Dptrs
  to address the same area of memory.

  This class does NOT provide a copy constructor or assignment operator.
 */

template <class T>
class AliasedDPtr
    : public DPtrRep<T>
{
public:

    /*@ManDoc: Construct a new internal LnClassPtr<SimpleDPtr<T>> of
               the specifed size.
    */
    AliasedDPtr (size_t size);

    /*@ManDoc: Construct a new object referencing the same internal
               representation as the old one.  Store size as the current
               size regardless of the size of *dptr, and have this object
               point to a location in memory offset units (the size of the
               type T) from the location to which the internal object points.
    */
    AliasedDPtr (AliasedDPtr<T>* dptr,
                 size_t          size,
                 int             offset = 0);
    //
    //@ManDoc: The destructor.
    //
    virtual ~AliasedDPtr ();

    /*@ManDoc: Smart array access, taking into account the current memory
               offset.
    */
    virtual T& operator[] (long) const;

    /*@ManDoc: Decrement the count on the internal object, remove it if this
               is the last alias, and set this object to point to null.
    */
    virtual void clear ();

    /*@ManDoc: Disassociate this object from the current internal
               representation and then create a new one of the specified
               size.
    */
    virtual void resize (size_t newsize);

    /*@ManDoc: Create in initial internal object.  The size is specified
               as the number of objects of type T that the arena will hold.
    */
    virtual void define (size_t size);
    //
    //@ManDoc: Returns the current arena size, in units of the size of type T.
    //
    virtual size_t size () const;
protected:
    LnClassPtr< SimpleDPtr<T> > dptr;
    int                         offset;
    size_t                      currentsize;
private:
    //
    // These are disallowed.
    //
    AliasedDPtr (const AliasedDPtr<T> &);
    AliasedDPtr<T>& operator= (const AliasedDPtr<T>&);
};


template <class T>
AliasedDPtr<T>::AliasedDPtr (size_t _size)
{
    dptr        = new SimpleDPtr<T>(_size);
    currentsize = _size;
    offset      = 0;

    if (dptr.isNull())
        BoxLib::OutOfMemory(__FILE__, __LINE__);
}

template <class T>
AliasedDPtr<T>::AliasedDPtr (AliasedDPtr<T>* _asdptr,
                             size_t          _size,
                             int             _offset)
{
    dptr        = _asdptr->dptr;
    offset      = _offset;
    currentsize = _size;
}

template <class T>
AliasedDPtr<T>::~AliasedDPtr ()
{}

template <class T>
size_t
AliasedDPtr<T>::size () const
{
    return currentsize;
}

template <class T>
T&
AliasedDPtr<T>::operator[] (long n) const
{
    return (*dptr)[n + offset];
}

template <class T>
void
AliasedDPtr<T>::define (size_t _size)
{
    dptr        = new SimpleDPtr<T>(_size);
    currentsize = _size;
    offset      = 0;

    if (dptr.isNull())
        BoxLib::OutOfMemory(__FILE__, __LINE__);
}

template <class T>
void
AliasedDPtr<T>::resize (size_t _size)
{
    dptr        = new SimpleDPtr<T>(_size);
    currentsize = _size;
    offset      = 0;

    if (dptr.isNull())
        BoxLib::OutOfMemory(__FILE__, __LINE__);
}

template <class T>
void
AliasedDPtr<T>::clear ()
{
    dptr = 0;
}

#endif /*_ALIASEDDPTR_H_*/
