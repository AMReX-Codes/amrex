//BL_COPYRIGHT_NOTICE

#ifndef _DPTR_H_
#define _DPTR_H_

//
// $Id: DPtr.H,v 1.1 1997-09-12 18:00:06 lijewski Exp $
//

#include <stddef.h>

#include <Assert.H>
#include <Arena.H>
#include <Boolean.H>

template <class T> class DPtrRep;

//
//@Man:
//@Memo: A Smart Pointer to Type T
/*@Doc:

  This is a smart pointer class.  It contains a private pointer to a
  DPtrRep<T> object, which is a virtual base class for classes such as
  SimpleDPtr<T> and AliasedDPtr<T> that implement the internal details of
  smart pointer memory management, using various forms of the Arena class.
  DPtr<T> is meant to look much like a T*, but be made to work with
  various types of memory-management schemes.

  This class does NOT provide a copy constructor or assignment operator.
*/

template <class T>
class DPtr
{
public:
    //
    //@ManDoc: Construct a null DPtr.
    //
    inline DPtr ();
    //
    //@ManDoc: Construct a DPtr pointing to ptr.
    //
    inline DPtr (DPtrRep<T>* ptr);
    //
    //@ManDoc: Delete to what this object points.
    //
    inline ~DPtr ();
    //
    //@ManDoc: Returns the internal pointer.
    //
    inline DPtrRep<T>* getDPtr ();
 
    /*@ManDoc: Set the internal pointer.  This will abort if the internal
               pointer is already set to something.
    */
    inline void setDPtr (DPtrRep<T>* ptr);
    //
    //@ManDoc: Returns true if the internal pointer is not null.
    //
    inline bool ok () const;
    //
    //@ManDoc: Sets the internal pointer to null.
    //
    inline void clearDPtr ();
    //
    //@ManDoc: Smart pointer access.  Abort if the pointer is null.
    //
    inline DPtrRep<T>* operator-> ();

    /*@ManDoc: Smart array element access. This assumes that the Dptr<T> was
               initialized with a pointer to an array of T.
               Abort if the pointer is null.
    */
    inline T& operator[] (long n) const;
    //
    //@ManDoc: Smart pointer dereference.  Abort if the pointer is null.
    //
    inline T& operator* () const;
    //
    //@ManDoc: Smart type conversion.  Abort if the pointer is null.
    //
    inline operator T* ();
    //
    //@ManDoc: Smart const type conversion.  Abort if the pointer is null.
    //
    inline operator const T* () const;

private:
    //
    // The underlying representation.
    //
    DPtrRep<T>* dptr;

private:
    //
    // These are disallowed.
    //
    DPtr (const DPtr<T>&);
    DPtr<T>& operator= (const DPtr<T>&);
};

//@Man:
//@Memo: Smart Pointer Virtual Base Class
/*@Doc:

  This is a virtual base class for classes implementing memory management
  for DPtr; i.e. its only use is as a base class for implementing memory
  management for smart pointers.

  This class does not provide a copy constructor of assignment operator.
*/

template <class T>
class DPtrRep
{
public:
    //
    //@ManDoc: Construct a DPtrRep using the given Arena.
    //
    DPtrRep (Arena *_arena = 0)
        : arena(_arena)
    {}
    //
    //@ManDoc: Virtual destructor.
    //
    virtual ~DPtrRep ();
    //
    //@ManDoc: Pure virtual function intended to clear the memory arena.
    //
    virtual void clear () = 0;
 
    /*@ManDoc: Pure virtual function intended to initialize the memory
               arena to hold size objects of type T.
    */
    virtual void define (size_t size) = 0;
 
    /*@ManDoc: Pure virtual function intended to resize the memory arena to
               hold at least newsize objects of type T.
    */
    virtual void resize (size_t newsize) = 0;
    //
    //@ManDoc: Pure virtual function intended to return the size of the arena.
    //
    virtual size_t size () const = 0;
    //
    //@ManDoc: Pure virtual smart array access.
    //
    virtual T& operator[] (long n) const = 0;
    //
    //@ManDoc: Virtual function returning address of first T in the arena.
    //
    virtual T* dataPtr ();
    //
    //@ManDoc: Set the arena.
    //
    void setArena (Arena *_arena)
    {
        arena = _arena;
    }
    //
    //@ManDoc: Get the arena.
    //
    Arena *getArena () const
    {
        return arena;
    }
    //
    //@ManDoc: No-op virtual function for future expansion.
    //
    virtual void Lock ();
    //
    //@ManDoc:  No-op virtual function for future expansion.
    //
    virtual void UnLock ();
    //
    //@ManDoc:  No-op virtual function for future expansion.
    //
    virtual void Demand ();
    //
    //@ManDoc:  No-op virtual function for future expansion.
    //
    virtual void Queue ();

protected:
    //
    // The arena.
    //
    Arena* arena;

private:
    //
    // These are disallowed.
    //
    DPtrRep (const DPtrRep<T>&);
    DPtrRep<T>& operator= (const DPtrRep<T>&);
};

//
// Inlines.
//

template <class T>
inline
DPtr<T>::DPtr ()
    : dptr(0)
{}

template <class T>
inline
DPtr<T>::DPtr (DPtrRep<T>* _dptr)
    : dptr(_dptr)
{}

template <class T>
inline
DPtr<T>::~DPtr ()
{
    delete dptr;
    dptr = 0;
}

template <class T>
inline
DPtrRep<T>*
DPtr<T>::getDPtr ()
{
    return dptr;
}

template <class T>
inline
void
DPtr<T>::setDPtr(DPtrRep<T>* _dptr)
{
    assert(dptr == 0);
    dptr = _dptr;
}

template <class T>
inline
bool
DPtr<T>::ok () const
{
    return dptr != 0;
}

template <class T>
inline
void
DPtr<T>::clearDPtr ()
{
    dptr = 0;
}

template <class T>
inline
DPtrRep<T>*
DPtr<T>::operator-> ()
{
    assert(dptr != 0);
    return dptr;
}

template <class T>
inline
T&
DPtr<T>::operator[] (long n) const
{
    assert(dptr != 0);
    return (*dptr)[n];
}

template <class T>
inline
T&
DPtr<T>::operator* () const
{
    assert(dptr != 0);
    return (*dptr)[0];
}

template <class T>
inline
DPtr<T>::operator T* ()
{
    assert(dptr != 0);
    return dptr->dataPtr();
}

template <class T>
inline
DPtr<T>::operator const T* () const
{
    assert(dptr != 0);
    return dptr->dataPtr();
}


template <class T>
DPtrRep<T>::~DPtrRep ()
{}

template <class T>
T*
DPtrRep<T>::dataPtr ()
{
    return &(*this)[0];
}

template <class T>
void
DPtrRep<T>::Lock ()
{}

template <class T>
void
DPtrRep<T>::UnLock ()
{}

template <class T>
void
DPtrRep<T>::Demand ()
{}

template <class T>
void
DPtrRep<T>::Queue ()
{}

#endif /*_DPTR_H_*/
