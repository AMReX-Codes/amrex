//BL_COPYRIGHT_NOTICE

#ifndef BL_POINTERS_H
#define BL_POINTERS_H

//
// $Id: Pointers.H,v 1.8 1999-05-10 17:18:47 car Exp $
//

#include <BLassert.H>
#include <BoxLib.H>
#include <UseCount.H>

//
//@Man:
//@Memo: A Smart Pointer for Intrinsic or User-Defined Types
/*@Doc:

  The template class CpPtr<T> provides a simple wrapper around a pointer
  to type T (T*) that builds a copy of the pointed-to object when copied
  from one CpPtr<T> to another.  This is in contrast to a reference-counted
  pointer class that would maintain one pointed-to object with a reference
  count indicating the number of references.  Hence we call this a
  "copied" smart pointer class.  It is intended for use with any type
  type T, including the intrinsic types.  This class does not supply
  an operator->(), as such an operator on intrinsic types has only recently
  become a part of the C++ language, and many compilers do not yet implement
  it.
*/

template <class T>
class CpPtr
{
public:
    //
    //@ManDoc: The default constructor.  The wrapped pointer is null.
    //
    CpPtr ();
    //
    //@ManDoc: Construct a CpPtr<T> setting the wrapped pointer to rhs.
    //
    explicit CpPtr (T* rhs);
    //
    //@ManDoc: The destructor.  Deletes the wrapped pointer.
    //
    ~CpPtr ();

    /*@ManDoc: The copy constructor.  If the pointer wrapped by rhs is null,
               the wrapped pointer is null here as well.  Otherwise,
               the contained pointer here is set to a new'd copy of that
               wrapped by rhs, with the two pointed-to values being identical.
               This assumes that type T has a well-defined and accessible copy
               constructor.  T must also be a concrete type, not a abstract
               type.
    */
    CpPtr (const CpPtr<T>& rhs);

    /*@ManDoc: Sets the wrapped pointer to rhs.  Deletes the previously
               wrapped pointer.
    */
    CpPtr<T>& operator= (T* rhs);

    /*@ManDoc: The copy assignment operator.  If the pointer wrapped by rhs
               is null, the wrapped pointer is null here as well.  Otherwise,
               the contained pointer here is set to a new'd copy of that
               wrapped by rhs, with the two pointed-to values being identical.
               This assumes that type T has a well-defined and accessible copy
               constructor.  T must also be a concrete type, not a abstract
               type.
    */
    CpPtr<T>& operator= (const CpPtr<T>& rhs);

    /*@ManDoc: Returns a reference to the value pointed to by the wrapped
               pointer; i.e. dereferencing this CpPtr<T>, returns the
               dereferenced wrapped pointer.  It is an error if the wrapped
               pointer is null.
    */
    T& operator* () const;
    //
    //@ManDoc: Returns true if the wrapped pointer null.
    //
    bool isNull () const;
    //
    //@ManDoc: Sets the wrapped pointer to null and returns the previous value.
    //
    T* release ();
    //
    //@ManDoc: Are the two pointers (not the values to which they point) equal?
    //
    bool operator== (const CpPtr<T>& rhs) const;
    //
    //@ManDoc: Are the two pointers not equal?
    //
    bool operator!= (const CpPtr<T>& rhs) const;

protected:
    T* ptr;
};

//
//@Man:
//@Memo: A Smart Pointer for User-Defined Types
/*@Doc:

  The template class CpClassPtr<T> is derived from CpPtr<T>.  It provides a
  simple wrapper around a pointer to type T (a T*) that "does the right thing"
  when copied from one CpPtr<T> to another.  The type T MUST be a user-defined
  type, not an intrinsic type.  Given this restriction, we can supply
  an operator->().
*/

template<class T>
class CpClassPtr :
    public CpPtr<T>
{
public:
    //
    //@ManDoc: The default constructor.  The wrapped pointer is null.
    //
    CpClassPtr ();
    //
    //@ManDoc: Construct a CpPtr<T> setting the wrapped pointer to rhs.
    //
    explicit CpClassPtr (T* rhs);

    /*@ManDoc: The copy constructor.  If the pointer wrapped by rhs is null,
               the wrapped pointer is null here as well.  Otherwise,
               the contained pointer here is set to a new'd copy of that
               wrapped by rhs, with the two pointed-to values being identical.
               This assumes that type T has a well-defined and accessible copy
               constructor.  T must also be a concrete type, not a abstract
               type.
    */
    CpClassPtr (const CpClassPtr<T>& rhs);

    /*@ManDoc: Sets the wrapped pointer to rhs.  Deletes the previously
               wrapped pointer.
    */
    CpClassPtr<T>& operator= (T* rhs);

    /*@ManDoc: The copy assignment operator.  If the pointer wrapped by rhs
               is null, the wrapped pointer is null here as well.  Otherwise,
               the contained pointer here is set to a new'd copy of that
               wrapped by rhs, with the two pointed-to values being identical.
               This assumes that type T has a well-defined and accessible copy
               constructor.  T must also be a concrete type, not a abstract
               type.
    */
    CpClassPtr<T>& operator= (const CpClassPtr<T>& rhs);
    //
    //@ManDoc: Applies operator-> to the wrapped pointer.
    //
    T* operator-> () const;
};

//
//@Man:
//@Memo: A Reference Counted Smart Pointer for Intrinsic or User-Defined Types
/*@Doc:

  The template class LnPtr<T> provides a reference counted wrapper around a
  pointer to type T (a T*).  This "smart" pointer is intended for use with
  any type type T, including the intrinsic types.  For this reason, we do
  not supply an operator->(), as such an operator on intrinsic types has only
  recently become a part of the C++ language and many compilers do not yet
  implement it.
*/

template<class T>
class LnPtr
{
public:
    //
    //@ManDoc: The default constructor.  The wrapped pointer is null.
    //
    LnPtr ();
    //
    //@ManDoc: Construct a LnPtr<T> setting the wrapped pointer to rhs.
    //
    explicit LnPtr (T* rhs);

    /*@ManDoc: The copy assignment operator.  The contained pointer is set
               to the one wrapped by rhs.  The reference count is decremented
               on this object and the reference count is incremented for
               the newly wrapped pointer.
    */
    LnPtr<T>& operator= (const LnPtr<T>& rhs);

    /*@ManDoc: Sets the wrapped pointer to rhs.  Decrements the count
               on the previously wrapped pointer and deletes it if there
               was only one reference.
    */
    LnPtr<T>& operator= (T* rhs)
    {
        //
        // This is inlined here as OCC won't have it any other way :-(
        //
        if (unique())
            delete ptr;
        ptr = rhs;
        ucnt = UseCount();
        return *this;
    }

    /*@ManDoc: The destructor -- decrements the reference count and deletes
               the wrapped pointer if there is only one reference.
    */
    ~LnPtr ();
    //
    //@ManDoc: Returns true if only one reference to the wrapped pointer.
    //
    bool unique () const;
    //
    //@ManDoc: Returns the number of references to the wrapped pointer.
    //
    int linkCount () const;

    /*@ManDoc: Returns a reference to the value pointed to by the wrapped
               pointer; i.e. dereferencing this LnPtr<T>, returns the
               dereferenced wrapped pointer.  It is an error if the wrapped
               pointer is null.
    */
    T& operator* () const;
    //
    //@ManDoc: Returns true if the wrapped pointer is null.
    //
    bool isNull () const;
    //
    //@ManDoc: Are the two pointers (not the values to which they point) equal?
    //
    bool operator== (const LnPtr<T>& rhs) const;
    //
    //@ManDoc: Are the two pointers not equal?
    //
     bool operator!= (const LnPtr<T>& rhs) const;

protected:
    T*       ptr;

private:
    UseCount ucnt;
};

//
//@Man:
//@Memo: A Smart Reference Counted Pointer for User-Defined Types
/*@Doc:

  The template class LnClassPtr<T> is derived from LnPtr<T>.  It provides a
  reference counted wrapper around a pointer to type T (a T*).  The type T
  MUST be a user-defined type, not an intrinsic type.  Given this
  restriction, we can supply an operator->().
*/

template<class T>
class LnClassPtr
    : public LnPtr<T>
{
public:
    //
    //@ManDoc: The default constructor.  The wrapped pointer is null.
    //
    LnClassPtr ();
    //
    //@ManDoc: Construct a LnPtr<T> setting the wrapped pointer to rhs.
    //
    explicit LnClassPtr (T* rhs);

    /*@ManDoc: The copy assignment operator.  The contained pointer is set
               to the one wrapped by rhs.  The reference count is decremented
               on this object and the reference count is incremented for
               the newly wrapped pointer.
    */
    LnClassPtr<T>& operator= (const LnClassPtr<T>& rhs);

    /*@ManDoc: Sets the wrapped pointer to rhs.  Decrements the count
               on the previously wrapped pointer and deletes it if there
               was only one reference.
    */
    LnClassPtr<T>& operator= (T* rhs);
    //
    //@ManDoc: Applies operator-> to the wrapped pointer.
    //
    T* operator->() const;
};

//
// Inlines.
//

template <class T>
inline
CpPtr<T>::CpPtr ()
    : ptr(0)
{}

template <class T>
inline
CpPtr<T>::CpPtr(T* rhs)
    : ptr(rhs)
{}

template <class T>
inline
CpPtr<T>::~CpPtr()
{
    delete ptr;
}

template <class T>
inline
bool
CpPtr<T>::isNull () const
{
    return ptr == 0;
}

template <class T>
inline
CpPtr<T>::CpPtr (const CpPtr<T>& rhs)
{
    ptr = rhs.isNull() ?  0 : new T(*rhs.ptr);
}

template <class T>
inline
CpPtr<T>&
CpPtr<T>::operator= (const CpPtr<T>& rhs)
{
    if (!(ptr == rhs.ptr))
    {
        delete ptr;
        ptr = rhs.isNull() ? 0 : new T(*rhs.ptr);
    }
    return *this;
}

template <class T>
inline
CpPtr<T>&
CpPtr<T>::operator= (T* rhs)
{
    delete ptr;
    ptr = rhs;
    return *this;
}

template <class T>
inline
T&
CpPtr<T>::operator* () const
{
    BLassert(ptr != 0);
    return *ptr;
}

template <class T>
inline
T*
CpPtr<T>::release ()
{
    T* old = ptr;
    ptr = 0;
    return old;
}

template <class T>
inline
bool
CpPtr<T>::operator== (const CpPtr<T>& rhs) const
{
    return ptr == rhs.ptr;
}

template <class T>
inline
bool
CpPtr<T>::operator!= (const CpPtr<T>& rhs) const
{
    return ptr != rhs.ptr;
}

template <class T>
inline
CpClassPtr<T>::CpClassPtr ()
    : CpPtr<T>()
{}

template <class T>
inline
CpClassPtr<T>::CpClassPtr (T* rhs)
    : CpPtr<T>(rhs)
{}

template <class T>
inline
CpClassPtr<T>::CpClassPtr (const CpClassPtr<T>& rhs)
    : CpPtr<T>(rhs) {}

template <class T>
inline
CpClassPtr<T>&
CpClassPtr<T>::operator= (T* rhs)
{
    CpPtr<T>::operator= (rhs);
    return *this;
}

template <class T>
inline
CpClassPtr<T>&
CpClassPtr<T>::operator= (const CpClassPtr<T>& rhs)
{
    CpPtr<T>::operator= (rhs);
    return *this;
}

template <class T>
inline
T*
CpClassPtr<T>::operator-> () const
{
    return ptr;
}

template <class T>
inline
LnPtr<T>::LnPtr ()
    : ptr(0)
{}

template <class T>
inline
LnPtr<T>::LnPtr(T* rhs)
    : ptr(rhs)
{}

template <class T>
inline
bool
LnPtr<T>::unique () const
{ 
    return ucnt.unique();
}

template <class T>
inline
LnPtr<T>::~LnPtr ()
{ 
    if (ucnt.unique())
        delete ptr;
}

template <class T>
inline
LnPtr<T>&
LnPtr<T>::operator= (const LnPtr<T>& rhs)
{
    if (ptr != rhs.ptr)
    {
        if (unique())
            delete ptr;
        ptr = rhs.ptr;
        ucnt = rhs.ucnt;
    }
    return *this;
}

template <class T>
inline
int
LnPtr<T>::linkCount () const
{ 
    return ucnt.linkCount();
}

template <class T>
inline
T&
LnPtr<T>::operator* () const
{
    BLassert(ptr != 0);
    return *ptr;
}

template <class T>
inline
bool
LnPtr<T>::isNull () const
{
    return ptr == 0;
}

template <class T>
inline
bool
LnPtr<T>::operator== (const LnPtr<T>& rhs) const
{
    return ptr == rhs.ptr;
}

template <class T>
inline
bool
LnPtr<T>::operator!= (const LnPtr<T>& rhs) const
{
    return ptr != rhs.ptr;
}

template <class T>
inline
LnClassPtr<T>::LnClassPtr ()
{}

template <class T>
inline
LnClassPtr<T>::LnClassPtr (T* rhs)
    : LnPtr<T>(rhs)
{}

template <class T>
inline
LnClassPtr<T>&
LnClassPtr<T>::operator= (const LnClassPtr<T>& rhs)
{
    LnPtr<T>::operator=(rhs);
    return *this;
}

template <class T>
inline
LnClassPtr<T>&
LnClassPtr<T>::operator= (T* rhs)
{
    LnPtr<T>::operator=(rhs);
    return *this;
}

template <class T>
inline
T*
LnClassPtr<T>::operator->() const
{
    return ptr;
}

#endif /*BL_POINTERS_H*/
