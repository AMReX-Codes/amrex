#ifndef BL_PARRAY_H
#define BL_PARRAY_H
//
// $Id: PArray.H,v 1.24 2002-11-12 20:49:35 lijewski Exp $
//
#include <Array.H>

//
//@Man:
//@Memo:
/*@Doc: An enum that controls whether or not class PArray<T>
         manages the pointed to objects, as well as the pointer itself.
         Valid values are PArrayNoManage and PArrayManage.
*/

enum PArrayPolicy
{
    PArrayNoManage = 0,
    PArrayManage   = 1
};

//
//@Man:
//@Memo: An Array of Pointers to Objects of Type T
/*@Doc:

  The PArray<T> class is similar to the Array<T> class except that it
  implements an array of pointers to objects of type T rather than an array
  of the objects themselves.  As with Array<T> objects, PArray<T> objects
  know about their length, can be dynamically resized and have built-in
  automatic bounds checking.

  When a PArray<T> object is destructed, should only the pointers be
  destructed or should the objects to which they point also be deleted?
  The PArray<T> class allows for both scenarios.  By default, it will
  only remove the pointers.  If the array is constructed with the
  PArrayManage argument it will also delete all the objects pointed to
  in the array.

  This is a concrete class, not a polymorphic one.

  PArray<T>s do NOT provide a copy constructor or an assignment operator.
*/

template <class T>
class PArray
{
public:

    /*@ManDoc: Constructs a PArray<T> with no elements.  The default memory
               management policy is not to delete the underlying objects
               upon destruction.
    */
    explicit PArray (PArrayPolicy _managed = PArrayNoManage);

    /*@ManDoc: Constructs a PArray<T> of length len of null pointers.  
               The default memory management policy is not to delete the
               underlying objects upon destruction.
    */
    explicit PArray (int          len,
                     PArrayPolicy _managed = PArrayNoManage);

    /*@ManDoc: The destructor.  If the memory management policy was
               set to Manage upon construction, this destructor will
               also delete all objects referenced by the pointers in
               this array.
    */
    ~PArray ();
    //
    //@ManDoc: Returns true if the i'th element is not the null pointer.
    //
    bool defined (int i) const;

    /*@ManDoc: Returns a constant reference to the object pointed to by the
               i'th element in this PArray<T>.  Performs bounds checking
               when the library is compiled in debug mode.
    */
    const T& operator[] (int i) const;

    /*@ManDoc: Returns a reference to the object pointed to by the
               i'th element in this PArray<T>.  Performs bounds checking
               when the library is compiled in debug mode.
    */
    T& operator[] (int i);

    /*@ManDoc: Returns a constant reference to the object pointed to by the
               i'th element in the PArray.  It is an error if the i'th element
               is the null pointer.
    */
    const T& get (int i) const;

    /*@ManDoc: Returns a reference to the object pointed to by the i'th
               element in the PArray.  It is an error if the i'th element
               is the null pointer.
    */
    T& get (int i);

    /*@ManDoc: Defines the i'th element of this PArray<T> to the pointer
               value elem.  It is an error if the i'th element previously
               had a value.
    */
    void set (int i,
              T*  elem);
    //
    //@ManDoc: Returns the number of elements.
    //
    int size () const;

    /*@ManDoc: This function sets each pointer in this PArray<T> to null.
               If the memory management policy was set to Manage
               upon construction, this function will also delete all objects
               referenced by the pointers in this array.
    */
    void clear ();

    /*@ManDoc: The i'th element in this PArray<T> is set to null.
               If the memory management policy was set to Manage
               upon construction, this function will also delete the object
               referenced by that pointer.
    */
    void clear (int i);

    /*@ManDoc: This function changes the size of the array to newsize
               while preserving the value of as many array elements as
               possible.  If newsize is greater than the current length,
               the array is grown and the new elements are null pointers.
               If newsize is less than the current length the array is
               cropped with the remaining elements retaining their
               original pointer values. If the memory management policy
               was set to Manage, any cropped pointers are deleted.
    */
    void resize (int newsize);

    /*@ManDoc: This function does a resize and changes the PArrayPolicy.
	       (Especially for arrays of PArray where the user wants
	       the default policy to be managed--the default constructor
	       makes the default unmanaged).  If the new managed policy is
               different from the old one, none of the elements of the
               original are retained.
    */
    void resize (int          newsize,
                 PArrayPolicy newmanagedpolicy);

    /*@ManDoc: Removes the i'th element from the PArray and returns the
               pointer to the calling function.  Nulls out the pointer
               in the PArray<T>.
    */
    T* remove (int i);

protected:
    //
    // The underlying representation.
    //
    Array<void*> vp;
    //
    // The memory management policy.
    //
    PArrayPolicy managed;

public:
    //
    // Disallow these. FIXME
    //
    PArray (const PArray<T>&);
    // PArray<T> operator= (const PArray<T>&);
};

template <class T>
PArray<T>::PArray (const PArray<T>& c)
{
    BL_ASSERT(c.size()==0);

    managed = c.managed;
}

template <class T>
PArray<T>::PArray (PArrayPolicy _managed)
    :
    managed(_managed)
{}

template <class T>
inline
bool
PArray<T>::defined (int i) const
{
    return vp[i] != 0;
}

template <class T>
inline
const T&
PArray<T>::operator[] (int i) const
{
    BL_ASSERT(vp[i] != 0);
    return *((T*)(vp[i]));
}

template <class T>
inline
T&
PArray<T>::operator[] (int i)
{
    BL_ASSERT(vp[i] != 0);
    return *((T*)(vp[i]));
}

template <class T>
inline
const T&
PArray<T>::get (int i) const
{
    BL_ASSERT(vp[i] != 0);
    return *((T*)(vp[i]));
}

template <class T>
inline
T&
PArray<T>::get (int i)
{
    BL_ASSERT(vp[i] != 0);
    return *((T*)(vp[i]));
}

template <class T>
inline
void
PArray<T>::set (int i,
                T*  elem)
{
    BL_ASSERT(vp[i] == 0);
    vp[i] = elem;
}

template <class T>
inline
int
PArray<T>::size () const
{
    return vp.size();
}

template <class T>
T*
PArray<T>::remove (int i)
{
    BL_ASSERT(vp[i] != 0);
    void* tmp = vp[i];
    vp[i]     = 0;
    return (T*) tmp;
}

template <class T>
void
PArray<T>::clear (int n)
{
    if (managed)
        delete ((T*)(vp[n]));
    vp[n] = 0;
}

template <class T>
PArray<T>::PArray (int          len,
                   PArrayPolicy _managed)
    :
    vp(len),
    managed(_managed)
{
    for (int i = 0; i < size(); i++)
        vp[i] = 0;
}

template <class T>
void
PArray<T>::clear ()
{
    if (managed)
    {
        for (int i = 0; i < size(); ++i)
        {
            delete ((T*)(vp[i]));
        }
    }
    for (int i = 0; i < size(); ++i)
        vp[i] = 0;
}

template <class T>
PArray<T>::~PArray ()
{
    clear();
}

template <class T>
void
PArray<T>::resize (int newlen)
{
    void** ovp = vp.size() ? &vp[0] : 0;

    int oldlen = vp.size();

    vp.resize(newlen);

    for (int i = oldlen; i < newlen; ++i)
        vp[i] = 0;

    if (managed)
    {
        for (int i = newlen; i < oldlen; i++)
        {
            delete ((T*)(ovp[i]));
            ovp[i]  = 0;
        }
    }
}

template <class T>
void
PArray<T>::resize (int          newlen,
                   PArrayPolicy newmanagedpolicy)
{
    if (newmanagedpolicy == managed)
    {
        resize(newlen);
    }
    else
    {
        clear();
        resize(newlen);
        managed = newmanagedpolicy;
    }
}

#endif /*BL_PARRAY_H*/
