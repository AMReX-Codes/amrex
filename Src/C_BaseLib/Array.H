#ifndef BL_ARRAY_H
#define BL_ARRAY_H
//
// $Id: Array.H,v 1.31 2004-03-08 23:06:36 car Exp $
//
#include <cstddef>
#include <vector>

#include <BLassert.H>
#include <BoxLib.H>

template <class T> class Array;

//
//@Man:
//@Memo: An Array of Objects of Type T
/*@Doc: 

   This class implements an array of objects of the parameterized type
   T.  In contrast with the predefined C++ array type, an `Array<T>'
   object knows its size, can be dynamically resized, and provides
   automatic bounds checking.  The bounds checking can be turned off by
   specifying the -DNDEBUG flag on the command line when compiling the
   BOXLIB library.  The main reason for defining the ARRAY class is that
   it is used, either by composition or inheritance, to implement many of
   the other classes in the BOXLIB library.

   The `Array<T>' class works by storing copies of the objects it
   contains.  If the objects are large, such as `FARRAYBOX's it is
   probably better to use the `PArray' class which is an array class that
   stores pointers to the objects (avoiding expensive copies).
   The `Array<T>' class destructs all objects in the array when it is 
   itself destructed.  If this is not the desired action, the `PArray' class
   should be used.

   In the Array<T> class, there are two different concepts of size: the amount
   of space allocated, and the amount of space actually used.  Obviously, the 
   allocated space must be larger than the used space.  We separate these two
   concepts to allow the user to optionally avoid memory allocation costs.
   Rather than construct and destroy a temporary Array<T> many times, it may 
   be less expensive to allocate it once with enough memory to handle all uses,
   and resize the Array<T> for each particular use.  See the member functions
   `reserve', `shrinkWrap', and `resize'. 

   Note that care must be taken when deriving classes from `Array<T>'.  It is
   a concrete class, not a polymorphic one.

   This class does NOT provide an assignment operator for assigning an integer
   to an Array<T>.
*/

template <class T>
class Array
    :
    public std::vector<T>
{
public:
    //
    //@ManDoc: Constructs an `Array<T>' with no elements
    //
    Array ();

    /*@ManDoc: Constructs an `Array<T>' of size len with the value of each
      element defined by the default constructor for `T'.
    */
    explicit Array (size_t len);

    /*@ManDoc: Constructs an `Array<T>' of size len with the value of each
               elements given by initialvalue.
    */
    Array (size_t   len,
           const T& initialvalue);

    /*@ManDoc: Constructs an `Array<T>' of size len in which the K'th
               value is a copy of vec[K].
    */
    Array (const T* vec,
           size_t   len);

    int size () const;
  
    /*@ManDoc: Returns a reference to the K'th element in this `Array<T>'.
               The element can be modified through this reference.  The
               result can be used as an L-value.
    */
    T& operator[] (size_t i);
    //
    //@ManDoc: Same as above, except acts on const Array's.
    //
    const T& operator[] (size_t i) const;

    //
    //@ManDoc: Different syntax for operator[] (long i).
    //
    T& get (size_t i);
    //
    //@ManDoc: Different syntax for const operator[] (long i).
    //
    const T& get (size_t i) const;
 
    /*@ManDoc: Returns pointer to vector of data.  This function breaks object
               encapsulation and should only be used for interfacing to
               Fortran subroutines.
    */
    T* dataPtr ();
    //
    //@ManDoc: Same as above for constant arrays.
    //
    const T* dataPtr () const;
    //
    //@ManDoc: Changes the i'th element of this `Array<T>' to elem.
    //
    void set (size_t   i,
              const T& elem);
private:
    //
    // This is disallowed.
    //
    Array<T>& operator= (int);
};

template <class T>
inline
int
Array<T>::size () const
{
    const std::vector<T>& a = *this;
    return a.size();
}

template <class T>
inline
T&
Array<T>::operator[] (size_t i)
{
    BL_ASSERT(i < size());
    std::vector<T>& a = *this;
    return a[i];
}

template <class T>
inline
const T&
Array<T>::operator[] (size_t i) const
{
    BL_ASSERT(i < size());
    const std::vector<T>& a = *this;
    return a[i];
}

template <class T>
inline
T&
Array<T>::get (size_t i)
{
    BL_ASSERT(i < size());
    return this->operator[](i);
}

template <class T>
inline
const T&
Array<T>::get (size_t i) const
{
    BL_ASSERT(i < size());
    return std::vector<T>::operator[](i);
}

template <class T>
inline
void
Array<T>::set (size_t     i,
               const T& elem)
{
    get(i) = elem;
}

template <class T>
inline
T*
Array<T>::dataPtr ()
{
    return &this->operator[](0);
}

template <class T>
inline
const T*
Array<T>::dataPtr () const
{
    return &this->operator[](0);
}

template <class T>
Array<T>::Array ()
    :
    std::vector<T>()
{}

template <class T>
Array<T>::Array (size_t len)
    :
    std::vector<T>(len)
{}

template <class T>
Array<T>::Array (size_t     len,
                 const T& initialValue)
    :
    std::vector<T>(len, initialValue)
{}

template <class T>
Array<T>::Array (const T* vec,
                 size_t     len)
    :
    std::vector<T>(len)
{
    BL_ASSERT(len == 0 || vec != 0);

    for (size_t i = 0; i < len; ++i)
    {
        this->operator[](i) = vec[i];
    }
}

#endif /*BL_ARRAY_H*/
