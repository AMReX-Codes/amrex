
#ifndef BL_TUPLE_H
#define BL_TUPLE_H

#include <cstdlib>

#include <BLassert.H>

//
// Ordered Tuples for Types T
//
// This class represents ordered tuples of some user-specified concrete
// type T for N > 0. The type T must have a default constructor.  If the
// non-default constructor, copy constructor, or copy assignment operator
// are used, T must also have a copy constructor.
//

template <class T, std::size_t N>
class Tuple
{
public:
    //
    // The default constructor.  For user-defined types T, the
    // default constructor for T will be run on each of the N
    // objects in the Tuple.  For builtin (intrinsic) types,
    // the values in the Tuple will be garbage.
    //
    Tuple () {}
    //
    // Constructs a Tuple, initializing the elements in the Tuple
    // with the corresponding elements in the vector v.  This assumes
    // that v contains at least N elements of type T -- an assumption
    // that is NOT checked.  For user-defined types, T must have a
    // well-defined and accessible copy constructor.
    //
    explicit Tuple (const T* v);
    //
    // The copy constructor.
    //
    Tuple (const Tuple& rhs);
    //
    // The copy assignment operator.
    //
    Tuple& operator= (const Tuple& rhs);
    //
    // Returns a reference to the i'th element in the Tuple,
    // counting from zero.  Performs range checking when the library
    // is compiled in debug mode.
    //
    T& operator[] (int i) { BL_ASSERT(0 <= i && i < N); return vect[i]; }
    //
    // Returns a constant reference to the i'th element in the Tuple,
    // counting from zero.  Performs range checking when the library
    // is compiled in debug mode.
    //
    const T& operator[] (int i) const { BL_ASSERT(0 <= i && i < N); return vect[i]; }
    //
    // Returns the address of the underlying vector of T
    // representation.  This should ONLY be used when interfacing
    // to Fortran as it breaks the encapsulation of the class.
    //
    operator const T* () const { return &vect[0]; }

protected:
    //
    // The underlying vector of T representing the Tuple.
    //
    T vect[N];
};

template <class T, std::size_t N>
Tuple<T,N>::Tuple (const T* v)
{
    BL_ASSERT(v != 0);
    for (std::size_t i = 0; i < N; ++i)
        vect[i] = v[i];
}

template <class T, std::size_t N>
Tuple<T,N>::Tuple (const Tuple<T,N>& rhs)
{
    for (std::size_t i = 0; i < N; ++i)
        vect[i] = rhs.vect[i];
}

template <class T, std::size_t N>
Tuple<T,N>&
Tuple<T,N>::operator= (const Tuple<T,N>& rhs)
{
    if (this != &rhs)
    {
        for (std::size_t i = 0; i < N; ++i)
            vect[i] = rhs.vect[i];
    }
    return *this;
}

#endif /*BL_TUPLE_H*/
