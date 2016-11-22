
#ifndef BL_ARRAY_H
#define BL_ARRAY_H

#include <vector>
#include <BLassert.H>

//
// This class is a thin wrapper around std::vector.  Unlike vector,
// Array::operator[] provides bound checking when compiled with
// DEBUG=TRUE.
//

template <class T>
class Array
    :
    public std::vector<T>
{
public:

    using std::vector<T>::vector;

    T& operator[] (std::size_t i)
    {
	BL_ASSERT(i < this->size());
	return this->std::vector<T>::operator[](i);
    }

    const T& operator[] (std::size_t i) const
    {
	BL_ASSERT(i < this->size());
	return this->std::vector<T>::operator[](i);
    }


    T* dataPtr () { return this->data(); }
    const T* dataPtr () const { return this->data(); }
};

#endif /*BL_ARRAY_H*/
