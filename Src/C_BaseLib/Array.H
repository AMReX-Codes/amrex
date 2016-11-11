
#ifndef BL_ARRAY_H
#define BL_ARRAY_H

#include <algorithm>
#include <vector>
#include <memory>
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

namespace BoxLib
{
    template <class T>
    Array<T*> GetArrOfPtrs (const Array<std::unique_ptr<T> >& a)
    {
	Array<T*> r(a.size());
	std::transform(a.begin(), a.end(), r.begin(),
		       [](const std::unique_ptr<T>& p) { return p.get(); });
	return r;
    }

    template <class T>
    Array<const T*> GetArrOfConstPtrs (const Array<std::unique_ptr<T> >& a)
    {
	Array<const T*> r(a.size());
	std::transform(a.begin(), a.end(), r.begin(),
		       [](const std::unique_ptr<T>& p) { return p.get(); });
	return r;
    }

    template <class T>
    Array<const T*> GetArrOfConstPtrs (const Array<T*>& a)
    {
	return {a.begin(), a.end()};
    }
}

#endif /*BL_ARRAY_H*/
