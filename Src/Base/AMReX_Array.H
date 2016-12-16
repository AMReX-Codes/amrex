
#ifndef BL_ARRAY_H
#define BL_ARRAY_H

#include <algorithm>
#include <vector>
#include <memory>
#include <AMReX_BLassert.H>

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

namespace amrex
{
    template <class T>
    Array<T*> GetArrOfPtrs (const Array<std::unique_ptr<T> >& a)
    {
	Array<T*> r;
	for (const auto& x : a)  r.push_back(x.get());
	return r;
    }

    template <class T>
    Array<const T*> GetArrOfConstPtrs (const Array<std::unique_ptr<T> >& a)
    {
	Array<const T*> r;
	for (const auto& x : a)  r.push_back(x.get());
	return r;
    }

    template <class T>
    Array<const T*> GetArrOfConstPtrs (const Array<T*>& a)
    {
	return {a.begin(), a.end()};
    }

    template <class T>
    Array<Array<T*> > GetArrOfArrOfPtrs (const Array<Array<std::unique_ptr<T> > >& a)
    {
	Array<Array<T*> > r;
	for (const auto& x : a)  r.push_back(GetArrOfPtrs(x));
	return r;
    }

    template <class T>
    void FillNull (Array<T*>& a)
    {
	std::for_each(a.begin(), a.end(), [](T*& p) { p = nullptr; });
    }

    template <class T>
    void FillNull (Array<std::unique_ptr<T> >& a)
    {
	std::for_each(a.begin(), a.end(), [](std::unique_ptr<T>& p) { p.reset(); });
    }
}

#endif /*BL_ARRAY_H*/
