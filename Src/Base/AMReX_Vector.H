#ifndef AMREX_VECTOR_H_
#define AMREX_VECTOR_H_

#include <algorithm>
#include <vector>
#include <memory>
#include <AMReX_BLassert.H>


namespace amrex {
/**
* \brief This class is a thin wrapper around std::vector.  Unlike vector,
* Vector::operator[] provides bound checking when compiled with
* DEBUG=TRUE.
*/

template <class T, class Allocator=std::allocator<T> >
class Vector
    :
        public std::vector<T, Allocator>
{
public:

    using std::vector<T, Allocator>::vector;

    T& operator[] (std::size_t i)
    {
	BL_ASSERT(static_cast<long>(i) < this->size());
	return this->std::vector<T, Allocator>::operator[](i);
    }

    const T& operator[] (std::size_t i) const
    {
	BL_ASSERT(static_cast<long>(i) < this->size());
	return this->std::vector<T, Allocator>::operator[](i);
    }

    //! get access to the underlying data pointer
    T* dataPtr () { return this->data(); }
    //! get access to the underlying data pointer
    const T* dataPtr () const { return this->data(); }

    long size() const {return static_cast<long>(std::vector<T, Allocator>::size());}

};

}

namespace amrex
{
    template <class T>
    Vector<const T*> GetVecOfConstPtrs (const Vector<T>& a)
    {
	Vector<const T*> r;
        r.reserve(a.size());
	for (const auto& x : a)  r.push_back(&x);
	return r;
    }

    template <class T>
    Vector<T*> GetVecOfPtrs (const Vector<std::unique_ptr<T> >& a)
    {
	Vector<T*> r;
        r.reserve(a.size());
	for (const auto& x : a)  r.push_back(x.get());
	return r;
    }

    template <class T>
    Vector<const T*> GetVecOfConstPtrs (const Vector<std::unique_ptr<T> >& a)
    {
	Vector<const T*> r;
        r.reserve(a.size());
	for (const auto& x : a)  r.push_back(x.get());
	return r;
    }

    template <class T>
    Vector<const T*> GetVecOfConstPtrs (const Vector<T*>& a)
    {
	return {a.begin(), a.end()};
    }

    template <class T>
    Vector<Vector<T*> > GetVecOfVecOfPtrs (const Vector<Vector<std::unique_ptr<T> > >& a)
    {
	Vector<Vector<T*> > r;
        r.reserve(a.size());
	for (const auto& x : a)  r.push_back(GetVecOfPtrs(x));
	return r;
    }

    template <class T>
    void FillNull (Vector<T*>& a)
    {
	std::for_each(a.begin(), a.end(), [](T*& p) { p = nullptr; });
    }

    template <class T>
    void FillNull (Vector<std::unique_ptr<T> >& a)
    {
	std::for_each(a.begin(), a.end(), [](std::unique_ptr<T>& p) { p.reset(); });
    }

    template <class T> 
    void RemoveDuplicates(Vector<T>& vec) {
        std::sort(vec.begin(), vec.end());
        auto it = std::unique(vec.begin(), vec.end());
        vec.erase(it, vec.end());
    }
}

#endif

