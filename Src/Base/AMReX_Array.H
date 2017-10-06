
#ifndef BL_ARRAY_H
#define BL_ARRAY_H

#include <algorithm>
#include <vector>
#include <memory>
#include <AMReX_BLassert.H>
#include <AMReX_Vector.H>

namespace amrex {

    template <class T> using Array = Vector<T>;
}

namespace amrex
{
    template <class T>
    Vector<T*> GetArrOfPtrs (const Vector<std::unique_ptr<T> >& a)
    {
	Vector<T*> r;
	for (const auto& x : a)  r.push_back(x.get());
	return r;
    }

    template <class T>
    Vector<const T*> GetArrOfConstPtrs (const Vector<std::unique_ptr<T> >& a)
    {
	Vector<const T*> r;
	for (const auto& x : a)  r.push_back(x.get());
	return r;
    }

    template <class T>
    Vector<const T*> GetArrOfConstPtrs (const Vector<T*>& a)
    {
	return {a.begin(), a.end()};
    }

    template <class T>
    Vector<Vector<T*> > GetArrOfArrOfPtrs (const Vector<Vector<std::unique_ptr<T> > >& a)
    {
	Vector<Vector<T*> > r;
	for (const auto& x : a)  r.push_back(GetArrOfPtrs(x));
	return r;
    }
}

#endif /*BL_ARRAY_H*/
