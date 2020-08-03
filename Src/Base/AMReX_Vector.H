#ifndef AMREX_VECTOR_H_
#define AMREX_VECTOR_H_

#include <algorithm>
#include <vector>
#include <memory>
#include <AMReX_BLassert.H>
#include <AMReX_INT.H>
#ifdef AMREX_SPACEDIM
#include <AMReX_Array.H>
#include <AMReX_TypeTraits.H>
#endif

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
    using typename std::vector<T, Allocator>::size_type;

    T& operator[] (size_type i) noexcept
    {
        BL_ASSERT( i < (this->std::vector<T, Allocator>::size()) );
	return this->std::vector<T, Allocator>::operator[](i);
    }

    const T& operator[] (size_type i) const noexcept
    {
        BL_ASSERT( i < (this->std::vector<T, Allocator>::size()) );
	return this->std::vector<T, Allocator>::operator[](i);
    }

    //! get access to the underlying data pointer
    T* dataPtr () noexcept { return this->data(); }
    //! get access to the underlying data pointer
    const T* dataPtr () const noexcept { return this->data(); }

    Long size () const noexcept {return static_cast<Long>(std::vector<T, Allocator>::size());}

};
    
}

namespace amrex
{
    /////////////////////////////////////////////////////////////

    template <class T, typename = typename T::FABType>
    Vector<T*> GetVecOfPtrs (Vector<T>& a)
    {
        Vector<T*> r;
        r.reserve(a.size());
        for (auto& x : a)  r.push_back(&x);
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

    /////////////////////////////////////////////////////////////

    template <class T, typename = typename T::FABType>
    Vector<const T*> GetVecOfConstPtrs (const Vector<T>& a)
    {
        Vector<const T*> r;
        r.reserve(a.size());
        for (const auto& x : a)  r.push_back(&x);
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

    template <class T, typename = typename T::FABType>
    Vector<const T*> GetVecOfConstPtrs (const Vector<T*>& a)
    {
        return {a.begin(), a.end()};
    }

    /////////////////////////////////////////////////////////////

    template <class T>
    Vector<Vector<T*> > GetVecOfVecOfPtrs (const Vector<Vector<std::unique_ptr<T> > >& a)
    {
        Vector<Vector<T*> > r;
        r.reserve(a.size());
        for (const auto& x : a)  r.push_back(GetVecOfPtrs(x));
        return r;
    }

    /////////////////////////////////////////////////////////////

#ifdef AMREX_SPACEDIM
    template <class T>
    Vector<std::array<T*,AMREX_SPACEDIM> >
    GetVecOfArrOfPtrs (const Vector<std::array<std::unique_ptr<T>,AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (const auto& x : a)  r.push_back(GetArrOfPtrs(x));
        return r;
    }

    template <class T>
    Vector<std::array<T const*,AMREX_SPACEDIM> >
    GetVecOfArrOfPtrsConst (const Vector<std::array<std::unique_ptr<T>,AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T const*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (const auto& x : a)  r.push_back(GetArrOfConstPtrs(x));
        return r;
    }

    template <class T>
    Vector<std::array<T const*,AMREX_SPACEDIM> >
    GetVecOfArrOfConstPtrs (const Vector<std::array<std::unique_ptr<T>,AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T const*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (const auto& x : a)  r.push_back(GetArrOfConstPtrs(x));
        return r;
    }

    template <class T, typename std::enable_if<IsFabArray<T>::value ||
                                               IsBaseFab<T>::value,
                                               int>::type = 0 >
    Vector<std::array<T const*,AMREX_SPACEDIM> >
    GetVecOfArrOfConstPtrs (const Vector<std::array<T,AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T const*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (const auto& x : a)  r.push_back(GetArrOfConstPtrs(x));
        return r;
    }

    template <class T, typename std::enable_if<IsFabArray<T>::value ||
                                               IsBaseFab<T>::value,
                                               int>::type = 0 >
    Vector<std::array<T*, AMREX_SPACEDIM> >
    GetVecOfArrOfPtrs(Vector<std::array<T, AMREX_SPACEDIM> >& a) 
    {
        Vector<std::array<T*, AMREX_SPACEDIM> > r; 
        r.reserve(a.size()); 
        for (auto &x: a) r.push_back(GetArrOfPtrs(x)); 
        return r; 
    }
#endif

    /////////////////////////////////////////////////////////////

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

    /////////////////////////////////////////////////////////////

    template <class T> 
    void RemoveDuplicates (Vector<T>& vec) {
        std::sort(vec.begin(), vec.end());
        auto it = std::unique(vec.begin(), vec.end());
        vec.erase(it, vec.end());
    }

    namespace detail {
        template <class T, class H>
        std::size_t removeDupDoit (Vector<T>& vec, std::size_t start, std::size_t stop)
        {
            std::size_t N = stop-start;
            if (N < 2) return stop;

            T* const data = vec.data() + start;
            T const sentinel = data[0]; // duplicates will be set to sentinel and removed later
            H const hasher;
            for (std::size_t i = 1; i < N; ) {
                if (data[i] == sentinel) {
                    ++i;
                    continue;
                }

                std::size_t const hash = hasher(data[i]) % N;
                if (i == hash) { // data[i] in correct hash position
                    ++i;
                    continue;
                }

                if (data[i] == data[hash]) {
                    data[i] = sentinel; // because it's a duplicate
                    ++i;
                    continue;
                }

                if (data[hash] == sentinel) {
                    std::swap(data[hash], data[i]);
                    // after swap, new data[i] holds sentinel
                    //   newdata[hash] in correct hash poitiion
                    ++i;
                    continue;
                }

                std::size_t const hashhash = hasher(data[hash]) % N;
                if (hashhash != hash) { // data[hash] not in correct has poision, thus will yield it's poistion
                    std::swap(data[i], data[hash]);
                    // after swap, new data[hash] in correct hash position
                    //   new data[i] not sure
                    if (hash < i) {  // we have seen new data[i]
                        ++i;
                    } // else next iteration we will work on data[i]
                } else { // data[hash] in correct hash position, but data[i] is not because of hash collision
                    ++i;
                }
            }

            // Now there are three types for data[i]
            //   (1) sentinel
            //   (2) data[i] != sentinel and hash(data[i]) == i
            //   (3) data[i] != sentinel and hash(data[i]) != i because of hash collision
            // All type 2s are unique, and all sentinels except one can be removed.
            // We will move all type 2s to the beginning, all sentinels to the end.
            // This will leave all type 3s in the middle.  Then we will work on the middle
            // part plus one sentinel.

            std::size_t swapPos = 0;
            for (std::size_t i = 0; i < N; ++i) {
                // move type 2 to the beginning pointed to by swapPos
                if (data[i] != sentinel and i == hasher(data[i]) % N) {
                    std::swap(data[i], data[swapPos++]);
                }
            }

            // Now we have moved all type 2 elements to the beginning, [0,swapPos)

            std::size_t sentinelPos = N;
            for (std::size_t i = swapPos; i < sentinelPos; ) {
                // move type 1 to the end
                if(data[i] == sentinel) {
                    std::swap(data[i], data[--sentinelPos]);
                } else {
                    ++i;
                }
            }

            // recursively work on the middle part
            return detail::removeDupDoit<T,H>(vec, start+swapPos, start+sentinelPos+1);
        }
    }

    template <class T, class H>
    void RemoveDuplicates (Vector<T>& vec) {
        // https://stackoverflow.com/questions/1532819/algorithm-efficient-way-to-remove-duplicate-integers-from-an-array
        std::size_t pos = detail::removeDupDoit<T,H>(vec, 0, vec.size());
        vec.erase(vec.begin()+pos, vec.end());
    }
}

#endif

