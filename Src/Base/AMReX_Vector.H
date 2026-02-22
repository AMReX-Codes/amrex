#ifndef AMREX_VECTOR_H_
#define AMREX_VECTOR_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Vector.H
 * \brief A thin wrapper around std::vector with bounds checking and AMReX utilities.
 */

#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_INT.H>
#ifdef AMREX_SPACEDIM
#include <AMReX_Array.H>
#include <AMReX_TypeTraits.H>
#endif

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace amrex {
/**
 * \class amrex::Vector
 * \ingroup amrex_linear_containers
 * \brief A thin wrapper around \c std::vector with bounds checking.
 *
 * Unlike \c std::vector, `Vector::operator[]` performs bounds checking
 * when compiled with \c DEBUG=TRUE.
 *
 * \tparam T         Element type.
 * \tparam Allocator Allocator type; defaults to \c std::allocator<T>.
 */
template <class T, class Allocator=std::allocator<T> >
class Vector
    :
        public std::vector<T, Allocator>
{
public:

    using std::vector<T, Allocator>::vector;
    using typename std::vector<T, Allocator>::size_type;

    //! Return a reference to element \p i, with bounds checking in debug builds.
    [[nodiscard]] T& operator[] (size_type i) noexcept
    {
        AMREX_ASSERT_WITH_MESSAGE(i < (this->std::vector<T, Allocator>::size()),
                                  "Out of bound error, index: " + std::to_string(i) + " size: " + std::to_string(size()));
        return this->std::vector<T, Allocator>::operator[](i);
    }

    //! Return a const reference to element \p i, with bounds checking in debug builds.
    [[nodiscard]] const T& operator[] (size_type i) const noexcept
    {
        AMREX_ASSERT_WITH_MESSAGE(i < (this->std::vector<T, Allocator>::size()),
                                  "Out of bound error, index: " + std::to_string(i) + " size: " + std::to_string(size()));
        return this->std::vector<T, Allocator>::operator[](i);
    }

    //! get access to the underlying data pointer
    [[nodiscard]] T* dataPtr () noexcept { return this->data(); }
    //! get access to the underlying data pointer
    [[nodiscard]] const T* dataPtr () const noexcept { return this->data(); }

    //! Return the number of elements as a \ref amrex::Long (signed, wider than \c int).
    [[nodiscard]] Long size () const noexcept {return static_cast<Long>(std::vector<T, Allocator>::size());}

};

}

namespace amrex
{
    /////////////////////////////////////////////////////////////

    /// \ingroup amrex_linear_containers
    //! Return a Vector of raw pointers to the elements of \p a (FAB types only).
    template <class T, typename = typename T::FABType>
    [[nodiscard]] Vector<T*> GetVecOfPtrs (Vector<T>& a)
    {
        Vector<T*> r;
        r.reserve(a.size());
        for (auto& x : a) { r.push_back(&x); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of raw pointers to the Array elements of \p a (FAB types only).
    template <class T, std::size_t N, typename = typename T::FABType>
    [[nodiscard]] Vector<Array<T,N>*> GetVecOfPtrs (Vector<Array<T,N>>& a)
    {
        Vector<Array<T,N>*> r;
        r.reserve(a.size());
        for (auto& x : a) { r.push_back(&x); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of raw pointers from a Vector of unique_ptr.
    template <class T>
    [[nodiscard]] Vector<T*> GetVecOfPtrs (const Vector<std::unique_ptr<T> >& a)
    {
        Vector<T*> r;
        r.reserve(a.size());
        for (const auto& x : a) { r.push_back(x.get()); }
        return r;
    }

    /////////////////////////////////////////////////////////////

    /// \ingroup amrex_linear_containers
    //! Return a Vector of const raw pointers to the elements of \p a (FAB types only).
    template <class T, typename = typename T::FABType>
    [[nodiscard]] Vector<const T*> GetVecOfConstPtrs (const Vector<T>& a)
    {
        Vector<const T*> r;
        r.reserve(a.size());
        for (const auto& x : a) { r.push_back(&x); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of const raw pointers to the Array elements of \p a (FAB types only).
    template <class T, std::size_t N, typename = typename T::FABType>
    [[nodiscard]] Vector<Array<T,N> const*> GetVecOfConstPtrs (Vector<Array<T,N>> const& a)
    {
        Vector<Array<T,N> const*> r;
        r.reserve(a.size());
        for (auto& x : a) { r.push_back(&x); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of const raw pointers from a Vector of unique_ptr.
    template <class T>
    [[nodiscard]] Vector<const T*> GetVecOfConstPtrs (const Vector<std::unique_ptr<T> >& a)
    {
        Vector<const T*> r;
        r.reserve(a.size());
        for (const auto& x : a) { r.push_back(x.get()); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of const raw pointers from a Vector of raw pointers.
    template <class T, typename = typename T::FABType>
    [[nodiscard]] Vector<const T*> GetVecOfConstPtrs (const Vector<T*>& a)
    {
        return {a.begin(), a.end()};
    }

    /////////////////////////////////////////////////////////////

    /// \ingroup amrex_linear_containers
    //! Return a Vector of Vectors of raw pointers from a Vector of Vectors of unique_ptr.
    template <class T>
    [[nodiscard]] Vector<Vector<T*> > GetVecOfVecOfPtrs (const Vector<Vector<std::unique_ptr<T> > >& a)
    {
        Vector<Vector<T*> > r;
        r.reserve(a.size());
        for (const auto& x : a) { r.push_back(GetVecOfPtrs(x)); }
        return r;
    }

    /////////////////////////////////////////////////////////////

#ifdef AMREX_SPACEDIM
    /// \ingroup amrex_linear_containers
    //! Return a Vector of AMREX_SPACEDIM-arrays of raw pointers from a Vector of arrays of unique_ptr.
    template <class T>
    [[nodiscard]] Vector<std::array<T*,AMREX_SPACEDIM> >
    GetVecOfArrOfPtrs (const Vector<std::array<std::unique_ptr<T>,AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (const auto& x : a) { r.push_back(GetArrOfPtrs(x)); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of AMREX_SPACEDIM-arrays of const raw pointers from a Vector of arrays of unique_ptr.
    template <class T>
    [[nodiscard]] Vector<std::array<T const*,AMREX_SPACEDIM> >
    GetVecOfArrOfPtrsConst (const Vector<std::array<std::unique_ptr<T>,AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T const*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (const auto& x : a) { r.push_back(GetArrOfConstPtrs(x)); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of AMREX_SPACEDIM-arrays of const raw pointers from a Vector of arrays of unique_ptr.
    template <class T>
    [[nodiscard]] Vector<std::array<T const*,AMREX_SPACEDIM> >
    GetVecOfArrOfConstPtrs (const Vector<std::array<std::unique_ptr<T>,AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T const*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (const auto& x : a) { r.push_back(GetArrOfConstPtrs(x)); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of AMREX_SPACEDIM-arrays of const raw pointers from a Vector of arrays of FAB objects.
    template <class T, std::enable_if_t<IsFabArray<T>::value ||
                                               IsBaseFab<T>::value,
                                               int> = 0 >
    //! Convert the array-of-FAB container into raw pointer arrays.
    [[nodiscard]] Vector<std::array<T const*,AMREX_SPACEDIM> >
    GetVecOfArrOfConstPtrs (const Vector<std::array<T,AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T const*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (const auto& x : a) { r.push_back(GetArrOfConstPtrs(x)); }
        return r;
    }

    /// \ingroup amrex_linear_containers
    //! Return a Vector of AMREX_SPACEDIM-arrays of raw pointers from a Vector of arrays of FAB objects.
    template <class T, std::enable_if_t<IsFabArray<T>::value ||
                                               IsBaseFab<T>::value,
                                               int> = 0 >
    [[nodiscard]] Vector<std::array<T*, AMREX_SPACEDIM> >
    GetVecOfArrOfPtrs(Vector<std::array<T, AMREX_SPACEDIM> >& a)
    {
        Vector<std::array<T*, AMREX_SPACEDIM> > r;
        r.reserve(a.size());
        for (auto &x: a) { r.push_back(GetArrOfPtrs(x)); }
        return r;
    }
#endif

    /////////////////////////////////////////////////////////////

    /// \ingroup amrex_linear_containers
    //! Set all raw pointers in \p a to nullptr.
    template <class T>
    void FillNull (Vector<T*>& a)
    {
        std::for_each(a.begin(), a.end(), [](T*& p) { p = nullptr; });
    }

    /// \ingroup amrex_linear_containers
    //! Reset all unique_ptr elements of \p a to empty.
    template <class T>
    void FillNull (Vector<std::unique_ptr<T> >& a)
    {
        std::for_each(a.begin(), a.end(), [](std::unique_ptr<T>& p) { p.reset(); });
    }

    /////////////////////////////////////////////////////////////

    /// \ingroup amrex_linear_containers
    //! Sort \p vec and remove duplicate elements (requires operator< and operator==).
    template <class T>
    void RemoveDuplicates (Vector<T>& vec) {
        std::sort(vec.begin(), vec.end());
        auto it = std::unique(vec.begin(), vec.end());
        vec.erase(it, vec.end());
    }

    /// \cond DOXYGEN_IGNORE
    namespace detail {
        template <class T, class H>
        std::size_t removeDupDoit (Vector<T>& vec, std::size_t start, std::size_t stop)
        {
            std::size_t N = stop-start;
            if (N < 2) { return stop; }

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
                if (hashhash != hash) { // data[hash] not in correct has poision, thus will yield it's position
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
                if (data[i] != sentinel && i == hasher(data[i]) % N) {
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
    /// \endcond

    /// \ingroup amrex_linear_containers
    /**
     * \brief Remove duplicate elements from \p vec using a hash function.
     *
     * Uses a hash-based algorithm that does not require \p T to be sortable.
     * Requires only equality comparison and the provided hash functor \p H.
     *
     * \tparam T Element type.
     * \tparam H Hash functor type satisfying the \c std::hash interface.
     */
    template <class T, class H>
    void RemoveDuplicates (Vector<T>& vec) {
        // https://stackoverflow.com/questions/1532819/algorithm-efficient-way-to-remove-duplicate-integers-from-an-array
        std::size_t pos = detail::removeDupDoit<T,H>(vec, 0, vec.size());
        vec.erase(vec.begin()+pos, vec.end());
    }
}

#endif
