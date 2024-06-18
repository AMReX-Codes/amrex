#ifndef AMREX_INTVECT_H_
#define AMREX_INTVECT_H_
#include <AMReX_Config.H>

#include <AMReX_INT.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Dim3.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_Tuple.H>
#include <AMReX_TypeTraits.H>

#include <iosfwd>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <climits>
#include <algorithm>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int coarsen (int i, int ratio) noexcept
{
    switch (ratio) {
    case  1: return i;
    case  2: return (i<0) ? -std::abs(i+1)/2    -1 : i/2;
    case  4: return (i<0) ? -std::abs(i+1)/4    -1 : i/4;
    default: return (i<0) ? -std::abs(i+1)/ratio-1 : i/ratio;
    }
}

/**
* An Integer Vector in dim-Dimensional Space
*
* The class IntVectND is an implementation of an integer vector in a
* dim-dimensional space.  It represents a point in a discrete space.
* IntVectND values are accessed using the operator[] function, as for a normal
* C++ array.  In addition, the basic arithmetic operators have been overloaded
* to implement scaling and translation operations.
*/
template<int dim>
class IntVectND
{
public:
    static_assert(dim >= 1, "The number of dimensions of IntVectND must be positive");

    struct shift_hasher {
        std::size_t operator()(const IntVectND<dim>& vec) const noexcept
        {
            static constexpr unsigned shift1 = sizeof(size_t)>=8 ? 20 : 10;
            static constexpr unsigned shift2 = sizeof(size_t)>=8 ? 40 : 20;
            if constexpr (dim == 1) {
                return static_cast<std::size_t>(vec[0]);
            } else if constexpr (dim == 2) {
                return static_cast<std::size_t>(vec[0]) ^
                      (static_cast<std::size_t>(vec[1]) << shift1);
            } else if constexpr (dim == 3) {
                return static_cast<std::size_t>(vec[0]) ^
                      (static_cast<std::size_t>(vec[1]) << shift1) ^
                      (static_cast<std::size_t>(vec[2]) << shift2);
            } else {
                std::size_t seed = dim;
                // hash function from
                // https://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
                for (int i=0; i<dim; ++i) {
                    auto x = static_cast<unsigned int>(vec[i]);
                    x = ((x >> 16) ^ x) * 0x45d9f3b;
                    x = ((x >> 16) ^ x) * 0x45d9f3b;
                    x = (x >> 16) ^ x;
                    seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                }
                return seed;
            }
        }
    };

    using hasher = shift_hasher;

    ///
    /**
    * \brief Construct an IntVectND whose components are all zero.
    */
    constexpr IntVectND () noexcept {} // cannot use = default due to Clang bug // NOLINT

    /**
    * \brief Construct an IntVectND given the specific values for its
    * coordinates. The inputs for this constructor are N integers,
    * where N is equal to the number of dimensions of the IntVectND.
    */
    template <class...Args,
        std::enable_if_t<
            (sizeof...(Args)+2 == dim) &&
            IsConvertible_v<int, Args...>,
        int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr IntVectND (int i, int j, Args...ks) noexcept : vect{i, j, static_cast<int>(ks)...} {}

    /**
    * \brief Construct an IntVectND whose components are all the same.
    */
    AMREX_GPU_HOST_DEVICE
    explicit constexpr IntVectND (int s) noexcept {
        for (int i=0; i<dim; ++i) {
            vect[i] = s;
        }
    }

    /**
    * \brief Construct an IntVectND setting the coordinates to the
    * corresponding values in the integer array a.
    */
    AMREX_GPU_HOST_DEVICE
    explicit IntVectND (const int* a) noexcept {
        for (int i=0; i<dim; ++i) {
            vect[i] = a[i];
        }
    }

    /**
    * \brief Construct an IntVectND from an Vector<int>.  It is an error if
    * the Vector<int> doesn't have the same dimension as this
    * IntVectND.
    */
    explicit IntVectND (const Vector<int>& a) noexcept {
        BL_ASSERT(a.size() == dim);
        for (int i=0; i<dim; ++i) {
            vect[i] = a[i];
        }
    }

    /**
    * \brief Construct an IntVectND from an Array<int,dim>.
    */
    explicit IntVectND (const Array<int,dim>& a) noexcept {
        for (int i=0; i<dim; ++i) {
            vect[i] = a[i];
        }
    }

    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    explicit constexpr IntVectND (Dim3 const& a) noexcept {
        vect[0] = a.x;
        if constexpr (dim >= 2) {
            vect[1] = a.y;
        }
        if constexpr (dim == 3) {
            vect[2] = a.z;
        }
    }

    // dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

    template <int N=dim, std::enable_if_t<( 1<=N && N<=3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 dim3 () const noexcept {
        if constexpr (dim == 1) {
            return Dim3{vect[0],0,0};
        } else if constexpr (dim == 2) {
            return Dim3{vect[0],vect[1],0};
        } else {
            return Dim3{vect[0],vect[1],vect[2]};
        }
    }

#if __cplusplus >= 201402L
    template< typename T = int >
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Array<T, dim>
    toArray () const noexcept {
        Array<T, dim> ret {};
        for (int i=0; i<dim; ++i) {
            ret[i] = T(vect[i]);
        }
        return ret;
    }
#endif

  ///
  /**
     Sum of all components of this IntVectND.
  */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int sum () const noexcept
    {
        int retval = vect[0];
        for (int i=1; i<dim; ++i) {
            retval += vect[i];
        }
        return retval;
    }

    //! maximum (no absolute values) value
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int max () const noexcept
    {
        int retval = vect[0];
        for (int i=1; i<dim; ++i) {
            retval = retval > vect[i] ? retval : vect[i];
        }
        return retval;
    }

    //! minimum (no absolute values) value
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int min () const noexcept
    {
        int retval = vect[0];
        for (int i=1; i<dim; ++i) {
            retval = retval < vect[i] ? retval : vect[i];
        }
        return retval;
    }

    //return coordinate with largest value
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int maxDir(bool a_doAbsValue) const noexcept;

    //! Returns a reference to the i'th coordinate of the IntVectND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    int& operator[] (int i) noexcept { BL_ASSERT(i>=0 && i < dim); return vect[i]; }

    //! Returns the i'th coordinate of the IntVectND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    const int& operator[] (int i) const noexcept { BL_ASSERT(i>=0 && i < dim); return vect[i]; }

    //! Returns a reference to the i'th coordinate of the IntVectND. Used by structured bindings.
    template<std::size_t i>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    int& get () noexcept {static_assert(0<=i && i<dim); return vect[i];}

    //! Returns a reference to the i'th coordinate of the IntVectND. Used by structured bindings.
    template<std::size_t i>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    const int& get () const noexcept {static_assert(0<=i && i<dim); return vect[i];}

    //! Returns a pointer to the first element of the IntVectND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    int* begin () noexcept { return &vect[0]; }

    //! Returns a pointer to the first element of the IntVectND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    const int* begin () const noexcept { return &vect[0]; }

    //! Returns a pointer to the (last+1) element of the IntVectND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    int* end () noexcept { return &vect[dim]; }

    //! Returns a pointer to the (last+1) element of the IntVectND.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    const int* end () const noexcept { return &vect[dim]; }

    //! Set i'th coordinate of IntVectND to val.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND& setVal (int i, int val) noexcept
    {
        BL_ASSERT(i>=0 && i<dim); vect[i] = val; return *this;
    }

    /**
    * \brief Returns a const pointer to an array of coordinates of the
    * IntVectND.  Useful for arguments to FORTRAN calls.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const int* getVect () const& noexcept { return vect; }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int*       getVect () &      noexcept { return vect; }
    AMREX_GPU_HOST_DEVICE
    int*       getVect () && = delete;

    //! Returns true if all components are equal to the argument val.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator== (int val) const noexcept
    {
        bool retval = vect[0] == val;
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] == val;
        }
        return retval;
    }

    //! Returns true if any component is not equal to the argument val.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator!= (int val) const noexcept
    {
        bool retval = vect[0] != val;
        for (int i=1; i<dim; ++i) {
            retval = retval || vect[i] != val;
        }
        return retval;
    }

    //! Returns true if this is equivalent to rhs.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator== (const IntVectND<dim>& rhs) const noexcept
    {
        bool retval = vect[0] == rhs[0];
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] == rhs[i];
        }
        return retval;
    }
    //! Returns true if this is different from rhs.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator!= (const IntVectND<dim>& rhs) const noexcept
    {
        bool retval = vect[0] != rhs[0];
        for (int i=1; i<dim; ++i) {
            retval = retval || vect[i] != rhs[i];
        }
        return retval;
    }
    //! Return true if this is lexicographically less than rhs.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator< (const IntVectND<dim>& rhs) const noexcept
    {
        for (int i=dim-1; i>=0; --i) {
            if (vect[i] < rhs[i]) {
                return true;
            } else if (vect[i] > rhs[i]) {
                return false;
            }
        }
        return false;
    }
    //! Return true if this is lexicographically less than or equal to rhs.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator<= (const IntVectND<dim>& rhs) const noexcept
    {
        return !(rhs < *this);
    }
    //! Return true if this is lexicographically greater than rhs.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator> (const IntVectND<dim>& rhs) const noexcept
    {
        return rhs < *this;
    }
    //! Return true if this is lexicographically greater than or equal to rhs.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator>= (const IntVectND<dim>& rhs) const noexcept
    {
        return !(*this < rhs);
    }
    /**
    * \brief Returns true if this is less than argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool allLT (const IntVectND<dim>& rhs) const noexcept
    {
        bool retval = vect[0] < rhs[0];
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] < rhs[i];
        }
        return retval;
    }
    /**
    * \brief Returns true if this is less than argument for all components.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool allLT (int rhs) const noexcept
    {
        bool retval = vect[0] < rhs;
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] < rhs;
        }
        return retval;
    }
    /**
    * \brief Returns true if this is less than or equal to argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool allLE (const IntVectND<dim>& rhs) const noexcept
    {
        bool retval = vect[0] <= rhs[0];
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] <= rhs[i];
        }
        return retval;
    }
    /**
    * \brief Returns true if this is less than or equal to argument for all components.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool allLE (int rhs) const noexcept
    {
        bool retval = vect[0] <= rhs;
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] <= rhs;
        }
        return retval;
    }
    /**
    * \brief Returns true if this is greater than argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool allGT (const IntVectND<dim>& rhs) const noexcept
    {
        bool retval = vect[0] > rhs[0];
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] > rhs[i];
        }
        return retval;
    }
    /**
    * \brief Returns true if this is greater than argument for all components.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool allGT (int rhs) const noexcept
    {
        bool retval = vect[0] > rhs;
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] > rhs;
        }
        return retval;
    }
    /**
    * \brief Returns true if this is greater than or equal to argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool allGE (const IntVectND<dim>& rhs) const noexcept
    {
        bool retval = vect[0] >= rhs[0];
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] >= rhs[i];
        }
        return retval;
    }
    /**
    * \brief Returns true if this is greater than or equal to argument for all components.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool allGE (int rhs) const noexcept
    {
        bool retval = vect[0] >= rhs;
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] >= rhs;
        }
        return retval;
    }
    //! Unary plus -- for completeness.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator+ () const noexcept { return *this; }
    //! Unary Minus -- negates all components.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator- () const noexcept {
        IntVectND<dim> retval(0);
        for (int i=0; i<dim; ++i) {
            retval[i] = -vect[i];
        }
        return retval;
    }
    //! Modifies IntVectND by addition of a scalar to each component.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& operator+= (int s) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] += s;
        }
        return *this;
    }
    //! Modifies IntVectND by component-wise addition with argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& operator+= (const IntVectND<dim>& p) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] += p[i];
        }
        return *this;
    }
    //! Modifies IntVectND by multiplication of a scalar to each component.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& operator*= (int s) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] *= s;
        }
        return *this;
    }
    //! Modifies IntVectND by component-wise multiplication with argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& operator*= (const IntVectND<dim>& p) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] *= p[i];
        }
        return *this;
    }
    //! Modifies IntVectND by division by a scalar to each component.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& operator/= (int s) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] /= s;
        }
        return *this;
    }
    //! Modifies IntVectND by component-wise division with argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& operator/= (const IntVectND<dim>& p) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] /= p[i];
        }
        return *this;
    }
    //! Modifies IntVectND by subtraction of a scalar to each component.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& operator-= (int s) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] -= s;
        }
        return *this;
    }
    //! Modifies IntVectND by component-wise subtraction with argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& operator-= (const IntVectND<dim>& p) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] -= p[i];
        }
        return *this;
    }
    //! Returns component-wise sum of IntVectND<and argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator+ (const IntVectND<dim>& p) const noexcept
    {
        IntVectND<dim> retval = *this;
        return retval += p;
    }
    //! Return an IntVectND that is this IntVectND + s.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator+ (int s) const noexcept
    {
        IntVectND<dim> retval = *this;
        return retval += s;
    }
    //! Returns component-wise difference of IntVectND and argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator- (const IntVectND<dim>& p) const noexcept
    {
        IntVectND<dim> retval = *this;
        return retval -= p;
    }
    //! Return an IntVectND that is this IntVectND - s.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator- (int s) const noexcept
    {
        IntVectND<dim> retval = *this;
        return retval -= s;
    }
    //! Returns component-wise product of IntVectND and argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator* (const IntVectND<dim>& p) const noexcept
    {
        IntVectND<dim> retval = *this;
        return retval *= p;
    }
    //! Returns component-wise product of IntVectND and s.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator* (int s) const noexcept
    {
        IntVectND<dim> retval = *this;
        return retval *= s;
    }
    //! Returns component-wise division of IntVectND by argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator/ (const IntVectND<dim>& p) const noexcept
    {
        IntVectND<dim> retval = *this;
        return retval /= p;
    }
    //! Returns component-wise division of IntVectND by s.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim> operator/ (int s) const noexcept
    {
        IntVectND<dim> retval = *this;
        return retval /= s;
    }
    //! Modifies IntVectND by taking component-wise min with argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& min (const IntVectND<dim>& p) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] = (vect[i] < p.vect[i] ? vect[i] : p.vect[i]);
        }
        return *this;
    }
    //! Modifies IntVectND by taking component-wise max with argument.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& max (const IntVectND<dim>& p) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] = (vect[i] > p.vect[i] ? vect[i] : p.vect[i]);
        }
        return *this;
    }
    //! Modify IntVectND by multiplying each coordinate by s.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& scale (int s) noexcept {
        for (int i=0; i<dim; ++i) {
            vect[i] *= s;
        }
        return *this;
    }
    /**
    * \brief Modify IntVectND by reflecting it in the plane defined by
    * the index ref_ix and with normal in the direction of idir.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& reflect (int ref_ix, int idir) noexcept
    {
        BL_ASSERT(idir >= 0 && idir < dim);
        vect[idir] = -vect[idir] + 2*ref_ix;
        return *this;
    }
    //! Modify IntVectND by adding s to given coordinate.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& shift (int coord, int s) noexcept
    {
        BL_ASSERT(coord >= 0 && coord < dim); vect[coord] += s; return *this;
    }
    //! Equivalent to shift(0,iv[0]).shift(1,iv[1]) ...
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& shift (const IntVectND<dim>& iv) noexcept { *this += iv; return *this; }
    //! Modify IntVectND by adding s to each coordinate.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& diagShift (int s) noexcept
    {
        for (int i=0; i<dim; ++i) {
            vect[i] += s;
        }
        return *this;
    }
    //! Modify IntVectND<dim> by component-wise integer projection.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& coarsen (const IntVectND<dim>& p) noexcept;
    //! Modify IntVectND<dim> by component-wise integer projection.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<dim>& coarsen (int p) noexcept;

    /**
    * \brief This static member function returns a reference to a constant IntVectND
    * object, all of whose dim arguments are set to zero (0).
    * Figuratively, it is the zero vector in dim-dimensional space.
    * It is provided as a convenient way to specify the zero vector.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IntVectND<dim> TheZeroVector () noexcept {
        return IntVectND<dim>(0);
    }
    /**
    * \brief This static member function returns a reference to a constant IntVectND
    * object, all of whose dim arguments are set to one (1).
    * Figuratively, it is the unit vector in dim-dimensional space.
    * It is provided as a convenient way to specify the unit vector.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IntVectND<dim> TheUnitVector () noexcept {
        return IntVectND<dim>(1);
    }
    /**
    * \brief This static member function returns a reference to a constant IntVectND
    * object, all of whose dim arguments are set to zero except that
    * the d-direction is set to one.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IntVectND<dim> TheDimensionVector (int d) noexcept {
        IntVectND<dim> retval(0);
        retval[d] = 1;
        return retval;
    }
    /**
    * \brief This static member function returns a reference to a constant IntVectND
    * object, all of whose dim arguments are set to IndexType::NODE.
    * It is provided as a convenience to our users when defining Boxes.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IntVectND<dim> TheNodeVector () noexcept {
        return IntVectND<dim>(1);
    }
    /**
    * \brief This static member function returns a reference to a constant IntVectND
    * object, all of whose dim arguments are set to IndexType::CELL.
    * It is provided as a convenience to our users when defining Boxes.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IntVectND<dim> TheCellVector () noexcept {
        return IntVectND<dim>(0);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IntVectND<dim> TheMaxVector () noexcept {
        return IntVectND<dim>(std::numeric_limits<int>::max());
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr IntVectND<dim> TheMinVector () noexcept {
        return IntVectND<dim>(std::numeric_limits<int>::lowest());
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr std::size_t size () noexcept {
        return static_cast<std::size_t>(dim);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int isize () noexcept {
        return dim;
    }

    using value_type = int;

    /**
    * \brief Returns a new IntVectND of size new_dim and
    * assigns the first new_dim values of this IntVectND to it.
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<new_dim> shrink () const noexcept {
        static_assert(new_dim <= dim);
        return IntVectND<new_dim>(this->begin());
    }

    /**
    * \brief Returns a new IntVectND of size new_dim and
    * assigns all values of this IntVectND to it and fill_extra to the remaining elements.
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<new_dim> expand (int fill_extra=0) const noexcept {
        static_assert(new_dim >= dim);
        IntVectND<new_dim> retval(fill_extra);
        for (int i=0; i<dim; ++i) {
            retval[i] = vect[i];
        }
        return retval;
    }

    /**
    * \brief Returns a new IntVectND of size new_dim
    * by either shrinking or expanding this IntVectND
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    IntVectND<new_dim> resize (int fill_extra=0) const noexcept {
        if constexpr (new_dim > dim) {
            return expand<new_dim>(fill_extra);
        } else {
            return shrink<new_dim>();
        }
    }

    /**
       This is an IntVect all of whose components are equal to zero.
    */
    static const IntVectND<dim> Zero;

    /**
       This is an IntVect all of whose components are equal to one.
    */
    static const IntVectND<dim> Unit;

private:

    int vect[dim] = {};
};

template <int dim>
inline constexpr const IntVectND<dim> IntVectND<dim>::Zero{0};

template <int dim>
inline constexpr const IntVectND<dim> IntVectND<dim>::Unit{1};

// Template deduction guide for IntVectND
template<std::size_t dim>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
IntVectND(const Array<int, dim>&) -> IntVectND<dim>;

// Template deduction guide for IntVectND
template <class...Args,
    std::enable_if_t<
        IsConvertible_v<int, Args...>,
    int> = 0>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
IntVectND(int, int, Args...) -> IntVectND<sizeof...(Args)+2>;

using IntVect = IntVectND<AMREX_SPACEDIM>;

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>&
IntVectND<dim>::coarsen (int s) noexcept
{
    BL_ASSERT(s > 0);
    switch (s) {
    case 1:
        break;
    case 2:
        for (int i=0; i<dim; ++i) {
            vect[i] = (vect[i]<0) ? -std::abs(vect[i]+1)/2-1 : vect[i]/2;
        }
        break;
    case 4:
        for (int i=0; i<dim; ++i) {
            vect[i] = (vect[i]<0) ? -std::abs(vect[i]+1)/4-1 : vect[i]/4;
        }
        break;
    default:
        for (int i=0; i<dim; ++i) {
            vect[i] = (vect[i]<0) ? -std::abs(vect[i]+1)/s-1 : vect[i]/s;
        }
    }
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>&
IntVectND<dim>::coarsen (const IntVectND<dim>& p) noexcept
{
    BL_ASSERT(p.allGT(0));
    for (int i=0; i<dim; ++i) {
        vect[i] = amrex::coarsen(vect[i], p.vect[i]);
    }
    return *this;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
IntVectND<dim>::maxDir(bool a_doAbsValue) const noexcept
{
    int retval = 0;
    if(a_doAbsValue)
    {
        int maxval = std::abs((*this)[0]);
        for(int idir = 1; idir < dim; idir++)
        {
            int curval = std::abs((*this)[idir]);
            if(curval > maxval)
            {
                maxval = curval;
                retval = idir;
            }
        }
    }
    else
    {
        int maxval = (*this)[0];
        for(int idir = 1; idir < dim; idir++)
        {
            int curval = (*this)[idir];
            if(curval > maxval)
            {
                maxval = curval;
                retval = idir;
            }
        }
    }
    return retval;
}

//! Returns p + s.
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> operator+ (int s, const IntVectND<dim>& p) noexcept
{
    IntVectND<dim> retval = p;
    for (int i=0; i<dim; ++i) {
        retval[i] = s + retval[i];
    }
    return retval;
}

//! Returns -p + s.
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
IntVectND<dim> operator- (int s, const IntVectND<dim>& p) noexcept
{
    IntVectND<dim> retval = p;
    for (int i=0; i<dim; ++i) {
        retval[i] = s - retval[i];
    }
    return retval;
}

//! Returns p * s.
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim> operator* (int s, const IntVectND<dim>& p) noexcept
{
    IntVectND<dim> retval = p;
    for (int i=0; i<dim; ++i) {
        retval[i] = s * retval[i];
    }
    return retval;
}

/**
 * \brief Returns the IntVectND that is the component-wise minimum of two
 * argument IntVectNDs.
 */
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
min (const IntVectND<dim>& p1, const IntVectND<dim>& p2) noexcept
{
    IntVectND<dim> p(p1);
    p.min(p2);
    return p;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
elemwiseMin (const IntVectND<dim>& p1, const IntVectND<dim>& p2) noexcept
{
    IntVectND<dim> p(p1);
    p.min(p2);
    return p;
}

/**
 * \brief Returns the IntVectND that is the component-wise maximum of two
 * argument IntVectNDs.
 */
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
max (const IntVectND<dim>& p1, const IntVectND<dim>& p2) noexcept
{
    IntVectND<dim> p(p1);
    p.max(p2);
    return p;
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
elemwiseMax (const IntVectND<dim>& p1, const IntVectND<dim>& p2) noexcept
{
    IntVectND<dim> p(p1);
    p.max(p2);
    return p;
}

    /**
    * \brief Returns a basis vector in the given coordinate direction;
    * eg.  IntVectND<3> BASISV<3>(1) == (0,1,0).  Note that the coordinate
    * directions are zero based.
    */
template<int dim = AMREX_SPACEDIM>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
BASISV (int dir) noexcept
{
    BL_ASSERT(dir >= 0 && dir < dim);
    IntVectND<dim> tmp(0);
    tmp[dir] = 1;
    return tmp;
}

    /**
    * \brief Returns a IntVectND obtained by multiplying each of the
    * components of this IntVectND by s.
    */
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
scale (const IntVectND<dim>& p, int s) noexcept
{
    return IntVectND<dim>(p).scale(s);
}

    /**
    * \brief Returns an IntVectND that is the reflection of input in the
    * plane which passes through ref_ix and normal to the
    * coordinate direction idir.
    */
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
reflect (const IntVectND<dim>& a, int ref_ix, int idir) noexcept
{
    return IntVectND<dim>(a).reflect(ref_ix, idir);
}

    /**
    * \brief Returns IntVectND obtained by adding s to each of the
    * components of this IntVectND.
    */
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
diagShift (const IntVectND<dim>& p, int s) noexcept
{
    return p + s;
}

    /**
    * \brief Returns an IntVectND that is the component-wise integer
    * projection of p by s.
    */
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
coarsen (const IntVectND<dim>& p, int s) noexcept
{
    BL_ASSERT(s > 0);
    IntVectND<dim> v = p;
    v.coarsen(s);
    return v;
}

    /**
    * \brief Returns an IntVectND which is the component-wise integer
    * projection of IntVectND p1 by IntVectND p2.
    */
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
IntVectND<dim>
coarsen (const IntVectND<dim>& p1, const IntVectND<dim>& p2) noexcept
{
    IntVectND<dim> v = p1;
    v.coarsen(p2);
    return v;
}

template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Dim3 refine (Dim3 const& coarse, IntVectND<dim> const& ratio) noexcept
{
    if constexpr (dim == 1) {
        return Dim3{coarse.x*ratio[0], coarse.y, coarse.z};
    } else if constexpr (dim == 2) {
        return Dim3{coarse.x*ratio[0], coarse.y*ratio[1], coarse.z};
    } else {
        return Dim3{coarse.x*ratio[0], coarse.y*ratio[1], coarse.z*ratio[2]};
    }
}

template<int dim, std::enable_if_t<( 1<=dim && dim<=3 ), int> = 0>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Dim3 coarsen (Dim3 const& fine, IntVectND<dim> const& ratio) noexcept
{
    if constexpr (dim == 1) {
        return Dim3{amrex::coarsen(fine.x, ratio[0]),
                    fine.y,
                    fine.z};
    } else if constexpr (dim == 2) {
        return Dim3{amrex::coarsen(fine.x, ratio[0]),
                    amrex::coarsen(fine.y, ratio[1]),
                    fine.z};
    } else {
        return Dim3{amrex::coarsen(fine.x, ratio[0]),
                    amrex::coarsen(fine.y, ratio[1]),
                    amrex::coarsen(fine.z, ratio[2])};
    }
}

namespace detail {
    std::ostream& int_vector_write (std::ostream& os, const int* iv, int dim);
    std::istream& int_vector_read (std::istream& is, int* iv, int dim);

    template<int dim>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    void IntVectCat_imp (int*& dst, const IntVectND<dim>& src) noexcept {
        for (int i=0; i<dim; ++i) {
            dst[i] = src[i];
        }
        dst += dim;
    }

    template<int dim>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    void IntVectSplit_imp2 (IntVectND<dim>& dst, const int*& src) noexcept {
        for (int i=0; i<dim; ++i) {
            dst[i] = src[i];
        }
        src += dim;
    }

    template<class T, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T IntVectSplit_imp (T& retval, std::index_sequence<Ns...>, const int * src) noexcept {
        (IntVectSplit_imp2(std::get<Ns>(retval), src), ...);
        return retval;
    }

    template<int...dims>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    int get_sum () {
        return (0 + ... + dims);
    }
}

template<int dim>
std::ostream&
operator<< (std::ostream& os, const IntVectND<dim>& iv)
{
    return detail::int_vector_write(os, iv.begin(), dim);
}

template<int dim>
std::istream&
operator>> (std::istream& is, IntVectND<dim>& iv)
{
    return detail::int_vector_read(is, iv.begin(), dim);
}

/**
* \brief Returns a IntVectND obtained by concatenating the input IntVectNDs.
* The dimension of the return value equals the sum of the dimensions of the inputted IntVectNDs.
*/
template<int d, int...dims>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr IntVectND<detail::get_sum<d, dims...>()>
IntVectCat (const IntVectND<d>& v, const IntVectND<dims>&...vects) noexcept {
    IntVectND<detail::get_sum<d, dims...>()> retval (0);
    int* dst = retval.begin();
    detail::IntVectCat_imp(dst, v);
    (detail::IntVectCat_imp(dst, vects), ...);
    return retval;
}

/**
* \brief Returns a tuple of IntVectND obtained by splitting the input IntVectND
* according to the dimensions specified by the template arguments.
*/
template<int d, int...dims>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr GpuTuple<IntVectND<d>, IntVectND<dims>...>
IntVectSplit (const IntVectND<detail::get_sum<d, dims...>()>& v) noexcept {
    GpuTuple<IntVectND<d>, IntVectND<dims>...> retval;
    return detail::IntVectSplit_imp(retval,
                                    std::make_index_sequence<1 + sizeof...(dims)>(),
                                    v.begin());
}

/**
* \brief Returns a new IntVectND of size new_dim and
* assigns the first new_dim values of iv to it.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr IntVectND<new_dim>
IntVectShrink (const IntVectND<old_dim>& iv) noexcept {
    return iv.template shrink<new_dim>();
}

/**
* \brief Returns a new IntVectND of size new_dim and
* assigns all values of iv to it and fill_extra to the remaining elements.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr IntVectND<new_dim>
IntVectExpand (const IntVectND<old_dim>& iv, int fill_extra=0) noexcept {
    return iv.template expand<new_dim>(fill_extra);
}

/**
* \brief Returns a new IntVectND of size new_dim
* by either shrinking or expanding iv.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr IntVectND<new_dim>
IntVectResize (const IntVectND<old_dim>& iv, int fill_extra=0) noexcept {
    return iv.template resize<new_dim>(fill_extra);
}

} // namespace amrex

// Spcialize std::tuple_size for IntVectND. Used by structured bindings.
template<int dim>
struct std::tuple_size<amrex::IntVectND<dim>> {
    static constexpr std::size_t value = dim;
};

// Spcialize std::tuple_element for IntVectND. Used by structured bindings.
template<std::size_t s, int dim>
struct std::tuple_element<s, amrex::IntVectND<dim>> {
    using type = int;
};

#endif /*AMREX_INTVECT_H*/
