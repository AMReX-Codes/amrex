#ifndef AMREX_REALVECT_H_
#define AMREX_REALVECT_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_IntVect.H>
#include <AMReX_Utility.H>
#include <AMReX_Math.H>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iosfwd>
#include <vector>

namespace amrex
{

/**
    \brief A Real vector in dim-dimensional space

    The class RealVectND is an implementation of a Real vector in a
    dim-dimensional space.
    RealVectND values are accessed using the operator[] function, as for a normal
    C++ array.  In addition, the basic arithmetic operators have been overloaded
    to implement scaling and translation operations.
*/
template <int dim>
class RealVectND {
    public:
    /**
        \name Constructors and Accessors
    */
    /*@{*/

    ///
    /**
        Construct a RealVectND whose components are zero.
    */
    constexpr RealVectND () noexcept {}  // NOLINT

    explicit RealVectND (const std::vector<Real>& vr) noexcept {
        BL_ASSERT(vr.size() == dim);
        for (int i = 0; i < dim; ++i) {
            vect[i] = vr[i];
        }
    }

    ///
    /**
        Construct a RealVectND given the specific values for its
        coordinates.
    */
    template <class... Args,
        std::enable_if_t<
            (sizeof...(Args) + 2 == dim) &&
            IsConvertible_v<Real, Args...>,
        int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr RealVectND (Real i, Real j, Args... ks) noexcept
        : vect{i, j, static_cast<Real>(ks)...} {}

    AMREX_GPU_HOST_DEVICE
    explicit constexpr RealVectND (Real s) noexcept {
        for (int i = 0; i < dim; ++i) {
            vect[i] = s;
        }
    }

    ///
    /**
        Construct a RealVectND setting the coordinates to the corresponding values in
        the Real array a.
    */
    AMREX_GPU_HOST_DEVICE
    explicit RealVectND (const Real* a) noexcept {
        for (int i = 0; i < dim; ++i) {
            vect[i] = a[i];
        }
    }

    ///
    /**
        Construct a RealVectND from an IntVectND by coercing each component
        from <tt>int</tt> to Real.
    */
    AMREX_GPU_HOST_DEVICE constexpr
    RealVectND (const IntVectND<dim>& iv) noexcept {
        for (int i = 0; i < dim; ++i) {
            vect[i] = static_cast<Real>(iv[i]);
        }
    }

    ///
    /**
        Construct a RealVectND from a GpuArray
    */
   AMREX_GPU_HOST_DEVICE constexpr
   RealVectND (const GpuArray<Real, dim> ga) noexcept {
       for (int i = 0; i < dim; ++i) {
           vect[i] = ga[i];
       }
   }

    ///
    /**
        Prevent rvalue assignment to the <i>i</i>'th coordinate of the
        RealVectND.
    */
    Real& operator[] (int i) && = delete;

    ///
    /**
        Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
        RealVectND.
    */
    AMREX_GPU_HOST_DEVICE inline constexpr
    Real& operator[] (int i) & noexcept {
        AMREX_ASSERT(i >= 0 && i < dim);
        return vect[i];
    }

    ///
    /**
        Returns the <i>i</i>'th coordinate of the RealVectND.
    */
    AMREX_GPU_HOST_DEVICE inline constexpr
    const Real& operator[] (int i) const& noexcept {
        AMREX_ASSERT(i >= 0 && i < dim);
        return vect[i];
    }

    //! Returns a reference to the i'th coordinate of the RealVectND. Used by structured bindings.
    template<std::size_t i>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    Real& get () noexcept {static_assert(0<=i && i<dim); return vect[i];}

    //! Returns a reference to the i'th coordinate of the RealVectND. Used by structured bindings.
    template<std::size_t i>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    const Real& get () const noexcept {static_assert(0<=i && i<dim); return vect[i];}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline constexpr
    GpuArray<Real, dim> to_array () const noexcept {
        GpuArray<Real, dim> retval{};
        for (int i = 0; i < dim; ++i) {
            retval[i] = vect[i];
        }
        return retval;
    }

    /*@}*/

    /**
        \name Iterators
    */
    /*@{*/

    ///
    /**
        Return pointer to the first component.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline constexpr
    Real* begin () noexcept { return vect; }

    ///
    /**
        Return pointer to the first component.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline constexpr
    const Real* begin () const noexcept { return vect; }

    ///
    /**
        Return pointer to the one past last component
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline constexpr
    Real* end () noexcept { return vect + dim; }

    ///
    /**
        Return pointer to the one past last component
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline constexpr
    const Real* end () const noexcept {
        return vect + dim;
    }

    /*@}*/

    /**
        \name Comparison Operators
    */
    /*@{*/

    ///
    /**
        Returns true if this RealVectND is equivalent to argument RealVectND.  All
        comparisons between analogous components must be satisfied.
    */
    AMREX_GPU_HOST_DEVICE constexpr
    bool operator== (const RealVectND& p) const noexcept {
        bool retval = vect[0] == p[0];
        for (int i=1; i<dim; ++i) {
            retval = retval && vect[i] == p[i];
        }
        return retval;
    }

    ///
    /**
        Returns true if this RealVectND is different from argument RealVectND.
        All comparisons between analogous components must be satisfied.
    */
    AMREX_GPU_HOST_DEVICE constexpr
    bool operator!= (const RealVectND& p) const noexcept {
        bool retval = vect[0] != p[0];
        for (int i=1; i<dim; ++i) {
            retval = retval || vect[i] != p[i];
        }
        return retval;
    }

    ///
    /**
        Returns true if this RealVectND is less than argument RealVectND.  All
        comparisons between analogous components must be satisfied.  Note
        that, since the comparison is component-wise, it is possible for
        an RealVectND to be neither greater than, less than, nor equal to
        another.
    */
    AMREX_GPU_HOST_DEVICE inline bool operator< (const RealVectND& p) const noexcept;

    ///
    /**
        Returns true if this RealVectND is less than or equal to argument
        RealVectND.  All comparisons between analogous components must be
        satisfied.  Note that, since the comparison is component-wise, it
        is possible for an RealVectND to be neither greater than or equal
        to, less than or equal to, nor equal to another.
    */
    AMREX_GPU_HOST_DEVICE inline bool operator<= (const RealVectND& p) const noexcept;

    ///
    /**
        Returns true if this RealVectND is greater than argument RealVectND.
        All comparisons between analogous components must be satisfied.
        Note that, since the comparison is component-wise, it is possible
        for an RealVectND to be neither greater than, less than, nor equal
        to another.
    */
    AMREX_GPU_HOST_DEVICE inline bool operator> (const RealVectND& p) const noexcept;

    ///
    /**
        Returns true if this RealVectND is greater than or equal to argument
        RealVectND.  All comparisons between analogous components must be
        satisfied.  Note that, since the comparison is component-wise, it
        is possible for an RealVectND to be neither greater than or equal
        to, less than or equal to, nor equal to another.
    */

    AMREX_GPU_HOST_DEVICE inline bool operator>= (const RealVectND& p) const noexcept;

    /*@}*/

    /**
        \name Arithmetic Operators
    */
    /*@{*/

    ///
    /**
        Modifies this RealVectND by addition of a scalar to each component.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& operator+= (Real s) noexcept;

    ///
    /**
        Returns a RealVectND that is this RealVectND with a scalar s added
        to each component.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND operator+ (Real s) const noexcept;

    ///
    /**
        Modifies this RealVectND by component-wise addition by argument.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& operator+= (const RealVectND& p) noexcept;

    ///
    /**
        Modifies this RealVectND by subtraction of a scalar from each
        component.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& operator-= (Real s) noexcept;

    ///
    /**
        Modifies this RealVectND by component-wise subtraction by argument.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& operator-= (const RealVectND& p) noexcept;

    ///
    /**
        Returns a RealVectND that is this RealVectND with a scalar s subtracted
        from each component.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND operator- (Real s) const noexcept;

    ///
    /**
        Modifies this RealVectND by multiplying each component by a scalar.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& operator*= (Real s) noexcept;

    ///
    /**
        Returns the result of the scalar product with another RealVectND
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline Real dotProduct (
        const RealVectND& a_rhs) const noexcept;

    ///
    /**
        Returns the result of the cross product with another RealVectND
    */
    template <int N=dim, std::enable_if_t<( N==3 ), int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline RealVectND crossProduct (
        const RealVectND& a_rhs) const noexcept;

    ///
    /**
        Modifies this RealVectND by component-wise multiplication by argument.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& operator*= (const RealVectND& p) noexcept;

    ///
    /**
        Returns a RealVectND that is this RealVectND with each component
        multiplied by a scalar.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND operator* (Real s) const noexcept;

    ///
    /**
        Modifies this RealVectND by dividing each component by a scalar.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& operator/= (Real s) noexcept;

    ///
    /**
        Modifies this RealVectND by component-wise division by argument.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& operator/= (const RealVectND& p) noexcept;

    ///
    /**
        Returns a RealVectND that is this RealVectND with each component
        divided by a scalar.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND operator/ (Real s) const noexcept;

    ///
    /**
        Modifies this RealVectND by multiplying each component by a scalar.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& scale (Real s) noexcept;

    ///
    /**
        Returns an IntVectND<dim> whose components are the std::floor of the RealVectND
        components.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline IntVectND<dim> floor () const noexcept;

    ///
    /**
        Returns an IntVectND<dim> whose components are the std::floor of the RealVectND
        components.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline IntVectND<dim> ceil () const noexcept;

    ///
    /**
        Returns an IntVectND<dim> whose components are the std::floor of the RealVectND
        components.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline IntVectND<dim> round () const noexcept;

    /*@}*/

    /**
        \name Other arithmetic operators
    */
    /*@{*/

    ///
    /**
        Modifies this RealVectND by taking component-wise min with RealVectND
        argument.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& min (const RealVectND& p) noexcept;

    ///
    /**
        Modifies this RealVectND by taking component-wise max with RealVectND
        argument.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND& max (const RealVectND& p) noexcept;

    /*@}*/

    /**
        \name Unary operators
    */
    /*@{*/

    ///
    /**
        Unary plus -- for completeness.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND operator+ () const noexcept;

    ///
    /**
        Unary minus -- negates all components of this RealVectND.
    */
    AMREX_GPU_HOST_DEVICE inline RealVectND operator- () const noexcept;

    ///
    /**
        Sum of all components of this RealVectND.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline Real sum () const noexcept;

    ///
    /**
        sqrt(sum squares)
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline Real vectorLength () const noexcept;

    ///
    /**
        sum squares--no square root
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline Real radSquared () const noexcept;

    ///
    /**
        Product of all components of this RealVectND.
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline Real product () const noexcept;

    ///
    /**
        Component with the minimum value of this RealVectND (returns 0 if they are all the same).
        a_doAbs : if true then take the absolute value before comparing
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline int minDir (const bool& a_doAbs) const noexcept;

    ///
    /**
        Component with the maximum value of this RealVectND (returns 0 if they are all the same).
        a_doAbs : if true then take the absolute value before comparing
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline int maxDir (const bool& a_doAbs) const noexcept;

    /*@}*/

    /**
        \name Data pointer functions
    */
    /*@{*/

    ///
    /**
        Only for sending stuff to Fortran
    */
    [[nodiscard]] AMREX_GPU_HOST_DEVICE const Real* dataPtr () const noexcept { return vect; }

    ///
    /**
        Only for sending stuff to Fortran
    */
    AMREX_GPU_HOST_DEVICE
    Real* dataPtr () noexcept { return vect; }

    /*@}*/

    /**
        \name Constants
    */
    /*@{*/

    ///
    /**
        This is a RealVectND all of whose components are equal to zero.
    */
    static const RealVectND Zero;

    ///
    /**
        This is a RealVectND all of whose components are equal to one.
    */
    static const RealVectND Unit;

    /*@}*/

    /**
        \name Static public members functions
    */
    /*@{*/

    ///
    /**
        The zero vector in dim-dimensional space.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr RealVectND TheZeroVector () noexcept {
        return RealVectND{0.};
    }

    /**
        The unit vector in dim-dimensional space.
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr RealVectND TheUnitVector () noexcept {
        return RealVectND{1.};
    }

    /*@}*/

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr std::size_t size () noexcept {
        return static_cast<std::size_t>(dim);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr int isize () noexcept {
        return dim;
    }

    using value_type = Real;

    /**
    * \brief Returns a new RealVectND of size new_dim and
    * assigns the first new_dim values of this RealVectND to it.
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVectND<new_dim> shrink () const noexcept {
        static_assert(new_dim <= dim);
        return RealVectND<new_dim>(this->begin());
    }

    /**
    * \brief Returns a new RealVectND of size new_dim and
    * assigns all values of this RealVectND to it and fill_extra to the remaining elements.
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVectND<new_dim> expand (Real fill_extra=0.) const noexcept {
        static_assert(new_dim >= dim);
        RealVectND<new_dim> retval(fill_extra);
        for (int i=0; i<dim; ++i) {
            retval[i] = vect[i];
        }
        return retval;
    }

    /**
    * \brief Returns a new RealVectND of size new_dim
    * by either shrinking or expanding this RealVectND
    */
    template<int new_dim>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    RealVectND<new_dim> resize (Real fill_extra=0.) const noexcept {
        if constexpr (new_dim > dim) {
            return expand<new_dim>(fill_extra);
        } else {
            return shrink<new_dim>();
        }
    }

protected:
    /**
        The individual components of this RealVectND.
    */
    Real vect[dim] = {};
};

// Template deduction guide for RealVectND
template <class...Args,
    std::enable_if_t<
        IsConvertible_v<Real, Args...>,
    int> = 0>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
RealVectND(Real, Real, Args...) -> RealVectND<sizeof...(Args)+2>;

// Template deduction guide for RealVectND
template<int dim>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
RealVectND(const IntVectND<dim>&) -> RealVectND<dim>;

// Template deduction guide for RealVectND
template<int dim>
AMREX_GPU_HOST_DEVICE // __device__ for HIP
RealVectND(const GpuArray<Real, dim>&) -> RealVectND<dim>;

template <int dim>
inline constexpr const RealVectND<dim> RealVectND<dim>::Zero{Real(0)};

template <int dim>
inline constexpr const RealVectND<dim> RealVectND<dim>::Unit{Real(1)};

using RealVect = RealVectND<AMREX_SPACEDIM>;

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator-= (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] -= s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator*= (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] *= s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator-= (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] -= p[i];
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator+ () const noexcept {
    return RealVectND<dim>(*this);
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator- () const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = -vect[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::scale (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] *= s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
IntVectND<dim>
RealVectND<dim>::floor () const noexcept {
    IntVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = static_cast<int>(std::floor(vect[i]));
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
IntVectND<dim>
RealVectND<dim>::ceil () const noexcept {
    IntVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = static_cast<int>(std::ceil(vect[i]));
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
IntVectND<dim>
RealVectND<dim>::round () const noexcept {
    IntVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = static_cast<int>(std::round(vect[i]));
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::sum () const noexcept {
    Real retval = vect[0];
    for (int i=1; i<dim; ++i) {
        retval += vect[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::vectorLength () const noexcept {
    Real len = this->radSquared();
    len = std::sqrt(len);
    return len;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::radSquared () const noexcept {
    Real retval = vect[0] * vect[0];
    for (int i=1; i<dim; ++i) {
        retval += vect[i] * vect[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::product () const noexcept {
    Real retval = vect[0];
    for (int i=1; i<dim; ++i) {
        retval *= vect[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
bool
RealVectND<dim>::operator< (const RealVectND<dim>& p) const noexcept {
    bool retval = vect[0] < p[0];
    for (int i=1; i<dim; ++i) {
        retval = retval && vect[i] < p[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
bool
RealVectND<dim>::operator<= (const RealVectND<dim>& p) const noexcept {
    bool retval = vect[0] <= p[0];
    for (int i=1; i<dim; ++i) {
        retval = retval && vect[i] <= p[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
bool
RealVectND<dim>::operator> (const RealVectND<dim>& p) const noexcept {
    bool retval = vect[0] > p[0];
    for (int i=1; i<dim; ++i) {
        retval = retval && vect[i] > p[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
bool
RealVectND<dim>::operator>= (const RealVectND<dim>& p) const noexcept {
    bool retval = vect[0] >= p[0];
    for (int i=1; i<dim; ++i) {
        retval = retval && vect[i] >= p[i];
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::min (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] = amrex::min(vect[i], p.vect[i]);
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::max (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] = amrex::max(vect[i], p.vect[i]);
    }
    return *this;
}

/**
    Returns the RealVectND that is the component-wise minimum of two
    argument RealVectNDs.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> min (const RealVectND<dim>& p1, const RealVectND<dim>& p2) noexcept {
    RealVectND<dim> p(p1);
    return p.min(p2);
}

/**
    Returns the RealVectND that is the component-wise maximum of two
    argument RealVectNDs.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> max (const RealVectND<dim>& p1, const RealVectND<dim>& p2) noexcept {
    RealVectND<dim> p(p1);
    return p.max(p2);
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
Real
RealVectND<dim>::dotProduct (const RealVectND<dim>& a_rhs) const noexcept {
    Real retval  = vect[0] * a_rhs.vect[0];
    for (int i=1; i<dim; ++i) {
        retval += vect[i] * a_rhs.vect[i];
    }
    return retval;
}

template <int dim>
template <int N, std::enable_if_t<( N==3 ), int>>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::crossProduct (const RealVectND<dim>& a_rhs) const noexcept {
    RealVectND<dim> tmp(vect[1] * a_rhs[2] - vect[2] * a_rhs[1],
                        vect[2] * a_rhs[0] - vect[0] * a_rhs[2],
                        vect[0] * a_rhs[1] - vect[1] * a_rhs[0]);
    return tmp;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator+= (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] += s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator+= (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] += p[i];
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator*= (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] *= p[i];
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator* (Real s) const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = vect[i] * s;
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator- (Real s) const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = vect[i] - s;
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator+ (Real s) const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = vect[i] + s;
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator/= (Real s) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] /= s;
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>&
RealVectND<dim>::operator/= (const RealVectND<dim>& p) noexcept {
    for (int i=0; i<dim; ++i) {
        vect[i] /= p[i];
    }
    return *this;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim>
RealVectND<dim>::operator/ (Real s) const noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = vect[i] / s;
    }
    return retval;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
int
RealVectND<dim>::minDir (const bool& a_doAbs) const noexcept {
    int mDir = 0;
    for (int idir = 0; idir < dim; idir++) {
        if (a_doAbs) {
            if (std::abs(vect[idir]) < std::abs(vect[mDir])) {
                mDir = idir;
            }
        } else {
            if (vect[idir] < vect[mDir]) {
                mDir = idir;
            }
        }
    }
    return mDir;
}

template <int dim>
AMREX_GPU_HOST_DEVICE inline
int
RealVectND<dim>::maxDir (const bool& a_doAbs) const noexcept {
    int mDir = 0;
    for (int idir = 0; idir < dim; idir++) {
        if (a_doAbs) {
            if (std::abs(vect[idir]) > std::abs(vect[mDir])) {
                mDir = idir;
            }
        } else {
            if (vect[idir] > vect[mDir]) {
                mDir = idir;
            }
        }
    }
    return mDir;
}

/**
    Returns a basis vector in the given coordinate direction.<br>
    In 2-D:<br>
    BASISREALV(0) == (1.,0.);
    BASISREALV(1) == (0.,1.).<br>
    In 3-D:<br>
    BASISREALV(0) == (1.,0.,0.);
    BASISREALV(1) == (0.,1.,0.);
    BASISREALV(2) == (0.,0.,1.).<br>
    Note that the coordinate directions are based at zero.
*/
template <int dim=AMREX_SPACEDIM>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> BASISREALV (int dir) noexcept {
    AMREX_ASSERT(dir >= 0 && dir < dim);
    RealVectND<dim> tmp{0.};
    tmp[dir] = 1;
    return tmp;
}

/**
    \name Arithmetic functions
*/
/*@{*/

/**
     Returns a RealVectND that is a RealVectND <i>p</i> with each component
    divided by a scalar <i>s</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator/ (Real s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s / p[i];
    }
    return retval;
}

/**
    Returns a RealVectND that is a RealVectND <i>p</i> with
    a scalar <i>s</i> added to each component.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator+ (Real s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = p[i] + s;
    }
    return retval;
}

/**
     Returns <i>s - p</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator- (Real s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s - p[i];
    }
    return retval;
}

/**
     Returns a RealVectND that is a RealVectND <i>p</i> with each component
    multiplied by a scalar <i>s</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator* (Real s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s * p[i];
    }
    return retval;
}

/**
     Returns component-wise quotient <i>p / s</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator/ (const RealVectND<dim>& s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s[i] / p[i];
    }
    return retval;
}

/**
     Returns component-wise sum of RealVectNDs <i>s</i> and <i>p</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator+ (const RealVectND<dim>& s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = p[i] + s[i];
    }
    return retval;
}

/**
     Returns <i>s - p</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator- (const RealVectND<dim>& s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s[i] - p[i];
    }
    return retval;
}

/**
     Returns component-wise product of <i>s</i> and <i>p</i>.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> operator* (const RealVectND<dim>& s, const RealVectND<dim>& p) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = p[i] * s[i];
    }
    return retval;
}

/**
    Returns a RealVectND obtained by multiplying each of the components
    of the given RealVectND by a scalar.
*/
template <int dim>
AMREX_GPU_HOST_DEVICE inline
RealVectND<dim> scale (const RealVectND<dim>& p, Real s) noexcept {
    RealVectND<dim> retval;
    for (int i=0; i<dim; ++i) {
        retval[i] = s * p[i];
    }
    return retval;
}

namespace detail {
    std::ostream& real_vector_write (std::ostream& os, const Real* p, int dim);
    std::istream& real_vector_read (std::istream& is, Real* p, int dim);

    template<int dim>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    void RealVectCat_imp (Real*& dst, const RealVectND<dim>& src) noexcept {
        for (int i=0; i<dim; ++i) {
            dst[i] = src[i];
        }
        dst += dim;
    }

    template<int dim>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    void RealVectSplit_imp2 (RealVectND<dim>& dst, const Real*& src) noexcept {
        for (int i=0; i<dim; ++i) {
            dst[i] = src[i];
        }
        src += dim;
    }

    template<class T, std::size_t...Ns>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr
    T RealVectSplit_imp (T& retval, std::index_sequence<Ns...>, const Real * src) noexcept {
        (RealVectSplit_imp2(amrex::get<Ns>(retval), src), ...);
        return retval;
    }
}

template<int dim>
std::ostream&
operator<< (std::ostream& os, const RealVectND<dim>& p)
{
    return detail::real_vector_write(os, p.begin(), dim);
}

template<int dim>
std::istream&
operator>> (std::istream& is, RealVectND<dim>& p)
{
    return detail::real_vector_read(is, p.begin(), dim);
}

/**
* \brief Returns a RealVectND obtained by concatenating the input RealVectNDs.
* The dimension of the return value equals the sum of the dimensions of the inputted RealVectNDs.
*/
template<int d, int...dims>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr RealVectND<detail::get_sum<d, dims...>()>
RealVectCat (const RealVectND<d>& v, const RealVectND<dims>&...vects) noexcept {
    RealVectND<detail::get_sum<d, dims...>()> retval {0.};
    Real* dst = retval.begin();
    detail::RealVectCat_imp(dst, v);
    (detail::RealVectCat_imp(dst, vects), ...);
    return retval;
}

/**
* \brief Returns a tuple of RealVectND obtained by splitting the input RealVectND
* according to the dimensions specified by the template arguments.
*/
template<int d, int...dims>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr GpuTuple<RealVectND<d>, RealVectND<dims>...>
RealVectSplit (const RealVectND<detail::get_sum<d, dims...>()>& v) noexcept {
    GpuTuple<RealVectND<d>, RealVectND<dims>...> retval;
    return detail::RealVectSplit_imp(retval,
                                    std::make_index_sequence<1 + sizeof...(dims)>(),
                                    v.begin());
}

/**
* \brief Returns a new RealVectND of size new_dim and
* assigns the first new_dim values of iv to it.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr RealVectND<new_dim>
RealVectShrink (const RealVectND<old_dim>& iv) noexcept {
    return iv.template shrink<new_dim>();
}

/**
* \brief Returns a new RealVectND of size new_dim and
* assigns all values of iv to it and fill_extra to the remaining elements.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr RealVectND<new_dim>
RealVectExpand (const RealVectND<old_dim>& iv, Real fill_extra=0) noexcept {
    return iv.template expand<new_dim>(fill_extra);
}

/**
* \brief Returns a new RealVectND of size new_dim
* by either shrinking or expanding iv.
*/
template<int new_dim, int old_dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
constexpr RealVectND<new_dim>
RealVectResize (const RealVectND<old_dim>& iv, Real fill_extra=0) noexcept {
    return iv.template resize<new_dim>(fill_extra);
}

/*@}*/

}  // namespace amrex

// Spcialize std::tuple_size for RealVectND. Used by structured bindings.
template<int dim>
struct std::tuple_size<amrex::RealVectND<dim>> {
    static constexpr std::size_t value = dim;
};

// Spcialize std::tuple_element for RealVectND. Used by structured bindings.
template<std::size_t s, int dim>
struct std::tuple_element<s, amrex::RealVectND<dim>> {
    using type = amrex::Real;
};

#endif
