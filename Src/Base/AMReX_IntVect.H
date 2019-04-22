
#ifndef BL_INTVECT_H
#define BL_INTVECT_H

#include <iosfwd>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <climits>
#include <algorithm>

#include <AMReX_ccse-mpi.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

namespace amrex {

/**
* An Integer Vector in SPACEDIM-Dimensional Space
*
* The class IntVect is an implementation of an integer vector in a
* SPACEDIM-dimensional space.  It represents a point in a discrete space.
* IntVect values are accessed using the operator[] function, as for a normal
* C++ array.  In addition, the basic arithmetic operators have been overloaded
* to implement scaling and translation operations.
*/

class IntVect
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<IntVect>::type();
public:

    class shift_hasher {
    private:
	static constexpr unsigned shift1 = sizeof(size_t)>=8 ? 20 : 10;
	static constexpr unsigned shift2 = sizeof(size_t)>=8 ? 40 : 20;
    public:
	std::size_t operator()(const IntVect& vec) const noexcept
        {
	    AMREX_D_DECL(std::size_t ret0 = vec[0], ret1 = vec[1], ret2 = vec[2]);
#if AMREX_SPACEDIM == 1
	    return ret0;
#elif AMREX_SPACEDIM == 2
	    return ret0 ^ (ret1 << shift1);
#else
	    return ret0 ^ (ret1 << shift1) ^ (ret2 << shift2);
#endif
	}
    };

    using hasher = shift_hasher;

    ///
    /**
       This is an IntVect all of whose components are equal to zero.
    */
    static const IntVect Zero;

    ///
    /**
       This is an IntVect all of whose components are equal to one.
    */
    static const IntVect Unit;

    ///
    /**
    * \brief Construct an IntVect whose components are all zero.
    */
    AMREX_GPU_HOST_DEVICE
    constexpr IntVect () noexcept : vect{AMREX_D_DECL(0,0,0)} {}

    /**
    * \brief Construct an IntVect given the specific values for its
    * coordinates.  AMREX_D_DECL is a macro that sets the constructor
    * to take AMREX_SPACEDIM arguments.
    */
    AMREX_GPU_HOST_DEVICE
    constexpr IntVect (AMREX_D_DECL(int i, int j, int k)) noexcept : vect{AMREX_D_DECL(i,j,k)} {}

#if (AMREX_SPACEDIM > 1)
    AMREX_GPU_HOST_DEVICE
    explicit constexpr IntVect (int i) noexcept : vect{AMREX_D_DECL(i,i,i)} {}
#endif

    /**
    * \brief Construct an IntVect setting the coordinates to the
    * corresponding values in the integer array a.
    */
    AMREX_GPU_HOST_DEVICE
    explicit IntVect (const int* a) noexcept : vect{AMREX_D_DECL(a[0],a[1],a[2])} {}

    /**
    * \brief Construct an IntVect from an Vector<int>.  It is an error if
    * the Vector<int> doesn' t have the same dimension as this
    * IntVect.
    */
    explicit IntVect (const Vector<int>& a) noexcept : vect{AMREX_D_DECL(a[0],a[1],a[2])} {
        BL_ASSERT(a.size() == AMREX_SPACEDIM);
    }

    // dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Dim3 dim3 () const noexcept {
#if (AMREX_SPACEDIM == 1)
        return Dim3{vect[0],0,0};
#elif (AMREX_SPACEDIM == 2)
        return Dim3{vect[0],vect[1],0};
#else
        return Dim3{vect[0],vect[1],vect[2]};
#endif
    }

  ///
  /**
     Sum of all components of this IntVect.
  */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int sum () const noexcept
    {
        int retval = AMREX_D_TERM((*this)[0], + (*this)[1], +(*this)[2]);
        return retval;
    }

    //! maximum (no absolute values) value
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int max () const noexcept
    {
#if (AMREX_SPACEDIM == 1)
        return vect[0];
#elif (AMREX_SPACEDIM == 2)
        return (vect[0] > vect[1] ? vect[0] : vect[1]);
#else
        int r = vect[0] > vect[1] ? vect[0] : vect[1];
        return (r > vect[2]) ? r : vect[2];
#endif
    }

    //! minimum (no absolute values) value
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int min () const noexcept
    {
#if (AMREX_SPACEDIM == 1)
        return vect[0];
#elif (AMREX_SPACEDIM == 2)
        return (vect[0] < vect[1] ? vect[0] : vect[1]);
#else
        int r = vect[0] < vect[1] ? vect[0] : vect[1];
        return (r < vect[2]) ? r : vect[2];
#endif
    }

    //return coordinate with largest value
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int maxDir(bool a_doAbsValue) const noexcept;

    //! Returns a reference to the i'th coordinate of the IntVect.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int& operator[] (int i) noexcept { BL_ASSERT(i>=0 && i < AMREX_SPACEDIM); return vect[i]; }

    //! Returns the i'th coordinate of the IntVect.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    const int& operator[] (int i) const noexcept { BL_ASSERT(i>=0 && i < AMREX_SPACEDIM); return vect[i]; }

    //! Set i'th coordinate of IntVect to val.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& setVal (int i, int val) noexcept
    {
        BL_ASSERT(i >=0 && i < AMREX_SPACEDIM); vect[i] = val; return *this;
    }

    /**
    * \brief Returns a const pointer to an array of coordinates of the
    * IntVect.  Useful for arguments to FORTRAN calls.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    const int* getVect () const& noexcept { return vect; }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int*       getVect () &      noexcept { return vect; }
    AMREX_GPU_HOST_DEVICE
    int*       getVect () && = delete;

    //! Returns true if all components are equal to the argument val.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator== (int val) const noexcept
    {
        return AMREX_D_TERM(vect[0] == val, && vect[1] == val, && vect[2] == val);
    }

    //! Returns true if any component is not equal to the argument val.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator!= (int val) const noexcept
    {
        return AMREX_D_TERM(vect[0] != val, || vect[1] != val, || vect[2] != val);
    }

    //! Returns true if this is equivalent to rhs.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator== (const IntVect& rhs) const noexcept
    {
        return AMREX_D_TERM(vect[0] == rhs[0], && vect[1] == rhs[1], && vect[2] == rhs[2]);
    }
    //! Returns true if this is different from rhs.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator!= (const IntVect& rhs) const noexcept
    {
        return AMREX_D_TERM(vect[0] != rhs[0], || vect[1] != rhs[1], || vect[2] != rhs[2]);
    }
    //! Return ture if this is lexicographically less than rhs.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator< (const IntVect& rhs) const noexcept
    {
#if (AMREX_SPACEDIM == 1)
        return vect[0] < rhs[0];
#elif (AMREX_SPACEDIM == 2)
        return (vect[1] < rhs[1]) || ((vect[1] == rhs[1]) && (vect[0] < rhs[0]));
#else
        return (vect[2] < rhs[2]) || ((vect[2] == rhs[2]) &&
            (  (vect[1] < rhs[1]) || ((vect[1] == rhs[1]) && (vect[0] < rhs[0])) ));
#endif
    }
    //! Return true if this is lexicographically less than or equal to rhs.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator<= (const IntVect& rhs) const noexcept
    {
	return !(rhs < *this);
    }
    //! Retrun true if this is lexicographically greater than rhs.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator> (const IntVect& rhs) const noexcept
    {
	return rhs < *this;
    }
    //! Retrun true if this is lexicographically greater than or equal to rhs.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool operator>= (const IntVect& rhs) const noexcept
    {
	return !(*this < rhs);
    }
    /**
    * \brief Returns true if this is less than argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool allLT (const IntVect& rhs) const noexcept
    {
        return AMREX_D_TERM(vect[0] < rhs[0], && vect[1] < rhs[1], && vect[2] < rhs[2]);
    }
    /**
    * \brief Returns true if this is less than or equal to argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool allLE (const IntVect& rhs) const noexcept
    {
        return AMREX_D_TERM(vect[0] <= rhs[0], && vect[1] <= rhs[1], && vect[2] <= rhs[2]);
    }
    /**
    * \brief Returns true if this is greater than argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool allGT (const IntVect& rhs) const noexcept
    {
        return AMREX_D_TERM(vect[0] > rhs[0], && vect[1] > rhs[1], && vect[2] > rhs[2]);
    }
    /**
    * \brief Returns true if this is greater than or equal to argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    bool allGE (const IntVect& rhs) const noexcept
    {
        return AMREX_D_TERM(vect[0] >= rhs[0], && vect[1] >= rhs[1], && vect[2] >= rhs[2]);
    }
    //! Unary plus -- for completeness.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    const IntVect operator+ () const noexcept { return *this; }
    //! Unary Minus -- negates all components.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator- () const noexcept { return IntVect(AMREX_D_DECL(-vect[0], -vect[1], -vect[2] )); }
    //! Modifies IntVect by addition of a scalar to each component.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& operator+= (int s) noexcept
    {
        AMREX_D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s); return *this;
    }
    //! Modifies IntVect by component-wise addition with argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& operator+= (const IntVect& p) noexcept
    {
        AMREX_D_EXPR(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2]); return *this;
    }
    //! Modifies IntVect by multiplication of a scalar to each component.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& operator*= (int s) noexcept
    {
        AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s); return *this;
    }
    //! Modifies IntVect by component-wise multiplication with argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& operator*= (const IntVect& p) noexcept
    {
        AMREX_D_EXPR(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2]); return *this;
    }
    //! Modifies IntVect by division by a scalar to each component.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& operator/= (int s) noexcept
    {
        AMREX_D_EXPR(vect[0] /= s, vect[1] /= s, vect[2] /= s); return *this;
    }
    //! Modifies IntVect by component-wise division with argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& operator/= (const IntVect& p) noexcept
    {
        AMREX_D_EXPR(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2]); return *this;
    }
    //! Modifies IntVect by subtraction of a scalar to each component.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& operator-= (int s) noexcept
    {
        AMREX_D_EXPR(vect[0] -= s, vect[1] -= s, vect[2] -= s); return *this;
    }
    //! Modifies IntVect by component-wise subtraction with argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& operator-= (const IntVect& p) noexcept
    {
        AMREX_D_EXPR(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2]); return *this;
    }
    //! Returns component-wise sum of IntVect and argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator+ (const IntVect& p) const noexcept
    {
        return IntVect(AMREX_D_DECL(vect[0] + p[0], vect[1] + p[1], vect[2] + p[2]));
    }
    //! Return an IntVect that is this IntVect + s.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator+ (int s) const noexcept
    {
        return IntVect(AMREX_D_DECL(vect[0] + s, vect[1] + s, vect[2] + s));
    }
    //! Returns component-wise difference of IntVect and argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator- (const IntVect& p) const noexcept
    {
        return IntVect(AMREX_D_DECL(vect[0] - p[0], vect[1] - p[1], vect[2] - p[2]));
    }
    //! Return an IntVect that is this IntVect - s.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator- (int s) const noexcept
    {
        return IntVect(AMREX_D_DECL(vect[0] - s, vect[1] - s, vect[2] - s));
    }
    //! Returns component-wise product of IntVect and argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator* (const IntVect& p) const noexcept
    {
        return IntVect(AMREX_D_DECL(vect[0] * p[0], vect[1] * p[1], vect[2] * p[2]));
    }
    //! Returns component-wise product of IntVect and s.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator* (int s) const noexcept
    {
        return IntVect(AMREX_D_DECL(vect[0] * s, vect[1] * s, vect[2] * s));
    }
    //! Returns component-wise division of IntVect by argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator/ (const IntVect& p) const noexcept
    {
        return IntVect(AMREX_D_DECL(vect[0] / p[0], vect[1] / p[1], vect[2] / p[2]));
    }
    //! Returns component-wise division of IntVect by s.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect operator/ (int s) const noexcept
    {
        return IntVect(AMREX_D_DECL(vect[0] / s, vect[1] / s, vect[2] / s));
    }
    //! Modifies IntVect by taking component-wise min with argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& min (const IntVect& p) noexcept
    {
        AMREX_D_EXPR(vect[0] = (vect[0] < p.vect[0] ? vect[0] : p.vect[0]),
                     vect[1] = (vect[1] < p.vect[1] ? vect[1] : p.vect[1]),
                     vect[2] = (vect[2] < p.vect[2] ? vect[2] : p.vect[2]));
        return *this;
    }
    //! Modifies IntVect by taking component-wise max with argument.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& max (const IntVect& p) noexcept
    {
        AMREX_D_EXPR(vect[0] = (vect[0] > p.vect[0] ? vect[0] : p.vect[0]),
                     vect[1] = (vect[1] > p.vect[1] ? vect[1] : p.vect[1]),
                     vect[2] = (vect[2] > p.vect[2] ? vect[2] : p.vect[2]));
        return *this;
    }
    //! Modify IntVect by multiplying each coordinate by s.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& scale (int s) noexcept { AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s); return *this; }
    /**
    * \brief Modify IntVect by reflecting it in the plane defined by
    * the index ref_ix and with normal in the direction of idir.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& reflect (int ref_ix, int idir) noexcept
    {
        BL_ASSERT(idir >= 0 && idir < AMREX_SPACEDIM);
        vect[idir] = -vect[idir] + 2*ref_ix;
        return *this;
    }
    //! Modify IntVect by adding s to given coordinate.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& shift (int coord, int s) noexcept
    {
        BL_ASSERT(coord >= 0 && coord < AMREX_SPACEDIM); vect[coord] += s; return *this;
    }
    //! Equivalent to shift(0,iv[0]).shift(1,iv[1]) ...
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& shift (const IntVect& iv) noexcept { *this += iv; return *this; }
    //! Modify IntVect by adding s to each coordinate.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& diagShift (int s) noexcept
    {
        AMREX_D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s); return *this;
    }
    //! Modify IntVect by component-wise integer projection.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& coarsen (const IntVect& p) noexcept;
    //! Modify IntVect by component-wise integer projection.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    IntVect& coarsen (int p) noexcept;

    //! for serialization
    static size_t AMREX_INLINE
    linearSize() noexcept
    {
        size_t retval = AMREX_SPACEDIM*sizeof(int);
        return retval;
    }

    //! for serialization
    void linearOut(void* a_buffer ) const noexcept;
    void linearIn(void* a_buffer ) noexcept;

    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose AMREX_SPACEDIM arguments are set to zero (0).
    * Figuratively, it is the zero vector in AMREX_SPACEDIM-dimensional space.
    * It is provided as a convenient way to specify the zero vector.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IntVect TheZeroVector () noexcept {
        return IntVect(0);
    }
    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose AMREX_SPACEDIM arguments are set to one (1).
    * Figuratively, it is the unit vector in AMREX_SPACEDIM-dimensional space.
    * It is provided as a convenient way to specify the unit vector.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IntVect TheUnitVector () noexcept {
        return IntVect(1);
    }
    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose AMREX_SPACEDIM arguments are set to zero except that
    * the d-direction is set to one.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IntVect TheDimensionVector (int d) noexcept {
        return  (d==0) ? IntVect(AMREX_D_DECL(1,0,0)) :
               ((d==1) ? IntVect(AMREX_D_DECL(0,1,0)) :
                         IntVect(AMREX_D_DECL(0,0,1)));
    }
    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose AMREX_SPACEDIM arguments are set to IndexType::NODE.
    * It is provided as a convenience to our users when defining Boxes.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IntVect TheNodeVector () noexcept {
        return IntVect(1);
    }
    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose AMREX_SPACEDIM arguments are set to IndexType::CELL.
    * It is provided as a convenience to our users when defining Boxes.
    */
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IntVect TheCellVector () noexcept {
        return IntVect(0);
    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IntVect TheMaxVector () noexcept {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        return IntVect(std::numeric_limits<int>::max());
#else
        return IntVect(INT_MAX);
#endif
    }
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr IntVect TheMinVector () noexcept {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
        return IntVect(std::numeric_limits<int>::lowest());
#else
        return IntVect(INT_MIN);
#endif
    }

private:

    friend std::ostream& operator<< (std::ostream& os, const IntVect& iv);
    friend std::istream& operator>> (std::istream& is, IntVect& iv);

    int vect[AMREX_SPACEDIM];
};

AMREX_INLINE
IntVect&
IntVect::coarsen (int s) noexcept
{
    BL_ASSERT(s > 0);
    return this->coarsen(IntVect(AMREX_D_DECL(s,s,s)));
}

AMREX_INLINE
IntVect&
IntVect::coarsen (const IntVect& p) noexcept
{
    BL_ASSERT(p.allGT(IntVect::TheZeroVector()));
    if (p != 1) {
        for (int i = 0; i <AMREX_SPACEDIM; ++i)
        {
            const int s = p.vect[i];
            vect[i] = ((vect[i]<0) ? -std::abs(vect[i]+1)/s-1 : vect[i]/s);
        }
    }
    return *this;
}

AMREX_INLINE
int
IntVect::maxDir(bool a_doAbsValue) const noexcept
{
    int retval = 0;
    if(a_doAbsValue)
    {
        int maxval = std::abs((*this)[0]);
        for(int idir = 1; idir < SpaceDim; idir++)
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
        for(int idir = 1; idir < SpaceDim; idir++)
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

AMREX_INLINE
void
IntVect::linearOut(void* a_buffer ) const noexcept
{
    int* intbuf = (int *) a_buffer;
    const IntVect& iv = *this;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
        *intbuf = iv[idir];
        intbuf++;
    }
}

AMREX_INLINE
void
IntVect::linearIn(void* a_buffer ) noexcept
{
    int* intbuf = (int *) a_buffer;
    IntVect& iv = *this;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
        iv[idir] = *intbuf;
        intbuf++;
    }
}

//! Returns p + s.
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect operator+ (int s, const IntVect& p) noexcept
{
    return IntVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));
}
//! Returns -p + s.
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
AMREX_GPU_HOST_DEVICE
IntVect operator- (int s, const IntVect& p) noexcept
{
    return IntVect(AMREX_D_DECL(s - p[0], s - p[1], s - p[2]));
}
//! Returns p * s.
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect operator* (int s, const IntVect& p) noexcept
{
    return IntVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
}

/**
 * \brief Returns the IntVect that is the component-wise minimum of two
 * argument IntVects.
 */
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect
min (const IntVect& p1, const IntVect& p2) noexcept
{
    IntVect p(p1);
    p.min(p2);
    return p;
}

/**
 * \brief Returns the IntVect that is the component-wise maximum of two
 * argument IntVects.
 */
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect
max (const IntVect& p1, const IntVect& p2) noexcept
{
    IntVect p(p1);
    p.max(p2);
    return p;
}

    /**
    * \brief Returns a basis vector in the given coordinate direction;
    * eg.  IntVect BASISV(1) == (0,1,0).  Note that the coordinate
    * directions are zero based.
    */
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect
BASISV (int dir) noexcept
{
    BL_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM);
    IntVect tmp;
    tmp[dir] = 1;
    return tmp;
}

    /**
    * \brief Returns a IntVect obtained by multiplying each of the
    * components of this IntVect by s.
    */
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect
scale (const IntVect& p, int s) noexcept
{
    return IntVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
}

    /**
    * \brief Returns an IntVect that is the reflection of input in the
    * plane which passes through ref_ix and normal to the
    * coordinate direction idir.
    */
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect
reflect (const IntVect& a, int ref_ix, int idir) noexcept
{
    BL_ASSERT(idir >= 0 && idir < AMREX_SPACEDIM);
    IntVect b(a);
    b[idir] = -b[idir] + 2*ref_ix;
    return b;
}

    /**
    * \brief Returns IntVect obtained by adding s to each of the
    * components of this IntVect.
    */
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect
diagShift (const IntVect& p, int s) noexcept
{
    return IntVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));
}

    /**
    * \brief Returns an IntVect that is the component-wise integer
    * projection of p by s.
    */
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect
coarsen (const IntVect& p, int s) noexcept
{
    BL_ASSERT(s > 0);
    IntVect v = p;
    v.coarsen(IntVect(AMREX_D_DECL(s,s,s)));
    return v;
}

    /**
    * \brief Returns an IntVect which is the component-wise integer
    * projection of IntVect p1 by IntVect p2.
    */
AMREX_GPU_HOST_DEVICE
AMREX_INLINE
IntVect
coarsen (const IntVect& p1, const IntVect& p2) noexcept
{
    IntVect v = p1;
    v.coarsen(p2);
    return v;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Dim3 refine (Dim3 const& coarse, IntVect const& ratio) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return {coarse.x*ratio[0], coarse.y, coarse.z};
#elif (AMREX_SPACEDIM == 2)
    return {coarse.x*ratio[0], coarse.y*ratio[1], coarse.z};
#else
    return {coarse.x*ratio[0], coarse.y*ratio[1], coarse.z*ratio[2]};
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
int coarsen (int i, int ratio) noexcept
{
    return (i<0) ? -std::abs(i+1)/ratio-1 : i/ratio;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Dim3 coarsen (Dim3 const& fine, IntVect const& ratio) noexcept
{
#if (AMREX_SPACEDIM == 1)
    return {amrex::coarsen(fine.x,ratio[0]), fine.y, fine.z};
#elif (AMREX_SPACEDIM == 2)
    return {amrex::coarsen(fine.x,ratio[0]), amrex::coarsen(fine.y, ratio[1]), fine.z};
#else
    return {amrex::coarsen(fine.x, ratio[0]),
            amrex::coarsen(fine.y, ratio[1]),
            amrex::coarsen(fine.z, ratio[2])};
#endif
}

}

#endif /*BL_INTVECT_H*/
