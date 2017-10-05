
#ifndef BL_INTVECT_H
#define BL_INTVECT_H

#include <iosfwd>

#include <AMReX_ccse-mpi.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BLassert.H>

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
	std::size_t operator()(const IntVect& vec) const
        {
	    AMREX_D_DECL(std::size_t ret0 = vec[0], ret1 = vec[1], ret2 = vec[2]);
#if BL_SPACEDIM == 1
	    return ret0;
#elif BL_SPACEDIM == 2
	    return ret0 ^ (ret1 << shift1);
#else
	    return ret0 ^ (ret1 << shift1) ^ (ret2 << shift2);	    
#endif
	}
    };

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

    /**
     Initializes Zero and Unit.
    */
    static int InitStatics();

    ///
    /**
    * \brief Construct an IntVect whose components are all zero.
    */
    IntVect () { AMREX_D_EXPR(vect[0] = 0, vect[1] = 0, vect[2] = 0); }

    /**
    * \brief Construct an IntVect given the specific values for its
    * coordinates.  AMREX_D_DECL is a macro that sets the constructor
    * to take BL_SPACEDIM arguments.
    */
    IntVect (AMREX_D_DECL(int i, int j, int k)) { AMREX_D_EXPR(vect[0] = i, vect[1] = j, vect[2] = k); }

    /**
    * \brief Construct an IntVect setting the coordinates to the
    * corresponding values in the integer array a.
    */
    explicit IntVect (const int* a);

    /**
    * \brief Construct an IntVect from an Array<int>.  It is an error if
    * the Array<int> doesn' t have the same dimension as this
    * IntVect.
    */
    explicit IntVect (const Array<int>& a);

    //! Returns a reference to the i'th coordinate of the IntVect.
    // dtor, copy-ctor, copy-op=, move-ctor, and move-op= are compiler generated.

    int& operator[] (int i) { BL_ASSERT(i>=0 && i < BL_SPACEDIM); return vect[i]; }

    //! Returns the i'th coordinate of the IntVect.
    const int& operator[] (int i) const { BL_ASSERT(i>=0 && i < BL_SPACEDIM); return vect[i]; }

    //! Set i'th coordinate of IntVect to val.
    IntVect& setVal (int i, int val)
    {
        BL_ASSERT(i >=0 && i < BL_SPACEDIM); vect[i] = val; return *this;
    }

    /**
    * \brief Returns a const pointer to an array of coordinates of the
    * IntVect.  Useful for arguments to FORTRAN calls.
    */
    const int* getVect () const& { return vect; }
    int*       getVect () &      { return vect; }
    int*       getVect () && = delete;

    //! Returns true if all components are equal to the argument val.
    bool operator== (int val) const
    {
        return AMREX_D_TERM(vect[0] == val, && vect[1] == val, && vect[2] == val);
    }

    //! Returns true if any component is not equal to the argument val.
    bool operator!= (int val) const
    {
        return AMREX_D_TERM(vect[0] != val, || vect[1] != val, || vect[2] != val);
    }

    //! Returns true if this is equivalent to rhs.
    bool operator== (const IntVect& rhs) const
    {
        return AMREX_D_TERM(vect[0] == rhs[0], && vect[1] == rhs[1], && vect[2] == rhs[2]);
    }
    //! Returns true if this is different from rhs.
    bool operator!= (const IntVect& rhs) const
    {
        return AMREX_D_TERM(vect[0] != rhs[0], || vect[1] != rhs[1], || vect[2] != rhs[2]);
    }
    //! Return ture if this is lexicographically less than rhs.
    bool operator< (const IntVect& rhs) const
    {
#if (BL_SPACEDIM == 1)
    return vect[0] < rhs[0];
#elif (BL_SPACEDIM == 2)
    return (vect[1] < rhs[1]) || ((vect[1] == rhs[1]) && (vect[0] < rhs[0]));
#else
    return (vect[2] < rhs[2]) || ((vect[2] == rhs[2]) &&
	(  (vect[1] < rhs[1]) || ((vect[1] == rhs[1]) && (vect[0] < rhs[0])) ));
#endif
    }
    //! Return true if this is lexicographically less than or equal to rhs.
    bool operator<= (const IntVect& rhs) const
    {
	return !(rhs < *this);
    }
    //! Retrun true if this is lexicographically greater than rhs.
    bool operator> (const IntVect& rhs) const
    {
	return rhs < *this;
    }
    //! Retrun true if this is lexicographically greater than or equal to rhs.
    bool operator>= (const IntVect& rhs) const
    {
	return !(*this < rhs);
    }
    /**
    * \brief Returns true if this is less than argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */	
    bool allLT (const IntVect& rhs) const
    {
        return AMREX_D_TERM(vect[0] < rhs[0], && vect[1] < rhs[1], && vect[2] < rhs[2]);
    }
    /**
    * \brief Returns true if this is less than or equal to argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    bool allLE (const IntVect& rhs) const
    {
        return AMREX_D_TERM(vect[0] <= rhs[0], && vect[1] <= rhs[1], && vect[2] <= rhs[2]);
    }
    /**
    * \brief Returns true if this is greater than argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    bool allGT (const IntVect& rhs) const
    {
        return AMREX_D_TERM(vect[0] > rhs[0], && vect[1] > rhs[1], && vect[2] > rhs[2]);
    }
    /**
    * \brief Returns true if this is greater than or equal to argument for all components.
    * NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    */
    bool allGE (const IntVect& rhs) const
    {
        return AMREX_D_TERM(vect[0] >= rhs[0], && vect[1] >= rhs[1], && vect[2] >= rhs[2]);
    }
    //! Unary plus -- for completeness.
    const IntVect operator+ () const { return *this; }
    //! Unary Minus -- negates all components.
    IntVect operator- () const { return IntVect(AMREX_D_DECL(-vect[0], -vect[1], -vect[2] )); }
    //! Modifies IntVect by addition of a scalar to each component.
    IntVect& operator+= (int s)
    {
        AMREX_D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s); return *this;
    }
    //! Modifies IntVect by component-wise addition with argument.
    IntVect& operator+= (const IntVect& p)
    {
        AMREX_D_EXPR(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2]); return *this;
    }
    //! Modifies IntVect by multiplication of a scalar to each component.
    IntVect& operator*= (int s)
    {
        AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s); return *this;
    }
    //! Modifies IntVect by component-wise multiplication with argument.
    IntVect& operator*= (const IntVect& p)
    {
        AMREX_D_EXPR(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2]); return *this;
    }
    //! Modifies IntVect by division by a scalar to each component.
    IntVect& operator/= (int s)
    {
        AMREX_D_EXPR(vect[0] /= s, vect[1] /= s, vect[2] /= s); return *this;
    }
    //! Modifies IntVect by component-wise division with argument.
    IntVect& operator/= (const IntVect& p)
    {
        AMREX_D_EXPR(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2]); return *this;
    }
    //! Modifies IntVect by subtraction of a scalar to each component.
    IntVect& operator-= (int s)
    {
        AMREX_D_EXPR(vect[0] -= s, vect[1] -= s, vect[2] -= s); return *this;
    }
    //! Modifies IntVect by component-wise subtraction with argument.
    IntVect& operator-= (const IntVect& p)
    {
        AMREX_D_EXPR(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2]); return *this;
    }
    //! Returns component-wise sum of IntVect and argument.
    IntVect operator+ (const IntVect& p) const
    {
        return IntVect(AMREX_D_DECL(vect[0] + p[0], vect[1] + p[1], vect[2] + p[2]));
    }
    //! Return an IntVect that is this IntVect + s.
    IntVect operator+ (int s) const
    {
        return IntVect(AMREX_D_DECL(vect[0] + s, vect[1] + s, vect[2] + s));
    }
    //! Returns component-wise difference of IntVect and argument.
    IntVect operator- (const IntVect& p) const
    {
        return IntVect(AMREX_D_DECL(vect[0] - p[0], vect[1] - p[1], vect[2] - p[2]));
    }
    //! Return an IntVect that is this IntVect - s.
    IntVect operator- (int s) const
    {
        return IntVect(AMREX_D_DECL(vect[0] - s, vect[1] - s, vect[2] - s));
    }
    //! Returns component-wise product of IntVect and argument.
    IntVect operator* (const IntVect& p) const
    {
        return IntVect(AMREX_D_DECL(vect[0] * p[0], vect[1] * p[1], vect[2] * p[2]));
    }
    //! Returns component-wise product of IntVect and s.
    IntVect operator* (int s) const
    {
        return IntVect(AMREX_D_DECL(vect[0] * s, vect[1] * s, vect[2] * s));
    }
    //! Returns component-wise division of IntVect by argument.
    IntVect operator/ (const IntVect& p) const
    {
        return IntVect(AMREX_D_DECL(vect[0] / p[0], vect[1] / p[1], vect[2] / p[2]));
    }
    //! Returns component-wise division of IntVect by s.
    IntVect operator/ (int s) const
    {
        return IntVect(AMREX_D_DECL(vect[0] / s, vect[1] / s, vect[2] / s));
    }
    //! Modifies IntVect by taking component-wise min with argument.
    IntVect& min (const IntVect& p)
    {
        AMREX_D_EXPR(vect[0] = std::min(vect[0], p.vect[0]),
               vect[1] = std::min(vect[1], p.vect[1]),
               vect[2] = std::min(vect[2], p.vect[2]));
        return *this;
    }
    //! Modifies IntVect by taking component-wise max with argument.
    IntVect& max (const IntVect& p)
    {
        AMREX_D_EXPR(vect[0] = std::max(vect[0], p.vect[0]),
               vect[1] = std::max(vect[1], p.vect[1]),
               vect[2] = std::max(vect[2], p.vect[2]));
        return *this;
    }
    //! Modify IntVect by multiplying each coordinate by s.
    IntVect& scale (int s) { AMREX_D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s); return *this; }
    /**
    * \brief Modify IntVect by reflecting it in the plane defined by
    * the index ref_ix and with normal in the direction of idir.
    */
    IntVect& reflect (int ref_ix, int idir)
    {
        BL_ASSERT(idir >= 0 && idir < BL_SPACEDIM);
        vect[idir] = -vect[idir] + 2*ref_ix;
        return *this;
    }
    //! Modify IntVect by adding s to given coordinate.
    IntVect& shift (int coord, int s)
    {
        BL_ASSERT(coord >= 0 && coord < BL_SPACEDIM); vect[coord] += s; return *this;
    }
    //! Equivalent to shift(0,iv[0]).shift(1,iv[1]) ...
    IntVect& shift (const IntVect& iv) { *this += iv; return *this; }
    //! Modify IntVect by adding s to each coordinate.
    IntVect& diagShift (int s)
    {
        AMREX_D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s); return *this;
    }
    //! Modify IntVect by component-wise integer projection.
    IntVect& coarsen (const IntVect& p);
    //! Modify IntVect by component-wise integer projection.
    IntVect& coarsen (int p);

    /// for serialization
    static size_t 
    linearSize()
    {
      size_t retval = SpaceDim*sizeof(int);
      return retval;
    }

    /// for serialization
    void 
    linearOut(void* a_buffer ) const
    {
      int* intbuf = (int *) a_buffer;
      const IntVect& iv = *this;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        *intbuf = iv[idir];
        intbuf++;
      }
    }

    /// for serialization
    void 
    linearIn(void* a_buffer )
    {
      int* intbuf = (int *) a_buffer;
      IntVect& iv = *this;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        iv[idir] = *intbuf;
        intbuf++;
      }
    }

    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose BL_SPACEDIM arguments are set to zero (0).
    * Figuratively, it is the zero vector in BL_SPACEDIM-dimensional space.
    * It is provided as a convenient way to specify the zero vector.
    */
    static const IntVect& TheZeroVector ();
    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose BL_SPACEDIM arguments are set to one (1).
    * Figuratively, it is the unit vector in BL_SPACEDIM-dimensional space.
    * It is provided as a convenient way to specify the unit vector.
    */
    static const IntVect& TheUnitVector ();
    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose BL_SPACEDIM arguments are set to zero except that
    * the d-direction is set to one. 
    */
    static const IntVect& TheDimensionVector (int d);
    /**
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose BL_SPACEDIM arguments are set to IndexType::NODE.
    * It is provided as a convenience to our users when defining Boxes.
    */
    static const IntVect& TheNodeVector ();
    /** 
    * \brief This static member function returns a reference to a constant IntVect
    * object, all of whose BL_SPACEDIM arguments are set to IndexType::CELL.
    * It is provided as a convenience to our users when defining Boxes.
    */
    static const IntVect& TheCellVector ();

    static const IntVect& TheMaxVector ();
    static const IntVect& TheMinVector ();

private:

    friend std::ostream& operator<< (std::ostream& os, const IntVect& iv);
    friend std::istream& operator>> (std::istream& is, IntVect& iv);

    int vect[BL_SPACEDIM];
};

//
// Static initialization.  Gotta make sure there are no static object
// definitions above here (except possibly stuff in headers).  Otherwise,
// the danger is that some static object's constructor calls IntVect::Zero or
// IntVect::Unit -- the very things the following definition is supposed to
// initialize.
//
static int s_dummyForIntVectH = IntVect::InitStatics();

    /**
    * \brief Returns a basis vector in the given coordinate direction;
    * eg.  IntVect BASISV(1) == (0,1,0).  Note that the coordinate
    * directions are zero based.
    */
    IntVect BASISV (int dir);
    /**
    * \brief Returns the IntVect that is the component-wise minimum of two
    * argument IntVects.
    */
    IntVect min (const IntVect& p1, const IntVect& p2);
    /**
    * \brief Returns the IntVect that is the component-wise maximum of two
    * argument IntVects.
    */
    IntVect max (const IntVect& p1, const IntVect& p2);
    /**
    * \brief Returns a IntVect obtained by multiplying each of the
    * components of this IntVect by s.
    */
    IntVect scale (const IntVect& p, int s);
    /**
    * \brief Returns an IntVect that is the reflection of input in the
    * plane which passes through ref_ix and normal to the
    * coordinate direction idir.
    */
    IntVect reflect(const IntVect& a, int ref_ix, int idir);
    /**
    * \brief Returns IntVect obtained by adding s to each of the
    * components of this IntVect.
    */
    IntVect diagShift (const IntVect& p, int s);
    /**
    * \brief Returns an IntVect that is the component-wise integer
    * projection of p by s.
    */
    IntVect coarsen (const IntVect& p, int s);
    /**
    * \brief Returns an IntVect which is the component-wise integer
    * projection of IntVect p1 by IntVect p2.
    */
    IntVect coarsen (const IntVect& p1, const IntVect& p2);


//! Returns p + s.
inline
IntVect operator+ (int s, const IntVect& p)
{
    return IntVect(AMREX_D_DECL(p[0] + s, p[1] + s, p[2] + s));
}
//! Returns -p + s.
inline
IntVect operator- (int s, const IntVect& p)
{
    return IntVect(AMREX_D_DECL(s - p[0], s - p[1], s - p[2]));
}
//! Returns p * s.
inline
IntVect operator* (int s, const IntVect& p)
{
    return IntVect(AMREX_D_DECL(s * p[0], s * p[1], s * p[2]));
}

}

#endif /*BL_INTVECT_H*/
