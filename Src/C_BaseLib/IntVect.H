
#ifndef BL_INTVECT_H
#define BL_INTVECT_H

#include <winstd.H>
#include <iosfwd>

#include <ccse-mpi.H>
#include <SPACE.H>
#include <Array.H>
#include <BLassert.H>

//
// An Integer Vector in SPACEDIM-Dimensional Space
//
// The class IntVect is an implementation of an integer vector in a
// SPACEDIM-dimensional space.  It represents a point in a discrete space.
// IntVect values are accessed using the operator[] function, as for a normal
// C++ array.  In addition, the basic arithmetic operators have been overloaded
// to implement scaling and translation operations.
//
	
class IntVect
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<IntVect>::type();
public:
    //
    // An IntVect comparison class.  Needed for STL associative
    // containers indexed by IntVects.  Uses lexicographical comparison.
    //
    class Compare
    {
    public:
        bool operator () (const IntVect& lhs,
                          const IntVect& rhs) const
        {
            return lhs.lexLT(rhs);
        }
    };
    //
    // Construct an IntVect whose components are all zero.
    //
    IntVect () { D_EXPR(vect[0] = 0, vect[1] = 0, vect[2] = 0); }
    //
    // Construct an IntVect given the specific values for its
    // coordinates.  D_DECL is a macro that sets the constructor
    // to take BL_SPACEDIM arguments.
    //
    IntVect (D_DECL(int i, int j, int k)) { D_EXPR(vect[0] = i, vect[1] = j, vect[2] = k); }
    //
    // Construct an IntVect setting the coordinates to the
    // corresponding values in the integer array a.
    //
    explicit IntVect (const int* a);
    //
    // Construct an IntVect from an Array<int>.  It is an error if
    // the Array<int> doesn' t have the same dimension as this
    // IntVect.
    //
    explicit IntVect (const Array<int>& a);
    //
    // Returns a reference to the i'th coordinate of the IntVect.
    //
    int& operator[] (int i) { BL_ASSERT(i>=0 && i < BL_SPACEDIM); return vect[i]; }
    //
    // Returns the i'th coordinate of the IntVect.
    //
    int operator[] (int i) const { BL_ASSERT(i>=0 && i < BL_SPACEDIM); return vect[i]; }
    //
    // Set i'th coordinate of IntVect to val.
    //
    IntVect& setVal (int i, int val)
    {
        BL_ASSERT(i >=0 && i < BL_SPACEDIM); vect[i] = val; return *this;
    }
    //
    // Returns a const pointer to an array of coordinates of the
    // IntVect.  Useful for arguments to FORTRAN calls.
    //
    const int* getVect () const& { return vect; }
    int*       getVect () &      { return vect; }
    int*       getVect () && = delete;
    //
    // Returns true of IntVect is equivalent to argument IntVect.
    // All comparisons between analogous components must be satisfied.
    //
    bool operator== (const IntVect& p) const
    {
        return D_TERM(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2]);
    }
    //
    // Returns true of IntVect is different from argument IntVect.
    // All comparisons between analogous components must be satisfied.
    //
    bool operator!= (const IntVect& p) const
    {
        return D_TERM(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2]);
    }
    //
    // Returns true if IntVect is less than argument IntVect.
    // All comparisons between analogous components must be satisfied.
    //
    // NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    //
    bool operator< (const IntVect& p) const
    {
        return D_TERM(vect[0] < p[0], && vect[1] < p[1], && vect[2] < p[2]);
    }
    //
    // Returns true if IntVect is less than or equal to argument
    // IntVect.  All comparisons between analogous components must
    // be satisfied.
    //
    // NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    //
    bool operator<= (const IntVect& p) const
    {
        return D_TERM(vect[0] <= p[0], && vect[1] <= p[1], && vect[2] <= p[2]);
    }
    //
    // Returns true if IntVect is greater than argument IntVect.
    // All comparisons between analogous components must be satisfied.
    //
    // NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    //
    bool operator> (const IntVect& p) const
    {
        return D_TERM(vect[0] > p[0], && vect[1] > p[1], && vect[2] > p[2]);
    }
    //
    // Returns true if IntVect is greater than or equal to argument
    // IntVect.  All comparisons between analogous components must
    // be satisfied.
    //
    // NOTE: This is NOT a strict weak ordering usable by STL sorting algorithms.
    //
    bool operator>= (const IntVect& p) const
    {
        return D_TERM(vect[0] >= p[0], && vect[1] >= p[1], && vect[2] >= p[2]);
    }
    //
    // Returns true if this IntVect is lexically less than the argument.
    // This is a strict weak ordering usable by STL sorting algorithms.
    //
    bool lexLT (const IntVect& s) const;
    //
    // Returns true if this IntVect is lexically greater than the argument.
    // This is a strict weak ordering usable by STL sorting algorithms.
    //
    bool lexGT (const IntVect& s) const;
    //
    // Unary plus -- for completeness.
    //
    const IntVect operator+ () const { return *this; }
    //
    // Unary Minus -- negates all components.
    //
    IntVect operator- () const { return IntVect(D_DECL(-vect[0], -vect[1], -vect[2] )); }
    //
    // Modifies IntVect by addition of a scalar to each component.
    //
    IntVect& operator+= (int s)
    {
        D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s); return *this;
    }
    //
    // Modifies IntVect by component-wise addition with argument.
    //
    IntVect& operator+= (const IntVect& p)
    {
        D_EXPR(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2]); return *this;
    }
    //
    // Modifies IntVect by multiplication of a scalar to each component.
    //
    IntVect& operator*= (int s)
    {
        D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s); return *this;
    }
    //
    // Modifies IntVect by component-wise multiplication with argument.
    //
    IntVect& operator*= (const IntVect& p)
    {
        D_EXPR(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2]); return *this;
    }
    //
    // Modifies IntVect by division by a scalar to each component.
    //
    IntVect& operator/= (int s)
    {
        D_EXPR(vect[0] /= s, vect[1] /= s, vect[2] /= s); return *this;
    }
    //
    // Modifies IntVect by component-wise division with argument.
    //
    IntVect& operator/= (const IntVect& p)
    {
        D_EXPR(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2]); return *this;
    }
    //
    // Modifies IntVect by subtraction of a scalar to each component.
    //
    IntVect& operator-= (int s)
    {
        D_EXPR(vect[0] -= s, vect[1] -= s, vect[2] -= s); return *this;
    }
    //
    // Modifies IntVect by component-wise subtraction with argument.
    //
    IntVect& operator-= (const IntVect& p)
    {
        D_EXPR(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2]); return *this;
    }
    //
    // Returns component-wise sum of IntVect and argument.
    //
    IntVect operator+ (const IntVect& p) const
    {
        return IntVect(D_DECL(vect[0] + p[0], vect[1] + p[1], vect[2] + p[2]));
    }
    //
    // Return an IntVect that is this IntVect + s.
    //
    IntVect operator+ (int s) const
    {
        return IntVect(D_DECL(vect[0] + s, vect[1] + s, vect[2] + s));
    }
    //
    // Returns component-wise difference of IntVect and argument.
    //
    IntVect operator- (const IntVect& p) const
    {
        return IntVect(D_DECL(vect[0] - p[0], vect[1] - p[1], vect[2] - p[2]));
    }
    //
    // Return an IntVect that is this IntVect - s.
    //
    IntVect operator- (int s) const
    {
        return IntVect(D_DECL(vect[0] - s, vect[1] - s, vect[2] - s));
    }
    //
    // Returns component-wise product of IntVect and argument.
    //
    IntVect operator* (const IntVect& p) const
    {
        return IntVect(D_DECL(vect[0] * p[0], vect[1] * p[1], vect[2] * p[2]));
    }
    //
    // Returns component-wise product of IntVect and s.
    //
    IntVect operator* (int s) const
    {
        return IntVect(D_DECL(vect[0] * s, vect[1] * s, vect[2] * s));
    }
    //
    // Returns component-wise division of IntVect by argument.
    //
    IntVect operator/ (const IntVect& p) const
    {
        return IntVect(D_DECL(vect[0] / p[0], vect[1] / p[1], vect[2] / p[2]));
    }
    //
    // Returns component-wise division of IntVect by s.
    //
    IntVect operator/ (int s) const
    {
        return IntVect(D_DECL(vect[0] / s, vect[1] / s, vect[2] / s));
    }
    //
    // Modifies IntVect by taking component-wise min with argument.
    //
    IntVect& min (const IntVect& p)
    {
        D_EXPR(vect[0] = std::min(vect[0], p.vect[0]),
               vect[1] = std::min(vect[1], p.vect[1]),
               vect[2] = std::min(vect[2], p.vect[2]));
        return *this;
    }
    //
    // Modifies IntVect by taking component-wise max with argument.
    //
    IntVect& max (const IntVect& p)
    {
        D_EXPR(vect[0] = std::max(vect[0], p.vect[0]),
               vect[1] = std::max(vect[1], p.vect[1]),
               vect[2] = std::max(vect[2], p.vect[2]));
        return *this;
    }
    //
    // Modify IntVect by multiplying each coordinate by s.
    //
    IntVect& scale (int s) { D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s); return *this; }
    //
    // Modify IntVect by reflecting it in the plane defined by
    // the index ref_ix and with normal in the direction of idir.
    //
    IntVect& reflect (int ref_ix, int idir)
    {
        BL_ASSERT(idir >= 0 && idir < BL_SPACEDIM);
        vect[idir] = -vect[idir] + 2*ref_ix;
        return *this;
    }
    //
    // Modify IntVect by adding s to given coordinate.
    //
    IntVect& shift (int coord, int s)
    {
        BL_ASSERT(coord >= 0 && coord < BL_SPACEDIM); vect[coord] += s; return *this;
    }
    //
    // Equivalent to shift(0,iv[0]).shift(1,iv[1]) ...
    //
    IntVect& shift (const IntVect& iv) { *this += iv; return *this; }
    //
    // Modify IntVect by adding s to each coordinate.
    //
    IntVect& diagShift (int s)
    {
        D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s); return *this;
    }
    //
    // Modify IntVect by component-wise integer projection.
    //
    IntVect& coarsen (const IntVect& p);
    //
    // Modify IntVect by component-wise integer projection.
    //
    IntVect& coarsen (int p);
    //
    // This static member function returns a reference to a constant IntVect
    // object, all of whose BL_SPACEDIM arguments are set to zero (0).
    // Figuratively, it is the zero vector in BL_SPACEDIM-dimensional space.
    // It is provided as a convenient way to specify the zero vector.
    //
    static const IntVect& TheZeroVector ();
    //
    // This static member function returns a reference to a constant IntVect
    // object, all of whose BL_SPACEDIM arguments are set to one (1).
    // Figuratively, it is the unit vector in BL_SPACEDIM-dimensional space.
    // It is provided as a convenient way to specify the unit vector.
    //
    static const IntVect& TheUnitVector ();
    //
    // This static member function returns a reference to a constant IntVect
    // object, all of whose BL_SPACEDIM arguments are set to zero except that
    // the d-direction is set to one. 
    //
    static const IntVect& TheDimensionVector (int d);
    //
    // This static member function returns a reference to a constant IntVect
    // object, all of whose BL_SPACEDIM arguments are set to IndexType::NODE.
    // It is provided as a convenience to our users when defining Boxes.
    //
    static const IntVect& TheNodeVector ();
    // 
    // This static member function returns a reference to a constant IntVect
    // object, all of whose BL_SPACEDIM arguments are set to IndexType::CELL.
    // It is provided as a convenience to our users when defining Boxes.
    //
    static const IntVect& TheCellVector ();

    static const IntVect& TheMaxVector ();
    static const IntVect& TheMinVector ();

private:

    int vect[BL_SPACEDIM];
};

namespace BoxLib
{
    //
    // Returns a basis vector in the given coordinate direction;
    // eg.  IntVect BASISV(1) == (0,1,0).  Note that the coordinate
    // directions are zero based.
    //
    IntVect BASISV (int dir);
    //
    // Returns the IntVect that is the component-wise minimum of two
    // argument IntVects.
    //
    IntVect min (const IntVect& p1, const IntVect& p2);
    //
    // Returns the IntVect that is the component-wise maximum of two
    // argument IntVects.
    //
    IntVect max (const IntVect& p1, const IntVect& p2);
    //
    // Returns a IntVect obtained by multiplying each of the
    // components of this IntVect by s.
    //
    IntVect scale (const IntVect& p, int s);
    //
    // Returns an IntVect that is the reflection of input in the
    // plane which passes through ref_ix and normal to the
    // coordinate direction idir.
    //
    IntVect reflect(const IntVect& a, int ref_ix, int idir);
    //
    // Returns IntVect obtained by adding s to each of the
    // components of this IntVect.
    //
    IntVect diagShift (const IntVect& p, int s);
    //
    // Returns an IntVect that is the component-wise integer
    // projection of p by s.
    //
    IntVect coarsen (const IntVect& p, int s);
    //
    // Returns an IntVect which is the component-wise integer
    // projection of IntVect p1 by IntVect p2.
    //
    IntVect coarsen (const IntVect& p1, const IntVect& p2);
}
//
// Returns p + s.
//
inline
IntVect operator+ (int s, const IntVect& p)
{
    return IntVect(D_DECL(p[0] + s, p[1] + s, p[2] + s));
}
//
// Returns -p + s.
//
inline
IntVect operator- (int s, const IntVect& p)
{
    return IntVect(D_DECL(s - p[0], s - p[1], s - p[2]));
}
//
// Returns p * s.
//
inline
IntVect operator* (int s, const IntVect& p)
{
    return IntVect(D_DECL(s * p[0], s * p[1], s * p[2]));
}
//
// Print the IntVect to given output stream in ASCII.
//
std::ostream& operator<< (std::ostream& os, const IntVect& iv);
//
// Read next IntVect from given input stream.
//
std::istream& operator>> (std::istream& is, IntVect& iv);

#endif /*BL_INTVECT_H*/
