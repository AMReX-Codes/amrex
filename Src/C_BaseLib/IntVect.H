#ifndef BL_INTVECT_H
#define BL_INTVECT_H
//
// $Id: IntVect.H,v 1.21 2001-10-10 20:12:44 car Exp $
//
#include <iosfwd>

#include <ccse-mpi.H>
#include <SPACE.H>
#include <Array.H>
#include <BLassert.H>

//
//@Man:
//@Memo: An Integer Vector in SPACEDIM-Dimensional Space
/*@Doc:

  The class IntVect is an implementation of an integer vector in a
  SPACEDIM-dimensional space.  It represents a point in a discrete space.
  IntVect values are accessed using the operator[] function, as for a normal
  C++ array.  In addition, the basic arithmetic operators have been overloaded
  to implement scaling and translation operations.
*/

class IntVect;

namespace BoxLib
{
    /*@ManDoc: Returns a basis vector in the given coordinate direction;
      eg.  IntVect3d BASISV(1) == (0,1,0).  Note that the coordinate
      directions are zero based.
    */
    IntVect BASISV (int dir);

    /*@ManDoc: Returns the IntVect that is the component-wise minimum of two
               argument IntVects.
    */
    IntVect min (const IntVect& p1, const IntVect& p2);

    /*@ManDoc: Returns the IntVect that is the component-wise minimum of two
      argument IntVects.
    */
    IntVect max (const IntVect& p1, const IntVect& p2);

    /*@ManDoc: Returns a IntVect obtained by multiplying each of the
      components of this IntVect by s.
    */
    IntVect scale (const IntVect& p, int s);

    /*@ManDoc: Returns an IntVect that is the reflection of input in the
      plane which passes through ref\_ix and normal to the
      coordinate direction idir.
    */
    IntVect reflect(const IntVect& a, int ref_ix, int idir);

    /*@ManDoc: Returns IntVect obtained by adding s to each of the
      components of this IntVect.
    */
    IntVect diagShift (const IntVect& p, int s);

    /*@ManDoc: Returns an IntVect that is the component-wise integer
      projection of p by s.
    */
    IntVect coarsen (const IntVect& p, int s);

    /*@ManDoc: Returns an IntVect which is the component-wise integer
      projection of IntVect p1 by IntVect p2.
    */
    IntVect coarsen (const IntVect& p1, const IntVect& p2);
}
//
//@ManDoc: Returns p + s.
//
IntVect operator+ (int s, const IntVect& p);
//
//@ManDoc: Returns -p + s.
//
IntVect operator- (int s, const IntVect& p);
//
//@ManDoc: Returns p * s.
//
IntVect operator* (int s, const IntVect& p);
//
//@ManDoc: Print the IntVect to given output stream in ASCII.
//
std::ostream& operator<< (std::ostream&  os, const IntVect& iv);
//
//@ManDoc: Read next IntVect from given input stream.
//
std::istream& operator>> (std::istream& os, IntVect& iv);
	

class IntVect
{
    friend MPI_Datatype ParallelDescriptor::Mpi_typemap<IntVect>::type();
public:
    //
    //@ManDoc: Construct an IntVect whose components are all zero.
    //
    IntVect ();

    /*@ManDoc: Construct an IntVect given the specific values for its
               coordinates.  D\_DECL is a macro that sets the constructor
               to take BL\_SPACEDIM arguments.
    */
    IntVect (D_DECL(int i, int j, int k));

    /*@ManDoc: Construct an IntVect setting the coordinates to the
               corresponding values in the integer array a.
    */
    explicit IntVect (const int* a);

    /*@ManDoc: Construct an IntVect from an Array<int>.  It is an error if
               the Array<int> doesn' t have the same dimension as this
               IntVect.
    */
    explicit IntVect (const Array<int>& a);
    //
    //@ManDoc: The copy constructor.
    //
    IntVect (const IntVect& rhs);
    //
    //@ManDoc: The assignment operator.
    //
    IntVect& operator= (const IntVect& rhs);
    //
    //@ManDoc: Returns a reference to the i'th coordinate of the IntVect.
    //
    int& operator[] (int i);
    //
    //@ManDoc: Returns the i'th coordinate of the IntVect.
    //
    int operator[] (int i) const;
    //
    //@ManDoc: Set i'th coordinate of IntVect to val.
    //
    IntVect& setVal (int i,
                     int val);

    /*@ManDoc: Returns a const pointer to an array of coordinates of the
               IntVect.  Useful for arguments to FORTRAN calls.
    */
    const int* getVect () const;

    /*@ManDoc: Returns true of IntVect is equivalent to argument IntVect.
               All comparisons between analogous components must be satisfied.
    */
    bool operator== (const IntVect& p) const;

    /*@ManDoc: Returns true of IntVect is different from argument IntVect.
               All comparisons between analogous components must be satisfied.
    */
    bool operator!= (const IntVect& p) const;

    /*@ManDoc: Returns true if IntVect is less than argument IntVect.
               All comparisons between analogous components must be satisfied.
    */
    bool operator< (const IntVect& p) const;

    /*@ManDoc: Returns true if IntVect is less than or equal to argument
               IntVect.  All comparisons between analogous components must
               be satisfied.
    */
    bool operator<= (const IntVect& p) const;

    /*@ManDoc: Returns true if IntVect is greater than argument IntVect.
               All comparisons between analogous components must be satisfied.
    */
    bool operator> (const IntVect& p) const;

    /*@ManDoc: Returns true if IntVect is greater than or equal to argument
               IntVect.  All comparisons between analogous components must
               be satisfied.
    */
    bool operator>= (const IntVect& p) const;

    /*@ManDoc: Returns true if this IntVect is lexically less than the
               argument.
    */
    bool lexLT (const IntVect& s) const;

    /*@ManDoc: Returns true if this IntVect is lexically greater than the
               argument.
    */
    bool lexGT (const IntVect& s) const;
    //
    //@ManDoc: Unary plus -- for completeness.
    //
    IntVect operator+ () const;
    //
    //@ManDoc: Unary Minus -- negates all components.
    //
    IntVect operator- () const;
    //
    //@ManDoc: Modifies IntVect by addition of a scaler to each component.
    //
    IntVect& operator+= (int s);
    //
    //@ManDoc: Modifies IntVect by component-wise addition with argument.
    //
    IntVect& operator+= (const IntVect& p);

    /*@ManDoc: Modifies IntVect by multiplication of a scaler to each
               component.
    */
    IntVect& operator*= (int s);
    //
    //@ManDoc: Modifies IntVect by component-wise multiplication with argument.
    //
    IntVect& operator*= (const IntVect& p);
    //
    //@ManDoc: Modifies IntVect by division by a scaler to each component.
    //
    IntVect& operator/= (int s);
    //
    //@ManDoc: Modifies IntVect by component-wise division with argument.
    //
    IntVect& operator/= (const IntVect& p);
    //
    //@ManDoc: Modifies IntVect by subtraction of a scaler to each component.
    //
    IntVect& operator-= (int s);
    //
    //@ManDoc: Modifies IntVect by component-wise subtraction with argument.
    //
    IntVect& operator-= (const IntVect& p);
    //
    //@ManDoc: Returns component-wise sum of IntVect and argument.
    //
    IntVect operator+ (const IntVect& p) const;
    //
    //@ManDoc: Return an IntVect that is this IntVect + s.
    //
    IntVect operator+ (int s) const;
    //
    //@ManDoc: Returns component-wise difference of IntVect and argument.
    //
    IntVect operator- (const IntVect& p) const;
    //
    //@ManDoc: Return an IntVect that is this IntVect - s.
    //
    IntVect operator- (int s) const;
    //
    //@ManDoc: Returns component-wise product of IntVect and argument.
    //
    IntVect operator* (const IntVect& p) const;
    //
    //@ManDoc: Returns component-wise product of IntVect and s.
    //
    IntVect operator* (int s) const;
    //
    //@ManDoc: Returns component-wise division of IntVect by argument.
    //
    IntVect operator/ (const IntVect& p) const;
    //
    //@ManDoc: Returns component-wise division of IntVect by s.
    //
    IntVect operator/ (int s) const;
    //
    //@ManDoc: Modifies IntVect by taking component-wise min with argument.
    //
    IntVect& min (const IntVect& p);
    //
    //@ManDoc: Modifies IntVect by taking component-wise max with argument.
    //
    IntVect& max (const IntVect& p);
    //
    //@ManDoc: Modify IntVect by multiplying each coordinate by s.
    //
    IntVect& scale (int s);

    /*@ManDoc: Modify IntVect by reflecting it in the plane defined by
               the index ref\_ix and with normal in the direction of idir.
    */
    IntVect& reflect (int ref_ix,
                      int idir);
    //
    //@ManDoc: Modify IntVect by adding s to given coordinate.
    //
    IntVect& shift (int coord,
                    int s);
    //
    //@ManDoc: Equivalent to shift(0,iv[0]).shift(1,iv[1]) ...
    //
    IntVect& shift (const IntVect& iv);
    //
    //@ManDoc: Modify IntVect by adding s to each coordinate.
    //
    IntVect& diagShift (int s);
    //
    //@ManDoc: Modify IntVect by component-wise integer projection.
    //
    IntVect& coarsen (const IntVect& p);
    //
    //@ManDoc: Modify IntVect by component-wise integer projection.
    //
    IntVect& coarsen (int p);

    /*@ManDoc:
      This static member function returns a reference to a constant IntVect
      object, all of whose BL\_SPACEDIM arguments are set to zero (0).
      Figuratively, it is the zero vector in BL\_SPACEDIM-dimensional space.
      It is provided as a convenient way to specify the zero vector.
      */
    static const IntVect& TheZeroVector ();

    /*@ManDoc:
      This static member function returns a reference to a constant IntVect
      object, all of whose BL\_SPACEDIM arguments are set to one (1).
      Figuratively, it is the unit vector in BL\_SPACEDIM-dimensional space.
      It is provided as a convenient way to specify the unit vector.
      */
    static const IntVect& TheUnitVector ();

    /*@ManDoc:
      This static member function returns a reference to a constant IntVect
      object, all of whose BL\_SPACEDIM arguments are set to IndexType::NODE.
      It is provided as a convenience to our users when defining Boxes.
      */
    static const IntVect& TheNodeVector ();

    /*@ManDoc: 
      This static member function returns a reference to a constant IntVect
      object, all of whose BL\_SPACEDIM arguments are set to IndexType::CELL.
      It is provided as a convenience to our users when defining Boxes.
      */
    static const IntVect& TheCellVector ();

protected:

    int vect[BL_SPACEDIM];
};

inline
IntVect::IntVect ()
{
    D_EXPR(vect[0] = 0, vect[1] = 0, vect[2] = 0);
}

inline
IntVect::IntVect (D_DECL(int i, int j, int k))
{
    D_EXPR(vect[0] = i, vect[1] = j, vect[2] = k);
}

inline
IntVect::IntVect (const IntVect &iv)
{
    D_EXPR(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2]);
}

inline
IntVect&
IntVect::operator= (const IntVect &iv)
{
    D_EXPR(vect[0]=iv.vect[0], vect[1]=iv.vect[1], vect[2]=iv.vect[2]);
    return *this;
}

inline
int&
IntVect::operator[] (int i)
{
    BL_ASSERT(i>=0 && i < BL_SPACEDIM);
    return vect[i];
}

inline
int
IntVect::operator[] (int i) const
{
    BL_ASSERT(i>=0 && i < BL_SPACEDIM);
    return vect[i];
}

inline
const int*
IntVect::getVect () const
{
    return vect;
}

inline
IntVect&
IntVect::setVal (int i,
                 int val)
{
    BL_ASSERT(i >=0 && i < BL_SPACEDIM);
    vect[i] = val;
    return *this;
}

inline
bool
IntVect::operator== (const IntVect& p) const
{
    return D_TERM(vect[0] == p[0], && vect[1] == p[1], && vect[2] == p[2]);
}

inline
bool
IntVect::operator!= (const IntVect& p) const
{
    return D_TERM(vect[0] != p[0], || vect[1] != p[1], || vect[2] != p[2]);
}

inline
bool
IntVect::operator< (const IntVect& p) const
{
    return D_TERM(vect[0] < p[0], && vect[1] < p[1], && vect[2] < p[2]);
}

inline
bool
IntVect::operator<= (const IntVect& p) const
{
    return D_TERM(vect[0] <= p[0], && vect[1] <= p[1], && vect[2] <= p[2]);
}

inline
bool
IntVect::operator> (const IntVect& p) const
{
    return D_TERM(vect[0] > p[0], && vect[1] > p[1], && vect[2] > p[2]);
}

inline
bool
IntVect::operator>= (const IntVect& p) const
{
    return D_TERM(vect[0] >= p[0], && vect[1] >= p[1], && vect[2] >= p[2]);
}

inline
IntVect
IntVect::operator+ () const
{
    return *this;
}

inline
IntVect
IntVect::operator- () const
{
    return IntVect(D_DECL(-vect[0], -vect[1], -vect[2] ));
}

inline
IntVect&
IntVect::operator+= (int s)
{
    D_EXPR(vect[0] += s, vect[1] += s, vect[2] += s);
    return *this;
}

inline
IntVect&
IntVect::operator+= (const IntVect& p)
{
    D_EXPR(vect[0] += p[0], vect[1] += p[1], vect[2] += p[2]);
    return *this;
}

inline
IntVect&
IntVect::operator*= (int s)
{
    D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
    return *this;
}

inline
IntVect&
IntVect::operator*= (const IntVect &p)
{
    D_EXPR(vect[0] *= p[0], vect[1] *= p[1], vect[2] *= p[2]);
    return *this;
}

inline
IntVect&
IntVect::operator/= (int s)
{
    D_EXPR(vect[0] /= s, vect[1] /= s, vect[2] /= s);
    return *this;
}

inline
IntVect&
IntVect::operator/= (const IntVect& p)
{
    D_EXPR(vect[0] /= p[0], vect[1] /= p[1], vect[2] /= p[2]);
    return *this;
}

inline
IntVect&
IntVect::operator-= (int s)
{
    D_EXPR(vect[0] -= s, vect[1] -= s, vect[2] -= s);
    return *this;
}

inline
IntVect&
IntVect::operator-= (const IntVect& p)
{
    D_EXPR(vect[0] -= p[0], vect[1] -= p[1], vect[2] -= p[2]);
    return *this;
}

inline
IntVect
IntVect::operator+ (const IntVect& p) const
{
    return IntVect(D_DECL(vect[0] + p[0], vect[1] + p[1], vect[2] + p[2]));
}

inline
IntVect
IntVect::operator+ (int s) const
{
    return IntVect(D_DECL(vect[0] + s, vect[1] + s, vect[2] + s));
}

inline
IntVect
IntVect::operator- (const IntVect& p) const
{
    return IntVect(D_DECL(vect[0] - p[0], vect[1] - p[1], vect[2] - p[2]));
}

inline
IntVect
IntVect::operator- (int s) const
{
    return IntVect(D_DECL(vect[0] - s, vect[1] - s, vect[2] - s));
}

inline
IntVect
IntVect::operator* (const IntVect& p) const
{
    return IntVect(D_DECL(vect[0] * p[0], vect[1] * p[1], vect[2] * p[2]));
}

inline
IntVect
IntVect::operator* (int s) const
{
    return IntVect(D_DECL(vect[0] * s, vect[1] * s, vect[2] * s));
}

inline
IntVect
IntVect::operator/ (const IntVect& p) const
{
    return IntVect(D_DECL(vect[0] / p[0], vect[1] / p[1], vect[2] / p[2]));
}

inline
IntVect
IntVect::operator/ (int s) const
{
    return IntVect(D_DECL(vect[0] / s, vect[1] / s, vect[2] / s));
}

#endif /*BL_INTVECT_H*/
