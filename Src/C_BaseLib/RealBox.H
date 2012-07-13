
#ifndef _RealBox_H_
#define _RealBox_H_ 

#include <iosfwd>

#include <Array.H>
#include <REAL.H>
#include <Box.H>

//
// A Real Box
//
//  A Box with real dimensions.  A RealBox is OK iff volume > 0.
//

class RealBox
{
public:
    //
    // The default constructor.  Builds invalid RealBox.
    //
    RealBox ();
    //
    // Construct region given diagonal points.
    //
    RealBox (const Real* lo,
             const Real* hi);
    //
    // Construct region given index box, cell spacing
    // and physical location of index (0,0,0).
    //
    RealBox (const Box&  bx,
             const Real* dx,
             const Real* base);
    //
    // Explicit dimension specific constructors.
    //
    RealBox (D_DECL(Real x0, Real y0, Real z0),
             D_DECL(Real x1, Real y1, Real z1));
    //
    // Set epsilon used by contains() function when comparing floating point numbers.
    //
    static void setEpsilon (Real epsilon);
    //
    // Returns epsilon.
    //
    static Real epsilon ();
    //
    // Returns lo side.
    //
    const Real* lo () const;
    //
    // Returns hide side.
    //
    const Real* hi () const;
    //
    // Returns length in specified direction.
    //
    Real lo (int dir) const;
    //
    // Returns hi side in specified direction.
    //
    Real hi (int dir) const;
    //
    // Returns length in specified direction.
    //
    Real length (int dir) const;
    //
    // Sets lo side.
    //
    void setLo (const Real* lo);
    //
    // Sets lo side.
    //
    void setLo (const Array<Real>& lo);
    //
    // Sets lo side in specified direction.
    //
    void setLo (int dir, Real lo);
    //
    // Sets hi side.
    //
    void setHi (const Real* hi);
    //
    // Sets hi side.
    //
    void setHi (const Array<Real>& hi);
    //
    // Sets hi side in specified direction.
    //
    void setHi (int dir, Real hi);
    //
    // Is the RealBox OK; i.e. does it have positive volume?
    //
    bool ok () const;
    //
    // Is the specified point contained in the RealBox?
    //
    bool contains (const Real* point) const;
    //
    // Is the specified RealBox contained in this RealBox?
    //
    bool contains (const RealBox& bx) const;

private:
    //
    // Static data.
    //
    static Real eps;
    //
    // The data.
    //
    Real xlo[BL_SPACEDIM];
    Real xhi[BL_SPACEDIM];
};
//
// Nice ASCII output.
//
std::ostream& operator<< (std::ostream&, const RealBox&);
//
// Nice ASCII input.
//
std::istream& operator>> (std::istream&, RealBox&);

inline
void
RealBox::setEpsilon (Real epsilon)
{
    eps = epsilon;
}

inline
Real
RealBox::epsilon ()
{
    return eps;
}

inline
const Real*
RealBox::lo () const
{
    return xlo;
}

inline
const Real*
RealBox::hi () const
{
    return xhi;
}

inline
Real
RealBox::lo (int dir) const
{
    return xlo[dir];
}

inline
Real
RealBox::hi (int dir) const
{
    return xhi[dir];
}

inline
Real
RealBox::length (int dir) const
{
    return xhi[dir]-xlo[dir];
}

inline
void
RealBox::setLo (const Real* lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
}

inline
void
RealBox::setLo (const Array<Real> &lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
}

inline
void
RealBox::setHi (const Real* hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
}

inline
void
RealBox::setHi (const Array<Real>& hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
}

inline
void
RealBox::setLo (int  indx,
                Real lo)
{
   BL_ASSERT(indx >= 0 && indx < BL_SPACEDIM);
   xlo[indx] = lo;
}

inline
void
RealBox::setHi (int  indx,
                Real hi)
{
    BL_ASSERT(indx >= 0 && indx < BL_SPACEDIM);
    xhi[indx] = hi;
}

#endif /*_RealBox_H_*/
