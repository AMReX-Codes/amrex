//BL_COPYRIGHT_NOTICE

#ifndef _RealBox_H_
#define _RealBox_H_ 

//
// $Id: RealBox.H,v 1.3 1997-12-11 05:01:10 lijewski Exp $
//

#include <Array.H>
#include <REAL.H>
#include <Box.H>

// rules:
//  ok if and only if (volume > 0)
//  is zero volume intersection valid?
//

class RealBox
{
public:

#ifdef __GNUG__
     bool operator== (const RealBox& rhs) const { return true; }
#endif 

    friend ostream& operator<< (ostream&, const RealBox& );
    friend istream& operator>> (istream&, RealBox& );
    //
    // Construct invalid region.
    //
    RealBox ();
    //
    // Construct region given diagional points.
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

    static void setEpsilon (Real epsilon);
    static Real epsilon ();

    const Real* lo () const;
    const Real* hi () const;

    const Real* length () const;

    Real lo (int dir) const;
    Real hi (int dir) const;

    Real length (int dir) const;

    void setLo (const Real* lo);
    void setLo (const Array<Real>& lo);
    void setLo (int indx, Real lo);

    void setHi (const Real* hi);
    void setHi (const Array<Real>& hi);
    void setHi (int indx, Real hi);

    bool ok () const;
    bool contains (const Real* point);
    bool contains (const RealBox &bx);

protected:
    //
    // Helper function.
    //
    void computeBoxLen ();
    //
    // The data.
    //
    Real xlo[BL_SPACEDIM];
    Real xhi[BL_SPACEDIM];
    Real len[BL_SPACEDIM];

    static Real eps;
};

//
// Inlines.
//

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
const Real*
RealBox::length () const
{
    return len;
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
    return len[dir];
}

inline
void
RealBox::computeBoxLen ()
{
    D_EXPR(len[0] = xhi[0]-xlo[0],
           len[1] = xhi[1]-xlo[1],
           len[2] = xhi[2]-xlo[2]);
}

inline
RealBox::RealBox ()
{
    D_TERM(xlo[0] , = xlo[1] , = xlo[2] ) = 0.;
    D_TERM(xhi[0] , = xhi[1] , = xhi[2] ) = -1.;
    computeBoxLen();
}

inline
RealBox::RealBox (const Real* lo,
                  const Real* hi)
{
    D_EXPR(xlo[0] = lo[0] , xlo[1] = lo[1] , xlo[2] = lo[2]);
    D_EXPR(xhi[0] = hi[0] , xhi[1] = hi[1] , xhi[2] = hi[2]);
    computeBoxLen() ;
}

inline
RealBox::RealBox (D_DECL(Real x0, Real y0, Real z0),
                  D_DECL(Real x1, Real y1, Real z1))
{
    D_EXPR(xlo[0] = x0 , xlo[1] = y0 , xlo[2] = z0);
    D_EXPR(xhi[0] = x1 , xhi[1] = y1 , xhi[2] = z1);
    computeBoxLen() ;
}

inline
void
RealBox::setLo (const Real* lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
    computeBoxLen();
}

inline
void
RealBox::setLo (const Array<Real> &lo)
{
    D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]);
    computeBoxLen();
}

inline
void
RealBox::setHi (const Real* hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
    computeBoxLen();
}

inline
void
RealBox::setHi (const Array<Real>& hi)
{
    D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]);
    computeBoxLen();
}

inline
void
RealBox::setLo (int  indx,
                Real lo)
{
   assert(indx >= 0 && indx < BL_SPACEDIM);
   xlo[indx] = lo;
   computeBoxLen();
}

inline
void
RealBox::setHi (int  indx,
                Real hi)
{
    assert(indx >= 0 && indx < BL_SPACEDIM);
    xhi[indx] = hi;
    computeBoxLen();
}

inline
bool
RealBox::contains (const RealBox& rb)
{
    return contains(rb.xlo) && contains(rb.xhi);
}

inline
bool
RealBox::ok () const
{
    return (len[0] > eps)
#if (BL_SPACEDIM > 1)
	&& (len[1] > eps)
#endif   
#if (BL_SPACEDIM > 2)
	&& (len[2] > eps)
#endif
   ;
}

#endif
