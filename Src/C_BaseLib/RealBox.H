//BL_COPYRIGHT_NOTICE

#ifndef _RealBox_H_
#define _RealBox_H_ 

//
// $Id: RealBox.H,v 1.9 1999-05-10 18:54:09 car Exp $
//

#include <Array.H>
#include <REAL.H>
#include <Box.H>

//
//@Man:
//@Memo: A Real Box
/*@Doc:

  A Box with real dimensions.  A RealBox is OK iff volume > 0.
*/

class RealBox
{
public:
    //
    //@ManDoc: Nice ASCII output.
    //
    friend ostream& operator<< (ostream&, const RealBox& );
    //
    //@ManDoc: Nice ASCII input.
    //
    friend istream& operator>> (istream&, RealBox& );
    //
    //@ManDoc: The default constructor.  Builds invalid RealBox.
    //
    RealBox ();
    //
    //@ManDoc: Construct region given diagonal points.
    //
    RealBox (const Real* lo,
             const Real* hi);

    /*@ManDoc: Construct region given index box, cell spacing
               and physical location of index (0,0,0).
    */
    RealBox (const Box&  bx,
             const Real* dx,
             const Real* base);
    //
    //@ManDoc: Explicit dimension specific constructors.
    //
    RealBox (D_DECL(Real x0, Real y0, Real z0),
             D_DECL(Real x1, Real y1, Real z1));
    //
    //@ManDoc: Set epsilon used by contains() function when comparing floating point numbers.
    //
    static void setEpsilon (Real epsilon);
    //
    //@ManDoc: Returns epsilon.
    //
    static Real epsilon ();
    //
    //@ManDoc: Returns lo side.
    //
    const Real* lo () const;
    //
    //@ManDoc: Returns hide side.
    //
    const Real* hi () const;
    //
    //@ManDoc: Returns length in each direction.
    //
    const Real* length () const;
    //
    //@ManDoc: Returns length in specified direction.
    //
    Real lo (int dir) const;
    //
    //@ManDoc: Returns hi side in specified direction.
    //
    Real hi (int dir) const;
    //
    //@ManDoc: Returns length in specified direction.
    //
    Real length (int dir) const;
    //
    //@ManDoc: Sets lo side.
    //
    void setLo (const Real* lo);
    //
    //@ManDoc: Sets lo side.
    //
    void setLo (const Array<Real>& lo);
    //
    //@ManDoc: Sets lo side in specified direction.
    //
    void setLo (int dir, Real lo);
    //
    //@ManDoc: Sets hi side.
    //
    void setHi (const Real* hi);
    //
    //@ManDoc: Sets hi side.
    //
    void setHi (const Array<Real>& hi);
    //
    //@ManDoc: Sets hi side in specified direction.
    //
    void setHi (int dir, Real hi);
    //
    //@ManDoc: Is the RealBox OK; i.e. does it have positive volume?
    //
    bool ok () const;
    //
    //@ManDoc: Is the specified point contained in the RealBox?
    //
    bool contains (const Real* point) const;
    //
    //@ManDoc: Is the specified RealBox contained in this RealBox?
    //
    bool contains (const RealBox& bx) const;

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
    //
    // Static data.
    //
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
   BL_ASSERT(indx >= 0 && indx < BL_SPACEDIM);
   xlo[indx] = lo;
   computeBoxLen();
}

inline
void
RealBox::setHi (int  indx,
                Real hi)
{
    BL_ASSERT(indx >= 0 && indx < BL_SPACEDIM);
    xhi[indx] = hi;
    computeBoxLen();
}

inline
bool
RealBox::contains (const RealBox& rb) const
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

#endif /*_RealBox_H_*/
