#ifndef _RealBox_H_
#define _RealBox_H_ 
//
// $Id: RealBox.H,v 1.11 2001-08-01 21:50:50 lijewski Exp $
//
#include <iosfwd>

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
//@ManDoc: Nice ASCII output.
//
std::ostream& operator<< (std::ostream&, const RealBox&);
//
//@ManDoc: Nice ASCII input.
//
std::istream& operator>> (std::istream&, RealBox&);

#endif /*_RealBox_H_*/
