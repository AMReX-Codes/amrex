
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
    static void setEpsilon (Real epsilon) { eps = epsilon; }
    //
    // Returns epsilon.
    //
    static Real epsilon () { return eps; }
    //
    // Returns lo side.
    //
    const Real* lo () const& { return xlo; }
    const Real* lo () && = delete;
    //
    // Returns hide side.
    //
    const Real* hi () const& { return xhi; }
    const Real* hi () && = delete;
    //
    // Returns length in specified direction.
    //
    Real lo (int dir) const { return xlo[dir]; }
    //
    // Returns hi side in specified direction.
    //
    Real hi (int dir) const { return xhi[dir]; }
    //
    // Returns length in specified direction.
    //
    Real length (int dir) const { return xhi[dir]-xlo[dir]; }
    //
    // Sets lo side.
    //
    void setLo (const Real* lo) { D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]); }
    //
    // Sets lo side.
    //
    void setLo (const Array<Real>& lo) { D_EXPR(xlo[0] = lo[0], xlo[1] = lo[1], xlo[2] = lo[2]); }
    //
    // Sets lo side in specified direction.
    //
    void setLo (int dir, Real lo) { BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM); xlo[dir] = lo; }
    //
    // Sets hi side.
    //
    void setHi (const Real* hi) { D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]); }
    //
    // Sets hi side.
    //
    void setHi (const Array<Real>& hi) { D_EXPR(xhi[0] = hi[0], xhi[1] = hi[1], xhi[2] = hi[2]); }
    //
    // Sets hi side in specified direction.
    //
    void setHi (int dir, Real hi) { BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM); xhi[dir] = hi; }
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

#endif /*_RealBox_H_*/
