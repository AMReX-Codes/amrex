//BL_COPYRIGHT_NOTICE

#ifndef _RealBox_H_
#define _RealBox_H_ 

//
// $Id: RealBox.H,v 1.1 1997-12-10 19:07:42 lijewski Exp $
//

#include <Array.H>
#include <REAL.H>
#include <Box.H>

// rules
// ok  if and only if (volume > 0)
// is zero volume intersection valid?

class RealBox
{
protected:
    Real xlo[BL_SPACEDIM], xhi[BL_SPACEDIM];
    Real len[BL_SPACEDIM];
    static Real eps;

public:
#ifdef __GNUG__
     bool operator== (const RealBox& rhs) const { return true; }
#endif 
      // -------------------------------------------
      // ---------  I/O friend functions  ----------
      // -------------------------------------------
    friend ostream& operator << (ostream&, const RealBox& );
    friend istream& operator >> (istream&, RealBox& );

      // -------------------------------------------
      // ------------  constructors  ---------------
      // -------------------------------------------

      // construct invalid region
    RealBox();

      // construct region given diagional points
    RealBox(const Real *lo, const Real* hi);
      
      // construct region given index box, cell spacing
      // and physical location of index (0,0,0)
    RealBox(const Box& bx, const Real* dx, const Real* base);

      // read region from input file
    RealBox(istream &is);

      // explicit dimension specific constructors
    RealBox(D_DECL(Real x0, Real y0, Real z0),
	    D_DECL(Real x1, Real y1, Real z1));

      // -------------------------------------------
      // -----------  static functions  ------------
      // -------------------------------------------
      // set epsilon for determining intersections
    static void setEpsilon(Real epsilon) { eps = epsilon; }
    static Real epsilon() { return eps; }

      // -------------------------------------------
      // ------------  access functions  -----------
      // -------------------------------------------
    const Real* lo() const {return xlo;}
    const Real* hi() const {return xhi;}
    const Real* length() const {return len;}
    Real  lo(int dir) const {return xlo[dir];}
    Real  hi(int dir) const {return xhi[dir];}
    Real  length(int dir) const {return len[dir];}

      // -------------------------------------------
      // ----------  Modify functions  -------------
      // -------------------------------------------
    void setLo(const Real* lo);
    void setLo(const Array<Real> &lo);
    void setHi(const Real* hi);
    void setHi(const Array<Real> &hi);
    void setLo(int indx, Real lo);
    void setHi(int indx, Real hi);

      // -------------------------------------------
      // ----------  Boolean functions  ------------
      // -------------------------------------------
    int ok() const;
    int contains(const Real* point);
    int contains(const RealBox &bx);
    int intersects(const RealBox &bx);

      // -------------------------------------------
      // -------  intersection functions  ----------
      // -------------------------------------------
    RealBox& operator &= (const RealBox& bx);
    RealBox  operator &  (const RealBox& bx);

      // -------------------------------------------
      // ---------  conversion functions  ----------
      // -------------------------------------------
    Box  CCBox(const Real *dx, const Real *base);
    Box  NodeBox(const Real *dx, const Real* base);

protected:
    void computeBoxLen() {
	D_EXPR(len[0] = xhi[0]-xlo[0],
	       len[1] = xhi[1]-xlo[1],
	       len[2] = xhi[2]-xlo[2]);
    }
};

#endif
