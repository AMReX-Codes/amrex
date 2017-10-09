#ifndef _INTERPOLATER_H_
#define _INTERPOLATER_H_ 

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_REAL.H>

namespace amrex {

class Geometry;
class FArrayBox;
class InterpolaterBoxCoarsener;

//
// Virtual base class for interpolaters.
//
// Specifies interpolater interface for coarse-to-fine interpolation in space.
//

class Interpolater
{
public:
    //
    // The destructor.
    //
    virtual ~Interpolater ();
    //
    // Returns coarsened box given fine box and refinement ratio.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) = 0;
    //
    // Returns coarsened box given fine box and refinement ratio.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) = 0;
    //
    // Coarse to fine interpolation in space.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec>&    bcr,
                         int              actual_comp,
                         int              actual_state) = 0;
    //
    // Re-visit the interpolation to protect against under- or overshoots.
    //
    virtual void protect (const FArrayBox& crse,
                          int              crse_comp,
                          FArrayBox&       fine,
                          int              fine_comp,
                          FArrayBox&       fine_state,
                          int              state_comp,
                          int              ncomp,
                          const Box&       fine_region,
                          const IntVect&   ratio,
                          const Geometry&  crse_geom,
                          const Geometry&  fine_geom,
                          Vector<BCRec>&    bcr) {};

    virtual InterpolaterBoxCoarsener BoxCoarsener (const IntVect& ratio);

    static Vector<int> GetBCArray (const Vector<BCRec>& bcr);
};

//
// Bilinear interpolation on node centered data.
//
// Bilinear interpolation on node centered data.
//

class NodeBilinear
    :
    public Interpolater
{
public:
    //
    // The destructor.
    //
    virtual ~NodeBilinear () override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;
    //
    // Coarse to fine interpolation in space.
    //
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec>&    bcr,
                         int              actual_comp,
                         int              actual_state) override;
};

//
// Bilinear interpolation on cell centered data.
//
// Bilinear interpolation on cell centered data.
//

class CellBilinear
    :
    public Interpolater
{
public:
    //
    // The destructor.
    //
    virtual ~CellBilinear () override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;
    //
    // Coarse to fine interpolation in space.
    //
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec>&    bcr,
                         int              actual_comp,
                         int              actual_state) override;
};

//
// Linear conservative interpolation on cell centered data
//
//
//  Linear conservative interpolation on cell centered data
//  I.e, conservative interpolation with a limiting scheme that 
//  preserves the value of any linear combination  of the
//  fab components; e.g.,
//
//  if sum_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
//  sum_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
//  in all fine cells if,jf covering coarse cell ic,jc.
//
//  If do_linear_limiting = 0, the interpolation scheme is identical to
//  the used in CellConservative for do_limited_slopes=1.
//

class CellConservativeLinear
    :
    public Interpolater
{
public:
    //
    // The constructor.
    //
    CellConservativeLinear (bool do_linear_limiting_= true);
    //
    // The destructor.
    //
    virtual ~CellConservativeLinear () override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;
    //
    // Coarse to fine interpolation in space.
    //
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec>&    bcr,
                         int              actual_comp,
                         int              actual_state) override;

protected:

    bool do_linear_limiting;
};

//
// Lin. cons. interp. on cc data with protection against under/over-shoots.
//
// Linear conservative interpolation on cell centered data
// but with protection against undershoots or overshoots.
//

class CellConservativeProtected
    :
    public Interpolater
{
public:
    //
    // The constructor.
    //
    CellConservativeProtected ();
    //
    // The destructor.
    //
    virtual ~CellConservativeProtected () override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;
    //
    // Coarse to fine interpolation in space.
    //
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec>&    bcr,
                         int              actual_comp,
                         int              actual_state) override;

    //
    // Re-visit the interpolation to protect against under- or overshoots.
    //
    virtual void protect (const FArrayBox& crse,
                          int              crse_comp,
                          FArrayBox&       fine,
                          int              fine_comp,
                          FArrayBox&       fine_state,
                          int              state_comp,
                          int              ncomp,
                          const Box&       fine_region,
                          const IntVect&   ratio,
                          const Geometry&  crse_geom,
                          const Geometry&  fine_geom,
                          Vector<BCRec>&    bcr) override;
};

//
// Quadratic interpolation on cell centered data.
//
// Quadratic interpolation on cell centered data.
//

class CellQuadratic
    :
    public Interpolater
{
public:
    //
    // The constructor.
    //
    CellQuadratic (bool limit = true);
    //
    // The destructor.
    //
    virtual ~CellQuadratic () override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;
    //
    // Coarse to fine interpolation in space.
    //
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec>&    bcr,
                         int              actual_comp,
                         int              actual_state) override;
private:

    bool  do_limited_slope;
};

//
// Piecewise Constant interpolation on cell centered data.
//

class PCInterp
    :
    public Interpolater
{
public:
    //
    // The destructor.
    //
    virtual ~PCInterp () override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;
    //
    // Coarse to fine interpolation in space.
    //
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec>&    bcr,
                         int              actual_comp,
                         int              actual_state) override;
};

//
// Conservative quartic interpolation on cell averaged data.
//
// An order 4 polynomial is used to fit the data.  For each cell involved
// in constructing the polynomial, the average of the polynomial inside that 
// cell is equal to the cell averaged value of the original data.
//

class CellConservativeQuartic
    :
    public Interpolater
{
public:
    //
    // The destructor.
    //
    virtual ~CellConservativeQuartic () override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) override;
    //
    // Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) override;
    //
    // Coarse to fine interpolation in space.
    //
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Vector<BCRec>&    bcr,
                         int              actual_comp,
                         int              actual_state) override;
};


//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
extern PCInterp                  pc_interp;
extern NodeBilinear              node_bilinear_interp;
extern CellBilinear              cell_bilinear_interp;
extern CellQuadratic             quadratic_interp;
extern CellConservativeLinear    lincc_interp;
extern CellConservativeLinear    cell_cons_interp;
extern CellConservativeProtected protected_interp;
extern CellConservativeQuartic   quartic_interp;

class InterpolaterBoxCoarsener
    : public BoxConverter
{
public:
    InterpolaterBoxCoarsener (Interpolater* mapper_, const IntVect& ratio_) 
	: mapper (mapper_), ratio(ratio_) { ; }
    virtual Box doit (const Box& fine) const;
    virtual BoxConverter* clone () const;
private:
    Interpolater* mapper;
    IntVect ratio;
};

}

#endif /*_INTERPOLATER_H_*/
