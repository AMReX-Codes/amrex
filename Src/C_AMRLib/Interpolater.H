#ifndef _INTERPOLATER_H_
#define _INTERPOLATER_H_ 
//
// $Id: Interpolater.H,v 1.10 2002-09-27 19:20:00 almgren Exp $
//
#include <Box.H>
#include <BCRec.H>
#include <REAL.H>

class Geometry;
class FArrayBox;

//
//@Man:
//@Memo: Virtual base class for interpolaters.
/*@Doc:

  Specifies interpolater interface for coarse-to-fine interpolation in space.
*/

class Interpolater
{
public:
    //
    //@ManDoc: The destructor.
    //
    virtual ~Interpolater ();

    /*@ManDoc: Returns coarsened box given fine box and refinement ratio.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual Box CoarseBox (const Box& fine,
                           int        ratio) = 0;

    /*@ManDoc: Returns coarsened box given fine box and refinement ratio.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio) = 0;

    /*@ManDoc: Coarse to fine interpolation in space.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual void interp (const FArrayBox& crse,
                         int              crse_comp,
                         FArrayBox&       fine,
                         int              fine_comp,
                         int              ncomp,
                         const Box&       fine_region,
                         const IntVect&   ratio,
                         const Geometry&  crse_geom,
                         const Geometry&  fine_geom,
                         Array<BCRec>&    bcr) = 0;

    //
    //@ManDoc: Re-visit the interpolation to protect against under- or overshoots.
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
                          Array<BCRec>&    bcr) {};
};

//
//@Man:
//@Memo: Bilinear interpolation on node centered data.
/*@Doc:

  Bilinear interpolation on node centered data.
*/

class NodeBilinear
    :
    public Interpolater
{
public:
    //
    //@ManDoc: The destructor.
    //
    virtual ~NodeBilinear ();
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio);
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio);
    //
    //@ManDoc: Coarse to fine interpolation in space.
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
                         Array<BCRec>&    bcr);
};

//
//@Man:
//@Memo: Bilinear interpolation on cell centered data.
/*@Doc:

  Bilinear interpolation on cell centered data.
*/

class CellBilinear
    :
    public Interpolater
{
public:
    //
    //@ManDoc: The destructor.
    //
    virtual ~CellBilinear ();
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio);
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio);
    //
    //@ManDoc: Coarse to fine interpolation in space.
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
                         Array<BCRec>&    bcr);
};

//
//@Man:
//@Memo: Conservative interpolation on cell centered data
/*@Doc:

  Conservative interpolation on cell centered data
*/

class CellConservative
    :
    public Interpolater
{
public:
    //
    //@ManDoc: The constructor.
    //
    CellConservative (bool limit = true);
    //
    //@ManDoc: The destructor.
    //
    virtual ~CellConservative ();
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio);
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio);
    //
    //@ManDoc: Coarse to fine interpolation in space.
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
                         Array<BCRec>&    bcr);
private:

    bool   do_limited_slope;
};

//
//@Man:
//@Memo: Linear conservative interpolation on cell centered data
/*@Doc:

  Linear conservative interpolation on cell centered data
  I.e, conservative interpolation with a limiting scheme that 
  preserves the value of any linear combination  of the
  fab components; e.g.,

  if sum\_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
  sum\_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
  in all fine cells if,jf covering coarse cell ic,jc.

  If do\_linear\_limiting = 0, the interpolation scheme is identical to
  the used in CellConservative for do\_limited\_slopes=1.
*/

class CellConservativeLinear
    :
    public Interpolater
{
public:
    //
    //@ManDoc: The constructor.
    //
    CellConservativeLinear (bool do_linear_limiting_= true);
    //
    //@ManDoc: The destructor.
    //
    virtual ~CellConservativeLinear ();
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio);
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio);
    //
    //@ManDoc: Coarse to fine interpolation in space.
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
                         Array<BCRec>&    bcr);

private:
    bool do_linear_limiting;
};

//
//@Man:
//@Memo: Lin. cons. interp. on cc data with protection against under/over-shoots.
/*@Doc:

  Linear conservative interpolation on cell centered data
  but with protection against undershoots or overshoots.
*/

class CellConservativeProtected
    :
    public Interpolater
{
public:
    //
    //@ManDoc: The constructor.
    //
    CellConservativeProtected ();
    //
    //@ManDoc: The destructor.
    //
    virtual ~CellConservativeProtected ();
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio);
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio);
    //
    //@ManDoc: Coarse to fine interpolation in space.
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
                         Array<BCRec>&    bcr);

    //
    //@ManDoc: Re-visit the interpolation to protect against under- or overshoots.
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
                          Array<BCRec>&    bcr);
};

//
//@Man:
//@Memo: Quadratic interpolation on cell centered data.
/*@Doc:

  Quadratic interpolation on cell centered data.
*/

class CellQuadratic
    :
    public Interpolater
{
public:
    //
    //@ManDoc: The constructor.
    //
    CellQuadratic (bool limit = true);
    //
    //@ManDoc: The destructor.
    //
    virtual ~CellQuadratic ();
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio);
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio);
    //
    //@ManDoc: Coarse to fine interpolation in space.
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
                         Array<BCRec>&    bcr);
private:

    bool  do_limited_slope;
};

//
//@Man:
//@Memo: Piecewise Constant interpolation on cell centered data.
/*@Doc:

  Piecewise Constant interpolation on cell centered data.
*/

class PCInterp
    :
    public Interpolater
{
public:
    //
    //@ManDoc: The destructor.
    //
    virtual ~PCInterp ();
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box& fine,
                           int        ratio);
    //
    //@ManDoc: Returns coarsened box given fine box and refinement ratio.
    //
    virtual Box CoarseBox (const Box&     fine,
                           const IntVect& ratio);
    //
    //@ManDoc: Coarse to fine interpolation in space.
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
                         Array<BCRec>&    bcr);
};

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//

extern NodeBilinear              node_bilinear_interp;
extern CellBilinear              cell_bilinear_interp;
extern CellConservative          cell_cons_interp;
extern CellQuadratic             quadratic_interp;
extern CellConservative          unlimited_cc_interp;
extern PCInterp                  pc_interp;
extern CellConservativeLinear    lincc_interp;
extern CellConservativeLinear    nonlincc_interp;
extern CellConservativeProtected protected_interp;

#endif /*_INTERPOLATER_H_*/
