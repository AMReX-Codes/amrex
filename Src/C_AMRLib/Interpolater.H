#ifndef _INTERPOLATER_H_
#define _INTERPOLATER_H_ 

//
// $Id: Interpolater.H,v 1.4 1997-12-02 21:15:48 lijewski Exp $
//

#include <Box.H>
#include <Geometry.H>
#include <FArrayBox.H>
#include <BCRec.H>

#define UNLIMITED_SLOPE 0
#define LIMITED_SLOPE 1

//
// Virtual base class.
//

class Interpolater
{
public:
    Interpolater () {}
    virtual ~Interpolater ();
    virtual Box  CoarseBox (const Box& fine, int ratio) = 0;
    virtual Box  CoarseBox (const Box& fine, const IntVect& ratio) = 0;
    virtual void interp (FArrayBox& crse, int crse_comp,
                         FArrayBox& fine, int fine_comp,
                         int ncomp,
                         const Box& fine_region, const IntVect & ratio,
                         const Geometry& crse_geom,
                         const Geometry& fine_geom,
                         Array<BCRec> &bcr) = 0;
};

//
// Bilinear interpolation on node centered data.
//

class NodeBilinear
    :
    public Interpolater
{
public:
    NodeBilinear ();
    virtual ~NodeBilinear ();
    virtual Box CoarseBox (const Box& fine, int ratio);
    virtual Box CoarseBox (const Box& fine, const IntVect& ratio);
    virtual void interp (FArrayBox& crse, int crse_comp,
                        FArrayBox& fine, int fine_comp,
                        int ncomp,
                        const Box& fine_region, const IntVect & ratio,
                        const Geometry& crse_geom,
                        const Geometry& fine_geom,
                        Array<BCRec> &bcr);
private:
    Real* strip;       // temp space passed to FORTRAN
    int   strip_len;
    Real* slope;
    int   slope_len;
};

//
// Bilinear interpolation on cell centered data.
//

class CellBilinear
    :
    public Interpolater
{
public:
    CellBilinear ();
    virtual ~CellBilinear ();
    virtual Box CoarseBox (const Box& fine, int ratio);
    virtual Box CoarseBox (const Box& fine, const IntVect& ratio);
    virtual void interp (FArrayBox& crse, int crse_comp,
                         FArrayBox& fine, int fine_comp,
                         int ncomp,
                         const Box& fine_region, const IntVect & ratio,
                         const Geometry& crse_geom,
                         const Geometry& fine_geom,
                         Array<BCRec> &bcr);
private:
    Real* strip;       // temp space passed to FORTRAN
    int   strip_len;
    Real* slope;
    int   slope_len;
};

//
// Conservative interpolation on cell centered data ---
//
class CellConservative
    :
    public Interpolater
{
public:
    CellConservative (bool limit = true);
    virtual ~CellConservative ();
    virtual Box CoarseBox (const Box& fine, int ratio);
    virtual Box CoarseBox (const Box& fine, const IntVect& ratio);
    virtual void interp (FArrayBox& crse, int crse_comp,
                         FArrayBox& fine, int fine_comp,
                         int ncomp,
                         const Box& fine_region, const IntVect & ratio,
                         const Geometry& crse_geom,
                         const Geometry& fine_geom,
                         Array<BCRec> &bcr);
private:
    Real*  strip;       // temp space passed to FORTRAN
    int    strip_len;
    Real*  cslope;
    int    slope_len;
    bool   do_limited_slope;
};

//
// "Linear" conservative interpolation on cell centered data
// I.e, conservative interpolation with a limiting scheme that 
// preserves the value of any linear combination  of the
// fab components--e.g.,
// if sum_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
//    sum_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
// in all fine cells if,jf covering coarse cell ic,jc.
//
// If do_linear_limiting = 0, the interpolation scheme is identical to
// the used in CellConservative for do_limited_slopes=1.
//

class CellConservativeLinear
    :
    public Interpolater
{
public:
    CellConservativeLinear (bool do_linear_limiting_= true);
    virtual ~CellConservativeLinear ();
    virtual Box CoarseBox (const Box& fine, int ratio);
    virtual Box CoarseBox (const Box& fine, const IntVect& ratio);
    virtual void interp (FArrayBox& crse, int crse_comp,
                         FArrayBox& fine, int fine_comp,
                         int ncomp,
                         const Box& fine_region, const IntVect & ratio,
                         const Geometry& crse_geom,
                         const Geometry& fine_geom,
                         Array<BCRec> &bcr);
private:
    bool do_linear_limiting;
};

//
// Quadratic interpolation on cell centered data.
//

class CellQuadratic
    :
    public Interpolater
{
public:
    CellQuadratic (bool limit = true);
    virtual ~CellQuadratic ();
    virtual Box CoarseBox (const Box& fine, int ratio);
    virtual Box CoarseBox (const Box& fine, const IntVect& ratio);
    virtual void interp (FArrayBox& crse, int crse_comp,
                         FArrayBox& fine, int fine_comp,
                         int ncomp,
                         const Box& fine_region, const IntVect & ratio,
                         const Geometry& crse_geom,
                         const Geometry& fine_geom,
                         Array<BCRec> &bcr);
private:
    Real* strip;       // temp space passed to FORTRAN
    int   strip_len;
    Real* cslope;
    int   slope_len;
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
    PCInterp ();
    virtual ~PCInterp ();
    virtual Box CoarseBox (const Box& fine, int ratio);
    virtual Box CoarseBox (const Box& fine, const IntVect& ratio);
    virtual void interp (FArrayBox& crse, int crse_comp,
                         FArrayBox& fine, int fine_comp,
                         int ncomp,
                         const Box& fine_region, const IntVect & ratio,
                         const Geometry& crse_geom,
                         const Geometry& fine_geom,
                         Array<BCRec> &bcr);
private:
    Real* strip;       // temp space passed to FORTRAN
    int   strip_len;
};

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//

extern NodeBilinear           node_bilinear_interp;
extern CellBilinear           cell_bilinear_interp;
extern CellConservative       cell_cons_interp;
extern CellQuadratic          quadratic_interp;
extern CellConservative       unlimited_cc_interp;
extern PCInterp               pc_interp;
extern CellConservativeLinear lincc_interp;
extern CellConservativeLinear nonlincc_interp;

#endif

