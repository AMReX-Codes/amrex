//BL_COPYRIGHT_NOTICE

#ifndef _INTERPOLATER_H_
#define _INTERPOLATER_H_ 

//
// $Id: Interpolater.H,v 1.2 1997-11-24 18:52:30 lijewski Exp $
//

#define _NodeBilinear_H_      "%W%  %G%"
#define _CellBilinear_H_      "%W%  %G%"
#define _CellConservative_H_  "%W%  %G%"
#define _PCInterp_H_	      "%W%  %G%"

#include <Box.H>
#include <Geometry.H>
#include <FArrayBox.H>
#include <BCRec.H>

#define UNLIMITED_SLOPE 0
#define LIMITED_SLOPE 1

// virtual base class
class Interpolater {
public:
    Interpolater() {};
    virtual ~Interpolater();
    virtual BOX  CoarseBox(const BOX& fine, int ratio) = 0;
    virtual BOX  CoarseBox(const BOX& fine, const IntVect& ratio) = 0;
    virtual void interp(FARRAYBOX& crse, int crse_comp,
                        FARRAYBOX& fine, int fine_comp,
                        int ncomp,
                        const BOX& fine_region, const IntVect & ratio,
                        const Geometry& crse_geom,
                        const Geometry& fine_geom,
                        Array<BCRec> &bcr) = 0;
};

// ------------------------------------------------------------
// ------  Bilinear interpolation on node centered data  ------
// ------------------------------------------------------------
class NodeBilinear : public Interpolater {
private:
    REAL *strip;       // temp space passed to FORTRAN
    int   strip_len;
    REAL *slope;
    int   slope_len;
public:
    NodeBilinear();
    virtual ~NodeBilinear();
    virtual BOX CoarseBox(const BOX& fine, int ratio)
    { return ::coarsen(fine,ratio); }
    virtual BOX CoarseBox(const BOX& fine, const IntVect& ratio)
    { return ::coarsen(fine,ratio); }
    virtual void interp(FARRAYBOX& crse, int crse_comp,
                        FARRAYBOX& fine, int fine_comp,
                        int ncomp,
                        const BOX& fine_region, const IntVect & ratio,
                        const Geometry& crse_geom,
                        const Geometry& fine_geom,
                        Array<BCRec> &bcr);
public:
};

// ------------------------------------------------------------
// ------  Bilinear interpolation on cell centered data  ------
// ------------------------------------------------------------
class CellBilinear : public Interpolater {
private:
    REAL *strip;       // temp space passed to FORTRAN
    int   strip_len;
    REAL *slope;
    int   slope_len;
public:
    CellBilinear();
    virtual ~CellBilinear();
    virtual BOX CoarseBox(const BOX& fine, int ratio);
    virtual BOX CoarseBox(const BOX& fine, const IntVect& ratio);
    virtual void interp(FARRAYBOX& crse, int crse_comp,
                        FARRAYBOX& fine, int fine_comp,
                        int ncomp,
                        const BOX& fine_region, const IntVect & ratio,
                        const Geometry& crse_geom,
                        const Geometry& fine_geom,
                        Array<BCRec> &bcr);
public:
};

// ------------------------------------------------------------
// ------  Conservative interpolation on cell centered data ---
// ------------------------------------------------------------
class CellConservative : public Interpolater {
private:
    int do_limited_slope;
    REAL *strip;       // temp space passed to FORTRAN
    int   strip_len;
    REAL *cslope;
    int   slope_len;
public:
    CellConservative(int limit=1);
    virtual ~CellConservative();
    virtual BOX CoarseBox(const BOX& fine, int ratio);
    virtual BOX CoarseBox(const BOX& fine, const IntVect& ratio);
    virtual void interp(FARRAYBOX& crse, int crse_comp,
                        FARRAYBOX& fine, int fine_comp,
                        int ncomp,
                        const BOX& fine_region, const IntVect & ratio,
                        const Geometry& crse_geom,
                        const Geometry& fine_geom,
                        Array<BCRec> &bcr);
public:
};

// -------------------------------------------------------------------
// ------  Piecewise Constant interpolation on cell centered data  ---
// -------------------------------------------------------------------
class PCInterp : public Interpolater {
private:
    REAL *strip;       // temp space passed to FORTRAN
    int   strip_len;
public:
    PCInterp();
    virtual ~PCInterp();
    virtual BOX CoarseBox(const BOX& fine, int ratio)
    { return ::coarsen(fine,ratio); }
    virtual BOX CoarseBox(const BOX& fine, const IntVect& ratio)
    { return ::coarsen(fine,ratio); }
    virtual void interp(FARRAYBOX& crse, int crse_comp,
                        FARRAYBOX& fine, int fine_comp,
                        int ncomp,
                        const BOX& fine_region, const IntVect & ratio,
                        const Geometry& crse_geom,
                        const Geometry& fine_geom,
                        Array<BCRec> &bcr);
public:
};


// ------------------------------------------------------------
// ------  CONSTRUCT A GLOBAL OBJECT OF EACH VERSION  ---------
// ------------------------------------------------------------
extern NodeBilinear     node_bilinear_interp;
extern CellBilinear     cell_bilinear_interp;
extern CellConservative cell_cons_interp;
extern CellConservative unlimited_cc_interp;
extern PCInterp         pc_interp;

#endif

