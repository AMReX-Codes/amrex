//BL_COPYRIGHT_NOTICE

//
// $Id: Interpolater.cpp,v 1.14 1998-11-03 18:16:38 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cmath>
#include <climits>
#else
#include <math.h>
#include <limits.h>
#endif

#include <Interpolater.H>
#include <INTERP_F.H>

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//

NodeBilinear           node_bilinear_interp;
CellBilinear           cell_bilinear_interp;
CellConservative       cell_cons_interp;
CellQuadratic          quadratic_interp;
CellConservative       unlimited_cc_interp(0);
PCInterp               pc_interp;
CellConservativeLinear lincc_interp;
CellConservativeLinear nonlincc_interp(0);

Interpolater::~Interpolater () {}

NodeBilinear::NodeBilinear ()
{
    strip = slope = 0;
    strip_len = slope_len = 0;
}

NodeBilinear::~NodeBilinear ()
{
    delete [] strip;
    delete [] slope;
}

Box
NodeBilinear::CoarseBox (const Box& fine,
                         int        ratio)
{
    Box b = ::coarsen(fine,ratio);

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        if (b.length(i) < 2)
        {
            //
            // Don't want degenerate boxes.
            //
            b.growHi(i,1);
        }
    }

    return b;
}

Box
NodeBilinear::CoarseBox (const Box&     fine,
                         const IntVect& ratio)
{
    Box b = ::coarsen(fine,ratio);

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        if (b.length(i) < 2)
        {
            //
            // Don't want degenerate boxes.
            //
            b.growHi(i,1);
        }
    }

    return b;
}

void
NodeBilinear::interp (const FArrayBox& crse,
                      int              crse_comp,
                      FArrayBox&       fine,
                      int              fine_comp,
                      int              ncomp,
                      const Box&       fine_region,
                      const IntVect&   ratio,
                      const Geometry& /* crse_geom */,
                      const Geometry& /* fine_geom */,
                      Array<BCRec>&   /*bcr*/)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo = crse.box().loVect();
    const int* chi = crse.box().hiVect();
    const int* flo = fine.loVect();
    const int* fhi = fine.hiVect();
    const int* lo  = fine_region.loVect();
    const int* hi  = fine_region.hiVect();
    int num_slope  = (int) pow(2.0,BL_SPACEDIM)-1;
    int len0       = crse.box().length()[0];
    int slp_len    = num_slope*len0;
    if (slope_len < slp_len)
    {
        delete [] slope;
        slope_len = slp_len;
        slope = new Real[slope_len];
    }
    int strp_len = len0*ratio[0];
    for (int n = 1; n < BL_SPACEDIM; n++ )
    {
      strp_len *= ratio[n]+1;
    }
    if (strip_len < strp_len)
    {
        delete [] strip;
        strip_len = strp_len;
        strip     = new Real[strip_len];
    }
    int strip_lo      = ratio[0] * clo[0];
    int strip_hi      = ratio[0] * chi[0];
    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);
    const int* ratioV = ratio.getVect();

    FORT_NBINTERP (cdat,ARLIM(clo),ARLIM(chi),ARLIM(clo),ARLIM(chi),
                   fdat,ARLIM(flo),ARLIM(fhi),ARLIM(lo),ARLIM(hi),
                   D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),&ncomp,
                   slope,&num_slope,strip,&strip_lo,&strip_hi);
}

CellBilinear::CellBilinear ()
{
    strip = slope = 0;
    strip_len = slope_len = 0;
}

CellBilinear::~CellBilinear ()
{
    delete [] strip;
    delete [] slope;
}

Box
CellBilinear::CoarseBox (const Box& fine,
                         int        ratio)
{
    return CoarseBox(fine, ratio*IntVect::TheUnitVector());
}

Box
CellBilinear::CoarseBox (const Box&     fine,
                         const IntVect& ratio)
{
    const int* lo = fine.loVect();
    const int* hi = fine.hiVect();

    Box crse(::coarsen(fine,ratio));
    const int* clo = crse.loVect();
    const int* chi = crse.hiVect();

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        int iratio = ratio[i];
        int hrat   = iratio/2;
        if (lo[i] <  clo[i]*ratio[i] + hrat)
            crse.growLo(i,1);
        if (hi[i] >= chi[i]*ratio[i] + hrat)
            crse.growHi(i,1);
    }
    return crse;
}

void
CellBilinear::interp (const FArrayBox& crse, int crse_comp,
                      FArrayBox& fine, int fine_comp,
                      int ncomp,
                      const Box& fine_region, const IntVect & ratio,
                      const Geometry& /* crse_geom */,
                      const Geometry& /* fine_geom */,
                      Array<BCRec>& /*bcr*/)
{
    BoxLib::Error("interp: not implemented");
    //
    // Set up to call FORTRAN.
    //
    const int* clo = crse.box().loVect();
    const int* chi = crse.box().hiVect();
    const int* flo = fine.loVect();
    const int* fhi = fine.hiVect();
    const int* lo  = fine_region.loVect();
    const int* hi  = fine_region.hiVect();
    int num_slope  = (int) pow(2.0,BL_SPACEDIM)-1;
    int len0       = crse.box().length()[0];
    int slp_len    = num_slope*len0;
    if (slope_len < slp_len)
    {
        delete [] slope;
        slope_len = slp_len;
        slope     = new Real[slope_len];
    }
    int strp_len = len0*ratio[0];
    if (strip_len < strp_len)
    {
        delete [] strip;
        strip_len = strp_len;
        strip     = new Real[strip_len];
    }
    int strip_lo = ratio[0] * clo[0];
    int strip_hi = ratio[0] * chi[0];

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);
    const int* ratioV = ratio.getVect();

    FORT_CBINTERP (cdat,ARLIM(clo),ARLIM(chi),ARLIM(clo),ARLIM(chi),
                   fdat,ARLIM(flo),ARLIM(fhi),ARLIM(lo),ARLIM(hi),
                   D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),&ncomp,
                   slope,&num_slope,strip,&strip_lo,&strip_hi);
}

CellConservative::CellConservative (bool limit)
{
    strip = cslope = 0;
    strip_len = slope_len = 0;
    do_limited_slope = limit;
}

CellConservative::~CellConservative ()
{
    delete [] strip;
    delete [] cslope;
}

Box
CellConservative::CoarseBox (const Box&     fine,
                             const IntVect& ratio)
{
    Box crse = ::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellConservative::CoarseBox (const Box& fine,
                             int        ratio)
{
    Box crse = ::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

static
Array<int>
GetBCArray (const Array<BCRec>& bcr)
{
    Array<int> bc(2*BL_SPACEDIM*bcr.length());

    for (int n = 0; n < bcr.length(); n++)
    {
        const int* b_rec = bcr[n].vect();

        for (int m = 0; m < 2*BL_SPACEDIM; m++)
        {
            bc[2*BL_SPACEDIM*n + m] = b_rec[m];
        }
    }

    return bc;
}

void
CellConservative::interp (const FArrayBox& crse,
                          int              crse_comp,
                          FArrayBox&       fine,
                          int              fine_comp,
                          int              ncomp,
                          const Box&       fine_region,
                          const IntVect &  ratio,
                          const Geometry&  crse_geom,
                          const Geometry&  fine_geom,
                          Array<BCRec>&    bcr)
{
    assert(bcr.length() >= ncomp);
    assert(fine_geom.Domain().contains(fine_region));
    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();
    Box crse_bx            = ::coarsen(target_fine_region,ratio);
    Box fslope_bx          = ::refine(crse_bx,ratio);
    Box cslope_bx          = crse_bx;
    cslope_bx.grow(1);
    assert(crse.box().contains(cslope_bx));
    //
    // Alloc temp space for coarse grid slopes.
    //
    long t_long = cslope_bx.numPts();
    assert(t_long < INT_MAX);
    int c_len = int(t_long);
    if (slope_len < BL_SPACEDIM*c_len)
    {
        slope_len = BL_SPACEDIM*c_len;
        delete [] cslope;
        cslope = new Real[slope_len];
    }

    int loslp = cslope_bx.index(crse_bx.smallEnd());
    int hislp = cslope_bx.index(crse_bx.bigEnd());

    t_long = cslope_bx.numPts();
    assert(t_long < INT_MAX);
    int cslope_vol = int(t_long);
    int clo        = 1 - loslp;
    int chi        = clo + cslope_vol - 1;
    c_len          = hislp - loslp + 1;
    //
    // Alloc temp space for one strip of fine grid slopes.
    //
    int dir;
    int f_len = fslope_bx.longside(dir);
    if (strip_len < (BL_SPACEDIM+2)*f_len)
    {
        strip_len = (BL_SPACEDIM+2)*f_len;
        delete [] strip;
        strip = new Real[strip_len];
    }

    Real* fstrip = strip;
    Real* foff   = fstrip + f_len;
    Real* fslope = foff + f_len;
    //
    // Get coarse in fine edge centered volume coordinates.
    //
    Array<Real> fvc[BL_SPACEDIM];
    Array<Real> cvc[BL_SPACEDIM];
    for (dir = 0; dir < BL_SPACEDIM; dir++)
    {
        fine_geom.GetEdgeVolCoord(fvc[dir],target_fine_region,dir);
        crse_geom.GetEdgeVolCoord(cvc[dir],crse_bx,dir);
    }
    //
    // Alloc tmp space for slope calc and to allow for vectorization.
    //
    Real* fdat        = fine.dataPtr(fine_comp);
    const Real* cdat  = crse.dataPtr(crse_comp);
    const int* flo    = fine.loVect();
    const int* fhi    = fine.hiVect();
    const int* fblo   = target_fine_region.loVect();
    const int* fbhi   = target_fine_region.hiVect();
    const int* cblo   = crse_bx.loVect();
    const int* cbhi   = crse_bx.hiVect();
    const int* cflo   = crse.loVect();
    const int* cfhi   = crse.hiVect();
    const int* fslo   = fslope_bx.loVect();
    const int* fshi   = fslope_bx.hiVect();
    int slope_flag    = (do_limited_slope ? 1 : 0);
    Array<int> bc     = GetBCArray(bcr);
    const int* ratioV = ratio.getVect();

    FORT_CCINTERP (fdat,ARLIM(flo),ARLIM(fhi),
                   ARLIM(fblo), ARLIM(fbhi),
                   &ncomp,D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),
                   cdat,&clo,&chi,
                   ARLIM(cblo), ARLIM(cbhi),
                   fslo,fshi,
                   cslope,&c_len,fslope,fstrip,&f_len,foff,
                   bc.dataPtr(), &slope_flag,
                   D_DECL(fvc[0].dataPtr(),fvc[1].dataPtr(),fvc[2].dataPtr()),
                   D_DECL(cvc[0].dataPtr(),cvc[1].dataPtr(),cvc[2].dataPtr()));
}

CellConservativeLinear::CellConservativeLinear (bool do_linear_limiting_)
{
    do_linear_limiting = do_linear_limiting_;
}

CellConservativeLinear::~CellConservativeLinear ()
{}

Box
CellConservativeLinear::CoarseBox (const Box&     fine,
                                   const IntVect& ratio)
{
    Box crse = ::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellConservativeLinear::CoarseBox (const Box& fine,
                                   int        ratio)
{
    Box crse(::coarsen(fine,ratio));
    crse.grow(1);
    return crse;
}

void
CellConservativeLinear::interp (const FArrayBox& crse, int crse_comp,
                         FArrayBox& fine, int fine_comp,
                         int ncomp,
                         const Box& fine_region, const IntVect & ratio,
                         const Geometry& crse_geom,
                         const Geometry& fine_geom,
                         Array<BCRec>& bcr)
{
    assert(bcr.length() >= ncomp);
    assert(fine_geom.Domain().contains(fine_region));
    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();
    //
    // crse_bx is coarsening of target_fine_region, grown by 1.
    //
    Box crse_bx = CoarseBox(target_fine_region,ratio);
    //
    // Slopes are needed only on coarsening of target_fine_region.
    //
    Box cslope_bx(crse_bx);
    cslope_bx.grow(-1);
    //
    // Get coarse in fine edge centered volume coordinates.
    //
    Array<Real> fvc[BL_SPACEDIM];
    Array<Real> cvc[BL_SPACEDIM];
    int dir;
    for (dir = 0; dir < BL_SPACEDIM; dir++)
    {
        fine_geom.GetEdgeVolCoord(fvc[dir],target_fine_region,dir);
        crse_geom.GetEdgeVolCoord(cvc[dir],crse_bx,dir);
    }
    //
    // alloc tmp space for slope calc.
    //
    // In ucc_slopes and lcc_slopes , there is a slight abuse of 
    // the number of compenents argument
    // --> there is a slope for each component in each coordinate 
    //     direction
    //
    FArrayBox ucc_slopes(cslope_bx,ncomp*BL_SPACEDIM);
    FArrayBox lcc_slopes(cslope_bx,ncomp*BL_SPACEDIM);
    FArrayBox slope_factors(cslope_bx,BL_SPACEDIM);

    Real* fdat       = fine.dataPtr(fine_comp);
    const Real* cdat = crse.dataPtr(crse_comp);
    Real* ucc_xsldat = ucc_slopes.dataPtr(0);
    Real* lcc_xsldat = lcc_slopes.dataPtr(0);
    Real* xslfac_dat = slope_factors.dataPtr(0);
    Real* ucc_ysldat = ucc_slopes.dataPtr(ncomp);
    Real* lcc_ysldat = lcc_slopes.dataPtr(ncomp);
    Real* yslfac_dat = slope_factors.dataPtr(1);
#if (BL_SPACEDIM==3)
    Real* ucc_zsldat = ucc_slopes.dataPtr(2*ncomp);
    Real* lcc_zsldat = lcc_slopes.dataPtr(2*ncomp);
    Real* zslfac_dat = slope_factors.dataPtr(2);
#endif
    
    const int* flo    = fine.loVect();
    const int* fhi    = fine.hiVect();
    const int* clo    = crse.loVect();
    const int* chi    = crse.hiVect();
    const int* fblo   = target_fine_region.loVect();
    const int* fbhi   = target_fine_region.hiVect();
    const int* csbhi  = cslope_bx.hiVect();
    const int* csblo  = cslope_bx.loVect();
    int lin_limit     = (do_linear_limiting ? 1 : 0);
    const int* cvcblo = crse_bx.loVect();
    const int* fvcblo = target_fine_region.loVect();

    int cvcbhi[BL_SPACEDIM];
    int fvcbhi[BL_SPACEDIM];

    for (dir=0; dir<BL_SPACEDIM; dir++)
    {
        cvcbhi[dir] = cvcblo[dir] + cvc[dir].length() - 1;
        fvcbhi[dir] = fvcblo[dir] + fvc[dir].length() - 1;
    }

    Real* voffx = new Real[fvc[0].length()];
    Real* voffy = new Real[fvc[1].length()];

#if (BL_SPACEDIM==3)
    Real* voffz = new Real[fvc[2].length()];
#endif

    Array<int> bc     = GetBCArray(bcr);
    const int* ratioV = ratio.getVect();

    FORT_LINCCINTERP (fdat,ARLIM(flo),ARLIM(fhi),
                      fblo, fbhi,
                      ARLIM(fvcblo), ARLIM(fvcbhi),
                      cdat,ARLIM(clo),ARLIM(chi),
                      ARLIM(cvcblo), ARLIM(cvcbhi),
                      ucc_xsldat, lcc_xsldat, xslfac_dat,
                      ucc_ysldat, lcc_ysldat, yslfac_dat,
#if (BL_SPACEDIM==3)
                      ucc_zsldat, lcc_zsldat, zslfac_dat,
#endif
                      ARLIM(csblo), ARLIM(csbhi),
                      csblo, csbhi,
                      &ncomp,D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),
                      bc.dataPtr(), &lin_limit,
                      D_DECL(fvc[0].dataPtr(),fvc[1].dataPtr(),fvc[2].dataPtr()),
                      D_DECL(cvc[0].dataPtr(),cvc[1].dataPtr(),cvc[2].dataPtr()),
                      D_DECL(voffx,voffy,voffz)
                      );
    delete [] voffx;
    delete [] voffy;
#if (BL_SPACEDIM==3)
    delete [] voffz;
#endif

}

CellQuadratic::CellQuadratic (bool limit)
{
    strip = cslope = 0;
    strip_len = slope_len = 0;
    do_limited_slope = limit;
}

CellQuadratic::~CellQuadratic ()
{
    delete [] strip;
    delete [] cslope;
}

Box
CellQuadratic::CoarseBox (const Box&     fine,
                          const IntVect& ratio)
{
    Box crse = ::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellQuadratic::CoarseBox (const Box& fine,
                          int        ratio)
{
    Box crse = ::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

void
CellQuadratic::interp (const FArrayBox& crse, int crse_comp,
                       FArrayBox& fine, int fine_comp,
                       int ncomp,
                       const Box& fine_region, const IntVect & ratio,
                       const Geometry& crse_geom,
                       const Geometry& fine_geom,
                       Array<BCRec>& bcr)
{
    assert(bcr.length() >= ncomp);
    assert(fine_geom.Domain().contains(fine_region));
    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    Box crse_bx(::coarsen(target_fine_region,ratio));
    Box fslope_bx(::refine(crse_bx,ratio));
    Box cslope_bx(crse_bx);
    cslope_bx.grow(1);
    assert(crse.box().contains(cslope_bx));
    //
    // Alloc temp space for coarse grid slopes: here we use 5 
    // instead of BL_SPACEDIM because of the x^2, y^2 and xy terms
    //
    long t_long = cslope_bx.numPts();
    assert(t_long < INT_MAX);
    int c_len = int(t_long);
    if (slope_len < 5*c_len)
    {
        slope_len = 5*c_len;
        delete [] cslope;
        cslope = new Real[slope_len];
    }
    int loslp = cslope_bx.index(crse_bx.smallEnd());
    int hislp = cslope_bx.index(crse_bx.bigEnd());

    t_long = cslope_bx.numPts();
    assert(t_long < INT_MAX);
    int cslope_vol = int(t_long);
    int clo        = 1 - loslp;
    int chi        = clo + cslope_vol - 1;
    c_len          = hislp - loslp + 1;
    //
    // Alloc temp space for one strip of fine grid slopes: here we use 5 
    // instead of BL_SPACEDIM because of the x^2, y^2 and xy terms.
    //
    int dir;
    int f_len = fslope_bx.longside(dir);
    if (strip_len < (5+2)*f_len)
    {
        strip_len = (5+2)*f_len;
        delete [] strip;
        strip = new Real[strip_len];
    }
    Real* fstrip = strip;
    Real* foff   = fstrip + f_len;
    Real* fslope = foff + f_len;
    //
    // Get coarse in fine edge centered volume coordinates.
    //
    Array<Real> fvc[BL_SPACEDIM];
    Array<Real> cvc[BL_SPACEDIM];
    for (dir = 0; dir < BL_SPACEDIM; dir++)
    {
        fine_geom.GetEdgeVolCoord(fvc[dir],target_fine_region,dir);
        crse_geom.GetEdgeVolCoord(cvc[dir],crse_bx,dir);
    }
    //
    // Alloc tmp space for slope calc and to allow for vectorization.
    //
    Real* fdat        = fine.dataPtr(fine_comp);
    const Real* cdat  = crse.dataPtr(crse_comp);
    const int* flo    = fine.loVect();
    const int* fhi    = fine.hiVect();
    const int* fblo   = target_fine_region.loVect();
    const int* fbhi   = target_fine_region.hiVect();
    const int* cblo   = crse_bx.loVect();
    const int* cbhi   = crse_bx.hiVect();
    const int* cflo   = crse.loVect();
    const int* cfhi   = crse.hiVect();
    const int* fslo   = fslope_bx.loVect();
    const int* fshi   = fslope_bx.hiVect();
    int slope_flag    = (do_limited_slope ? 1 : 0);
    Array<int> bc     = GetBCArray(bcr);
    const int* ratioV = ratio.getVect();

    FORT_CQINTERP (fdat,ARLIM(flo),ARLIM(fhi),
                   ARLIM(fblo), ARLIM(fbhi),
                   &ncomp,D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),
                   cdat,&clo,&chi,
                   ARLIM(cblo), ARLIM(cbhi),
                   fslo,fshi,
                   cslope,&c_len,fslope,fstrip,&f_len,foff,
                   bc.dataPtr(), &slope_flag,
                   D_DECL(fvc[0].dataPtr(),fvc[1].dataPtr(),fvc[2].dataPtr()),
                   D_DECL(cvc[0].dataPtr(),cvc[1].dataPtr(),cvc[2].dataPtr()));
}

PCInterp::PCInterp ()
{
    strip = 0;
    strip_len = 0;
}

PCInterp::~PCInterp ()
{
    delete [] strip;
}

Box
PCInterp::CoarseBox (const Box& fine,
                     int        ratio)
{
    return ::coarsen(fine,ratio);
}

Box
PCInterp::CoarseBox (const Box&     fine,
                     const IntVect& ratio)
{
    return ::coarsen(fine,ratio);
}

void
PCInterp::interp (const FArrayBox& crse, int crse_comp,
                  FArrayBox& fine, int fine_comp,
                  int ncomp,
                  const Box& fine_region, const IntVect & ratio,
                  const Geometry& /* crse_geom */,
                  const Geometry& /* fine_geom */,
                  Array<BCRec>& /* bcr */)
{
    //
    // Set up to call FORTRAN.
    //
    const int* clo  = crse.box().loVect();
    const int* chi  = crse.box().hiVect();
    const int* flo  = fine.loVect();
    const int* fhi  = fine.hiVect();
    const int* fblo = fine_region.loVect();
    const int* fbhi = fine_region.hiVect();

    Box cregion(::coarsen(fine_region,ratio));

    const int* cblo = cregion.loVect();
    const int* cbhi = cregion.hiVect();
    int long_dir;
    int long_len = cregion.longside(long_dir);
    int s_len = long_len*ratio[long_dir];
    if (strip_len < s_len)
    {
        delete [] strip;
        strip_len = s_len;
        strip     = new Real[strip_len];
    }
    int strip_lo = ratio[long_dir] * cblo[long_dir];
    int strip_hi = ratio[long_dir] * (cbhi[long_dir]+1) - 1;
    //
    // Convert long_dir to FORTRAN (1 based) index.
    //
    long_dir++;
    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);
    const int* ratioV = ratio.getVect();


    FORT_PCINTERP (cdat,ARLIM(clo),ARLIM(chi),cblo,cbhi,
                   fdat,ARLIM(flo),ARLIM(fhi),fblo,fbhi,
                   &long_dir,D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),
                   &ncomp,strip,&strip_lo,&strip_hi);
}

