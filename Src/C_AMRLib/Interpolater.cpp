//
// $Id: Interpolater.cpp,v 1.31 2003-09-15 21:17:03 lijewski Exp $
//
#include <winstd.H>

#include <cmath>
#include <climits>

#include <FArrayBox.H>
#include <Geometry.H>
#include <Interpolater.H>
#include <INTERP_F.H>
#include <Profiler.H>

//
// Note that in 1D, CellConservativeLinear and CellQuadratic
// interpolation are turned off in a hardwired way.
//

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
PCInterp                  pc_interp;
NodeBilinear              node_bilinear_interp;
CellBilinear              cell_bilinear_interp;
CellQuadratic             quadratic_interp;
CellConservativeLinear    lincc_interp;
CellConservativeLinear    cell_cons_interp(0);
CellConservativeProtected protected_interp;

Interpolater::~Interpolater () {}

NodeBilinear::~NodeBilinear () {}

Box
NodeBilinear::CoarseBox (const Box& fine,
                         int        ratio)
{
    Box b = BoxLib::coarsen(fine,ratio);

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
    Box b = BoxLib::coarsen(fine,ratio);

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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::interp");
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

    Array<Real> strip(slp_len);

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);
    const int* ratioV = ratio.getVect();

    FORT_NBINTERP (cdat,ARLIM(clo),ARLIM(chi),ARLIM(clo),ARLIM(chi),
                   fdat,ARLIM(flo),ARLIM(fhi),ARLIM(lo),ARLIM(hi),
                   D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),&ncomp,
                   strip.dataPtr(),&num_slope);
}

CellBilinear::~CellBilinear () {}

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

    Box crse(BoxLib::coarsen(fine,ratio));
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
CellBilinear::interp (const FArrayBox& crse,
                      int              crse_comp,
                      FArrayBox&       fine,
                      int              fine_comp,
                      int              ncomp,
                      const Box&       fine_region,
                      const IntVect &  ratio,
                      const Geometry&  /* crse_geom */,
                      const Geometry&  /* fine_geom */,
                      Array<BCRec>&    /*bcr*/)
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

    Array<Real> slope(slp_len);

    int strp_len = len0*ratio[0];

    Array<Real> strip(strp_len);

    int strip_lo = ratio[0] * clo[0];
    int strip_hi = ratio[0] * chi[0];

    const Real* cdat  = crse.dataPtr(crse_comp);
    Real*       fdat  = fine.dataPtr(fine_comp);
    const int* ratioV = ratio.getVect();

    FORT_CBINTERP (cdat,ARLIM(clo),ARLIM(chi),ARLIM(clo),ARLIM(chi),
                   fdat,ARLIM(flo),ARLIM(fhi),ARLIM(lo),ARLIM(hi),
                   D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),&ncomp,
                   slope.dataPtr(),&num_slope,strip.dataPtr(),&strip_lo,&strip_hi);
}

static
Array<int>
GetBCArray (const Array<BCRec>& bcr)
{
    Array<int> bc(2*BL_SPACEDIM*bcr.size());

    for (int n = 0; n < bcr.size(); n++)
    {
        const int* b_rec = bcr[n].vect();

        for (int m = 0; m < 2*BL_SPACEDIM; m++)
        {
            bc[2*BL_SPACEDIM*n + m] = b_rec[m];
        }
    }

    return bc;
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
    Box crse = BoxLib::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellConservativeLinear::CoarseBox (const Box& fine,
                                   int        ratio)
{
    Box crse(BoxLib::coarsen(fine,ratio));
    crse.grow(1);
    return crse;
}

void
CellConservativeLinear::interp (const FArrayBox& crse,
                                int              crse_comp,
                                FArrayBox&       fine,
                                int              fine_comp,
                                int              ncomp,
                                const Box&       fine_region,
                                const IntVect&   ratio,
                                const Geometry&  crse_geom,
                                const Geometry&  fine_geom,
                                Array<BCRec>& bcr)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::interp");

    BL_ASSERT(bcr.size() >= ncomp);
    BL_ASSERT(fine_geom.Domain().contains(fine_region));
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
    // Get coarse and fine edge-centered volume coordinates.
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

    FArrayBox  cmax(cslope_bx,ncomp);
    FArrayBox  cmin(cslope_bx,ncomp);
    FArrayBox alpha(cslope_bx,ncomp);

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
    int slope_flag    = 1;

    int cvcbhi[BL_SPACEDIM];
    int fvcbhi[BL_SPACEDIM];

    for (dir=0; dir<BL_SPACEDIM; dir++)
    {
        cvcbhi[dir] = cvcblo[dir] + cvc[dir].size() - 1;
        fvcbhi[dir] = fvcblo[dir] + fvc[dir].size() - 1;
    }

    D_TERM(Real* voffx = new Real[fvc[0].size()];,
           Real* voffy = new Real[fvc[1].size()];,
           Real* voffz = new Real[fvc[2].size()];);

    Array<int> bc     = GetBCArray(bcr);
    const int* ratioV = ratio.getVect();

#if (BL_SPACEDIM > 1)

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
                      bc.dataPtr(), &slope_flag, &lin_limit,
                      D_DECL(fvc[0].dataPtr(),fvc[1].dataPtr(),fvc[2].dataPtr()),
                      D_DECL(cvc[0].dataPtr(),cvc[1].dataPtr(),cvc[2].dataPtr()),
                      D_DECL(voffx,voffy,voffz),
                      alpha.dataPtr(),cmax.dataPtr(),cmin.dataPtr());

    D_TERM(delete [] voffx;, delete [] voffy;, delete [] voffz;);

#endif /*(BL_SPACEDIM > 1)*/
}

CellQuadratic::CellQuadratic (bool limit)
{
    do_limited_slope = limit;
}

CellQuadratic::~CellQuadratic () {}

Box
CellQuadratic::CoarseBox (const Box&     fine,
                          const IntVect& ratio)
{
    Box crse = BoxLib::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellQuadratic::CoarseBox (const Box& fine,
                          int        ratio)
{
    Box crse = BoxLib::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

void
CellQuadratic::interp (const FArrayBox& crse,
                       int              crse_comp,
                       FArrayBox&       fine,
                       int              fine_comp,
                       int              ncomp,
                       const Box&       fine_region,
                       const IntVect&   ratio,
                       const Geometry&  crse_geom,
                       const Geometry&  fine_geom,
                       Array<BCRec>&    bcr)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::interp");

    BL_ASSERT(bcr.size() >= ncomp);
    BL_ASSERT(fine_geom.Domain().contains(fine_region));
    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    Box crse_bx(BoxLib::coarsen(target_fine_region,ratio));
    Box fslope_bx(BoxLib::refine(crse_bx,ratio));
    Box cslope_bx(crse_bx);
    cslope_bx.grow(1);
    BL_ASSERT(crse.box().contains(cslope_bx));
    //
    // Alloc temp space for coarse grid slopes: here we use 5 
    // instead of BL_SPACEDIM because of the x^2, y^2 and xy terms
    //
    long t_long = cslope_bx.numPts();
    BL_ASSERT(t_long < INT_MAX);
    int c_len = int(t_long);

    Array<Real> cslope(5*c_len);

    int loslp = cslope_bx.index(crse_bx.smallEnd());
    int hislp = cslope_bx.index(crse_bx.bigEnd());

    t_long = cslope_bx.numPts();
    BL_ASSERT(t_long < INT_MAX);
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

    Array<Real> strip((5+2)*f_len);

    Real* fstrip = strip.dataPtr();
    Real* foff   = fstrip + f_len;
    Real* fslope = foff + f_len;
    //
    // Get coarse and fine edge-centered volume coordinates.
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

#if (BL_SPACEDIM > 1)

    FORT_CQINTERP (fdat,ARLIM(flo),ARLIM(fhi),
                   ARLIM(fblo), ARLIM(fbhi),
                   &ncomp,D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),
                   cdat,&clo,&chi,
                   ARLIM(cblo), ARLIM(cbhi),
                   fslo,fshi,
                   cslope.dataPtr(),&c_len,fslope,fstrip,&f_len,foff,
                   bc.dataPtr(), &slope_flag,
                   D_DECL(fvc[0].dataPtr(),fvc[1].dataPtr(),fvc[2].dataPtr()),
                   D_DECL(cvc[0].dataPtr(),cvc[1].dataPtr(),cvc[2].dataPtr()));

#endif /*(BL_SPACEDIM > 1)*/
}

PCInterp::~PCInterp () {}

Box
PCInterp::CoarseBox (const Box& fine,
                     int        ratio)
{
    return BoxLib::coarsen(fine,ratio);
}

Box
PCInterp::CoarseBox (const Box&     fine,
                     const IntVect& ratio)
{
    return BoxLib::coarsen(fine,ratio);
}

void
PCInterp::interp (const FArrayBox& crse,
                  int              crse_comp,
                  FArrayBox&       fine,
                  int              fine_comp,
                  int              ncomp,
                  const Box&       fine_region,
                  const IntVect&   ratio,
                  const Geometry&  /*crse_geom*/,
                  const Geometry&  /*fine_geom*/,
                  Array<BCRec>&    /*bcr*/)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::interp");
    //
    // Set up to call FORTRAN.
    //
    const int* clo  = crse.box().loVect();
    const int* chi  = crse.box().hiVect();
    const int* flo  = fine.loVect();
    const int* fhi  = fine.hiVect();
    const int* fblo = fine_region.loVect();
    const int* fbhi = fine_region.hiVect();

    Box cregion(BoxLib::coarsen(fine_region,ratio));

    const int* cblo = cregion.loVect();
    const int* cbhi = cregion.hiVect();

    int long_dir;
    int long_len = cregion.longside(long_dir);
    int s_len    = long_len*ratio[long_dir];

    Array<Real> strip(s_len);

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
                   &ncomp,strip.dataPtr(),&strip_lo,&strip_hi);
}

CellConservativeProtected::CellConservativeProtected () {}

CellConservativeProtected::~CellConservativeProtected () {}

Box
CellConservativeProtected::CoarseBox (const Box&     fine,
                                      const IntVect& ratio)
{
    Box crse = BoxLib::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellConservativeProtected::CoarseBox (const Box& fine,
                                      int        ratio)
{
    Box crse(BoxLib::coarsen(fine,ratio));
    crse.grow(1);
    return crse;
}

void
CellConservativeProtected::interp (const FArrayBox& crse,
                                   int              crse_comp,
                                   FArrayBox&       fine,
                                   int              fine_comp,
                                   int              ncomp,
                                   const Box&       fine_region,
                                   const IntVect&   ratio,
                                   const Geometry&  crse_geom,
                                   const Geometry&  fine_geom,
                                   Array<BCRec>& bcr)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::interp");

    BL_ASSERT(bcr.size() >= ncomp);
    BL_ASSERT(fine_geom.Domain().contains(fine_region));
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
    // Get coarse and fine edge-centered volume coordinates.
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

    FArrayBox  cmax(cslope_bx,ncomp);
    FArrayBox  cmin(cslope_bx,ncomp);
    FArrayBox alpha(cslope_bx,ncomp);

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
    int lin_limit     = 1;
    const int* cvcblo = crse_bx.loVect();
    const int* fvcblo = target_fine_region.loVect();

    int cvcbhi[BL_SPACEDIM];
    int fvcbhi[BL_SPACEDIM];

    for (dir=0; dir<BL_SPACEDIM; dir++)
    {
        cvcbhi[dir] = cvcblo[dir] + cvc[dir].size() - 1;
        fvcbhi[dir] = fvcblo[dir] + fvc[dir].size() - 1;
    }

    D_TERM(Real* voffx = new Real[fvc[0].size()];,
           Real* voffy = new Real[fvc[1].size()];,
           Real* voffz = new Real[fvc[2].size()];);

    Array<int> bc     = GetBCArray(bcr);
    const int* ratioV = ratio.getVect();
    int slope_flag    = 1;

#if (BL_SPACEDIM > 1)

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
                      bc.dataPtr(), &slope_flag, &lin_limit,
                      D_DECL(fvc[0].dataPtr(),fvc[1].dataPtr(),fvc[2].dataPtr()),
                      D_DECL(cvc[0].dataPtr(),cvc[1].dataPtr(),cvc[2].dataPtr()),
                      D_DECL(voffx,voffy,voffz),
                      alpha.dataPtr(),cmax.dataPtr(),cmin.dataPtr());

    D_TERM(delete [] voffx;, delete [] voffy;, delete [] voffz;);

#endif /*(BL_SPACEDIM > 1)*/
}

void
CellConservativeProtected::protect (const FArrayBox& crse,
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
                                    Array<BCRec>& bcr)
{
    BL_ASSERT(bcr.size() >= ncomp);
    BL_ASSERT(fine_geom.Domain().contains(fine_region));

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    //
    // crse_bx is coarsening of target_fine_region, grown by 1.
    //
    Box crse_bx = CoarseBox(target_fine_region,ratio);

    //
    // cs_bx is coarsening of target_fine_region.
    //
    Box cs_bx(crse_bx);
    cs_bx.grow(-1);

    //
    // Get coarse and fine edge-centered volume coordinates.
    //
    int dir;
    Array<Real> fvc[BL_SPACEDIM];
    Array<Real> cvc[BL_SPACEDIM];
    for (dir = 0; dir < BL_SPACEDIM; dir++)
    {
        fine_geom.GetEdgeVolCoord(fvc[dir],target_fine_region,dir);
        crse_geom.GetEdgeVolCoord(cvc[dir],crse_bx,dir);
    }

#if (BL_SPACEDIM == 2)
    const int* cvcblo = crse_bx.loVect();
    const int* fvcblo = target_fine_region.loVect();

    int cvcbhi[BL_SPACEDIM];
    int fvcbhi[BL_SPACEDIM];

    for (dir=0; dir<BL_SPACEDIM; dir++)
    {
        cvcbhi[dir] = cvcblo[dir] + cvc[dir].size() - 1;
        fvcbhi[dir] = fvcblo[dir] + fvc[dir].size() - 1;
    }
#endif

    Real* fdat       = fine.dataPtr(fine_comp);
    Real* state_dat  = fine_state.dataPtr(state_comp);
    const Real* cdat = crse.dataPtr(crse_comp);
    
    const int* flo    = fine.loVect();
    const int* fhi    = fine.hiVect();
    const int* slo    = fine_state.loVect();
    const int* shi    = fine_state.hiVect();
    const int* clo    = crse.loVect();
    const int* chi    = crse.hiVect();
    const int* fblo   = target_fine_region.loVect();
    const int* fbhi   = target_fine_region.hiVect();
    const int* csbhi  = cs_bx.hiVect();
    const int* csblo  = cs_bx.loVect();

    Array<int> bc     = GetBCArray(bcr);
    const int* ratioV = ratio.getVect();

#if (BL_SPACEDIM > 1)

    FORT_PROTECT_INTERP (fdat,ARLIM(flo),ARLIM(fhi),
                         fblo, fbhi,
                         cdat,ARLIM(clo),ARLIM(chi),
                         csblo, csbhi,
#if (BL_SPACEDIM == 2)
                         fvc[0].dataPtr(),fvc[1].dataPtr(),
                         ARLIM(fvcblo), ARLIM(fvcbhi),
                         cvc[0].dataPtr(),cvc[1].dataPtr(),
                         ARLIM(cvcblo), ARLIM(cvcbhi),
#endif
                         state_dat, ARLIM(slo), ARLIM(shi),
                         &ncomp,D_DECL(&ratioV[0],&ratioV[1],&ratioV[2]),
                         bc.dataPtr());

#endif /*(BL_SPACEDIM > 1)*/
 
}
