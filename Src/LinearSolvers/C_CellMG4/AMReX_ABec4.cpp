#include <algorithm>
#include <AMReX_ABec4.H>
#include <AMReX_ABec4_F.H>
#include <AMReX_ParallelDescriptor.H>

#include <AMReX_LO_BCTYPES.H>
#include <AMReX_LO_F.H>

namespace amrex {

Real ABec4::a_def     = 0.0;
Real ABec4::b_def     = 1.0;
Real ABec4::alpha_def = 1.0;
Real ABec4::beta_def  = 1.0;

ABec4::ABec4 (const BndryData& _bd,
	      const Real*      _h)
    :
    LinOp(_bd,_h)
{
    LO_Op = new ABec2(_bd,_h);

    buildWorkSpace();

    initCoefficients(_bd.boxes(),_bd.DistributionMap());
}

ABec4::~ABec4 ()
{
    clearToLevel(-1);
    delete LO_Op; LO_Op = 0;
}

void
ABec4::buildWorkSpace()
{
  const BoxArray& ba = boxArray();
  const DistributionMapping& dm = DistributionMapping();
  BL_ASSERT(resL.size()==0);
  resL.define(ba, dm, 1, 0);
}

int
ABec4::NumGrow(int level) const
{
  if (level == 0) {
    return 2;
  }
  
  BL_ASSERT(LO_Op != 0);
  return LO_Op->NumGrow(level);
}

int
ABec4::numLevels () const
{
  BL_ASSERT(LO_Op != 0);
  return LO_Op->numLevels();
}

int
ABec4::numLevelsHO () const
{
  return acoefs.size();
}

const BoxArray&
ABec4::boxArray (int level) const
{
  if (level == 0) {
    return gbox[level];
  }
  
  BL_ASSERT(LO_Op != 0);
  return LO_Op->boxArray(level);
}

void
ABec4::applyBC (MultiFab&     inout,
                int            src_comp,
                int            num_comp,
                int            level,
                LinOp::BC_Mode bc_mode,
                bool           local,
		int            bndry_comp)
{
    //
    // The inout MultiFab needs enough ghost cells for applyBC.
    //
    BL_ASSERT(inout.nGrow() >= NumGrow(level));
    //
    // No coarsened boundary values, cannot apply inhomog at lev>0.
    //
    BL_ASSERT(level < numLevelsHO());
    BL_ASSERT(!(level > 0 && bc_mode == Inhomogeneous_BC));

    int flagden = 1; // Fill in undrrelxr.
    int flagbc  = 1; // Fill boundary data.

    if (bc_mode == LinOp::Homogeneous_BC)
        //
        // No data if homogeneous.
        //
        flagbc = 0;

    prepareForLevel(level);

    const bool cross = false;
    inout.FillBoundary(src_comp,num_comp,geomarray[level].periodicity(),cross);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(inout); mfi.isValid(); ++mfi)
    {
        const int gn = mfi.index();

        BL_ASSERT(gbox[level][gn] == inout.box(gn));

        BL_ASSERT(level<undrrelxr.size());

        const BndryData::RealTuple&      bdl = bgb->bndryLocs(gn);
        const Vector< Vector<BoundCond> >& bdc = bgb->bndryConds(gn);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation o = oitr();

            FabSet&       f   = undrrelxr[level][o];
            int           cdr = o;
            const FabSet& fs  = bgb->bndryValues(o);
            const Mask&   m   = local ? lmaskvals[level][o][mfi] : maskvals[level][o][mfi];
            Real          bcl = bdl[o];
            BL_ASSERT(bdc[o].size()>bndry_comp);
            int           bct = bdc[o][bndry_comp];

            const Box&       vbx   = inout.box(gn);
            FArrayBox&       iofab = inout[gn];
            BL_ASSERT(f.size()>gn);
            BL_ASSERT(fs.size()>gn);

            FArrayBox&       ffab  = f[gn];
            const FArrayBox& fsfab = fs[gn];

            amrex_ab4_applybc4(&flagden, &flagbc, &maxorder,
			  iofab.dataPtr(src_comp),
			  ARLIM(iofab.loVect()), ARLIM(iofab.hiVect()),
			  &cdr, &bct, &bcl,
			  fsfab.dataPtr(bndry_comp), 
			  ARLIM(fsfab.loVect()), ARLIM(fsfab.hiVect()),
			  m.dataPtr(),
			  ARLIM(m.loVect()), ARLIM(m.hiVect()),
			  ffab.dataPtr(),
			  ARLIM(ffab.loVect()), ARLIM(ffab.hiVect()),
			  vbx.loVect(),
			  vbx.hiVect(), &num_comp, h[level].data());
        }
    }

  // Clean up corners:
  // The problem here is that APPLYBC fills only grow cells normal to the boundary.
  // As a result, any corner cell on the boundary (either coarse-fine or fine-fine)
  // is not filled.  For coarse-fine, the operator adjusts itself, sliding away from
  // the box edge to avoid referencing that corner point.  On the physical boundary
  // though, the corner point is needed.  Particularly if a fine-fine boundary intersects
  // the physical boundary, since we want the stencil to be independent of the box
  // blocking. 
    inout.EnforcePeriodicity(geomarray[level].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(inout); mfi.isValid(); ++mfi) {

    const int gn = mfi.index();

    BL_ASSERT(gbox[level][gn] == inout.box(gn));

    const Box& vbx   = inout.box(gn);
    FArrayBox& iofab = inout[gn];

    amrex_ab4_applybc4_touchup(
      iofab.dataPtr(src_comp),ARLIM(iofab.loVect()), ARLIM(iofab.hiVect()),
      vbx.loVect(), vbx.hiVect(), &num_comp);
  }
}

void
ABec4::ca2cc(const MultiFab& ca, MultiFab& cc,
             int sComp, int dComp, int nComp)
{
  const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(ca,tiling); mfi.isValid(); ++mfi) {
    const FArrayBox& caf = ca[mfi];
    FArrayBox& ccf = cc[mfi];
    const Box& box = mfi.tilebox();
    BL_ASSERT(caf.box().contains(amrex::grow(box,1)));
    amrex_ab4_ca2cc(box.loVect(), box.hiVect(),
               caf.dataPtr(sComp), ARLIM(caf.box().loVect()), ARLIM(caf.box().hiVect()),
               ccf.dataPtr(dComp), ARLIM(ccf.box().loVect()), ARLIM(ccf.box().hiVect()),
               &nComp);
  }
}

void
ABec4::cc2ca(const MultiFab& cc, MultiFab& ca,
             int sComp, int dComp, int nComp)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(ca,true); mfi.isValid(); ++mfi) {
    const FArrayBox& ccf = cc[mfi];
    FArrayBox& caf = ca[mfi];
    const Box& box = mfi.growntilebox();
    BL_ASSERT(ccf.box().contains(amrex::grow(box,1)));
    amrex_ab4_cc2ca(box.loVect(), box.hiVect(),
               ccf.dataPtr(sComp), ARLIM(ccf.box().loVect()), ARLIM(ccf.box().hiVect()),
               caf.dataPtr(dComp), ARLIM(caf.box().loVect()), ARLIM(caf.box().hiVect()),
               &nComp);
  }
}

void
ABec4::lo_cc2ec(const MultiFab& cc, MultiFab& ec,
                int sComp, int dComp, int nComp, int dir, bool do_harm)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(ec,true); mfi.isValid(); ++mfi) {
    const FArrayBox& ccf = cc[mfi];
    FArrayBox& ecf = ec[mfi];
    const Box& box = mfi.growntilebox();
    BL_ASSERT(ccf.box().contains(Box(box).enclosedCells().grow(dir,1)));
    int iharm = (int)do_harm;
    amrex_ab4_lo_cc2ec(box.loVect(), box.hiVect(),
                  ccf.dataPtr(sComp), ARLIM(ccf.box().loVect()), ARLIM(ccf.box().hiVect()),
                  ecf.dataPtr(dComp), ARLIM(ecf.box().loVect()), ARLIM(ecf.box().hiVect()),
                  &nComp,&dir,&iharm);
  }
}

Real
ABec4::norm (int nm, int level, const bool local)
{
  return 1;
}

void
ABec4::clearToLevel (int level)
{
  BL_ASSERT(level >= -1);
  
  for (int i = level+1; i < numLevelsHO(); ++i)
  {
    BL_ASSERT(acoefs.size()>=i);
    delete acoefs[i];
    a_valid[i] = false;
    BL_ASSERT(bcoefs.size()>=i);
    delete bcoefs[i];
    b_valid[i] = false;
  }

  BL_ASSERT(LO_Op != 0);
  LO_Op->clearToLevel(level);
}

void
ABec4::prepareForLevel (int level)
{
    LinOp::prepareForLevel(level);

    BL_ASSERT(LO_Op != 0);
    LO_Op->prepareForLevel(level);

    if (level == 0 )
        return;

    LO_Op->prepareForLevel(level-1);
}

void
ABec4::initCoefficients (const BoxArray& _ba, const DistributionMapping& _dm)
{
    const int nComp=1;
    const int nGrow=2;
    acoefs.resize(1);
    bcoefs.resize(1);
    acoefs[0] = new MultiFab(_ba, _dm, nComp, nGrow);
    acoefs[0]->setVal(a_def);
    a_valid.resize(1);
    a_valid[0] = true;

    const int nGrowb=2;
    bcoefs[0] = new MultiFab(_ba, _dm, nComp, nGrowb);
    bcoefs[0]->setVal(b_def);
    b_valid.resize(1);
    b_valid[0] = true;
}

void
ABec4::aCoefficients (const MultiFab& _a)
{
    BL_ASSERT(_a.ok());
    BL_ASSERT(_a.boxArray() == (acoefs[0])->boxArray());
    invalidate_a_to_level(0);
    MultiFab::Copy(*acoefs[0],_a,0,0,acoefs[0]->nComp(),acoefs[0]->nGrow());
}

void
ABec4::ZeroACoefficients ()
{
    invalidate_a_to_level(0);
    (*acoefs[0]).setVal(0,0,acoefs[0]->nComp(),acoefs[0]->nGrow());
}

void
ABec4::bCoefficients (const MultiFab& _b)
{
    BL_ASSERT(_b.ok());
    BL_ASSERT(_b.boxArray() == (bcoefs[0])->boxArray());
    invalidate_b_to_level(0);
    MultiFab::Copy(*bcoefs[0],_b,0,0,bcoefs[0]->nComp(),bcoefs[0]->nGrow());
}

void
ABec4::bCoefficients (const FArrayBox& _b,
		      int              gridno)
{
    BL_ASSERT(_b.box().contains((bcoefs[0])->boxArray()[gridno]));
    invalidate_b_to_level(0);
    (*bcoefs[0])[gridno].copy(_b,0,0,bcoefs[0]->nComp());
}

const MultiFab&
ABec4::aCoefficients (int level)
{
    prepareForLevel(level);
    return *acoefs[level];
}

const MultiFab&
ABec4::bCoefficients (int level)
{
    prepareForLevel(level);
    return *bcoefs[level];
}

void
ABec4::setCoefficients (const MultiFab &_a,
			const MultiFab &_b)
{
  aCoefficients(_a);
  bCoefficients(_b);

  if (LO_Op) {
    int level = 0;
    const BoxArray& cba = boxArray(level);
    const DistributionMapping& dm = DistributionMap();
    LO_Op->aCoefficients(_a);
    bool do_harm = true;
    for (int d=0; d<BL_SPACEDIM; ++d) {
      BoxArray eba = BoxArray(cba).surroundingNodes(d);
      MultiFab btmp(eba,dm,1,0);
      lo_cc2ec(_b,btmp,0,0,1,d,do_harm);
      LO_Op->bCoefficients(btmp,d);
    }
  }
}

void
ABec4::invalidate_a_to_level (int lev)
{
    lev = (lev >= 0 ? lev : 0);
    for (int i = lev; i < numLevelsHO(); i++) {
        a_valid[i] = false;
    }
    LO_Op->invalidate_a_to_level(lev);
}

void
ABec4::invalidate_b_to_level (int lev)
{
    lev = (lev >= 0 ? lev : 0);
    for (int i = lev; i < numLevelsHO(); i++) {
        b_valid[i] = false;
    }
    LO_Op->invalidate_b_to_level(lev);
}

void
ABec4::compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
		 MultiFab& in, const BC_Mode& bc_mode,
		 int src_comp, int dst_comp, int num_comp, int bnd_comp)
{
  compFlux(AMREX_D_DECL(xflux, yflux, zflux), in, true, bc_mode, src_comp, dst_comp, num_comp, bnd_comp);
}

void
ABec4::compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
		 MultiFab& in, bool do_ApplyBC, const BC_Mode& bc_mode,
		 int src_comp, int dst_comp, int num_comp, int bnd_comp)
{
    const int level = 0;
    BL_ASSERT(num_comp==1);

    if (do_ApplyBC)
      applyBC(in,src_comp,num_comp,level,bc_mode,bnd_comp);

    const MultiFab& a = aCoefficients(level);
    const MultiFab& b = bCoefficients(level);

    const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter inmfi(in,tiling); inmfi.isValid(); ++inmfi)
    {
        AMREX_D_TERM(const Box& xbx   = inmfi.nodaltilebox(0);,
	       const Box& ybx   = inmfi.nodaltilebox(1);,
	       const Box& zbx   = inmfi.nodaltilebox(2););

        FArrayBox& infab = in[inmfi];
	const FArrayBox& bfab = b[inmfi];

        AMREX_D_TERM(FArrayBox& xfluxfab = xflux[inmfi];,
               FArrayBox& yfluxfab = yflux[inmfi];,
               FArrayBox& zfluxfab = zflux[inmfi];);

        amrex_ab4_flux(infab.dataPtr(src_comp),
		  ARLIM(infab.loVect()), ARLIM(infab.hiVect()),
		  &alpha, &beta, a[inmfi].dataPtr(), 
		  ARLIM(a[inmfi].loVect()), ARLIM(a[inmfi].hiVect()),
		  bfab.dataPtr(), 
		  ARLIM(bfab.loVect()), ARLIM(bfab.hiVect()),
		  &num_comp,
		  h[level].data(),
		  xbx.loVect(), xbx.hiVect(), 
		  xfluxfab.dataPtr(dst_comp),
		  ARLIM(xfluxfab.loVect()), ARLIM(xfluxfab.hiVect())
#if (BL_SPACEDIM >= 2)
		  ,ybx.loVect(), ybx.hiVect()
		  ,yfluxfab.dataPtr(dst_comp),
		  ARLIM(yfluxfab.loVect()), ARLIM(yfluxfab.hiVect())
#endif
#if (BL_SPACEDIM == 3)
		  ,zbx.loVect(), zbx.hiVect()
		  ,zfluxfab.dataPtr(dst_comp),
		  ARLIM(zfluxfab.loVect()), ARLIM(zfluxfab.hiVect())
#endif
		  );
    }
}
        
//
// Must be defined for MultiGrid/CGSolver to work.
//

void
ABec4::Fsmooth (MultiFab&       solnL,
		const MultiFab& rhsL,
		int             level,
		int             redBlackFlag)
{
  amrex::Abort("ABec4 does not surrport Fsmooth at level = 0");
}

void
ABec4::Fsmooth_jacobi (MultiFab&       solnL,
		       const MultiFab& rhsL,
		       int             level)
{
  amrex::Abort("ABec4 does not surrport Fsmooth_jacobi");
}

void
ABec4::smooth (MultiFab&       solnL,
               const MultiFab& rhsL,
               int             level,
               LinOp::BC_Mode  bc_mode)
{
  BL_ASSERT(LO_Op != 0);

  if (level == 0)
  {
    bool local = false;
    for (int redBlackFlag = 0; redBlackFlag < 2; redBlackFlag++)
    {
      residual(resL,rhsL,solnL,level,bc_mode,local);
      LO_Op->altSmooth(solnL, resL, level, redBlackFlag);
    }
  }
  else {
    LO_Op->smooth(solnL,rhsL,level,bc_mode);
  }
}

void
ABec4::jacobi_smooth (MultiFab&       solnL,
                      const MultiFab& rhsL,
                      int             level,
                      LinOp::BC_Mode  bc_mode)
{
    BL_ASSERT(LO_Op != 0);

    LO_Op->jacobi_smooth(solnL,rhsL,level,bc_mode);
}

void
ABec4::Fapply (MultiFab&       y,
	       const MultiFab& x,
	       int             level)
{
  int num_comp = 1;
  int src_comp = 0;
  int dst_comp = 0;
  Fapply(y,dst_comp,x,src_comp,num_comp,level);
}

void
ABec4::Fapply (MultiFab&       y,
	       int             dst_comp,
	       const MultiFab& x,
	       int             src_comp,
	       int             num_comp,
	       int             level)
{
  if (level == 0) {

    BL_ASSERT(y.nComp()>=dst_comp+num_comp);
    BL_ASSERT(x.nComp()>=src_comp+num_comp);

    const MultiFab& a = aCoefficients(level);
    const MultiFab& b = bCoefficients(level);

    prepareForLevel(level);

    const bool cross = false;
    const_cast<MultiFab&>(b).FillBoundary(src_comp,num_comp,geomarray[level].periodicity(),cross);

    const bool tiling = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter ymfi(y,tiling); ymfi.isValid(); ++ymfi)
    {
        const Box&       tbx  = ymfi.tilebox();
        FArrayBox&       yfab = y[ymfi];
        const FArrayBox& xfab = x[ymfi];
        const FArrayBox& afab = a[ymfi];
	const FArrayBox& bfab = b[ymfi];

        amrex_ab4_adotx(yfab.dataPtr(dst_comp),
                   ARLIM(yfab.loVect()),ARLIM(yfab.hiVect()),
                   xfab.dataPtr(src_comp),
                   ARLIM(xfab.loVect()), ARLIM(xfab.hiVect()),
                   &alpha, &beta, afab.dataPtr(), 
                   ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
                   bfab.dataPtr(), 
                   ARLIM(bfab.loVect()), ARLIM(bfab.hiVect()),
                   tbx.loVect(), tbx.hiVect(), &num_comp,
                   h[level].data());
    }
  }
  else {
    amrex::Abort("ABec4 cannot do Fapply on level != 0");
  }
}

void
ABec4::apply (MultiFab&      out,
              MultiFab&      in,
              int            level,
              LinOp::BC_Mode bc_mode,
              bool           local,
	      int            src_comp,
	      int            dst_comp,
	      int            num_comp,
	      int            bndry_comp)
{
  if (level == 0) {
    applyBC(in,src_comp,num_comp,level,bc_mode,local,bndry_comp);
    Fapply(out,dst_comp,in,src_comp,num_comp,level);
  }
  else {

    BL_ASSERT(LO_Op != 0);

    LO_Op->apply(out,in,level,bc_mode,local,src_comp,dst_comp,num_comp,bndry_comp);
  }
}

void
ABec4::residual (MultiFab&       residL,
                 const MultiFab& rhsL,
                 MultiFab&       solnL,
                 int             level,
                 LinOp::BC_Mode  bc_mode,
                 bool            local)
{
  if (level == 0) {
      apply(residL, solnL, level, bc_mode, local);
      MultiFab::Xpay(residL, -1.0, rhsL, 0, 0, residL.nComp(), 0);
  }
  else {
      BL_ASSERT(LO_Op != 0);
      LO_Op->residual(residL,rhsL,solnL,level,bc_mode,local);
  }
}

}
