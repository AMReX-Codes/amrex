/*
 * MLSDCAmrEncap encapsulation for SDCLib.
 *
 * Notes:
 *   - State/solution encaps are created with grow/ghost cells.
 *   - Function evaluation encaps are created without grow/ghost cells.
 *   - Integral encaps are created without grow/ghost cells.
 *
 * XXX: Note that the FEVAL encapsulations have flux registers, and
 * since we're using the IMEX sweeper, both the "explicit" feval and
 * "implicit" feval will have flux registers, but this isn't
 * necessary.  Matt should clean this up sometime.
 */

#include <MultiFab.H>
#include <MLSDCAmr.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>

#include <cassert>

BEGIN_EXTERN_C

void mf_encap_setval(void *Qptr, sdc_dtype val, const int flags);

void *mf_encap_create(int type, void *encap_ctx)
{
  MLSDCAmrEncapCtx* ctx   = (MLSDCAmrEncapCtx*) encap_ctx;
  MLSDCAmrEncap*    encap = new MLSDCAmrEncap;

  encap->amrlevel  = ctx->amrlevel;
  encap->type      = type;
  encap->fine_flux = 0;
  encap->crse_flux = 0;

  switch (type) {
  case SDC_SOLUTION:
  case SDC_WORK:
    encap->U = new MultiFab(*ctx->ba, ctx->ncomp, ctx->ngrow);
    break;
  case SDC_FEVAL:
  case SDC_INTEGRAL:
  case SDC_TAU:
    encap->U = new MultiFab(*ctx->ba, ctx->ncomp, 0);
    if (ctx->level > 0)
      encap->fine_flux = new FluxRegister(*ctx->ba, ctx->crse_ratio, ctx->level, ctx->ncomp);
    if (ctx->level < ctx->finest) {
      MLSDCAmr& amr  = *encap->amrlevel->getMLSDCAmr();
      AmrLevel& lvl = amr.getLevel(ctx->level+1);
      encap->crse_flux = new FluxRegister(lvl.boxArray(), amr.refRatio(ctx->level), lvl.Level(), ctx->ncomp);
    }
    break;
  }

  mf_encap_setval(encap, 0.0, SDC_ENCAP_ALL);
  return encap;
}

void mf_encap_destroy(void *Qptr)
{
  MLSDCAmrEncap* Q = (MLSDCAmrEncap*) Qptr;
  delete Q->U;
  if (Q->fine_flux != NULL) delete Q->fine_flux;
  if (Q->crse_flux != NULL) delete Q->crse_flux;
  delete Q;
}

void mf_encap_setval_flux(FluxRegister& dst, sdc_dtype val)
{
  for (OrientationIter face; face; ++face)
    for (FabSetIter bfsi(dst[face()]); bfsi.isValid(); ++bfsi)
      dst[face()][bfsi].setVal(val);
}

void mf_encap_setval(void *Qptr, sdc_dtype val, const int flags)
{
  MLSDCAmrEncap& Q = *((MLSDCAmrEncap*) Qptr);
  MultiFab& U = *Q.U;

  if ((flags & SDC_ENCAP_INTERIOR) && (flags & SDC_ENCAP_GHOST))
    U.setVal(val, U.nGrow());
  else
    U.setVal(val, 0);

  if (Q.fine_flux) mf_encap_setval_flux(*Q.fine_flux, val);
  if (Q.crse_flux) mf_encap_setval_flux(*Q.crse_flux, val);
}

void mf_encap_copy_flux(FluxRegister& dst, FluxRegister& src)
{
  for (OrientationIter face; face; ++face)
    for (FabSetIter bfsi(dst[face()]); bfsi.isValid(); ++bfsi)
      dst[face()][bfsi].copy(src[face()][bfsi]);
}

void mf_encap_copy(void *dstp, const void *srcp, int flags)
{
  MLSDCAmrEncap& Qdst = *((MLSDCAmrEncap*) dstp);
  MLSDCAmrEncap& Qsrc = *((MLSDCAmrEncap*) srcp);
  MultiFab& Udst = *Qdst.U;
  MultiFab& Usrc = *Qsrc.U;

  if ((flags & SDC_ENCAP_INTERIOR) && (flags & SDC_ENCAP_GHOST)) {
    int ngsrc = Usrc.nGrow();
    int ngdst = Udst.nGrow();
    int nghost = (ngdst < ngsrc) ? ngdst : ngsrc;
    MultiFab::Copy(Udst, Usrc, 0, 0, Usrc.nComp(), nghost);
  } else {
    MultiFab::Copy(Udst, Usrc, 0, 0, Usrc.nComp(), 0);
  }

  if (Qdst.fine_flux && Qsrc.fine_flux) mf_encap_copy_flux(*Qdst.fine_flux, *Qsrc.fine_flux);
  if (Qdst.crse_flux && Qsrc.crse_flux) mf_encap_copy_flux(*Qdst.crse_flux, *Qsrc.crse_flux);

#ifndef NDEBUG
  BL_ASSERT(Usrc.contains_nan() == false);
  BL_ASSERT(Udst.contains_nan() == false);
#endif
}

void mf_encap_saxpy_flux(FluxRegister& y, sdc_dtype a, FluxRegister& x)
{
  for (OrientationIter face; face; ++face)
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FabSetIter bfsi(y[face()]); bfsi.isValid(); ++bfsi)
      y[face()][bfsi].saxpy(a, x[face()][bfsi]);
}

void mf_encap_saxpy(void *yp, sdc_dtype a, void *xp, int flags)
{
  MLSDCAmrEncap& Qy = *((MLSDCAmrEncap*) yp);
  MLSDCAmrEncap& Qx = *((MLSDCAmrEncap*) xp);
  MultiFab& Uy = *Qy.U;
  MultiFab& Ux = *Qx.U;
  int ncomp = Uy.nComp();
  int ngrow = std::min(Ux.nGrow(), Uy.nGrow());

  BL_ASSERT(Uy.boxArray() == Ux.boxArray());

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(Uy,true); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.growntilebox(ngrow);
      Uy[mfi].saxpy(a, Ux[mfi], bx, bx, 0, 0, ncomp);
  }

  if ((Qy.type==SDC_TAU) && (Qx.fine_flux!=NULL)) mf_encap_saxpy_flux(*Qy.fine_flux, a, *Qx.fine_flux);
  if ((Qy.type==SDC_TAU) && (Qx.crse_flux!=NULL)) mf_encap_saxpy_flux(*Qy.crse_flux, a, *Qx.crse_flux);
}

void mf_encap_norm(void *qp, double *n)
{
  MLSDCAmrEncap& Q = *((MLSDCAmrEncap*) qp);
  MultiFab& U = *Q.U;

  double m = 0.0;
  for (int c=0; c<U.nComp(); c++)
    m += U.norm0(c);
  *n = m / U.nComp();
}

END_EXTERN_C

sdc_encap* MLSDCAmr::BuildEncap(int lev)
{
  const DescriptorList& dl = getLevel(lev).get_desc_lst();
  assert(dl.size() == 1);

  MLSDCAmrEncapCtx* ctx = new MLSDCAmrEncapCtx;
  ctx->level    = lev;
  ctx->ba       = &boxArray(lev);
  ctx->amrlevel = dynamic_cast<MLSDCAmrLevel*>(&getLevel(lev));
  ctx->finest   = finest_level;
  ctx->ncomp    = dl[0].nComp();
  ctx->ngrow    = dl[0].nExtra();
  if (lev > 0)
    ctx->crse_ratio = refRatio(lev-1);

  sdc_encap* encap = new sdc_encap;
  encap->create  = mf_encap_create;
  encap->destroy = mf_encap_destroy;
  encap->setval  = mf_encap_setval;
  encap->copy    = mf_encap_copy;
  encap->saxpy   = mf_encap_saxpy;
  encap->norm    = mf_encap_norm;
  encap->ctx     = ctx;

  return encap;
}
