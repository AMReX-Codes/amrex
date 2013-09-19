
#include <MultiFab.H>
#include <MLSDCAmr.H>

BEGIN_EXTERN_C

void *mf_encap_create(int type, void *encap_ctx)
{
  throw "SDCLib should not call mf_encap_create.";
  // mf_encap_t* ctx = (mf_encap_t*) encap_ctx;
  // return new MultiFab(ctx->ba, ctx->ncomp, ctx->ngrow);
}

void mf_encap_destroy(void *sol)
{
  throw "SDCLib should not call mf_encap_destroy.";
  // MultiFab* mf = (MultiFab*) sol;
  // delete mf;
}

void mf_encap_setval(void *sol, sdc_dtype_t val)
{
  MultiFab& mf = *((MultiFab*) sol);

  // XXX: does this set ghost cells too?
  mf.setVal(val);
}

void mf_encap_copy(void *dstp, void *srcp)
{
  MultiFab& dst = *((MultiFab*) dstp);
  MultiFab& src = *((MultiFab*) srcp);

  // XXX: does this copy ghost cells too?
  for (MFIter mfi(dst); mfi.isValid(); ++mfi)
    dst[mfi].copy(src[mfi]);
}

void mf_encap_saxpy(void *yp, sdc_dtype_t a, void *xp)
{
  MultiFab& y = *((MultiFab*) yp);
  MultiFab& x = *((MultiFab*) xp);

  // XXX: does this do ghost cells too?
  for (MFIter mfi(y); mfi.isValid(); ++mfi)
    y[mfi].saxpy(a, x[mfi]);
}

END_EXTERN_C


void MLSDCAmr::build_encap()
{
  encap.create  = mf_encap_create;
  encap.destroy = mf_encap_destroy;
  encap.setval  = mf_encap_setval;
  encap.copy    = mf_encap_copy;
  encap.saxpy   = mf_encap_saxpy;
  encap.ctx     = &encap_ctx;
}
