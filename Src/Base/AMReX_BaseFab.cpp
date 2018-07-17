#include <cstring>
#include <cstdlib>

#include <AMReX_BaseFab.H>
#include <AMReX_BArena.H>
#include <AMReX_CArena.H>

#if !defined(BL_NO_FORT)
#include <AMReX_BaseFab_f.H>
#endif

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

namespace amrex {

long private_total_bytes_allocated_in_fabs     = 0L;
long private_total_bytes_allocated_in_fabs_hwm = 0L;
long private_total_cells_allocated_in_fabs     = 0L;
long private_total_cells_allocated_in_fabs_hwm = 0L;

int BF_init::m_cnt = 0;

namespace
{
    Arena* the_arena = 0;
#ifdef AMREX_USE_CUDA
    Arena* the_nvar_arena = 0;
#endif
}

BF_init::BF_init ()
{
    if (m_cnt++ == 0)
    {
        BL_ASSERT(the_arena == 0);

#if defined(BL_COALESCE_FABS)
        the_arena = new CArena;
#else
        the_arena = new BArena;
#endif

#ifdef AMREX_USE_CUDA
        the_arena->SetPreferred();
#endif

#ifdef AMREX_USE_CUDA
        const std::size_t hunk_size = 64 * 1024;
        the_nvar_arena = new CArena(hunk_size);
        the_nvar_arena->SetHostAlloc();
#endif

#ifdef _OPENMP
#pragma omp parallel
	{
	    amrex::private_total_bytes_allocated_in_fabs     = 0;
	    amrex::private_total_bytes_allocated_in_fabs_hwm = 0;
	    amrex::private_total_cells_allocated_in_fabs     = 0;
	    amrex::private_total_cells_allocated_in_fabs_hwm = 0;
	}
#endif

#ifdef BL_MEM_PROFILING
	MemProfiler::add("Fab", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {amrex::TotalBytesAllocatedInFabs(),
				     amrex::TotalBytesAllocatedInFabsHWM()};
			 }));
#endif
    }
}

BF_init::~BF_init ()
{
    if (--m_cnt == 0) {
        delete the_arena;
#ifdef AMREX_USE_CUDA
        delete the_nvar_arena;
#endif
    }
}

long 
TotalBytesAllocatedInFabs()
{
#ifdef _OPENMP
    long r=0;
#pragma omp parallel reduction(+:r)
    {
	r += private_total_bytes_allocated_in_fabs;
    }
    return r;
#else
    return private_total_bytes_allocated_in_fabs;
#endif
}

long 
TotalBytesAllocatedInFabsHWM()
{
#ifdef _OPENMP
    long r=0;
#pragma omp parallel reduction(+:r)
    {
	r += private_total_bytes_allocated_in_fabs_hwm;
    }
    return r;
#else
    return private_total_bytes_allocated_in_fabs_hwm;
#endif
}

long 
TotalCellsAllocatedInFabs()
{
#ifdef _OPENMP
    long r=0;
#pragma omp parallel reduction(+:r)
    {
	r += private_total_cells_allocated_in_fabs;
    }
    return r;
#else
    return private_total_cells_allocated_in_fabs;
#endif
}

long 
TotalCellsAllocatedInFabsHWM()
{
#ifdef _OPENMP
    long r=0;
#pragma omp parallel reduction(+:r)
    {
	r += private_total_cells_allocated_in_fabs_hwm;
    }
    return r;
#else
    return private_total_cells_allocated_in_fabs_hwm;
#endif
}

void 
ResetTotalBytesAllocatedInFabsHWM()
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	private_total_bytes_allocated_in_fabs_hwm = 0;
    }
}

void
update_fab_stats (long n, long s, size_t szt)
{
    long tst = s*szt;
    amrex::private_total_bytes_allocated_in_fabs += tst;
    amrex::private_total_bytes_allocated_in_fabs_hwm 
	= std::max(amrex::private_total_bytes_allocated_in_fabs_hwm,
		   amrex::private_total_bytes_allocated_in_fabs);
	
    if(szt == sizeof(Real)) {
	amrex::private_total_cells_allocated_in_fabs += n;
	amrex::private_total_cells_allocated_in_fabs_hwm 
	    = std::max(amrex::private_total_cells_allocated_in_fabs_hwm,
		       amrex::private_total_cells_allocated_in_fabs);
    }
}

Arena*
The_Arena ()
{
    BL_ASSERT(the_arena != 0);

    return the_arena;
}

#ifdef AMREX_USE_CUDA
Arena*
The_Nvar_Arena ()
{
    BL_ASSERT(the_nvar_arena != 0);

    return the_nvar_arena;
}
#endif

#if !defined(BL_NO_FORT)
template<>
void
BaseFab<Real>::performCopy (const BaseFab<Real>& src,
                            const Box&           srcbox,
                            int                  srccomp,
                            const Box&           destbox,
                            int                  destcomp,
                            int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

#pragma gpu
    amrex_fort_fab_copy
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         BL_TO_FORTRAN_N_ANYD(src,srccomp),
         AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(destbox.loVect()),
         numcomp);
}

template <>
std::size_t
BaseFab<Real>::copyToMem (const Box& srcbox,
                          int        srccomp,
                          int        numcomp,
                          void*      dst) const
{
    BL_ASSERT(box().contains(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= nComp());

    if (srcbox.ok())
    {
#pragma gpu
	amrex_fort_fab_copytomem
            (AMREX_ARLIM_ARG(srcbox.loVect()), AMREX_ARLIM_ARG(srcbox.hiVect()),
             AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(srcbox.hiVect()),
             static_cast<Real*>(dst),
             BL_TO_FORTRAN_N_ANYD(*this,srccomp),
             numcomp);
        return sizeof(Real) * srcbox.numPts() * numcomp;
    }
    else
    {
        return 0;
    }
}

template <>
std::size_t
BaseFab<Real>::copyFromMem (const Box&  dstbox,
                            int         dstcomp,
                            int         numcomp,
                            const void* src)
{
    BL_ASSERT(box().contains(dstbox));
    BL_ASSERT(dstcomp >= 0 && dstcomp+numcomp <= nComp());

    if (dstbox.ok()) 
    {
#pragma gpu
	amrex_fort_fab_copyfrommem
            (AMREX_ARLIM_ARG(dstbox.loVect()), AMREX_ARLIM_ARG(dstbox.hiVect()),
             AMREX_ARLIM_3D(dstbox.loVect()), AMREX_ARLIM_3D(dstbox.hiVect()),
             BL_TO_FORTRAN_N_ANYD(*this,dstcomp), numcomp,
             static_cast<const Real*>(src));
        return sizeof(Real) * dstbox.numPts() * numcomp;
    }
    else
    {
        return 0;
    }
}

template<>
void
BaseFab<Real>::performSetVal (Real       val,
                              const Box& bx,
                              int        comp,
                              int        ncomp)
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);

#pragma gpu
    amrex_fort_fab_setval
        (AMREX_ARLIM_ARG(bx.loVect()), AMREX_ARLIM_ARG(bx.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,comp), ncomp,
         val);
}

template<>
BaseFab<Real>&
BaseFab<Real>::invert (Real       val,
                       const Box& bx,
                       int        comp,
                       int        ncomp)
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);

#pragma gpu
    amrex_fort_fab_invert
        (AMREX_ARLIM_ARG(bx.loVect()), AMREX_ARLIM_ARG(bx.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,comp), ncomp,
         val);
    return *this;
}

template <>
Real
BaseFab<Real>::norminfmask (const Box& bx, const BaseFab<int>& mask, int comp, int ncomp) const
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);
    
    return amrex_fort_fab_norminfmask(AMREX_ARLIM_3D(bx.loVect()), AMREX_ARLIM_3D(bx.hiVect()),
                                BL_TO_FORTRAN_ANYD(mask),
                                BL_TO_FORTRAN_N_3D(*this,comp), &ncomp);
}

template<>
Real
BaseFab<Real>::norm (const Box& bx,
                     int        p,
                     int        comp,
                     int        ncomp) const
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);
    
    Real nrm = 0.0;
    
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
    std::shared_ptr<Real> sptr = Device::create_device_pointer<Real>();
    Real* nrm_f = sptr.get();
    CudaAPICheck(cudaMemset(nrm_f, 0, sizeof(Real)));
#else
    Real* nrm_f = &nrm;
#endif
    
    if (p == 0 || p == 1)
    {
#pragma gpu
        amrex_fort_fab_norm
            (AMREX_ARLIM_ARG(bx.loVect()), AMREX_ARLIM_ARG(bx.hiVect()),
             BL_TO_FORTRAN_N_ANYD(*this,comp), ncomp,
             p, nrm_f);
    }
    else
        {
            amrex::Error("BaseFab<Real>::norm(): only p == 0 or p == 1 are supported");
        }
    
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
    CudaAPICheck(cudaMemcpy(&nrm, nrm_f, sizeof(Real), cudaMemcpyDeviceToHost));
#endif
    
    return nrm;
}

template<>
Real
BaseFab<Real>::sum (const Box& bx,
                    int        comp,
                    int        ncomp) const
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);
    
    Real sm = 0.0;
    
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
    std::shared_ptr<Real> sptr = Device::create_device_pointer<Real>();
    Real* sm_f = sptr.get();
    CudaAPICheck(cudaMemset(sm_f, 0, sizeof(Real)));
#else
    Real* sm_f = &sm;
#endif

#pragma gpu
    amrex_fort_fab_sum
        (AMREX_ARLIM_ARG(bx.loVect()), AMREX_ARLIM_ARG(bx.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,comp), ncomp, sm_f);
    
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
    CudaAPICheck(cudaMemcpy(&sm, sm_f, sizeof(Real), cudaMemcpyDeviceToHost));
#endif
    
    return sm;
}

template<>
BaseFab<Real>&
BaseFab<Real>::plus (const BaseFab<Real>& src,
                     const Box&           srcbox,
                     const Box&           destbox,
                     int                  srccomp,
                     int                  destcomp,
                     int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

#pragma gpu
    amrex_fort_fab_plus
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         BL_TO_FORTRAN_N_ANYD(src,srccomp),
         AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(destbox.loVect()),
         numcomp);

    return *this;
}

template<>
BaseFab<Real>&
BaseFab<Real>::mult (const BaseFab<Real>& src,
                     const Box&           srcbox,
                     const Box&           destbox,
                     int                  srccomp,
                     int                  destcomp,
                     int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

#pragma gpu
    amrex_fort_fab_mult
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         BL_TO_FORTRAN_N_ANYD(src,srccomp),
         AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(destbox.loVect()),
         numcomp);
    return *this;
}

template <>
BaseFab<Real>&
BaseFab<Real>::saxpy (Real a, const BaseFab<Real>& src,
                      const Box&        srcbox,
                      const Box&        destbox,
                      int               srccomp,
                      int               destcomp,
                      int               numcomp)
{
    BL_ASSERT(srcbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(destbox.ok());
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT( srccomp >= 0 &&  srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <=     nComp());

#pragma gpu
    amrex_fort_fab_saxpy
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         a,
         BL_TO_FORTRAN_N_ANYD(src,srccomp),
         AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(destbox.loVect()),
         numcomp);
    return *this;
}

template <>
BaseFab<Real>&
BaseFab<Real>::xpay (Real a, const BaseFab<Real>& src,
		     const Box&        srcbox,
		     const Box&        destbox,
		     int               srccomp,
		     int               destcomp,
		     int               numcomp)
{
    BL_ASSERT(srcbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(destbox.ok());
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT( srccomp >= 0 &&  srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <=     nComp());

#pragma gpu
    amrex_fort_fab_xpay
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         a,
         BL_TO_FORTRAN_N_ANYD(src,srccomp),
         AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(destbox.loVect()),
         numcomp);
    return *this;
}

template <>
BaseFab<Real>&
BaseFab<Real>::addproduct (const Box&           destbox,
			   int                  destcomp,
			   int                  numcomp,
			   const BaseFab<Real>& src1,
			   int                  comp1,
			   const BaseFab<Real>& src2,
			   int                  comp2)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(   comp1 >= 0 &&    comp1+numcomp <= src1.nComp());
    BL_ASSERT(   comp2 >= 0 &&    comp2+numcomp <= src2.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <=      nComp());

#pragma gpu
    amrex_fort_fab_addproduct
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         BL_TO_FORTRAN_N_ANYD(src1,comp1),
         BL_TO_FORTRAN_N_ANYD(src2,comp2),
         numcomp);
    return *this;
}

template<>
BaseFab<Real>&
BaseFab<Real>::minus (const BaseFab<Real>& src,
                      const Box&           srcbox,
                      const Box&           destbox,
                      int                  srccomp,
                      int                  destcomp,
                      int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

#pragma gpu
    amrex_fort_fab_minus
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         BL_TO_FORTRAN_N_ANYD(src,srccomp),
         AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(destbox.loVect()),
         numcomp);
    return *this;
}

template<>
BaseFab<Real>&
BaseFab<Real>::divide (const BaseFab<Real>& src,
                       const Box&           srcbox,
                       const Box&           destbox,
                       int                  srccomp,
                       int                  destcomp,
                       int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

#pragma gpu
    amrex_fort_fab_divide
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         BL_TO_FORTRAN_N_ANYD(src,srccomp),
         AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(destbox.loVect()),
         numcomp);
    return *this;
}

template<>
BaseFab<Real>&
BaseFab<Real>::protected_divide (const BaseFab<Real>& src,
                                 const Box&           srcbox,
                                 const Box&           destbox,
                                 int                  srccomp,
                                 int                  destcomp,
                                 int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

#pragma gpu
    amrex_fort_fab_protdivide
        (AMREX_ARLIM_ARG(destbox.loVect()), AMREX_ARLIM_ARG(destbox.hiVect()),
         BL_TO_FORTRAN_N_ANYD(*this,destcomp),
         BL_TO_FORTRAN_N_ANYD(src,srccomp),
         AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(destbox.loVect()),
         numcomp);
    return *this;
}

template <>
BaseFab<Real>&
BaseFab<Real>::linComb (const BaseFab<Real>& f1,
			const Box&           b1,
			int                  comp1,
			const BaseFab<Real>& f2,
			const Box&           b2,
			int                  comp2,
			Real                 alpha,
			Real                 beta,
			const Box&           b,
			int                  comp,
			int                  numcomp)
{
    BL_ASSERT(b1.ok());
    BL_ASSERT(f1.box().contains(b1));
    BL_ASSERT(b2.ok());
    BL_ASSERT(f2.box().contains(b2));
    BL_ASSERT(b.ok());
    BL_ASSERT(box().contains(b));
    BL_ASSERT(b.sameSize(b1));
    BL_ASSERT(b.sameSize(b2));
    BL_ASSERT(comp1 >= 0 && comp1+numcomp <= f1.nComp());
    BL_ASSERT(comp2 >= 0 && comp2+numcomp <= f2.nComp());
    BL_ASSERT(comp  >= 0 && comp +numcomp <=    nComp());

#pragma gpu
    amrex_fort_fab_lincomb
        (AMREX_ARLIM_ARG(b.loVect()), AMREX_ARLIM_ARG(b.hiVect()),
         AMREX_ARLIM_3D(b.loVect()),
         BL_TO_FORTRAN_N_ANYD(*this,comp),
         alpha, BL_TO_FORTRAN_N_ANYD(f1,comp1), AMREX_ARLIM_3D(b1.loVect()),
         beta,  BL_TO_FORTRAN_N_ANYD(f2,comp2), AMREX_ARLIM_3D(b2.loVect()),
         numcomp);
    return *this;
}

template <>
Real
BaseFab<Real>::dot (const Box& xbx, int xcomp, 
                    const BaseFab<Real>& y, const Box& ybx, int ycomp,
                    int numcomp) const
{
    BL_ASSERT(xbx.ok());
    BL_ASSERT(box().contains(xbx));
    BL_ASSERT(y.box().contains(ybx));
    BL_ASSERT(xbx.sameSize(ybx));
    BL_ASSERT(xcomp >= 0 && xcomp+numcomp <=   nComp());
    BL_ASSERT(ycomp >= 0 && ycomp+numcomp <= y.nComp());
    
    Real dp = 0.0;
    
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
    std::shared_ptr<Real> sptr = Device::create_device_pointer<Real>();
    Real* dp_f = sptr.get();
    CudaAPICheck(cudaMemset(dp_f, 0, sizeof(Real)));
#else
    Real* dp_f = &dp;
#endif

#pragma gpu
    amrex_fort_fab_dot
        (AMREX_ARLIM_ARG(xbx.loVect()), AMREX_ARLIM_ARG(xbx.hiVect()),
         AMREX_ARLIM_3D(xbx.loVect()),
         BL_TO_FORTRAN_N_ANYD(*this,xcomp),
         BL_TO_FORTRAN_N_ANYD(y,ycomp), AMREX_ARLIM_3D(ybx.loVect()),
         numcomp, dp_f);
    
#if (defined(AMREX_USE_CUDA) && !defined(AMREX_NO_DEVICE_LAUNCH))
    CudaAPICheck(cudaMemcpy(&dp, dp_f, sizeof(Real), cudaMemcpyDeviceToHost));
#endif
    
    return dp;
}

template <>
Real
BaseFab<Real>::dotmask (const BaseFab<int>& mask, const Box& xbx, int xcomp, 
                        const BaseFab<Real>& y, const Box& ybx, int ycomp,
                        int numcomp) const
{
    BL_ASSERT(xbx.ok());
    BL_ASSERT(box().contains(xbx));
    BL_ASSERT(y.box().contains(ybx));
    BL_ASSERT(xbx.sameSize(ybx));
    BL_ASSERT(xcomp >= 0 && xcomp+numcomp <=   nComp());
    BL_ASSERT(ycomp >= 0 && ycomp+numcomp <= y.nComp());

    return amrex_fort_fab_dot_mask(AMREX_ARLIM_3D(xbx.loVect()), AMREX_ARLIM_3D(xbx.hiVect()),
                             BL_TO_FORTRAN_N_3D(*this,xcomp),
                             BL_TO_FORTRAN_N_3D(y,ycomp), AMREX_ARLIM_3D(ybx.loVect()),
                             BL_TO_FORTRAN_ANYD(mask),
                             &numcomp);
}

template<>
void
BaseFab<int>::performCopy (const BaseFab<int>& src,
                           const Box&          srcbox,
                           int                 srccomp,
                           const Box&          destbox,
                           int                 destcomp,
                           int                 numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    amrex_fort_ifab_copy(AMREX_ARLIM_3D(destbox.loVect()), AMREX_ARLIM_3D(destbox.hiVect()),
                   BL_TO_FORTRAN_N_3D(*this,destcomp),
                   BL_TO_FORTRAN_N_3D(src,srccomp), AMREX_ARLIM_3D(srcbox.loVect()),
                   &numcomp);
}

template <>
std::size_t
BaseFab<int>::copyToMem (const Box& srcbox,
                         int        srccomp,
                         int        numcomp,
                         void*      dst) const
{
    BL_ASSERT(box().contains(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= nComp());

    if (srcbox.ok())
    {
	long nints =  amrex_fort_ifab_copytomem(AMREX_ARLIM_3D(srcbox.loVect()), AMREX_ARLIM_3D(srcbox.hiVect()),
                                          static_cast<int*>(dst),
                                          BL_TO_FORTRAN_N_3D(*this,srccomp),
                                          &numcomp);
        return sizeof(int) * nints;
    }
    else
    {
        return 0;
    }
}

template <>
std::size_t
BaseFab<int>::copyFromMem (const Box&  dstbox,
                           int         dstcomp,
                           int         numcomp,
                           const void* src)
{
    BL_ASSERT(box().contains(dstbox));
    BL_ASSERT(dstcomp >= 0 && dstcomp+numcomp <= nComp());

    if (dstbox.ok()) 
    {
	long nints = amrex_fort_ifab_copyfrommem(AMREX_ARLIM_3D(dstbox.loVect()), AMREX_ARLIM_3D(dstbox.hiVect()),
                                           BL_TO_FORTRAN_N_3D(*this,dstcomp), &numcomp,
                                           static_cast<const int*>(src));
        return sizeof(int) * nints;
    }
    else
    {
        return 0;
    }
}

template<>
void
BaseFab<int>::performSetVal (int        val,
                             const Box& bx,
                             int        comp,
                             int        ncomp)
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);

    amrex_fort_ifab_setval(AMREX_ARLIM_3D(bx.loVect()), AMREX_ARLIM_3D(bx.hiVect()),
                     BL_TO_FORTRAN_N_3D(*this,comp), &ncomp,
                     &val);
}

template<>
BaseFab<int>&
BaseFab<int>::plus (const BaseFab<int>& src,
                    const Box&           srcbox,
                    const Box&           destbox,
                    int                  srccomp,
                    int                  destcomp,
                    int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    amrex_fort_ifab_plus(AMREX_ARLIM_3D(destbox.loVect()), AMREX_ARLIM_3D(destbox.hiVect()),
                   BL_TO_FORTRAN_N_3D(*this,destcomp),
                   BL_TO_FORTRAN_N_3D(src,srccomp), AMREX_ARLIM_3D(srcbox.loVect()),
                   &numcomp);
    
    return *this;
}

template<>
BaseFab<int>&
BaseFab<int>::minus (const BaseFab<int>& src,
                     const Box&           srcbox,
                     const Box&           destbox,
                     int                  srccomp,
                     int                  destcomp,
                     int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    amrex_fort_ifab_minus(AMREX_ARLIM_3D(destbox.loVect()), AMREX_ARLIM_3D(destbox.hiVect()),
		   BL_TO_FORTRAN_N_3D(*this,destcomp),
		   BL_TO_FORTRAN_N_3D(src,srccomp), AMREX_ARLIM_3D(srcbox.loVect()),
		   &numcomp);
    return *this;
}

#endif

}
