
#include <AMReX_BArena.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_FluxReg_C.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_iMultiFab.H>

#ifndef BL_NO_FORT
#include <AMReX_FLUXREG_F.H>
#endif

#include <vector>

namespace amrex {

FluxRegister::FluxRegister ()
{
    fine_level = ncomp = -1;
    ratio = IntVect::TheUnitVector();
    ratio.scale(-1);
}

FluxRegister::FluxRegister (const BoxArray&            fine_boxes, 
                            const DistributionMapping& dm,
                            const IntVect&             ref_ratio,
                            int                        fine_lev,
                            int                        nvar)
{
    define(fine_boxes,dm,ref_ratio,fine_lev,nvar);
}

const IntVect&
FluxRegister::refRatio () const noexcept
{
    return ratio;
}

int
FluxRegister::fineLevel () const noexcept
{
    return fine_level;
}

int
FluxRegister::crseLevel () const noexcept
{
    return fine_level-1;
}

int
FluxRegister::nComp () const noexcept
{
    return ncomp;
}

const BoxArray&
FluxRegister::coarsenedBoxes () const noexcept
{
    return grids;
}

void
FluxRegister::define (const BoxArray&            fine_boxes, 
                      const DistributionMapping& dm,
                      const IntVect&             ref_ratio,
                      int                        fine_lev,
                      int                        nvar)
{
    BL_ASSERT(fine_boxes.isDisjoint());
    BL_ASSERT(grids.size() == 0);

    ratio      = ref_ratio;
    fine_level = fine_lev;
    ncomp      = nvar;

    grids = fine_boxes;
    grids.coarsen(ratio);

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        const Orientation lo_face(dir,Orientation::low);
        const Orientation hi_face(dir,Orientation::high);

        IndexType typ(IndexType::TheCellType());

        typ.setType(dir,IndexType::NODE);

        BndryRegister::define(lo_face,typ,0,1,0,nvar,dm);
        BndryRegister::define(hi_face,typ,0,1,0,nvar,dm);
    }
}

void
FluxRegister::clear ()
{
    BndryRegister::clear();
}

FluxRegister::~FluxRegister () {}

Real
FluxRegister::SumReg (int comp) const
{
    Real sum = 0.0;

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        const FabSet& lofabs = bndry[Orientation(dir,Orientation::low) ];
        const FabSet& hifabs = bndry[Orientation(dir,Orientation::high)];

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            sum += amrex::ReduceSum(lofabs.m_mf, 0,
                   [=] AMREX_GPU_HOST_DEVICE (Box const& tbx,
                                              Array4<Real const> const& lofab) -> Real
                   {
                       Real r = 0.0;
                       AMREX_LOOP_3D(tbx, i, j, k,
                       {
                           r += lofab(i,j,k,comp);
                       });
                       return r;
                   });
            sum += amrex::ReduceSum(hifabs.m_mf, 0,
                   [=] AMREX_GPU_HOST_DEVICE (Box const& tbx,
                                              Array4<Real const> const& hifab) -> Real
                   {
                       Real r = 0.0;
                       AMREX_LOOP_3D(tbx, i, j, k,
                       {
                           r -= hifab(i,j,k,comp);
                       });
                       return r;
                   });
        }
        else
#endif
        {
#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
            for (FabSetIter fsi(lofabs); fsi.isValid(); ++fsi)
            {
                Array4<Real const> const& lofab = lofabs.const_array(fsi);
                Box lobx(lofab);
                AMREX_LOOP_3D(lobx, i, j, k,
                {
                    sum += lofab(i,j,k,comp);
                });
                Array4<Real const> const& hifab = hifabs.const_array(fsi);
                Box hibx(hifab);
                AMREX_LOOP_3D(hibx, i, j, k,
                {
                    sum -= hifab(i,j,k,comp);
                });
            }
        }
    }

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        const MultiFab& area,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult,
                        FrOp            op)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Orientation face_lo(dir,Orientation::low);
    const Orientation face_hi(dir,Orientation::high);
 
    MultiFab mf(mflx.boxArray(),mflx.DistributionMap(),numcomp,0,
                MFInfo(), mflx.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif    
    for (MFIter mfi(mflx,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();
        auto       dfab =   mf.array(mfi);
        auto const sfab = mflx.const_array(mfi);
        auto const afab = area.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, numcomp, i, j, k, n,
        {
            dfab(i,j,k,n) = sfab(i,j,k,n+srccomp)*mult*afab(i,j,k);
        });
    }

    for (int pass = 0; pass < 2; pass++)
    {
        const Orientation face = ((pass == 0) ? face_lo : face_hi);

        if (op == FluxRegister::COPY)
        {
            bndry[face].copyFrom(mf,0,0,destcomp,numcomp);
        }
        else
        {
            FabSet fs(bndry[face].boxArray(),bndry[face].DistributionMap(),numcomp);

            fs.setVal(0);

            fs.copyFrom(mf,0,0,0,numcomp);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (FabSetIter mfi(fs); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                auto const sfab = fs.const_array(mfi);
                auto       dfab = bndry[face].array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (bx, numcomp, i, j, k, n,
                {
                    dfab(i,j,k,n+destcomp) += sfab(i,j,k,n);
                });
            }
        }
    }
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult,
                        FrOp            op)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    MultiFab area(mflx.boxArray(), mflx.DistributionMap(), 1, 0,
                  MFInfo(), mflx.Factory());

    area.setVal(1, 0, 1, 0);

    CrseInit(mflx,area,dir,srccomp,destcomp,numcomp,mult,op);
}

void
FluxRegister::CrseAdd (const MultiFab& mflx,
                       const MultiFab& area,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult,
                       const Geometry& geom)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Orientation face_lo(dir,Orientation::low);
    const Orientation face_hi(dir,Orientation::high);
 
    MultiFab mf(mflx.boxArray(),mflx.DistributionMap(),numcomp,0,
                MFInfo(), mflx.Factory());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mflx,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();
        auto       dfab =   mf.array(mfi);
        auto const sfab = mflx.const_array(mfi);
        auto const afab = area.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, numcomp, i, j, k, n,
        {
            dfab(i,j,k,n) = sfab(i,j,k,n+srccomp)*mult*afab(i,j,k);
        });
    }

    for (int pass = 0; pass < 2; pass++)
    {
        const Orientation face = ((pass == 0) ? face_lo : face_hi);
        bndry[face].plusFrom(mf,0,0,destcomp,numcomp,geom.periodicity());
    }
}

void
FluxRegister::CrseAdd (const MultiFab& mflx,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult,
                       const Geometry& geom)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    MultiFab area(mflx.boxArray(), mflx.DistributionMap(), 1, 0,
                  MFInfo(), mflx.Factory());

    area.setVal(1, 0, 1, 0);

    CrseAdd(mflx,area,dir,srccomp,destcomp,numcomp,mult,geom);
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],dir,k,srccomp,destcomp,numcomp,mult,RunOn::Gpu);
    }
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       const MultiFab& area,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],area[mfi],dir,k,srccomp,destcomp,numcomp,mult,RunOn::Gpu);
    }
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult,
                       RunOn            runon) noexcept
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];
    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];
    const Box& lobox = loreg.box();
    const Box& hibox = hireg.box();

    Array4<Real> loarr = loreg.array();
    Array4<Real> hiarr = hireg.array();
    Array4<Real const> farr = flux.const_array();
    const Dim3 local_ratio = ratio.dim3();

    if ((runon == RunOn::Gpu) && Gpu::inLaunchRegion())
    {
        AMREX_LAUNCH_DEVICE_LAMBDA
        ( lobox, tlobx,
          {
              fluxreg_fineadd(tlobx, loarr, destcomp, farr, srccomp,
                              numcomp, dir, local_ratio, mult);
          },
          hibox, thibx,
          {
              fluxreg_fineadd(thibx, hiarr, destcomp, farr, srccomp,
                              numcomp, dir, local_ratio, mult);
          }
        );
    }
    else
    {
        fluxreg_fineadd(lobox, loarr, destcomp, farr, srccomp,
                        numcomp, dir, local_ratio, mult);
        fluxreg_fineadd(hibox, hiarr, destcomp, farr, srccomp,
                        numcomp, dir, local_ratio, mult);
    }
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       const FArrayBox& area,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult,
                       RunOn            runon) noexcept
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];
    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];
    const Box& lobox = loreg.box();
    const Box& hibox = hireg.box();

    Array4<Real> loarr = loreg.array();
    Array4<Real> hiarr = hireg.array();
    Array4<Real const> farr = flux.const_array();
    Array4<Real const> aarr = area.const_array();
    const Dim3 local_ratio = ratio.dim3();

    if ((runon == RunOn::Gpu) && Gpu::inLaunchRegion())
    {
        AMREX_LAUNCH_DEVICE_LAMBDA
        ( lobox, tlobx,
          {
              fluxreg_fineareaadd(tlobx, loarr, destcomp,
                                  aarr, farr, srccomp,
                                  numcomp, dir, local_ratio, mult);
          },
          hibox, thibx,
          {
              fluxreg_fineareaadd(thibx, hiarr, destcomp,
                                  aarr, farr, srccomp,
                                  numcomp, dir, local_ratio, mult);
          }
        );
    }
    else
    {
        fluxreg_fineareaadd(lobox, loarr, destcomp,
                            aarr, farr, srccomp,
                            numcomp, dir, local_ratio, mult);
        fluxreg_fineareaadd(hibox, hiarr, destcomp,
                            aarr, farr, srccomp,
                            numcomp, dir, local_ratio, mult);
    }
}

void
FluxRegister::FineSetVal (int              dir,
                          int              boxno,
                          int              destcomp,
                          int              numcomp,
                          Real             val,
                          RunOn           /*runon*/) noexcept
{
    Gpu::LaunchSafeGuard lsg(false); // xxxxx gpu todo

    // This routine used by FLASH does NOT run on gpu for safety.

    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];
    BL_ASSERT(numcomp <= loreg.nComp());
    loreg.setVal<RunOn::Host>(val, loreg.box(), destcomp, numcomp);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];
    BL_ASSERT(numcomp <= hireg.nComp());
    hireg.setVal<RunOn::Host>(val, hireg.box(), destcomp, numcomp);
}

void
FluxRegister::Reflux (MultiFab&       mf,
		      const MultiFab& volume,
		      Real            scale,
		      int             scomp,
		      int             dcomp,
		      int             nc,
		      const Geometry& geom)
{
    for (OrientationIter fi; fi; ++fi)
    {
	const Orientation& face = fi();
        Reflux(mf, volume, face, scale, scomp, dcomp, nc, geom);
    }
}

void
FluxRegister::Reflux (MultiFab&       mf,
                      const MultiFab& volume,
                      int             dir,
		      Real            scale,
		      int             scomp,
		      int             dcomp,
		      int             nc,
		      const Geometry& geom)
{
    for (int s = 0; s < 2; ++s)
    {
        Orientation::Side side = (s==0) ? Orientation::low : Orientation::high;
	Orientation face(dir, side);
        Reflux(mf, volume, face, scale, scomp, dcomp, nc, geom);
    }
}

void
FluxRegister::Reflux (MultiFab&       mf,
		      Real            scale,
		      int             scomp,
		      int             dcomp,
		      int             nc,
		      const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    MultiFab volume(mf.boxArray(), mf.DistributionMap(), 1, 0,
                    MFInfo(), mf.Factory());

    volume.setVal(AMREX_D_TERM(dx[0],*dx[1],*dx[2]), 0, 1, 0);

    Reflux(mf,volume,scale,scomp,dcomp,nc,geom);
}

void 
FluxRegister::Reflux (MultiFab&       mf,
                      int             dir,
		      Real            scale,
		      int             scomp,
		      int             dcomp,
		      int             nc,
		      const Geometry& geom)
{
    const Real* dx = geom.CellSize();
    
    MultiFab volume(mf.boxArray(), mf.DistributionMap(), 1, 0,
                    MFInfo(), mf.Factory());
    
    volume.setVal(AMREX_D_TERM(dx[0],*dx[1],*dx[2]), 0, 1, 0);

    Reflux(mf,volume,dir,scale,scomp,dcomp,nc,geom);
}

void
FluxRegister::Reflux (MultiFab& mf, const MultiFab& volume, Orientation face,
                      Real scale, int scomp, int dcomp, int nc, const Geometry& geom)
{
    BL_PROFILE("FluxRegister::Reflux()");

    int idir = face.coordDir();

    MultiFab flux(amrex::convert(mf.boxArray(), IntVect::TheDimensionVector(idir)),
                  mf.DistributionMap(), nc, 0, MFInfo(), mf.Factory());
    flux.setVal(0.0);

    bndry[face].copyTo(flux, 0, scomp, 0, nc, geom.periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& sfab = mf.array(mfi);
        Array4<Real const> const& ffab = flux.const_array(mfi);
        Array4<Real const> const& vfab = volume.const_array(mfi);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA (bx, tbx,
        {
            fluxreg_reflux(tbx, sfab, dcomp, ffab, vfab, nc, scale, face);
        });
    }
}

void
FluxRegister::ClearInternalBorders (const Geometry& geom)
{
    Gpu::LaunchSafeGuard lsg(false); // xxxxx gpu todo

    int nc = this->nComp();
    const Box& domain = geom.Domain();
    
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
	Orientation lo(dir, Orientation::low);
	Orientation hi(dir, Orientation::high);
	
	FabSet& frlo = (*this)[lo];
	FabSet& frhi = (*this)[hi];
	
	const BoxArray& balo = frlo.boxArray();
	const BoxArray& bahi = frhi.boxArray();
	
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
	    for (FabSetIter fsi(frlo); fsi.isValid(); ++fsi) {
		const Box& bx = fsi.validbox();
		const std::vector< std::pair<int,Box> >& isects = bahi.intersections(bx);
		for (int ii = 0; ii < static_cast<int>(isects.size()); ++ii) {
		    frlo[fsi].setVal<RunOn::Host>(0.0, isects[ii].second, 0, nc);
		}
		if (geom.isPeriodic(dir)) {
		    if (bx.smallEnd(dir) == domain.smallEnd(dir)) {
			const Box& sbx = amrex::shift(bx, dir, domain.length(dir));
			const std::vector<std::pair<int,Box> >& isects2 = bahi.intersections(sbx);
			for (int ii = 0; ii < static_cast<int>(isects2.size()); ++ii) {
			    const Box& bx2 = amrex::shift(isects2[ii].second, dir, -domain.length(dir));
			    frlo[fsi].setVal<RunOn::Host>(0.0, bx2, 0, nc);
			}		      
		    }
                }
	    }
	    
	    for (FabSetIter fsi(frhi); fsi.isValid(); ++fsi) {
		const Box& bx = fsi.validbox();
		const std::vector< std::pair<int,Box> >& isects = balo.intersections(bx);
		for (int ii = 0; ii < static_cast<int>(isects.size()); ++ii) {
		    frhi[fsi].setVal<RunOn::Host>(0.0, isects[ii].second, 0, nc);
		}
		if (geom.isPeriodic(dir)) {
		    if (bx.bigEnd(dir) == domain.bigEnd(dir)) {
			const Box& sbx = amrex::shift(bx, dir, -domain.length(dir));
			const std::vector<std::pair<int,Box> >& isects2 = balo.intersections(sbx);
			for (int ii = 0; ii < static_cast<int>(isects2.size()); ++ii) {
			    const Box& bx2 = amrex::shift(isects2[ii].second, dir, domain.length(dir));
			    frhi[fsi].setVal<RunOn::Host>(0.0, bx2, 0, nc);
			}		      
		    }
		}
	    }
	}
    }
}

#ifndef BL_NO_FORT
void
FluxRegister::OverwriteFlux (Array<MultiFab*,AMREX_SPACEDIM> const& crse_fluxes,
                             Real scale, int srccomp, int destcomp, int numcomp,
                             const Geometry& crse_geom)
{
    Gpu::LaunchSafeGuard lsg(false); // xxxxx gpu todo

    BL_PROFILE("FluxRegister::OverwriteFlux()");

    const auto& cperiod = crse_geom.periodicity();

    Box cdomain = crse_geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (crse_geom.isPeriodic(idim)) cdomain.grow(idim, 1);
    }

    // cell-centered mask: 0: coarse, 1: covered by fine, 2: phys bc
    const BoxArray& cba = amrex::convert(crse_fluxes[0]->boxArray(), IntVect::TheCellVector());
    iMultiFab cc_mask(cba, crse_fluxes[0]->DistributionMap(), 1, 1);
    {
        const std::vector<IntVect>& pshifts = cperiod.shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::vector< std::pair<int,Box> > isects;

            for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
            {
                IArrayBox& fab = cc_mask[mfi];
                fab.setVal<RunOn::Host>(0);  // coarse by default
                fab.setComplement<RunOn::Host>(2, cdomain, 0, 1); // phys bc

                const Box& bx = fab.box();
                for (const auto& iv : pshifts)
                {
                    grids.intersections(bx+iv, isects);
                    for (const auto& is : isects)
                    {
                        fab.setVal<RunOn::Host>(1, is.second-iv, 0, 1);
                    }
                }
            }
        }
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        MultiFab& crse_flux = *crse_fluxes[idim];
        MultiFab fine_flux(crse_flux.boxArray(), crse_flux.DistributionMap(),
                           numcomp, 0);
        fine_flux.setVal(0.0);

        Orientation lo_face(idim, Orientation::low);
        Orientation hi_face(idim, Orientation::high);
        
        bndry[lo_face].copyTo(fine_flux, 0, srccomp, 0, numcomp, cperiod);
        bndry[hi_face].plusTo(fine_flux, 0, srccomp, 0, numcomp, cperiod);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse_flux,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            amrex_froverwrite_cfb(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_N_ANYD(crse_flux[mfi],destcomp),
                                  BL_TO_FORTRAN_ANYD(fine_flux[mfi]),
                                  BL_TO_FORTRAN_ANYD(cc_mask[mfi]),
                                  &numcomp, &idim, &scale);
        }
    }
}
#endif

void
FluxRegister::write (const std::string& name, std::ostream& os) const
{
    if (ParallelDescriptor::IOProcessor())
    {
        os << ratio      << '\n';
        os << fine_level << '\n';
        os << ncomp      << '\n';
    }

    const BndryRegister* br = this;

    br->write(name,os);
}


void
FluxRegister::read (const std::string& name, std::istream& is)
{
    if (ncomp < 0) {
	amrex::Abort("FluxRegister::read: FluxRegister not defined");
    }

    IntVect ratio_in;
    int fine_level_in;
    int ncomp_in;

    is >> ratio_in;
    is >> fine_level_in;
    is >> ncomp_in;

    if (ratio_in != ratio || fine_level_in != fine_level || ncomp_in != ncomp) {
	amrex::Abort("FluxRegister::read: predefined FluxRegister does not match the one in istream");
    }

    BndryRegister* br = this;

    br->read(name,is);
}

}
