
#include <AMReX_BArena.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_FLUXREG_F.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ccse-mpi.H>

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
FluxRegister::refRatio () const
{
    return ratio;
}

int
FluxRegister::fineLevel () const
{
    return fine_level;
}

int
FluxRegister::crseLevel () const
{
    return fine_level-1;
}

int
FluxRegister::nComp () const
{
    return ncomp;
}

const BoxArray&
FluxRegister::coarsenedBoxes () const
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

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
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

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        const FabSet& lofabs = bndry[Orientation(dir,Orientation::low) ];
        const FabSet& hifabs = bndry[Orientation(dir,Orientation::high)];

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
        for (FabSetIter fsi(lofabs); fsi.isValid(); ++fsi)
        {
            sum += (lofabs[fsi].sum(comp) - hifabs[fsi].sum(comp));
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
#pragma omp parallel
#endif    
    for (MFIter mfi(mflx,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();
	
        const FArrayBox& mflxFAB = mflx[mfi];
        mf[mfi].copy(mflxFAB,bx,srccomp,bx,0,numcomp);

        mf[mfi].mult(mult,bx,0,numcomp);

        for (int i = 0; i < numcomp; i++)
            mf[mfi].mult(area[mfi],bx,bx,0,i,1);
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
#pragma omp parallel
#endif
            for (FabSetIter mfi(fs); mfi.isValid(); ++mfi)
                bndry[face][mfi].plus(fs[mfi],0,destcomp,numcomp);
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

    MultiFab area(mflx.boxArray(), mflx.DistributionMap(), 1, mflx.nGrow(),
                  MFInfo(), mflx.Factory());

    area.setVal(1, 0, 1, area.nGrow());

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
#pragma omp parallel
#endif    
    for (MFIter mfi(mflx,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();
	
        mf[mfi].copy(mflx[mfi],bx,srccomp,bx,0,numcomp);

        mf[mfi].mult(mult,bx,0,numcomp);

        for (int i = 0; i < numcomp; i++)
            mf[mfi].mult(area[mfi],bx,bx,0,i,1);
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

    MultiFab area(mflx.boxArray(), mflx.DistributionMap(), 1, mflx.nGrow(),
                  MFInfo(), mflx.Factory());

    area.setVal(1, 0, 1, area.nGrow());

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
#pragma omp parallel
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],dir,k,srccomp,destcomp,numcomp,mult);
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
#pragma omp parallel
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],area[k],dir,k,srccomp,destcomp,numcomp,mult);
    }
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Box&  flxbox = flux.box();
    const int*  flo    = flxbox.loVect();
    const int*  fhi    = flxbox.hiVect();
    const Real* flxdat = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

#ifndef NDEBUG
    Box cbox = amrex::coarsen(flux.box(),ratio);
    BL_ASSERT(cbox.contains(loreg.box()));
#endif
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFINEADD(lodat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

#ifndef NDEBUG
    BL_ASSERT(cbox.contains(hireg.box()));
#endif
    rlo = hireg.box().loVect();
    rhi = hireg.box().hiVect();
    Real* hidat = hireg.dataPtr(destcomp);
    FORT_FRFINEADD(hidat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       const FArrayBox& area,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Real* area_dat = area.dataPtr();
    const int*  alo      = area.loVect();
    const int*  ahi      = area.hiVect();
    const Box&  flxbox   = flux.box();
    const int*  flo      = flxbox.loVect();
    const int*  fhi      = flxbox.hiVect();
    const Real* flxdat   = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

#ifndef NDEBUG
    Box cbox = amrex::coarsen(flux.box(),ratio);
    BL_ASSERT(cbox.contains(loreg.box()));
#endif
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFAADD(lodat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

#ifndef NDEBUG
    BL_ASSERT(cbox.contains(hireg.box()));
#endif
    rlo = hireg.box().loVect();
    rhi = hireg.box().hiVect();
    Real* hidat = hireg.dataPtr(destcomp);
    FORT_FRFAADD(hidat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);
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
    BL_PROFILE("FluxRegister::Reflux()");

    for (OrientationIter fi; fi; ++fi)
    {
	const Orientation& face = fi();
	int idir = face.coordDir();
	int islo = face.isLow();

        MultiFab flux(amrex::convert(mf.boxArray(), IntVect::TheDimensionVector(idir)),
                      mf.DistributionMap(), nc, 0, MFInfo(), mf.Factory());
	flux.setVal(0.0);

	bndry[face].copyTo(flux, 0, scomp, 0, nc, geom.periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    FArrayBox& sfab = mf[mfi];
	    const Box& sbox = sfab.box();

	    const FArrayBox& ffab = flux[mfi];
	    const Box& fbox = ffab.box();
	    
	    const FArrayBox& vfab = volume[mfi];
	    const Box& vbox = vfab.box();

	    FORT_FRREFLUX(bx.loVect(), bx.hiVect(),
			  sfab.dataPtr(dcomp), sbox.loVect(), sbox.hiVect(),
			  ffab.dataPtr(     ), fbox.loVect(), fbox.hiVect(),
			  vfab.dataPtr(     ), vfab.loVect(), vbox.hiVect(),
			  &nc, &scale, &idir, &islo);
			  
	}
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
    
    MultiFab volume(mf.boxArray(), mf.DistributionMap(), 1, mf.nGrow(),
                    MFInfo(), mf.Factory());
    
    volume.setVal(AMREX_D_TERM(dx[0],*dx[1],*dx[2]), 0, 1, mf.nGrow());

    Reflux(mf,volume,scale,scomp,dcomp,nc,geom);
}

void
FluxRegister::ClearInternalBorders (const Geometry& geom)
{
    int nc = this->nComp();
    const Box& domain = geom.Domain();
    
    for (int dir = 0; dir < BL_SPACEDIM; dir++) {
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
		    frlo[fsi].setVal(0.0, isects[ii].second, 0, nc);
		}
		if (geom.isPeriodic(dir)) {
		    if (bx.smallEnd(dir) == domain.smallEnd(dir)) {
			const Box& sbx = amrex::shift(bx, dir, domain.length(dir));
			const std::vector<std::pair<int,Box> >& isects2 = bahi.intersections(sbx);
			for (int ii = 0; ii < static_cast<int>(isects2.size()); ++ii) {
			    const Box& bx2 = amrex::shift(isects2[ii].second, dir, -domain.length(dir));
			    frlo[fsi].setVal(0.0, bx2, 0, nc);
			}		      
		    }
		}
	    }
	    
	    for (FabSetIter fsi(frhi); fsi.isValid(); ++fsi) {
		const Box& bx = fsi.validbox();
		const std::vector< std::pair<int,Box> >& isects = balo.intersections(bx);
		for (int ii = 0; ii < static_cast<int>(isects.size()); ++ii) {
		    frhi[fsi].setVal(0.0, isects[ii].second, 0, nc);
		}
		if (geom.isPeriodic(dir)) {
		    if (bx.bigEnd(dir) == domain.bigEnd(dir)) {
			const Box& sbx = amrex::shift(bx, dir, -domain.length(dir));
			const std::vector<std::pair<int,Box> >& isects2 = balo.intersections(sbx);
			for (int ii = 0; ii < static_cast<int>(isects2.size()); ++ii) {
			    const Box& bx2 = amrex::shift(isects2[ii].second, dir, domain.length(dir));
			    frhi[fsi].setVal(0.0, bx2, 0, nc);
			}		      
		    }
		}
	    }
	}
    }
}

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

void
FluxRegister::AddProcsToComp(int ioProcNumSCS, int ioProcNumAll,
                             int scsMyId, MPI_Comm scsComm)
{
  // ---- ints
  ParallelDescriptor::Bcast(&fine_level, 1, ioProcNumSCS, scsComm);
  ParallelDescriptor::Bcast(&ncomp, 1, ioProcNumSCS, scsComm);

  // ---- IntVects
  Vector<int> iv(BL_SPACEDIM, -1);
  if(scsMyId == ioProcNumSCS) {
    for(int i(0); i < BL_SPACEDIM; ++i) { iv[i] = ratio[i]; }
  }
  ParallelDescriptor::Bcast(iv.dataPtr(), iv.size(), ioProcNumSCS, scsComm);
  if(scsMyId != ioProcNumSCS) {
    for(int i(0); i < BL_SPACEDIM; ++i) { ratio[i] = iv[i]; }
  }
  BndryRegister::AddProcsToComp(ioProcNumSCS, ioProcNumAll,
                                scsMyId, scsComm);
}

}
