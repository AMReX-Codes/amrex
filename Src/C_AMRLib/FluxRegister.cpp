
//
// $Id: FluxRegister.cpp,v 1.2 1997-11-20 00:49:43 lijewski Exp $
//

#include <FluxRegister.H>
#include <Geometry.H>

#include <FLUXREG_F.H>
#include <ParallelDescriptor.H>

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();          \
const int* fabhi = (fab).hiVect();          \
const Real* fabdat = (fab).dataPtr();

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();          \
const int* fabhi = (fab).hiVect();          \
Real* fabdat = (fab).dataPtr();

// -------------------------------------------------------------
FluxRegister::FluxRegister() : BndryRegister()
{
    fine_level = ncomp = -1;

    ratio = IntVect::TheUnitVector(); ratio.scale(-1);
}

// -------------------------------------------------------------
FluxRegister::FluxRegister(const BoxArray& fine_boxes, 
			   const IntVect & ref_ratio, int fine_lev, int nvar)
    : BndryRegister()
{
    define(fine_boxes,ref_ratio,fine_lev,nvar);
}

// -------------------------------------------------------------
void
FluxRegister::define(const BoxArray& fine_boxes, 
		     const IntVect & ref_ratio, int fine_lev, int nvar)
{
    assert(fine_boxes.isDisjoint());
    assert( ! grids.ready());

    ratio = ref_ratio;
    fine_level = fine_lev;
    ncomp = nvar;
    grids.define(fine_boxes);
    grids.coarsen(ratio);
    int dir;
    for(dir = 0; dir < BL_SPACEDIM; dir++) {
        Orientation lo_face(dir,Orientation::low);
        Orientation hi_face(dir,Orientation::high);
        IndexType typ(IndexType::TheCellType());
        typ.setType(dir,IndexType::NODE);
        BndryRegister::define(lo_face,typ,0,1,0,nvar);
        BndryRegister::define(hi_face,typ,0,1,0,nvar);
    }
}

// -------------------------------------------------------------
FluxRegister::~FluxRegister() {
}

// -------------------------------------------------------------
Real FluxRegister::SumReg(int comp) const {
    Real sum = 0.0;
    for(int dir = 0; dir < BL_SPACEDIM; dir++) {
        Orientation lo_face(dir,Orientation::low);
        Orientation hi_face(dir,Orientation::high);
        const FabSet &lofabs = bndry[lo_face];
        const FabSet &hifabs = bndry[hi_face];
        for(ConstFabSetIterator fsi(lofabs); fsi.isValid(); ++fsi) {
          ConstDependentFabSetIterator dfsi(fsi, hifabs);
            sum += fsi().sum(comp);
            sum -= dfsi().sum(comp);
        }
    }
    ParallelDescriptor::ReduceRealSum(sum);
    return sum;
}

// -------------------------------------------------------------
void FluxRegister::copyTo(MultiFab& mflx, int dir,
                          int src_comp, int dest_comp, int num_comp)
{
    assert( dir >= 0 && dir < BL_SPACEDIM);

    Orientation lo_face(dir,Orientation::low);
    const FabSet &lofabs = bndry[lo_face];
    lofabs.copyTo(mflx,0,src_comp,dest_comp,num_comp);

    Orientation hi_face(dir,Orientation::high);
    const FabSet &hifabs = bndry[hi_face];
    hifabs.copyTo(mflx,0,src_comp,dest_comp,num_comp);
}

// -------------------------------------------------------------
void
FluxRegister::copyTo(FARRAYBOX& flx, int dir,
                     int src_comp, int dest_comp, int num_comp)
{
    assert( dir >= 0 && dir < BL_SPACEDIM);

    Orientation lo_face(dir,Orientation::low);
    const FabSet &lofabs = bndry[lo_face];
    lofabs.copyTo(flx,src_comp,dest_comp,num_comp);

    Orientation hi_face(dir,Orientation::high);
    const FabSet &hifabs = bndry[hi_face];
    hifabs.copyTo(flx,src_comp,dest_comp,num_comp);
}


/*
// -------------------------------------------------------------
void
FluxRegister::addTo(MultiFab& mflx, int dir, REAL mult,
		    int src_comp, int dest_comp, int num_comp)
{
    assert( dir >= 0 && dir < BL_SPACEDIM);
    for(MultiFabIterator mfi(mflx); mfi.isValid(); ++mfi) {
	addTo(mfi(),dir,mult,src_comp,dest_comp,num_comp);
    }
}

// -------------------------------------------------------------
void
FluxRegister::addTo(FARRAYBOX& flx, int dir, REAL mult,
		     int src_comp, int dest_comp, int num_comp)
{
      // flx = flx + mult*register
    int ngrds = grids.length();

    const BOX& bx = flx.box();
    for(FabSetIterator fsi(bndry[face_lo]); fsi.isValid(); ++fsi) {
      DependentFabSetIterator dfsi(fsi, fndry[face_hi]);

        Orientation face_lo(dir,Orientation::low);
        Orientation face_hi(dir,Orientation::high);
        //FARRAYBOX& lo_fab = bndry[face_lo][k];
        //FARRAYBOX& hi_fab = bndry[face_hi][k];
        FARRAYBOX& lo_fab = fsi();
        FARRAYBOX& hi_fab = dfsi();
        BOX lo_box(lo_fab.box());
        BOX hi_box(hi_fab.box());
        lo_box &= bx;
        hi_box &= bx;
        if (lo_box.ok()) {
            FARRAYBOX tmp(lo_box,num_comp);
            tmp.copy(lo_fab,src_comp,0,num_comp);
            tmp.mult(mult);
            flx.plus(tmp,0,dest_comp,num_comp);
        }
        if (hi_box.ok()) {
            FARRAYBOX tmp(hi_box,num_comp);
            tmp.copy(hi_fab,src_comp,0,num_comp);
            tmp.mult(mult);
            flx.plus(tmp,0,dest_comp,num_comp);
        }
    }
}
*/

/*
// -------------------------------------------------------------
void
FluxRegister::scaleAddTo(MultiFab& mflx, const MultiFab& area,
                         int dir, REAL mult,
                         int src_comp, int dest_comp, int num_comp)
{
    assert( dir >= 0 && dir < BL_SPACEDIM);
    for(MultiFabIterator mfi(mflx); mfi.isValid(); ++mfi) {
        DependentMultiFabIterator dmfi(mfi, area);
	scaleAddTo(mfi(),dmfi(),dir,mult,src_comp,dest_comp,num_comp);
    }
}

// -------------------------------------------------------------
void
FluxRegister::scaleAddTo(FARRAYBOX& flx, const FARRAYBOX& area,
                         int dir, REAL mult,
                         int src_comp, int dest_comp, int num_comp)
{
      // flx = flx + mult*register
    int ngrds = grids.length();


    const BOX& flux_bx = flx.box();
    const int * flo = flux_bx.loVect();
    const int * fhi = flux_bx.hiVect();

    const BOX& areabox = area.box();
    const REAL* area_dat = area.dataPtr();
    const int * alo = areabox.loVect();
    const int * ahi = areabox.hiVect();

    for(FabSetIterator fsi(bndry[face_lo]); fsi.isValid(); ++fsi) {
      DependentFabSetIterator dfsi(fsi, fndry[face_hi]);

        Orientation lo_face(dir,Orientation::low);
        //FARRAYBOX& lo_fab = bndry[lo_face][k];
        FARRAYBOX& lo_fab = fsi();

        Orientation hi_face(dir,Orientation::high);
        //FARRAYBOX& hi_fab = bndry[hi_face][k];
        FARRAYBOX& hi_fab = dfsi();

        BOX rlo_box(lo_fab.box());
        BOX rhi_box(hi_fab.box());

        BOX lo_box(lo_fab.box());
        BOX hi_box(hi_fab.box());

        lo_box &= flux_bx;
        hi_box &= flux_bx;

        if (lo_box.ok()) {

            const int * rlo = rlo_box.loVect();
            const int * rhi = rlo_box.hiVect();

            const int * lo = lo_box.loVect();
            const int * hi = lo_box.hiVect();

            FORT_SCALADDTO(flx.dataPtr(dest_comp),ARLIM(flo),ARLIM(fhi),
                           area_dat,ARLIM(alo),ARLIM(ahi),
                           lo_fab.dataPtr(src_comp),ARLIM(rlo),ARLIM(rhi),
                           lo,hi,&num_comp,&mult);
        }

        if (hi_box.ok()) {

            const int * rlo = rhi_box.loVect();
            const int * rhi = rhi_box.hiVect();

            const int * lo = hi_box.loVect();
            const int * hi = hi_box.hiVect();

            FORT_SCALADDTO(flx.dataPtr(dest_comp),ARLIM(flo),ARLIM(fhi),
                           area_dat,ARLIM(alo),ARLIM(ahi),
                           hi_fab.dataPtr(src_comp),ARLIM(rlo),ARLIM(rhi),
                           lo,hi,&num_comp,&mult);


        }
    }
}
*/





// -------------------------------------------------------------
void FluxRegister::Reflux(MultiFab &S, const MultiFab &volume, Real scale,
                          int src_comp, int dest_comp, int num_comp, 
                          const Geometry&geom)
{
    BoxLib::Error("FluxRegister::Reflux(MultiFab&, const MultiFab&, ...) not implemented");
/*
    int nreg = grids.length();
    const BoxArray& grd_boxes = S.boxArray();
    int ngrd = grd_boxes.length();

    int grd;
    for (grd = 0; grd < ngrd; grd++) {
        FARRAYBOX& s = S[grd];
        const BOX& s_box = grd_boxes[grd];
        REAL* s_dat = s.dataPtr(dest_comp);
        const int* slo = s.loVect();
        const int* shi = s.hiVect();
        const FARRAYBOX& vol = volume[grd];
        const REAL* vol_dat = vol.dataPtr();
        const int* vlo = vol.loVect();
        const int* vhi = vol.hiVect();

          // find flux register that intersect with this grid
        int k;
        for (k = 0; k < nreg; k++) {
            const BOX& reg_box = grids[k];
            BOX bx(grow(reg_box,1));
            if (bx.intersects(s_box)) {
                for (OrientationIter fi; fi; ++fi) {
                    Orientation face = fi();
                    BOX fine_face(adjCell(reg_box,face));
                      // low(hight)  face of fine grid => high (low)
                      // face of the exterior crarse grid cell updated.
                      // adjust sign of scale accordingly.
                    REAL mult = (face.isLow() ? -scale : scale);
                    BOX ovlp(s_box);
                    ovlp &= fine_face;
                    if (ovlp.ok()) {
                        const FARRAYBOX& reg = bndry[face][k];
                        const REAL* reg_dat = reg.dataPtr(src_comp);
                        const int* rlo = fine_face.loVect();
                        const int* rhi = fine_face.hiVect();
                        const int* lo = ovlp.loVect();
                        const int* hi = ovlp.hiVect();
                        FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                                      vol_dat,ARLIM(vlo),ARLIM(vhi),
                                      reg_dat,ARLIM(rlo),ARLIM(rhi),
                                      lo,hi,&num_comp,&mult);
                    }
                }
            }
#if 1
	    // add periodic possibilities
	    if( geom.isAnyPeriodic() && !geom.Domain().contains(bx)){
	      Array<IntVect> pshifts(27);
	      geom.periodicShift(bx,s_box,pshifts);
	      for(int iiv=0; iiv<pshifts.length(); iiv++){
		IntVect iv = pshifts[iiv];
		s.shift(iv);
		const int* slo = s.loVect();
		const int* shi = s.hiVect();
		// this is a funny situation.  I don't want to permanently
		// change vol, but I need to do a shift on it.  I'll shift
		// it back later, so the overall change is nil.  But to do
		// this, I have to cheat and do a cast.  This is pretty 
		// disgusting.
		FArrayBox& cheatvol = *(FArrayBox *)&vol;
		cheatvol.shift(iv);
		const int* vlo = cheatvol.loVect();
		const int* vhi = cheatvol.hiVect();
		Box s_box = grd_boxes[grd];
		D_TERM( s_box.shift(0,iv[0]);,
			s_box.shift(1,iv[1]);,
			s_box.shift(2,iv[2]); )
		if( !bx.intersects(s_box) ){
		  cerr << "FluxRegister::Reflux logic error"<<'\n';
		  exit(1);
		}

		for (OrientationIter fi; fi; ++fi) {
		  Orientation face = fi();
		  Box fine_face(adjCell(reg_box,face));
		  // low(hight)  face of fine grid => high (low)
		  // face of the exterior crarse grid cell updated.
		  // adjust sign of scale accordingly.
		  Real mult = (face.isLow() ? -scale : scale);
		  Box ovlp(s_box);
		  ovlp &= fine_face;
		  if (ovlp.ok()) {
		    const FArrayBox& reg = bndry[face][k];
		    const Real* reg_dat = reg.dataPtr(src_comp);
		    const int* rlo = fine_face.loVect();
		    const int* rhi = fine_face.hiVect();
		    const int* lo = ovlp.loVect();
		    const int* hi = ovlp.hiVect();
		    FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
				  vol_dat,ARLIM(vlo),ARLIM(vhi),
				  reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
				  &num_comp,&mult);
		  }
                }

		s.shift(-iv);
		cheatvol.shift(-iv);
	      }
	    }
#endif
        }
    }
*/
}

// -------------------------------------------------------------
void FluxRegister::Reflux(MultiFab &S, Real scale,
                     int src_comp, int dest_comp, int num_comp, 
		     const Geometry &geom)
{
/*
    int nreg = grids.length();
    const BoxArray& grd_boxes = S.boxArray();
    int ngrd = grd_boxes.length();

    const REAL* dx = geom.CellSize();

    for (int grd = 0; grd < ngrd; grd++) {
        FARRAYBOX& s = S[grd];
        const BOX& s_box = grd_boxes[grd];
        REAL* s_dat = s.dataPtr(dest_comp);
        const int* slo = s.loVect();
        const int* shi = s.hiVect();

          // find flux register that intersect with this grid
        int k;
        for (k = 0; k < nreg; k++) {
            const BOX& reg_box = grids[k];
            BOX bx(grow(reg_box,1));
            if (bx.intersects(s_box)) {
                for (OrientationIter fi; fi; ++fi) {
                    Orientation face = fi();
                    BOX fine_face(adjCell(reg_box,face));
                      // low(hight)  face of fine grid => high (low)
                      // face of the exterior crarse grid cell updated.
                      // adjust sign of scale accordingly.
                    REAL mult = (face.isLow() ? -scale : scale);
                    BOX ovlp(s_box);
                    ovlp &= fine_face;
                    if (ovlp.ok()) {
                        const FARRAYBOX& reg = bndry[face][k];
                        const REAL* reg_dat = reg.dataPtr(src_comp);
                        const int* rlo = fine_face.loVect();
                        const int* rhi = fine_face.hiVect();
                        const int* lo = ovlp.loVect();
                        const int* hi = ovlp.hiVect();
                        FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                                        reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                                        &num_comp,&mult);
                    }
                }
            }
#if 1
	    // add periodic possibilities
	    if( geom.isAnyPeriodic() && !geom.Domain().contains(bx)){
	      Array<IntVect> pshifts(27);
	      geom.periodicShift(bx,s_box,pshifts);
	      for(int iiv=0; iiv<pshifts.length(); iiv++){
		IntVect iv = pshifts[iiv];
		s.shift(iv);
		const int* slo = s.loVect();
		const int* shi = s.hiVect();
		Box s_box = grd_boxes[grd];
		D_TERM( s_box.shift(0,iv[0]);,
			s_box.shift(1,iv[1]);,
			s_box.shift(2,iv[2]); )
		if( !bx.intersects(s_box) ){
		  cerr << "FluxRegister::Reflux logic error"<<'\n';
		  exit(1);
		}

		for (OrientationIter fi; fi; ++fi) {
		  Orientation face = fi();
		  BOX fine_face(adjCell(reg_box,face));
		  // low(hight)  face of fine grid => high (low)
		  // face of the exterior crarse grid cell updated.
		  // adjust sign of scale accordingly.
		  REAL mult = (face.isLow() ? -scale : scale);
		  BOX ovlp(s_box);
		  ovlp &= fine_face;
		  if (ovlp.ok()) {
		    const FARRAYBOX& reg = bndry[face][k];
		    const REAL* reg_dat = reg.dataPtr(src_comp);
		    const int* rlo = fine_face.loVect();
		    const int* rhi = fine_face.hiVect();
		    const int* lo = ovlp.loVect();
		    const int* hi = ovlp.hiVect();
		    FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
				    reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
				    &num_comp,&mult);
		  }
                }

		s.shift(-iv);
	      }
	    }
#endif
        }
    }
*/
    const REAL* dx = geom.CellSize();

    FabSetCopyDescriptor fscd(true);
    FabSetId fsid[2*BL_SPACEDIM];
    for(OrientationIter fi; fi; ++fi) {
      fsid[fi()] = fscd.RegisterFabSet(&(bndry[fi()]));
    }
    List<FillBoxId> fillBoxIdList;
    FillBoxId tempFillBoxId;

    for(MultiFabIterator mfi(S); mfi.isValid(); ++mfi) {
        const Box &s_box = mfi.validbox();
          // find flux register that intersect with this grid
        for(int k = 0; k < grids.length(); k++) {
            const Box &reg_box = grids[k];
            Box bx(grow(reg_box,1));
            if(bx.intersects(s_box)) {
                for(OrientationIter fi; fi; ++fi) {
                    Orientation face = fi();
                    Box fine_face(adjCell(reg_box,face));
                    Box ovlp(s_box);
                    ovlp &= fine_face;
                    if(ovlp.ok()) {
                        //Box regBox(bndry[face][k].box());
                        Box regBox(bndry[face].box(k));
                        BoxList unfilledBoxes(regBox.ixType());

                        tempFillBoxId = fscd.AddBox(fsid[face],regBox,unfilledBoxes,
                                                    src_comp, dest_comp, num_comp);
                        //assert(unfilledBoxes.isEmpty());
                        fillBoxIdList.append(tempFillBoxId);
                    }
                }
            }  // end if(bx.intersects(s_box))

#if 1
            // add periodic possibilities
            if( geom.isAnyPeriodic() && !geom.Domain().contains(bx)){
              FArrayBox &s = mfi();
              Array<IntVect> pshifts(27);
              geom.periodicShift(bx,s_box,pshifts);
              for(int iiv=0; iiv<pshifts.length(); iiv++){
                IntVect iv = pshifts[iiv];
                s.shift(iv);
                const int *slo = s.loVect();
                const int *shi = s.hiVect();
                Box s_box = mfi.validbox();
                D_TERM( s_box.shift(0,iv[0]);,
                        s_box.shift(1,iv[1]);,
                        s_box.shift(2,iv[2]); )
                assert(bx.intersects(s_box));

                for(OrientationIter fi; fi; ++fi) {
                  Orientation face = fi();
                  BOX fine_face(adjCell(reg_box,face));
                  // low(hight)  face of fine grid => high (low)
                  // face of the exterior coarse grid cell updated.
                  // adjust sign of scale accordingly.
                  REAL mult = (face.isLow() ? -scale : scale);
                  BOX ovlp(s_box);
                  ovlp &= fine_face;
                  if(ovlp.ok()) {
                    Box regBox(bndry[face].box(k));
                    BoxList unfilledBoxes(regBox.ixType());

                    tempFillBoxId = fscd.AddBox(fsid[face],regBox,unfilledBoxes,
                                                src_comp, dest_comp, num_comp);
                    //assert(unfilledBoxes.isEmpty());
                    fillBoxIdList.append(tempFillBoxId);

                    //const FARRAYBOX &reg = bndry[face][k];
                    //const REAL *reg_dat = reg.dataPtr(src_comp);
                    //const int *rlo = fine_face.loVect();
                    //const int *rhi = fine_face.hiVect();
                    //const int *lo = ovlp.loVect();
                    //const int *hi = ovlp.hiVect();
                    //FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                                    //reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                                    //&num_comp,&mult);
                  }
                }  // end for(fi...)
                s.shift(-iv);
              }
            }  // end if(periodic...)
#endif

        }  // end for(k...)
    }  // end for(MultiFabIterator...)

    Array<FillBoxId> fillBoxId(fillBoxIdList.length());
    int ifbi = 0;
    for(ListIterator<FillBoxId> li(fillBoxIdList); li; ++li) {
     fillBoxId[ifbi] = li();
     ++ifbi;
    }
    fillBoxIdList.clear();

    fscd.CollectData();

    int overlapId = 0;
    for(MultiFabIterator mfi(S); mfi.isValid(); ++mfi) {
        const Box &s_box = mfi.validbox();
          // find flux register that intersect with this grid
        for(int k = 0; k < grids.length(); k++) {
            const Box &reg_box = grids[k];
            Box bx(grow(reg_box,1));
            if(bx.intersects(s_box)) {
                for(OrientationIter fi; fi; ++fi) {
                    Orientation face = fi();
                    Box fine_face(adjCell(reg_box,face));
                      // low(hight)  face of fine grid => high (low)
                      // face of the exterior coarse grid cell updated.
                      // adjust sign of scale accordingly.
                    Real mult = (face.isLow() ? -scale : scale);
                    Box ovlp(s_box);
                    ovlp &= fine_face;
                    if(ovlp.ok()) {
                      FArrayBox &sfab = mfi();
                      Real *s_dat = sfab.dataPtr(dest_comp);
                      const int *slo = sfab.loVect();
                      const int *shi = sfab.hiVect();
                      FillBoxId fbid = fillBoxId[overlapId];
                      FArrayBox reg(fbid.box(), num_comp);
                      fscd.FillFab(fsid[face], fbid, reg);
                      //const FArrayBox &reg = bndry[face][k];
                      //const Real *reg_dat = reg.dataPtr(src_comp);
                       const Real *reg_dat = reg.dataPtr(0);
                      const int *rlo = fine_face.loVect();
                      const int *rhi = fine_face.hiVect();
                      const int *lo  = ovlp.loVect();
                      const int *hi  = ovlp.hiVect();
                      FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                                      reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                                      &num_comp,&mult);
                      ++overlapId;
                    }
                }  // end for(fi...)
            }  // end if(bx.intersects(s_box))

#if 1
            // add periodic possibilities
            if( geom.isAnyPeriodic() && !geom.Domain().contains(bx)){
              FArrayBox &s = mfi();
              Array<IntVect> pshifts(27);
              geom.periodicShift(bx,s_box,pshifts);
              for(int iiv=0; iiv<pshifts.length(); iiv++){
                IntVect iv = pshifts[iiv];
                s.shift(iv);
                const int *slo = s.loVect();
                const int *shi = s.hiVect();
                Box s_box = mfi.validbox();
                D_TERM( s_box.shift(0,iv[0]);,
                        s_box.shift(1,iv[1]);,
                        s_box.shift(2,iv[2]); )
                assert(bx.intersects(s_box));

                for(OrientationIter fi; fi; ++fi) {
                  Orientation face = fi();
                  BOX fine_face(adjCell(reg_box,face));
                  REAL mult = (face.isLow() ? -scale : scale);
                  BOX ovlp(s_box);
                  ovlp &= fine_face;
                  if(ovlp.ok()) {
                    FArrayBox &sfab = mfi();
                    Real *s_dat = sfab.dataPtr(dest_comp);
                    const int *slo = sfab.loVect();
                    const int *shi = sfab.hiVect();
                    FillBoxId fbid = fillBoxId[overlapId];
                    Box regBox(bndry[face].box(k));
                    assert(regBox == fbid.box());
                    FArrayBox reg(fbid.box(), num_comp);
                    fscd.FillFab(fsid[face], fbid, reg);

                    //const FARRAYBOX &reg = bndry[face][k];
                    //const REAL *reg_dat = reg.dataPtr(src_comp);
                    const REAL *reg_dat = reg.dataPtr(0);
                    const int *rlo = fine_face.loVect();
                    const int *rhi = fine_face.hiVect();
                    const int *lo = ovlp.loVect();
                    const int *hi = ovlp.hiVect();
                    FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                                    reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                                    &num_comp,&mult);
                    ++overlapId;
                  }
                }  // end for(fi...)
                s.shift(-iv);
              }
            }  // end if(periodic...)
#endif

        }  // end for(k...)
    }  // end for(MultiFabIterator...)
}  // end FluxRegister::Reflux(...)


// -------------------------------------------------------------
void
FluxRegister::CrseInit(const MultiFab& mflx, int dir,
		       int srccomp, int destcomp, int numcomp, REAL mult)
{
    BoxLib::Error("CrseInit(multifab, ...) not implemented");
/*
    const BoxArray& bxa = mflx.boxArray();
    for(ConstMultiFabIterator mfi(mflx); mfi.isValid(); ++mfi) {
        assert(mfi.box() == bxa[mfi.index()]);
	CrseInit(mfi(),mfi.box(),dir,srccomp,destcomp,numcomp,mult);
    }
*/
}

// -------------------------------------------------------------
void
FluxRegister::CrseInit(const MultiFab& mflx, const MultiFab& area,
		       int dir, int srccomp, int destcomp,
		       int numcomp, REAL mult)
{
    BoxLib::Error("CrseInit(multifab, multifab, ...) not implemented");
/*
    const BoxArray& bxa = mflx.boxArray();
    for(ConstMultiFabIterator mfi(mflx); mfi.isValid(); ++mfi) {
        ConstDependentMultiFabIterator dmfi(area); dmfi.isValid(); ++dmfi);
        assert(mfi.box() == bxa[mfi.index()]);
	CrseInit(mfi(),dmfi(),mfi.box(),dir,srccomp,destcomp,numcomp,mult);
    }
*/
}




/*
// working nonparallel code vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
// -------------------------------------------------------------
void
FluxRegister::CrseInit(const FARRAYBOX& flux, const BOX& subbox, int dir,
		       int srccomp, int destcomp, int numcomp, REAL mult)
{
    int nvf = flux.nComp();
    assert(srccomp >= 0 && srccomp+numcomp <= nvf);
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const BOX& flxbox = flux.box();
    assert(flxbox.contains(subbox));
    const int* flo = flxbox.loVect();
    const int* fhi = flxbox.hiVect();
    const REAL* flxdat = flux.dataPtr(srccomp);

    int nreg = grids.length();
    int k;
    for (k = 0; k < nreg; k++) {
	Orientation face_lo(dir,Orientation::low);
	FARRAYBOX& loreg = bndry[face_lo][k];
	BOX lobox(loreg.box());
	lobox &= subbox;
	if (lobox.ok()) {
	    const int* rlo = loreg.loVect();
	    const int* rhi = loreg.hiVect();
	    REAL* lodat = loreg.dataPtr(destcomp);
	    const int* lo = lobox.loVect();
	    const int* hi = lobox.hiVect();
	    FORT_FRCRSEINIT(lodat,ARLIM(rlo),ARLIM(rhi),
                            flxdat,ARLIM(flo),ARLIM(fhi),
                            lo,hi,&numcomp,&dir,&mult);
	  // should be able to replace FORT_FRCRSEINIT with
	  // loreg.copy(flux, lobox, srccomp, lobox, destcomp, numcomp);
	  // loreg.mult(mult, lobox, destcomp, numcomp);
	}
	Orientation face_hi(dir,Orientation::high);
	FARRAYBOX& hireg = bndry[face_hi][k];
	BOX hibox(hireg.box());
	hibox &= subbox;
	if (hibox.ok()) {
	    const int* rlo = hireg.loVect();
	    const int* rhi = hireg.hiVect();
	    REAL* hidat = hireg.dataPtr(destcomp);
	    const int* lo = hibox.loVect();
	    const int* hi = hibox.hiVect();
	    FORT_FRCRSEINIT(hidat,ARLIM(rlo),ARLIM(rhi),
                            flxdat,ARLIM(flo),ARLIM(fhi),
                            lo,hi,&numcomp,&dir,&mult);
	  // should be able to replace FORT_FRCRSEINIT with
	  // hireg.copy(flux, hibox, srccomp, hibox, destcomp, numcomp);
          // hireg.mult(mult, hibox, destcomp, numcomp);
	}
    }
}


void FluxRegister::CrseInitFinish() {
}
*/



// -------------------------------------------------------------
void FluxRegister::CrseInit(const FArrayBox &flux, const Box &subbox, int dir,
                            int srccomp, int destcomp, int numcomp, Real mult)
{
    assert(srccomp  >= 0 && srccomp+numcomp  <= flux.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);

    int myproc = ParallelDescriptor::MyProc();
    FabComTag fabComTag;

    const Box &flxbox = flux.box();
    assert(flxbox.contains(subbox));
    const int  *flo    = flxbox.loVect();
    const int  *fhi    = flxbox.hiVect();
    const Real *flxdat = flux.dataPtr(srccomp);

    for(int k = 0; k < grids.length(); k++) {
        Orientation face_lo(dir,Orientation::low);
        Box lobox(bndry[face_lo].box(k));
        lobox &= subbox;
        if(lobox.ok()) {
          const DistributionMapping &distributionMap = bndry[face_lo].DistributionMap();
          if(myproc == distributionMap[k]) {  // local
            FArrayBox &loreg = bndry[face_lo][k];
            // replaced FORT_FRCRSEINIT(...)
            loreg.copy(flux, lobox, srccomp, lobox, destcomp, numcomp);
            loreg.mult(mult, lobox, destcomp, numcomp);
          } else {
            FArrayBox fabCom(lobox, numcomp);
            int fabComDestComp = 0;
            fabCom.copy(flux, lobox, srccomp, lobox, fabComDestComp, numcomp);
            fabCom.mult(mult, lobox, fabComDestComp, numcomp);
            fabComTag.fromProc = myproc;
            fabComTag.toProc   = distributionMap[k];
            fabComTag.fabIndex = k;
            fabComTag.destComp = destcomp;
            fabComTag.nComp    = fabCom.nComp();
            fabComTag.box      = fabCom.box();
            fabComTag.face     = face_lo;
            ParallelDescriptor::SendData(fabComTag.toProc, &fabComTag,
                                         fabCom.dataPtr(),
                                         fabComTag.box.numPts() *
                                         fabComTag.nComp * sizeof(Real));
          }
        }
        Orientation face_hi(dir,Orientation::high);
        Box hibox(bndry[face_hi].box(k));
        hibox &= subbox;
        if(hibox.ok()) {
          const DistributionMapping &distributionMap = bndry[face_hi].DistributionMap();
          if(myproc == distributionMap[k]) {  // local
            FArrayBox &hireg = bndry[face_hi][k];
            // replaced FORT_FRCRSEINIT(...)
            hireg.copy(flux, hibox, srccomp, hibox, destcomp, numcomp);
            hireg.mult(mult, hibox, destcomp, numcomp);
          } else {
            FArrayBox fabCom(hibox, numcomp);
            int fabComDestComp = 0;
            fabCom.copy(flux, hibox, srccomp, hibox, fabComDestComp, numcomp);
            fabCom.mult(mult, hibox, fabComDestComp, numcomp);
            fabComTag.fromProc = myproc;
            fabComTag.toProc   = distributionMap[k];
            fabComTag.fabIndex = k;
            fabComTag.destComp = destcomp;
            fabComTag.nComp    = fabCom.nComp();
            fabComTag.box      = fabCom.box();
            fabComTag.face     = face_hi;
            ParallelDescriptor::SendData(fabComTag.toProc, &fabComTag,
                                         fabCom.dataPtr(),
                                         fabComTag.box.numPts() *
                                         fabComTag.nComp * sizeof(Real));
          }
        }
    }
}



// -------------------------------------------------------------
void FluxRegister::CrseInitFinish() {

    FabComTag fabComTag;
    ParallelDescriptor::SetMessageHeaderSize(sizeof(FabComTag));

    int dataWaitingSize;

    while (ParallelDescriptor::GetMessageHeader(dataWaitingSize, &fabComTag))
    {
        //
        // Data was sent to this processor.
        //
        long t_long = fabComTag.box.numPts() * fabComTag.nComp * sizeof(Real);
        assert(t_long < INT_MAX);
        int shouldReceiveBytes = int(t_long);

      if (dataWaitingSize != shouldReceiveBytes)
      {
        cerr << "Error in FluxRegister::CrseInitFinish():  "
             << "dataWaitingSize != shouldReceiveBytes:  = "
             << dataWaitingSize << " != " << shouldReceiveBytes << '\n';
        BoxLib::Error("Bad received nbytes");
      }
      if (!fabComTag.box.ok())
      {
          BoxLib::Error("FluxRegister::CrseInitFinish(): bad fabComTag.box");
      }


      FArrayBox tempFab(fabComTag.box, fabComTag.nComp);
      ParallelDescriptor::ReceiveData(tempFab.dataPtr(),
               fabComTag.box.numPts() * fabComTag.nComp * sizeof(Real));
      int srcComp = 0;
      bndry[fabComTag.face][fabComTag.fabIndex].copy(tempFab, fabComTag.box,
                                       srcComp, fabComTag.box,
                                       fabComTag.destComp, fabComTag.nComp);
    }
}



// -------------------------------------------------------------
void
FluxRegister::CrseInit(const FARRAYBOX& flux, const FARRAYBOX& area,
		       const BOX& subbox, int dir,
		       int srccomp, int destcomp, int numcomp, REAL mult)
{
    BoxLib::Error("CrseInit(fab, fab, ...) not implemented");
/*
    int nvf = flux.nComp();
    assert(srccomp >= 0 && srccomp+numcomp <= nvf);
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const BOX& flxbox = flux.box();
    assert(flxbox.contains(subbox));
    const int* flo = flxbox.loVect();
    const int* fhi = flxbox.hiVect();
    const REAL* flx_dat = flux.dataPtr(srccomp);

    const BOX& areabox = area.box();
    assert(areabox.contains(subbox));
    const int* alo = areabox.loVect();
    const int* ahi = areabox.hiVect();
    const REAL* area_dat = area.dataPtr();

    int nreg = grids.length();
    int k;
    for (k = 0; k < nreg; k++) {
	Orientation face_lo(dir,Orientation::low);
	FARRAYBOX& loreg = bndry[face_lo][k];
	BOX lobox(loreg.box());
	lobox &= subbox;
	if (lobox.ok()) {
	    const int* rlo = loreg.loVect();
	    const int* rhi = loreg.hiVect();
	    REAL* lodat = loreg.dataPtr(destcomp);
	    const int* lo = lobox.loVect();
	    const int* hi = lobox.hiVect();
	    FORT_FRCAINIT(lodat,ARLIM(rlo),ARLIM(rhi),
                          flx_dat,ARLIM(flo),ARLIM(fhi),
			  area_dat,ARLIM(alo),ARLIM(ahi),         
		          lo,hi,&numcomp,&dir,&mult);
	}
	Orientation face_hi(dir,Orientation::high);
	FARRAYBOX& hireg = bndry[face_hi][k];
	BOX hibox(hireg.box());
	hibox &= subbox;
	if (hibox.ok()) {
	    const int* rlo = hireg.loVect();
	    const int* rhi = hireg.hiVect();
	    REAL* hidat = hireg.dataPtr(destcomp);
	    const int* lo = hibox.loVect();
	    const int* hi = hibox.hiVect();
	    FORT_FRCAINIT(hidat,ARLIM(rlo),ARLIM(rhi), 
                          flx_dat,ARLIM(flo),ARLIM(fhi),
			  area_dat,ARLIM(alo),ARLIM(ahi),lo,hi,&numcomp,
			  &dir,&mult);
	}
    }
*/
}

// -------------------------------------------------------------
void
FluxRegister::FineAdd(const MultiFab& mflx, int dir,
		      int srccomp, int destcomp, int numcomp, REAL mult)
{
    BoxLib::Error("FineAdd(multifab, ...) not implemented");
/*
    const BoxArray& bxa = mflx.boxArray();
    int ngrd = bxa.length();
    int k;
    for (k = 0; k < ngrd; k++) {
	FineAdd(mflx[k],dir,k,srccomp,destcomp,numcomp,mult);
    }
*/
}

// -------------------------------------------------------------
void
FluxRegister::FineAdd(const MultiFab& mflx, const MultiFab& area, int dir,
		      int srccomp, int destcomp, int numcomp, REAL mult)
{
    BoxLib::Error("FineAdd(multifab, multifab, ...) not implemented");
/*
    const BoxArray& bxa = mflx.boxArray();
    int ngrd = bxa.length();
    int k;
    for (k = 0; k < ngrd; k++) {
	FineAdd(mflx[k],area[k],dir,k,srccomp,destcomp,numcomp,mult);
    }
*/
}

// -------------------------------------------------------------
void
FluxRegister::FineAdd(const FARRAYBOX& flux, int dir, int boxno,
		      int srccomp, int destcomp, int numcomp, REAL mult)
{
    assert(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);
    Box cbox(flux.box());
    cbox.coarsen(ratio);
    
    const Box  &flxbox = flux.box();
    const int  *flo    = flxbox.loVect();
    const int  *fhi    = flxbox.hiVect();
    const Real *flxdat = flux.dataPtr(srccomp);

    Orientation face_lo(dir,Orientation::low);
    FArrayBox& loreg = bndry[face_lo][boxno];
    const Box& lobox = loreg.box();
    assert(cbox.contains(lobox));
    const int* rlo = lobox.loVect();
    const int* rhi = lobox.hiVect();
    Real *lodat = loreg.dataPtr(destcomp);
    FORT_FRFINEADD(lodat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);

    Orientation face_hi(dir,Orientation::high);
    FArrayBox &hireg = bndry[face_hi][boxno];
    const Box &hibox = hireg.box();
    assert(cbox.contains(hibox));
    rlo = hibox.loVect();
    rhi = hibox.hiVect();
    Real *hidat = hireg.dataPtr(destcomp);
    FORT_FRFINEADD(hidat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);
}

// -------------------------------------------------------------
void
FluxRegister::FineAdd(const FARRAYBOX& flux, const FARRAYBOX& area,
                      int dir, int boxno,
		      int srccomp, int destcomp, int numcomp, REAL mult)
{
    BoxLib::Error("FineAdd(flux, area, ..., boxno, ...) not implemented");
/*
    int nvf = flux.nComp();
    assert(srccomp >= 0 && srccomp+numcomp <= nvf);
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);
    BOX cbox(flux.box());
    cbox.coarsen(ratio);
    
    const REAL* area_dat = area.dataPtr();
    const int* alo = area.loVect();
    const int* ahi = area.hiVect();
    
    const BOX& flxbox = flux.box();
    const int* flo = flxbox.loVect();
    const int* fhi = flxbox.hiVect();
    const REAL* flxdat = flux.dataPtr(srccomp);

    Orientation face_lo(dir,Orientation::low);
    FARRAYBOX& loreg = bndry[face_lo][boxno];
    const BOX& lobox = loreg.box();
    assert(cbox.contains(lobox));
    const int* rlo = lobox.loVect();
    const int* rhi = lobox.hiVect();
    REAL *lodat = loreg.dataPtr(destcomp);
    FORT_FRFAADD(lodat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);

    Orientation face_hi(dir,Orientation::high);
    FARRAYBOX& hireg = bndry[face_hi][boxno];
    const BOX& hibox = hireg.box();
    assert(cbox.contains(hibox));
    rlo = hibox.loVect();
    rhi = hibox.hiVect();
    REAL *hidat = hireg.dataPtr(destcomp);
    FORT_FRFAADD(hidat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);
*/
}

// -------------------------------------------------------------
/*
static void printFAB(ostream& os, const FARRAYBOX& f, int comp)
{

    const BOX& bx = f.box();
    BOX subbox(bx);
    os << "[box = " << subbox << ", comp = "
         << comp << ']' << '\n';
    const int* len = bx.length();
    const int* lo = bx.loVect();
    const int* s_len = subbox.length();
    const int* s_lo = subbox.loVect();
    const REAL* d = f.dataPtr(comp);
    char str[80];
    int j;
    for (j = 0; j < s_len[1]; j++) {
        int jrow = s_lo[1] + s_len[1]-1-j;
        const REAL* d_x = d + (jrow - lo[1])*len[0] + s_lo[0]-lo[0];
        sprintf(str,"%04d : ",jrow);
        os << str;
        int i;
        for (i = 0; i < s_len[0]; i++) {
            sprintf(str,"%18.12f ",d_x[i]);
            os << str;
        }
        os << '\n';
    }
}
*/

// -------------------------------------------------------------
void
FluxRegister::print(ostream &os)
{
    BoxLib::Error("FluxRegister::print() not implemented");
/*
    int ngrd = grids.length();
    os << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    os << "FluxRegister with coarse level = " << crseLevel() << '\n';
    int k;
    for (k = 0; k < ngrd; k++) {
        os << "  Registers surrounding coarsened box " << grids[k] << '\n';
        int comp;
        for (comp = 0; comp < ncomp; comp++) {
            for (OrientationIter face; face; ++face) {
                const FARRAYBOX& reg = bndry[face()][k];
                printFAB(os,reg,comp);
            }
        }
    }
    os << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << '\n';
*/
}
    

