//BL_COPYRIGHT_NOTICE

//
// $Id: FluxRegister.cpp,v 1.9 1998-03-30 20:45:51 lijewski Exp $
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

FluxRegister::FluxRegister ()
    :
    BndryRegister()
{
    fine_level = ncomp = -1;
    ratio = IntVect::TheUnitVector();
    ratio.scale(-1);
}

FluxRegister::FluxRegister (const BoxArray& fine_boxes, 
                            const IntVect&  ref_ratio,
                            int             fine_lev,
                            int             nvar)
    :
    BndryRegister()
{
    define(fine_boxes,ref_ratio,fine_lev,nvar);
}

void
FluxRegister::define (const BoxArray& fine_boxes, 
                      const IntVect&  ref_ratio,
                      int             fine_lev,
                      int             nvar)
{
    assert(fine_boxes.isDisjoint());
    assert(!grids.ready());

    ratio = ref_ratio;
    fine_level = fine_lev;
    ncomp = nvar;
    grids.define(fine_boxes);
    grids.coarsen(ratio);
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        Orientation lo_face(dir,Orientation::low);
        Orientation hi_face(dir,Orientation::high);
        IndexType typ(IndexType::TheCellType());
        typ.setType(dir,IndexType::NODE);
        BndryRegister::define(lo_face,typ,0,1,0,nvar);
        BndryRegister::define(hi_face,typ,0,1,0,nvar);
    }
}

FluxRegister::~FluxRegister () {}

Real
FluxRegister::SumReg (int comp) const
{
    Real sum = 0.0;
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        Orientation lo_face(dir,Orientation::low);
        Orientation hi_face(dir,Orientation::high);
        const FabSet &lofabs = bndry[lo_face];
        const FabSet &hifabs = bndry[hi_face];
        for (ConstFabSetIterator fsi(lofabs); fsi.isValid(false); ++fsi)
        {
            ConstDependentFabSetIterator dfsi(fsi, hifabs);
            sum += fsi().sum(comp);
            sum -= dfsi().sum(comp);
        }
    }
    ParallelDescriptor::ReduceRealSum(sum);
    return sum;
}

void
FluxRegister::copyTo (FArrayBox& flx,
                      int        dir,
                      int        src_comp,
                      int        dest_comp,
                      int        num_comp)
{
    assert( dir >= 0 && dir < BL_SPACEDIM);

    Orientation lo_face(dir,Orientation::low);
    const FabSet &lofabs = bndry[lo_face];
    lofabs.copyTo(flx,src_comp,dest_comp,num_comp);

    Orientation hi_face(dir,Orientation::high);
    const FabSet &hifabs = bndry[hi_face];
    hifabs.copyTo(flx,src_comp,dest_comp,num_comp);
}

void
FluxRegister::Reflux (MultiFab&       S,
                      const MultiFab& volume,
                      Real            scale,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp, 
                      const Geometry& geom)
{
    FabSetCopyDescriptor fscd(true);
    FabSetId fsid[2*BL_SPACEDIM];
    for (OrientationIter fi; fi; ++fi)
    {
      fsid[fi()] = fscd.RegisterFabSet(&(bndry[fi()]));
    }
    List<FillBoxId> fillBoxIdList;
    FillBoxId tempFillBoxId;

    const BoxArray& grd_boxes = S.boxArray();

    for (MultiFabIterator mfi(S); mfi.isValid(false); ++mfi)
    {
        int grd = mfi.index();  // Punt for now.
        DependentMultiFabIterator mfi_volume(mfi, volume);
        FArrayBox &s = mfi();
        assert(grd_boxes[mfi.index()] == mfi.validbox());
        const Box &s_box = mfi.validbox();
        Real *s_dat = s.dataPtr(dest_comp);
        const int *slo = s.loVect();
        const int *shi = s.hiVect();
        const FArrayBox &vol = mfi_volume();
        const Real *vol_dat = vol.dataPtr();
        const int *vlo = vol.loVect();
        const int *vhi = vol.hiVect();
        //
        // Find flux register that intersect with this grid.
        //
        for (int k = 0; k < grids.length(); k++)
        {
            const Box &reg_box = grids[k];
            Box bx(grow(reg_box,1));
            if (bx.intersects(s_box))
            {
                for (OrientationIter fi; fi; ++fi)
                {
                    Orientation face = fi();
                    Box fine_face(adjCell(reg_box,face));
                    //
                    // low(hight)  face of fine grid => high (low)
                    // face of the exterior crarse grid cell updated.
                    // adjust sign of scale accordingly.
                    //
                    Real mult = (face.isLow() ? -scale : scale);
                    Box ovlp(s_box);
                    ovlp &= fine_face;
                    if (ovlp.ok())
                    {
                        Box regBox(bndry[face].box(k));
                        BoxList unfilledBoxes(regBox.ixType());

                        tempFillBoxId = fscd.AddBox(fsid[face],regBox,
                                                    unfilledBoxes, src_comp,
                                                    dest_comp, num_comp);
                        fillBoxIdList.append(tempFillBoxId);
                    }
                }
            }
            //
            // Add periodic possibilities.
            //
            if (geom.isAnyPeriodic() && !geom.Domain().contains(bx))
            {
              Array<IntVect> pshifts(27);
              geom.periodicShift(bx,s_box,pshifts);
              for (int iiv=0; iiv<pshifts.length(); iiv++)
              {
                IntVect iv = pshifts[iiv];
                s.shift(iv);
                const int *slo = s.loVect();
                const int *shi = s.hiVect();
                //
                // This is a funny situation.  I don't want to permanently
                // change vol, but I need to do a shift on it.  I'll shift
                // it back later, so the overall change is nil.  But to do
                // this, I have to cheat and do a cast.  This is pretty 
                // disgusting.
                //
                FArrayBox &cheatvol = *(FArrayBox *)&vol;
                cheatvol.shift(iv);
                const int *vlo = cheatvol.loVect();
                const int *vhi = cheatvol.hiVect();
                Box s_box = grd_boxes[grd];
                D_TERM( s_box.shift(0,iv[0]);,
                        s_box.shift(1,iv[1]);,
                        s_box.shift(2,iv[2]); )
                if (!bx.intersects(s_box))
                    BoxLib::Error("FluxRegister::Reflux logic error");

                for (OrientationIter fi; fi; ++fi)
                {
                    Orientation face = fi();
                    Box fine_face(adjCell(reg_box,face));
                    //
                    // low(hight)  face of fine grid => high (low)
                    // face of the exterior crarse grid cell updated.
                    // adjust sign of scale accordingly.
                    //
                    Real mult = (face.isLow() ? -scale : scale);
                    Box ovlp(s_box);
                    ovlp &= fine_face;
                    if (ovlp.ok())
                    {
                        Box regBox(bndry[face].box(k));
                        BoxList unfilledBoxes(regBox.ixType());

                        tempFillBoxId = fscd.AddBox(fsid[face], regBox,
                                                    unfilledBoxes, src_comp,
                                                    dest_comp, num_comp);
                        fillBoxIdList.append(tempFillBoxId);
                    }
                }
                s.shift(-iv);
                cheatvol.shift(-iv);
              }
            }
        }
    }

    Array<FillBoxId> fillBoxId(fillBoxIdList.length());
    int ifbi = 0;
    for (ListIterator<FillBoxId> li(fillBoxIdList); li; ++li)
    {
        fillBoxId[ifbi] = li();
        ++ifbi;
    }
    fillBoxIdList.clear();

    fscd.CollectData();

    int overlapId = 0;

    for (MultiFabIterator mfi(S); mfi.isValid(false); ++mfi)
    {
        int grd = mfi.index();  // Punt for now.
        DependentMultiFabIterator mfi_volume(mfi, volume);
        FArrayBox &s = mfi();
        assert(grd_boxes[mfi.index()] == mfi.validbox());
        const Box &s_box = mfi.validbox();
        Real *s_dat = s.dataPtr(dest_comp);
        const int *slo = s.loVect();
        const int *shi = s.hiVect();
        const FArrayBox &vol = mfi_volume();
        const Real *vol_dat = vol.dataPtr();
        const int *vlo = vol.loVect();
        const int *vhi = vol.hiVect();
        //
        // Find flux register that intersect with this grid.
        //
        for (int k = 0; k < grids.length(); k++)
        {
            const Box &reg_box = grids[k];
            Box bx(grow(reg_box,1));
            if (bx.intersects(s_box))
            {
                for (OrientationIter fi; fi; ++fi)
                {
                    Orientation face = fi();
                    Box fine_face(adjCell(reg_box,face));
                    //
                    // low(hight)  face of fine grid => high (low)
                    // face of the exterior crarse grid cell updated.
                    // adjust sign of scale accordingly.
                    //
                    Real mult = (face.isLow() ? -scale : scale);
                    Box ovlp(s_box);
                    ovlp &= fine_face;
                    if (ovlp.ok())
                    {
                        FillBoxId fbid = fillBoxId[overlapId];
                        FArrayBox reg(fbid.box(), num_comp);
                        fscd.FillFab(fsid[face], fbid, reg);
                        if (src_comp != 0)
                        {
                            cerr << endl << "Check me:  FluxRegister::Reflux(MultiFab &S, const MultiFab,...)"
                                 << endl << endl;
                            // check that reg.dataPtr(0) is correct (should it be reg.dataPtr(src_comp))?
                            // vince
                        }
                        const Real *reg_dat = reg.dataPtr(0);

                        const int *rlo = fine_face.loVect();
                        const int *rhi = fine_face.hiVect();
                        const int *lo = ovlp.loVect();
                        const int *hi = ovlp.hiVect();
                        FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                                      vol_dat,ARLIM(vlo),ARLIM(vhi),
                                      reg_dat,ARLIM(rlo),ARLIM(rhi),
                                      lo,hi,&num_comp,&mult);
                        ++overlapId;
                    }
                }
            }
            //
            // Add periodic possibilities.
            //
            if (geom.isAnyPeriodic() && !geom.Domain().contains(bx))
            {
                Array<IntVect> pshifts(27);
                geom.periodicShift(bx,s_box,pshifts);
                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    IntVect iv = pshifts[iiv];
                    s.shift(iv);
                    const int* slo = s.loVect();
                    const int* shi = s.hiVect();
                    //
                    // This is a funny situation.  I don't want to permanently
                    // change vol, but I need to do a shift on it.  I'll shift
                    // it back later, so the overall change is nil.  But to do
                    // this, I have to cheat and do a cast.  This is pretty 
                    // disgusting.
                    //
                    FArrayBox &cheatvol = *(FArrayBox *)&vol;
                    cheatvol.shift(iv);
                    const int *vlo = cheatvol.loVect();
                    const int *vhi = cheatvol.hiVect();
                    Box s_box = grd_boxes[grd];
                    D_TERM( s_box.shift(0,iv[0]);,
                            s_box.shift(1,iv[1]);,
                            s_box.shift(2,iv[2]); )
                    if( ! bx.intersects(s_box) )
                        BoxLib::Error("FluxRegister::Reflux logic error");

                    for(OrientationIter fi; fi; ++fi) {
                        Orientation face = fi();
                        Box fine_face(adjCell(reg_box,face));
                        //
                        // low(hight)  face of fine grid => high (low)
                        // face of the exterior crarse grid cell updated.
                        // adjust sign of scale accordingly.
                        //
                        Real mult = (face.isLow() ? -scale : scale);
                        Box ovlp(s_box);
                        ovlp &= fine_face;
                        if (ovlp.ok())
                        {
                            FillBoxId fbid = fillBoxId[overlapId];
                            Box regBox(bndry[face].box(k));
                            assert(regBox == fbid.box());
                            FArrayBox reg(fbid.box(), num_comp);
                            fscd.FillFab(fsid[face], fbid, reg);
                            const Real *reg_dat = reg.dataPtr(0);

                            const int *rlo = fine_face.loVect();
                            const int *rhi = fine_face.hiVect();
                            const int *lo = ovlp.loVect();
                            const int *hi = ovlp.hiVect();
                            FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                                          vol_dat,ARLIM(vlo),ARLIM(vhi),
                                          reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                                          &num_comp,&mult);
                            ++overlapId;
                        }
                    }

                    s.shift(-iv);
                    cheatvol.shift(-iv);
                }
            }
        }
    }
}

void
FluxRegister::Reflux (MultiFab&       S,
                      Real            scale,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp, 
                      const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    FabSetCopyDescriptor fscd(true);
    FabSetId fsid[2*BL_SPACEDIM];
    for (OrientationIter fi; fi; ++fi)
    {
        fsid[fi()] = fscd.RegisterFabSet(&(bndry[fi()]));
    }
    List<FillBoxId> fillBoxIdList;
    FillBoxId tempFillBoxId;

    for (MultiFabIterator mfi(S); mfi.isValid(false); ++mfi)
    {
        const Box &s_box = mfi.validbox();
        //
        // Find flux register that intersect with this grid.
        //
        for (int k = 0; k < grids.length(); k++)
        {
            const Box &reg_box = grids[k];
            Box bx(grow(reg_box,1));
            if (bx.intersects(s_box))
            {
                for (OrientationIter fi; fi; ++fi)
                {
                    Orientation face = fi();
                    Box fine_face(adjCell(reg_box,face));
                    Box ovlp(s_box);
                    ovlp &= fine_face;
                    if (ovlp.ok())
                    {
                        Box regBox(bndry[face].box(k));
                        BoxList unfilledBoxes(regBox.ixType());

                        tempFillBoxId = fscd.AddBox(fsid[face],regBox,
                                                    unfilledBoxes, src_comp,
                                                    dest_comp, num_comp);
                        fillBoxIdList.append(tempFillBoxId);
                    }
                }
            }
            //
            // Add periodic possibilities.
            //
            if (geom.isAnyPeriodic() && !geom.Domain().contains(bx))
            {
                FArrayBox &s = mfi();
                Array<IntVect> pshifts(27);
                geom.periodicShift(bx,s_box,pshifts);
                for (int iiv=0; iiv<pshifts.length(); iiv++)
                {
                    IntVect iv = pshifts[iiv];
                    s.shift(iv);
                    const int *slo = s.loVect();
                    const int *shi = s.hiVect();
                    Box s_box = mfi.validbox();
                    D_TERM( s_box.shift(0,iv[0]);,
                            s_box.shift(1,iv[1]);,
                            s_box.shift(2,iv[2]); )
                        assert(bx.intersects(s_box));

                for (OrientationIter fi; fi; ++fi)
                {
                    Orientation face = fi();
                    Box fine_face(adjCell(reg_box,face));
                    //
                    // low(hight)  face of fine grid => high (low)
                    // face of the exterior coarse grid cell updated.
                    // adjust sign of scale accordingly.
                    //
                    Real mult = (face.isLow() ? -scale : scale);
                    Box ovlp(s_box);
                    ovlp &= fine_face;
                    if (ovlp.ok())
                    {
                        Box regBox(bndry[face].box(k));
                        BoxList unfilledBoxes(regBox.ixType());

                        tempFillBoxId = fscd.AddBox(fsid[face],regBox,unfilledBoxes,
                                                    src_comp, dest_comp, num_comp);
                        //assert(unfilledBoxes.isEmpty());
                        fillBoxIdList.append(tempFillBoxId);

                        //const FArrayBox &reg = bndry[face][k];
                        //const Real *reg_dat = reg.dataPtr(src_comp);
                        //const int *rlo = fine_face.loVect();
                        //const int *rhi = fine_face.hiVect();
                        //const int *lo = ovlp.loVect();
                        //const int *hi = ovlp.hiVect();
                        //FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                        //reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                        //&num_comp,&mult);
                    }
                }
                s.shift(-iv);
                }
            }
        }
    }

    Array<FillBoxId> fillBoxId(fillBoxIdList.length());
    int ifbi = 0;
    for (ListIterator<FillBoxId> li(fillBoxIdList); li; ++li)
    {
        fillBoxId[ifbi] = li();
        ++ifbi;
    }
    fillBoxIdList.clear();

    fscd.CollectData();

    int overlapId = 0;
    for (MultiFabIterator mfi(S); mfi.isValid(false); ++mfi)
    {
        const Box &s_box = mfi.validbox();
        //
        // Find flux register that intersect with this grid.
        //
        for (int k = 0; k < grids.length(); k++)
        {
            const Box &reg_box = grids[k];
            Box bx(grow(reg_box,1));
            if (bx.intersects(s_box))
            {
                for (OrientationIter fi; fi; ++fi)
                {
                    Orientation face = fi();
                    Box fine_face(adjCell(reg_box,face));
                    //
                    // low(hight)  face of fine grid => high (low)
                    // face of the exterior coarse grid cell updated.
                    // adjust sign of scale accordingly.
                    //
                    Real mult = (face.isLow() ? -scale : scale);
                    Box ovlp(s_box);
                    ovlp &= fine_face;
                    if (ovlp.ok())
                    {
                        FArrayBox &sfab = mfi();
                        Real *s_dat = sfab.dataPtr(dest_comp);
                        const int *slo = sfab.loVect();
                        const int *shi = sfab.hiVect();
                        FillBoxId fbid = fillBoxId[overlapId];
                        FArrayBox reg(fbid.box(), num_comp);
                        fscd.FillFab(fsid[face], fbid, reg);
if(src_comp != 0) {
  cerr << endl << "Check me:  FluxRegister::Reflux(MultiFab &S, const MultiFab,...)"
       << endl << endl;
  // check that reg.dataPtr(0) is correct (should it be reg.dataPtr(src_comp))?
  // vince
}
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
                }
            }
            //
            // Add periodic possibilities.
            //
            if (geom.isAnyPeriodic() && !geom.Domain().contains(bx))
            {
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

                    for (OrientationIter fi; fi; ++fi)
                    {
                        Orientation face = fi();
                        Box fine_face(adjCell(reg_box,face));
                        Real mult = (face.isLow() ? -scale : scale);
                        Box ovlp(s_box);
                        ovlp &= fine_face;
                        if (ovlp.ok())
                        {
                            FArrayBox &sfab = mfi();
                            Real *s_dat = sfab.dataPtr(dest_comp);
                            const int *slo = sfab.loVect();
                            const int *shi = sfab.hiVect();
                            FillBoxId fbid = fillBoxId[overlapId];
                            Box regBox(bndry[face].box(k));
                            assert(regBox == fbid.box());
                            FArrayBox reg(fbid.box(), num_comp);
                            fscd.FillFab(fsid[face], fbid, reg);
                            const Real *reg_dat = reg.dataPtr(0);
                            const int *rlo = fine_face.loVect();
                            const int *rhi = fine_face.hiVect();
                            const int *lo = ovlp.loVect();
                            const int *hi = ovlp.hiVect();
                            FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                                            reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                                            &num_comp,&mult);
                            ++overlapId;
                        }
                    }
                    s.shift(-iv);
                }
            }
        }
    }
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        const MultiFab& area,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult)
{
    assert(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);

    Orientation face_lo(dir,Orientation::low);
    Orientation face_hi(dir,Orientation::high);

    MultiFabCopyDescriptor mfcd(true);
    MultiFabId mfid_mflx = mfcd.RegisterFabArray((MultiFab *) &mflx);
    MultiFabId mfid_area = mfcd.RegisterFabArray((MultiFab *) &area);

    List<FillBoxId> fillBoxIdList_mflx;
    List<FillBoxId> fillBoxIdList_area;

    const BoxArray &bxa = mflx.boxArray();
    for (FabSetIterator mfi_bndry_lo(bndry[face_lo]); mfi_bndry_lo.isValid(false);
        ++mfi_bndry_lo)
    {
      DependentFabSetIterator mfi_bndry_hi(mfi_bndry_lo, bndry[face_hi]);

      for (int k = 0; k < bxa.length(); k++)
      {
        const Box subbox(bxa[k]);
        Box lobox(mfi_bndry_lo.fabbox());
        lobox &= subbox;
        if (lobox.ok())
        {
              BoxList unfilledBoxes(lobox.ixType());  // unused here
              FillBoxId fbid_mflx;
              fbid_mflx = mfcd.AddBox(mfid_mflx, mflx.fabbox(k),
                                      unfilledBoxes, 0, 0,
                                      mflx.nComp());
              fillBoxIdList_mflx.append(fbid_mflx);

              FillBoxId fbid_area;
              fbid_area = mfcd.AddBox(mfid_area, area.fabbox(k),
                                      unfilledBoxes, 0, 0,
                                      area.nComp());
              fillBoxIdList_area.append(fbid_area);
        }
        Box hibox(mfi_bndry_hi.fabbox());
        hibox &= subbox;
        if (hibox.ok())
        {
              BoxList unfilledBoxes(hibox.ixType());  // unused here
              FillBoxId fbid_mflx;
              fbid_mflx = mfcd.AddBox(mfid_mflx, mflx.fabbox(k),
                                      unfilledBoxes, 0, 0,
                                      mflx.nComp());
              fillBoxIdList_mflx.append(fbid_mflx);

              FillBoxId fbid_area;
              fbid_area = mfcd.AddBox(mfid_area, area.fabbox(k),
                                      unfilledBoxes, 0, 0,
                                      area.nComp());
              fillBoxIdList_area.append(fbid_area);
        }
      }
    }

    mfcd.CollectData();

    ListIterator<FillBoxId> fbidli_mflx(fillBoxIdList_mflx);
    ListIterator<FillBoxId> fbidli_area(fillBoxIdList_area);

    for (FabSetIterator mfi_bndry_lo(bndry[face_lo]); mfi_bndry_lo.isValid(false);
        ++mfi_bndry_lo)
    {
      DependentFabSetIterator mfi_bndry_hi(mfi_bndry_lo, bndry[face_hi]);

      for (int k = 0; k < bxa.length(); k++)
      {
        const Box subbox(bxa[k]);
        Box lobox(mfi_bndry_lo.fabbox());
        lobox &= subbox;
        if (lobox.ok())
        {
          assert(fbidli_mflx);
          FillBoxId fbid_mflx = fbidli_mflx();
          ++fbidli_mflx;
          FArrayBox mflx_fab(fbid_mflx.box(), mflx.nComp());
          mfcd.FillFab(mfid_mflx,  fbid_mflx, mflx_fab);

          assert(fbidli_area);
          FillBoxId fbid_area = fbidli_area();
          ++fbidli_area;
          FArrayBox area_fab(fbid_area.box(), area.nComp());
          mfcd.FillFab(mfid_area,  fbid_area, area_fab);

          const Box &flxbox = mflx_fab.box();
          const int *flo = flxbox.loVect();
          const int *fhi = flxbox.hiVect();
          const Real *flx_dat = mflx_fab.dataPtr(srccomp);

          const Box &areabox = area_fab.box();
          const int *alo = areabox.loVect();
          const int *ahi = areabox.hiVect();
          const Real *area_dat = area_fab.dataPtr();

          FArrayBox &loreg = mfi_bndry_lo();
          const int *rlo = loreg.loVect();
          const int *rhi = loreg.hiVect();
          Real *lodat = loreg.dataPtr(destcomp);
          const int *lo = lobox.loVect();
          const int *hi = lobox.hiVect();
          FORT_FRCAINIT(lodat,ARLIM(rlo),ARLIM(rhi),
                        flx_dat,ARLIM(flo),ARLIM(fhi),
                        area_dat,ARLIM(alo),ARLIM(ahi),
                        lo,hi,&numcomp,&dir,&mult);
        }
        Box hibox(mfi_bndry_hi.fabbox());
        hibox &= subbox;
        if (hibox.ok())
        {
          assert(fbidli_mflx);
          FillBoxId fbid_mflx = fbidli_mflx();
          ++fbidli_mflx;
          FArrayBox mflx_fab(fbid_mflx.box(), mflx.nComp());
          mfcd.FillFab(mfid_mflx,  fbid_mflx, mflx_fab);

          assert(fbidli_area);
          FillBoxId fbid_area = fbidli_area();
          ++fbidli_area;
          FArrayBox area_fab(fbid_area.box(), area.nComp());
          mfcd.FillFab(mfid_area,  fbid_area, area_fab);

          const Box &flxbox = mflx_fab.box();
          const int *flo = flxbox.loVect();
          const int *fhi = flxbox.hiVect();
          const Real *flx_dat = mflx_fab.dataPtr(srccomp);

          const Box &areabox = area_fab.box();
          const int *alo = areabox.loVect();
          const int *ahi = areabox.hiVect();
          const Real *area_dat = area_fab.dataPtr();

          FArrayBox &hireg = mfi_bndry_hi();
          const int *rlo = hireg.loVect();
          const int *rhi = hireg.hiVect();
          Real *hidat = hireg.dataPtr(destcomp);
          const int *lo = hibox.loVect();
          const int *hi = hibox.hiVect();
          FORT_FRCAINIT(hidat,ARLIM(rlo),ARLIM(rhi),
                        flx_dat,ARLIM(flo),ARLIM(fhi),
                        area_dat,ARLIM(alo),ARLIM(ahi),lo,hi,&numcomp,
                        &dir,&mult);
        }
      }
    }
}

void
FluxRegister::CrseInit (const FArrayBox& flux,
                        const Box&       subbox,
                        int              dir,
                        int              srccomp,
                        int              destcomp,
                        int              numcomp,
                        Real             mult)
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

    for (int k = 0; k < grids.length(); k++)
    {
        Orientation face_lo(dir,Orientation::low);
        Box lobox(bndry[face_lo].box(k));
        lobox &= subbox;
        if (lobox.ok())
        {
            const DistributionMapping &distributionMap = bndry[face_lo].DistributionMap();
            if (myproc == distributionMap[k])
            {
                //
                // Local
                //
                FArrayBox &loreg = bndry[face_lo][k];
                loreg.copy(flux, lobox, srccomp, lobox, destcomp, numcomp);
                loreg.mult(mult, lobox, destcomp, numcomp);
            }
            else
            {
                FArrayBox fabCom(lobox, numcomp);
                int fabComDestComp = 0;
                fabCom.copy(flux, lobox, srccomp, lobox, fabComDestComp,numcomp);
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
        if (hibox.ok())
        {
            const DistributionMapping &distributionMap = bndry[face_hi].DistributionMap();
            if (myproc == distributionMap[k])
            {
                //
                // Local.
                //
                FArrayBox &hireg = bndry[face_hi][k];
                hireg.copy(flux, hibox, srccomp, hibox, destcomp, numcomp);
                hireg.mult(mult, hibox, destcomp, numcomp);
            }
            else
            {
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

void
FluxRegister::CrseInitFinish ()
{
    FabComTag fabComTag;
    int tagSize = sizeof(FabComTag);
    ParallelDescriptor::SetMessageHeaderSize(tagSize);

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
                 << dataWaitingSize << " != " << shouldReceiveBytes
                 << '\n';
            BoxLib::Abort("Bad received nbytes");
        }
        if (!fabComTag.box.ok())
        {
            cerr << "Error in FluxRegister::CrseInitFinish():  "
                 << "bad fabComTag.box"
                 << '\n';
            BoxLib::Abort("Bad received box");
        }

        FArrayBox tempFab(fabComTag.box, fabComTag.nComp);
        ParallelDescriptor::ReceiveData(tempFab.dataPtr(),
                                        fabComTag.box.numPts()*fabComTag.nComp
                                        * sizeof(Real)); 
        int srcComp = 0;
        bndry[fabComTag.face][fabComTag.fabIndex].copy(tempFab, fabComTag.box,
                                                       srcComp, fabComTag.box,
                                                       fabComTag.destComp,
                                                       fabComTag.nComp);
    }
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
    for (ConstMultiFabIterator mflxmfi(mflx); mflxmfi.isValid(false); ++mflxmfi)
    {
        FineAdd(mflxmfi(),dir,mflxmfi.index(),srccomp,destcomp,numcomp,mult);
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
    for (ConstMultiFabIterator mflxmfi(mflx); mflxmfi.isValid(false); ++mflxmfi)
    {
        ConstDependentMultiFabIterator areamfi(mflxmfi, area);
        FineAdd(mflxmfi(),areamfi(),dir,mflxmfi.index(),
                srccomp,destcomp,numcomp,mult);
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
    int nvf = flux.nComp();
    assert(srccomp >= 0 && srccomp+numcomp <= nvf);
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);
    Box cbox(flux.box());
    cbox.coarsen(ratio);
    
    const Real* area_dat = area.dataPtr();
    const int* alo = area.loVect();
    const int* ahi = area.hiVect();
    
    const Box& flxbox = flux.box();
    const int* flo = flxbox.loVect();
    const int* fhi = flxbox.hiVect();
    const Real* flxdat = flux.dataPtr(srccomp);

    Orientation face_lo(dir,Orientation::low);
    FArrayBox& loreg = bndry[face_lo][boxno];
    const Box& lobox = loreg.box();
    assert(cbox.contains(lobox));
    const int* rlo = lobox.loVect();
    const int* rhi = lobox.hiVect();
    Real *lodat = loreg.dataPtr(destcomp);
    FORT_FRFAADD(lodat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);

    Orientation face_hi(dir,Orientation::high);
    FArrayBox& hireg = bndry[face_hi][boxno];
    const Box& hibox = hireg.box();
    assert(cbox.contains(hibox));
    rlo = hibox.loVect();
    rhi = hibox.hiVect();
    Real *hidat = hireg.dataPtr(destcomp);
    FORT_FRFAADD(hidat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);
}
