//
// $Id: FluxRegister.cpp,v 1.62 2001-08-01 21:50:45 lijewski Exp $
//

#include <FluxRegister.H>
#include <Geometry.H>
#include <FLUXREG_F.H>
#include <ParallelDescriptor.H>
#include <Profiler.H>
#include <ccse-mpi.H>

#include <vector>

FluxRegister::FluxRegister ()
{
    fine_level = ncomp = -1;
    ratio = IntVect::TheUnitVector();
    ratio.scale(-1);
}

FluxRegister::FluxRegister (const BoxArray& fine_boxes, 
                            const IntVect&  ref_ratio,
                            int             fine_lev,
                            int             nvar)
{
    define(fine_boxes,ref_ratio,fine_lev,nvar);
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
FluxRegister::define (const BoxArray& fine_boxes, 
                      const IntVect&  ref_ratio,
                      int             fine_lev,
                      int             nvar)
{
    BL_ASSERT(fine_boxes.isDisjoint());
    BL_ASSERT(grids.size() == 0);

    ratio      = ref_ratio;
    fine_level = fine_lev;
    ncomp      = nvar;

    grids.define(fine_boxes);
    grids.coarsen(ratio);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        const Orientation lo_face(dir,Orientation::low);
        const Orientation hi_face(dir,Orientation::high);

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
        const FabSet& lofabs = bndry[Orientation(dir,Orientation::low)];
        const FabSet& hifabs = bndry[Orientation(dir,Orientation::high)];

        for (FabSetIter fsi(lofabs); fsi.isValid(); ++fsi)
        {
            sum += lofabs[fsi].sum(comp);
            sum -= hifabs[fsi].sum(comp);
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
    BL_ASSERT(dir >= 0 && dir < BL_SPACEDIM);

    const FabSet& lofabs = bndry[Orientation(dir,Orientation::low)];
    const FabSet& hifabs = bndry[Orientation(dir,Orientation::high)];

    lofabs.copyTo(flx,src_comp,dest_comp,num_comp);
    hifabs.copyTo(flx,src_comp,dest_comp,num_comp);
}

//
// Structure used by Reflux()s.
//

struct RF
{
    RF ()
        :
        m_fabidx(-1),
        m_fridx(-1),
        m_shifted(false) {}

    RF (int         fabidx,
        int         fridx,
        Orientation face,
        FillBoxId   fbid)
        :
        m_fabidx(fabidx),
        m_fridx(fridx),
        m_face(face),
        m_fbid(fbid),
        m_shifted(false) {}

    RF (const IntVect& iv,
        int            fabidx,
        int            fridx,
        Orientation    face,
        FillBoxId      fbid)
        :
        m_iv(iv),
        m_fabidx(fabidx),
        m_fridx(fridx),
        m_face(face),
        m_fbid(fbid),
        m_shifted(true) {}

    RF (const RF& rhs)
        :
        m_iv(rhs.m_iv),
        m_fabidx(rhs.m_fabidx),
        m_fridx(rhs.m_fridx),
        m_face(rhs.m_face),
        m_fbid(rhs.m_fbid),
        m_shifted(rhs.m_shifted) {}

    RF& operator= (const RF& rhs)
    {
        m_iv      = rhs.m_iv;
        m_fabidx  = rhs.m_fabidx;
        m_fridx   = rhs.m_fridx;
        m_face    = rhs.m_face;
        m_fbid    = rhs.m_fbid;
        m_shifted = rhs.m_shifted;
        return *this;
    }

    IntVect     m_iv;
    int         m_fabidx;
    int         m_fridx;
    Orientation m_face;
    FillBoxId   m_fbid;
    bool        m_shifted;
};

void
FluxRegister::Reflux (MultiFab&       S,
                      const MultiFab& volume,
                      Real            scale,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp, 
                      const Geometry& geom)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::Reflux(MultiFab&,...)");

    FabSetCopyDescriptor fscd;

    FabSetId fsid[2*BL_SPACEDIM];

    for (OrientationIter fi; fi; ++fi)
    {
        fsid[fi()] = fscd.RegisterFabSet(&bndry[fi()]);
    }

    std::vector<RF> RFs;
    Array<IntVect>  pshifts(27);

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        Real*       s_dat   = S[mfi].dataPtr(dest_comp);
        const int*  slo     = S[mfi].loVect();
        const int*  shi     = S[mfi].hiVect();
        const Real* vol_dat = volume[mfi].dataPtr();
        const int*  vlo     = volume[mfi].loVect();
        const int*  vhi     = volume[mfi].hiVect();
        //
        // Find flux register that intersect with this grid.
        //
        for (int k = 0; k < grids.size(); k++)
        {
            Box bx = BoxLib::grow(grids[k],1);

            if (bx.intersects(mfi.validbox()))
            {
                for (OrientationIter fi; fi; ++fi)
                {
                    //
                    // low(high) face of fine grid => high (low)
                    // face of the exterior coarse grid cell updated.
                    //
                    Box ovlp = mfi.validbox() & BoxLib::adjCell(grids[k],fi());

                    if (ovlp.ok())
                    {
                        FillBoxId fbid = fscd.AddBox(fsid[fi()],
                                                     bndry[fi()].box(k),
                                                     0,
                                                     k,
                                                     src_comp,
                                                     0,
                                                     num_comp);

                        RFs.push_back(RF(mfi.index(),k,fi(),fbid));
                    }
                }
            }
            //
            // Add periodic possibilities.
            //
            if (geom.isAnyPeriodic() && !geom.Domain().contains(bx))
            {
                geom.periodicShift(bx,mfi.validbox(),pshifts);

                for (int iiv = 0; iiv < pshifts.size(); iiv++)
                {
                    const IntVect& iv = pshifts[iiv];
                    S[mfi].shift(iv);
                    const int* slo = S[mfi].loVect();
                    const int* shi = S[mfi].hiVect();
                    //
                    // This is a funny situation.  I don't want to permanently
                    // change vol, but I need to do a shift on it.  I'll shift
                    // it back later, so the overall change is nil.  But to do
                    // this, I have to cheat and do a cast.  This is pretty 
                    // disgusting.
                    //
                    FArrayBox* cheatvol = const_cast<FArrayBox*>(&volume[mfi]);
                    BL_ASSERT(cheatvol != 0);
                    cheatvol->shift(iv);
                    const int* vlo = cheatvol->loVect();
                    const int* vhi = cheatvol->hiVect();
                    Box sftbox = mfi.validbox();
                    sftbox.shift(iv);
                    BL_ASSERT(bx.intersects(sftbox));

                    for (OrientationIter fi; fi; ++fi)
                    {
                        //
                        // low(high)  face of fine grid => high (low)
                        // face of the exterior coarse grid cell updated.
                        //
                        Box ovlp = sftbox & BoxLib::adjCell(grids[k],fi());

                        if (ovlp.ok())
                        {
                            FillBoxId fbid = fscd.AddBox(fsid[fi()],
                                                         bndry[fi()].box(k),
                                                         0,
                                                         k,
                                                         src_comp,
                                                         0,
                                                         num_comp);

                            RFs.push_back(RF(iv,mfi.index(),k,fi(),fbid));
                        }
                    }
                    S[mfi].shift(-iv);
                    cheatvol->shift(-iv);
                }
            }
        }
    }

    fscd.CollectData();

    const int MyProc = ParallelDescriptor::MyProc();

    FArrayBox reg;

    for (int i = 0; i < RFs.size(); i++)
    {
        const RF&        rf   = RFs[i];
        const FillBoxId& fbid = rf.m_fbid;

        BL_ASSERT(bndry[rf.m_face].box(rf.m_fridx) == fbid.box());
        BL_ASSERT(S.DistributionMap()[rf.m_fabidx] == MyProc);
        BL_ASSERT(volume.DistributionMap()[rf.m_fabidx] == MyProc);

        FArrayBox&       fab_S      = S[rf.m_fabidx];
        const FArrayBox& fab_volume = volume[rf.m_fabidx];
        Real*            s_dat      = fab_S.dataPtr(dest_comp);
        const int*       slo        = fab_S.loVect();
        const int*       shi        = fab_S.hiVect();
        const Real*      vol_dat    = fab_volume.dataPtr();
        Box              fine_face  = BoxLib::adjCell(grids[rf.m_fridx],rf.m_face);
        Real             mult       = rf.m_face.isLow() ? -scale : scale;
        const int*       rlo        = fine_face.loVect();
        const int*       rhi        = fine_face.hiVect();

        if (!rf.m_shifted)
        {
            Box ovlp = S.box(rf.m_fabidx) & fine_face;

            BL_ASSERT(ovlp.ok());

            reg.resize(fbid.box(), num_comp);
            fscd.FillFab(fsid[rf.m_face], fbid, reg);

            const Real* reg_dat = reg.dataPtr(0);
            const int*  vlo     = fab_volume.loVect();
            const int*  vhi     = fab_volume.hiVect();
            const int*  lo      = ovlp.loVect();
            const int*  hi      = ovlp.hiVect();

            FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                          vol_dat,ARLIM(vlo),ARLIM(vhi),
                          reg_dat,ARLIM(rlo),ARLIM(rhi),
                          lo,hi,&num_comp,&mult);
        }
        else
        {
            fab_S.shift(rf.m_iv);
            //
            // This is a funny situation.  I don't want to permanently
            // change vol, but I need to do a shift on it.  I'll shift
            // it back later, so the overall change is nil.  But to do
            // this, I have to cheat and do a cast.  This is pretty 
            // disgusting.
            //
            FArrayBox* cheatvol = const_cast<FArrayBox*>(&fab_volume);
            BL_ASSERT(cheatvol != 0);
            cheatvol->shift(rf.m_iv);
            Box sftbox = S.box(rf.m_fabidx);
            sftbox.shift(rf.m_iv);
            Box ovlp = sftbox & fine_face;

            BL_ASSERT(ovlp.ok());

            reg.resize(fbid.box(), num_comp);
            fscd.FillFab(fsid[rf.m_face], fbid, reg);

            const Real* reg_dat = reg.dataPtr(0);
            const int*  vlo      = cheatvol->loVect();
            const int*  vhi      = cheatvol->hiVect();
            const int*  lo       = ovlp.loVect();
            const int*  hi       = ovlp.hiVect();

            FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                          vol_dat,ARLIM(vlo),ARLIM(vhi),
                          reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                          &num_comp,&mult);
            fab_S.shift(-rf.m_iv);
            cheatvol->shift(-rf.m_iv);
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
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::Reflux(MultiFab&, Real,...)");

    const Real* dx = geom.CellSize();

    FabSetCopyDescriptor fscd;

    FabSetId fsid[2*BL_SPACEDIM];

    for (OrientationIter fi; fi; ++fi)
    {
        fsid[fi()] = fscd.RegisterFabSet(&bndry[fi()]);
    }

    std::vector<RF> RFs;
    Array<IntVect>  pshifts(27);

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        //
        // Find flux register that intersects with this grid.
        //
        for (int k = 0; k < grids.size(); k++)
        {
            Box bx = BoxLib::grow(grids[k],1);

            if (bx.intersects(mfi.validbox()))
            {
                for (OrientationIter fi; fi; ++fi)
                {
                    Box ovlp = mfi.validbox() & BoxLib::adjCell(grids[k],fi());

                    if (ovlp.ok())
                    {
                        FillBoxId fbid = fscd.AddBox(fsid[fi()],
                                                     bndry[fi()].box(k),
                                                     0,
                                                     k,
                                                     src_comp,
                                                     0,
                                                     num_comp);

                        RFs.push_back(RF(mfi.index(),k,fi(),fbid));
                    }
                }
            }
            //
            // Add periodic possibilities.
            //
            if (geom.isAnyPeriodic() && !geom.Domain().contains(bx))
            {
                geom.periodicShift(bx,mfi.validbox(),pshifts);

                for (int iiv = 0; iiv < pshifts.size(); iiv++)
                {
                    const IntVect& iv = pshifts[iiv];
                    S[mfi].shift(iv);
                    const int* slo = S[mfi].loVect();
                    const int* shi = S[mfi].hiVect();
                    Box sftbox     = mfi.validbox();
                    sftbox.shift(iv);
                    BL_ASSERT(bx.intersects(sftbox));

                    for (OrientationIter fi; fi; ++fi)
                    {
                        //
                        // low(high) face of fine grid => high (low)
                        // face of the exterior coarse grid cell updated.
                        //
                        Box ovlp = sftbox & BoxLib::adjCell(grids[k],fi());

                        if (ovlp.ok())
                        {
                            FillBoxId fbid = fscd.AddBox(fsid[fi()],
                                                         bndry[fi()].box(k),
                                                         0,
                                                         k,
                                                         src_comp,
                                                         0,
                                                         num_comp);

                            RFs.push_back(RF(iv,mfi.index(),k,fi(),fbid));
                        }
                    }
                    S[mfi].shift(-iv);
                }
            }
        }
    }

    fscd.CollectData();

    const int MyProc = ParallelDescriptor::MyProc();

    FArrayBox reg;

    for (int i = 0; i < RFs.size(); i++)
    {
        const RF& rf          = RFs[i];
        const FillBoxId& fbid = rf.m_fbid;

        BL_ASSERT(bndry[rf.m_face].box(rf.m_fridx) == fbid.box());
        BL_ASSERT(S.DistributionMap()[rf.m_fabidx] == MyProc);

        FArrayBox& fab_S     = S[rf.m_fabidx];
        Box        fine_face = BoxLib::adjCell(grids[rf.m_fridx],rf.m_face);
        Real       mult      = rf.m_face.isLow() ? -scale : scale;
        const int* rlo       = fine_face.loVect();
        const int* rhi       = fine_face.hiVect();
        Real*      s_dat     = fab_S.dataPtr(dest_comp);
        const int* slo       = fab_S.loVect();
        const int* shi       = fab_S.hiVect();

        if (!rf.m_shifted)
        {
            Box ovlp = S.box(rf.m_fabidx) & fine_face;

            BL_ASSERT(ovlp.ok());

            reg.resize(fbid.box(), num_comp);
            fscd.FillFab(fsid[rf.m_face], fbid, reg);

            const Real* reg_dat = reg.dataPtr(0);
            const int*  lo      = ovlp.loVect();
            const int*  hi      = ovlp.hiVect();

            FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                            reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                            &num_comp,&mult);
        }
        else
        {
            fab_S.shift(rf.m_iv);
            Box sftbox = S.box(rf.m_fabidx);
            sftbox.shift(rf.m_iv);
            Box ovlp = sftbox & fine_face;

            BL_ASSERT(ovlp.ok());

            reg.resize(fbid.box(), num_comp);
            fscd.FillFab(fsid[rf.m_face], fbid, reg);

            const Real* reg_dat = reg.dataPtr(0);
            const int*  lo      = ovlp.loVect();
            const int*  hi      = ovlp.hiVect();

            FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                            reg_dat,ARLIM(rlo),ARLIM(rhi),
                            lo,hi,&num_comp,&mult);

            fab_S.shift(-rf.m_iv);
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
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Orientation face_lo(dir,Orientation::low);
    const Orientation face_hi(dir,Orientation::high);

    MultiFabCopyDescriptor mfcd;

    MultiFabId mfid_mflx = mfcd.RegisterFabArray(const_cast<MultiFab*>(&mflx));
    MultiFabId mfid_area = mfcd.RegisterFabArray(const_cast<MultiFab*>(&area));

    std::vector<FillBoxId> fillBoxId_mflx, fillBoxId_area;

    for (FabSetIter mfi_lo(bndry[face_lo]); mfi_lo.isValid(); ++mfi_lo)
    {
        for (int k = 0; k < mflx.boxArray().size(); k++)
        {
            if (mfi_lo.fabbox().intersects(mflx.boxArray()[k]))
            {
                Box lobox = mfi_lo.fabbox() & mflx.boxArray()[k];

                fillBoxId_mflx.push_back(mfcd.AddBox(mfid_mflx,
                                                     lobox,
                                                     0,
                                                     k,
                                                     srccomp,
                                                     0,
                                                     numcomp));

                BL_ASSERT(fillBoxId_mflx.back().box() == lobox);
                //
                // Here we'll save the index into the FabSet.
                //
                fillBoxId_mflx.back().FabIndex(mfi_lo.index());

                fillBoxId_area.push_back(mfcd.AddBox(mfid_area,
                                                     lobox,
                                                     0,
                                                     k,
                                                     0,
                                                     0,
                                                     1));

                BL_ASSERT(fillBoxId_area.back().box() == lobox);
                //
                // Here we'll save the direction.
                //
                fillBoxId_area.back().FabIndex(Orientation::low);
            }

            if (bndry[face_hi].fabbox(mfi_lo.index()).intersects(mflx.boxArray()[k]))
            {
                Box hibox = bndry[face_hi].fabbox(mfi_lo.index()) & mflx.boxArray()[k];

                fillBoxId_mflx.push_back(mfcd.AddBox(mfid_mflx,
                                                     hibox,
                                                     0,
                                                     k,
                                                     srccomp,
                                                     0,
                                                     numcomp));

                BL_ASSERT(fillBoxId_mflx.back().box() == hibox);
                //
                // Here we'll save the index into the FabSet.
                //
                fillBoxId_mflx.back().FabIndex(mfi_lo.index());

                fillBoxId_area.push_back(mfcd.AddBox(mfid_area,
                                                     hibox,
                                                     0,
                                                     k,
                                                     0,
                                                     0,
                                                     1));

                BL_ASSERT(fillBoxId_area.back().box() == hibox);
                //
                // Here we'll save the direction.
                //
                fillBoxId_area.back().FabIndex(Orientation::high);
            }
        }
    }

    mfcd.CollectData();

    BL_ASSERT(fillBoxId_mflx.size() == fillBoxId_area.size());

    const int MyProc = ParallelDescriptor::MyProc();

    FArrayBox mflx_fab, area_fab;

    for (int i = 0; i < fillBoxId_mflx.size(); i++)
    {
        const FillBoxId& fbid_mflx = fillBoxId_mflx[i];
        const FillBoxId& fbid_area = fillBoxId_area[i];

        BL_ASSERT(fbid_mflx.box() == fbid_area.box());

        Orientation the_face(dir,Orientation::Side(fbid_area.FabIndex()));

        BL_ASSERT(the_face == face_lo || the_face == face_hi);

        mflx_fab.resize(fbid_mflx.box(), numcomp);
        mfcd.FillFab(mfid_mflx, fbid_mflx, mflx_fab);
        area_fab.resize(fbid_mflx.box(), 1);
        mfcd.FillFab(mfid_area, fbid_area, area_fab);

        FabSet&   fabset   = bndry[the_face];
        const int fabindex = fbid_mflx.FabIndex();

        BL_ASSERT(fabset.DistributionMap()[fabindex] == MyProc);

        FArrayBox&  fab      = fabset[fabindex];
        const int*  flo      = mflx_fab.box().loVect();
        const int*  fhi      = mflx_fab.box().hiVect();
        const Real* flx_dat  = mflx_fab.dataPtr();
        const int*  alo      = area_fab.box().loVect();
        const int*  ahi      = area_fab.box().hiVect();
        const Real* area_dat = area_fab.dataPtr();
        const int*  rlo      = fab.loVect();
        const int*  rhi      = fab.hiVect();
        Real*       lodat    = fab.dataPtr(destcomp);
        const int*  lo       = fbid_mflx.box().loVect();
        const int*  hi       = fbid_mflx.box().hiVect();
        FORT_FRCAINIT(lodat,ARLIM(rlo),ARLIM(rhi),
                      flx_dat,ARLIM(flo),ARLIM(fhi),
                      area_dat,ARLIM(alo),ARLIM(ahi),
                      lo,hi,&numcomp,&dir,&mult);
    }
}

//
// Helper function and data for CrseInit()/CrseInitFinish().
//

static Array<int>              CIMsgs;
static std::vector<FabComTag>  CITags;
static std::vector<FArrayBox*> CIFabs;

static
void
DoIt (Orientation        face,
      int                k,
      FabSet*            bndry,
      const Box&         bx,
      const FArrayBox&   flux,
      int                srccomp,
      int                destcomp,
      int                numcomp,
      Real               mult)
{
    const DistributionMapping& dMap = bndry[face].DistributionMap();

    if (ParallelDescriptor::MyProc() == dMap[k])
    {
        //
        // Local data.
        //
        bndry[face][k].copy(flux, bx, srccomp, bx, destcomp, numcomp);
        bndry[face][k].mult(mult, bx, destcomp, numcomp);
    }
#ifdef BL_USE_MPI
    else
    {
        BL_ASSERT(CIMsgs.size() == ParallelDescriptor::NProcs());

        FabComTag tag;

        tag.toProc   = dMap[k];
        tag.fabIndex = k;
        tag.box      = bx;
        tag.face     = face;
        tag.destComp = destcomp;
        tag.nComp    = numcomp;

        FArrayBox* fab = new FArrayBox(bx, numcomp);

        fab->copy(flux, bx, srccomp, bx, 0, numcomp);
        fab->mult(mult, bx, 0, numcomp);

        CITags.push_back(tag);
        CIFabs.push_back(fab);

        CIMsgs[dMap[k]]++;
    }
#endif /*BL_USE_MPI*/
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
    BL_ASSERT(flux.box().contains(subbox));
    BL_ASSERT(srccomp  >= 0 && srccomp+numcomp  <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    if (CIMsgs.size() == 0)
        CIMsgs.resize(ParallelDescriptor::NProcs(), 0);

    for (int k = 0; k < grids.size(); k++)
    {
        const Orientation lo(dir,Orientation::low);

        if (subbox.intersects(bndry[lo].box(k)))
        {
            Box lobox = bndry[lo].box(k) & subbox;

            DoIt(lo,k,bndry,lobox,flux,srccomp,destcomp,numcomp,mult);
        }
        const Orientation hi(dir,Orientation::high);

        if (subbox.intersects(bndry[hi].box(k)))
        {
            Box hibox = bndry[hi].box(k) & subbox;

            DoIt(hi,k,bndry,hibox,flux,srccomp,destcomp,numcomp,mult);
        }
    }
}

void
FluxRegister::CrseInitFinish ()
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::CrseInitFinish()");

    if (ParallelDescriptor::NProcs() == 1) return;

#ifdef BL_USE_MPI

    const int seqno_1 = ParallelDescriptor::SeqNum();
    const int seqno_2 = ParallelDescriptor::SeqNum();
    const int MyProc  = ParallelDescriptor::MyProc();
    const int NProcs  = ParallelDescriptor::NProcs();
    const int IOProc  = ParallelDescriptor::IOProcessorNumber();

    BL_ASSERT(CITags.size() == CIFabs.size());

    if (CIMsgs.size() == 0)
        CIMsgs.resize(ParallelDescriptor::NProcs(),0);

    BL_ASSERT(CIMsgs[MyProc] == 0);

    Array<int>         Rcvs(NProcs,0);
    Array<int>         indx(NProcs);
    Array<Real*>       fab_data(NProcs);
    Array<CommData>    senddata;
    Array<MPI_Status>  status(NProcs);
    Array<MPI_Request> req_cd(NProcs,MPI_REQUEST_NULL);
    Array<MPI_Request> req_data(NProcs,MPI_REQUEST_NULL);

    int rc, NumRcvs = 0, idx = 0, NWaits = 0;
    //
    // Set Rcvs[i] to # of blocks we expect to get from CPU i ...
    //
    for (int i = 0; i < NProcs; i++)
    {
        if ((rc = MPI_Gather(&CIMsgs[i],
                             1,
                             MPI_INT,
                             Rcvs.dataPtr(),
                             1,
                             MPI_INT,
                             i,
                             ParallelDescriptor::Communicator())) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }

    BL_ASSERT(Rcvs[MyProc] == 0);

    for (int i = 0; i < NProcs; i++)
        NumRcvs += Rcvs[i];

    Array<CommData> recvdata(NumRcvs);
    //
    // Make sure we can treat CommData as a stream of integers.
    //
    BL_ASSERT(sizeof(CommData) == CommData::DIM*sizeof(int));
    //
    // Post one receive for each chunk being sent by other CPUs.
    // This is the CommData describing the FAB data that will be sent.
    //
    for (int i = 0; i < NProcs; i++)
    {
        if (Rcvs[i] > 0)
        {
            NWaits++;

            if ((rc = MPI_Irecv(&recvdata[idx],
                                Rcvs[i] * CommData::DIM,
                                MPI_INT,
                                i,
                                seqno_1,
                                ParallelDescriptor::Communicator(),
                                &req_cd[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

            idx += Rcvs[i];
        }
    }

    BL_ASSERT(idx == NumRcvs);
    //
    // Now send the CommData.
    //
    for (int i = 0; i < NProcs; i++)
    {
        if (CIMsgs[i] > 0)
        {
            senddata.resize(CIMsgs[i]);

            int Processed = 0;

            for (int j = 0; j < CITags.size(); j++)
            {
                if (CITags[j].toProc == i)
                {
                    CommData data(CITags[j].face,
                                  CITags[j].fabIndex,
                                  MyProc,
                                  0,
                                  CITags[j].nComp,
                                  CITags[j].destComp,   // Store as srcComp()
                                  0,                    // Not used.
                                  CITags[j].box);

                    senddata[Processed++] = data;
                }
            }

            BL_ASSERT(Processed == CIMsgs[i]);

            if ((rc = MPI_Send(senddata.dataPtr(),
                               senddata.size() * CommData::DIM,
                               MPI_INT,
                               i,
                               seqno_1,
                               ParallelDescriptor::Communicator())) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
        }
    }
    //
    // Post one receive for data being sent by CPU i ...
    //
    for (int completed; NWaits > 0; NWaits -= completed)
    {
        if ((rc = MPI_Waitsome(NProcs,
                               req_cd.dataPtr(),
                               &completed,
                               indx.dataPtr(),
                               status.dataPtr())) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);

        for (int k = 0; k < completed; k++)
        {
            //
            // Got to figure out # of Reals to expect from this CPU.
            //
            const int Ncpu = indx[k];

            idx = 0;
            for (int j = 0; j < Ncpu; j++)
                idx += Rcvs[j];

            size_t N = 0;
            for (int j = 0; j < Rcvs[Ncpu]; j++)
                N += recvdata[idx+j].box().numPts() * recvdata[idx+j].nComp();

            BL_ASSERT(N < INT_MAX);
            BL_ASSERT(!(The_FAB_Arena == 0));

            fab_data[Ncpu] = static_cast<Real*>(The_FAB_Arena->alloc(N*sizeof(Real)));

            if ((rc = MPI_Irecv(fab_data[Ncpu],
                                int(N),
                                mpi_data_type(fab_data[Ncpu]),
                                Ncpu,
                                seqno_2,
                                ParallelDescriptor::Communicator(),
                                &req_data[Ncpu])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
        }
    }
    //
    // Send the agglomerated FAB data.
    //
    for (int i = 0; i < NProcs; i++)
    {
        if (CIMsgs[i] > 0)
        {
            size_t N = 0;

            for (int j = 0; j < CITags.size(); j++)
                if (CITags[j].toProc == i)
                    N += CITags[j].box.numPts() * CITags[j].nComp;

            BL_ASSERT(N < INT_MAX);
            BL_ASSERT(!(The_FAB_Arena == 0));

            Real* data = static_cast<Real*>(The_FAB_Arena->alloc(N*sizeof(Real)));
            Real* dptr = data;

            for (int j = 0; j < CITags.size(); j++)
            {
                if (CITags[j].toProc == i)
                {
                    BL_ASSERT(CITags[j].box == CIFabs[j]->box());
                    BL_ASSERT(CITags[j].nComp == CIFabs[j]->nComp());
                    int count = CITags[j].box.numPts() * CITags[j].nComp;
                    memcpy(dptr, CIFabs[j]->dataPtr(), count * sizeof(Real));
                    delete CIFabs[j];
                    CIFabs[j] = 0;
                    dptr += count;
                }
            }

            BL_ASSERT(data + N == dptr);

            if ((rc = MPI_Send(data,
                               int(N),
                               mpi_data_type(data),
                               i,
                               seqno_2,
                               ParallelDescriptor::Communicator())) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

            The_FAB_Arena->free(data);
        }
    }
    //
    // Now receive and unpack FAB data.
    //
    FArrayBox fab;

    NWaits = 0;
    for (int i = 0; i < NProcs; i++)
        if (req_data[i] != MPI_REQUEST_NULL)
            NWaits++;

    for (int completed; NWaits > 0; NWaits -= completed)
    {
        if ((rc = MPI_Waitsome(NProcs,
                               req_data.dataPtr(),
                               &completed,
                               indx.dataPtr(),
                               status.dataPtr())) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);

        for (int k = 0; k < completed; k++)
        {
            const int Ncpu = indx[k];

            idx = 0;
            for (int j = 0; j < Ncpu; j++)
                idx += Rcvs[j];

            int Processed = 0;

            Real* dptr = fab_data[Ncpu];

            BL_ASSERT(!(dptr == 0));

            for (int j = 0; j < Rcvs[Ncpu]; j++)
            {
                const CommData& cd = recvdata[idx+j];
                fab.resize(cd.box(),cd.nComp());
                int N = fab.box().numPts() * fab.nComp();
                BL_ASSERT(N < INT_MAX);
                memcpy(fab.dataPtr(), dptr, N * sizeof(Real));
                bndry[cd.face()][cd.fabindex()].copy(fab,
                                                     fab.box(),
                                                     0,
                                                     fab.box(),
                                                     cd.srcComp(),
                                                     cd.nComp());
                dptr += N;
                Processed++;
            }

            BL_ASSERT(Processed == Rcvs[Ncpu]);

            The_FAB_Arena->free(fab_data[Ncpu]);
        }
    }
    //
    // Null out vectors.
    //
    CIFabs.erase(CIFabs.begin(), CIFabs.end());
    CITags.erase(CITags.begin(), CITags.end());
    //
    // Zero out CIMsgs.
    //
    for (int i = 0; i < NProcs; i++) CIMsgs[i] = 0;
#endif
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
    for (MFIter mflxmfi(mflx); mflxmfi.isValid(); ++mflxmfi)
    {
        FineAdd(mflx[mflxmfi],dir,mflxmfi.index(),srccomp,destcomp,numcomp,mult);
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
    for (MFIter mflxmfi(mflx); mflxmfi.isValid(); ++mflxmfi)
    {
        FineAdd(mflx[mflxmfi],area[mflxmfi],dir,mflxmfi.index(),srccomp,destcomp,numcomp,mult);
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
#ifndef NDEBUG
    Box cbox = BoxLib::coarsen(flux.box(),ratio);
#endif
    const Box&  flxbox = flux.box();
    const int*  flo    = flxbox.loVect();
    const int*  fhi    = flxbox.hiVect();
    const Real* flxdat = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

    BL_ASSERT(cbox.contains(loreg.box()));
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFINEADD(lodat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

    BL_ASSERT(cbox.contains(hireg.box()));
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
#ifndef NDEBUG
    Box cbox = BoxLib::coarsen(flux.box(),ratio);
#endif
    const Real* area_dat = area.dataPtr();
    const int*  alo      = area.loVect();
    const int*  ahi      = area.hiVect();
    const Box&  flxbox   = flux.box();
    const int*  flo      = flxbox.loVect();
    const int*  fhi      = flxbox.hiVect();
    const Real* flxdat   = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

    BL_ASSERT(cbox.contains(loreg.box()));
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFAADD(lodat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

    BL_ASSERT(cbox.contains(hireg.box()));
    rlo = hireg.box().loVect();
    rhi = hireg.box().hiVect();
    Real* hidat = hireg.dataPtr(destcomp);
    FORT_FRFAADD(hidat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);
}
