
//
// $Id: FluxRegister.cpp,v 1.59 2000-10-02 20:48:42 lijewski Exp $
//

#include <FluxRegister.H>
#include <Geometry.H>
#include <FLUXREG_F.H>
#include <ParallelDescriptor.H>
#include <RunStats.H>

#ifdef BL_USE_NEW_HFILES
#include <vector>
using std::vector;
#else
#include <vector.h>
#endif

#ifdef BL3_PROFILING
#include <BoxLib3/Profiler.H>
#endif

#include <ccse-mpi.H>

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

void
FluxRegister::define (const BoxArray& fine_boxes, 
                      const IntVect&  ref_ratio,
                      int             fine_lev,
                      int             nvar)
{
    BL_ASSERT(fine_boxes.isDisjoint());
    BL_ASSERT(!grids.ready());

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

        for (ConstFabSetIterator fsi(lofabs); fsi.isValid(); ++fsi)
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
#ifdef BL3_PROFILING
  BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::Reflux(MultiFab&,...)");
#endif
    static RunStats stats("reflux");

    stats.start();

    FabSetCopyDescriptor fscd;

    FabSetId fsid[2*BL_SPACEDIM];

    for (OrientationIter fi; fi; ++fi)
    {
        fsid[fi()] = fscd.RegisterFabSet(&bndry[fi()]);
    }

    vector<RF>     RFs;
    Array<IntVect> pshifts(27);

    for (MultiFabIterator mfi(S); mfi.isValid(); ++mfi)
    {
        DependentMultiFabIterator mfi_volume(mfi, volume);

        Real*       s_dat   = mfi().dataPtr(dest_comp);
        const int*  slo     = mfi().loVect();
        const int*  shi     = mfi().hiVect();
        const Real* vol_dat = mfi_volume().dataPtr();
        const int*  vlo     = mfi_volume().loVect();
        const int*  vhi     = mfi_volume().hiVect();
        //
        // Find flux register that intersect with this grid.
        //
        for (int k = 0; k < grids.length(); k++)
        {
            Box bx = ::grow(grids[k],1);

            if (bx.intersects(mfi.validbox()))
            {
                for (OrientationIter fi; fi; ++fi)
                {
                    //
                    // low(high) face of fine grid => high (low)
                    // face of the exterior coarse grid cell updated.
                    //
                    Box ovlp = mfi.validbox() & ::adjCell(grids[k],fi());

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

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    const IntVect& iv = pshifts[iiv];
                    mfi().shift(iv);
                    const int* slo = mfi().loVect();
                    const int* shi = mfi().hiVect();
                    //
                    // This is a funny situation.  I don't want to permanently
                    // change vol, but I need to do a shift on it.  I'll shift
                    // it back later, so the overall change is nil.  But to do
                    // this, I have to cheat and do a cast.  This is pretty 
                    // disgusting.
                    //
                    FArrayBox* cheatvol = const_cast<FArrayBox*>(&mfi_volume());
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
                        Box ovlp = sftbox & ::adjCell(grids[k],fi());

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
                    mfi().shift(-iv);
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
        Box              fine_face  = ::adjCell(grids[rf.m_fridx],rf.m_face);
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

    stats.end();
}

void
FluxRegister::Reflux (MultiFab&       S,
                      Real            scale,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp, 
                      const Geometry& geom)
{
#ifdef BL3_PROFILING
  BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::Reflux(MultiFab&, Real,...)");
#endif
    static RunStats stats("reflux");

    stats.start();

    const Real* dx = geom.CellSize();

    FabSetCopyDescriptor fscd;

    FabSetId fsid[2*BL_SPACEDIM];

    for (OrientationIter fi; fi; ++fi)
    {
        fsid[fi()] = fscd.RegisterFabSet(&bndry[fi()]);
    }

    vector<RF>     RFs;
    Array<IntVect> pshifts(27);

    for (MultiFabIterator mfi(S); mfi.isValid(); ++mfi)
    {
        //
        // Find flux register that intersects with this grid.
        //
        for (int k = 0; k < grids.length(); k++)
        {
            Box bx = ::grow(grids[k],1);

            if (bx.intersects(mfi.validbox()))
            {
                for (OrientationIter fi; fi; ++fi)
                {
                    Box ovlp = mfi.validbox() & ::adjCell(grids[k],fi());

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

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    const IntVect& iv = pshifts[iiv];
                    mfi().shift(iv);
                    const int* slo = mfi().loVect();
                    const int* shi = mfi().hiVect();
                    Box sftbox     = mfi.validbox();
                    sftbox.shift(iv);
                    BL_ASSERT(bx.intersects(sftbox));

                    for (OrientationIter fi; fi; ++fi)
                    {
                        //
                        // low(high) face of fine grid => high (low)
                        // face of the exterior coarse grid cell updated.
                        //
                        Box ovlp = sftbox & ::adjCell(grids[k],fi());

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
                    mfi().shift(-iv);
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
        Box        fine_face = ::adjCell(grids[rf.m_fridx],rf.m_face);
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

    stats.end();
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

    vector<FillBoxId> fillBoxId_mflx, fillBoxId_area;

    for (FabSetIterator mfi_bndry_lo(bndry[face_lo]);
         mfi_bndry_lo.isValid(); ++mfi_bndry_lo)
    {
        DependentFabSetIterator mfi_bndry_hi(mfi_bndry_lo, bndry[face_hi]);

        for (int k = 0; k < mflx.boxArray().length(); k++)
        {
            if (mfi_bndry_lo.fabbox().intersects(mflx.boxArray()[k]))
            {
                Box lobox = mfi_bndry_lo.fabbox() & mflx.boxArray()[k];

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
                fillBoxId_mflx.back().FabIndex(mfi_bndry_lo.index());

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
            if (mfi_bndry_hi.fabbox().intersects(mflx.boxArray()[k]))
            {
                Box hibox = mfi_bndry_hi.fabbox() & mflx.boxArray()[k];

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
                fillBoxId_mflx.back().FabIndex(mfi_bndry_hi.index());

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

static Array<int>         CIMsgs;
static vector<FabComTag>  CITags;
static vector<FArrayBox*> CIFabs;

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
        BL_ASSERT(CIMsgs.length() == ParallelDescriptor::NProcs());

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

    if (CIMsgs.length() == 0)
        CIMsgs.resize(ParallelDescriptor::NProcs(), 0);

    for (int k = 0; k < grids.length(); k++)
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
#ifdef BL3_PROFILING
  BL3_PROFILE(BL3_PROFILE_THIS_NAME() + "::CrseInitFinish()");
#endif
    if (ParallelDescriptor::NProcs() == 1)
        return;

#ifdef BL_USE_MPI
    static RunStats mpi_recv("mpi_recv");
    static RunStats mpi_send("mpi_send");
    static RunStats mpi_wait("mpi_wait");
    static RunStats mpi_gath("mpi_gather");
    static RunStats mpi_redu("mpi_reduce");
    static RunStats mpi_stat("crse_init_finish");

    mpi_stat.start();

    const int MyProc = ParallelDescriptor::MyProc();

    BL_ASSERT(CITags.size() == CIFabs.size());

    if (CIMsgs.length() == 0)
        CIMsgs.resize(ParallelDescriptor::NProcs(),0);

    BL_ASSERT(CIMsgs[MyProc] == 0);

    int rc;

    Array<int> Rcvs(ParallelDescriptor::NProcs(),0);
    //
    // Set Rcvs[i] to # of blocks we expect to get from CPU i ...
    //
    mpi_gath.start();
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
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
    mpi_gath.end();

    BL_ASSERT(Rcvs[MyProc] == 0);

    int NumRcvs = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        NumRcvs += Rcvs[i];

    Array<MPI_Request> req_cd(ParallelDescriptor::NProcs());
    Array<CommData>    rcv_cd(NumRcvs);
    //
    // Make sure we can treat CommData as a stream of integers.
    //
    BL_ASSERT(sizeof(CommData) == CommData::DIM*sizeof(int));
    //
    // Post one receive for each chunk being sent by other CPUs.
    // This is the CommData describing the FAB data that will be sent.
    //
    int idx = 0;

    mpi_recv.start();
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Rcvs[i] > 0)
        {
            if ((rc = MPI_Irecv(&rcv_cd[idx],
                                Rcvs[i] * CommData::DIM,
                                MPI_INT,
                                i,
                                741,
                                ParallelDescriptor::Communicator(),
                                &req_cd[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);

            idx += Rcvs[i];
        }
    }
    mpi_recv.end();

    BL_ASSERT(idx == NumRcvs);
    //
    // Now send the CommData.
    //
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (CIMsgs[i] > 0)
        {
            Array<CommData> senddata(CIMsgs[i]);

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
            //
            // Use MPI_Ssend() to try and force the system not to buffer.
            //
            mpi_send.start();
            if ((rc = MPI_Ssend(senddata.dataPtr(),
                                senddata.length() * CommData::DIM,
                                MPI_INT,
                                i,
                                741,
                                ParallelDescriptor::Communicator())) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_send.end();
        }
    }
    //
    // Post one receive for data being sent by CPU i ...
    //
    MPI_Status status;

    Array<Real*>       fab_data(ParallelDescriptor::NProcs());
    Array<MPI_Request> req_data(ParallelDescriptor::NProcs());

    idx = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Rcvs[i] > 0)
        {
            mpi_wait.start();
            if ((rc = MPI_Wait(&req_cd[i], &status)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_wait.end();
            //
            // Got to figure out # of Reals to expect from this CPU.
            //
            size_t N = 0;

            for (int j = 0; j < Rcvs[i]; j++)
                N += rcv_cd[idx+j].box().numPts() * rcv_cd[idx+j].nComp();

            BL_ASSERT(N < INT_MAX);
            BL_ASSERT(!(The_FAB_Arena == 0));

            fab_data[i] = static_cast<Real*>(The_FAB_Arena->alloc(N*sizeof(Real)));

            mpi_recv.start();
            if ((rc = MPI_Irecv(fab_data[i],
                                int(N),
                                mpi_data_type(fab_data[i]),
                                i,
                                719,
                                ParallelDescriptor::Communicator(),
                                &req_data[i])) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_recv.end();

            idx += Rcvs[i];
        }
    }

    BL_ASSERT(idx == NumRcvs);
    //
    // Send the agglomerated FAB data.
    //
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
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
            //
            // Use MPI_Ssend() to try and force the system not to buffer.
            //
            mpi_send.start();
            if ((rc = MPI_Ssend(data,
                                int(N),
                                mpi_data_type(data),
                                i,
                                719,
                                ParallelDescriptor::Communicator())) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_send.end();

            BL_ASSERT(!(The_FAB_Arena == 0));

            The_FAB_Arena->free(data);
        }
    }
    //
    // Now receive and unpack FAB data.
    //
    FArrayBox fab;

    idx = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (Rcvs[i] > 0)
        {
            mpi_wait.start();
            if ((rc = MPI_Wait(&req_data[i], &status)) != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
            mpi_wait.end();

            int Processed = 0;

            Real* dptr = fab_data[i];

            BL_ASSERT(!(dptr == 0));

            for (int j = 0; j < Rcvs[i]; j++)
            {
                const CommData& cd = rcv_cd[idx+j];

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

            BL_ASSERT(Processed == Rcvs[i]);
            BL_ASSERT(!(The_FAB_Arena == 0));

            The_FAB_Arena->free(fab_data[i]);

            fab_data[i] = 0;

            idx += Rcvs[i];
        }
    }

    BL_ASSERT(idx == NumRcvs);
    //
    // Null out vectors.
    //
    CIFabs.erase(CIFabs.begin(), CIFabs.end());
    CITags.erase(CITags.begin(), CITags.end());
    //
    // Zero out CIMsgs.
    //
    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        CIMsgs[i] = 0;

    mpi_stat.end();
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
    for (ConstMultiFabIterator mflxmfi(mflx); mflxmfi.isValid(); ++mflxmfi)
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
    for (ConstMultiFabIterator mflxmfi(mflx); mflxmfi.isValid(); ++mflxmfi)
    {
        ConstDependentMultiFabIterator areamfi(mflxmfi, area);

        FineAdd(mflxmfi(),
                areamfi(),
                dir,
                mflxmfi.index(),
                srccomp,
                destcomp,
                numcomp,
                mult);
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
    Box cbox = ::coarsen(flux.box(),ratio);
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
    Box cbox = ::coarsen(flux.box(),ratio);
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
