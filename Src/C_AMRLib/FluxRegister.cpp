//BL_COPYRIGHT_NOTICE

//
// $Id: FluxRegister.cpp,v 1.40 1998-06-21 17:56:23 lijewski Exp $
//

#include <FluxRegister.H>
#include <Geometry.H>
#include <FLUXREG_F.H>
#include <ParallelDescriptor.H>
#include <RunStats.H>
#include <Tracer.H>

#ifdef BL_USE_NEW_HFILES
#include <vector>
using std::vector;
#else
#include <vector.h>
#endif

#ifdef BL_USE_MPI
#include <mpi.h>
#endif

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
    assert(fine_boxes.isDisjoint());
    assert(!grids.ready());

    ratio      = ref_ratio;
    fine_level = fine_lev;
    ncomp      = nvar;

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
        const FabSet& lofabs = bndry[lo_face];
        const FabSet& hifabs = bndry[hi_face];
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
    assert(dir >= 0 && dir < BL_SPACEDIM);

    Orientation lo_face(dir,Orientation::low);
    const FabSet& lofabs = bndry[lo_face];
    lofabs.copyTo(flx,src_comp,dest_comp,num_comp);

    Orientation hi_face(dir,Orientation::high);
    const FabSet& hifabs = bndry[hi_face];
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
        Orientation face)
        :
        m_fabidx(fabidx),
        m_fridx(fridx),
        m_face(face),
        m_shifted(false) {}

    RF (const IntVect& iv,
        int            fabidx,
        int            fridx,
        Orientation    face)
        :
        m_iv(iv),
        m_fabidx(fabidx),
        m_fridx(fridx),
        m_face(face),
        m_shifted(true) {}

    IntVect     m_iv;
    int         m_fabidx;
    int         m_fridx;
    Orientation m_face;
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
    RunStats stats("reflux");

    stats.start();

    FabSetCopyDescriptor fscd;

    FabSetId fsid[2*BL_SPACEDIM];

    for (OrientationIter fi; fi; ++fi)
    {
        fsid[fi()] = fscd.RegisterFabSet(&bndry[fi()]);
    }

    vector<FillBoxId> fillBoxId;
    vector<RF>        RFs;
    Array<IntVect>    pshifts(27);

    for (MultiFabIterator mfi(S); mfi.isValid(false); ++mfi)
    {
        DependentMultiFabIterator mfi_volume(mfi, volume);

        Real* s_dat         = mfi().dataPtr(dest_comp);
        const int* slo      = mfi().loVect();
        const int* shi      = mfi().hiVect();
        const Real* vol_dat = mfi_volume().dataPtr();
        const int* vlo      = mfi_volume().loVect();
        const int* vhi      = mfi_volume().hiVect();
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
                        fillBoxId.push_back(fscd.AddBox(fsid[fi()],
                                                        bndry[fi()].box(k),
                                                        0,
                                                        k,
                                                        src_comp,
                                                        0,
                                                        num_comp));
                        //
                        // Push back a parallel RF for later use.
                        //
                        RFs.push_back(RF(mfi.index(),k,fi()));
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
                    IntVect iv = pshifts[iiv];
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
                    assert(cheatvol != 0);
                    cheatvol->shift(iv);
                    const int* vlo = cheatvol->loVect();
                    const int* vhi = cheatvol->hiVect();
                    Box sftbox = mfi.validbox();
                    sftbox.shift(iv);
                    assert(bx.intersects(sftbox));

                    for (OrientationIter fi; fi; ++fi)
                    {
                        //
                        // low(high)  face of fine grid => high (low)
                        // face of the exterior coarse grid cell updated.
                        //
                        Box ovlp = sftbox & ::adjCell(grids[k],fi());

                        if (ovlp.ok())
                        {
                            fillBoxId.push_back(fscd.AddBox(fsid[fi()],
                                                            bndry[fi()].box(k),
                                                            0,
                                                            k,
                                                            src_comp,
                                                            0,
                                                            num_comp));
                            //
                            // Push back a parallel RF for later use.
                            //
                            RFs.push_back(RF(iv,mfi.index(),k,fi()));
                        }
                    }
                    mfi().shift(-iv);
                    cheatvol->shift(-iv);
                }
            }
        }
    }

    fscd.CollectData();

    assert(fillBoxId.size() == RFs.size());

    const int MyProc = ParallelDescriptor::MyProc();

    FArrayBox reg;

    for (int i = 0; i < fillBoxId.size(); i++)
    {
        const FillBoxId& fbid = fillBoxId[i];
        const RF& rf          = RFs[i];

        assert(bndry[rf.m_face].box(rf.m_fridx) == fbid.box());
        assert(S.DistributionMap().ProcessorMap()[rf.m_fabidx] == MyProc);
        assert(volume.DistributionMap().ProcessorMap()[rf.m_fabidx] == MyProc);

        FArrayBox& fab_S            = S[rf.m_fabidx];
        const FArrayBox& fab_volume = volume[rf.m_fabidx];
        Real* s_dat                 = fab_S.dataPtr(dest_comp);
        const int* slo              = fab_S.loVect();
        const int* shi              = fab_S.hiVect();
        const Real* vol_dat         = fab_volume.dataPtr();
        Box fine_face               = ::adjCell(grids[rf.m_fridx],rf.m_face);
        Real mult                   = rf.m_face.isLow() ? -scale : scale;
        const int* rlo              = fine_face.loVect();
        const int* rhi              = fine_face.hiVect();

        if (!rf.m_shifted)
        {
            Box ovlp = S.box(rf.m_fabidx) & fine_face;

            assert(ovlp.ok());

            reg.resize(fbid.box(), num_comp);
            fscd.FillFab(fsid[rf.m_face], fbid, reg);

            const Real* reg_dat = reg.dataPtr(0);
            const int* vlo      = fab_volume.loVect();
            const int* vhi      = fab_volume.hiVect();
            const int* lo       = ovlp.loVect();
            const int* hi       = ovlp.hiVect();

            FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                          vol_dat,ARLIM(vlo),ARLIM(vhi),
                          reg_dat,ARLIM(rlo),ARLIM(rhi),
                          lo,hi,&num_comp,&mult);
        }
        else
        {
            IntVect iv = rf.m_iv;
            fab_S.shift(iv);
            //
            // This is a funny situation.  I don't want to permanently
            // change vol, but I need to do a shift on it.  I'll shift
            // it back later, so the overall change is nil.  But to do
            // this, I have to cheat and do a cast.  This is pretty 
            // disgusting.
            //
            FArrayBox* cheatvol = const_cast<FArrayBox*>(&fab_volume);
            assert(cheatvol != 0);
            cheatvol->shift(iv);
            Box sftbox = S.box(rf.m_fabidx);
            sftbox.shift(iv);
            Box ovlp = sftbox & fine_face;

            assert(ovlp.ok());

            reg.resize(fbid.box(), num_comp);
            fscd.FillFab(fsid[rf.m_face], fbid, reg);

            const Real* reg_dat = reg.dataPtr(0);
            const int* vlo      = cheatvol->loVect();
            const int* vhi      = cheatvol->hiVect();
            const int* lo       = ovlp.loVect();
            const int* hi       = ovlp.hiVect();

            FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                          vol_dat,ARLIM(vlo),ARLIM(vhi),
                          reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                          &num_comp,&mult);
            fab_S.shift(-iv);
            cheatvol->shift(-iv);
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
    RunStats stats("reflux");

    stats.start();

    const Real* dx = geom.CellSize();

    FabSetCopyDescriptor fscd;

    FabSetId fsid[2*BL_SPACEDIM];

    for (OrientationIter fi; fi; ++fi)
    {
        fsid[fi()] = fscd.RegisterFabSet(&bndry[fi()]);
    }

    vector<FillBoxId> fillBoxId;
    vector<RF>        RFs;
    Array<IntVect>    pshifts(27);

    for (MultiFabIterator mfi(S); mfi.isValid(false); ++mfi)
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
                        fillBoxId.push_back(fscd.AddBox(fsid[fi()],
                                                        bndry[fi()].box(k),
                                                        0,
                                                        k,
                                                        src_comp,
                                                        0,
                                                        num_comp));
                        //
                        // Push back a parallel RF for later use.
                        //
                        RFs.push_back(RF(mfi.index(),k,fi()));
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
                    IntVect iv = pshifts[iiv];
                    mfi().shift(iv);
                    const int* slo = mfi().loVect();
                    const int* shi = mfi().hiVect();
                    Box sftbox     = mfi.validbox();
                    sftbox.shift(iv);
                    assert(bx.intersects(sftbox));

                    for (OrientationIter fi; fi; ++fi)
                    {
                        //
                        // low(high) face of fine grid => high (low)
                        // face of the exterior coarse grid cell updated.
                        //
                        Box ovlp = sftbox & ::adjCell(grids[k],fi());

                        if (ovlp.ok())
                        {
                            fillBoxId.push_back(fscd.AddBox(fsid[fi()],
                                                            bndry[fi()].box(k),
                                                            0,
                                                            k,
                                                            src_comp,
                                                            0,
                                                            num_comp));
                            //
                            // Push back a parallel RF for later use.
                            //
                            RFs.push_back(RF(iv,mfi.index(),k,fi()));
                        }
                    }
                    mfi().shift(-iv);
                }
            }
        }
    }

    fscd.CollectData();

    assert(fillBoxId.size() == RFs.size());

    const int MyProc = ParallelDescriptor::MyProc();

    FArrayBox reg;

    for (int i = 0; i < fillBoxId.size(); i++)
    {
        const FillBoxId& fbid = fillBoxId[i];
        const RF& rf          = RFs[i];

        assert(bndry[rf.m_face].box(rf.m_fridx) == fbid.box());
        assert(S.DistributionMap().ProcessorMap()[rf.m_fabidx] == MyProc);

        FArrayBox& fab_S = S[rf.m_fabidx];
        Box fine_face    = ::adjCell(grids[rf.m_fridx],rf.m_face);
        Real mult        = rf.m_face.isLow() ? -scale : scale;
        const int* rlo   = fine_face.loVect();
        const int* rhi   = fine_face.hiVect();
        Real* s_dat      = fab_S.dataPtr(dest_comp);
        const int* slo   = fab_S.loVect();
        const int* shi   = fab_S.hiVect();

        if (!rf.m_shifted)
        {
            Box ovlp = S.box(rf.m_fabidx) & fine_face;

            assert(ovlp.ok());

            reg.resize(fbid.box(), num_comp);
            fscd.FillFab(fsid[rf.m_face], fbid, reg);

            const Real* reg_dat = reg.dataPtr(0);
            const int* lo       = ovlp.loVect();
            const int* hi       = ovlp.hiVect();

            FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                            reg_dat,ARLIM(rlo),ARLIM(rhi),lo,hi,
                            &num_comp,&mult);
        }
        else
        {
            IntVect iv = rf.m_iv;
            fab_S.shift(iv);
            Box sftbox = S.box(rf.m_fabidx);
            sftbox.shift(iv);
            Box ovlp = sftbox & fine_face;

            assert(ovlp.ok());

            reg.resize(fbid.box(), num_comp);
            fscd.FillFab(fsid[rf.m_face], fbid, reg);

            const Real* reg_dat = reg.dataPtr(0);
            const int* lo       = ovlp.loVect();
            const int* hi       = ovlp.hiVect();

            FORT_FRCVREFLUX(s_dat,ARLIM(slo),ARLIM(shi),dx,
                            reg_dat,ARLIM(rlo),ARLIM(rhi),
                            lo,hi,&num_comp,&mult);

            fab_S.shift(-iv);
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
    assert(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Orientation face_lo(dir,Orientation::low);
    const Orientation face_hi(dir,Orientation::high);

    MultiFabCopyDescriptor mfcd;

    MultiFabId mfid_mflx = mfcd.RegisterFabArray(const_cast<MultiFab*>(&mflx));
    MultiFabId mfid_area = mfcd.RegisterFabArray(const_cast<MultiFab*>(&area));

    vector<FillBoxId> fillBoxId_mflx, fillBoxId_area;

    for (FabSetIterator mfi_bndry_lo(bndry[face_lo]);
         mfi_bndry_lo.isValid(false); ++mfi_bndry_lo)
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

                assert(fillBoxId_mflx.back().box() == lobox);
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

                assert(fillBoxId_area.back().box() == lobox);
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

                assert(fillBoxId_mflx.back().box() == hibox);
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

                assert(fillBoxId_area.back().box() == hibox);
                //
                // Here we'll save the direction.
                //
                fillBoxId_area.back().FabIndex(Orientation::high);
            }
        }
    }

    mfcd.CollectData();

    assert(fillBoxId_mflx.size() == fillBoxId_area.size());

    const int MyProc = ParallelDescriptor::MyProc();

    FArrayBox mflx_fab, area_fab;

    for (int i = 0; i < fillBoxId_mflx.size(); i++)
    {
        const FillBoxId& fbid_mflx = fillBoxId_mflx[i];
        const FillBoxId& fbid_area = fillBoxId_area[i];
        assert(fbid_mflx.box() == fbid_area.box());

        Orientation the_face(dir,Orientation::Side(fbid_area.FabIndex()));
        assert(the_face == face_lo || the_face == face_hi);

        mflx_fab.resize(fbid_mflx.box(), numcomp);
        mfcd.FillFab(mfid_mflx, fbid_mflx, mflx_fab);
        area_fab.resize(fbid_mflx.box(), 1);
        mfcd.FillFab(mfid_area, fbid_area, area_fab);

        FabSet& fabset = bndry[the_face];
        int fabindex   = fbid_mflx.FabIndex();

        assert(fabset.DistributionMap().ProcessorMap()[fabindex] == MyProc);

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

#ifdef BL_USE_MPI
static Array<int>         CIMsgs;
static vector<FabComTag>  CITags;
static vector<FArrayBox*> CIFabs;
#endif

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
    else
    {
        FabComTag tag;

        tag.toProc   = dMap[k];
        tag.fabIndex = k;
        tag.box      = bx;
        tag.face     = face;
        tag.destComp = destcomp;
        tag.nComp    = numcomp;

#ifdef BL_USE_MPI
        assert(CIMsgs.length() == ParallelDescriptor::NProcs());

        FArrayBox* fabCom = new FArrayBox(bx, numcomp);

        fabCom->copy(flux, bx, srccomp, bx, 0, numcomp);
        fabCom->mult(mult, bx, 0, numcomp);

        CITags.push_back(tag);
        CIFabs.push_back(fabCom);

        CIMsgs[dMap[k]]++;
#else
        FArrayBox fabCom(bx, numcomp);

        fabCom.copy(flux, bx, srccomp, bx, 0, numcomp);
        fabCom.mult(mult, bx, 0, numcomp);

        ParallelDescriptor::SendData(dMap[k],
                                     &tag,
                                     fabCom.dataPtr(),
                                     bx.numPts() * numcomp * sizeof(Real));
#endif /*BL_USE_MPI*/
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
    TRACER("FluxRegister::CrseInit()");

    assert(flux.box().contains(subbox));
    assert(srccomp  >= 0 && srccomp+numcomp  <= flux.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);

#ifdef BL_USE_MPI
    if (CIMsgs.length() == 0)
    {
        CIMsgs.resize(ParallelDescriptor::NProcs(), 0);
    }
#endif

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
#ifdef BL_USE_MPI
    //
    // Pass each processor # of IRecv()s it'll need to post.
    //
    const int NProcs = ParallelDescriptor::NProcs();
    const int MyProc = ParallelDescriptor::MyProc();

    assert(CITags.size() == CIFabs.size());

    if (CIMsgs.length() == 0)
    {
        CIMsgs.resize(NProcs, 0);
    }

    int rc;

    Array<int> nrcv(NProcs, 0);

    for (int i = 0; i < NProcs; i++)
    {
        if ((rc = MPI_Reduce(&CIMsgs[i],
                             &nrcv[i],
                             1,
                             MPI_INT,
                             MPI_SUM,
                             i,
                             MPI_COMM_WORLD)) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }

    const int NumRecv = nrcv[MyProc];

    Array<MPI_Request> reqs(NumRecv);
    Array<MPI_Status>  stat(NumRecv);
    Array<CommData>    recv(NumRecv);
    PArray<FArrayBox>  fabs(NumRecv,PArrayManage);
    //
    // First receive/send the box information.
    // I'll receive the NumRecv boxes in any order.
    //
    for (int i = 0; i < NumRecv; i++)
    {
        if ((rc = MPI_Irecv(recv[i].dataPtr(),
                            recv[i].length(),
                            MPI_INT,
                            MPI_ANY_SOURCE,
                            711,
                            MPI_COMM_WORLD,
                            &reqs[i])) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }

    for (int i = 0; i < CITags.size(); i++)
    {
        CommData senddata(CITags[i].face,
                          CITags[i].fabIndex,
                          MyProc,
                          //
                          // We use the index into loop over CITags as the ID.
                          // The combination of the loop index and the
                          // processor from which the message was sent forms
                          // a unique identifier.  We'll later use the
                          // combination of fromproc() and id() to match up
                          // the box()s being sent now with the FAB data on
                          // those box()s to be sent next.
                          //
                          i,
                          CITags[i].nComp,
                          CITags[i].destComp, // Store as srcComp() component.
                          0,                // Not Used.
                          CITags[i].box);

        if ((rc = MPI_Ssend(senddata.dataPtr(),
                            senddata.length(),
                            MPI_INT,
                            CITags[i].toProc,
                            711,
                            MPI_COMM_WORLD)) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }

    if ((rc = MPI_Waitall(NumRecv,
                          reqs.dataPtr(),
                          stat.dataPtr())) != MPI_SUCCESS)
        ParallelDescriptor::Abort(rc);
    //
    // Now the FAB data itself.
    //
    for (int i = 0; i < NumRecv; i++)
    {
        fabs.set(i, new FArrayBox(recv[i].box(), recv[i].nComp()));

        if ((rc = MPI_Irecv(fabs[i].dataPtr(),
                            fabs[i].box().numPts() * recv[i].nComp(),
                            mpi_data_type(fabs[i].dataPtr()),
                            recv[i].fromproc(),
                            recv[i].id(),
                            MPI_COMM_WORLD,
                            &reqs[i])) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }

    for (int i = 0; i < CITags.size(); i++)
    {
        long count = CITags[i].box.numPts() * CITags[i].nComp;

        assert(count < INT_MAX);
        assert(CITags[i].box == CIFabs[i]->box());
        assert(CITags[i].nComp == CIFabs[i]->nComp());
        //
        // Use MPI_Ssend() to try and force the system not to buffer.
        //
        if ((rc = MPI_Ssend(CIFabs[i]->dataPtr(),
                            int(count),
                            mpi_data_type(CIFabs[i]->dataPtr()),
                            CITags[i].toProc,
                            //
                            // We use the index into loop over CITags as ID.
                            // The combination of the loop index and the
                            // processor from which the message was sent forms
                            // a unique identifier.
                            //
                            // Note that the form of this MPI_Ssend() MUST
                            // match the MPI_Send() of the box()es
                            // corresponding to this FAB above.
                            //
                            i,
                            MPI_COMM_WORLD)) != MPI_SUCCESS)
            ParallelDescriptor::Abort(rc);
    }

    if ((rc = MPI_Waitall(NumRecv,
                          reqs.dataPtr(),
                          stat.dataPtr())) != MPI_SUCCESS)
        ParallelDescriptor::Abort(rc);

    for (int i = 0; i < NumRecv; i++)
    {
        bndry[recv[i].face()][recv[i].fabindex()].copy(fabs[i],
                                                       fabs[i].box(),
                                                       0,
                                                       fabs[i].box(),
                                                       recv[i].srcComp(),
                                                       recv[i].nComp());
    }
    //
    // Delete buffered FABs.
    //
    for (int i = 0; i < CIFabs.size(); i++)
    {
        delete CIFabs[i];
    }
    //
    // Null out vectors.
    //
    CIFabs.erase(CIFabs.begin(), CIFabs.end());
    CITags.erase(CITags.begin(), CITags.end());
    //
    // Zero out CIMsgs.  It's size will not need to change.
    //
    for (int i = 0; i < NProcs; i++)
    {
        CIMsgs[i] = 0;
    }
#else
    FabComTag tag;

    ParallelDescriptor::SetMessageHeaderSize(sizeof(FabComTag));

    int dataWaitingSize;
    while (ParallelDescriptor::GetMessageHeader(dataWaitingSize, &tag))
    {
        long t_long = tag.box.numPts() * tag.nComp * sizeof(Real);

        assert(t_long < INT_MAX);
        assert(dataWaitingSize == int(t_long));
        assert(tag.box.ok());

        FArrayBox tempFab(tag.box, tag.nComp);

        ParallelDescriptor::ReceiveData(tempFab.dataPtr(), int(t_long));

        bndry[tag.face][tag.fabIndex].copy(tempFab,
                                           tag.box,
                                           0,
                                           tag.box,
                                           tag.destComp,
                                           tag.nComp);
    }
#endif /*BL_USE_MPI*/
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
#ifndef NDEBUG
    Box cbox = ::coarsen(flux.box(),ratio);
#endif
    const Box&  flxbox = flux.box();
    const int*  flo    = flxbox.loVect();
    const int*  fhi    = flxbox.hiVect();
    const Real* flxdat = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

    assert(cbox.contains(loreg.box()));
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFINEADD(lodat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

    assert(cbox.contains(hireg.box()));
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
    assert(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    assert(destcomp >= 0 && destcomp+numcomp <= ncomp);
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

    assert(cbox.contains(loreg.box()));
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFAADD(lodat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

    assert(cbox.contains(hireg.box()));
    rlo = hireg.box().loVect();
    rhi = hireg.box().hiVect();
    Real* hidat = hireg.dataPtr(destcomp);
    FORT_FRFAADD(hidat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);
}
