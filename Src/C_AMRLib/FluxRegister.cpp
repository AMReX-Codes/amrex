
#include <winstd.H>
#include <BArena.H>
#include <FluxRegister.H>
#include <Geometry.H>
#include <FLUXREG_F.H>
#include <ParallelDescriptor.H>
#include <ccse-mpi.H>

#include <deque>
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

FluxRegister::FluxRegister (const BoxArray&            fine_boxes, 
                            const IntVect&             ref_ratio,
                            int                        fine_lev,
                            int                        nvar,
                            const DistributionMapping& dm)
{
    define(fine_boxes,ref_ratio,fine_lev,nvar,dm);
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

void
FluxRegister::define (const BoxArray&            fine_boxes, 
                      const IntVect&             ref_ratio,
                      int                        fine_lev,
                      int                        nvar,
                      const DistributionMapping& dm)
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

        BndryRegister::define(lo_face,typ,0,1,0,nvar,dm);
        BndryRegister::define(hi_face,typ,0,1,0,nvar,dm);
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

void
FluxRegister::Reflux (MultiFab&       S,
                      const MultiFab& volume,
                      Real            scale,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp, 
                      const Geometry& geom,
		      const Real*     multf)
{
    BoxArray ba = grids; ba.grow(1);

    FabSetId                          fsid[2*BL_SPACEDIM];
    FabSetCopyDescriptor              fscd;
    std::deque<FluxRegister::Rec>     Recs;
    std::vector< std::pair<int,Box> > isects;

    for (OrientationIter fi; fi; ++fi)
        fsid[fi()] = fscd.RegisterFabSet(&bndry[fi()]);

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int  idx = mfi.index();
        const Box& vbx = mfi.validbox();
        //
        // Find flux register that intersects with this grid.
        //
        ba.intersections(vbx,isects);

        for (int i = 0, N = isects.size(); i < N; i++)
        {
            const int k = isects[i].first;

            for (OrientationIter fi; fi; ++fi)
            {
                //
                // low (high) face of fine grid => high (low)
                // face of the exterior coarse grid cell updated.
                //
                const Box ovlp = vbx & BoxLib::adjCell(grids[k],fi());

                if (ovlp.ok())
                {
                    FillBoxId fbid = fscd.AddBox(fsid[fi()],
                                                 bndry[fi()].box(k),
                                                 0,
                                                 k,
                                                 src_comp,
                                                 0,
                                                 num_comp);

                    Recs.push_back(Rec(idx,k,fi(),fbid));
                }
            }
        }
    }
    //
    // Add periodic possibilities.
    //
    if (geom.isAnyPeriodic())
    {
        Array<IntVect> pshifts(27);

        for (MFIter mfi(S); mfi.isValid(); ++mfi)
        {
            const int  idx  = mfi.index();
            const Box& vbx  = mfi.validbox();

            for (int k = 0, N = grids.size(); k < N; k++)
            {
                const Box& bx = ba[k];

                if (!geom.Domain().contains(bx))
                {
                    geom.periodicShift(bx,vbx,pshifts);

                    const Box& kgrid = grids[k];

                    for (int iiv = 0, M = pshifts.size(); iiv < M; iiv++)
                    {
                        const Box sftbox = vbx + pshifts[iiv];

                        BL_ASSERT(bx.intersects(sftbox));

                        for (OrientationIter fi; fi; ++fi)
                        {
                            //
                            // low (high)  face of fine grid => high (low)
                            // face of the exterior coarse grid cell updated.
                            //
                            const Box ovlp = sftbox & BoxLib::adjCell(kgrid,fi());

                            if (ovlp.ok())
                            {
                                FillBoxId fbid = fscd.AddBox(fsid[fi()],
                                                             bndry[fi()].box(k),
                                                             0,
                                                             k,
                                                             src_comp,
                                                             0,
                                                             num_comp);

                                Recs.push_back(Rec(pshifts[iiv],idx,k,fi(),fbid));
                            }
                        }
                    }
                }
            }
        }
    }

    fscd.CollectData();

    FArrayBox reg;

    for (std::deque<FluxRegister::Rec>::const_iterator it = Recs.begin(),
             End = Recs.end();
         it != End;
         ++it)
    {
        const Rec&       rf   = *it;
        const FillBoxId& fbid = rf.m_fbid;

        BL_ASSERT(bndry[rf.m_face].box(rf.m_idx) == fbid.box());
        BL_ASSERT(S.DistributionMap()[rf.m_fabidx] == ParallelDescriptor::MyProc());
        BL_ASSERT(volume.DistributionMap()[rf.m_fabidx] == ParallelDescriptor::MyProc());

	Real mult; 
	if (multf == 0)
	  mult = rf.m_face.isLow() ? -scale : scale;
	else
	  mult = (*multf)*scale;

        FArrayBox&       fab_S      = S[rf.m_fabidx];
        const FArrayBox& fab_volume = volume[rf.m_fabidx];
        Real*            s_dat      = fab_S.dataPtr(dest_comp);
        const int*       slo        = fab_S.loVect();
        const int*       shi        = fab_S.hiVect();
        const Real*      vol_dat    = fab_volume.dataPtr();
        const Box        fine_face  = BoxLib::adjCell(grids[rf.m_idx],rf.m_face);
        const Box        sftbox     = S.box(rf.m_fabidx) + rf.m_shift;
        const Box        ovlp       = sftbox & fine_face;
        const int*       lo         = ovlp.loVect();
        const int*       hi         = ovlp.hiVect();
        const int*       rlo        = fine_face.loVect();
        const int*       rhi        = fine_face.hiVect();
        const int*       shft       = rf.m_shift.getVect();
        const int*       vlo        = fab_volume.loVect();
        const int*       vhi        = fab_volume.hiVect();

        reg.resize(fbid.box(),num_comp);
        fscd.FillFab(fsid[rf.m_face],fbid,reg);

        const Real* reg_dat = reg.dataPtr(0);

        BL_ASSERT(ovlp.ok());

        FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                      vol_dat,ARLIM(vlo),ARLIM(vhi),
                      reg_dat,ARLIM(rlo),ARLIM(rhi),
                      lo,hi,shft,&num_comp,&mult);
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

    MultiFab volume(S.boxArray(), 1, S.nGrow());

    volume.setVal(D_TERM(dx[0],*dx[1],*dx[2]), 0, 1, S.nGrow());

    Reflux(S,volume,scale,src_comp,dest_comp,num_comp,geom);
}

//
// Some useful typedefs.
//
typedef FabArrayBase::CopyComTag::CopyComTagsContainer CopyComTagsContainer;

typedef FabArrayBase::CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;

void
FluxRegister::CrseInitDoit (const MultiFab& mflx,
                            const MultiFab* area,
                            int             dir,
                            int             srccomp,
                            int             destcomp,
                            int             numcomp,
                            Real            mult,
                            FrOp            op)
{
    const int MyProc = ParallelDescriptor::MyProc();

    FArrayBox                         fab;
    FabArrayBase::CopyComTag          tag;
    MapOfCopyComTagContainers         m_SndTags, m_RcvTags;
    std::map<int,int>                 m_SndVols, m_RcvVols;
    std::vector< std::pair<int,Box> > isects;
    const Orientation                 face_lo(dir,Orientation::low);
    const Orientation                 face_hi(dir,Orientation::high);

    for (int pass = 0; pass < 2; pass++)
    {
        const int face = ((pass == 0) ? face_lo : face_hi);

        FabSet&                    fabset  = bndry[face];
        const BoxArray&            ba      = mflx.boxArray();
        const DistributionMapping& srcDMap = mflx.DistributionMap();
        const DistributionMapping& dstDMap = fabset.DistributionMap();

        for (int i = 0, N = fabset.size(); i < N; i++)
        {
            ba.intersections(fabset.fabbox(i),isects);

            const int dst_owner = dstDMap[i];

            for (int j = 0, M = isects.size(); j < M; j++)
            {
                const Box& bx        = isects[j].second;
                const int  k         = isects[j].first;
                const int  src_owner = srcDMap[k];

                if (dst_owner != MyProc && src_owner != MyProc) continue;

                tag.box      = bx;
                tag.srcIndex = face;

                const int vol = bx.numPts();

                if (dst_owner == MyProc)
                {
                    tag.fabIndex = i;

                    if (src_owner == MyProc)
                    {
                        //
                        // Do the local work right here.
                        //
                        if (op == COPY)
                        {
                            fabset[i].copy(mflx[k],bx,srccomp,bx,destcomp,numcomp);
                            fabset[i].mult(mult,bx,destcomp,numcomp);
                            if (area)
                            {
                                for (int n = 0; n < numcomp; n++)
                                    fabset[i].mult((*area)[k],bx,bx,0,destcomp+n,1);
                            }
                        }
                        else
                        {
                            fab.resize(bx,numcomp);
                            fab.copy(mflx[k],bx,srccomp,bx,0,numcomp);
                            fab.mult(mult);
                            if (area)
                            {
                                for (int n = 0; n < numcomp; n++)
                                    fab.mult((*area)[k],bx,bx,0,n,1);
                            }
                            fabset[i].plus(fab,bx,bx,0,destcomp,numcomp);
                        }
                    }
                    else
                    {
                        FabArrayBase::SetRecvTag(m_RcvTags,src_owner,tag,m_RcvVols,vol);
                    }
                }
                else if (src_owner == MyProc)
                {
                    tag.fabIndex = k;

                    FabArrayBase::SetSendTag(m_SndTags,dst_owner,tag,m_SndVols,vol);
                }
            }
        }
    }

#ifdef BL_USE_MPI
    if (ParallelDescriptor::NProcs() == 1) return;
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    const int SeqNum = ParallelDescriptor::SeqNum();

    if (m_SndTags.empty() && m_RcvTags.empty())
        //
        // No parallel work for this MPI process to do.
        //
        return;
    //
    // If area is defined we'll have one component of area appended to numcomp data.
    //
    const int NumCompArea = (area ? numcomp+1 : numcomp);

    Array<MPI_Status>  stats;
    Array<int>         recv_from, index;
    Array<double*>     recv_data, send_data;
    Array<MPI_Request> recv_reqs, send_reqs;
    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    double* the_recv_data = 0;

    FabArrayBase::PostRcvs(m_RcvTags,m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,NumCompArea,SeqNum);
    //
    // Send the data.
    //
    for (MapOfCopyComTagContainers::const_iterator m_it = m_SndTags.begin(),
             m_End = m_SndTags.end();
         m_it != m_End;
         ++m_it)
    {
        std::map<int,int>::const_iterator vol_it = m_SndVols.find(m_it->first);

        BL_ASSERT(vol_it != m_SndVols.end());

        const int N = vol_it->second*NumCompArea;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        double* data = static_cast<double*>(BoxLib::The_Arena()->alloc(N*sizeof(double)));
        double* dptr = data;

        for (CopyComTagsContainer::const_iterator it = m_it->second.begin(),
                 End = m_it->second.end();
             it != End;
             ++it)
        {
            const Box& bx = it->box;
            fab.resize(bx,numcomp);
            fab.copy(mflx[it->fabIndex],bx,srccomp,bx,0,numcomp);
            const int NumPts = bx.numPts(), Cnt = NumPts*numcomp;
            memcpy(dptr,fab.dataPtr(),Cnt*sizeof(double));
            dptr += Cnt;
            if (area)
            {
                //
                // For every box of fab data append one component of area data.
                //
                fab.copy((*area)[it->fabIndex],bx,0,bx,0,1);
                memcpy(dptr,fab.dataPtr(),NumPts*sizeof(double));
                dptr += NumPts;
            }
        }
        BL_ASSERT(data+N == dptr);

        if (FabArrayBase::do_async_sends)
        {
            send_data.push_back(data);
            send_reqs.push_back(ParallelDescriptor::Asend(data,N,m_it->first,SeqNum).req());
        }
        else
        {
            ParallelDescriptor::Send(data,N,m_it->first,SeqNum);
            BoxLib::The_Arena()->free(data);
        }
    }
    //
    // Now receive and unpack FAB data as it becomes available.
    //
    const int N_rcvs = m_RcvTags.size();

    index.resize(N_rcvs);
    stats.resize(N_rcvs);

    for (int NWaits = N_rcvs, completed; NWaits > 0; NWaits -= completed)
    {
        ParallelDescriptor::Waitsome(recv_reqs, completed, index, stats);

        for (int k = 0; k < completed; k++)
        {
            const double* dptr = recv_data[index[k]];

            BL_ASSERT(dptr != 0);

            MapOfCopyComTagContainers::const_iterator m_it = m_RcvTags.find(recv_from[index[k]]);

            BL_ASSERT(m_it != m_RcvTags.end());

            for (CopyComTagsContainer::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                //
                // Area will be the "numcomp" component in our data.
                //
                const Box& bx = it->box;
                fab.resize(bx,NumCompArea);
                const int Cnt = bx.numPts()*NumCompArea;
                memcpy(fab.dataPtr(),dptr,Cnt*sizeof(double));
                fab.mult(mult,0,numcomp);
                if (area)
                {
                    for (int n = 0; n < numcomp; n++)
                        fab.mult(fab,bx,bx,numcomp,n,1);
                }

                if (op == COPY)
                {
                    bndry[it->srcIndex][it->fabIndex].copy(fab,bx,0,bx,destcomp,numcomp);
                }
                else
                {
                    bndry[it->srcIndex][it->fabIndex].plus(fab,bx,bx,0,destcomp,numcomp);
                }
                dptr += Cnt;
            }
        }
    }

    BoxLib::The_Arena()->free(the_recv_data);

    if (FabArrayBase::do_async_sends && !m_SndTags.empty())
        FabArrayBase::GrokAsyncSends(m_SndTags.size(),send_reqs,send_data,stats);
#endif /*BL_USE_MPI*/
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
    BL_ASSERT(mflx.boxArray() == area.boxArray());
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    CrseInitDoit(mflx,&area,dir,srccomp,destcomp,numcomp,mult,op);
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

    CrseInitDoit(mflx,0,dir,srccomp,destcomp,numcomp,mult,op);
}

//
// Helper function and data for CrseInit()/CrseInitFinish().
//

static Array<int>                           CIMsgs;
static std::vector<FabArrayBase::FabComTag> CITags;
static std::vector<FArrayBox*>              CIFabs;
static BArena                               CIArena;

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
      Real               mult,
      FluxRegister::FrOp op = FluxRegister::COPY)
{
    const DistributionMapping& dMap = bndry[face].DistributionMap();

    FArrayBox tmp;

    if (ParallelDescriptor::MyProc() == dMap[k])
    {
        //
        // Local data.
        //
        if (op == FluxRegister::COPY) 
        {
            bndry[face][k].copy(flux, bx, srccomp, bx, destcomp, numcomp);
            bndry[face][k].mult(mult, bx, destcomp, numcomp);    
        }
        else
        {
            tmp.resize(bx, numcomp);
            tmp.copy(flux, bx, srccomp, bx, 0, numcomp);
            tmp.mult(mult);
            bndry[face][k].plus(tmp, bx, bx, 0, destcomp, numcomp);
        }
    }
    else
    {
        FabArrayBase::FabComTag tag;

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

        if (CIMsgs.size() == 0)
            CIMsgs.resize(ParallelDescriptor::NProcs(), 0);

        CIMsgs[dMap[k]]++;
    }
}

void
FluxRegister::CrseInit (const FArrayBox& flux,
                        const Box&       subbox,
                        int              dir,
                        int              srccomp,
                        int              destcomp,
                        int              numcomp,
                        Real             mult,
                        FrOp             op)
{
    BL_ASSERT(flux.box().contains(subbox));
    BL_ASSERT(srccomp  >= 0 && srccomp+numcomp  <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    if (ParallelDescriptor::IOProcessor())
        BoxLib::Warning("\n*** FluxRegister::CrseInit(const FArrayBox&,...) is deprecated; please use CrseInit(MultiFab&,...) instead!!");
    
    const Orientation lo(dir,Orientation::low);

    std::vector< std::pair<int,Box> > isects;

    bndry[lo].boxArray().intersections(subbox,isects);

    for (int i = 0, N = isects.size(); i < N; i++)
    {
        DoIt(lo,isects[i].first,bndry,isects[i].second,flux,srccomp,destcomp,numcomp,mult,op);
    }

    const Orientation hi(dir,Orientation::high);

    bndry[hi].boxArray().intersections(subbox,isects);

    for (int i = 0, N = isects.size(); i < N; i++)
    {
        DoIt(hi,isects[i].first,bndry,isects[i].second,flux,srccomp,destcomp,numcomp,mult,op);
    }
}

void
FluxRegister::CrseInitFinish (FrOp op)
{
    if (ParallelDescriptor::IOProcessor())
        BoxLib::Warning("\n*** FluxRegister::CrseInitFinish() is deprecated; please use CrseInit(MultiFab&,...) instead!!");

    if (ParallelDescriptor::NProcs() == 1) return;

#if BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();

    const bool verbose = false;

    if (verbose)
    {
        long count = 0;

        for (int i = 0, N = CIFabs.size(); i < N; i++)
        {
            count += CIFabs[i]->box().numPts()*CIFabs[i]->nComp()*sizeof(Real);
        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceLongMax(count,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "FluxRegister::CrseInitFinish(): HWM = " << count << std::endl;
    }

    BL_ASSERT(CITags.size() == CIFabs.size());

    if (CIMsgs.size() == 0)
        CIMsgs.resize(ParallelDescriptor::NProcs(),0);

    BL_ASSERT(CIMsgs[MyProc] == 0);

    Array<int> Rcvs(NProcs,0);
    //
    // Set Rcvs[i] to # of blocks we expect to get from CPU i ...
    //
    BL_MPI_REQUIRE( MPI_Alltoall(CIMsgs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    int NumRcvs = 0;
    for (int i = 0; i < NProcs; i++)
        NumRcvs += Rcvs[i];
    if (NumRcvs == 0) NumRcvs = 1;
    Array<ParallelDescriptor::CommData> recvdata(NumRcvs);

    int NumSnds = 0;
    for (int i = 0; i < NProcs; i++)
        NumSnds += CIMsgs[i];
    if (NumSnds == 0) NumSnds = 1;
    Array<ParallelDescriptor::CommData> senddata(NumSnds);
    //
    // Make sure we can treat CommData as a stream of integers.
    //
    BL_ASSERT(sizeof(ParallelDescriptor::CommData) == ParallelDescriptor::CommData::DIM*sizeof(int));
    {
        Array<int> sendcnts(NProcs,0), sdispls(NProcs,0);
        Array<int> recvcnts(NProcs,0), rdispls(NProcs,0), offset(NProcs,0);

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i]   * ParallelDescriptor::CommData::DIM;
            sendcnts[i] = CIMsgs[i] * ParallelDescriptor::CommData::DIM;

            if (i < NProcs-1)
            {
                rdispls[i+1] = rdispls[i] + recvcnts[i];
                sdispls[i+1] = sdispls[i] + sendcnts[i];
            }
        }

        for (int i = 1; i < NProcs; i++)
            offset[i] = offset[i-1] + CIMsgs[i-1];

        for (int j = 0, N = CITags.size(); j < N; j++)
        {
            ParallelDescriptor::CommData data(CITags[j].face,
                                              CITags[j].fabIndex,
                                              MyProc,
                                              0,
                                              CITags[j].nComp,
                                              CITags[j].destComp,   // Store as srcComp()
                                              0,                    // Not used.
                                              CITags[j].box);

            senddata[offset[CITags[j].toProc]++] = data;
        }

        BL_MPI_REQUIRE( MPI_Alltoallv(senddata.dataPtr(),
                                      sendcnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      recvdata.dataPtr(),
                                      recvcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      ParallelDescriptor::Communicator()) );
    }
    Array<int> sendcnts(NProcs,0), sdispls(NProcs,0);
    Array<int> recvcnts(NProcs,0), rdispls(NProcs,0);

    int send_sz = 0, recv_sz = 0, roffset = 0, soffset = 0;

    for (int i = 0; i < NProcs; i++)
    {
        size_t recv_N = 0;
        for (int j = 0; j < Rcvs[i]; j++)
            recv_N += recvdata[roffset+j].box().numPts() * recvdata[roffset+j].nComp();
        recv_sz    += recv_N;
        recvcnts[i] = recv_N;
        roffset    += Rcvs[i];

        size_t send_N = 0;
        for (int j = 0; j < CIMsgs[i]; j++)
            send_N += senddata[soffset+j].box().numPts() * senddata[soffset+j].nComp();
        send_sz    += send_N;
        sendcnts[i] = send_N;
        soffset    += CIMsgs[i];

        if (i < NProcs-1)
        {
            rdispls[i+1] = rdispls[i] + recvcnts[i];
            sdispls[i+1] = sdispls[i] + sendcnts[i];
        }
    }

    BL_ASSERT((send_sz*sizeof(Real)) < std::numeric_limits<size_t>::max());

    Real* sendbuf = static_cast<Real*>(BoxLib::The_Arena()->alloc(send_sz*sizeof(Real)));

    Array<int> offset = sdispls;

    for (int j = 0; j < CITags.size(); j++)
    {
        BL_ASSERT(CITags[j].box == CIFabs[j]->box());
        BL_ASSERT(CITags[j].nComp == CIFabs[j]->nComp());
        const int N = CITags[j].box.numPts() * CITags[j].nComp;
        memcpy(&sendbuf[offset[CITags[j].toProc]], CIFabs[j]->dataPtr(), N * sizeof(Real));
        delete CIFabs[j];
        CIFabs[j] = 0;
        offset[CITags[j].toProc] += N;
    }

    BL_ASSERT((recv_sz*sizeof(Real)) < std::numeric_limits<size_t>::max());

    Real* recvbuf = static_cast<Real*>(BoxLib::The_Arena()->alloc(recv_sz*sizeof(Real)));

    BL_MPI_REQUIRE( MPI_Alltoallv(sendbuf,
                                  sendcnts.dataPtr(),
                                  sdispls.dataPtr(),
                                  ParallelDescriptor::Mpi_typemap<Real>::type(),
                                  recvbuf,
                                  recvcnts.dataPtr(),
                                  rdispls.dataPtr(),
                                  ParallelDescriptor::Mpi_typemap<Real>::type(),
                                  ParallelDescriptor::Communicator()) );

    BoxLib::The_Arena()->free(sendbuf);

    FArrayBox fab;

    roffset = 0;

    for (int i = 0; i < NProcs; i++)
    {
        const Real* dptr = &recvbuf[rdispls[i]];

        for (int j = 0; j < Rcvs[i]; j++)
        {
            const ParallelDescriptor::CommData& cd = recvdata[roffset+j];
            fab.resize(cd.box(),cd.nComp());
            const int N = fab.box().numPts() * fab.nComp();
            memcpy(fab.dataPtr(), dptr, N * sizeof(Real));
            if (op == COPY)
            {
                bndry[cd.face()][cd.fabindex()].copy(fab, fab.box(), 0, fab.box(), cd.srcComp(), cd.nComp());
            }
            else
            {
                bndry[cd.face()][cd.fabindex()].plus(fab, fab.box(), fab.box(), 0, cd.srcComp(), cd.nComp());
            }
            dptr += N;
        }

        roffset += Rcvs[i];
    }

    BoxLib::The_Arena()->free(recvbuf);

    CIFabs.erase(CIFabs.begin(), CIFabs.end());
    CITags.erase(CITags.begin(), CITags.end());

    for (int i = 0; i < NProcs; i++) CIMsgs[i] = 0;
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
    const int N = mflx.IndexMap().size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
    {
        const int k = mflx.IndexMap()[i];
        FineAdd(mflx[k],dir,k,srccomp,destcomp,numcomp,mult);
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
    const int N = mflx.IndexMap().size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
    {
        const int k = mflx.IndexMap()[i];
        FineAdd(mflx[k],area[k],dir,k,srccomp,destcomp,numcomp,mult);
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

    is >> ratio;
    is >> fine_level;
    is >> ncomp;

    BndryRegister* br = this;

    br->read(name,is);
}
