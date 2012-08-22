
#include <winstd.H>
#include <map>

#include <FabSet.H>
#include <ParallelDescriptor.H>
#include <VisMF.H>

FabSetCopyDescriptor::FabSetCopyDescriptor ()
    :
    MultiFabCopyDescriptor() {}

FabSetCopyDescriptor::~FabSetCopyDescriptor () {}

FabSet::FabSet () {}

FabSet::~FabSet () {}

FabSet::FabSet (const BoxArray& grids, int ncomp)
    :
    MultiFab(grids,ncomp,0,Fab_allocate)
{}

void
FabSet::define (const BoxArray& grids, int ncomp)
{
    MultiFab* tmp = this;

    tmp->define(grids, ncomp, 0, Fab_allocate);
}

void
FabSet::define (const BoxArray&            grids,
                int                        ncomp,
                const DistributionMapping& dm)
{
    MultiFab* tmp = this;

    tmp->define(grids, ncomp, 0, dm, Fab_allocate);
}

const FabSet&
FabSet::copyTo (FArrayBox& dest) const
{
    copy(dest);
    return *this;
}

const FabSet&
FabSet::copyTo (FArrayBox& dest,
                int        src_comp,
                int        dest_comp,
                int        num_comp) const
{
    copy(dest,src_comp,dest_comp,num_comp);
    return *this;
}

const FabSet&
FabSet::copyTo (FArrayBox& dest,
                const Box& subbox,
                int        src_comp,
                int        dest_comp,
                int        num_comp) const
{
    copy(dest,subbox,src_comp,dest_comp,num_comp);
    return *this;
}

void
FabSet::copyTo (MultiFab& dest) const
{
    dest.copy(*this);
}

FabSet&
FabSet::copyFrom (const FabSet& src)
{
    copy(src);
    return *this;
}

FabSet&
FabSet::copyFrom (const FabSet& src,
                  int           src_comp,
                  int           dest_comp,
                  int           num_comp)
{
    copy(src,src_comp,dest_comp,num_comp);
    return *this;
}

//
// The following are different from MultiFab only in the return value
//

FabSet&
FabSet::plus (Real v,
              int  comp,
              int  num_comp)
{
    MultiFab* tmp = this;
    tmp->plus(v, comp, num_comp);
    return *this;
}

FabSet&
FabSet::plus (Real       v,
              const Box& subreg,
              int        comp,
              int        num_comp)
{
    MultiFab* tmp = this;
    tmp->plus(v, subreg, comp, num_comp);
    return *this;
}

FabSet&
FabSet::mult (Real v,
              int  comp,
              int  num_comp)
{
    MultiFab* tmp = this;
    tmp->mult(v, comp, num_comp);
    return *this;
}

FabSet&
FabSet::mult (Real       v,
              const Box& subreg,
              int        comp,
              int        num_comp)
{
    MultiFab* tmp = this;
    tmp->mult(v, subreg, comp, num_comp);
    return *this;
}


FabSet&
FabSet::copyFrom (const FArrayBox& src)
{
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        get(fsi).copy(src);
    }
    return *this;
}

FabSet&
FabSet::copyFrom (const FArrayBox& src,
                  int              src_comp,
                  int              dest_comp,
                  int              num_comp)
{
    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        get(fsi).copy(src,src_comp,dest_comp,num_comp);
    }
    return *this;
}

FabSet&
FabSet::copyFrom (const FArrayBox& src,
                  const Box&       subbox,
                  int              src_comp,
                  int              dest_comp,
                  int              num_comp)
{
    BL_ASSERT(src.box().contains(subbox));

    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        Box dbox = get(fsi).box() & subbox;

        if (dbox.ok())
        {
            get(fsi).copy(src,dbox,src_comp,dbox,dest_comp,num_comp);
        }
    }

    return *this;
}

void
FabSet::DoIt (const MultiFab& src,
              int             ngrow,
              int             scomp,
              int             dcomp,
              int             ncomp,
              How             how)
{
    BL_ASSERT(ngrow <= src.nGrow());
    BL_ASSERT((dcomp+ncomp) <= nComp());
    BL_ASSERT((scomp+ncomp) <= src.nComp());
    BL_ASSERT(how == FabSet::COPYFROM || how == FabSet::PLUSFROM);
    //
    // Some shorthand.
    //
    typedef FabArrayBase::CopyComTag::CopyComTagsContainer CopyComTagsContainer;

    typedef std::map<int,CopyComTagsContainer> MapOfCopyComTagContainers;

    const int                  MyProc  = ParallelDescriptor::MyProc();
    const DistributionMapping& srcDMap = src.DistributionMap();
    const DistributionMapping& dstDMap = DistributionMap();

    BoxArray ba_src = src.boxArray();

    ba_src.grow(ngrow);

    FArrayBox                         fab;
    FabArrayBase::CopyComTag          tag;
    MapOfCopyComTagContainers         m_SndTags, m_RcvTags;
    std::map<int,int>                 m_SndVols, m_RcvVols;
    std::vector< std::pair<int,Box> > isects;

    for (int i = 0, N = size(); i < N; i++)
    {
        ba_src.intersections(this->fabbox(i),isects);

        const int dst_owner = dstDMap[i];

        for (int j = 0, M = isects.size(); j < M; j++)
        {
            const Box& bx        = isects[j].second;
            const int  k         = isects[j].first;
            const int  src_owner = srcDMap[k];

            if (dst_owner != MyProc && src_owner != MyProc) continue;

            tag.box = bx;

            const int vol = bx.numPts();

            if (dst_owner == MyProc)
            {
                tag.fabIndex = i;

                if (src_owner == MyProc)
                {
                    //
                    // Do the local work right here.
                    //
                    if (how == COPYFROM)
                    {
                        (*this)[i].copy(src[k],bx,scomp,bx,dcomp,ncomp);
                    }
                    else
                    {
                        fab.resize(bx,ncomp);
                        fab.copy(src[k],bx,scomp,bx,0,ncomp);
                        (*this)[i].plus(fab,bx,0,dcomp,ncomp);
                    }
                }
                else
                {
                    FabArrayBase::CopyComTag::SetRecvTag(m_RcvTags,src_owner,tag,m_RcvVols,vol);
                }
            }
            else if (src_owner == MyProc)
            {
                tag.fabIndex = k;

                FabArrayBase::CopyComTag::SetSendTag(m_SndTags,dst_owner,tag,m_SndVols,vol);
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

    Array<MPI_Status>  stats;
    Array<int>         recv_from, index;
    Array<double*>     recv_data, send_data;
    Array<MPI_Request> recv_reqs, send_reqs;
    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    double* the_recv_data = 0;

    FabArrayBase::CopyComTag::PostRcvs(m_RcvTags,m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);
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

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        double* data = static_cast<double*>(BoxLib::The_Arena()->alloc(N*sizeof(double)));
        double* dptr = data;

        for (CopyComTagsContainer::const_iterator it = m_it->second.begin(),
                 End = m_it->second.end();
             it != End;
             ++it)
        {
            const Box& bx = it->box;
            fab.resize(bx, ncomp);
            fab.copy(src[it->fabIndex],bx,scomp,bx,0,ncomp);
            const int Cnt = bx.numPts()*ncomp;
            memcpy(dptr,fab.dataPtr(),Cnt*sizeof(double));
            dptr += Cnt;
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
                const Box& bx = it->box;
                fab.resize(bx,ncomp);
                const int Cnt = bx.numPts()*ncomp;
                memcpy(fab.dataPtr(),dptr,Cnt*sizeof(double));

                if (how == COPYFROM)
                {
                    (*this)[it->fabIndex].copy(fab,bx,0,bx,dcomp,ncomp);
                }
                else
                {
                    (*this)[it->fabIndex].plus(fab,bx,0,dcomp,ncomp);
                }

                dptr += Cnt;
            }
        }
    }

    BoxLib::The_Arena()->free(the_recv_data);

    if (FabArrayBase::do_async_sends && !m_SndTags.empty())
        FabArrayBase::CopyComTag::GrokAsyncSends(m_SndTags,send_reqs,send_data,stats);
#endif /*BL_USE_MPI*/
}

FabSet&
FabSet::copyFrom (const MultiFab& src,
                  int             ngrow,
                  int             scomp,
                  int             dcomp,
                  int             ncomp)
{
    DoIt(src,ngrow,scomp,dcomp,ncomp,FabSet::COPYFROM);

    return *this;
}

FabSet&
FabSet::plusFrom (const MultiFab& src,
                  int             ngrow,
                  int             scomp,
                  int             dcomp,
                  int             ncomp)
{
    DoIt(src,ngrow,scomp,dcomp,ncomp,FabSet::PLUSFROM);

    return *this;
}

//
// Linear combination this := a*this + b*src
// Note: corresponding fabsets must be commensurate.
//
FabSet&
FabSet::linComb (Real          a,
                 Real          b,
                 const FabSet& src,
                 int           scomp,
                 int           dcomp,
                 int           ncomp)
{
    BL_ASSERT(size() == src.size());

    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        BL_ASSERT(get(fsi).box() == src[fsi].box());
        //
        // WARNING: same fab used as src and dest here.
        //
        get(fsi).linComb(get(fsi),
                      get(fsi).box(),
                      dcomp,
                      src[fsi],
                      src[fsi].box(),
                      scomp,
                      a,
                      b,
                      get(fsi).box(),
                      dcomp,
                      ncomp);
    }
    return *this;
}

FabSet&
FabSet::linComb (Real            a,
                 const MultiFab& mfa,
                 int             a_comp,
                 Real            b,
                 const MultiFab& mfb,
                 int             b_comp,
                 int             dcomp,
                 int             ncomp,
                 int             ngrow)
{
    BL_ASSERT(ngrow <= mfa.nGrow());
    BL_ASSERT(ngrow <= mfb.nGrow());

    const BoxArray& bxa = mfa.boxArray();

    BL_ASSERT(bxa == mfb.boxArray());

    MultiFabCopyDescriptor mfcd;

    MultiFabId mfid_mfa = mfcd.RegisterFabArray(const_cast<MultiFab*>(&mfa));
    MultiFabId mfid_mfb = mfcd.RegisterFabArray(const_cast<MultiFab*>(&mfb));

    std::vector<FillBoxId> fbids_mfa, fbids_mfb;

    fbids_mfa.reserve(16);
    fbids_mfb.reserve(16);

    BoxArray ba_isects = bxa;
    ba_isects.grow(ngrow);

    std::vector< std::pair<int,Box> > isects;

    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        ba_isects.intersections(get(fsi).box(),isects);

        const int index = fsi.index();

        for (int j = 0, N = isects.size(); j < N; j++)
        {
            const int  grd  = isects[j].first;
            const Box& ovlp = isects[j].second;

            fbids_mfa.push_back(mfcd.AddBox(mfid_mfa,
                                            ovlp,
                                            0,
                                            grd,
                                            a_comp,
                                            0,
                                            ncomp,
                                            false));

            BL_ASSERT(fbids_mfa.back().box() == ovlp);
            //
            // Also save the index of the FAB in the FabSet.
            //
            fbids_mfa.back().FabIndex(index);

            fbids_mfb.push_back(mfcd.AddBox(mfid_mfb,
                                            ovlp,
                                            0,
                                            grd,
                                            b_comp,
                                            0,
                                            ncomp,
                                            false));

            BL_ASSERT(fbids_mfb.back().box() == ovlp);
        }
    }

    mfcd.CollectData();

    FArrayBox a_fab, b_fab;

    BL_ASSERT(fbids_mfa.size() == fbids_mfb.size());

    for (int i = 0, N = fbids_mfa.size(); i < N; i++)
    {
        a_fab.resize(fbids_mfa[i].box(), ncomp);
        b_fab.resize(fbids_mfb[i].box(), ncomp);

        mfcd.FillFab(mfid_mfa, fbids_mfa[i], a_fab);
        mfcd.FillFab(mfid_mfb, fbids_mfb[i], b_fab);

        BL_ASSERT(DistributionMap()[fbids_mfa[i].FabIndex()] == ParallelDescriptor::MyProc());

        (*this)[fbids_mfa[i].FabIndex()].linComb(a_fab,
                                                 fbids_mfa[i].box(),
                                                 0,
                                                 b_fab,
                                                 fbids_mfa[i].box(),
                                                 0,
                                                 a,
                                                 b,
                                                 fbids_mfa[i].box(),
                                                 dcomp,
                                                 ncomp);
    }

    return *this;
}

void
FabSet::write(const std::string& name) const
{
    VisMF::Write(*this,name);
}

void
FabSet::read(const std::string& name)
{
    VisMF::Read(*this,name);
}
