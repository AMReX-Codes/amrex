//
// $Id: FabSet.cpp,v 1.43 2002-12-11 17:05:25 lijewski Exp $
//
#include <winstd.H>

#include <list>

#include <FabSet.H>
#include <ParallelDescriptor.H>

FabSetIter::FabSetIter (const FabSet& fabset)
    :
    MFIter(fabset)
{}

FabSetIter::~FabSetIter () {}

FabSetCopyDescriptor::FabSetCopyDescriptor ()
    :
    MultiFabCopyDescriptor() {}

FabSetCopyDescriptor::~FabSetCopyDescriptor () {}

FabSetId
FabSetCopyDescriptor::RegisterFabSet (FabSet* fabset)
{
    return RegisterMultiFab(fabset);
}

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

//
// Used in caching CollectData() stuff for copyFrom() and plusFrom().
//
struct FSRec
{
    FSRec ();

    FSRec (const BoxArray& src,
           const BoxArray& dst,
           int             ngrow);

    FSRec (const FSRec& rhs);

    ~FSRec ();

    bool operator== (const FSRec& rhs) const;
    bool operator!= (const FSRec& rhs) const;

    std::vector<Box> m_box;
    std::vector<int> m_mfidx;
    std::vector<int> m_fsidx;
    Array<int>       m_snds;
    CommDataCache    m_commdata;
    BoxArray         m_src;
    BoxArray         m_dst;
    int              m_ngrow;
};

FSRec::FSRec ()
    :
    m_ngrow(-1)
{}

FSRec::FSRec (const BoxArray& src,
              const BoxArray& dst,
              int             ngrow)
    :
    m_src(src),
    m_dst(dst),
    m_ngrow(ngrow)
{
    BL_ASSERT(ngrow >= 0);
}

FSRec::FSRec (const FSRec& rhs)
    :
    m_box(rhs.m_box),
    m_mfidx(rhs.m_mfidx),
    m_fsidx(rhs.m_fsidx),
    m_snds(rhs.m_snds),
    m_commdata(rhs.m_commdata),
    m_src(rhs.m_src),
    m_dst(rhs.m_dst),
    m_ngrow(rhs.m_ngrow)
{}

FSRec::~FSRec () {}

bool
FSRec::operator== (const FSRec& rhs) const
{
    return
        m_ngrow == rhs.m_ngrow &&
        m_src   == rhs.m_src   &&
        m_dst   == rhs.m_dst;
}

bool
FSRec::operator!= (const FSRec& rhs) const
{
    return !operator==(rhs);
}

//
// A useful typedef.
//
typedef std::list<FSRec> FSRecList;

//
// Cache of FSRec info.
//
static FSRecList TheCache;

void
FabSet::FlushCache ()
{
    TheCache.clear();
}

static
FSRec&
TheFSRec (const MultiFab& src,
          const FabSet&   dst,
          int             ngrow,
          int             scomp,
          int             ncomp)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(scomp >= 0);
    BL_ASSERT(ncomp >  0);

    const FSRec rec(src.boxArray(),dst.boxArray(),ngrow);

    for (FSRecList::iterator it = TheCache.begin(); it != TheCache.end(); ++it)
    {
        if (*it == rec)
        {
            //
            // Adjust the ncomp & scomp in CommData.
            //
            Array<CommData>& cd = (*it).m_commdata.theCommData();

            for (int i = 0; i < cd.size(); i++)
            {
                cd[i].nComp(ncomp);
                cd[i].srcComp(scomp);
            }

            return *it;
        }
    }

    TheCache.push_front(rec);
    //
    // Calculate and cache intersection info.
    //
    for (FabSetIter fsi(dst); fsi.isValid(); ++fsi)
    {
        for (int i = 0; i < src.size(); i++)
        {
            Box ovlp = dst[fsi].box() & BoxLib::grow(src.boxArray()[i],ngrow);

            if (ovlp.ok())
            {
                TheCache.front().m_box.push_back(ovlp);
                //
                // Maintain parallel array of indices into MultiFab.
                //
                TheCache.front().m_mfidx.push_back(i);
                //
                // Maintain parallel array of indices into FabSet.
                //
                TheCache.front().m_fsidx.push_back(fsi.index());
            }
        }
    }

    return TheCache.front();
}

void
FabSet::DoIt (const MultiFab& src,
              int             ngrow,
              int             scomp,
              int             dcomp,
              int             ncomp,
              How             how)
{
    BL_ASSERT((dcomp+ncomp) <= nComp());
    BL_ASSERT((scomp+ncomp) <= src.nComp());

    BL_ASSERT(how == FabSet::COPYFROM || how == FabSet::PLUSFROM);

    FArrayBox              tmp;
    FabSetCopyDescriptor   fscd;
    std::vector<FillBoxId> fbids;

    const int  MyProc = ParallelDescriptor::MyProc();
    FSRec&     fsrec  = TheFSRec(src,*this,ngrow,scomp,ncomp);
    MultiFabId mfid   = fscd.RegisterFabArray(const_cast<MultiFab*>(&src));

    BL_ASSERT(fsrec.m_box.size() == fsrec.m_mfidx.size());
    BL_ASSERT(fsrec.m_box.size() == fsrec.m_fsidx.size());

    for (int i = 0; i < fsrec.m_box.size(); i++)
    {
        fbids.push_back(fscd.AddBox(mfid,
                                    fsrec.m_box[i],
                                    0,
                                    fsrec.m_mfidx[i],
                                    scomp,
                                    how == COPYFROM ? dcomp : 0,
                                    ncomp,
                                    false));

        BL_ASSERT(fbids.back().box() == fsrec.m_box[i]);
        //
        // Also save the index of our FAB needing filling.
        //
        fbids.back().FabIndex(fsrec.m_fsidx[i]);
    }

    fscd.CollectData(&fsrec.m_snds, &fsrec.m_commdata);

    for (int i = 0; i < fbids.size(); i++)
    {
        BL_ASSERT(DistributionMap()[fbids[i].FabIndex()] == MyProc);

        if (how == COPYFROM)
        {
            fscd.FillFab(mfid,fbids[i],(*this)[fbids[i].FabIndex()]);
        }
        else
        {
            tmp.resize(fbids[i].box(), ncomp);

            fscd.FillFab(mfid, fbids[i], tmp);

            (*this)[fbids[i].FabIndex()].plus(tmp,tmp.box(),0,dcomp,ncomp);
        }
    }
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
    const BoxArray& bxb = mfb.boxArray();

    BL_ASSERT(bxa == bxb);

    MultiFabCopyDescriptor mfcd;

    MultiFabId mfid_mfa = mfcd.RegisterFabArray(const_cast<MultiFab*>(&mfa));
    MultiFabId mfid_mfb = mfcd.RegisterFabArray(const_cast<MultiFab*>(&mfb));

    std::vector<FillBoxId> fbids_mfa, fbids_mfb;

    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        for (int grd = 0; grd < bxa.size(); grd++)
        {
            Box ovlp = get(fsi).box() & BoxLib::grow(bxa[grd],ngrow);

            if (ovlp.ok())
            {
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
                fbids_mfa.back().FabIndex(fsi.index());

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
    }

    mfcd.CollectData();

    FArrayBox a_fab, b_fab;

    const int MyProc = ParallelDescriptor::MyProc();

    BL_ASSERT(fbids_mfa.size() == fbids_mfb.size());

    for (int i = 0; i < fbids_mfa.size(); i++)
    {
        a_fab.resize(fbids_mfa[i].box(), ncomp);
        b_fab.resize(fbids_mfb[i].box(), ncomp);

        mfcd.FillFab(mfid_mfa, fbids_mfa[i], a_fab);
        mfcd.FillFab(mfid_mfb, fbids_mfb[i], b_fab);

        BL_ASSERT(DistributionMap()[fbids_mfa[i].FabIndex()] == MyProc);

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
