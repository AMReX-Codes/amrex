
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
    BL_ASSERT((dcomp+ncomp) <= nComp());
    BL_ASSERT((scomp+ncomp) <= src.nComp());

    BL_ASSERT(how == FabSet::COPYFROM || how == FabSet::PLUSFROM);

    std::vector<Box> boxes;
    std::vector<int> mfidx;
    std::vector<int> fsidx;

    boxes.reserve(32);
    mfidx.reserve(32);
    fsidx.reserve(32);
    //
    // Calculate and cache intersection info.
    //
    BoxArray ba_src = src.boxArray();

    ba_src.grow(ngrow);

    std::vector< std::pair<int,Box> > isects;

    for (FabSetIter fsi(*this); fsi.isValid(); ++fsi)
    {
        ba_src.intersections((*this)[fsi].box(),isects);

        for (int j = 0, N = isects.size(); j < N; j++)
        {
            boxes.push_back(isects[j].second);
            //
            // Maintain parallel array of indices into MultiFab.
            //
            mfidx.push_back(isects[j].first);
            //
            // Maintain parallel array of indices into FabSet.
            //
            fsidx.push_back(fsi.index());
        }
    }

    FabSetCopyDescriptor fscd;

    MultiFabId mfid = fscd.RegisterFabArray(const_cast<MultiFab*>(&src));

    BL_ASSERT(boxes.size() == mfidx.size());
    BL_ASSERT(boxes.size() == fsidx.size());

    std::vector<FillBoxId> fbids;

    fbids.reserve(boxes.size());

    const int DCOMP = (how == COPYFROM) ? dcomp : 0;

    for (int i = 0, N = boxes.size(); i < N; i++)
    {
        fbids.push_back(fscd.AddBox(mfid,
                                    boxes[i],
                                    0,
                                    mfidx[i],
                                    scomp,
                                    DCOMP,
                                    ncomp,
                                    false));

        BL_ASSERT(fbids.back().box() == boxes[i]);
        //
        // Also save the index of our FAB needing filling.
        //
        fbids.back().FabIndex(fsidx[i]);
    }

    fscd.CollectData();

    FArrayBox tmp;

    for (std::vector<FillBoxId>::const_iterator it = fbids.begin(), End = fbids.end();
         it != End;
         ++it)
    {
        const FillBoxId& fbid = *it;

        BL_ASSERT(DistributionMap()[fbid.FabIndex()] == ParallelDescriptor::MyProc());

        if (how == COPYFROM)
        {
            fscd.FillFab(mfid,fbid,(*this)[fbid.FabIndex()]);
        }
        else
        {
            tmp.resize(fbid.box(), ncomp);

            fscd.FillFab(mfid, fbid, tmp);

            (*this)[fbid.FabIndex()].plus(tmp,tmp.box(),0,dcomp,ncomp);
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
