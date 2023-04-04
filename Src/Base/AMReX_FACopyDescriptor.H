
#ifndef BL_FACOPYDESCRIPTOR_H_
#define BL_FACOPYDESCRIPTOR_H_
#include <AMReX_Config.H>

#include <AMReX_FabArray.H>
#include <string>

namespace amrex {

//!
//! This enum and the FabCopyDescriptor class should really be nested
//! in FabArrayCopyDescriptor (not done for portability reasons).
//!

enum FillType { FillLocally, FillRemotely, Unfillable };

class FillBoxId
{
  public:

    FillBoxId () = default;

    FillBoxId (int newid, const Box& fillbox)
        :
        m_fillBox(fillbox),
        m_fillBoxId(newid)
        {}

    [[nodiscard]] int Id () const              { return m_fillBoxId;    }
    [[nodiscard]] int FabIndex () const        { return m_fabIndex;     }
    void FabIndex (int fabindex) { m_fabIndex = fabindex; }
    [[nodiscard]] const Box& box () const      { return m_fillBox;      }

private:

    Box m_fillBox;
    int m_fillBoxId{-1};
    int m_fabIndex{-1};
};

class FabArrayId
{
public:

    explicit FabArrayId (int newid = -1)
        :
        fabArrayId(newid) {}

    [[nodiscard]] int Id () const { return fabArrayId; }

    bool operator== (const FabArrayId& rhs) const
    {
        return fabArrayId == rhs.fabArrayId;
    }

private:

    int fabArrayId;
};

template <class FAB>
struct FabCopyDescriptor
{
    FabCopyDescriptor ();

    ~FabCopyDescriptor ();

    FabCopyDescriptor (const FabCopyDescriptor&) = delete;
    FabCopyDescriptor (FabCopyDescriptor&&) = delete;
    FabCopyDescriptor& operator= (const FabCopyDescriptor&) = delete;
    FabCopyDescriptor& operator= (FabCopyDescriptor&&) = delete;

    FAB*     localFabSource;
    Box      subBox;
    int      myProc{-1};
    int      copyFromProc{-1};
    int      copyFromIndex{-1};
    int      fillBoxId{-1};
    int      srcComp{-1};
    int      destComp{-1};
    int      nComp{-1};
    FillType fillType{Unfillable};
    bool     cacheDataAllocated{false};
};

template <class FAB>
FabCopyDescriptor<FAB>::FabCopyDescriptor ()
    :
    localFabSource(nullptr)

{}

template <class FAB>
FabCopyDescriptor<FAB>::~FabCopyDescriptor ()
{
    if (cacheDataAllocated)
        delete localFabSource;
}

/**
* \brief This class orchestrates filling a destination fab of size destFabBox
* from fabarray on the local processor (myProc).
*/

template <class FAB>
class FabArrayCopyDescriptor
{
  typedef std::multimap<int,FabCopyDescriptor<FAB>*> FCDMap;

  typedef typename FCDMap::value_type     FCDMapValueType;
  typedef typename FCDMap::iterator       FCDMapIter;
  typedef typename FCDMap::const_iterator FCDMapConstIter;

  public:

    FabArrayCopyDescriptor () = default;

    ~FabArrayCopyDescriptor ();

    FabArrayCopyDescriptor (const FabArrayCopyDescriptor<FAB>&) = delete;
    FabArrayCopyDescriptor (FabArrayCopyDescriptor<FAB>&&) = delete;
    FabArrayCopyDescriptor<FAB>& operator= (const FabArrayCopyDescriptor<FAB> &) = delete;
    FabArrayCopyDescriptor<FAB>& operator= (FabArrayCopyDescriptor<FAB> &&) = delete;

    FabArrayId RegisterFabArray(FabArray<FAB> *fabarray);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes);

    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp);
    //!
    //! Add a box but only from FabArray[fabarrayindex].
    //!
    FillBoxId AddBox (FabArrayId fabarrayid,
                      const Box& destFabBox,
                      BoxList*   unfilledBoxes,
                      int        fabarrayindex,
                      int        srccomp,
                      int        destcomp,
                      int        numcomp,
                      bool       bUseValidBox = true);

    void CollectData ();

    void FillFab (FabArrayId       faid,
                  const FillBoxId& fillboxid,
                  FAB&             destFab);

    void FillFab (FabArrayId       faid,
                  const FillBoxId& fillboxid,
                  FAB&             destFab,
                  const Box&       destBox);

    void PrintStats () const;

    [[nodiscard]] bool DataAvailable () const { return dataAvailable; }

    void clear ();

    [[nodiscard]] int CurrentNFabArrays () const { return fabArrays.size(); }

    [[nodiscard]] int nFabComTags () const { return fabComTagList.size(); }

    [[nodiscard]] int nFabCopyDescs () const { return fabCopyDescList.size(); }

private:
    //!
    //! Helper function for AddBox() routines.
    //!
    void AddBoxDoIt (FabArrayId fabarrayid,
                     const Box& destFabBox,
                     BoxList*   returnedUnfilledBoxes,
                     int        faindex,
                     int        srccomp,
                     int        destcomp,
                     int        numcomp,
                     bool       bUseValidBox,
                     BoxDomain& unfilledBoxDomain);
    //!
    //! Some useful typedefs.
    //!
    typedef std::vector<FabArrayBase::FabComTag> FabComTagContainer;

    typedef std::vector<FabComTagContainer::const_iterator> FabComTagIterContainer;
    //!
    //! The data.
    //!
    std::vector<FabArray<FAB>*> fabArrays;
    std::vector<FCDMap>         fabCopyDescList;
    FabComTagContainer          fabComTagList;
    int                         nextFillBoxId{0};
    bool                        dataAvailable{false};
};

template <class FAB>
FabArrayId
FabArrayCopyDescriptor<FAB>::RegisterFabArray(FabArray<FAB>* fabarray)
{
    BL_ASSERT(fabArrays.size() == fabCopyDescList.size());

    FabArrayId result(fabArrays.size());

    fabArrays.push_back(fabarray);  /* Bump size() by one */

    fabCopyDescList.push_back(FCDMap());

    return result;
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::AddBoxDoIt (FabArrayId fabarrayid,
                                         const Box& destFabBox,
                                         BoxList*   returnedUnfilledBoxes,
                                         int        faindex,
                                         int        srccomp,
                                         int        destcomp,
                                         int        numcomp,
                                         bool       bUseValidBox,
                                         BoxDomain& unfilledBoxDomain)
{
    const int myProc = ParallelDescriptor::MyProc();

    FabArray<FAB>* fabArray = fabArrays[fabarrayid.Id()];

    BL_ASSERT(faindex >= 0 && faindex < fabArray->size());

    Box intersect = destFabBox;

    if (bUseValidBox)
    {
        intersect &= fabArray->box(faindex);
    }
    else
    {
        intersect &= fabArray->fabbox(faindex);
    }

    if (intersect.ok())
    {
        auto* fcd = new FabCopyDescriptor<FAB>;

        int remoteProc     = fabArray->DistributionMap()[faindex];
        if(remoteProc >= ParallelDescriptor::NProcs()) {
            amrex::Abort("Bad remoteProc: "
                         + std::to_string(ParallelDescriptor::MyProc())
                         + ":: _in AddBoxDoIt:  nProcs remoteProc = "
                         + std::to_string(ParallelDescriptor::NProcs())
                         + "  " + std::to_string(remoteProc) + "\n");
        }
        fcd->fillBoxId     = nextFillBoxId;
        fcd->subBox        = intersect;
        fcd->myProc        = myProc;
        fcd->copyFromProc  = remoteProc;
        fcd->copyFromIndex = faindex;
        fcd->srcComp       = srccomp;
        fcd->destComp      = destcomp;
        fcd->nComp         = numcomp;

        if (ParallelDescriptor::sameTeam(remoteProc))
        {
            //
            // Data is local.
            //
            fcd->fillType       = FillLocally;
            fcd->localFabSource = &(*fabArray)[faindex];
        }
        else
        {
            //
            // Data is remote.
            //
            FabArrayBase::FabComTag fabComTag;

            dataAvailable               = false;
            fcd->fillType               = FillRemotely;
            fcd->localFabSource         = new FAB(intersect, numcomp);
            fcd->cacheDataAllocated     = true;
            fabComTag.fabArrayId        = fabarrayid.Id();
            fabComTag.fillBoxId         = nextFillBoxId;
            fabComTag.fabIndex          = faindex;
            fabComTag.procThatNeedsData = myProc;
            fabComTag.procThatHasData   = remoteProc;
            fabComTag.box               = intersect;
            fabComTag.srcComp           = srccomp;
            fabComTag.destComp          = destcomp;
            fabComTag.nComp             = numcomp;
            //
            // Do not send the data yet.
            //
            fabComTagList.push_back(fabComTag);
        }

        fabCopyDescList[fabarrayid.Id()].insert(FCDMapValueType(fcd->fillBoxId,fcd));

        if (returnedUnfilledBoxes != nullptr)
        {
            unfilledBoxDomain.rmBox(intersect);
        }
    }
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes,
                                     int        srccomp,
                                     int        destcomp,
                                     int        numcomp)
{
    BoxDomain unfilledBoxDomain(destFabBox.ixType());

    if (returnedUnfilledBoxes != nullptr)
    {
        unfilledBoxDomain.add(destFabBox);
    }

    std::vector< std::pair<int,Box> > isects;

    fabArrays[fabarrayid.Id()]->boxArray().intersections(destFabBox,isects);

    for (auto & isect : isects)
    {
        AddBoxDoIt(fabarrayid,
                   destFabBox,
                   returnedUnfilledBoxes,
                   isect.first,
                   srccomp,
                   destcomp,
                   numcomp,
                   true,
                   unfilledBoxDomain);
    }

    if (returnedUnfilledBoxes != nullptr)
    {
        returnedUnfilledBoxes->clear();
        (*returnedUnfilledBoxes) = unfilledBoxDomain.boxList();
    }

    return FillBoxId(nextFillBoxId++, destFabBox);
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes,
                                     int        fabarrayindex,
                                     int        srccomp,
                                     int        destcomp,
                                     int        numcomp,
                                     bool       bUseValidBox)
{
    BoxDomain unfilledBoxDomain(destFabBox.ixType());

    if (returnedUnfilledBoxes != nullptr)
    {
        unfilledBoxDomain.add(destFabBox);
    }

    AddBoxDoIt(fabarrayid,
               destFabBox,
               returnedUnfilledBoxes,
               fabarrayindex,
               srccomp,
               destcomp,
               numcomp,
               bUseValidBox,
               unfilledBoxDomain);

    if (returnedUnfilledBoxes != nullptr)
    {
        returnedUnfilledBoxes->clear();
        (*returnedUnfilledBoxes) = unfilledBoxDomain.boxList();
    }

    return FillBoxId(nextFillBoxId++, destFabBox);
}

template <class FAB>
FillBoxId
FabArrayCopyDescriptor<FAB>::AddBox (FabArrayId fabarrayid,
                                     const Box& destFabBox,
                                     BoxList*   returnedUnfilledBoxes)
{
    return AddBox(fabarrayid,
                  destFabBox,
                  returnedUnfilledBoxes,
                  0,
                  0,
                  fabArrays[fabarrayid.Id()]->nComp(),
                  true);
}

template <class FAB>
FabArrayCopyDescriptor<FAB>::~FabArrayCopyDescriptor()
{
   clear();
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::clear ()
{
    for (unsigned int i = 0, N = fabCopyDescList.size(); i < N; ++i)
    {
        for (auto fmi = fabCopyDescList[i].begin(), End = fabCopyDescList[i].end();
             fmi != End;
             ++fmi)
        {
            delete (*fmi).second;
        }
    }

    fabArrays.clear();
    fabCopyDescList.clear();
    fabComTagList.clear();

    nextFillBoxId = 0;
    dataAvailable = false;
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::CollectData ()
{
    dataAvailable = true;

    if (ParallelDescriptor::NProcs() == 1) return;

#ifdef BL_USE_MPI
    typedef typename FAB::value_type value_type;

    BL_PROFILE("FabArrayCopyDescriptor::CollectData()");

    const int MyProc = ParallelDescriptor::MyProc();
    amrex::ignore_unused(MyProc);

    int Total_Rcvs_Size = 0;
    //
    // We use this to make finding matching FabComTags more efficient.
    //
    std::map< int,FabComTagIterContainer > RcvTags;

    std::map<int,int> Snds, Rcvs, Npts;
    //
    // Set Rcvs[i] to # of blocks needed from CPU i
    //
    for (auto it = fabComTagList.begin(),
             End = fabComTagList.end();
         it != End;
         ++it)
    {
        BL_ASSERT(it->box.ok());
        BL_ASSERT(it->procThatNeedsData == MyProc);
        BL_ASSERT(it->procThatHasData   != MyProc);

        const int Who = it->procThatHasData;
        const auto Cnt = static_cast<int>((it->box.numPts())*(it->nComp));

        RcvTags[Who].push_back(it);

        Total_Rcvs_Size += Cnt;

        if (Rcvs.count(Who) > 0)
        {
            Rcvs[Who] += 1;
        }
        else
        {
            Rcvs[Who] = 1;
        }

        if (Npts.count(Who) > 0)
        {
            Npts[Who] += Cnt;
        }
        else
        {
            Npts[Who] = Cnt;
        }
    }
    BL_ASSERT(Rcvs.count(MyProc) == 0);

    BL_ASSERT((Total_Rcvs_Size*sizeof(value_type))
              < static_cast<std::size_t>(std::numeric_limits<int>::max()));

    const int NProcs = ParallelDescriptor::NProcs();

    {
        Vector<int> SndsArray(NProcs,0), RcvsArray(NProcs,0);

        for (auto Rcv : Rcvs)
        {
            RcvsArray[Rcv.first] = Rcv.second;
        }

        {
            BL_PROFILE_VAR("CollectData_Alltoall()", blpvCDATA);
            BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int), ParallelDescriptor::MyProc(),
                            BLProfiler::BeforeCall());

            BL_MPI_REQUIRE( MPI_Alltoall(RcvsArray.dataPtr(),
                                         1,
                                         ParallelDescriptor::Mpi_typemap<int>::type(),
                                         SndsArray.dataPtr(),
                                         1,
                                         ParallelDescriptor::Mpi_typemap<int>::type(),
                                         ParallelDescriptor::Communicator()) );

            BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int), ParallelDescriptor::MyProc(),
                            BLProfiler::AfterCall());

            BL_PROFILE_VAR_STOP(blpvCDATA);
        }
        BL_ASSERT(SndsArray[MyProc] == 0);

        for (int i = 0; i < NProcs; i++)
            if (SndsArray[i] > 0)
                Snds[i] = SndsArray[i];
    }

    // There are two rounds of send and recv.
    // First, the data receivers need to send the data senders meta-data (e.g., boxes).
    // Then, the senders know what data to send and perform send.
    const int SeqNum_md   = ParallelDescriptor::SeqNum();
    const int SeqNum_data = ParallelDescriptor::SeqNum();

    const auto N_snds = static_cast<int>(Snds.size());
    const auto N_rcvs = static_cast<int>(Rcvs.size());

    if ( N_snds == 0 && N_rcvs == 0 ) return;

    const int Nints = 4 + 3*AMREX_SPACEDIM;  // # of ints in a meta-data

    // for meta-data
    Vector<int> md_sender, md_offset, md_icnts, md_bcnts;
    int* md_recv_data = nullptr;
    Vector<int*> md_send_data;
    Vector<MPI_Request> md_recv_reqs, md_send_reqs;

    // for data
    Vector<int> data_sender, data_offset;
    value_type* recv_data = nullptr;
    Vector<value_type*> send_data;
    Vector<MPI_Request> data_recv_reqs, data_send_reqs;

    if (N_snds > 0)
    {
        // Recv meta-data

        int N = 0;
        for (auto Snd : Snds)
        {
            md_sender.push_back(Snd.first);
            md_bcnts.push_back(Snd.second);
            int cnt = Snd.second * Nints;
            md_icnts.push_back(cnt);
            md_offset.push_back(N);
            N += cnt;
        }

        md_recv_data = static_cast<int*>(amrex::The_Arena()->alloc(N*sizeof(int)));

        for (int i = 0; i < N_snds; ++i)
        {
            md_recv_reqs.push_back(ParallelDescriptor::Arecv(&md_recv_data[md_offset[i]],
                                                             md_icnts[i], md_sender[i],
                                                             SeqNum_md).req());
        }
    }

    if (N_rcvs > 0)
    {
        // Send meta-data
        for (auto Rcv : Rcvs)
        {
            int rank = Rcv.first;
            int Nmds = Rcv.second;
            int cnt = Nmds * Nints;

            int* p = static_cast<int*>(amrex::The_Arena()->alloc(cnt*sizeof(int)));
            md_send_data.push_back(p);

            const FabComTagIterContainer& tags = RcvTags[rank];

            // initialized the data
            int * md = p;
            for (int i = 0; i < Nmds; ++i, md += Nints)
            {
                md[0] = tags[i]->fabArrayId;
                md[1] = tags[i]->fabIndex;
                md[2] = tags[i]->srcComp;
                md[3] = tags[i]->nComp;
                const int* lo = tags[i]->box.loVect();
                const int* hi = tags[i]->box.hiVect();
                const IntVect& bxtyp = tags[i]->box.type();
                const int* tp = bxtyp.getVect();
                AMREX_D_EXPR(md[4] = lo[0],
                             md[5] = lo[1],
                             md[6] = lo[2]);
                AMREX_D_EXPR(md[4+  AMREX_SPACEDIM] = hi[0],
                             md[5+  AMREX_SPACEDIM] = hi[1],
                             md[6+  AMREX_SPACEDIM] = hi[2]);
                AMREX_D_EXPR(md[4+2*AMREX_SPACEDIM] = tp[0],
                             md[5+2*AMREX_SPACEDIM] = tp[1],
                             md[6+2*AMREX_SPACEDIM] = tp[2]);
            }

            md_send_reqs.push_back(ParallelDescriptor::Asend(p,cnt,rank,SeqNum_md).req());
        }
    }

    if (N_rcvs > 0)
    {
        recv_data = static_cast<value_type*>(amrex::The_Arena()->alloc(Total_Rcvs_Size*sizeof(value_type)));

        // Post receives for data
        int Idx = 0;
        for (auto & Npt : Npts)
        {
            int Who = Npt.first;
            int Cnt = Npt.second;
            BL_ASSERT(Cnt > 0);
            BL_ASSERT(Cnt < std::numeric_limits<int>::max());
            data_sender.push_back(Who);
            data_recv_reqs.push_back(ParallelDescriptor::Arecv(&recv_data[Idx],
                                                               Cnt,Who,SeqNum_data).req());
            data_offset.push_back(Idx);
            Idx += Cnt;
        }
    }

    // Wait on meta-data and do send
    if (N_snds > 0)
    {
        int send_counter = 0;
        while (send_counter++ < N_snds)
        {
            MPI_Status status;
            int index;
            ParallelDescriptor::Waitany(md_recv_reqs, index, status);

            int rank = status.MPI_SOURCE;
            BL_ASSERT(status.MPI_TAG == SeqNum_md);
            BL_ASSERT(rank == md_sender[index]);

            const int* p = &md_recv_data[md_offset[index]];
            int numboxes = md_bcnts[index];
            Vector<int> faid(numboxes);
            Vector<int> fidx(numboxes);
            Vector<int> scomp(numboxes);
            Vector<int> ncomp(numboxes);
            Vector<int> npts(numboxes);
            Vector<Box> bxs;
            int N = 0;
            const int * md = p;
            for (int i = 0; i < numboxes; ++i, md += Nints)
            {
                faid[i] = md[0];
                fidx[i] = md[1];
                scomp[i] = md[2];
                ncomp[i] = md[3];
                bxs.push_back(Box(IntVect(&md[4]),
                                  IntVect(&md[4+AMREX_SPACEDIM]),
                                  IntVect(&md[4+AMREX_SPACEDIM*2])));
                npts[i] = static_cast<int>(bxs.back().numPts()*ncomp[i]);
                N += npts[i];
            }

            BL_ASSERT(N < std::numeric_limits<int>::max());

            auto* data = static_cast<value_type*>(amrex::The_Arena()->alloc(N*sizeof(value_type)));
            value_type* dptr = data;
            send_data.push_back(data);

            for (int i = 0; i < numboxes; ++i)
            {
                (*fabArrays[faid[i]])[fidx[i]].template copyToMem<RunOn::Host>(bxs[i],scomp[i],ncomp[i],dptr);
                dptr += npts[i];
            }

            data_send_reqs.push_back(ParallelDescriptor::Asend(data,N,rank,SeqNum_data).req());
        }

        amrex::The_Arena()->free(md_recv_data);
    }

    // Wait and upack data
    if (N_rcvs > 0)
    {
        Vector<MPI_Status> stats(N_rcvs);

        ParallelDescriptor::Waitall(md_send_reqs, stats);
        for (int i = 0; i < N_rcvs; ++i) {
            amrex::The_Arena()->free(md_send_data[i]);
        }

        ParallelDescriptor::Waitall(data_recv_reqs, stats);

        std::pair<FCDMapIter,FCDMapIter> match;
        std::map< int,FabComTagIterContainer >::const_iterator found;

        for (int k = 0; k < N_rcvs; k++)
        {
            const int         Who  = data_sender[k];
            const value_type* dptr = &recv_data[data_offset[k]];

            BL_ASSERT(dptr != nullptr);

            found = RcvTags.find(Who);

            BL_ASSERT(found != RcvTags.end());

            const FabComTagIterContainer& tags = found->second;

            for (auto it : tags)
            {
                const FabArrayBase::FabComTag& tag = *it;

                BL_ASSERT(tag.procThatHasData == Who);

                match = fabCopyDescList[tag.fabArrayId].equal_range(tag.fillBoxId);

                for (auto fmi = match.first; fmi != match.second; ++fmi)
                {
                    FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

                    BL_ASSERT(fcdp->fillBoxId == tag.fillBoxId);

                    if (fcdp->subBox == tag.box)
                    {
                        BL_ASSERT(fcdp->localFabSource->dataPtr() != nullptr);
                        BL_ASSERT(fcdp->localFabSource->box() == tag.box);
                        auto Cnt = static_cast<int>(tag.box.numPts()*tag.nComp);
                        fcdp->localFabSource->template copyFromMem<RunOn::Host>(tag.box,0,tag.nComp,dptr);
                        dptr += Cnt;
                        break;
                    }
                }
            }
        }

        amrex::The_Arena()->free(recv_data);
    }

    // Finished send
    if (N_snds > 0)
    {
        Vector<MPI_Status> stats(N_snds);
        ParallelDescriptor::Waitall(data_send_reqs, stats);

        for (int i = 0; i < N_snds; ++i) {
            amrex::The_Arena()->free(send_data[i]);
        }
    }

#endif /*BL_USE_MPI*/
}


template <class FAB>
void
FabArrayCopyDescriptor<FAB>::FillFab (FabArrayId       faid,
                                      const FillBoxId& fillboxid,
                                      FAB&             destFab)
{
    BL_ASSERT(dataAvailable);

    std::pair<FCDMapIter,FCDMapIter> match = fabCopyDescList[faid.Id()].equal_range(fillboxid.Id());

    for (auto fmi = match.first; fmi != match.second; ++fmi)
    {
        FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

        BL_ASSERT(fcdp->fillBoxId == fillboxid.Id());

        destFab.template copy<RunOn::Host>
                   (*fcdp->localFabSource,
                     fcdp->subBox,
                     fcdp->fillType == FillLocally ? fcdp->srcComp : 0,
                     fcdp->subBox,
                     fcdp->destComp,
                     fcdp->nComp);
    }
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::FillFab (FabArrayId       faid,
                                        const FillBoxId& fillboxid,
                                        FAB&             destFab,
                                        const Box&       destBox)
{
    BL_ASSERT(dataAvailable);

    FCDMapIter fmi = fabCopyDescList[faid.Id()].lower_bound(fillboxid.Id());

    BL_ASSERT(fmi != fabCopyDescList[faid.Id()].end());

    FabCopyDescriptor<FAB>* fcdp = (*fmi).second;

    BL_ASSERT(fcdp->fillBoxId == fillboxid.Id());

    BL_ASSERT(fcdp->subBox.sameSize(destBox));

    destFab.ParallelCopy(*fcdp->localFabSource,
                         fcdp->subBox,
                         fcdp->fillType == FillLocally ? fcdp->srcComp : 0,
                         destBox,
                         fcdp->destComp,
                         fcdp->nComp);

    BL_ASSERT(++fmi == fabCopyDescList[faid.Id()].upper_bound(fillboxid.Id()));
}

template <class FAB>
void
FabArrayCopyDescriptor<FAB>::PrintStats () const
{
    const int MyProc = ParallelDescriptor::MyProc();

    amrex::AllPrint() << "----- "
         << MyProc
         << ":  Parallel stats for FabArrayCopyDescriptor:" << '\n';

    for (int fa = 0; fa < fabArrays.size(); ++fa)
    {
      amrex::AllPrint() << "fabArrays["
             << fa
             << "]->boxArray() = "
             << fabArrays[fa]->boxArray()
             << '\n';
    }
}

}

#endif
