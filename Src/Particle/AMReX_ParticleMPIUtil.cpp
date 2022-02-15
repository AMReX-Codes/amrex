#include <AMReX_ParticleMPIUtil.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_BLProfiler.H>

namespace amrex {

#ifdef AMREX_USE_MPI

    Long CountSnds(const std::map<int, Vector<char> >& not_ours, Vector<Long>& Snds)
    {
        Long NumSnds = 0;
        for (const auto& kv : not_ours)
        {
            NumSnds       += kv.second.size();
            Snds[kv.first] = kv.second.size();
        }

        ParallelAllReduce::Max(NumSnds, ParallelContext::CommunicatorSub());

        return NumSnds;
    }

    Long doHandShake(const std::map<int, Vector<char> >& not_ours,
                     Vector<Long>& Snds, Vector<Long>& Rcvs)
    {
        Long NumSnds = CountSnds(not_ours, Snds);
        if (NumSnds == 0) return NumSnds;

        BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(Long),
                        ParallelContext::MyProcSub(), BLProfiler::BeforeCall());

        BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                     1,
                                     ParallelDescriptor::Mpi_typemap<Long>::type(),
                                     Rcvs.dataPtr(),
                                     1,
                                     ParallelDescriptor::Mpi_typemap<Long>::type(),
                                     ParallelContext::CommunicatorSub()) );

        AMREX_ASSERT(Rcvs[ParallelContext::MyProcSub()] == 0);

        BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(Long),
                        ParallelContext::MyProcSub(), BLProfiler::AfterCall());

        return NumSnds;
    }

    Long doHandShakeLocal(const std::map<int, Vector<char> >& not_ours,
                          const Vector<int>& neighbor_procs, Vector<Long>& Snds, Vector<Long>& Rcvs)
    {

        Long NumSnds = 0;
        for (const auto& kv : not_ours)
        {
            NumSnds       += kv.second.size();
            Snds[kv.first] = kv.second.size();
        }

        const int SeqNum = ParallelDescriptor::SeqNum();

        const int num_rcvs = neighbor_procs.size();
        Vector<MPI_Status>  stats(num_rcvs);
        Vector<MPI_Request> rreqs(num_rcvs);

        // Post receives
        for (int i = 0; i < num_rcvs; ++i) {
            const int Who = neighbor_procs[i];
            const Long Cnt = 1;

            AMREX_ASSERT(Who >= 0 && Who < ParallelContext::NProcsSub());

            rreqs[i] = ParallelDescriptor::Arecv(&Rcvs[Who], Cnt, Who, SeqNum,
                                                 ParallelContext::CommunicatorSub()).req();
        }

        // Send.
        for (int i = 0; i < num_rcvs; ++i) {
            const int Who = neighbor_procs[i];
            const Long Cnt = 1;

            AMREX_ASSERT(Who >= 0 && Who < ParallelContext::NProcsSub());

            ParallelDescriptor::Send(&Snds[Who], Cnt, Who, SeqNum,
                                     ParallelContext::CommunicatorSub());
        }

        if (num_rcvs > 0) {
            ParallelDescriptor::Waitall(rreqs, stats);
        }

        return NumSnds;
    }
#endif  // AMREX_USE_MPI

}
