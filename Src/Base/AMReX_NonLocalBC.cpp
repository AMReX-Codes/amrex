#include "AMReX_NonLocalBC.H"

namespace amrex {
namespace NonLocalBC {
#ifdef AMREX_USE_MPI
// Note, this is copied and modified from PrepareSendBuffers and PostRcvs
void PrepareCommBuffers(CommData& comm,
                        const FabArrayBase::MapOfCopyComTagContainers& cctc, int n_components,
                        std::size_t object_size, std::size_t align)
{
    comm.data.clear();
    comm.size.clear();
    comm.rank.clear();
    comm.request.clear();
    comm.offset.clear();
    comm.cctc.clear();
    comm.stats.clear();

    const int N_comms = cctc.size();
    if (N_comms == 0) return;
    // reserve for upcominf push_backs
    comm.data.reserve(N_comms);
    comm.size.reserve(N_comms);
    comm.rank.reserve(N_comms);
    comm.request.reserve(N_comms);
    comm.offset.reserve(N_comms);
    comm.cctc.reserve(N_comms);
    // resize to provide buffer for later
    comm.stats.resize(N_comms);

    std::size_t total_volume = 0;
    for (const auto& kv : cctc)
    {
        std::size_t nbytes = 0;
        for (auto const& cct : kv.second)
        {
            // Note: Does this hold for all FAB types?
            // This nBytes() implementation is currently also assumed in unpack_recv_buffers
            nbytes += cct.sbox.numPts() * object_size * n_components;
        }

        std::size_t acd = ParallelDescriptor::alignof_comm_data(nbytes);
        nbytes = amrex::aligned_size(acd, nbytes);  // so that nbytes are aligned

        // Also need to align the offset properly
        total_volume = amrex::aligned_size(std::max(align, acd), total_volume);

        comm.offset.push_back(total_volume);
        comm.data.push_back(nullptr);
        comm.size.push_back(nbytes);
        comm.rank.push_back(kv.first);
        comm.request.push_back(MPI_REQUEST_NULL);
        comm.cctc.push_back(&kv.second);

        total_volume += nbytes;
    }

    if (total_volume == 0)
    {
        comm.the_data.reset();
    }
    else
    {
        comm.the_data.reset(static_cast<char*>(amrex::The_FA_Arena()->alloc(total_volume)));
        for (int i = 0; i < N_comms; ++i) {
            comm.data[i] = comm.the_data.get() + comm.offset[i];
        }
    }
}

void PostRecvs(CommData& recv, int mpi_tag) {
    AMREX_ASSERT(recv.data.size() == recv.offset.size());
    AMREX_ASSERT(recv.data.size() == recv.size.size());
    AMREX_ASSERT(recv.data.size() == recv.rank.size());
    AMREX_ASSERT(recv.data.size() == recv.request.size());
    MPI_Comm comm = ParallelContext::CommunicatorSub();
    for (int i = 0; i < recv.data.size(); ++i) {
        if (recv.size[i] > 0) {
            const int rank = ParallelContext::global_to_local_rank(recv.rank[i]);
            AMREX_ASSERT(recv.data[i] != nullptr);
            recv.request[i] =
                ParallelDescriptor::Arecv(recv.data[i], recv.size[i], rank, mpi_tag, comm).req();
        }
    }
}

void PostSends(CommData& send, int mpi_tag) {
    const int n_sends = send.data.size();
    AMREX_ASSERT(n_sends == send.size.size());
    AMREX_ASSERT(n_sends == send.rank.size());
    AMREX_ASSERT(n_sends == send.request.size());
    MPI_Comm comm = ParallelContext::CommunicatorSub();
    for (int j = 0; j < n_sends; ++j) {
        if (send.size[j] > 0) {
            const int rank = ParallelContext::global_to_local_rank(send.rank[j]);
            AMREX_ASSERT(send.data[j] != nullptr);
            send.request[j] =
                ParallelDescriptor::Asend(send.data[j], send.size[j], rank, mpi_tag, comm).req();
        }
    }
}
#endif

template MultiBlockCommMetaData ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox,
                                             const FabArray<FArrayBox>& src, int destcomp,
                                             int srccomp, int numcomp, const IntVect& ngrow,
                                             MultiBlockIndexMapping, Identity);

} // namespace NonLocalBC
} // namespace amrex
