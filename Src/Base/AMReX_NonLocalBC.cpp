#include "AMReX_NonLocalBC.H"

namespace amrex {
namespace NonLocalBC {
#ifdef AMREX_USE_MPI
// Note, this is copied and modified from PrepareSendBuffers and PostRcvs
void PrepareCommBuffers(CommData& comm, const PackComponents& components, 
                        const FabArrayBase::MapOfCopyComTagContainers& cctc,
                        std::size_t object_size, std::size_t align)
{
    comm.data.clear();
    comm.size.clear();
    comm.rank.clear();
    comm.request.clear();
    comm.offset.clear();
    comm.cctc.clear();

    std::size_t total_volume = 0;
    for (const auto& kv : cctc)
    {
        std::size_t nbytes = 0;
        for (auto const& cct : kv.second)
        {
            // Note: Does this hold for all FAB types? 
            // This nBytes() implementation is currently also assumed in unpack_recv_buffers
            nbytes += cct.sbox.numPts() * object_size * components.n_components;
        }

        std::size_t acd = ParallelDescriptor::alignof_comm_data(nbytes);
        nbytes = amrex::aligned_size(acd, nbytes);  // so that nbytes are aligned

        // Also need to align the offset properly
        total_volume = amrex::aligned_size(std::max(align, acd), total_volume);

        comm.offset.push_back(total_volume);
        total_volume += nbytes;

        comm.data.push_back(nullptr);
        comm.size.push_back(nbytes);
        comm.rank.push_back(kv.first);
        comm.request.push_back(MPI_REQUEST_NULL);
        comm.cctc.push_back(&kv.second);
    }

    const int N_recvs = comm.data.size();
    if (total_volume == 0)
    {
        comm.the_data = nullptr;
    }
    else
    {
        comm.the_data.reset(static_cast<char*>(amrex::The_FA_Arena()->alloc(total_volume)));
        for (int i = 0; i < N_recvs; ++i) {
            comm.data[i] = comm.the_data.get() + comm.offset[i];
        }
    }
    comm.stats.resize(N_recvs);
}

void PostRecvs(CommData& recv, int mpi_tag) {
    const int n_recv = recv.data.size();
    AMREX_ASSERT(n_recv == recv.offset.size());
    AMREX_ASSERT(n_recv == recv.size.size());
    AMREX_ASSERT(n_recv == recv.rank.size());
    AMREX_ASSERT(n_recv == recv.request.size());
    MPI_Comm comm = ParallelContext::CommunicatorSub();
    for (int i = 0; i < recv.data.size(); ++i) {
        if (recv.size[i] > 0) {
            const int rank = ParallelContext::global_to_local_rank(recv.rank[i]);
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
            send.request[j] =
                ParallelDescriptor::Asend(send.data[j], send.size[j], rank, mpi_tag, comm).req();
        }
    }
}
#endif

template void Rotate90(FabArray<FArrayBox>& mf, int scomp, int ncomp, IntVect const& nghost,
                       Box const& domain);

template void Rotate90(FabArray<FArrayBox>& mf, Box const& domain);

template void Rotate180(FabArray<FArrayBox>& mf, int scomp, int ncomp, IntVect const& nghost,
                        Box const& domain);

template void Rotate180(FabArray<FArrayBox>& mf, Box const& domain);

template void FillPolar(FabArray<FArrayBox>& mf, int scomp, int ncomp, IntVect const& nghost,
                        Box const& domain);

template void FillPolar(FabArray<FArrayBox>& mf, Box const& domain);

template MultiBlockCommMetaData ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox, const FabArray<FArrayBox>& src, int destcomp,
             int srccomp, int numcomp, const IntVect& ngrow, Identity, Identity);

template MultiBlockCommMetaData ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox, const FabArray<FArrayBox>& src, int destcomp,
             int srccomp, int numcomp, const IntVect& ngrow, MultiBlockIndexMapping, Identity);

} // namespace NonLocalBC
} // namespace amrex