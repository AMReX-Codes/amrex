#include "AMReX_NonLocalBC.H"

namespace amrex {
namespace NonLocalBC {
#ifdef AMREX_USE_MPI
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

    std::size_t TotalRcvsVolume = 0;
    for (const auto& kv : cctc) // loop over senders
    {
        std::size_t nbytes = 0;
        for (auto const& cct : kv.second)
        {
            // the followng does not work since src[cct.srcIndex] will throw an assertion
            //   nbytes += src[cct.srcIndex].nBytes(cct.sbox, components.n_components);
            // Can we have a static function FAB::nBytes(Box, int) ?
            nbytes += cct.sbox.numPts() * object_size * components.n_components;
        }

        std::size_t acd = ParallelDescriptor::alignof_comm_data(nbytes);
        nbytes = amrex::aligned_size(acd, nbytes);  // so that nbytes are aligned

        // Also need to align the offset properly
        TotalRcvsVolume = amrex::aligned_size(std::max(align, acd), TotalRcvsVolume);

        comm.offset.push_back(TotalRcvsVolume);
        TotalRcvsVolume += nbytes;

        comm.data.push_back(nullptr);
        comm.size.push_back(nbytes);
        comm.rank.push_back(kv.first);
        comm.request.push_back(MPI_REQUEST_NULL);
        comm.cctc.push_back(&kv.second);
    }

    const int N_recvs = comm.data.size();
    if (TotalRcvsVolume == 0)
    {
        comm.the_data = nullptr;
    }
    else
    {
        comm.the_data.reset(static_cast<char*>(amrex::The_FA_Arena()->alloc(TotalRcvsVolume)));
        for (int i = 0; i < N_recvs; ++i) {
            comm.data[i] = comm.the_data.get() + comm.offset[i];
        }
    }
    comm.stats.resize(N_recvs);
}

void PostRecvs(CommData& recv, int mpi_tag) {
    const int n_recv = recv.data.size();
    BL_ASSERT(n_recv == recv.offset.size());
    BL_ASSERT(n_recv == recv.size.size());
    BL_ASSERT(n_recv == recv.rank.size());
    BL_ASSERT(n_recv == recv.request.size());
    MPI_Comm comm = ParallelContext::CommunicatorSub();
    char* const the_recv_data = recv.the_data.get();
    for (int i = 0; i < recv.data.size(); ++i) {
        recv.data[i] = the_recv_data + recv.offset[i];
        if (recv.size[i] > 0) {
            const int rank = ParallelContext::global_to_local_rank(recv.rank[i]);
            recv.request[i] =
                ParallelDescriptor::Arecv(recv.data[i], recv.size[i], rank, mpi_tag, comm).req();
        }
    }
}

void PostSends(CommData& send, int mpi_tag) {
    MPI_Comm comm = ParallelContext::CommunicatorSub();
    const int n_sends = send.data.size();
    BL_ASSERT(n_sends == send.size.size());
    BL_ASSERT(n_sends == send.rank.size());
    BL_ASSERT(n_sends == send.request.size());
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

template void ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox, const FabArray<FArrayBox>& src, int destcomp,
             int srccomp, int numcomp, const IntVect& ngrow, Identity, Identity);

} // namespace NonLocalBC
} // namespace amrex