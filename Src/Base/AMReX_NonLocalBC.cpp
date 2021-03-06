#include "AMReX_NonLocalBC.H"

namespace amrex {
namespace NonLocalBC {
#ifdef AMREX_USE_MPI
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