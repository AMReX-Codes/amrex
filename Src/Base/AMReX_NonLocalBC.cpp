#include "AMReX_NonLocalBC.H"

#include <limits>

/// \cond DOXYGEN_IGNORE
namespace amrex::NonLocalBC::detail {
void split_boxes (BoxList& bl, Box const& domain)
{
    BoxList bltmp(bl.ixType());
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (auto& b : bl) {
            if (b.smallEnd(idim) < domain.smallEnd(idim) &&
                b.bigEnd(idim) >= domain.smallEnd(idim))
            {
                auto btmp = b;
                b.setBig(idim,domain.smallEnd(idim)-1);
                btmp.setSmall(idim,domain.smallEnd(idim));
                bltmp.push_back(btmp);
            }
        }
        bl.join(bltmp);
        bltmp.clear();

        for (auto& b : bl) {
            if (b.smallEnd(idim) <= domain.bigEnd(idim) &&
                b.bigEnd(idim) > domain.bigEnd(idim))
            {
                auto btmp = b;
                b.setBig(idim,domain.bigEnd(idim));
                btmp.setSmall(idim,domain.bigEnd(idim)+1);
                bltmp.push_back(btmp);
            }
        }
        bl.join(bltmp);
        bltmp.clear();
    }
}
}
/// \endcond

namespace amrex::NonLocalBC {

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

    comm.id = FabArrayBase::getNextCommMetaDataId();

    const auto N_comms = static_cast<int>(cctc.size());
    if (N_comms == 0) {
        comm.total_bytes = 0;
        comm.the_data.reset();
        return;
    }
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

        std::size_t acd = ParallelDescriptor::sizeof_selected_comm_data_type(nbytes);
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

    comm.total_bytes = total_volume;
    comm.the_data.reset();
    AcquireCommBuffer(comm);
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
    const auto n_sends = send.data.size();
    AMREX_ASSERT(n_sends == send.size.size());
    AMREX_ASSERT(n_sends == send.rank.size());
    AMREX_ASSERT(n_sends == send.request.size());
    MPI_Comm comm = ParallelContext::CommunicatorSub();
    for (int j = 0; j < static_cast<int>(n_sends); ++j) {
        if (send.size[j] > 0) {
            const int rank = ParallelContext::global_to_local_rank(send.rank[j]);
            AMREX_ASSERT(send.data[j] != nullptr);
            send.request[j] =
                ParallelDescriptor::Asend(send.data[j], send.size[j], rank, mpi_tag, comm).req();
        }
    }
}

void TranslateRanks (CommData& comm, MPI_Comm c)
{
    const int n = static_cast<int>(comm.rank.size());
    comm.comm_rank.resize(n);
    if (n == 0) { return; }
    if (c == ParallelContext::CommunicatorSub()) {
        ParallelContext::global_to_local_rank(comm.comm_rank.data(), comm.rank.data(), n);
    } else {
        MPI_Group group;
        MPI_Comm_group(c, &group);
        MPI_Group_translate_ranks(ParallelContext::GroupAll(), n, comm.rank.data(),
                                  group, comm.comm_rank.data());
        MPI_Group_free(&group);
        for (int i = 0; i < n; ++i) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(comm.comm_rank[i] != MPI_UNDEFINED,
                "NonLocalBC::TranslateRanks: process not in communicator");
        }
    }
}

void PrepareAlltoallvCounts (CommData& comm, MPI_Comm c)
{
    int nprocs;
    MPI_Comm_size(c, &nprocs);
    comm.counts.assign(nprocs, 0);
    comm.displs.assign(nprocs, 0);
    constexpr auto int_max = static_cast<std::size_t>(std::numeric_limits<int>::max());
    for (int i = 0, n = static_cast<int>(comm.rank.size()); i < n; ++i) {
        AMREX_ALWAYS_ASSERT(comm.offset[i] + comm.size[i] <= int_max);
        comm.counts[comm.comm_rank[i]] = static_cast<int>(comm.size[i]);
        comm.displs[comm.comm_rank[i]] = static_cast<int>(comm.offset[i]);
    }
}

void PrepareCommBuffersUniform (CommData& comm,
                                const FabArrayBase::MapOfCopyComTagContainers& cctc,
                                int n_components, std::size_t object_size,
                                std::size_t block_size, MPI_Comm c)
{
    comm.data.clear();
    comm.size.clear();
    comm.rank.clear();
    comm.request.clear();
    comm.offset.clear();
    comm.cctc.clear();
    comm.stats.clear();

    comm.id = FabArrayBase::getNextCommMetaDataId();

    int nprocs;
    MPI_Comm_size(c, &nprocs);
    const auto N_comms = static_cast<int>(cctc.size());
    comm.data.reserve(N_comms);
    comm.size.reserve(N_comms);
    comm.rank.reserve(N_comms);
    comm.request.reserve(N_comms);
    comm.offset.reserve(N_comms);
    comm.cctc.reserve(N_comms);
    comm.stats.resize(N_comms);

    for (const auto& kv : cctc)
    {
        std::size_t nbytes = 0;
        for (auto const& cct : kv.second) {
            nbytes += cct.sbox.numPts() * object_size * n_components;
        }
        AMREX_ALWAYS_ASSERT(nbytes <= block_size);

        comm.data.push_back(nullptr);
        comm.size.push_back(nbytes);
        comm.rank.push_back(kv.first);
        comm.request.push_back(MPI_REQUEST_NULL);
        comm.cctc.push_back(&kv.second);
    }

    // Block i of the buffer belongs to rank i of c.
    TranslateRanks(comm, c);
    for (int i = 0; i < N_comms; ++i) {
        comm.offset.push_back(block_size * std::size_t(comm.comm_rank[i]));
    }

    std::size_t const total_volume = block_size * std::size_t(nprocs);
    comm.total_bytes = total_volume;
    comm.the_data.reset();
    AcquireCommBuffer(comm);
}

void ReleaseCommBuffer(CommData& comm)
{
    comm.the_data.reset();
    for (auto& p : comm.data) { p = nullptr; }
}

void AcquireCommBuffer(CommData& comm)
{
    if (comm.the_data || comm.total_bytes == 0) { return; }
    comm.the_data.reset(static_cast<char*>(amrex::The_Comms_Arena()->alloc(comm.total_bytes)));
    for (Long i = 0; i < comm.data.size(); ++i) {
        comm.data[i] = comm.the_data.get() + comm.offset[i];
    }
}

namespace {
// Valid buffer pointers for zero counts. They are separate in case an MPI
// implementation rejects aliased send and receive buffers.
char s_dummy_send = 0;
char s_dummy_recv = 0;
}

void Alltoall (CommData const& send, CommData& recv, std::size_t block_size, MPI_Comm comm)
{
    BL_PROFILE("NonLocalBC::Alltoall");

    AMREX_ALWAYS_ASSERT(block_size <= std::size_t(std::numeric_limits<int>::max()));

    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    std::size_t const nbytes = block_size * std::size_t(nprocs);
    // Catches a CommHandler that was prepared by another ParallelCopy variant.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(send.total_bytes == nbytes && recv.total_bytes == nbytes,
        "NonLocalBC::Alltoall: buffers were not prepared by PrepareCommBuffersUniform");

    char const* sbuf = send.the_data ? send.the_data.get() : &s_dummy_send;
    char* rbuf = recv.the_data ? recv.the_data.get() : &s_dummy_recv;

    BL_COMM_PROFILE(BLProfiler::Alltoall, nbytes, ParallelDescriptor::MyProc(),
                    BLProfiler::BeforeCall());
    BL_MPI_REQUIRE( MPI_Alltoall(sbuf, static_cast<int>(block_size), MPI_BYTE,
                                 rbuf, static_cast<int>(block_size), MPI_BYTE, comm) );
    BL_COMM_PROFILE(BLProfiler::Alltoall, nbytes, ParallelDescriptor::MyProc(),
                    BLProfiler::AfterCall());
}

void Alltoallv (CommData const& send, CommData& recv, MPI_Comm comm)
{
    BL_PROFILE("NonLocalBC::Alltoallv");

    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    // Catches a CommHandler that was prepared by another ParallelCopy variant.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(send.counts.size() == Long(nprocs) &&
                                     recv.counts.size() == Long(nprocs),
        "NonLocalBC::Alltoallv: buffers were not prepared by PrepareAlltoallvCounts");

#ifdef BL_COMM_PROFILING
    std::size_t nbytes = 0;
    for (int i = 0, n = static_cast<int>(send.rank.size()); i < n; ++i) {
        nbytes += send.size[i];
    }
#endif

    char const* sbuf = send.the_data ? send.the_data.get() : &s_dummy_send;
    char* rbuf = recv.the_data ? recv.the_data.get() : &s_dummy_recv;

    BL_COMM_PROFILE(BLProfiler::Alltoallv, nbytes, ParallelDescriptor::MyProc(),
                    BLProfiler::BeforeCall());
    BL_MPI_REQUIRE( MPI_Alltoallv(sbuf, send.counts.data(), send.displs.data(), MPI_BYTE,
                                  rbuf, recv.counts.data(), recv.displs.data(), MPI_BYTE,
                                  comm) );
    BL_COMM_PROFILE(BLProfiler::Alltoallv, nbytes, ParallelDescriptor::MyProc(),
                    BLProfiler::AfterCall());
}
#endif

template MultiBlockCommMetaData ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox,
                                             const FabArray<FArrayBox>& src, int destcomp,
                                             int srccomp, int numcomp, const IntVect& ngrow,
                                             MultiBlockIndexMapping const&, Identity const&);

}
