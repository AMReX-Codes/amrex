#ifndef AMREX_MPMD_H_
#define AMREX_MPMD_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_MPI

#include <AMReX_FabArray.H>

#include <mpi.h>

namespace amrex::MPMD {

MPI_Comm Initialize (int argc, char* argv[]);

void Finalize ();

bool Initialized ();

int MyProc ();   //! Process ID in MPI_COMM_WORLD
int NProcs ();   //! Number of processes in MPI_COMM_WORLD
int MyProgId (); //! Program ID

class Copier
{
public:
    Copier (BoxArray const& ba, DistributionMapping const& dm);

    template <typename FAB>
    void send (FabArray<FAB> const& mf, int icomp, int ncomp) const;

    template <typename FAB>
    void recv (FabArray<FAB>& mf, int icomp, int ncomp) const;

private:
    std::map<int,FabArrayBase::CopyComTagsContainer> m_SndTags;
    std::map<int,FabArrayBase::CopyComTagsContainer> m_RcvTags;
};

template <typename FAB>
void Copier::send (FabArray<FAB> const& mf, int icomp, int ncomp) const
{
    const auto N_snds = static_cast<int>(m_SndTags.size());

    if (N_snds == 0) return;

    // Prepare buffer

    Vector<char*>       send_data;
    Vector<std::size_t> send_size;
    Vector<int>         send_rank;
    Vector<MPI_Request> send_reqs;
    Vector<FabArrayBase::CopyComTagsContainer const*> send_cctc;

    Vector<std::size_t> offset;
    std::size_t total_volume = 0;
    for (auto const& kv : m_SndTags) {
        auto const& cctc = kv.second;

        std::size_t nbytes = 0;
        for (auto const& cct : cctc) {
            nbytes += cct.sbox.numPts() * ncomp * sizeof(typename FAB::value_type);
        }

        std::size_t acd = ParallelDescriptor::alignof_comm_data(nbytes);
        nbytes = amrex::aligned_size(acd, nbytes); // so that bytes are aligned

        // Also need to align the offset properly
        total_volume = amrex::aligned_size(std::max(alignof(typename FAB::value_type),
                                                    acd), total_volume);

        offset.push_back(total_volume);
        total_volume += nbytes;

        send_data.push_back(nullptr);
        send_size.push_back(nbytes);
        send_rank.push_back(kv.first);
        send_reqs.push_back(MPI_REQUEST_NULL);
        send_cctc.push_back(&cctc);
    }

    Gpu::PinnedVector<char> send_buffer(total_volume);
    char* the_send_data = send_buffer.data();
    for (int i = 0; i < N_snds; ++i) {
        send_data[i] = the_send_data + offset[i];
    }

    // Pack buffer
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && (mf.arena()->isDevice() || mf.arena()->isManaged())) {
        mf.pack_send_buffer_gpu(mf, icomp, ncomp, send_data, send_size, send_cctc);
    } else
#endif
    {
        mf.pack_send_buffer_cpu(mf, icomp, ncomp, send_data, send_size, send_cctc);
    }

    // Send
    for (int i = 0; i < N_snds; ++i) {
        send_reqs[i] = ParallelDescriptor::Asend
            (send_data[i], send_size[i], send_rank[i], 100, MPI_COMM_WORLD).req();
    }
    Vector<MPI_Status> stats(N_snds);
    ParallelDescriptor::Waitall(send_reqs, stats);
}

template <typename FAB>
void Copier::recv (FabArray<FAB>& mf, int icomp, int ncomp) const
{
    const auto N_rcvs = static_cast<int>(m_RcvTags.size());

    if (N_rcvs == 0) return;

    // Prepare buffer

    Vector<char*>       recv_data;
    Vector<std::size_t> recv_size;
    Vector<int>         recv_from;
    Vector<MPI_Request> recv_reqs;

    Vector<std::size_t> offset;
    std::size_t TotalRcvsVolume = 0;
    for (auto const& kv : m_RcvTags) {
        std::size_t nbytes = 0;
        for (auto const& cct : kv.second) {
            nbytes += cct.dbox.numPts() * ncomp * sizeof(typename FAB::value_type);
        }

        std::size_t acd = ParallelDescriptor::alignof_comm_data(nbytes);
        nbytes = amrex::aligned_size(acd, nbytes); // so that nbytes are aligned

        // Also need to align the offset properly
        TotalRcvsVolume = amrex::aligned_size(std::max(alignof(typename FAB::value_type),
                                                       acd), TotalRcvsVolume);

        offset.push_back(TotalRcvsVolume);
        TotalRcvsVolume += nbytes;

        recv_data.push_back(nullptr);
        recv_size.push_back(nbytes);
        recv_from.push_back(kv.first);
        recv_reqs.push_back(MPI_REQUEST_NULL);
    }

    Gpu::PinnedVector<char> recv_buffer(TotalRcvsVolume);
    char* the_recv_data = recv_buffer.data();

    // Recv
    for (int i = 0; i < N_rcvs; ++i) {
        recv_data[i] = the_recv_data + offset[i];
        recv_reqs[i] = ParallelDescriptor::Arecv
            (recv_data[i], recv_size[i], recv_from[i], 100, MPI_COMM_WORLD).req();
    }

    Vector<FabArrayBase::CopyComTagsContainer const*> recv_cctc(N_rcvs, nullptr);
    for (int i = 0; i < N_rcvs; ++i) {
        recv_cctc[i] = &(m_RcvTags.at(recv_from[i]));
    }

    Vector<MPI_Status> stats(N_rcvs);
    ParallelDescriptor::Waitall(recv_reqs, stats);

    // Unpack buffer
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && (mf.arena()->isDevice() || mf.arena()->isManaged())) {
        mf.unpack_recv_buffer_gpu(mf, icomp, ncomp, recv_data, recv_size, recv_cctc,
                                  FabArrayBase::COPY, true);
    } else
#endif
    {
        mf.unpack_recv_buffer_cpu(mf, icomp, ncomp, recv_data, recv_size, recv_cctc,
                                  FabArrayBase::COPY, true);
    }
}

}

#endif
#endif
