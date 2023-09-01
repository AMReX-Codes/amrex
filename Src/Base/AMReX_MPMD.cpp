#include <AMReX_MPMD.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <utility>
#include <vector>

#ifdef AMREX_USE_MPI

namespace amrex::MPMD {

namespace {
    bool initialized = false;
    bool mpi_initialized_by_us = false;
    MPI_Comm app_comm = MPI_COMM_NULL;
    int myproc;
    int nprocs;
}

namespace {

template <typename T>
int num_unique_elements (std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto last = std::unique(v.begin(), v.end());
    return last - v.begin();
}

}

MPI_Comm Initialize (int argc, char* argv[])
{
    initialized = true;
    int flag;
    MPI_Initialized(&flag);
    if (!flag) {
        MPI_Init(&argc, &argv);
        mpi_initialized_by_us = true;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int* p;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_APPNUM, &p, &flag);
    int appnum = *p;

    std::vector<int> all_appnum(nprocs);
    MPI_Allgather(&appnum, 1, MPI_INT, all_appnum.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int napps = num_unique_elements(all_appnum);

    // MPI_APPNUM does not appear to work with slurm on some systems.
    if (napps != 2) {
        std::vector<int> all_argc(nprocs);
        MPI_Allgather(&argc, 1, MPI_INT, all_argc.data(), 1, MPI_INT, MPI_COMM_WORLD);
        napps = num_unique_elements(all_argc);
        if (napps == 2) {
            appnum = static_cast<int>(argc != all_argc[0]);
        }
    }

    if (napps != 2) {
        std::string exename;
        if (argc > 0) {
            exename = std::string(argv[0]);
        }
        unsigned long long hexe = std::hash<std::string>{}(exename);
        std::vector<unsigned long long> all_hexe(nprocs);
        MPI_Allgather(&hexe, 1, MPI_UNSIGNED_LONG_LONG,
                      all_hexe.data(), 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
        napps = num_unique_elements(all_hexe);
        if (napps == 2) {
            appnum = static_cast<int>(hexe != all_hexe[0]);
        }
    }

    if (napps == 2) {
        MPI_Comm_split(MPI_COMM_WORLD, appnum, myproc, &app_comm);
    } else {
        std::cout << "amrex::MPMD only supports two programs." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return app_comm;
}

void Finalize ()
{
    MPI_Comm_free(&app_comm);
    if (mpi_initialized_by_us) {
        MPI_Finalize();
        mpi_initialized_by_us = false;
    }
    initialized = false;
}

bool Initialized () { return initialized; }

int MyProc ()
{
    return myproc;
}

int NProcs ()
{
    return nprocs;
}

int MyProgId ()
{
    return (myproc == ParallelDescriptor::MyProc()) ? 0 : 1;
}

Copier::Copier (BoxArray const& ba, DistributionMapping const& dm)
{
    int rank_offset = myproc - ParallelDescriptor::MyProc();
    int this_root, other_root;
    if (rank_offset == 0) { // First program
        this_root = 0;
        other_root = ParallelDescriptor::NProcs();
    } else {
        this_root = rank_offset;
        other_root = 0;
    }

    Vector<Box> bv = ba.boxList().data();

    int this_nboxes = static_cast<int>(ba.size());
    Vector<int> procs = dm.ProcessorMap();
    if (rank_offset != 0) {
        for (int i = 0; i < this_nboxes; ++i) {
            procs[i] += rank_offset;
        }
    }

    Vector<Box> obv;
    Vector<int> oprocs;
    int other_nboxes;
    if (myproc == this_root) {
        if (rank_offset == 0) // the first program
        {
            MPI_Send(&this_nboxes, 1, MPI_INT, other_root, 0, MPI_COMM_WORLD);
            MPI_Recv(&other_nboxes, 1, MPI_INT, other_root, 1, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            obv.resize(other_nboxes);
            MPI_Send(bv.data(), this_nboxes,
                     ParallelDescriptor::Mpi_typemap<Box>::type(),
                     other_root, 2, MPI_COMM_WORLD);
            MPI_Recv(obv.data(), other_nboxes,
                     ParallelDescriptor::Mpi_typemap<Box>::type(),
                     other_root, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            oprocs.resize(other_nboxes);
            MPI_Send(procs.data(), this_nboxes, MPI_INT, other_root, 4, MPI_COMM_WORLD);
            MPI_Recv(oprocs.data(), other_nboxes, MPI_INT, other_root, 5, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }
        else // the second program
        {
            MPI_Recv(&other_nboxes, 1, MPI_INT, other_root, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Send(&this_nboxes, 1, MPI_INT, other_root, 1, MPI_COMM_WORLD);
            obv.resize(other_nboxes);
            MPI_Recv(obv.data(), other_nboxes,
                     ParallelDescriptor::Mpi_typemap<Box>::type(),
                     other_root, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(bv.data(), this_nboxes,
                     ParallelDescriptor::Mpi_typemap<Box>::type(),
                     other_root, 3, MPI_COMM_WORLD);
            oprocs.resize(other_nboxes);
            MPI_Recv(oprocs.data(), other_nboxes, MPI_INT, other_root, 4, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Send(procs.data(), this_nboxes, MPI_INT, other_root, 5, MPI_COMM_WORLD);
        }
    }

    ParallelDescriptor::Bcast(&other_nboxes, 1);
    if (obv.empty()) {
        obv.resize(other_nboxes);
        oprocs.resize(other_nboxes);
    }
    ParallelDescriptor::Bcast(obv.data(), obv.size());
    ParallelDescriptor::Bcast(oprocs.data(), oprocs.size());

    BoxArray oba(BoxList(std::move(obv)));

    // At this point, ba and bv hold our boxes, and oba holds the other
    // program's boxes. procs holds mpi ranks of our boxes, and oprocs holds
    // mpi ranks of the other program's boxes.  All mpi ranks are in
    // MPI_COMM_WORLD.

    // Build communication meta-data

    AMREX_ALWAYS_ASSERT(ba.ixType().cellCentered());

    std::vector<std::pair<int,Box> > isects;

    for (int i = 0; i < this_nboxes; ++i) {
        if (procs[i] == myproc) {
            oba.intersections(bv[i], isects);
            for (auto const& isec : isects) {
                const int oi = isec.first;
                const Box& bx = isec.second;
                const int orank = oprocs[oi];
                m_SndTags[orank].push_back
                    (FabArrayBase::CopyComTag(bx, bx, oi, i));
                m_RcvTags[orank].push_back
                    (FabArrayBase::CopyComTag(bx, bx, i, oi));
            }
        }
    }

    for (auto& kv : m_SndTags) {
        std::sort(kv.second.begin(), kv.second.end());
    }
    for (auto& kv : m_RcvTags) {
        std::sort(kv.second.begin(), kv.second.end());
    }
}

}

#endif
