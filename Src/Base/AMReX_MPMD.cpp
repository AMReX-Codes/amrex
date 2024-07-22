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
    int appnum;
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

/*
Initialize_without_split function assigns and checks the required
AMReX_MPMD variables. This function is internally leveraged by
Initialize function.

This function needs to be used EXPLICITLY ONLY with pyAMReX (python)
so that the communication split can be performed using a python
library, for example, mpi4py.
*/
void Initialize_without_split (int argc, char* argv[])
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
    appnum = *p;

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

    if (napps != 2) {
        std::cout << "amrex::MPMD only supports two programs." << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

}

MPI_Comm Initialize (int argc, char* argv[])
{
    Initialize_without_split(argc,argv);
    MPI_Comm_split(MPI_COMM_WORLD, appnum, myproc, &app_comm);

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

/*
AppNum function is provided so that appnum (color)
can be passed to python library (mpi4py) to perform
a pythonic version of MPI_Comm_split.
*/
int AppNum ()
{
    return appnum;
}

int MyProgId ()
{
    return (myproc == ParallelDescriptor::MyProc()) ? 0 : 1;
}

Copier::Copier (BoxArray const& ba, DistributionMapping const& dm,
        bool send_ba)
        : m_ba(ba), m_dm(dm)
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
    int other_nboxes = this_nboxes;
    if (myproc == this_root) {
        if (rank_offset == 0) // the first program
        {
            MPI_Send(&this_nboxes, 1, MPI_INT, other_root, 0, MPI_COMM_WORLD);
            if (!send_ba)
            {
                MPI_Recv(&other_nboxes, 1, MPI_INT, other_root, 1, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
                obv.resize(other_nboxes);
            }
            MPI_Send(bv.data(), this_nboxes,
                     ParallelDescriptor::Mpi_typemap<Box>::type(),
                     other_root, 2, MPI_COMM_WORLD);
            if (!send_ba)
            {
                MPI_Recv(obv.data(), other_nboxes,
                     ParallelDescriptor::Mpi_typemap<Box>::type(),
                     other_root, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Send(procs.data(), this_nboxes, MPI_INT, other_root, 4, MPI_COMM_WORLD);
            oprocs.resize(other_nboxes);
            MPI_Recv(oprocs.data(), other_nboxes, MPI_INT, other_root, 5, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }
        else // the second program
        {
            if (!send_ba)
            {
                MPI_Recv(&other_nboxes, 1, MPI_INT, other_root, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
                obv.resize(other_nboxes);
            }
            MPI_Send(&this_nboxes, 1, MPI_INT, other_root, 1, MPI_COMM_WORLD);
            if (!send_ba)
            {
                MPI_Recv(obv.data(), other_nboxes,
                     ParallelDescriptor::Mpi_typemap<Box>::type(),
                     other_root, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Send(bv.data(), this_nboxes,
                     ParallelDescriptor::Mpi_typemap<Box>::type(),
                     other_root, 3, MPI_COMM_WORLD);
            oprocs.resize(other_nboxes);
            MPI_Recv(oprocs.data(), other_nboxes, MPI_INT, other_root, 4, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Send(procs.data(), this_nboxes, MPI_INT, other_root, 5, MPI_COMM_WORLD);
        }
    }

    if (!send_ba) {
        ParallelDescriptor::Bcast(&other_nboxes, 1);
        if (obv.empty()){
            obv.resize(other_nboxes);
        }
        ParallelDescriptor::Bcast(obv.data(), obv.size());
    }

    if (oprocs.empty()) {
        oprocs.resize(other_nboxes);
    }
    ParallelDescriptor::Bcast(oprocs.data(), oprocs.size());

    BoxArray oba;
    if (!obv.empty()) {
        oba.define(BoxList(std::move(obv)));
    }

    // At this point, ba and bv hold our boxes, and oba holds the other
    // program's boxes. procs holds mpi ranks of our boxes, and oprocs holds
    // mpi ranks of the other program's boxes.  All mpi ranks are in
    // MPI_COMM_WORLD.

    // Build communication meta-data
    if (!send_ba){
        AMREX_ALWAYS_ASSERT(ba.ixType() == oba.ixType());
        m_is_thread_safe = ba.ixType().cellCentered();
    }else{
        m_is_thread_safe = true;
    }

    std::vector<std::pair<int,Box> > isects;

    for (int i = 0; i < this_nboxes; ++i) {
        if (procs[i] == myproc) {
            if (!send_ba){
                oba.intersections(bv[i], isects);
            }
            else{
                isects.resize(0);
                isects.emplace_back(i,bv[i]);
            }
            for (auto const& isec : isects) {
                const int oi = isec.first;
                const Box& bx = isec.second;
                const int orank = oprocs[oi];
                m_SndTags[orank].emplace_back(bx, bx, oi, i);
                m_RcvTags[orank].emplace_back(bx, bx, i, oi);
            }
        }
    }

    if (!send_ba){
        for (auto& kv : m_SndTags) {
            std::sort(kv.second.begin(), kv.second.end());
        }
        for (auto& kv : m_RcvTags) {
            std::sort(kv.second.begin(), kv.second.end());
        }
    }
}

Copier::Copier (bool)
    : m_is_thread_safe(true)
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

    Vector<Box> bv;
    int this_nboxes;

    if (myproc == this_root) {
        int tags[2];
        if (rank_offset == 0) // the first program
        {
            tags[0] = 1;
            tags[1] = 3;
        }
        else // the second program
        {
            tags[0] = 0;
            tags[1] = 2;
        }

        MPI_Recv(&this_nboxes, 1, MPI_INT, other_root, tags[0], MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        bv.resize(this_nboxes);
        MPI_Recv(bv.data(), this_nboxes,
                 ParallelDescriptor::Mpi_typemap<Box>::type(),
                 other_root, tags[1], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    ParallelDescriptor::Bcast(&this_nboxes, 1);
    if (bv.empty()) {
        bv.resize(this_nboxes);
    }

    ParallelDescriptor::Bcast(bv.data(), bv.size());
    m_ba.define(BoxList(std::move(bv)));
    m_dm.define(m_ba);
    Vector<int> procs = m_dm.ProcessorMap();
    if (rank_offset != 0) {
        for (int i = 0; i < this_nboxes; ++i) {
            procs[i] += rank_offset;
        }
    }

    Vector<int> oprocs(this_nboxes);
    if (myproc == this_root) {
        if (rank_offset == 0) // the first program
        {
            MPI_Send(procs.data(), this_nboxes, MPI_INT, other_root, 4, MPI_COMM_WORLD);
            MPI_Recv(oprocs.data(), this_nboxes, MPI_INT, other_root, 5, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }
        else // the second program
        {
            MPI_Recv(oprocs.data(), this_nboxes, MPI_INT, other_root, 4, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Send(procs.data(), this_nboxes, MPI_INT, other_root, 5, MPI_COMM_WORLD);
        }
    }

    ParallelDescriptor::Bcast(oprocs.data(), oprocs.size());

    // procs holds mpi ranks of our boxes, and oprocs holds
    // mpi ranks of the other program's boxes.  All mpi ranks are in
    // MPI_COMM_WORLD.

    // Build communication meta-data

    for (int i = 0; i < this_nboxes; ++i) {
        if (procs[i] == myproc) {
            const Box& bx = m_ba[i];
            const int orank = oprocs[i];
            m_SndTags[orank].emplace_back(bx, bx, i, i);
            m_RcvTags[orank].emplace_back(bx, bx, i, i);
        }
    }
}

BoxArray const& Copier::boxArray () const
{
    return m_ba;
}

DistributionMapping const& Copier::DistributionMap () const
{
    return m_dm;
}

}

#endif
