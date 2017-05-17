
#include <iostream>
#include <cassert>
#include <mpi.h>
#include <atomic>

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm team_comm;
    MPI_Group team_grp;

#ifdef MANUAL_SPLIT

    int team_size = nprocs/2;  // There will be two teams.
    assert(team_size*2 == nprocs);

    MPI_Group grp;
    MPI_Comm_group(MPI_COMM_WORLD, &grp);

    int team_ranks[team_size];
    if (rank < team_size) {  // Team 0
	for (int i = 0; i < team_size; ++i) {
	    team_ranks[i] = i;
	}
    } else {                 // Team 1
	for (int i = 0; i < team_size; ++i) {
	    team_ranks[i] = i + team_size;
	}
    }

    MPI_Group_incl(grp, team_size, team_ranks, &team_grp);

    MPI_Comm_create(MPI_COMM_WORLD, team_grp, &team_comm);

#else

    int r = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, 
                                MPI_INFO_NULL, &team_comm);

    if (r != MPI_SUCCESS) std::cout << "MPI_Comm_split_type failed" << std::endl;
#endif

    int real_team_size;
    int real_team_rank;
    MPI_Comm_size(team_comm, &real_team_size);
    MPI_Comm_rank(team_comm, &real_team_rank);

#if 0
    for (int i = 0; i < nprocs; ++i) {
        if ( i == rank) {
            std::cout << "rank " << rank << ", team_comm size " << real_team_size
		      << ", rank in team " << real_team_rank << std::endl;
        }
	MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    const int N = 8;
    const int NN = (rank == 0) ? N : 0;
    MPI_Win win_shared;
    double* p;
    MPI_Win_allocate_shared(NN*sizeof(double), sizeof(double),
			    MPI_INFO_NULL, team_comm, &p, &win_shared);

    if (rank != 0) {
	MPI_Aint sz;
	int disp;
	MPI_Win_shared_query(win_shared, MPI_PROC_NULL, &sz, &disp, &p);
    }

    std::atomic_thread_fence(std::memory_order_release);
    MPI_Barrier(team_comm);
    std::atomic_thread_fence(std::memory_order_acquire);

    if (rank == 1) {
	for (int i = 0; i < N; ++i) {
	    p[i] = (double) (i*(rank+1));
	}
    }

    std::atomic_thread_fence(std::memory_order_release);
    MPI_Barrier(team_comm);
    std::atomic_thread_fence(std::memory_order_acquire);    

    for (int i = 0; i < nprocs; ++i) {
        if ( i == rank) {
            std::cout << "rank " << rank << ", data =";
	    for (int j = 0; j < N; ++j) {
		std::cout << " " << p[j];
	    }
	    std::cout << std::endl;
        }
	MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Win_free(&win_shared);
#ifdef MANUAL_SPLIT
    MPI_Group_free(&team_grp);
#endif
    MPI_Comm_free(&team_comm);
    MPI_Finalize();
}
