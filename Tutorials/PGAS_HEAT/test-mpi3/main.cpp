
#include <iostream>
#include <cassert>
#include <mpi.h>

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

    for (int i = 0; i < nprocs; ++i) {
        if ( i == rank) {
            std::cout << "rank " << rank << ", team_comm size " << real_team_size
		      << ", rank in team " << real_team_rank << std::endl;
        }
	MPI_Barrier(MPI_COMM_WORLD);
    }

    const int N = 8;
    MPI_Win win_shared;
    void* baseptr;
    MPI_Win_allocate_shared(8*sizeof(double), sizeof(double),
			    MPI_INFO_NULL, team_comm, &baseptr, &win_shared);

    double* p = (double*) baseptr;
    for (int i = 0; i < N; ++i) {
	p[i] = (double) i;
    }

    MPI_Win_free(&win_shared);
#ifdef MANUAL_SPLIT
    MPI_Group_free(&team_grp);
#endif
    MPI_Comm_free(&team_comm);
    MPI_Finalize();
}
