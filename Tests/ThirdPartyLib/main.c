#include <stdio.h>
#include <mpi.h>

void foo(MPI_Comm comm);    // foo is a C++ function.
void bar(MPI_Comm comm);    // bar is a Fortran subroutine.

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank_world, np_world;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);
    MPI_Comm_size(MPI_COMM_WORLD, &np_world);

    int fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);

    if (rank_world == 0) {
        printf("\nCalling foo with MPI_COMM_WORLD\n");
    }
    foo(MPI_COMM_WORLD);

    if (rank_world == 0) {
        printf("\nCalling bar with MPI_COMM_WORLD\n");
    }
    bar(fcomm);

    int subcomm, subrank;
    int color = rank_world % 2;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank_world, &subcomm);
    MPI_Comm_rank(subcomm, &subrank);
    int fsubcomm = MPI_Comm_c2f(subcomm);

    if (color == 0) {
        if (subrank == 0) {
            printf("\nCalling bar with sub-communicator 0\n");
        }
        bar(subcomm);

        if (subrank == 0) {
            printf("\nCalling foo with sub-communicator 0\n");
        }
        foo(subcomm);
    } else {
        if (subrank == 0) {
            printf("\nCalling foo with sub-communicator 1\n");
        }
        foo(subcomm);

        if (subrank == 0) {
            printf("\nCalling bar with sub-communicator 1\n");
        }
        bar(subcomm);
    }

    if (rank_world == 0) {
        printf("\nCalling bar with MPI_COMM_WORLD\n");
    }
    bar(MPI_COMM_WORLD);

    if (rank_world == 0) {
        printf("\nCalling foo with MPI_COMM_WORLD\n");
    }
    foo(MPI_COMM_WORLD);

    MPI_Comm_free(&subcomm);

    MPI_Finalize();
}
