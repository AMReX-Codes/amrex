#ifdef USE_MPI
#include <mpi.h>
#endif

void hpgmg_setup (const int log2_box_dim,
                  const int target_boxes_per_rank,
                  const int OMP_Threads,
                  const int OMP_Nested,
                  const int requested_threading_model,
                  const int actual_threading_model);

int
main (int argc, char *argv[])
{

  const int log2_box_dim = 6;
  const int target_boxes_per_rank = 1;

  int OMP_Threads = 1;
  int OMP_Nested = 0;

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp master
    {
      OMP_Threads = omp_get_num_threads ();
      OMP_Nested = omp_get_nested ();
    }
  }
#endif

#ifdef USE_MPI
  int actual_threading_model = -1;
  int requested_threading_model = -1;
  requested_threading_model = MPI_THREAD_SINGLE;
  //requested_threading_model = MPI_THREAD_FUNNELED;
  //requested_threading_model = MPI_THREAD_SERIALIZED;
  //requested_threading_model = MPI_THREAD_MULTIPLE;
  //MPI_Init(&argc, &argv);
#ifdef _OPENMP
  requested_threading_model = MPI_THREAD_FUNNELED;
  //requested_threading_model = MPI_THREAD_SERIALIZED;
  //requested_threading_model = MPI_THREAD_MULTIPLE;
  //MPI_Init_thread(&argc, &argv, requested_threading_model, &actual_threading_model);
#endif
  MPI_Init_thread (&argc, &argv, requested_threading_model,
                   &actual_threading_model);
#ifdef USE_HPM                  // IBM HPM counters for BGQ...
  HPM_Init ();
#endif
#endif

  hpgmg_setup (log2_box_dim,
               target_boxes_per_rank,
               OMP_Threads,
               OMP_Nested, requested_threading_model, actual_threading_model);

  return 0;
}
