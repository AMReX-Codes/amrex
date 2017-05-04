/* 
   Contains architecture specific timer routines used by bl_timer.f90
*/

#if defined(_BL_ANSI_TIME)

#include <stdlib.h>
#include <time.h>

void
cpu_second(double* rslt)
{
  const double clock_rate = 1.0/(double)(CLOCKS_PER_SEC);
  static clock_t start;
  static int inited = 0;
  if ( inited == 0 )
    {
      inited = 1;
      start = clock();
    }
  *rslt = (double)(clock()-start)*clock_rate;
}

void
wall_second(double* rslt)
{
  static time_t start;
  static int inited = 0;
  if ( inited == 0 )
    {
      inited = 1;
      start = time(0);
    }
  *rslt = (double)(time(0)-start);
}

#else

#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

void
cpu_second(double* rslt)
{
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  *rslt = (double)(ru.ru_utime.tv_sec + ru.ru_stime.tv_sec)
    + (double)(ru.ru_utime.tv_usec + ru.ru_stime.tv_usec)*1.0e-6;
}

#if defined(_BL_USE_MPI_WTIME)
#include <mpi.h>

void
wall_second(double* rslt)
{
  int fr;
  assert(  MPI_Initialized(&fr) == 0 );
  *rslt = MPI_Wtime();
}
#else
void
wall_second(double* rslt)
{
  struct timeval tp;
  /* cannot fail, so why check */
  gettimeofday(&tp, 0);
  *rslt = (double)(tp.tv_sec) + (double)(tp.tv_usec)*1.0e-6;
}
#endif
#endif


void
wall_second_tick( double* r)
{
  static double tickval = -1.0;
  double t1, t2;
  int cnt, icnt;
  if ( tickval < 0.0 )
    {
      tickval = 1.0e6;
      for (icnt=0; icnt<100; icnt++)
	{
	  cnt = 1000;
	  wall_second( &t1 );
	  while (cnt--) 
	    {
	      wall_second( &t2 );
	      if (t2 > t1) break;
	    }
	  if (cnt && t2 > t1 && t2 - t1 < tickval)
	    {
	      tickval = t2 - t1;
	    }
	}
    }
  *r = tickval;
}

void
cpu_second_tick( double* r)
{
  static double tickval = -1.0;
  double t1, t2;
  int cnt, icnt;
  if ( tickval < 0.0 )
    {
      tickval = 1.0e6;
      for (icnt=0; icnt<1000; icnt++)
	{
	  cnt = 1000;
	  cpu_second( &t1 );
	  while (cnt--)
	    {
	      cpu_second( &t2 );
	      if (t2 > t1) break;
	    }
	  if (cnt && t2 > t1 && t2 - t1 < tickval)
	    {
	      tickval = t2 - t1;
	    }
	}
    }
  *r = tickval;
}

#include <stdlib.h>

void
sys_abort()
{
  abort();
}
