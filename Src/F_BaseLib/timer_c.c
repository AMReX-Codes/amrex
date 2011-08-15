/* 
   Contains architecture specific timer routines used by bl_timer.f90
*/
#if defined(BL_FORT_USE_UNDERSCORE)
#define CPU_SECOND cpu_second_
#define WALL_SECOND wall_second_
#define CPU_SECOND_TICK cpu_second_tick_
#define WALL_SECOND_TICK wall_second_tick_
#elif defined(BL_FORT_USE_DBL_UNDERSCORE)
#define CPU_SECOND cpu_second__
#define WALL_SECOND wall_second__
#define CPU_SECOND_TICK cpu_second_tick__
#define WALL_SECOND_TICK wall_second_tick__
#elif defined(BL_FORT_USE_LOWERCASE)
#define CPU_SECOND cpu_second
#define WALL_SECOND wall_second
#define CPU_SECOND_TICK cpu_second_tick
#define WALL_SECOND_TICK wall_second_tick
#endif

#if defined(_BL_ANSI_TIME)

#include <stdlib.h>
#include <time.h>

void
CPU_SECOND(double* rslt)
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
WALL_SECOND(double* rslt)
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

#elif defined(WIN32)

#include <windows.h>
#include <time.h>

void
WALL_SECOND(double* rslt)
{
  static int inited = 0;
  static LONGLONG llStart;
  static double rate;
  LARGE_INTEGER li;

  if ( inited == 0 )
    {
      LARGE_INTEGER lii;
      inited = 1;
      QueryPerformanceFrequency(&lii);
      rate = 1.0/lii.QuadPart;
      QueryPerformanceCounter(&lii);
      llStart = lii.QuadPart;
    }
  QueryPerformanceCounter(&li);
  *rslt = (double)((li.QuadPart-llStart)*rate);
}

void
CPU_SECOND (double* r)
{
  static int inited = 0;
  static clock_t start;

  if (inited == 0 )
    {
      inited = 1;
      start = clock();
    }

  *r = (double)(clock() - start)/CLOCKS_PER_SEC;
}

#else

#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

void
CPU_SECOND(double* rslt)
{
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  *rslt = (double)(ru.ru_utime.tv_sec + ru.ru_stime.tv_sec)
    + (double)(ru.ru_utime.tv_usec + ru.ru_stime.tv_usec)*1.0e-6;
}

#if defined(_BL_USE_MPI_WTIME)
#include <mpi.h>

void
WALL_SECOND(double* rslt)
{
  int fr;
  assert(  MPI_Initialized(&fr) == 0 );
  *rslt = MPI_Wtime();
}
#else
void
WALL_SECOND(double* rslt)
{
  struct timeval tp;
  /* cannot fail, so why check */
  gettimeofday(&tp, 0);
  *rslt = (double)(tp.tv_sec) + (double)(tp.tv_usec)*1.0e-6;
}
#endif
#endif


void
WALL_SECOND_TICK( double* r)
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
	  WALL_SECOND( &t1 );
	  while (cnt--) 
	    {
	      WALL_SECOND( &t2 );
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
CPU_SECOND_TICK( double* r)
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
	  CPU_SECOND( &t1 );
	  while (cnt--)
	    {
	      CPU_SECOND( &t2 );
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

#if defined(BL_FORT_USE_UNDERSCORE)
#define SYS_ABORT sys_abort_
#elif defined(BL_FORT_USE_DBL_UNDERSCORE)
#define SYS_ABORT sys_abort__
#elif defined(BL_FORT_USE_LOWERCASE)
#define SYS_ABORT sys_abort
#endif

void
SYS_ABORT()
{
  abort();
}
