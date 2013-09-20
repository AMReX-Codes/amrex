/* 
   Contains the IO routines for fabio module
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#if !defined(WIN32)
#include <unistd.h>
#else
#include <io.h>
typedef int mode_t;
#endif

#ifndef WIN32
#define FILE_MODE  (            S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)
#define DIR_MODE   (FILE_MODE | S_IXUSR |           S_IXGRP | S_IXOTH)
#define O_BINARY   0
#else
#define FILE_MODE (-1)
#define DIR_MODE  (-1)
#endif

#if defined(BL_FORT_USE_UNDERSCORE)
#define FABIO_UNLINK_IF_EMPTY_STR fabio_unlink_if_empty_str_
#define FABIO_OPEN_STR    fabio_open_str_
#define FABIO_MKDIR_STR   fabio_mkdir_str_
#define FABIO_READ_D      fabio_read_d_
#define FABIO_READ_S      fabio_read_s_
#define FABIO_WRITE_RAW_ARRAY_D fabio_write_raw_array_d_
#define FABIO_WRITE_RAW_ARRAY_I fabio_write_raw_array_i_
#define FABIO_READ_RAW_ARRAY_D fabio_read_raw_array_d_
#define FABIO_READ_RAW_ARRAY_I fabio_read_raw_array_i_
#define FABIO_WRITE_RAW_D fabio_write_raw_d_
#define FABIO_WRITE_RAW_S fabio_write_raw_s_
#define FABIO_READ_SKIP_D fabio_read_skip_d_
#define FABIO_READ_SKIP_S fabio_read_skip_s_
#define FABIO_CLOSE       fabio_close_
#define FAB_CONTAINS_NAN  fab_contains_nan_
#define FAB_CONTAINS_INF  fab_contains_inf_
#define VAL_IS_INF        val_is_inf_
#define VAL_IS_NAN        val_is_nan_
#elif defined(BL_FORT_USE_DBL_UNDERSCORE)
#define FABIO_UNLINK_IF_EMPTY_STR fabio_unlink_if_empty_str__
#define FABIO_OPEN_STR    fabio_open_str__
#define FABIO_MKDIR_STR   fabio_mkdir_str__
#define FABIO_READ_D      fabio_read_d__
#define FABIO_READ_S      fabio_read_s__
#define FABIO_WRITE_RAW_ARRAY_D fabio_write_raw_array_d__
#define FABIO_WRITE_RAW_ARRAY_I fabio_write_raw_array_i__
#define FABIO_READ_RAW_ARRAY_D fabio_read_raw_array_d__
#define FABIO_READ_RAW_ARRAY_I fabio_read_raw_array_i__
#define FABIO_WRITE_RAW_D fabio_write_raw_d__
#define FABIO_WRITE_RAW_S fabio_write_raw_s__
#define FABIO_READ_SKIP_D fabio_read_skip_d__
#define FABIO_READ_SKIP_S fabio_read_skip_s__
#define FABIO_CLOSE       fabio_close__
#define FAB_CONTAINS_NAN  fab_contains_nan__
#define FAB_CONTAINS_INF  fab_contains_inf__
#define VAL_IS_INF        val_is_inf__
#define VAL_IS_NAN        val_is_nan__
#elif defined(BL_FORT_USE_LOWERCASE)
#define FABIO_UNLINK_IF_EMPTY_STR fabio_unlink_if_empty_str
#define FABIO_OPEN_STR    fabio_open_str
#define FABIO_MKDIR_STR   fabio_mkdir_str
#define FABIO_READ_D      fabio_read_d
#define FABIO_READ_S      fabio_read_s
#define FABIO_WRITE_RAW_ARRAY_D fabio_write_raw_array_d
#define FABIO_WRITE_RAW_ARRAY_I fabio_write_raw_array_i
#define FABIO_READ_RAW_ARRAY_D fabio_read_raw_array_d
#define FABIO_READ_RAW_ARRAY_I fabio_read_raw_array_i
#define FABIO_WRITE_RAW_D fabio_write_raw_d
#define FABIO_WRITE_RAW_S fabio_write_raw_s
#define FABIO_READ_SKIP_D fabio_read_skip_d
#define FABIO_READ_SKIP_S fabio_read_skip_s
#define FABIO_CLOSE       fabio_close
#define FAB_CONTAINS_NAN  fab_contains_nan
#define FAB_CONTAINS_INF  fab_contains_inf
#define VAL_IS_INF        val_is_inf
#define VAL_IS_NAN        val_is_nan
#endif

static
void
int_2_str(char f[], int n, const int* fi)
{
  int i;
  for ( i = 0; i < n; ++i )
    {
      if ( fi[i] < 0 ) 
	{
	  f[i] = 0;
	  break;
	}
      f[i] = (char)fi[i];
    }
  if ( i == n )
    {
      fprintf(stderr, "name to long, probably not terminated ifilename\n");
      exit(1);
    }
}

void
FABIO_OPEN_STR(int* fdp, const int* ifilename, const int* flagp)
{
  int lflag;
  int lmode;
  char filename[128];

  int_2_str(filename, sizeof(filename), ifilename);
  switch ( *flagp ) 
    {
    case 0: 
      lflag = O_RDONLY; 
      break;
    case 1: 
      lflag = O_WRONLY | O_CREAT | O_TRUNC | O_BINARY; 
      lmode = FILE_MODE;
      break;
    case 2: 
      lflag = O_RDWR;
      break;
    case 3:
      lflag = O_RDWR | O_APPEND;
      break;
    default:
      fprintf(stderr, "FABIO_OPEN_STR: invalid flag, %d, must be <=0<=2", *flagp);
      exit(1);
    }
  *fdp = open(filename, lflag, lmode);
  if ( *fdp == -1 )
    {
      fprintf(stderr, "FABIO_OPEN_STR: failed to open \"%s\": %s\n",
	      filename, strerror(errno));
      exit(1);
    }
}

/* 
 * DOUBLE data
 * FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (1 2 3 4 5 6 7 8)))((0,0) (63,63) (0,0)) 27
 * FLOAT data
 * FAB ((8, (32 8 23 0 1 9 0 127)),(4, (1 2 3 4)        ))((0,0) (63,63) (0,0)) 27
 */

/*
 * NORDER_? : normal byte order floats(f), doubles(d) on this architecture
 */
static const char* str_ieee_d = "64 11 52 0 1 12 0 1023";
static const char* str_ieee_f = "32 8 23 0 1 9 0 127";
#if defined(BL_AIX) || defined(BL_Darwin)
static const int norder_d[8] = { 1, 2, 3, 4, 5, 6, 7, 8};
static const char* str_norder_d = "1 2 3 4 5 6 7 8";
static const int norder_f[4] = { 1, 2, 3, 4};
static const char* str_norder_f = "1 2 3 4";
#else
static const int norder_d[8] = { 8, 7, 6, 5, 4, 3, 2, 1};
static const char* str_norder_d = "8 7 6 5 4 3 2 1";
static const int norder_f[4] = { 4, 3, 2, 1 };
static const char* str_norder_f = "4 3 2 1";
#endif

enum
  {
    FABIO_ERR    = 0,
    /* cf. fabio.f90 */
    FABIO_SINGLE = 2, 
    FABIO_DOUBLE = 1
  };

static
int
scan_buffer(const char* buffer, int border[])
{
  int i;
  int bcount;
  char bstr[1024];

  /* first try for double data */
  i = sscanf(buffer,
	     "FAB ((8, (64 11 52 0 1 12 0 1023)),(%d, (%[^)])))",
	     &bcount,
	     bstr);
  if ( i == 2 )
    {
      i = sscanf(bstr, "%d %d %d %d %d %d %d %d",
		 border + 0, border + 1, border + 2, border + 3,
		 border + 4, border + 5, border + 6, border + 7
		 );
      if ( i != 8 )
	{
	  fprintf(stderr, "FABIO: scan_buffer failed to parse FAB border\n"
		  "Not double precision data\n");
	  exit(1);
	}
      return FABIO_DOUBLE;
    }

  /* second, try for float data */
  i = sscanf(buffer,
	     "FAB ((8, (32 8 23 0 1 9 0 127)),(%d, (%[^)])))",
	     &bcount,
	     bstr);
  if ( i == 2 )
    {
      i = sscanf(bstr, "%d %d %d %d",
		 border + 0, border + 1, border + 2, border + 3
		 );
      if ( i != 4)
	{
	  fprintf(stderr, "FABIO: scan_buffer failed to parse FAB border\n"
		  "Not double precision data\n");
	  exit(1);
	}
      return FABIO_SINGLE;
    }
  
  fprintf(stderr, "FABIO: scan_buffer failed to parse FAB header\n"
	  "Architecture difference for floating point format\n");
  exit(1);
  return FABIO_ERR;
}

void
FABIO_READ_SKIP_D(const int* fdp, const int* offsetp, const int* skipp, 
		  double dp[], const int* countp)
{
  int fd = *fdp;
  char c;
  size_t count = *countp;
  off_t offset = *offsetp;
  off_t skip = *skipp;
  int i,j;
  char buffer[1024];
  int border[8];
  int swap_bytes = 0;
  
  if ( lseek(fd, offset, SEEK_SET) < 0 )
    {
      fprintf(stderr, "FABIO_READ_SKIP_D: failed to seek to %ld: %s\n", 
	      offset, strerror(errno));
      exit(1);
    }
  for (i=0;;i++)
    {
      if ( read(fd, &c, 1) != 1 )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_D: failed to read a char: %s\n",
		  strerror(errno));
	  exit(1);
	}
      if ( c == '\n' ) break;
      if ( i == sizeof(buffer) )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_D: failed FAB header\n");
	  exit(1);
	}
      buffer[i] = c;
    }
  buffer[i] = 0;

  i = scan_buffer(buffer, border);
  if ( i == FABIO_DOUBLE ) 
    {
      /* should be positioned to read the doubles */
      if ( skip && lseek(fd, skip*sizeof(double), SEEK_CUR) < 0 )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_D: failed to seek to comp %ld: %s\n", 
		  offset, strerror(errno));
	  exit(1);
	}
      if ( count*sizeof(double) != read(fd, dp, count*sizeof(double)) )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_D: failed to read %ld doubles: %s\n", 
		  (long)count, strerror(errno));
	  exit(1);
	}
      for ( j = 0; j < 8; ++j )
	{
	  if (border[j] != norder_d[j] )
	    {
	      swap_bytes = 1;
	      break;
	    }
	}
      if ( swap_bytes )
	{
	  unsigned char* cdp = (unsigned char*)dp;
	  for ( i = 0; i < count; i++ )
	    {
	      unsigned char t[8];
	      for ( j = 0; j < 8; j++ )
		{
		  t[j] = cdp[border[j]-1];
		}
	      for ( j = 0; j < 8; j++ )
		{
		  cdp[j] = t[norder_d[j]-1];
		}
	      cdp += 8;
	    }
	}
    }
  else if ( i == FABIO_SINGLE )
    {
      float* fp;
      if ( (fp = (float *) malloc(count*sizeof(float))) == NULL)
	{
	  fprintf(stderr, "FABIO_READ_SKIP_D: failed to allocate fp\n");
	  exit(1);
	}
      /* should be positioned to read the doubles */
      if ( skip && lseek(fd, skip*sizeof(float), SEEK_CUR) < 0 )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_D: failed to seek to comp %ld: %s\n", 
		  offset, strerror(errno));
	  exit(1);
	}
      if ( count*sizeof(float) != read(fd, fp, count*sizeof(float)) )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_D: failed to read %ld doubles: %s\n", 
		  (long)count, strerror(errno));
	  exit(1);
	}
      for ( j = 0; j < 4; ++j )
	{
	  if (border[j] != norder_f[j] )
	    {
	      swap_bytes = 1;
	      break;
	    }
	}
      if ( swap_bytes )
	{
	  unsigned char* csp = (unsigned char*)fp;
	  for ( i = 0; i < count; i++ )
	    {
	      unsigned char t[4];
	      for ( j = 0; j < 4; j++ )
		{
		  t[j] = csp[border[j]-1];
		}
	      for ( j = 0; j < 4; j++ )
		{
		  csp[j] = t[norder_f[j]-1];
		}
	      csp += 4;
	    }
	}
      for ( i = 0; i < count; i++)
	{
	  dp[i] = (double)fp[i];
	}
      free(fp);
    }
}

void
FABIO_READ_SKIP_S(const int* fdp, const int* offsetp, const int* skipp, 
		  float sp[], const int* countp)
{
  int fd = *fdp;
  char c;
  size_t count = *countp;
  off_t offset = *offsetp;
  off_t skip = *skipp;
  int i,j;
  char buffer[1024];
  int border[8];
  int swap_bytes = 0;
  
  if ( lseek(fd, offset, SEEK_SET) < 0 )
    {
      fprintf(stderr, "FABIO_READ_SKIP_S: failed to seek to %ld: %s\n", 
	      offset, strerror(errno));
      exit(1);
    }
  for (i=0;;i++)
    {
      if ( read(fd, &c, 1) != 1 )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_S: failed to read a char: %s\n",
		  strerror(errno));
	  exit(1);
	}
      if ( c == '\n' ) break;
      if ( i == sizeof(buffer) )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_S: failed FAB header\n");
	  exit(1);
	}
      buffer[i] = c;
    }
  buffer[i] = 0;

  i = scan_buffer(buffer, border);
  if ( i == FABIO_DOUBLE ) 
    {
      double* dp;
      if ( (dp = (double *) malloc(count*sizeof(double))) == NULL)
	{
	  fprintf(stderr, "FABIO_READ_SKIP_S: failed to allocate sp\n");
	  exit(1);
	}
      /* should be positioned to read the doubles */
      if ( skip && lseek(fd, skip*sizeof(double), SEEK_CUR) < 0 )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_S: failed to seek to comp %ld: %s\n", 
		  offset, strerror(errno));
	  exit(1);
	}
      if ( count*sizeof(double) != read(fd, dp, count*sizeof(double)) )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_S: failed to read %ld doubles: %s\n", 
		  (long)count, strerror(errno));
	  exit(1);
	}
      for ( j = 0; j < 8; ++j )
	{
	  if (border[j] != norder_d[j] )
	    {
	      swap_bytes = 1;
	      break;
	    }
	}
      if ( swap_bytes )
	{
	  unsigned char* cdp = (unsigned char*)dp;
	  for ( i = 0; i < count; i++ )
	    {
	      unsigned char t[8];
	      for ( j = 0; j < 8; j++ )
		{
		  t[j] = cdp[border[j]-1];
		}
	      for ( j = 0; j < 8; j++ )
		{
		  cdp[j] = t[norder_d[j]-1];
		}
	      cdp += 8;
	    }
	}
      free(dp);
      for ( i = 0; i < count; i++ )
	{
	  if      ( dp[i] >  FLT_MAX )  sp[i] =  FLT_MAX;
	  else if ( dp[i] < -FLT_MAX )  sp[i] = -FLT_MAX;
	  else                          sp[i] = (float)dp[i];
	}
    }
  else if ( i == FABIO_SINGLE )
    {
      /* should be positioned to read the doubles */
      if ( skip && lseek(fd, skip*sizeof(float), SEEK_CUR) < 0 )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_S: failed to seek to comp %ld: %s\n", 
		  offset, strerror(errno));
	  exit(1);
	}
      if ( count*sizeof(float) != read(fd, sp, count*sizeof(float)) )
	{
	  fprintf(stderr, "FABIO_READ_SKIP_S: failed to read %ld doubles: %s\n", 
		  (long)count, strerror(errno));
	  exit(1);
	}
      for ( j = 0; j < 4; ++j )
	{
	  if (border[j] != norder_f[j] )
	    {
	      swap_bytes = 1;
	      break;
	    }
	}
      if ( swap_bytes )
	{
	  unsigned char* csp = (unsigned char*)sp;
	  for ( i = 0; i < count; i++ )
	    {
	      unsigned char t[4];
	      for ( j = 0; j < 4; j++ )
		{
		  t[j] = csp[border[j]-1];
		}
	      for ( j = 0; j < 4; j++ )
		{
		  csp[j] = t[norder_f[j]-1];
		}
	      csp += 4;
	    }
	}
    }
}

/*
** These four guys are used by the particle code.
*/

void
FABIO_WRITE_RAW_ARRAY_D(const int* fdp, const double* vp, const int* countp)
{
  int    fd    = *fdp;
  size_t count = *countp;
  int    ilen  = sizeof(double) * count;

  lseek(fd, 0, SEEK_END);

  if ( ilen != write(fd, vp, ilen) )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_ARRAY_D: failed to write %d bytes: %s\n", 
	      ilen, strerror(errno));
      exit(1);
    }
}

void
FABIO_WRITE_RAW_ARRAY_I(const int* fdp, const int* vp, const int* countp)
{
  int    fd    = *fdp;
  size_t count = *countp;
  int    ilen  = sizeof(int) * count;

  lseek(fd, 0, SEEK_END);

  if ( ilen != write(fd, vp, ilen) )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_ARRAY_I: failed to write %d bytes: %s\n", 
	      ilen, strerror(errno));
      exit(1);
    }
}

void
FABIO_READ_RAW_ARRAY_D(const int* fdp, double* vp, const int* countp)
{
  int    fd    = *fdp;
  size_t count = *countp;
  int    ilen  = sizeof(double) * count;

  if ( ilen != read(fd, vp, ilen) )
    {
      fprintf(stderr, "FABIO_READ_RAW_ARRAY_D: failed to write %d bytes: %s\n", 
	      ilen, strerror(errno));
      exit(1);
    }
}

void
FABIO_READ_RAW_ARRAY_I(const int* fdp, int* vp, const int* countp)
{
  int    fd    = *fdp;
  size_t count = *countp;
  int    ilen  = sizeof(int) * count;

  if ( ilen != read(fd, vp, ilen) )
    {
      fprintf(stderr, "FABIO_READ_RAW_ARRAY_I: failed to write %d bytes: %s\n", 
	      ilen, strerror(errno));
      exit(1);
    }
}

void
FABIO_WRITE_RAW_D(const int* fdp, int* offsetp, const double* vp, const int* countp, 
		  const int* dmp, const int lo[], const int hi[], const int nd[],
		  const int* ncp)
{
  int fd = *fdp;
  int dm = *dmp;
  int nc = *ncp;
  size_t count = *countp;
  off_t offset;
  char buffer[256];
  int ilen;
  double* dp = (double*)vp;

  offset = lseek(fd, 0, SEEK_END);
  if ( snprintf(buffer, 256, "FAB ((8, (%s)),(8, (%s)))", str_ieee_d, str_norder_d) >= 256 )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_D: buffer too small");
      exit(1);
    }
  ilen = strlen(buffer);
  if ( ilen != write(fd, buffer, ilen) )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_D: failed to write %d bytes: %s\n", 
	      ilen, strerror(errno));
      exit(1);
    }
  switch ( dm ) 
    {
    case 1:
      ilen = snprintf(buffer, 256, "((%d) (%d) (%d)) %d\n", 
                      lo[0], hi[0], nd[0], nc);
      break;
    case 2:
      ilen = snprintf(buffer, 256, "((%d,%d) (%d,%d) (%d,%d)) %d\n", 
                      lo[0], lo[1], hi[0], hi[1], nd[0], nd[1], nc);
      break;
    case 3:
      ilen = snprintf(buffer, 256, "((%d,%d,%d) (%d,%d,%d) (%d,%d,%d)) %d\n", 
                      lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], nd[0], nd[1], nd[2], nc);
      break;
    default:
      fprintf(stderr, "FABIO_WRITE_RAW_D: strange dimension = %d\n", dm);
      exit(1);
    }

  if ( ilen >= 256 )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_D: buffer too small");
      exit(1);
    }

  ilen = write(fd, buffer, strlen(buffer));

  if ( ilen !=  strlen(buffer) )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_D: write of buffer failed\n");
      exit(1);
    }
  
  ilen = nc*count*sizeof(double);
  if ( ilen != write(fd, dp, ilen) )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_D: failed to write %ld doubles: %s\n", 
	      (long)nc*count, strerror(errno));
      exit(1);
    }

  if ( offset > INT_MAX )
  {
      fprintf(stderr, "FABIO_WRITE_RAW_D: offset will overflow offsetp");
      exit(1);
  }

  *offsetp = offset;
}

void
FABIO_WRITE_RAW_S(const int* fdp, int* offsetp, const float* vp, const int* countp, 
		  const int* dmp, const int lo[], const int hi[], const int nd[], 
		  const int* ncp)
{
  int fd = *fdp;
  int dm = *dmp;
  int nc = *ncp;
  size_t count = *countp;
  off_t offset;
  char buffer[256];
  int ilen;
  float* sp = (float*)vp;

  offset = lseek(fd, 0, SEEK_END);
  if ( snprintf(buffer, 256, "FAB ((8, (%s)),(4, (%s)))", str_ieee_f, str_norder_f) >= 256 )
  {
      fprintf(stderr, "FABIO_WRITE_RAW_S: buffer too small");
      exit(1);
  }
  ilen = strlen(buffer);
  if ( ilen != write(fd, buffer, ilen) )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_S: failed to write %d bytes: %s\n", 
	      ilen, strerror(errno));
      exit(1);
    }
  switch ( dm ) 
    {
    case 1:
      ilen = snprintf(buffer, 256, "((%d) (%d) (%d)) %d\n", 
                      lo[0], hi[0], nd[0], nc);
      break;
    case 2:
      ilen = snprintf(buffer, 256, "((%d,%d) (%d,%d) (%d,%d)) %d\n", 
                      lo[0], lo[1], hi[0], hi[1], nd[0], nd[1], nc);
      break;
    case 3:
      ilen = snprintf(buffer, 256, "((%d,%d,%d) (%d,%d,%d) (%d,%d,%d)) %d\n", 
                      lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], nd[0], nd[1], nd[2], nc);
      break;
    default:
      fprintf(stderr, "FABIO_WRITE_RAW_S: strange dimension = %d\n", dm);
      exit(1);
    }

  if ( ilen >= 256 )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_S: buffer too small");
      exit(1);
    }

  ilen = write(fd, buffer, strlen(buffer));
  if ( ilen !=  strlen(buffer) )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_S: write of buffer failed\n");
      exit(1);
    }

  ilen = nc*count*sizeof(float);
  if ( ilen != write(fd, sp, ilen) )
    {
      fprintf(stderr, "FABIO_WRITE_RAW_S: failed to write %ld floats: %s\n", 
	      (long)nc*count, strerror(errno));
      exit(1);
    }

  if ( offset > INT_MAX )
  {
      fprintf(stderr, "FABIO_WRITE_RAW_S: offset will overflow offsetp");
      exit(1);
  }

  *offsetp = offset;
}


void
FABIO_READ_D(const int* fdp, const int* offsetp, double dp[], const int* countp)
{
  int skip = 0;
  FABIO_READ_SKIP_D(fdp, offsetp, &skip, dp, countp);
}

void
FABIO_READ_S(const int* fdp, const int* offsetp, float sp[], const int* countp)
{
  int skip = 0;
  FABIO_READ_SKIP_S(fdp, offsetp, &skip, sp, countp);
}

void
FABIO_CLOSE(const int* fdp)
{
  int fd = *fdp;
  if ( close(fd) < 0 )
    {
      fprintf(stderr, "FABIO_CLOSE: failed to close %d: %s\n", 
	      fd, strerror(errno));
      exit(1);
    }
}

#if defined(WIN32)
#include <direct.h>
#define mkdir(a,b) _mkdir((a))
/* static const char* path_sep_str = "\\"; */
#else
/* static const char* path_sep_str = "/"; */
#endif

void
FABIO_MKDIR_STR(const int* idirname, int* statp)
{
  mode_t mode = DIR_MODE;
  int st = *statp;
  char dirname[128];

  int_2_str(dirname, sizeof(dirname), idirname);

  *statp = 0;
  /* we allow the mkdir on an existing directory */
  if ( mkdir(dirname, mode) <0 && errno != EEXIST ) 
    {
      if ( st )
	{
	  *statp = 1;
	  return;
	}
      else
	{
	  fprintf(stderr, "FABIO_MKDIR_STR: mkdir(%s,%d): %s\n", 
		  dirname, mode, strerror(errno));
	  exit(1);
	}
    }
}

void
FABIO_UNLINK_IF_EMPTY_STR(const int* ifilename)
{
  int fd;
  char filename[128];
  int lmode = FILE_MODE;
  int pos;

  int_2_str(filename, sizeof(filename), ifilename);

  if ((fd = open(filename, O_RDONLY, lmode)) < 0)
    {
      fprintf(stderr, "FABIO_UNLINK_IF_EMPTY: open() failed: \"%s\": %s\n",
	      filename, strerror(errno));
      exit(1);
    }

  if ((pos = lseek(fd, 0, SEEK_END)) < 0)
    {
      fprintf(stderr, "FABIO_UNLINK_IF_EMPTY: lseek() failed: \"%s\": %s\n",
	      filename, strerror(errno));
      exit(1);
    }

  close(fd);

  if (pos == 0)
    {
      if (unlink(filename) < 0)
        {
          fprintf(stderr, "FABIO_UNLINK_IF_EMPTY: unlink() failed: \"%s\": %s\n",
                  filename, strerror(errno));
          exit(1);
        }
    }
}

void
FAB_CONTAINS_NAN (double dptr[], const int* countp, int* result)
{
#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(_AIX) || defined(__PATHSCALE__)
    int i, r;
    int rr=0;
#ifdef _OPENMP
#pragma omp parallel private(r) reduction(+:rr)
#endif
    {
      r = 0;
#ifdef _OPENMP
#pragma omp for private(i)
#endif
      for (i = 0; i < *countp; i++) {
	if (isnan(dptr[i])) {
	  r = 1;
	}
      }
      rr += r;
    }
    *result = (rr>0) ? 1 : 0;
#endif
}

void
FAB_CONTAINS_INF (double dptr[], const int* countp, int* result)
{
#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(_AIX) || defined(__PATHSCALE__)
    int i, r;
    int rr=0;
#ifdef _OPENMP
#pragma omp parallel private(r) reduction(+:rr)
#endif
    {
      r = 0;
#ifdef _OPENMP
#pragma omp for private(i)
#endif
      for (i = 0; i < *countp; i++) {
        if (isinf(dptr[i])) {
	  r = 1;
	}
      }
      rr += r;
    }
    *result = (rr>0) ? 1 : 0;
#endif
}

void
VAL_IS_INF (double* val, int* result)
{
#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(_AIX) || defined(__PATHSCALE__)

    *result = (isinf(*val) ? 1 : 0);

#endif
}

void
VAL_IS_NAN (double* val, int* result)
{
#if defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(_AIX) || defined(__PATHSCALE__)

    *result = (isnan(*val) ? 1 : 0);

#endif
}
