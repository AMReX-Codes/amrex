/* 
   Contains architecture specific PPM writing routines used by ppm_util.f90
*/
#include <stdlib.h>
#include <stdio.h>

static const char PGM_MAGIC1 = 'P';
static const char RPGM_MAGIC2 = '5';
static const char RPGM_MAGIC3 = '6';
#define NCOLOR 256

#if defined(BL_FORT_USE_UNDERSCORE)
#define STORE_PPM_STR store_ppm_str_
#define STORE_PGM_STR store_pgm_str_
#define LOAD_PALETTE_STR load_palette_str_
#elif defined(BL_FORT_USE_DBL_UNDERSCORE)
#define STORE_PPM_STR store_ppm_str__
#define STORE_PGM_STR store_pgm_str__
#define LOAD_PALETTE_STR load_palette_str__
#elif defined(BL_FORT_USE_LOWERCASE)
#define STORE_PPM_STR store_ppm_str
#define STORE_PGM_STR store_pgm_str
#define LOAD_PALETTE_STR load_palette_str
#endif

static void
int_2_str(char filename[], int n, const int ifilename[])
{
  int i;
  for (i = 0; i < n; i++) 
    {
      if ( ifilename[i] < 0 )
	{
	  filename[i] = 0;
	  break;
	}
      filename[i] = (char)ifilename[i];
    }
  if ( i == n )
    {
      fprintf(stderr, "name to long, probably not terminated ifilename\n");
      exit(1);
    }
}

void
STORE_PPM_STR (const int ifilename[], const int* width, const int* height, int iimage[],
	       const int r[], const int g[], const int b[])
{
  FILE *image_fp;	/* file descriptor for image (output) */
  int i, j, wid, hgt;
  char filename[128];
  unsigned char* image;
  
  /* create image file */
  int_2_str(filename, sizeof(filename), ifilename);
  if ((image_fp = fopen(filename, "w")) == NULL)
    {
      fprintf(stderr, "cannot open output file %s\n", filename);
      exit(1);
    }

  /* translate Fortran image data to chars */
  wid = *width;
  hgt = *height;
  image = (unsigned char*) malloc(3*wid*hgt*sizeof(unsigned char));
  if ( image == NULL ) 
    {
      fprintf(stderr, "STORE_PPM_STR: failed to allocate image buffer\n");
      exit(1);
    }
  for ( i = 0; i < wid*hgt; i++ )
    {
      j = iimage[i];
      if ( j < 0 || j > 255 ) 
	{
	  fprintf(stderr,"out of bounds on image[%d] = %d\n", i, j);
	  exit(1);
	}
      image[3*i+0] = (unsigned char)r[j];
      image[3*i+1] = (unsigned char)g[j];
      image[3*i+2] = (unsigned char)b[j];
    }
  fprintf(image_fp, "%c%c\n%d %d\n%d\n", PGM_MAGIC1, RPGM_MAGIC3,
	  wid, hgt, 255);
  fwrite(image, 1, 3*wid*hgt, image_fp);
  free(image);
  fclose(image_fp);
}

void
STORE_PGM_STR (const int ifilename[], const int* width, const int* height, int iimage[])
{
  FILE *image_fp;	/* file descriptor for image (output) */
  int i, wid, hgt;
  char filename[128];
  unsigned char* image;
  
  /* create image file */
  int_2_str(filename, sizeof(filename), ifilename);
  if ((image_fp = fopen(filename, "w")) == NULL)
    {
      fprintf(stderr, "cannot open output file %s\n", filename);
      exit(1);
    }

  /* translate Fortran image data to chars */
  wid = *width;
  hgt = *height;
  image = (unsigned char*)malloc(wid*hgt*sizeof(unsigned char));
  if ( image == NULL ) 
    {
      fprintf(stderr, "STORE_PGM_STR: failed to allocate image buffer\n");
      exit(1);
    }
  for ( i = 0; i < wid*hgt; ++i )
    {
      image[i] = (unsigned char)iimage[i];
    }
  fprintf(image_fp, "%c%c\n%d %d\n%d\n", PGM_MAGIC1, RPGM_MAGIC2,
	  wid, hgt, 255);
  fwrite(image, 1, wid*hgt, image_fp);
  free(image);
  fclose(image_fp);
}

void
LOAD_PALETTE_STR (const int ifilename[], int r[], int g[], int b[], int a[])
{
  char filename[128];
  FILE *pal_fp;	/* file descriptor for image (output) */
  unsigned char c[NCOLOR];
  int i;
  long length;
  int num_elements;

  int_2_str(filename, sizeof(filename), ifilename);
  if ((pal_fp = fopen(filename, "rb")) == NULL)
    {
      fprintf(stderr,
	      "cannot open palette file %s\n", 
	      filename);
      exit(1);
    }

  fseek(pal_fp, 0, SEEK_END);
  length = ftell(pal_fp);
  fseek(pal_fp, 0, SEEK_SET);
  
  /* check for RGB or RGBA palette */
  num_elements = length/(NCOLOR*sizeof(unsigned char));

  if ( num_elements != 3 && num_elements != 4 )
    {
      fprintf(stderr,
	      "cannot process palette file %s, num(r,g,b,a) = %d, must be 3 or 4\n",
	      filename, num_elements);
      exit(1);
    }

  if (fread(c, 1, NCOLOR, pal_fp) != NCOLOR)
    {
      fprintf(stderr, "LOAD_PALETTE_STR: fread() failed!\n");
      exit(1);
    }

  for ( i = 0; i < NCOLOR; ++i )
    {
      r[i] = c[i];
    }
  if (fread(c, 1, NCOLOR, pal_fp) != NCOLOR)
    {
      fprintf(stderr, "LOAD_PALETTE_STR: fread() failed!\n");
      exit(1);
    }
  for ( i = 0; i < NCOLOR; ++i )
    {
      g[i] = c[i];
    }
  if (fread(c, 1, NCOLOR, pal_fp) != NCOLOR)
    {
      fprintf(stderr, "LOAD_PALETTE_STR: fread() failed!\n");
      exit(1);
    }
  for ( i = 0; i < NCOLOR; ++i )
    {
      b[i] = c[i];
    }
  if ( num_elements == 4 ) 
    {
      if (fread(c, 1, NCOLOR, pal_fp) != NCOLOR)
        {
          fprintf(stderr, "LOAD_PALETTE_STR: fread() failed!\n");
          exit(1);
        }
      for ( i = 0; i < NCOLOR; ++i )
	{
	  a[i] = c[i];
	}
    }
  fclose(pal_fp);
}
