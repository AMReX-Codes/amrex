/* ------------------------------------------------------------------------- */
/* FakeWriteMultifab.c                                                       */
/*                                                                           */
/* this will only work in serial with one grid.                              */
/* ------------------------------------------------------------------------- */
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>

void FakeWriteMultifab(int lo[3], int hi[3], int fablo[3], int fabhi[3],
                       int ng, int ncomp, double ** __restrict__ _cons,
		       char *name)
{
  int c;
  int lo0=lo[0]; int hi0=hi[0];
  int lo1=lo[1]; int hi1=hi[1];
  int lo2=lo[2]; int hi2=hi[2];
  int pencilg =  (hi0-lo0+1)+2*ng;
  int planeg  = ((hi1-lo1+1)+2*ng)*pencilg;

  double * __restrict__ cons[5];
  for(c = 0; c < ncomp; ++c) {
    cons[c] = _cons[c];
  }

  int i, j, k;
  double rholoc, uvel, vvel, wvel, eloc;
  double dmin[ncomp], dmax[ncomp];

  for(c = 0; c < ncomp; ++c) {
    int ijkg = (lo[0]-fablo[0])+(lo[1]-fablo[1])*pencilg+(lo[2]-fablo[2])*planeg;
    dmin[c] = cons[c][ijkg];
    dmax[c] = cons[c][ijkg];
  }

  for(k=lo2;k<=hi2;k++){
    for(j=lo1;j<=hi1;j++){
      for(i=lo0;i<=hi0;i++){
        int ijkg = (i-fablo[0]) + (j-fablo[1])*pencilg + (k-fablo[2])*planeg;
        for(c = 0; c < ncomp; ++c) {
          dmin[c] = fmin(dmin[c], cons[c][ijkg]);
          dmax[c] = fmax(dmax[c], cons[c][ijkg]);
        }
      }
    }
  }

  FILE *outfileH, *outfileD;
  char fnameH[64], fnameD[64];
  sprintf(fnameH, "%s%s", name, "_H");
  sprintf(fnameD, "%s%s", name, "_D_0000");
  outfileH = fopen(fnameH, "w");

  fprintf(outfileH, "1\n0\n%d\n0\n", ncomp);
  fprintf(outfileH, "(1 0\n((%d,%d,%d) (%d,%d,%d) (0,0,0))\n)\n",
                    lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
  fprintf(outfileH, "1\nFabOnDisk: %s 0\n\n", fnameD);
  fprintf(outfileH, "1,%d\n", ncomp);
  for(c = 0; c < ncomp; ++c) {
    fprintf(outfileH, "%20.17e,",dmin[c]);
  }
  fprintf(outfileH, "\n");
  fprintf(outfileH, "1,%d\n", ncomp);
  for(c = 0; c < ncomp; ++c) {
    fprintf(outfileH, "%20.17e,",dmax[c]);
  }
  fprintf(outfileH, "\n");
  fclose(outfileH);

  outfileD = fopen(fnameD, "w");
  fprintf(outfileD, "FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))");
  fprintf(outfileD, "((%d,%d,%d) (%d,%d,%d) (0,0,0)) %d\n",
                     fablo[0], fablo[1], fablo[2], fabhi[0], fabhi[1], fabhi[2],
		     ncomp);
  int volume = (fabhi[0] - fablo[0] + 1) *
               (fabhi[1] - fablo[1] + 1) * 
               (fabhi[2] - fablo[2] + 1);
  for(c = 0; c < ncomp; ++c) {
    int ec = fwrite(cons[c], sizeof(double), volume, outfileD);
    if(ec == -1) {
      printf("errno = %s\n", strerror(errno));
    }
  }

  fclose(outfileD);
}

