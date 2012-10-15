/* ------------------------------------------------------------------------- */
/* this will only work in serial with one grid. */
/* ------------------------------------------------------------------------- */

void FakeWriteMultifab(int lo[3], int hi[3],  int fablo[3], int fabhi[3],
                       int ng, int ncomp, double ** __restrict__ _cons,
		       char *name);
