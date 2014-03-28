#ifndef _support_H_
#define _support_H_ 1

#include <FArrayBox.H>
#include <BoxArray.H>

void
SAXPYjpdf(FArrayBox&       dest,
          Real             scale,
          const FArrayBox& src);

std::vector<double>
condMean(FArrayBox& src,
         bool       cdir);

extern int _charArrayLen;

char *
makePGM( FArrayBox* fab, int icomp, double usemin, double usemax);

// helper functions for manipulating BoxArray's
BoxArray *
complementIn(const BoxArray& ba1, const BoxArray& ba2);
BoxArray *
join(const BoxArray &ba1, const BoxArray &ba2);
//BoxArray *
//intersect( const BoxArray &ba1, const BoxArray &ba2);

void
dummyFunc();

// transfers between coarse and fine grids
void tocoarse(int ratio,FArrayBox&fine,FArrayBox&crse,IntVect&iv);
void tofine(int ratio,FArrayBox&crse,FArrayBox&fine,IntVect&iv,Real defval);
void injectCoarse(FArrayBox&fine,int ratio, const IntVect&iv, 
		  FArrayBox&crse, const Box&tbox);
// conditional vector merges
void cvmgnfab(FArrayBox&res,FArrayBox&n,FArrayBox&p,FArrayBox&trg);
void cvmgpfab(FArrayBox&res,FArrayBox&p,FArrayBox&n,FArrayBox&trg);
// misc
void gradfillfab(FArrayBox&fab, const IntVect&iv);

#endif
