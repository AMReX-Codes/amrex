#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "support.H"

#include <Box.H>
#include <FArrayBox.H>
#include <BoxList.H>
#include <BoxDomain.H>
#include <BoxArray.H>

void
SAXPYjpdf(FArrayBox&       dest,
          Real             scale,
          const FArrayBox& src)
{
  if (dest.box() != src.box() || dest.nComp() != 3 || src.nComp() != 3) {
    BoxLib::Abort("Bad boxes");
  }

  const Box& box = src.box();
  if (scale < 0) {
    // Must "flip" data in src fab before adding
    // FIXME:  We simply assume that the bins were build symmetric around zero for this

    if (box.length(1)%2 != 0) {
      BoxLib::Abort("Box length in j must be even");
    }
    int jm = (box.smallEnd()[1] + box.bigEnd()[1] + 1)/2;
    Box sbox(box); sbox.setSmall(1,jm);
    for (IntVect iv=sbox.smallEnd(), End=sbox.bigEnd(); iv<=End; sbox.next(iv)) {
      IntVect ivr(D_DECL(iv[0],2*jm - iv[1] - 1,iv[2]));

      dest(iv,0) -= scale * src(ivr,0);
      dest(ivr,0) -= scale * src(iv,0);

      dest(iv,1) -= scale * src(ivr,1);
      dest(ivr,1) -= scale * src(iv,1);

      dest(iv,2)  = scale * src(ivr,2);
      dest(ivr,2)  = scale * src(iv,2);
    }
  }
  else {
    dest.saxpy(scale, src, box, box, 0, 0, 3);
  }
}

std::vector<double>
condMean(FArrayBox& src,
         bool       cdir)
{
  const Box& box = src.box();
  const IntVect se=box.smallEnd();
  const IntVect be=box.bigEnd();

  std::vector<double> cm(be[cdir]-se[cdir]+1);
  bool idir = !cdir;
  int jm = (se[idir] + be[idir] + 1)/2;
  for (int i=se[cdir]; i<=be[cdir]; ++i) {
    Real s, v; s = v = 0;
    for (int j=se[idir]; j<=be[idir]; ++j) {
      IntVect iv; iv[cdir]=i; iv[idir]=j;
      v += src(iv,0);
      s += src(iv,1+idir);
    }
    cm[i-se[cdir]] = (v != 0 ? s / v : 0);
  }
  return cm;
}

// fill fab with uniform gradient
// cell(i,j) -> i*iv[0]+j*iv[1]
// occasionally useful
#if (BL_SPACEDIM==2)
void gradfillfab(FArrayBox&fab, const IntVect&iv){
  int iv0 = iv[0];
  int iv1 = iv[1];
  Box bx = fab.box();
  int nc = fab.nComp();
  int i,j;
  int ilo = bx.loVect()[0];
  int jlo = bx.loVect()[1];
  int ihi = bx.hiVect()[0];
  int jhi = bx.hiVect()[1];
  int ilen = ihi-ilo+1;
  int jlen = jhi-jlo+1;
  Real *ptr = fab.dataPtr();
  for(int n=0; n<nc; n++){
    for(int j = jlo; j<=jhi; j++){
      for(int i=ilo; i<=ihi; i++){
	int offset = i-ilo+ilen*(j-jlo+jlen*n);
	ptr[offset] = i*iv0+j*iv1;
      }
    }
  }
}
#elif (BL_SPACEDIM==3)
void gradfillfab(FArrayBox&fab, const IntVect&iv){
  int iv0 = iv[0];
  int iv1 = iv[1];
  int iv2 = iv[2];
  Box bx = fab.box();
  int nc = fab.nComp();
  int ilo = bx.loVect()[0];
  int jlo = bx.loVect()[1];
  int klo = bx.loVect()[2];
  int ihi = bx.hiVect()[0];
  int jhi = bx.hiVect()[1];
  int khi = bx.hiVect()[2];
  int ilen = ihi-ilo+1;
  int jlen = jhi-jlo+1;
  int klen = khi-klo+1;
  Real *ptr = fab.dataPtr();
  for(int n=0; n<nc; n++){
    for(int k = klo; k<=khi; k++){
      for(int j = jlo; j<=jhi; j++){
	for(int i=ilo; i<=ihi; i++){
	  int offset = i-ilo+ilen*(j-jlo+jlen*(k-klo+klen*n));
	  ptr[offset] = i*iv0+j*iv1+k*iv2;
	}
      }
    }
  }
}
#else
whoa! this cant happen
#endif

// do a conditional vector merge operation on fab's and put into res
// res = p if trg>0              on a point by point basis
// res = n if trg<= 0            on a point by point basis
// is implicit assumption that all fab's are exactly the same size
#if (BL_SPACEDIM==2)
void cvmgpfab(FArrayBox&res,FArrayBox&p,FArrayBox&n,FArrayBox&trg)
{
  Box bx = res.box();
  int nc = res.nComp();
  int i,j;
  int ilo = bx.loVect()[0];
  int jlo = bx.loVect()[1];
  int ihi = bx.hiVect()[0];
  int jhi = bx.hiVect()[1];
  int ilen = ihi-ilo+1;
  int jlen = jhi-jlo+1;
  Real *rptr = res.dataPtr();
  Real *pptr = p.dataPtr();
  Real *nptr = n.dataPtr();
  Real *tptr = trg.dataPtr();
  int istart = 0;
  int iend = ihi-ilo+ilen*(jhi-jlo+jlen*(nc-1));
  for( i=istart; i<=iend; i++){
    if( tptr[i]>0. ){
      rptr[i] = pptr[i];
    }
    else{
      rptr[i] = nptr[i];
    }
  }
}
#elif (BL_SPACEDIM==3)
void cvmgpfab(FArrayBox&res,FArrayBox&p,FArrayBox&n,FArrayBox&trg)
{
  Box bx = res.box();
  int nc = res.nComp();
  int i,j,k;
  int ilo = bx.loVect()[0];
  int jlo = bx.loVect()[1];
  int klo = bx.loVect()[2];
  int ihi = bx.hiVect()[0];
  int jhi = bx.hiVect()[1];
  int khi = bx.hiVect()[2];
  int ilen = ihi-ilo+1;
  int jlen = jhi-jlo+1;
  int klen = khi-klo+1;
  Real *rptr = res.dataPtr();
  Real *pptr = p.dataPtr();
  Real *nptr = n.dataPtr();
  Real *tptr = trg.dataPtr();
  int istart = 0;
  int iend = ihi-ilo+ilen*(jhi-jlo+jlen*(khi-klo+klen*(nc-1)));
  for( i=istart; i<=iend; i++){
    if( tptr[i]>0. ){
      rptr[i] = pptr[i];
    }
    else{
      rptr[i] = nptr[i];
    }
  }
}
#else
whoa! this cant happen
#endif

// do a conditional vector merge operation on fab's and put into res
// res = n if trg<0              on a point by point basis
// res = p if trg>= 0            on a point by point basis
// is implicit assumption that all fab's are exactly the same size
#if (BL_SPACEDIM==2)
void cvmgnfab(FArrayBox&res,FArrayBox&n,FArrayBox&p,FArrayBox&trg)
{
  Box bx = res.box();
  int nc = res.nComp();
  int i,j;
  int ilo = bx.loVect()[0];
  int jlo = bx.loVect()[1];
  int ihi = bx.hiVect()[0];
  int jhi = bx.hiVect()[1];
  int ilen = ihi-ilo+1;
  int jlen = jhi-jlo+1;
  Real *rptr = res.dataPtr();
  Real *pptr = p.dataPtr();
  Real *nptr = n.dataPtr();
  Real *tptr = trg.dataPtr();
  int istart = 0;
  int iend = ihi-ilo+ilen*(jhi-jlo+jlen*(nc-1));
  for( i=istart; i<=iend; i++){
    if( tptr[i]<0. ){
      rptr[i] = nptr[i];
    }
    else{
      rptr[i] = pptr[i];
    }
  }
}
#elif (BL_SPACEDIM==3)
void cvmgnfab(FArrayBox&res,FArrayBox&n,FArrayBox&p,FArrayBox&trg)
{
  Box bx = res.box();
  int nc = res.nComp();
  int i,j,k;
  int ilo = bx.loVect()[0];
  int jlo = bx.loVect()[1];
  int klo = bx.loVect()[2];
  int ihi = bx.hiVect()[0];
  int jhi = bx.hiVect()[1];
  int khi = bx.hiVect()[2];
  int ilen = ihi-ilo+1;
  int jlen = jhi-jlo+1;
  int klen = khi-klo+1;
  Real *rptr = res.dataPtr();
  Real *pptr = p.dataPtr();
  Real *nptr = n.dataPtr();
  Real *tptr = trg.dataPtr();
  int istart = 0;
  int iend = ihi-ilo+ilen*(jhi-jlo+jlen*(khi-klo+klen*(nc-1)));
  for( i=istart; i<=iend; i++){
    if( tptr[i]<0. ){
      rptr[i] = nptr[i];
    }
    else{
      rptr[i] = pptr[i];
    }
  }
}
#else
whoa! this cant happen
#endif

#if (BL_SPACEDIM==2)
void tocoarse(int ratio,FArrayBox&fine,FArrayBox&crse,IntVect&iv)
{
  // error checking
  const Box & cbox = crse.box();
  const Box & fbox = fine.box();
  if( BoxLib::refine(cbox,ratio) != fbox){
    BoxLib::Error("tocoarse: coarse and fine boxes incommensurate");
  }
  int ioff = iv[0];
  int joff = iv[1];
  if( ioff < 0 || ioff >= ratio || joff < 0 || joff >= ratio ){
    BoxLib::Error("tocoarse: offset improper for ratio");
  }
  int nc = crse.nComp();
  if( nc != fine.nComp() ){
    BoxLib::Error("tocoarse: number of components not match");
  }
  int i,j;
  int cilo = cbox.loVect()[0];
  int cjlo = cbox.loVect()[1];
  int cihi = cbox.hiVect()[0];
  int cjhi = cbox.hiVect()[1];
  int cilen = cihi-cilo+1;
  int cjlen = cjhi-cjlo+1;
  Real *cptr = crse.dataPtr();
  int filo = fbox.loVect()[0];
  int fjlo = fbox.loVect()[1];
  int fihi = fbox.hiVect()[0];
  int fjhi = fbox.hiVect()[1];
  int filen = fihi-filo+1;
  int fjlen = fjhi-fjlo+1;
  Real *fptr = fine.dataPtr();
  for( int n=0; n<nc; n++){
    for( j=cjlo; j<=cjhi; j++){
      for( i=cilo; i<=cihi; i++){
	int coffset = i-cilo+cilen*(j-cjlo+cjlen*n);
	int foffset = ioff+i*ratio-filo+filen*(joff+j*ratio-fjlo+fjlen*n);
	cptr[coffset] = fptr[foffset];
      }
    }
  }
}
#elif (BL_SPACEDIM==3)
void tocoarse(int ratio,FArrayBox&fine,FArrayBox&crse,IntVect&iv)
{
  // error checking
  const Box & cbox = crse.box();
  const Box & fbox = fine.box();
  if( BoxLib::refine(cbox,ratio) != fbox){
    BoxLib::Error("tocoarse: coarse and fine boxes incommensurate");
  }
  int ioff = iv[0];
  int joff = iv[1];
  int koff = iv[2];
  if(ioff<0||ioff>=ratio||joff<0||joff>=ratio||koff<0||koff>=ratio){
    BoxLib::Error("tocoarse: offset improper for ratio");
  }
  int nc = crse.nComp();
  if( nc != fine.nComp() ){
    BoxLib::Error("tocoarse: number of components not match");
  }
  int i,j,k;
  int cilo = cbox.loVect()[0];
  int cjlo = cbox.loVect()[1];
  int cklo = cbox.loVect()[2];
  int cihi = cbox.hiVect()[0];
  int cjhi = cbox.hiVect()[1];
  int ckhi = cbox.hiVect()[2];
  int cilen = cihi-cilo+1;
  int cjlen = cjhi-cjlo+1;
  int cklen = ckhi-cklo+1;
  Real *cptr = crse.dataPtr();
  int filo = fbox.loVect()[0];
  int fjlo = fbox.loVect()[1];
  int fklo = fbox.loVect()[2];
  int fihi = fbox.hiVect()[0];
  int fjhi = fbox.hiVect()[1];
  int fkhi = fbox.hiVect()[2];
  int filen = fihi-filo+1;
  int fjlen = fjhi-fjlo+1;
  int fklen = fkhi-fklo+1;
  Real *fptr = fine.dataPtr();
  for( int n=0; n<nc; n++){
    for( k=cklo; k<=ckhi; k++){
      for( j=cjlo; j<=cjhi; j++){
	for( i=cilo; i<=cihi; i++){
	  int coffset = i-cilo+cilen*(j-cjlo+cjlen*(k-cklo+cklen*n));
	  int foffset = ioff+i*ratio-filo+filen*(joff+j*ratio-fjlo+fjlen*
						 (koff+k*ratio-fklo+fklen*n));
	  cptr[coffset] = fptr[foffset];
	}
      }
    }
  }
}
#else
whoa! this cant happen
#endif


#if (BL_SPACEDIM==2)
void tofine(int ratio,FArrayBox&crse,FArrayBox&fine,IntVect&iv,Real defval)
{
  // error checking
  const Box & cbox = crse.box();
  const Box & fbox = fine.box();
  if( BoxLib::refine(cbox,ratio) != fbox){
    BoxLib::Error("tofine: coarse and fine boxes incommensurate");
  }
  int ioff = iv[0];
  int joff = iv[1];
  if( ioff < 0 || ioff >= ratio || joff < 0 || joff >= ratio ){
    BoxLib::Error("tofine: offset improper for ratio");
  }
  int nc = crse.nComp();
  if( nc != fine.nComp() ){
    BoxLib::Error("tofine: number of components not match");
  }
  int i,j;
  int cilo = cbox.loVect()[0];
  int cjlo = cbox.loVect()[1];
  int cihi = cbox.hiVect()[0];
  int cjhi = cbox.hiVect()[1];
  int cilen = cihi-cilo+1;
  int cjlen = cjhi-cjlo+1;
  Real *cptr = crse.dataPtr();
  int filo = fbox.loVect()[0];
  int fjlo = fbox.loVect()[1];
  int fihi = fbox.hiVect()[0];
  int fjhi = fbox.hiVect()[1];
  int filen = fihi-filo+1;
  int fjlen = fjhi-fjlo+1;
  Real *fptr = fine.dataPtr();
  fine.setVal(defval);
  for( int n=0; n<nc; n++){
    for( j=cjlo; j<=cjhi; j++){
      for( i=cilo; i<=cihi; i++){
	int coffset = i-cilo+cilen*(j-cjlo+cjlen*n);
	int foffset = ioff+i*ratio-filo+filen*(joff+j*ratio-fjlo+fjlen*n);
	fptr[foffset] = cptr[coffset];
      }
    }
  }
}
#elif (BL_SPACEDIM==3)
void tofine(int ratio,FArrayBox&crse,FArrayBox&fine,IntVect&iv,Real defval)
{
  // error checking
  const Box & cbox = crse.box();
  const Box & fbox = fine.box();
  if( BoxLib::refine(cbox,ratio) != fbox){
    BoxLib::Error("tofine: coarse and fine boxes incommensurate");
  }
  int ioff = iv[0];
  int joff = iv[1];
  int koff = iv[2];
  if(ioff<0||ioff>=ratio||joff<0||joff>=ratio||koff<0||koff>=ratio){
    BoxLib::Error("tofine: offset improper for ratio");
  }
  int nc = crse.nComp();
  if( nc != fine.nComp() ){
    BoxLib::Error("tofine: number of components not match");
  }
  int i,j,k;
  int cilo = cbox.loVect()[0];
  int cjlo = cbox.loVect()[1];
  int cklo = cbox.loVect()[2];
  int cihi = cbox.hiVect()[0];
  int cjhi = cbox.hiVect()[1];
  int ckhi = cbox.hiVect()[2];
  int cilen = cihi-cilo+1;
  int cjlen = cjhi-cjlo+1;
  int cklen = ckhi-cklo+1;
  Real *cptr = crse.dataPtr();
  int filo = fbox.loVect()[0];
  int fjlo = fbox.loVect()[1];
  int fklo = fbox.loVect()[2];
  int fihi = fbox.hiVect()[0];
  int fjhi = fbox.hiVect()[1];
  int fkhi = fbox.hiVect()[2];
  int filen = fihi-filo+1;
  int fjlen = fjhi-fjlo+1;
  int fklen = fkhi-fklo+1;
  Real *fptr = fine.dataPtr();
  fine.setVal(defval);
  for( int n=0; n<nc; n++){
    for( k=cklo; k<=ckhi; k++){
      for( j=cjlo; j<=cjhi; j++){
	for( i=cilo; i<=cihi; i++){
	  int coffset = i-cilo+cilen*(j-cjlo+cjlen*(k-cklo+cklen*n));
	  int foffset = ioff+i*ratio-filo+filen*(joff+j*ratio-fjlo+
					 fjlen*(koff+k*ratio-fklo+fklen*n));
	  fptr[foffset] = cptr[coffset];
	}
      }
    }
  }
}
#else
whoa! this cant happen
#endif



#if (BL_SPACEDIM==2)
void injectCoarse(FArrayBox&fine,int ratio, const IntVect&iv, 
		  FArrayBox&crse,const Box&tbox)
{
  // error checking
  const Box & cbox = crse.box();
  const Box & fbox = fine.box();
  if( ! cbox.contains(tbox) ){
    BoxLib::Error("injectCoarse: coarse and target boxes incommensurate");
  }
  int ioff = iv[0];
  int joff = iv[1];
  if( ioff < 0 || ioff >= ratio || joff < 0 || joff >= ratio ){
    BoxLib::Error("injectCoarse: offset improper for ratio");
  }
  int nc = crse.nComp();
  if( nc != fine.nComp() ){
    BoxLib::Error("injectCoarse: number of components not match");
  }
  int i,j;
  int cilo = cbox.loVect()[0];
  int cjlo = cbox.loVect()[1];
  int cihi = cbox.hiVect()[0];
  int cjhi = cbox.hiVect()[1];
  int cilen = cihi-cilo+1;
  int cjlen = cjhi-cjlo+1;
  Real *cptr = crse.dataPtr();
  int filo = fbox.loVect()[0];
  int fjlo = fbox.loVect()[1];
  int fihi = fbox.hiVect()[0];
  int fjhi = fbox.hiVect()[1];
  int filen = fihi-filo+1;
  int fjlen = fjhi-fjlo+1;
  Real *fptr = fine.dataPtr();
  int tilo = tbox.loVect()[0];
  int tjlo = tbox.loVect()[1];
  int tihi = tbox.hiVect()[0];
  int tjhi = tbox.hiVect()[1];
  // run over coarse grid points in tbox and stuff it in fine
  for( int n=0; n<nc; n++){
    for( j=tjlo; j<=tjhi; j++){
      for( i=tilo; i<=tihi; i++){
	int coffset = i-cilo+cilen*(j-cjlo+cjlen*n);
	int fioff = i*ratio+ioff;
	int fjoff = j*ratio+joff;
	if( fioff >= filo && fioff <= fihi && fjoff >= fjlo && fjoff <= fjhi ){
	  int foffset = fioff-filo+filen*(fjoff-fjlo+fjlen*n);
	  fptr[foffset] = cptr[coffset];
	}
      }
    }
  }
}
#elif (BL_SPACEDIM==3)
void injectCoarse(FArrayBox&fine,int ratio, const IntVect&iv, 
		  FArrayBox&crse,const Box&tbox)
{
  // error checking
  if( !tbox.ok() ) return;
  const Box & cbox = crse.box();
  const Box & fbox = fine.box();
  if( ! cbox.contains(tbox) ){
    BoxLib::Error("tofine: coarse and target boxes incommensurate");
  }
  // test to make sure that tbox image + iv is in fbox
  Box timage = BoxLib::refine(tbox,ratio);
  const IntVect &tsm = timage.smallEnd();
  IntVect testiv = tsm+iv;
  if( !fbox.contains(testiv) ) return;

  int ioff = iv[0];
  int joff = iv[1];
  int koff = iv[2];
  if(ioff<0||ioff>=ratio||joff<0||joff>=ratio||koff<0||koff>=ratio){
    BoxLib::Error("tofine: offset improper for ratio");
  }
  int nc = crse.nComp();
  if( nc != fine.nComp() ){
    BoxLib::Error("tofine: number of components not match");
  }
  int i,j,k;
  int cilo = cbox.loVect()[0];
  int cjlo = cbox.loVect()[1];
  int cklo = cbox.loVect()[2];
  int cihi = cbox.hiVect()[0];
  int cjhi = cbox.hiVect()[1];
  int ckhi = cbox.hiVect()[2];
  int cilen = cihi-cilo+1;
  int cjlen = cjhi-cjlo+1;
  int cklen = ckhi-cklo+1;
  Real *cptr = crse.dataPtr();
  int filo = fbox.loVect()[0];
  int fjlo = fbox.loVect()[1];
  int fklo = fbox.loVect()[2];
  int fihi = fbox.hiVect()[0];
  int fjhi = fbox.hiVect()[1];
  int fkhi = fbox.hiVect()[2];
  int filen = fihi-filo+1;
  int fjlen = fjhi-fjlo+1;
  int fklen = fkhi-fklo+1;
  Real *fptr = fine.dataPtr();
  int tilo = tbox.loVect()[0];
  int tjlo = tbox.loVect()[1];
  int tklo = tbox.loVect()[2];
  int tihi = tbox.hiVect()[0];
  int tjhi = tbox.hiVect()[1];
  int tkhi = tbox.hiVect()[2];
  // run over coarse grid points in tbox and stuff it in fine
  for( int n=0; n<nc; n++){
    for( k=tklo; k<=tkhi; k++){
      for( j=tjlo; j<=tjhi; j++){
	for( i=tilo; i<=tihi; i++){
	  int coffset = i-cilo+cilen*(j-cjlo+cjlen*(k-cklo+cklen*n));
	  int fioff = i*ratio+ioff;
	  int fjoff = j*ratio+joff;
	  int fkoff = k*ratio+koff;
	  if( fioff >= filo && fioff <= fihi &&
	      fjoff >= fjlo && fjoff <= fjhi &&
	      fkoff >= fklo && fkoff <= fkhi ) {
	    int foffset = fioff-filo+filen*
	      (fjoff-fjlo+fjlen*(fkoff-fklo+fklen*n));
	    fptr[foffset] = cptr[coffset];
	  }
	}
      }
    }
  }
}
#else
whoa! this cant happen
#endif

void
dummyFunc()
{
  std::cerr << "the impossible has happened: dummyFunc has been called!"<<std::endl;
  exit(1);
}

// compute complement of two BoxArrays
// result will be complement of ba2 in ba1, i.e. ba1 minus intersection 
BoxArray *
complementIn(const BoxArray& ba1, const BoxArray& ba2)
{
  BoxList res;
  BoxList bl2(ba2);
  int len = ba1.size();
  for( int i=0; i<len; i++){
    BoxList tmp;
    tmp.complementIn(ba1[i],bl2);
    res.join(tmp);
  }
  res.minimize();
  return new BoxArray(res);
}

// compute union of two BoxArray's (output will be disjoint boxes)
BoxArray *
join(const BoxArray &ba1, const BoxArray &ba2)
{
  BoxDomain res;
  res.add(BoxList(ba1));
  res.add(BoxList(ba2));
  res.simplify();
  return new BoxArray(res.boxList());
}

#if 0
// recent version of boxlib doesn't allow this
// compute intersection of two BoxArray's
BoxArray *
intersect( const BoxArray &ba1, const BoxArray &ba2)
{
  BoxList accum;
  BoxList bl2(ba2);
  int len1 = ba1.size();
  for( int i=0; i<len1; i++){
    const Box &bx = ba1[i];
    accum.join(intersect(bl2,bx));
  }
  accum.simplify();
  return new BoxArray(accum);
}
#endif


// create block of memory containing portable gray map representation of
// FAB

int _charArrayLen;
//#include <malloc.h>

#if (BL_SPACEDIM==2)
char *
makePGM( FArrayBox* _fab, int icomp, double usemin, double usemax)
{
  FArrayBox &fab = *_fab;
  char tmp[512];  // use for temp storage of header
  Box bx = fab.box();
  int width = bx.length(0);
  int height = bx.length(1);
  int npts = bx.numPts();
  sprintf(tmp,"P5\n%d\n%d\n255\n",width,height);
  int len = std::strlen(tmp);
  // now allocate the big buffer
  //char *ptr = new char[len+npts];
  char *ptr = (char *)malloc(len+npts);
  strcpy(ptr,tmp);
  char *nextchar = ptr+len;
  Real *nextreal = fab.dataPtr(icomp);

  int i,j;
  int ilo = bx.loVect()[0];
  int jlo = bx.loVect()[1];
  int ihi = bx.hiVect()[0];
  int jhi = bx.hiVect()[1];
  int ilen = (ihi-ilo+1);

  // this reverses j direction to agree with graphics idiom
  for( j=jlo; j<=jhi; ++j){
    for(i=ilo; i<=ihi; ++i){
      int val = (int)
	((nextreal[i-ilo+(jhi-j)*ilen]-usemin)/(usemax-usemin)*255);
      if( val < 0 ) val = 0;
      if( val > 255 ) val = 255;
      nextchar[i-ilo+(jhi-j)*ilen] = (char)val;
    }
  }

  _charArrayLen = len+npts;
  return ptr;
}
#elif (BL_SPACEDIM==3)
char *
makePGM( FArrayBox* _fab, int icomp, double usemin, double usemax)
{
  FArrayBox &fab = *_fab;
  char tmp[512];  // use for temp storage of header
  Box bx = fab.box();
  int width,height;
  if( bx.length(0)==1 ){
    width = bx.length(1);
    height = bx.length(2);
  }
  else if(bx.length(1)==1){
    width = bx.length(0);
    height = bx.length(2);
  }
  else{
    width = bx.length(0);
    height = bx.length(1);
  }
  int npts = width*height;
  sprintf(tmp,"P5\n%d\n%d\n255\n",width,height);
  int len = strlen(tmp);
  // now allocate the big buffer
  char *ptr = (char *)malloc(len+npts);
  strcpy(ptr,tmp);
  char *nextchar = ptr+len;
  Real *nextreal = fab.dataPtr(icomp);

  int i,j,k;
  int ilo = bx.loVect()[0];
  int jlo = bx.loVect()[1];
  int klo = bx.loVect()[2];
  int ihi = bx.hiVect()[0];
  int jhi = bx.hiVect()[1];
  int khi = bx.hiVect()[2];
  int ilen = (ihi-ilo+1);
  int jlen = (jhi-jlo+1);

  if( bx.length(0)==1 ){
    for( k=klo; k<=khi; ++k){
      for(j=jlo; j<=jhi; ++j){
	int val = (int)
	  ((nextreal[j-jlo+(khi-k)*jlen]-usemin)/(usemax-usemin)*255);
	if( val < 0 ) val = 0;
	if( val > 255 ) val = 255;
	nextchar[j-jlo+(khi-k)*jlen] = (char)val;
      }
    }
  }
  else if( bx.length(1)==1 ){
    for( k=klo; k<=khi; ++k){
      for(i=ilo; i<=ihi; ++i){
	int val = (int)
	  ((nextreal[i-ilo+(khi-k)*ilen]-usemin)/(usemax-usemin)*255);
	if( val < 0 ) val = 0;
	if( val > 255 ) val = 255;
	nextchar[i-ilo+(khi-k)*ilen] = (char)val;
      }
    }
  }
  else {
    for( j=jlo; j<=jhi; ++j){
      for(i=ilo; i<=ihi; ++i){
	int val = (int)
	  ((nextreal[i-ilo+(jhi-j)*ilen]-usemin)/(usemax-usemin)*255);
	if( val < 0 ) val = 0;
	if( val > 255 ) val = 255;
	nextchar[i-ilo+(jhi-j)*ilen] = (char)val;
      }
    }
  }

  _charArrayLen = len+npts;
  return ptr;
}
#else
whoa! this cant happen
#endif


// this is to fix up some problems with tkintr and swig
extern "C" {
  char * getprogramname();
  char * Py_GetProgramName();
  char * PyString_FromStringAndSize();
};

static void dummy_stuff_()
{
  PyString_FromStringAndSize();
}

#if 0
char *
Py_GetProgramName()
{
  return getprogramname();
}
#endif
