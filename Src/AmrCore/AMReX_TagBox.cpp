#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <climits>

#include <AMReX_TagBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ccse-mpi.H>

namespace amrex {

TagBox::TagBox () noexcept {}

TagBox::TagBox (const Box& bx,
                int        n,
                bool       alloc,
		bool       shared)
    :
    BaseFab<TagBox::TagType>(bx,n,alloc,shared)
{
    if (alloc) setVal(TagBox::CLEAR);
}

TagBox::TagBox (const TagBox& rhs, MakeType make_type, int scomp, int ncomp)
    :
    BaseFab<TagBox::TagType>(rhs,make_type,scomp,ncomp)
{
}

void
TagBox::coarsen (const IntVect& ratio, bool owner) noexcept
{
    BL_ASSERT(nComp() == 1);

    TagType*   fdat     = dataPtr();
    IntVect    lov      = domain.smallEnd();
    IntVect    hiv      = domain.bigEnd();
    IntVect    d_length = domain.size();
    const int* flo      = lov.getVect();
    const int* fhi      = hiv.getVect();
    const int* flen     = d_length.getVect();

    const Box& cbox = amrex::coarsen(domain,ratio);

    this->nvar = 1;
    this->domain = cbox;

    if (!owner) {
        return;
    }

    const int* clo      = cbox.loVect();
    IntVect    cbox_len = cbox.size();
    const int* clen     = cbox_len.getVect();

    Box b1(amrex::refine(cbox,ratio));
    const int* lo       = b1.loVect();
    int        longlen  = b1.longside();

    long numpts = domain.numPts();
    Vector<TagType> cfab(numpts);
    TagType* cdat = cfab.dataPtr();

    Vector<TagType> t(longlen,TagBox::CLEAR);

    int klo = 0, khi = 0, jlo = 0, jhi = 0, ilo, ihi;
    AMREX_D_TERM(ilo=flo[0]; ihi=fhi[0]; ,
           jlo=flo[1]; jhi=fhi[1]; ,
           klo=flo[2]; khi=fhi[2];)

#define IXPROJ(i,r) (((i)+(r)*std::abs(i))/(r) - std::abs(i))
#define IOFF(j,k,lo,len) AMREX_D_TERM(0, +(j-lo[1])*len[0], +(k-lo[2])*len[0]*len[1])
   
   int ratiox = 1, ratioy = 1, ratioz = 1;
   AMREX_D_TERM(ratiox = ratio[0];,
          ratioy = ratio[1];,
          ratioz = ratio[2];)

   for (int k = klo; k <= khi; k++)
   {
       const int kc = IXPROJ(k,ratioz);
       for (int j = jlo; j <= jhi; j++)
       {
           const int     jc = IXPROJ(j,ratioy);
           TagType*       c = cdat + IOFF(jc,kc,clo,clen);
           const TagType* f = fdat + IOFF(j,k,flo,flen);
           //
           // Copy fine grid row of values into tmp array.
           //
           for (int i = ilo; i <= ihi; i++)
               t[i-lo[0]] = f[i-ilo];

           for (int off = 0; off < ratiox; off++)
           {
               for (int ic = 0; ic < clen[0]; ic++)
               {
                   const int i = ic*ratiox + off;
                   c[ic] = std::max(c[ic],t[i]);
               }
           }
       }
   }

#undef IXPROJ
#undef IOFF

   for (int i = 0; i < numpts; ++i) {
       fdat[i] = cdat[i];
   }
}

void 
TagBox::buffer (const IntVect& nbuff, const IntVect& nwid) noexcept
{
    //
    // Note: this routine assumes cell with TagBox::SET tag are in
    // interior of tagbox (region = grow(domain,-nwid)).
    //
    Box inside(domain);
    inside.grow(-nwid);
    const int* inlo = inside.loVect();
    const int* inhi = inside.hiVect();

    int klo = 0, khi = 0, jlo = 0, jhi = 0, ilo, ihi;
    AMREX_D_TERM(ilo=inlo[0]; ihi=inhi[0]; ,
                 jlo=inlo[1]; jhi=inhi[1]; ,
                 klo=inlo[2]; khi=inhi[2];)

    int ni = 0, nj = 0, nk = 0;
    AMREX_D_TERM(ni=nbuff[0];, nj=nbuff[1];, nk=nbuff[2];)

    IntVect d_length = domain.size();
    const int* len = d_length.getVect();
    const int* lo = domain.loVect();
    TagType* d = dataPtr();

#define OFF(i,j,k,lo,len) AMREX_D_TERM(i-lo[0], +(j-lo[1])*len[0] , +(k-lo[2])*len[0]*len[1])
   
    for (int k = klo; k <= khi; k++)
    {
        for (int j = jlo; j <= jhi; j++)
        {
            for (int i = ilo; i <= ihi; i++)
            {
                TagType* d_check = d + OFF(i,j,k,lo,len);
                if (*d_check == TagBox::SET)
                {
                    for (int kk = -nk; kk <= nk; kk++)
                    {
                        for (int jj = -nj; jj <= nj; jj++)
                        {
                            for (int ii = -ni; ii <= ni; ii++)
                            {
                                TagType* dn = d_check+ AMREX_D_TERM(ii, +jj*len[0], +kk*len[0]*len[1]);
                                if (*dn !=TagBox::SET)
                                    *dn = TagBox::BUF;
                            }
                        }
                    }
                }
            }
        }
    }
#undef OFF
}

void 
TagBox::merge (const TagBox& src) noexcept
{
    //
    // Compute intersections.
    //
    const Box& bx = domain & src.domain;

    if (bx.ok())
    {
        const int*     dlo        = domain.loVect();
        IntVect        d_length   = domain.size();
        const int*     dleng      = d_length.getVect();
        const int*     slo        = src.domain.loVect();
        IntVect        src_length = src.domain.size();
        const int*     sleng      = src_length.getVect();
        const int*     lo         = bx.loVect();
        const int*     hi         = bx.hiVect();
        const TagType* ds0        = src.dataPtr();
        TagType*       dd0        = dataPtr();

        int klo = 0, khi = 0, jlo = 0, jhi = 0, ilo, ihi;
        AMREX_D_TERM(ilo=lo[0]; ihi=hi[0]; ,
               jlo=lo[1]; jhi=hi[1]; ,
               klo=lo[2]; khi=hi[2];)

#define OFF(i,j,k,lo,len) AMREX_D_TERM(i-lo[0], +(j-lo[1])*len[0] , +(k-lo[2])*len[0]*len[1])
      
        for (int k = klo; k <= khi; k++)
        {
            for (int j = jlo; j <= jhi; j++)
            {
                for (int i = ilo; i <= ihi; i++)
                {
                    const TagType* ds = ds0 + OFF(i,j,k,slo,sleng);
                    if (*ds != TagBox::CLEAR)
                    {
                        TagType* dd = dd0 + OFF(i,j,k,dlo,dleng);
                        *dd = TagBox::SET;
                    }            
                }
            }
        }
    }
#undef OFF
}

long
TagBox::numTags () const noexcept
{
    long nt = 0L;
    long len = domain.numPts();
    const TagType* d = dataPtr();
    for (long n = 0; n < len; ++n)
    {
	if (d[n] != TagBox::CLEAR)
	    ++nt;
    }
    return nt;
}

long
TagBox::numTags (const Box& b) const noexcept
{
   TagBox tempTagBox(b,1);
   tempTagBox.copy(*this);
   return tempTagBox.numTags();
}

long
TagBox::collate (Vector<IntVect>& ar, int start) const noexcept
{
    BL_ASSERT(start >= 0);
    //
    // Starting at given offset of array ar, enter location (IntVect) of
    // each tagged cell in tagbox.
    //
    long count       = 0;
    IntVect d_length = domain.size();
    const int* len   = d_length.getVect();
    const int* lo    = domain.loVect();
    const TagType* d = dataPtr();
    int ni = 1, nj = 1, nk = 1;
    AMREX_D_TERM(ni = len[0]; , nj = len[1]; , nk = len[2];)

    for (int k = 0; k < nk; k++)
    {
        for (int j = 0; j < nj; j++)
        {
            for (int i = 0; i < ni; i++)
            {
                const TagType* dn = d + AMREX_D_TERM(i, +j*len[0], +k*len[0]*len[1]);
                if (*dn != TagBox::CLEAR)
                {
                    ar[start++] = IntVect(AMREX_D_DECL(lo[0]+i,lo[1]+j,lo[2]+k));
                    count++;
                }
            }
        }
    }
    return count;
}

Vector<int>
TagBox::tags () const noexcept
{
    Vector<int> ar(domain.numPts(), TagBox::CLEAR);

    const TagType* cptr = dataPtr();
    int*           iptr = ar.dataPtr();

    for (int i = 0; i < ar.size(); i++, cptr++, iptr++)
    {
        if (*cptr)
            *iptr = *cptr;
    }

    return ar;
}


// Set values as specified by the array -- this only tags.
// It's an error if ar.length() != domain.numPts().
void
TagBox::tags (const Vector<int>& ar) noexcept
{
    BL_ASSERT(ar.size() == domain.numPts());

    TagType*   cptr = dataPtr();
    const int* iptr = ar.dataPtr();

    for (int i = 0; i < ar.size(); i++, cptr++, iptr++)
    {
        if (*iptr)
            *cptr = *iptr;
    }
}

// Set values as specified by the array -- this tags and untags.
// It's an error if ar.length() != domain.numPts().
void
TagBox::tags_and_untags (const Vector<int>& ar) noexcept
{
    BL_ASSERT(ar.size() == domain.numPts());

    TagType*   cptr = dataPtr();
    const int* iptr = ar.dataPtr();

    // This clears as well as sets tags.
    for (int i = 0; i < ar.size(); i++, cptr++, iptr++)
    {
        *cptr = *iptr;
    }
}

// Since a TagBox is a BaseFab<char>, we can use this utility
// function to allocate an integer array to have the same number 
// of elements as cells in tilebx
void 
TagBox::get_itags(Vector<int>& ar, const Box& tilebx) const noexcept
{
    auto dlen = length();
    int Lbx[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
	Lbx[idim] = dlen[idim];
    }
    
    long stride[] = {1, Lbx[0], long(Lbx[0])*long(Lbx[1])};

    long Ntb = 1, stb=0;
    int Ltb[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
	Ltb[idim] = tilebx.length(idim);
	Ntb *= Ltb[idim];
	stb += stride[idim] * (tilebx.smallEnd(idim) - domain.smallEnd(idim));
    }
    
    if (ar.size() < Ntb) ar.resize(Ntb);
    
    const TagType* const p0   = dataPtr() + stb;  // +stb to the lower corner of tilebox
    int*                 iptr = ar.dataPtr();

    for (int k=0; k<Ltb[2]; k++) {
        for (int j=0; j<Ltb[1]; j++) {
	    const TagType* cptr = p0 + j*stride[1] + k*stride[2];
	    for (int i=0; i<Ltb[0]; i++, cptr++, iptr++) {
		if (*cptr) {
		    *iptr = *cptr;
		}
		else {
		    *iptr = TagBox::CLEAR;
		}
	    }
        }
    }
}

// Set values as specified by the array -- this only tags.
// only changes values in the tilebx region
void 
TagBox::tags (const Vector<int>& ar, const Box& tilebx) noexcept
{
    auto dlen = length();
    int Lbx[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
	Lbx[idim] = dlen[idim];
    }
    
    long stride[] = {1, Lbx[0], long(Lbx[0])*long(Lbx[1])};

    long stb=0;
    int Ltb[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
	Ltb[idim] = tilebx.length(idim);
	stb += stride[idim] * (tilebx.smallEnd(idim) - domain.smallEnd(idim));
    }
    
    TagType* const p0   = dataPtr() + stb;  // +stb to the lower corner of tilebox
    const int*     iptr = ar.dataPtr();

    for (int k=0; k<Ltb[2]; k++) {
        for (int j=0; j<Ltb[1]; j++) {
	    TagType* cptr = p0 + j*stride[1] + k*stride[2];
	    for (int i=0; i<Ltb[0]; i++, cptr++, iptr++) {
		if (*iptr) 
		    *cptr = *iptr;
	    }
        }
    }
}

// Set values as specified by the array -- this tags and untags.
// only changes values in the tilebx region
void 
TagBox::tags_and_untags (const Vector<int>& ar, const Box& tilebx) noexcept
{
    auto dlen = length();
    int Lbx[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
	Lbx[idim] = dlen[idim];
    }
    
    long stride[] = {1, Lbx[0], long(Lbx[0])*long(Lbx[1])};

    long stb=0;
    int Ltb[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
	Ltb[idim] = tilebx.length(idim);
	stb += stride[idim] * (tilebx.smallEnd(idim) - domain.smallEnd(idim));
    }
    
    TagType* const p0   = dataPtr() + stb;  // +stb to the lower corner of tilebox
    const int*     iptr = ar.dataPtr();

    for (int k=0; k<Ltb[2]; k++) {
        for (int j=0; j<Ltb[1]; j++) {
	    TagType* cptr = p0 + j*stride[1] + k*stride[2];
	    for (int i=0; i<Ltb[0]; i++, cptr++, iptr++) {
		*cptr = *iptr;
	    }
        }
    }
}

TagBoxArray::TagBoxArray (const BoxArray& ba,
			  const DistributionMapping& dm,
                          int             _ngrow)
    :
    FabArray<TagBox>(ba,dm,1,_ngrow,MFInfo(),DefaultFabFactory<TagBox>())
{
    if (SharedMemory()) setVal(TagBox::CLEAR);
}

TagBoxArray::TagBoxArray (const BoxArray& ba,
			  const DistributionMapping& dm,
                          const IntVect&  _ngrow)
    :
    FabArray<TagBox>(ba,dm,1,_ngrow,MFInfo(),DefaultFabFactory<TagBox>())
{
    if (SharedMemory()) setVal(TagBox::CLEAR);
}

IntVect
TagBoxArray::borderSize () const noexcept
{
    return n_grow;
}

void 
TagBoxArray::buffer (const IntVect& nbuf)
{
    AMREX_ASSERT(nbuf.allLE(n_grow));

    if (nbuf.max() > 0)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
       for (MFIter mfi(*this); mfi.isValid(); ++mfi)
           get(mfi).buffer(nbuf, n_grow);
    }
}

void
TagBoxArray::mapPeriodic (const Geometry& geom)
{
    if (!geom.isAnyPeriodic()) return;

    BL_PROFILE("TagBoxArray::mapPeriodic()");

    // This function is called after coarsening.
    // So we can assume that n_grow is 0.
    BL_ASSERT(n_grow[0] == 0);

    TagBoxArray tmp(boxArray(),DistributionMap()); // note that tmp is filled w/ CLEAR.

    tmp.copy(*this, geom.periodicity(), FabArrayBase::ADD);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
	get(mfi).merge(tmp[mfi]);
    }
}

long
TagBoxArray::numTags () const
{
    long ntag = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:ntag)
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
	ntag += get(mfi).numTags();
    }
    
    ParallelDescriptor::ReduceLongSum(ntag);
    
    return ntag;
}

void
TagBoxArray::collate (Vector<IntVect>& TheGlobalCollateSpace) const
{
    BL_PROFILE("TagBoxArray::collate()");

    long count = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:count)
#endif
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        count += get(fai).numTags();
    }

    //
    // Local space for holding just those tags we want to gather to the root cpu.
    //
    Vector<IntVect> TheLocalCollateSpace(count);

    count = 0;

    // unsafe to do OMP
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        count += get(fai).collate(TheLocalCollateSpace,count);
    }

    if (count > 0)
    {
        amrex::RemoveDuplicates(TheLocalCollateSpace);
	count = TheLocalCollateSpace.size();
    }
    //
    // The total number of tags system wide that must be collated.
    // This is really just an estimate of the upper bound due to duplicates.
    // While we've removed duplicates per MPI process there's still more systemwide.
    //
    long numtags = count;

    ParallelDescriptor::ReduceLongSum(numtags);

    if (numtags == 0) {
	TheGlobalCollateSpace.clear();
	return;
    }

    //
    // This holds all tags after they've been gather'd and unique'ified.
    //
    // Each CPU needs an identical copy since they all must go through grid_places() which isn't parallelized.

    TheGlobalCollateSpace.resize(numtags);

#ifdef BL_USE_MPI
    //
    // Tell root CPU how many tags each CPU will be sending.
    //
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    count *= AMREX_SPACEDIM;  // Convert from count of tags to count of integers to expect.
    const std::vector<long>& countvec = ParallelDescriptor::Gather(count, IOProcNumber);
    
    std::vector<long> offset(countvec.size(),0L);
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 1, N = offset.size(); i < N; i++) {
	    offset[i] = offset[i-1] + countvec[i-1];
	}
    }
    //
    // Gather all the tags to IOProcNumber into TheGlobalCollateSpace.
    //
    BL_ASSERT(sizeof(IntVect) == AMREX_SPACEDIM * sizeof(int));
    const int* psend = (count > 0) ? TheLocalCollateSpace[0].getVect() : 0;
    int* precv = TheGlobalCollateSpace[0].getVect();
    ParallelDescriptor::Gatherv(psend, count,
				precv, countvec, offset, IOProcNumber); 

    if (ParallelDescriptor::IOProcessor())
    {
        amrex::RemoveDuplicates(TheGlobalCollateSpace);
	numtags = TheGlobalCollateSpace.size();
    }

    //
    // Now broadcast them back to the other processors.
    //
    ParallelDescriptor::Bcast(&numtags, 1, IOProcNumber);
    ParallelDescriptor::Bcast(TheGlobalCollateSpace[0].getVect(), numtags*AMREX_SPACEDIM, IOProcNumber);
    TheGlobalCollateSpace.resize(numtags);

#else
    //
    // Copy TheLocalCollateSpace to TheGlobalCollateSpace.
    //
    TheGlobalCollateSpace = TheLocalCollateSpace;
#endif
}

void
TagBoxArray::setVal (const BoxList& bl,
                     TagBox::TagVal val)
{
    BoxArray ba(bl);
    setVal(ba,val);
}

void
TagBoxArray::setVal (const BoxDomain& bd,
                     TagBox::TagVal   val)
{
    setVal(bd.boxList(),val);
}

void
TagBoxArray::setVal (const BoxArray& ba,
                     TagBox::TagVal  val)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
	std::vector< std::pair<int,Box> > isects;

        ba.intersections(mfi.fabbox(),isects);

        TagBox& tags = get(mfi);

        for (int i = 0, N = isects.size(); i < N; i++)
        {
            tags.setVal(val,isects[i].second,0);
        }
    }
}

void
TagBoxArray::coarsen (const IntVect & ratio)
{
    // If team is used, all team workers need to go through all the fabs, including ones they don't own.
    int teamsize = ParallelDescriptor::TeamSize();
    unsigned char flags = (teamsize == 1) ? 0 : MFIter::AllBoxes;

#if defined(_OPENMP)
#pragma omp parallel if (teamsize == 1)
#endif
    for (MFIter mfi(*this,flags); mfi.isValid(); ++mfi)
    {
        this->fabHostPtr(mfi)->coarsen(ratio,isOwner(mfi.LocalIndex()));
#ifdef AMREX_USE_GPU
        this->fabDevicePtr(mfi)->coarsen(ratio,false);
#endif
    }

    boxarray.growcoarsen(n_grow,ratio);
    updateBDKey();  // because we just modify boxarray in-place.

    n_grow = IntVect::TheZeroVector();
}

}

