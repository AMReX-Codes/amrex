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

TagBox::TagBox (Arena* ar) noexcept
    : BaseFab<TagBox::TagType>(ar)
{}

TagBox::TagBox (const Box& bx, int n, Arena* ar)
    : BaseFab<TagBox::TagType>(bx,n,ar)
{
    setVal<RunOn::Host>(TagBox::CLEAR);
}

TagBox::TagBox (const Box& bx, int n, bool alloc, bool shared, Arena* ar)
    : BaseFab<TagBox::TagType>(bx,n,alloc,shared,ar)
{
    if (alloc) setVal<RunOn::Host>(TagBox::CLEAR);
}

TagBox::TagBox (const TagBox& rhs, MakeType make_type, int scomp, int ncomp)
    : BaseFab<TagBox::TagType>(rhs,make_type,scomp,ncomp)
{}

void
TagBox::coarsen (const IntVect& ratio, const Box& cbox) noexcept
{
    // xxxxx TODO: gpu

    BL_ASSERT(nComp() == 1);
    Array4<char const> const& farr = this->const_array();

    TagBox cfab(cbox, 1, The_Cpu_Arena());
    Array4<char> const& carr = cfab.array();

    const auto flo = amrex::lbound(domain);
    const auto fhi = amrex::ubound(domain);
    Dim3 r{1,1,1};
    AMREX_D_TERM(r.x = ratio[0];, r.y = ratio[1];, r.z = ratio[2]);

    for (int k = flo.z; k <= fhi.z; ++k) {
        int kc = amrex::coarsen(k,r.z);
        for (int j = flo.y; j <= fhi.y; ++j) {
            int jc = amrex::coarsen(j,r.y);
            for (int i = flo.x; i <= fhi.x; ++i) {
                int ic = amrex::coarsen(i,r.x);
                carr(ic,jc,kc) = carr(ic,jc,kc) || farr(i,j,k);
            }
        }
    }

    std::memcpy(this->dataPtr(), cfab.dataPtr(), sizeof(TagType)*cbox.numPts());
    this->domain = cbox;
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

Long
TagBox::numTags () const noexcept
{
    Long nt = 0L;
    Long len = domain.numPts();
    const TagType* d = dataPtr();
    for (Long n = 0; n < len; ++n)
    {
        if (d[n] != TagBox::CLEAR)
            ++nt;
    }
    return nt;
}

Long
TagBox::numTags (const Box& b) const noexcept
{
   TagBox tempTagBox(b,1);
   tempTagBox.copy<RunOn::Host>(*this);
   return tempTagBox.numTags();
}

Long
TagBox::collate (Vector<IntVect>& ar, int start) const noexcept
{
    BL_ASSERT(start >= 0);
    //
    // Starting at given offset of array ar, enter location (IntVect) of
    // each tagged cell in tagbox.
    //
    Long count       = 0;
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
    
    Long stride[] = {1, Lbx[0], Long(Lbx[0])*Long(Lbx[1])};

    Long Ntb = 1, stb=0;
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
    
    Long stride[] = {1, Lbx[0], Long(Lbx[0])*Long(Lbx[1])};

    Long stb=0;
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
                if (*iptr) *cptr = *iptr;
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
    
    Long stride[] = {1, Lbx[0], Long(Lbx[0])*Long(Lbx[1])};

    Long stb=0;
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

void 
TagBoxArray::buffer (const IntVect& nbuf)
{
    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

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
TagBoxArray::mapPeriodicRemoveDuplicates (const Geometry& geom)
{
    BL_PROFILE("TagBoxArray::mapPRD");

    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

    TagBoxArray tmp(boxArray(),DistributionMap(),0); // note that tmp is filled w/ CLEAR.

    tmp.ParallelAdd(*this, 0, 0, 1, nGrowVect(), IntVect{0}, geom.periodicity());

    std::swap(*this, tmp);
}

Long
TagBoxArray::numTags () const
{
    Long ntag = 0;

    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

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

    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

    Long count = 0;

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

    //
    // The total number of tags system wide that must be collated.
    //
    Long numtags = count;
    ParallelDescriptor::ReduceLongSum(numtags);

    if (numtags == 0) {
        TheGlobalCollateSpace.clear();
        return;
    } else if (numtags > static_cast<Long>(std::numeric_limits<int>::max())) {
        // xxxxx todo
        amrex::Abort("TagBoxArray::collate: Too many tags. Using a larger blocking factor might help. Please file an issue on github");
    }

#ifdef BL_USE_MPI
    //
    // On I/O proc. this holds all tags after they've been gather'd.
    // On other procs. non-mempty signals size is not zero.
    //
    if (ParallelDescriptor::IOProcessor()) {
        TheGlobalCollateSpace.resize(numtags);
    } else {
        TheGlobalCollateSpace.resize(1);
    }

    //
    // Tell root CPU how many tags each CPU will be sending.
    //
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    const std::vector<int>& countvec = ParallelDescriptor::Gather(static_cast<int>(count),
                                                                  IOProcNumber);
    std::vector<int> offset(countvec.size(),0);
    if (ParallelDescriptor::IOProcessor()) {
        for (int i = 1, N = offset.size(); i < N; i++) {
	    offset[i] = offset[i-1] + countvec[i-1];
	}
    }
    //
    // Gather all the tags to IOProcNumber into TheGlobalCollateSpace.
    //
    const IntVect* psend = (count > 0) ? TheLocalCollateSpace.data() : nullptr;
    IntVect* precv = TheGlobalCollateSpace.data();
    ParallelDescriptor::Gatherv(psend, count, precv, countvec, offset, IOProcNumber);

#else
    TheGlobalCollateSpace = std::move(TheLocalCollateSpace);
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
    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

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
            tags.setVal<RunOn::Host>(val,isects[i].second,0);
        }
    }
}

void
TagBoxArray::coarsen (const IntVect & ratio)
{
    // If team is used, all team workers need to go through all the fabs, including ones they don't own.
    int teamsize = ParallelDescriptor::TeamSize();
    unsigned char flags = (teamsize == 1) ? 0 : MFIter::AllBoxes;

    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

    IntVect new_n_grow;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        new_n_grow[idim] = (n_grow[idim]+ratio[idim]-1)/ratio[idim];
    }

#if defined(_OPENMP)
#pragma omp parallel if (teamsize == 1)
#endif
    for (MFIter mfi(*this,flags); mfi.isValid(); ++mfi)
    {
        Box const& cbox = amrex::grow(amrex::coarsen(mfi.validbox(),ratio),new_n_grow);
        this->fabPtr(mfi)->coarsen(ratio,cbox);
    }

    boxarray.coarsen(ratio);
    n_grow = new_n_grow;
}

}

