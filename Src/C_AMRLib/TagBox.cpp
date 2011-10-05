
#include <winstd.H>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <climits>

#include <TagBox.H>
#include <Geometry.H>
#include <ParallelDescriptor.H>
#include <ccse-mpi.H>

TagBox::TagBox () {}

TagBox::TagBox (const Box& bx,
                int        n)
    :
    BaseFab<TagBox::TagType>(bx,n)
{
    setVal(TagBox::CLEAR);
}

void
TagBox::resize (const Box& b, int ncomp)
{
    BaseFab<TagType>::resize(b,ncomp);
}

TagBox*
TagBox::coarsen (const IntVect& ratio)
{
    TagBox* crse = new TagBox(BoxLib::coarsen(domain,ratio));

    const Box& cbox = crse->box();

    Box b1(BoxLib::refine(cbox,ratio));

    const int* flo      = domain.loVect();
    const int* fhi      = domain.hiVect();
    IntVect    d_length = domain.size();
    const int* flen     = d_length.getVect();
    const int* clo      = cbox.loVect();
    IntVect    cbox_len = cbox.size();
    const int* clen     = cbox_len.getVect();
    const int* lo       = b1.loVect();
    int        longlen  = b1.longside();

    TagType* fdat = dataPtr();
    TagType* cdat = crse->dataPtr();

    Array<TagType> t(longlen,TagBox::CLEAR);

    int klo = 0, khi = 0, jlo = 0, jhi = 0, ilo, ihi;
    D_TERM(ilo=flo[0]; ihi=fhi[0]; ,
           jlo=flo[1]; jhi=fhi[1]; ,
           klo=flo[2]; khi=fhi[2];)

#define IXPROJ(i,r) (((i)+(r)*std::abs(i))/(r) - std::abs(i))
#define IOFF(j,k,lo,len) D_TERM(0, +(j-lo[1])*len[0], +(k-lo[2])*len[0]*len[1])
   
   int ratiox = 1, ratioy = 1, ratioz = 1;
   D_TERM(ratiox = ratio[0];,
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

   return crse;

#undef IXPROJ
#undef IOFF
}

void 
TagBox::buffer (int nbuff,
                int nwid)
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
    D_TERM(ilo=inlo[0]; ihi=inhi[0]; ,
           jlo=inlo[1]; jhi=inhi[1]; ,
           klo=inlo[2]; khi=inhi[2];)

    IntVect d_length = domain.size();
    const int* len = d_length.getVect();
    const int* lo = domain.loVect();
    int ni = 0, nj = 0, nk = 0;
    D_TERM(ni, =nj, =nk) = nbuff;
    TagType* d = dataPtr();

#define OFF(i,j,k,lo,len) D_TERM(i-lo[0], +(j-lo[1])*len[0] , +(k-lo[2])*len[0]*len[1])
   
    for (int k = klo; k <= khi; k++)
    {
        for (int j = jlo; j <= jhi; j++)
        {
            for (int i = ilo; i <= ihi; i++)
            {
                TagType* d_check = d + OFF(i,j,k,lo,len);
                if (*d_check == TagBox::SET)
                {
                    for (int k = -nk; k <= nk; k++)
                    {
                        for (int j = -nj; j <= nj; j++)
                        {
                            for (int i = -ni; i <= ni; i++)
                            {
                                TagType* dn = d_check+ D_TERM(i, +j*len[0], +k*len[0]*len[1]);
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
TagBox::merge (const TagBox& src)
{
    //
    // Compute intersections.
    //
    Box bx = domain & src.domain;

    if (bx.ok())
    {
        const int* dlo  = domain.loVect();
        IntVect d_length = domain.size();
        const int* dlen = d_length.getVect();
        const int* slo  = src.domain.loVect();
        IntVect src_length = src.domain.size();
        const int* slen = src_length.getVect();
        const int* lo   = bx.loVect();
        const int* hi   = bx.hiVect();

        const TagType* ds0 = src.dataPtr();
        TagType* dd0       = dataPtr();

        int klo = 0, khi = 0, jlo = 0, jhi = 0, ilo, ihi;
        D_TERM(ilo=lo[0]; ihi=hi[0]; ,
               jlo=lo[1]; jhi=hi[1]; ,
               klo=lo[2]; khi=hi[2];)

#define OFF(i,j,k,lo,len) D_TERM(i-lo[0], +(j-lo[1])*len[0] , +(k-lo[2])*len[0]*len[1])
      
        for (int k = klo; k <= khi; k++)
        {
            for (int j = jlo; j <= jhi; j++)
            {
                for (int i = ilo; i <= ihi; i++)
                {
                    const TagType* ds = ds0 + OFF(i,j,k,slo,slen);
                    if (*ds != TagBox::CLEAR)
                    {
                        TagType* dd = dd0 + OFF(i,j,k,dlo,dlen);
                        *dd = TagBox::SET;
                    }            
                }
            }
        }
    }
#undef OFF
}

int
TagBox::numTags () const
{
   int nt = 0;
   long t_long = domain.numPts();
   BL_ASSERT(t_long < INT_MAX);
   int len = static_cast<int>(t_long);
   const TagType* d = dataPtr();
   for (int n = 0; n < len; n++)
   {
      if (d[n] != TagBox::CLEAR)
          nt++;
   }
   return nt;
}

int
TagBox::numTags (const Box& b) const
{
   TagBox tempTagBox(b,1);
   tempTagBox.copy(*this);
   return tempTagBox.numTags();
}

int
TagBox::collate (IntVect* ar,
                 int      start) const
{
    BL_ASSERT(!(ar == 0));
    BL_ASSERT(start >= 0);
    //
    // Starting at given offset of array ar, enter location (IntVect) of
    // each tagged cell in tagbox.
    //
    int count        = 0;
    IntVect d_length = domain.size();
    const int* len   = d_length.getVect();
    const int* lo    = domain.loVect();
    const TagType* d = dataPtr();
    int ni = 1, nj = 1, nk = 1;
    D_TERM(ni = len[0]; , nj = len[1]; , nk = len[2];)

    for (int k = 0; k < nk; k++)
    {
        for (int j = 0; j < nj; j++)
        {
            for (int i = 0; i < ni; i++)
            {
                const TagType* dn = d + D_TERM(i, +j*len[0], +k*len[0]*len[1]);
                if (*dn != TagBox::CLEAR)
                {
                    ar[start++] = IntVect(D_DECL(lo[0]+i,lo[1]+j,lo[2]+k));
                    count++;
                }
            }
        }
    }
    return count;
}

Array<int>
TagBox::tags () const
{
    Array<int> ar(domain.numPts(), TagBox::CLEAR);

    const TagType* cptr = dataPtr();
    int*           iptr = ar.dataPtr();

    for (int i = 0; i < ar.size(); i++, cptr++, iptr++)
    {
        if (*cptr)
            *iptr = *cptr;
    }

    return ar;
}

void
TagBox::tags (const Array<int>& ar)
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

void
TagBox::tags_and_untags (const Array<int>& ar)
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

TagBoxArray::TagBoxArray (const BoxArray& ba,
                          int             ngrow)
    :
    m_border(ngrow)
{
    BoxArray grownBoxArray(ba);
    grownBoxArray.grow(ngrow);
    define(grownBoxArray, 1, 0, Fab_allocate);
}

int
TagBoxArray::borderSize () const
{
    return m_border;
}

void 
TagBoxArray::buffer (int nbuf)
{
    if (!(nbuf == 0))
    {
        BL_ASSERT(nbuf <= m_border);
        for (MFIter fai(*this); fai.isValid(); ++fai)
        {
            get(fai).buffer(nbuf, m_border);
        } 
    }
}

void
TagBoxArray::mapPeriodic (const Geometry& geom)
{
    if (!geom.isAnyPeriodic()) return;

    FabArrayCopyDescriptor<TagBox> facd;

    FabArrayId     faid   = facd.RegisterFabArray(this);
    const int      MyProc = ParallelDescriptor::MyProc();
    const Box&     domain = geom.Domain();
    Array<IntVect> pshifts(27);

    std::list< std::pair<FillBoxId,IntVect> > IDs;

    for (int i = 0, N = boxarray.size(); i < N; i++)
    {
        if (!domain.contains(boxarray[i]))
        {
            geom.periodicShift(domain, boxarray[i], pshifts);

            for (int iiv = 0, M = pshifts.size(); iiv < M; iiv++)
            {
                Box shiftbox(boxarray[i]);

                shiftbox.shift(pshifts[iiv]);

                for (int j = 0; j < N; j++)
                {
                    if (distributionMap[j] == MyProc)
                    {
                        Box isect = boxarray[j] & shiftbox;

                        if (isect.ok())
                        {
                            isect.shift(-pshifts[iiv]);

                            FillBoxId fbid = facd.AddBox(faid, isect, 0, i, 0, 0, n_comp);

                            BL_ASSERT(fbid.box() == isect);
                            //
                            // Here we'll save the fab index.
                            //
                            fbid.FabIndex(j);

                            IDs.push_back(std::pair<FillBoxId,IntVect>(fbid,pshifts[iiv]));
                        }
                    }
                }
            }
        }
    }

    facd.CollectData();

    TagBox src;

    for (std::list< std::pair<FillBoxId,IntVect> >::const_iterator it = IDs.begin(), End = IDs.end();
         it != End;
         ++it)
    {
        BL_ASSERT(distributionMap[it->first.FabIndex()] == MyProc);

        src.resize(it->first.box(), n_comp);

        facd.FillFab(faid, it->first, src);

        src.shift(it->second);

        get(it->first.FabIndex()).merge(src);
    }
}

long
TagBoxArray::numTags () const 
{
   long ntag = 0;

   for (MFIter fai(*this); fai.isValid(); ++fai)
   {
      ntag += get(fai).numTags();
   }

   ParallelDescriptor::ReduceLongSum(ntag);

   return ntag;
}

IntVect*
TagBoxArray::collate (long& numtags) const
{
    //
    // The total number of tags system wide that must be collated.
    //
    numtags = numTags();

    int count = 0;
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        count += get(fai).numTags();
    }
    //
    // Local space for holding just those tags we want to gather to the root cpu.
    //
    IntVect*  TheLocalCollateSpace  = new IntVect[count];
    //
    // This holds all tags after they've been gather'd and unique'ified.
    //
    // Each CPU needs an identical copy since they all must go through grid_places() which isn't parallelized.
    //
    // The caller of collate() is responsible for delete[]ing this space.
    //
    IntVect*  TheGlobalCollateSpace = new IntVect[numtags];

    count = 0;
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        get(fai).collate(TheLocalCollateSpace,count);
        count += get(fai).numTags();
    }

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

#if BL_USE_MPI
    Array<int> nmtags(ParallelDescriptor::NProcs(),0);
    Array<int> offset(ParallelDescriptor::NProcs(),0);
    //
    // Tell root CPU how many tags each CPU will be sending.
    //
    MPI_Gather(&count,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nmtags.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 0, N = nmtags.size(); i < N; i++)
            //
            // Convert from count of tags to count of integers to expect.
            //
            nmtags[i] *= BL_SPACEDIM;

        for (int i = 1, N = offset.size(); i < N; i++)
            offset[i] = offset[i-1] + nmtags[i-1];
    }
    //
    // Gather all the tags to IOProc into TheGlobalCollateSpace.
    //
    BL_ASSERT(sizeof(IntVect) == BL_SPACEDIM * sizeof(int));

    MPI_Gatherv(reinterpret_cast<int*>(TheLocalCollateSpace),
                count*BL_SPACEDIM,
                ParallelDescriptor::Mpi_typemap<int>::type(),
                reinterpret_cast<int*>(TheGlobalCollateSpace),
                nmtags.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<int>::type(),
                IOProc,
                ParallelDescriptor::Communicator());
#else
    //
    // Copy TheLocalCollateSpace to TheGlobalCollateSpace.
    //
    for (int i = 0; i < count; i++)
        TheGlobalCollateSpace[i] = TheLocalCollateSpace[i];
#endif

    delete [] TheLocalCollateSpace;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Remove duplicate IntVects.
        //
        std::sort(TheGlobalCollateSpace, TheGlobalCollateSpace+numtags, IntVect::Compare());
        IntVect* end = std::unique(TheGlobalCollateSpace,TheGlobalCollateSpace+numtags);
        ptrdiff_t duplicates = (TheGlobalCollateSpace+numtags) - end;
        BL_ASSERT(duplicates >= 0);
        numtags -= duplicates;
    }
    //
    // Now broadcast them back to the other processors.
    //
    ParallelDescriptor::Bcast(&numtags, 1, IOProc);
    ParallelDescriptor::Bcast(reinterpret_cast<int*>(TheGlobalCollateSpace), numtags*BL_SPACEDIM, IOProc);

    return TheGlobalCollateSpace;
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
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        std::vector< std::pair<int,Box> > isects = ba.intersections(fai.validbox());

        for (int i = 0, N = isects.size(); i < N; i++)
        {
            get(fai).setVal(val,isects[i].second,0);
        }
    } 
}

void
TagBoxArray::coarsen (const IntVect & ratio)
{
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        TagBox* tfine = m_fabs[fai.index()];
        m_fabs[fai.index()] = tfine->coarsen(ratio);
        delete tfine;
    }
    boxarray.coarsen(ratio);
    m_border = 0;
}
