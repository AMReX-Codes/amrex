//BL_COPYRIGHT_NOTICE

//
// $Id: TagBox.cpp,v 1.14 1998-02-05 17:02:28 lijewski Exp $
//

#include <TagBox.H>
#include <Misc.H>
#include <Geometry.H>
#include <ParallelDescriptor.H>

extern void inspectTAGArray (const TagBoxArray& tba);
extern void inspectTAG (const TagBox& tb, int n);
extern void inspectFAB (FArrayBox& unfab, int n);

//
// Number of IntVects that can fit into m_CollateSpace
//
long TagBoxArray::m_CollateCount = TagBoxArray::ChunkSize;

//
// Static space used by collate().
//
IntVect* TagBoxArray::m_CollateSpace = new IntVect[TagBoxArray::ChunkSize];

void
TagBoxArray::BumpCollateSpace (long numtags)
{
    assert(TagBoxArray::m_CollateCount < numtags);

    do
    {
        TagBoxArray::m_CollateCount += TagBoxArray::ChunkSize;
    }
    while (TagBoxArray::m_CollateCount < numtags);

    delete [] TagBoxArray::m_CollateSpace;

    TagBoxArray::m_CollateSpace = new IntVect[TagBoxArray::m_CollateCount];
}

struct TagBoxMergeDesc
{
  bool      destLocal;
  int       mergeIndexSrc;
  int       mergeIndexDest;
  int       nOverlap;
  Box       overlapBox;
  FillBoxId fillBoxId;
};

TagBox::TagBox () {}

TagBox::TagBox (const Box& bx,
                int        n)
    :
    BaseFab<TagBox::TagType>(bx,n)
{
    setVal(TagBox::CLEAR);
}

TagBox::~TagBox () {}

TagBox*
TagBox::coarsen (const IntVect& ratio)
{
    Box cbx(domain);
    cbx.coarsen(ratio);
    TagBox* crse = new TagBox(cbx);
    const Box& cbox = crse->box();
    Box b1(::refine(cbox,ratio));

    const int* flo  = domain.loVect();
    const int* fhi  = domain.hiVect();
    const int* flen = domain.length().getVect();

    const int* clo  = cbox.loVect();
    const int* chi  = cbox.hiVect();
    const int* clen = cbox.length().getVect();

    const int* lo  = b1.loVect();
    const int* hi  = b1.hiVect();
    const int* len = b1.length().getVect();

    int longlen, dir;
    longlen = b1.longside(dir);

    TagType* fdat = dataPtr();
    TagType* cdat = crse->dataPtr();

    TagType* t = new TagType[longlen];
    for (int i = 0; i < longlen; i++)
        t[i] = TagBox::CLEAR;

    int klo = 0, khi = 0, jlo = 0, jhi = 0, ilo, ihi;
    D_TERM(ilo=flo[0]; ihi=fhi[0]; ,
           jlo=flo[1]; jhi=fhi[1]; ,
           klo=flo[2]; khi=fhi[2];)

#define IXPROJ(i,r) (((i)+(r)*Abs(i))/(r) - Abs(i))
#define IOFF(j,k,lo,len) D_TERM(0, +(j-lo[1])*len[0], +(k-lo[2])*len[0]*len[1])
#define JOFF(i,k,lo,len) D_TERM(i-lo[0], +0, +(k-lo[2])*len[0]*len[1])
#define KOFF(i,j,lo,len) D_TERM(i-lo[0], +(j-lo[1])*len[0], +0)
   //
   // hack
   //
   dir = 0;
   
   int ratiox = 1, ratioy = 1, ratioz = 1;
   D_TERM(ratiox = ratio[0];,
          ratioy = ratio[1];,
          ratioz = ratio[2];)

   int dummy_ratio = 1;

   if (dir == 0)
   {
      for (int k = klo; k <= khi; k++)
      {
          int kc = IXPROJ(k,ratioz);
          for (int j = jlo; j <= jhi; j++)
          {
              int jc = IXPROJ(j,ratioy);
              TagType* c = cdat + IOFF(jc,kc,clo,clen);
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
                      int i = ic*ratiox + off;
                      c[ic] = Max(c[ic],t[i]);
                  }
              }
          }
      }
   }
   else if (dir == 1)
   {
      for (int k = klo; k <= khi; k++)
      {
          int kc = IXPROJ(k,dummy_ratio);
          for (int i = ilo; i <= ihi; i++)
          {
              int ic = IXPROJ(i,dummy_ratio);
              TagType* c = cdat + JOFF(ic,kc,clo,clen);
              const TagType* f = fdat + JOFF(i,k,flo,flen);
              //
              // Copy fine grid row of values into tmp array.
              //
              int strd = flen[0];
              for (int j = jlo; j <= jhi; j++)
                  t[j-lo[1]] = f[(j-jlo)*strd];

              for (int off = 0; off < dummy_ratio; off++)
              {
                  int jc = 0;
                  strd = clen[0];
                  for (int jcnt = 0; jcnt < clen[1]; jcnt++)
                  {
                      int j = jcnt*dummy_ratio + off;
                      c[jc] = Max(c[jc],t[j]);
                      jc += strd;
                  }
              }
          }
      }
   }
   else
   {
       for (int j = jlo; j <= jhi; j++)
       {
           int jc = IXPROJ(j,dummy_ratio);
           for (int i = ilo; i <= ihi; i++)
           {
               int ic = IXPROJ(i,dummy_ratio);
               TagType* c = cdat + KOFF(ic,jc,clo,clen);
               const TagType* f = fdat + KOFF(i,j,flo,flen);
               //
               // Copy fine grid row of values into tmp array.
               //
               int strd = flen[0]*flen[1];
               for (int k = klo; k <= khi; k++)
                   t[k-lo[2]] = f[(k-klo)*strd];

               for (int off = 0; off < dummy_ratio; off++)
               {
                   int kc = 0;
                   strd = clen[0]*clen[1];
                   for (int kcnt = 0; kcnt < clen[2]; kcnt++)
                   {
                       int k = kcnt*dummy_ratio + off;
                       c[kc] = Max(c[kc],t[k]);
                       kc += strd;
                   }
               }
           }
      }
   }

   delete [] t;

   return crse;

#undef ABS
#undef IXPROJ
#undef IOFF
#undef JOFF
#undef KOFF
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

    const int* len = domain.length().getVect();
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
    Box bx(domain);
    bx &= src.domain;
    if (bx.ok())
    {
        const int* dlo = domain.loVect();
        const int* dlen = domain.length().getVect();
        const int* slo = src.domain.loVect();
        const int* slen = src.domain.length().getVect();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        const TagType* ds0 = src.dataPtr();
        TagType* dd0 = dataPtr();
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
   assert(t_long < INT_MAX);
   int len = int(t_long);
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

void
TagBox::collate (IntVect* ar,
                 int      start) const
{
    assert(!(ar == 0));
    assert(start >= 0);
    //
    // Starting at given offset of array ar, enter location (IntVect) of
    // each tagged cell in tagbox.
    //
    const int* len = domain.length().getVect();
    const int* lo = domain.loVect();
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
                }
            }
        }
    }
}

Array<int>
TagBox::tags () const
{
    Array<int> ar(domain.numPts(), TagBox::CLEAR);

    const TagType* cptr = dataPtr();
    int*           iptr = ar.dataPtr();

    for (int i = 0; i < ar.length(); i++, cptr++, iptr++)
    {
        if (*cptr)
            *iptr = *cptr;
    }

    return ar;
}

void
TagBox::tags (const Array<int>& ar)
{
    assert(ar.length() == domain.numPts());

    TagType*   cptr = dataPtr();
    const int* iptr = ar.dataPtr();

    for (int i = 0; i < ar.length(); i++, cptr++, iptr++)
    {
        if (*iptr)
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

TagBoxArray::~TagBoxArray () {}

void 
TagBoxArray::buffer (int nbuf)
{
    if (!(nbuf == 0))
    {
        assert(nbuf <= m_border);
        for (FabArrayIterator<TagType,TagBox> fai(*this); fai.isValid(); ++fai)
        {
            fai().buffer(nbuf, m_border);
        } 
    }
}

void
TagBoxArray::mergeUnique ()
{
    FabArrayCopyDescriptor<TagType,TagBox> facd(true);
    FabArrayId faid = facd.RegisterFabArray(this);
    int nOverlap = 0;
    int myproc = ParallelDescriptor::MyProc();
    List<TagBoxMergeDesc> tbmdList;
    List<FabComTag> tbmdClearList;
    //
    // Don't use FabArrayIterator here.
    //
    BoxList unfilledBoxes;  // Returned by AddBox, not used.
    for (int idest = 0; idest < fabparray.length(); ++idest)
    {
        bool destLocal = (distributionMap[idest] == myproc);
        for (int isrc = idest + 1; isrc < fabparray.length(); ++isrc)
        {
            Box ovlp(boxarray[idest]);
            ovlp &= boxarray[isrc];
            if (ovlp.ok())
            {
                TagBoxMergeDesc tbmd;
                tbmd.destLocal      = destLocal;
                tbmd.mergeIndexSrc  = isrc;
                tbmd.mergeIndexDest = idest;
                tbmd.nOverlap       = nOverlap;
                tbmd.overlapBox     = ovlp;
                if (destLocal)
                {
                    tbmd.fillBoxId = facd.AddBox(faid,ovlp,unfilledBoxes,isrc,0,0,1);
                }
                tbmdList.append(tbmd);
                ++nOverlap;
            }
        }
    }
    facd.CollectData();

    int listIndex = 0;
    for (ListIterator<TagBoxMergeDesc> tbmdli(tbmdList); tbmdli; ++tbmdli)
    {
        const TagBoxMergeDesc& tbmd = tbmdli();
        if (tbmd.destLocal)
        {
            TagBox& dest = fabparray[tbmd.mergeIndexDest];
            TagBox src(tbmd.overlapBox, 1);
            facd.FillFab(faid, tbmd.fillBoxId, src);
            for (ListIterator<TagBoxMergeDesc> tbmdliprev(tbmdList);
                 tbmdliprev && tbmdliprev().nOverlap <= listIndex;
                 ++tbmdliprev)
            {
                Box ovlpBox(src.box());
                ovlpBox &= tbmdliprev().overlapBox;
                if (ovlpBox.ok() && tbmdliprev().mergeIndexSrc == tbmd.mergeIndexSrc)
                {
                    if (tbmdliprev().nOverlap < listIndex)
                    {
                        src.setVal(TagBox::CLEAR, ovlpBox, 0);
                    }
                    FabComTag tbmdClear;
                    tbmdClear.fabIndex = tbmd.mergeIndexSrc;
                    tbmdClear.ovlpBox = tbmd.overlapBox;
                    tbmdClearList.append(tbmdClear);
                }
            }
            dest.merge(src);
        }
        ++listIndex;
    }
    //
    // Now send the clear list elements to the processor they belong to.
    //
    ParallelDescriptor::SetMessageHeaderSize(sizeof(FabComTag));

    ListIterator<FabComTag> tbmdsendli(tbmdClearList);

    for ( ; tbmdsendli; ++tbmdsendli)
    {
        ParallelDescriptor::SendData(distributionMap[tbmdsendli().fabIndex],
                                     &tbmdsendli(), 0, 0);
    }

    ParallelDescriptor::Synchronize();  // To guarantee messages are sent.
    //
    // Now clear the overlaps in the TagBoxArray.
    //
    int dataWaitingSize;
    FabComTag tbmdClear;
    while (ParallelDescriptor::GetMessageHeader(dataWaitingSize, &tbmdClear))
    {
        //
        // Data was sent to this processor.
        //
        if (!(distributionMap[tbmdClear.fabIndex] == myproc))
            BoxLib::Error("tbmdClear.fabIndex is not local");

        TagBox& src = fabparray[tbmdClear.fabIndex];
        src.setVal(TagBox::CLEAR, tbmdClear.ovlpBox, 0);

        ParallelDescriptor::ReceiveData(0, 0);  // To advance message header.
    }
    ParallelDescriptor::Synchronize();
}

void
TagBoxArray::mapPeriodic (const Geometry& geom)
{
    FabArrayCopyDescriptor<TagType,TagBox> facd(true);
    FabArrayId faid = facd.RegisterFabArray(this);
    int myproc = ParallelDescriptor::MyProc();
    List<FillBoxId> fillBoxIdList;
    FillBoxId tempFillBoxId;

    Box domain(geom.Domain());
    TagBox tagtmp;
    int srcComp  = 0;
    int destComp = 0;
    int nComp    = n_comp;
    //
    // This logic needs to be turned inside out to use a FabArrayIterator
    //
    for (int i = 0; i < fabparray.length(); i++)
    {
        if (!domain.contains( boxarray[i]))
        {
            //
            // src is candidate for periodic mapping.
            //
            Array<IntVect> pshifts(27);
            geom.periodicShift( domain, boxarray[i], pshifts );
            for (int iiv = 0; iiv < pshifts.length(); iiv++)
            {
                IntVect iv = pshifts[iiv];
                Box shiftbox( boxarray[i] );
                D_TERM(shiftbox.shift(0,iv[0]);,
                       shiftbox.shift(1,iv[1]);,
                       shiftbox.shift(2,iv[2]);)
                //
                // Possible periodic remapping, try each tagbox.
                //
                for (int j = 0; j < fabparray.length(); j++)
                {
                    if (distributionMap[j] == myproc)
                    {
                        TagBox& dest = fabparray[j];
                        Box intbox = dest.box() & shiftbox;
                        if (intbox.ok())
                        {
                            BoxList unfilledBoxes(intbox.ixType());
                            //
                            // Ok, got a hit, but be careful if is same TagBox.
                            //
                            if (i != j)
                            {
                                tempFillBoxId = facd.AddBox(faid, intbox,
                                                            unfilledBoxes,
                                                            srcComp, destComp,
                                                            nComp);
                                fillBoxIdList.append(tempFillBoxId);
                            }
                            else
                            {
                                //
                                // Is same tagbox, must be careful.
                                //
                                Box shintbox(intbox);
                                IntVect tmpiv( -iv );
                                D_TERM(shintbox.shift(0,tmpiv[0]);,
                                       shintbox.shift(1,tmpiv[1]);,
                                       shintbox.shift(2,tmpiv[2]);)
                                tempFillBoxId = facd.AddBox(faid, shintbox,
                                                            unfilledBoxes,
                                                            srcComp, destComp,
                                                            nComp);
                                fillBoxIdList.append(tempFillBoxId);
                            }
                        }
                    }
                }
            }
        }
    }

    Array<FillBoxId> fillBoxId(fillBoxIdList.length());
    int ifbi = 0;
    for (ListIterator<FillBoxId> li(fillBoxIdList); li; ++li)
    {
        fillBoxId[ifbi] = li();
        ++ifbi;
    }
    fillBoxIdList.clear();

    facd.CollectData();

    int iFillBox = 0;
    //
    // This logic needs to be turned inside out to use a FabArrayIterator.
    //
    for (int i = 0; i < fabparray.length(); i++)
    {
        if (!domain.contains(boxarray[i]))
        {
            //
            // src is candidate for periodic mapping.
            //
            Array<IntVect> pshifts(27);
            geom.periodicShift( domain, boxarray[i], pshifts );
            for (int iiv = 0; iiv < pshifts.length(); iiv++)
            {
                IntVect iv = pshifts[iiv];
                Box shiftbox( boxarray[i] );
                D_TERM(shiftbox.shift(0,iv[0]);,
                       shiftbox.shift(1,iv[1]);,
                       shiftbox.shift(2,iv[2]);)
                //
                // Possible periodic remapping, try each tagbox.
                //
                for (int j = 0; j < fabparray.length(); j++)
                {
                    if (distributionMap[j] == myproc)
                    {
                        //
                        // Local dest fab.
                        //
                        TagBox& dest = fabparray[j];
                        Box intbox = dest.box() & shiftbox;
                        if (intbox.ok())
                        {
                            //
                            // Ok, got a hit, but be careful if is same TagBox.
                            //
                            if (i != j)
                            {
                                FillBoxId fillboxid = fillBoxId[iFillBox];
                                ++iFillBox;
                                TagBox src(fillboxid.box(), n_comp);
                                facd.FillFab(faid, fillboxid, src);
                                src.shift(iv);
                                dest.merge(src);
                                src.shift(-iv);
                            }
                            else
                            {
                                //
                                // Is same tagbox, must be careful.
                                //
                                FillBoxId fillboxid = fillBoxId[iFillBox];
                                ++iFillBox;
                                TagBox src(fillboxid.box(), n_comp);
                                facd.FillFab(faid, fillboxid, src);
                                tagtmp.resize(intbox);
                                Box shintbox(intbox);
                                IntVect tmpiv( -iv );
                                D_TERM(shintbox.shift(0,tmpiv[0]);,
                                       shintbox.shift(1,tmpiv[1]);,
                                       shintbox.shift(2,tmpiv[2]);)
                                assert(shintbox == fillboxid.box());
                                tagtmp.copy(src,shintbox,0,intbox,0,1);
                                dest.merge(tagtmp);
                            }
                        }
                    }
                }
            }
        }
    }
}

long
TagBoxArray::numTags () const 
{
   long ntag = 0;
   for (ConstFabArrayIterator<TagType,TagBox> fai(*this); fai.isValid(); ++fai)
   {
      ntag += fai().numTags();
   } 
   ParallelDescriptor::ReduceLongSum(ntag);
   return ntag;
}

IntVect*
TagBoxArray::collate (long& numtags) const
{
    const int myproc = ParallelDescriptor::MyProc();
    const int nGrids = fabparray.length();
    int* sharedNTags = new int[nGrids]; // Shared numTags per grid.
    for (int isn = 0; isn < nGrids; ++isn)
    {
        sharedNTags[isn] = -1;  // A bad value.
    }
    for (ConstFabArrayIterator<TagType,TagBox> fai(*this); fai.isValid();++fai)
    {
        sharedNTags[fai.index()] = fai().numTags();
    }
    //
    // Communicate number of local tags for each grid.
    //
    const int nProcs = ParallelDescriptor::NProcs();
    ParallelDescriptor::ShareVar(sharedNTags, nGrids * sizeof(int));
    ParallelDescriptor::Synchronize();  // For ShareVar.
    for (int iGrid = 0; iGrid < nGrids; ++iGrid)
    {
        if (sharedNTags[iGrid] != -1)
        {
            for (int iProc = 0; iProc < nProcs; ++iProc)
            {
                if (iProc != myproc)
                    ParallelDescriptor::WriteData(iProc,
                                                  &sharedNTags[iGrid],
                                                  sharedNTags,
                                                  iGrid * sizeof(int),
                                                  sizeof(int));
            }
        }
    }
    ParallelDescriptor::Synchronize();
    ParallelDescriptor::UnshareVar(sharedNTags);

    int* startOffset = new int[nGrids]; // Start locations per grid.
    startOffset[0] = 0;
    for (int iGrid = 1; iGrid < nGrids; ++iGrid)
    {
        startOffset[iGrid] = startOffset[iGrid-1] + sharedNTags[iGrid-1];
    }
    //
    // Communicate all local points so all procs have the same global set.
    //
    // Use our static 1D array for contiguous parallel copies.
    //
    numtags = numTags();

    if (TagBoxArray::m_CollateCount < numtags)
        TagBoxArray::BumpCollateSpace(numtags);

    for (ConstFabArrayIterator<TagType,TagBox> fai(*this); fai.isValid();++fai)
    {
        fai().collate(TagBoxArray::m_CollateSpace, startOffset[fai.index()]);
    }

    const size_t IVSize = sizeof(IntVect);
    //
    // Now copy the local IntVects to all other processors.
    //
    ParallelDescriptor::ShareVar(TagBoxArray::m_CollateSpace, numtags*IVSize);
    ParallelDescriptor::Synchronize();  // For ShareVar.

    for (ConstFabArrayIterator<TagType,TagBox> fai(*this); fai.isValid();++fai)
    {
        int idx = fai.index();

        IntVect* ivDestBase = TagBoxArray::m_CollateSpace + startOffset[idx];

        for (int iProc = 0; iProc < nProcs; ++iProc)
        {
            if (iProc != myproc)
            {
                if (sharedNTags[fai.index()] != 0)
                {
                    ParallelDescriptor::WriteData(iProc,
                                                  ivDestBase,
                                                  TagBoxArray::m_CollateSpace,
                                                  startOffset[idx]*IVSize,
                                                  sharedNTags[idx]*IVSize);
                }
            }
        }
    }
    ParallelDescriptor::Synchronize();  // Need this sync after the put.
    ParallelDescriptor::UnshareVar(TagBoxArray::m_CollateSpace);

    delete [] startOffset;
    delete [] sharedNTags;

    return TagBoxArray::m_CollateSpace;
}

void
TagBoxArray::setVal (BoxDomain&     bd,
                     TagBox::TagVal val)
{
    for (FabArrayIterator<TagType,TagBox> fai(*this); fai.isValid(); ++fai)
    {
        for (BoxDomainIterator bdi(bd); bdi; ++bdi)
        {
            Box bx(fai.validbox());
            bx &= bdi();
            if (bx.ok())
            {
                fai().setVal(val,bx,0);
            }
        }
    }
}

void
TagBoxArray::setVal (BoxArray&      ba,
                     TagBox::TagVal val)
{
    for (FabArrayIterator<TagType,TagBox> fai(*this); fai.isValid(); ++fai)
    {
        for (int j = 0; j < ba.length(); j++)
        {
            Box bx(fai.validbox());
            bx &= ba[j];
            if (bx.ok())
            {
                fai().setVal(val,bx,0);
            }
        }
    } 
}

void 
TagBoxArray::coarsen (const IntVect & ratio)
{
    for (FabArrayIterator<TagType,TagBox> fai(*this); fai.isValid(); ++fai)
    {
        TagBox* tfine = fabparray.remove(fai.index());
        TagBox* tcrse = tfine->coarsen(ratio);
        fabparray.set(fai.index(),tcrse);
        delete tfine;
    } 
    m_border = 0;
    n_grow   = 0;
}
