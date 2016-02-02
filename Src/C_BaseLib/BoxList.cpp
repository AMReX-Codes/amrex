
#include <winstd.H>

#include <algorithm>
#include <iostream>

#include <BoxArray.H>
#include <BoxList.H>
#include <BLProfiler.H>

void
BoxList::clear ()
{
    //
    // Really clear out the boxes.
    //
    std::list<Box>().swap(lbox);
}

void
BoxList::join (const BoxList& blist)
{
    BL_ASSERT(ixType() == blist.ixType());
    std::list<Box> lb = blist.lbox;
    lbox.splice(lbox.end(), lb);
}

void
BoxList::catenate (BoxList& blist)
{
    BL_ASSERT(ixType() == blist.ixType());
    lbox.splice(lbox.end(), blist.lbox);
    BL_ASSERT(blist.isEmpty());
}

BoxList&
BoxList::remove (const Box& bx)
{
    BL_ASSERT(ixType() == bx.ixType());
    lbox.remove(bx);
    return *this;
}

BoxList&
BoxList::remove (iterator bli)
{
    BL_ASSERT(ixType() == bli->ixType());
    lbox.erase(bli);
    return *this;
}

BoxList
BoxLib::intersect (const BoxList& bl,
		   const Box&     b)
{
    BL_ASSERT(bl.ixType() == b.ixType());
    BoxList newbl(bl);
    newbl.intersect(b);
    return newbl;
}

BoxList
BoxLib::intersect (const BoxList& bl,
                   const BoxList& br)
{
    BL_ASSERT(bl.ixType() == br.ixType());
    BoxList newbl(bl);
    newbl.intersect(br);
    return newbl;
}

BoxList
BoxLib::refine (const BoxList& bl,
		int            ratio)
{
    BoxList nbl(bl);
    nbl.refine(ratio);
    return nbl;
}

BoxList
BoxLib::coarsen (const BoxList& bl,
                 int            ratio)
{
    BoxList nbl(bl);
    nbl.coarsen(ratio);
    return nbl;
}

BoxList
BoxLib::accrete (const BoxList& bl,
                 int            sz)
{
    BoxList nbl(bl);
    nbl.accrete(sz);
    return nbl;
}

BoxList
BoxLib::removeOverlap (const BoxList& bl)
{
    BoxArray ba(bl);
    ba.removeOverlap();
    return ba.boxList();
}

bool
BoxList::operator!= (const BoxList& rhs) const
{
    return !operator==(rhs);
}

BoxList::BoxList ()
    :
    lbox(),
    btype(IndexType::TheCellType())
{}

BoxList::BoxList (const Box& bx)
    : btype(bx.ixType())
{
    push_back(bx);
}

BoxList::BoxList (IndexType _btype)
    :
    lbox(),
    btype(_btype)
{}

BoxList::BoxList (const BoxArray &ba)
    :
    lbox(),
    btype()
{
    if (ba.size() > 0)
        btype = ba.ixType();
    for (int i = 0, N = ba.size(); i < N; ++i)
        push_back(ba[i]);
}

BoxList::BoxList(const Box& bx, const IntVect& tilesize)
    : btype(bx.ixType())
{
    int ntiles = 1;
    IntVect nt;
    for (int d=0; d<BL_SPACEDIM; d++) {
	nt[d] = (bx.length(d)+tilesize[d]-1)/tilesize[d];
	ntiles *= nt[d];
    }

    IntVect small, big, ijk;  // note that the initial values are all zero.
    ijk[0] = -1;
    for (int t=0; t<ntiles; ++t) {
	for (int d=0; d<BL_SPACEDIM; d++) {
	    if (ijk[d]<nt[d]-1) {
		ijk[d]++;
		break;
	    } else {
		ijk[d] = 0;
	    }
	}

	for (int d=0; d<BL_SPACEDIM; d++) {
	    small[d] = ijk[d]*tilesize[d];
	    big[d] = std::min(small[d]+tilesize[d]-1, bx.length(d)-1);
	}

	Box tbx(small, big, btype);
	tbx.shift(bx.smallEnd());
	push_back(tbx);
    }
}

bool
BoxList::ok () const
{
    const_iterator bli = begin(), End = end();
    if ( bli != End )
    {
        for (Box b(*bli); bli != End; ++bli)
            if (!(bli->ok() && bli->sameType(b)))
                return false;
    }
    return true;
}

bool
BoxList::isDisjoint () const
{
    for (const_iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        const_iterator bli2 = bli;
        //
        // Skip the first element.
        //
        ++bli2;
        for (; bli2 != End; ++bli2)
            if (bli->intersects(*bli2))
                return false;
    }
    return true;
}

bool
BoxList::contains (const IntVect& v) const
{
    for (const_iterator bli = begin(), End = end(); bli != End; ++bli)
        if (bli->contains(v))
            return true;

    return false;
}

bool
BoxList::contains (const Box& b) const
{
    if (isEmpty()) return false;

    BL_ASSERT(ixType() == b.ixType());

    BoxList bnew = BoxLib::complementIn(b,*this);

    return bnew.isEmpty();
}

bool
BoxList::contains (const BoxList&  bl) const
{
    if (isEmpty() || bl.isEmpty()) return false;

    BL_ASSERT(ixType() == bl.ixType());

    if (!minimalBox().contains(bl.minimalBox())) return false;

    BoxArray ba(*this);

    for (const_iterator bli = bl.begin(), End = bl.end(); bli != End; ++bli)
        if (!ba.contains(*bli))
            return false;

    return true;
}

bool
BoxList::contains (const BoxArray&  ba) const
{
    BoxArray tba(*this);
    return ba.contains(tba);
}

BoxList&
BoxList::intersect (const Box& b)
{
    BL_ASSERT(ixType() == b.ixType());

    for (iterator bli = begin(), End = end(); bli != End; )
    {
        const Box& bx = *bli & b;

        if (bx.ok())
        {
            *bli = bx;
            ++bli;
        }
        else
        {
            lbox.erase(bli++);
        }
    }
    return *this;
}

BoxList&
BoxList::intersect (const BoxList& b)
{
    BL_ASSERT(ixType() == b.ixType());

    BoxList bl(b.ixType());

    for (iterator lhs = begin(); lhs != end(); ++lhs)
    {
        for (const_iterator rhs = b.begin(), End = b.end(); rhs != End; ++rhs)
        {
            const Box& bx = *lhs & *rhs;
            if (bx.ok())
                bl.push_back(bx);
        }
    }

    *this = bl;

    return *this;
}

BoxList
BoxLib::complementIn (const Box&     b,
                      const BoxList& bl)
{
    BL_ASSERT(bl.ixType() == b.ixType());
    BoxList newb(b.ixType());
    newb.complementIn(b,bl);
    return newb;
}

BoxList&
BoxList::complementIn (const Box&     b,
                       const BoxList& bl)
{
    BL_ASSERT(bl.ixType() == b.ixType());

    if (bl.size() == 1)
    {
        *this = BoxLib::boxDiff(b,bl.front());
    }
    else
    {
        clear();

        Box     mbox = bl.minimalBox();
        BoxList diff = BoxLib::boxDiff(b,mbox);

        catenate(diff);

        BoxArray ba(bl);

        BoxList mesh(b.ixType());
        if (mbox.ok())
            mesh.push_back(mbox);
        mesh.maxSize(BL_SPACEDIM == 3 ? 64 : 128);

        std::vector< std::pair<int,Box> > isects;

        for (BoxList::const_iterator bli = mesh.begin(), End = mesh.end(); bli != End; ++bli)
        {
            const Box& bx = *bli & b;

            if (!bx.ok()) continue;

            ba.intersections(bx,isects);

            if (isects.empty())
            {
                push_back(bx);
            }
            else
            {
                BoxList tm(b.ixType()), tmpbl(b.ixType());
                for (int i = 0, N = isects.size(); i < N; i++)
                    tmpbl.push_back(isects[i].second);
                tm.complementIn_base(bx,tmpbl);
                catenate(tm);
            }
        }
    }

    return *this;
}

BoxList&
BoxList::complementIn_base (const Box&     b,
                            const BoxList& bl)
{
    BL_ASSERT(bl.ixType() == b.ixType());

    clear();

    push_back(b);

    BoxList diff;

    for (const_iterator bli = bl.begin(), End = bl.end(); bli != End && isNotEmpty(); ++bli)
    {
        for (iterator newbli = lbox.begin(); newbli != lbox.end(); )
        {
            if (newbli->intersects(*bli))
            {
                diff = BoxLib::boxDiff(*newbli, *bli);
                lbox.splice(lbox.begin(), diff.lbox);
                lbox.erase(newbli++);
            }
            else
            {
                ++newbli;
            }
        }
    }

    return *this;
}

BoxList&
BoxList::refine (int ratio)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->refine(ratio);
    }
    return *this;
}

BoxList&
BoxList::refine (const IntVect& ratio)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->refine(ratio);
    }
    return *this;
}

BoxList&
BoxList::coarsen (int ratio)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->coarsen(ratio);
    }
    return *this;
}

BoxList&
BoxList::coarsen (const IntVect& ratio)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->coarsen(ratio);
    }
    return *this;
}

BoxList&
BoxList::accrete (int sz)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->grow(sz);
    }
    return *this;
}

BoxList&
BoxList::accrete (IntVect sz)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->grow(sz);
    }
    return *this;
}

BoxList&
BoxList::shift (int dir,
                int nzones)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->shift(dir, nzones);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (int dir,
                    int num_halfs)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->shiftHalf(dir, num_halfs);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (const IntVect& iv)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->shiftHalf(iv);
    }
    return *this;
}

//
// Returns a list of boxes defining the compliment of b2 in b1in.
//

BoxList
BoxLib::boxDiff (const Box& b1in,
		 const Box& b2)
{
   BL_ASSERT(b1in.sameType(b2));
  
   Box b1(b1in);
   BoxList b_list(b1.ixType());

   if ( !b2.contains(b1) )
   {
       if ( !b1.intersects(b2) )
       {
           b_list.push_back(b1);
       }
       else
       {
           const int* b2lo = b2.loVect();
           const int* b2hi = b2.hiVect();

           for (int i = BL_SPACEDIM-1; i >= 0; i--)
           {
               const int* b1lo = b1.loVect();
               const int* b1hi = b1.hiVect();

               if ((b1lo[i] < b2lo[i]) && (b2lo[i] <= b1hi[i]))
               {
                   Box bn(b1);
                   bn.setSmall(i,b1lo[i]);
                   bn.setBig(i,b2lo[i]-1);
                   b_list.push_back(bn);
                   b1.setSmall(i,b2lo[i]);
               }
               if ((b1lo[i] <= b2hi[i]) && (b2hi[i] < b1hi[i]))
               {
                   Box bn(b1);
                   bn.setSmall(i,b2hi[i]+1);
                   bn.setBig(i,b1hi[i]);
                   b_list.push_back(bn);
                   b1.setBig(i,b2hi[i]);
               }
           }
       }
   }
   return b_list;
}

namespace
{
    struct BoxCmp
    {
        bool operator () (const Box& lhs,
                          const Box& rhs) const
            {
                return lhs.smallEnd().lexLT(rhs.smallEnd());
            }
    };
}

int
BoxList::simplify (bool best)
{
    lbox.sort(BoxCmp());

    return simplify_doit(best);
}

int
BoxList::simplify_doit (bool best)
{
    //
    // Try to merge adjacent boxes.
    //
    int count = 0, lo[BL_SPACEDIM], hi[BL_SPACEDIM];

    for (iterator bla = begin(), End = end(); bla != End; )
    {
        const int* alo   = bla->loVect();
        const int* ahi   = bla->hiVect();
        bool       match = false;
        iterator   blb   = bla;
        ++blb;
        //
        // If we're not looking for the "best" we can do in one pass, we
        // limit how far afield we look for abutting boxes.  This greatly
        // speeds up this routine for large numbers of boxes.  It does not
        // do quite as good a job though as full brute force.
        //
        const int MaxCnt = (best ? size() : 100);

        for (int cnt = 0; blb != End && cnt < MaxCnt; cnt++)
        {
            const int* blo = blb->loVect();
            const int* bhi = blb->hiVect();
            //
            // Determine if a and b can be coalesced.
            // They must have equal extents in all index directions
            // except possibly one, and must abutt in that direction.
            //
            bool canjoin = true;
            int  joincnt = 0;
            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                if (alo[i]==blo[i] && ahi[i]==bhi[i])
                {
                    lo[i] = alo[i];
                    hi[i] = ahi[i];
                }
                else if (alo[i]<=blo[i] && blo[i]<=ahi[i]+1)
                {
                    lo[i] = alo[i];
                    hi[i] = std::max(ahi[i],bhi[i]);
                    joincnt++;
                }
                else if (blo[i]<=alo[i] && alo[i]<=bhi[i]+1)
                {
                    lo[i] = blo[i];
                    hi[i] = std::max(ahi[i],bhi[i]);
                    joincnt++;
                }
                else
                {
                    canjoin = false;
                    break;
                }
            }
            if (canjoin && (joincnt <= 1))
            {
                //
                // Modify b and remove a from the list.
                //
                blb->setSmall(IntVect(lo));
                blb->setBig(IntVect(hi));
                lbox.erase(bla++);
                count++;
                match = true;
                break;
            }
            else
            {
                //
                // No match found, try next element.
                //
                ++blb;
            }
        }
        //
        // If a match was found, a was already advanced in the list.
        //
        if (!match)
            ++bla;
    }
    return count;
}

int
BoxList::minimize ()
{
    int cnt = 0;
    for (int n; (n=simplify(true)) > 0; )
        cnt += n;
    return cnt;
}

Box
BoxList::minimalBox () const
{
    Box minbox(IntVect::TheUnitVector(), IntVect::TheZeroVector(), ixType());
    if ( !isEmpty() )
    {
        const_iterator bli = begin(), End = end();
        minbox = *bli;
        while ( bli != End )
	{
            minbox.minBox(*bli++);
	}
    }
    return minbox;
}

BoxList&
BoxList::maxSize (const IntVect& chunk)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        IntVect boxlen = bli->size();
        const int* len = boxlen.getVect();

        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            if (len[i] > chunk[i])
            {
                //
                // Reduce by powers of 2.
                //
                int ratio = 1;
                int bs    = chunk[i];
                int nlen  = len[i];
                while ((bs%2 == 0) && (nlen%2 == 0))
                {
                    ratio *= 2;
                    bs    /= 2;
                    nlen  /= 2;
                }
                //
                // Determine number and size of (coarsened) cuts.
                //
                const int numblk = nlen/bs + (nlen%bs ? 1 : 0);
                const int size   = nlen/numblk;
                const int extra  = nlen%numblk;
                //
                // Number of cuts = number of blocks - 1.
                //
                for (int k = 0; k < numblk-1; k++)
                {
                    //
                    // Compute size of this chunk, expand by power of 2.
                    //
                    const int ksize = (k < extra ? size+1 : size) * ratio;
                    //
                    // Chop from high end.
                    //
                    const int pos = bli->bigEnd(i) - ksize + 1;

                    push_back(bli->chop(i,pos));
                }
            }
        }
        //
        // b has been chopped down to size and pieces split off
        // have been added to the end of the list so that they
        // can be checked for splitting (in other directions) later.
        //
    }
    return *this;
}

BoxList&
BoxList::maxSize (int chunk)
{
    return maxSize(IntVect(D_DECL(chunk,chunk,chunk)));
}

BoxList&
BoxList::surroundingNodes ()
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->surroundingNodes();
    }
    return *this;
}

BoxList&
BoxList::surroundingNodes (int dir)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->surroundingNodes(dir);
    }
    return *this;
}

BoxList&
BoxList::enclosedCells ()
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->enclosedCells();
    }
    return *this;
}

BoxList&
BoxList::enclosedCells (int dir)
{
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->enclosedCells(dir);
    }
    return *this;
}

BoxList&
BoxList::convert (IndexType typ)
{
    btype = typ;
    for (iterator bli = begin(), End = end(); bli != End; ++bli)
    {
        bli->convert(typ);
    }
    return *this;
}

std::ostream&
operator<< (std::ostream&  os,
            const BoxList& blist)
{
    BoxList::const_iterator bli = blist.begin(), End = blist.end();
    os << "(BoxList " << blist.size() << ' ' << blist.ixType() << '\n';
    for (int count = 1; bli != End; ++bli, ++count)
    {
        os << count << " : " << *bli << '\n';
    }
    os << ')' << '\n';

    if (os.fail())
        BoxLib::Error("operator<<(ostream&,BoxList&) failed");

    return os;
}

bool
BoxList::operator== (const BoxList& rhs) const
{
    if ( !(size() == rhs.size()) ) return false;

    BoxList::const_iterator liter = begin(), riter = rhs.begin(), End = end();
    for (; liter != End; ++liter, ++riter)
        if ( !( *liter == *riter) )
            return false;
    return true;
}
