//
// $Id: BoxList.cpp,v 1.17 2001-07-23 21:15:16 car Exp $
//
#include <winstd.H>

#include <algorithm>

#include <BoxArray.H>
#include <BoxList.H>

BoxList::Iterator
BoxList::begin()
{
    return lbox.begin();
}

BoxList::ConstIterator
BoxList::begin() const
{
    return lbox.begin();
}

BoxList::Iterator
BoxList::end()
{
    return lbox.end();
}

BoxList::ConstIterator
BoxList::end() const
{
    return lbox.end();
}

BoxList::~BoxList()
{}

IndexType
BoxList::ixType () const
{
    return btype;
}

void
BoxList::append (const Box& bn)
{
    BL_ASSERT(ixType() == bn.ixType());
    lbox.push_back(bn);
}

void
BoxList::add (const Box& bn)
{
    append(bn);
}

void
BoxList::prepend (const Box& bn)
{
    BL_ASSERT(ixType() == bn.ixType());
    lbox.push_front(bn);
}

void
BoxList::join (const BoxList& blist)
{
    BL_ASSERT(ixType() == blist.ixType());
    std::list<Box> lb = blist.lbox;
    lbox.splice(lbox.end(), lb);
}

bool
BoxList::isEmpty () const
{
    return lbox.empty();
}

void
BoxList::catenate (BoxList& blist)
{
    BL_ASSERT(ixType() == blist.ixType());
    lbox.splice(lbox.end(), blist.lbox);
    BL_ASSERT(blist.isEmpty());
}

void
BoxList::clear ()
{
    lbox.clear();
}

int
BoxList::length () const
{
    return lbox.size();
}

bool
BoxList::isNotEmpty () const
{
    return !lbox.empty();
}

bool
BoxList::contains (const Box& b) const
{
    BoxList bnew = BoxLib::complementIn(b,*this);

    return bnew.isEmpty();
}

BoxList&
BoxList::remove (const Box& bx)
{
    BL_ASSERT(ixType() == bx.ixType());
    lbox.remove(bx);
    return *this;
}

BoxList&
BoxList::remove (Iterator bli)
{
    BL_ASSERT(ixType() == bli->ixType());
    lbox.erase(bli);
    return *this;
}

BoxList
BoxLib::intersect (const BoxList& bl,
		   const Box&     b)
{
    BoxList newbl(bl);
    return newbl.intersect(b);
}

BoxList
BoxLib::intersect (const BoxList& bl,
           const BoxList& br)
{
    BoxList newbl(bl);
    return newbl.intersect(br);
}

BoxList
BoxLib::refine (const BoxList& bl,
		int            ratio)
{
    BoxList nbl(bl);
    return nbl.refine(ratio);
}

BoxList
BoxLib::coarsen (const BoxList& bl,
         int            ratio)
{
    BoxList nbl(bl);
    return nbl.coarsen(ratio);
}

BoxList
accrete (const BoxList& bl,
         int            sz)
{
    BoxList nbl(bl);
    return nbl.accrete(sz);
}

std::list<Box>&
BoxList::listBox()
{
    return lbox;
}

const std::list<Box>&
BoxList::listBox() const
{
    return lbox;
}

bool
BoxList::operator!= (const BoxList& rhs) const
{
    return !operator==(rhs);
}

BoxList::BoxList ()
    : lbox(), btype(IndexType::TheCellType())
{}

BoxList::BoxList (IndexType _btype)
    : lbox(), btype(_btype)
{}

BoxList::BoxList (const BoxList &_blst)
    : lbox(_blst.lbox),
      btype(_blst.btype)
{}

BoxList&
BoxList::operator= (const BoxList& rhs)
{
    lbox  = rhs.lbox;
    btype = rhs.btype;
    return *this;
}

BoxList::BoxList (const BoxArray &ba)
    : lbox(),
      btype()
{
    if (ba.length() > 0)
        btype = ba[0].ixType();
    for (int i = 0; i < ba.length(); ++i)
        append(ba[i]);
}

bool
BoxList::ok () const
{
    bool isok = true;
    ConstIterator bli = begin();
    if ( bli != end() )
    {
        for (Box b(*bli); bli != end() && isok; ++bli)
	{
            isok = bli->ok() && bli->sameType(b);
	}
    }
    return isok;
}

bool
BoxList::isDisjoint () const
{
    bool isdisjoint = true;
    for (ConstIterator bli = begin(); bli != end() && isdisjoint; ++bli)
    {
        ConstIterator bli2 = bli;
        //
        // Skip the first element.
        //
        ++bli2; 
        for (; bli2 != end() && isdisjoint; ++bli2)
	{
            if (bli->intersects(*bli2))
	    {
                isdisjoint = false;
	    }
	}
    }
    return isdisjoint;
}

bool
BoxList::contains (const IntVect& v) const
{
    bool contained = false;
    for (ConstIterator bli = begin(); bli != end() && !contained; ++bli)
    {
        if (bli->contains(v))
	{
            contained = true;
	}
    }
    return contained;
}

bool
BoxList::contains (const BoxList&  bl) const
{
    for (ConstIterator bli = bl.begin(); bli != bl.end(); ++bli)
    {
	if ( !contains(*bli) )
	{
	    return false;
	}
    }
    return true;
}

bool
BoxList::contains (const BoxArray&  ba) const
{
    for (int i = 0; i < ba.length(); i++)
        if (!contains(ba[i]))
            return false;
    return true;
}

BoxList&
BoxList::intersect (const Box& b)
{
    for (Iterator bli= begin(); bli != end(); )
    {
        if (bli->intersects(b))
        {
            *bli &= b;
            ++bli;
        }
        else
        {
            bli = lbox.erase(bli);
        }
    }
    return *this;
}

BoxList&
BoxList::intersect (const BoxList& b)
{
    BoxList bl;

    for (Iterator lhs = begin(); lhs != end(); ++lhs)
    {
        for (ConstIterator rhs = b.begin(); rhs != b.end(); ++rhs)
        {
            if ( lhs->intersects(*rhs) )
            {
                bl.add(*lhs & *rhs);
            }
        }
    }

    *this = bl;

    return *this;
}

BoxList
BoxLib::complementIn (const Box&     b,
                      const BoxList& bl)
{
    BoxList newb(b.ixType());
    newb.append(b);
    for (BoxList::ConstIterator bli = bl.begin(); bli != bl.end() && newb.isNotEmpty(); ++bli)
    {
        for (BoxList::Iterator newbli = newb.begin(); newbli != newb.end(); )
        {
            if ( newbli->intersects(*bli))
            {
                BoxList tm = BoxLib::boxDiff(*newbli, *bli);
                newb.catenate(tm);
                newb.remove(newbli++);
            }
            else
            {
                ++newbli;
            }
        }
    }
    return newb;
}

BoxList&
BoxList::complementIn (const Box&     b,
                       const BoxList& bl)
{
    clear();
    append(b);
    for (ConstIterator bli = bl.begin(); bli != bl.end() && isNotEmpty(); ++bli)
    {
        for (Iterator newbli = lbox.begin(); newbli != lbox.end(); )
        {
            if ( newbli->intersects(*bli) )
            {
                BoxList tm = BoxLib::boxDiff(*newbli, *bli);
                lbox.splice(lbox.end(), tm.lbox);
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
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->refine(ratio);
    }
    return *this;
}

BoxList&
BoxList::refine (const IntVect& ratio)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->refine(ratio);
    }
    return *this;
}

BoxList&
BoxList::coarsen (int ratio)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->coarsen(ratio);
    }
    return *this;
}

BoxList&
BoxList::coarsen (const IntVect& ratio)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->coarsen(ratio);
    }
    return *this;
}

BoxList&
BoxList::accrete (int sz)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->grow(sz);
    }
    return *this;
}

BoxList&
BoxList::shift (int dir,
                int nzones)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->shift(dir, nzones);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (int dir,
                    int num_halfs)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->shiftHalf(dir, num_halfs);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (const IntVect& iv)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
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
   Box b1(b1in);
   BoxList b_list(b1.ixType());

   if ( !b2.contains(b1) )
   {
       if ( !b1.intersects(b2) )
       {
           b_list.append(b1);
       }
       else
       {
           const int* b2lo = b2.loVect();
           const int* b2hi = b2.hiVect();

           for (int i = 0; i < BL_SPACEDIM; i++)
           {
               const int* b1lo = b1.loVect();
               const int* b1hi = b1.hiVect();

               if ((b1lo[i] < b2lo[i]) && (b2lo[i] <= b1hi[i]))
               {
                   Box bn(b1);
                   bn.setSmall(i,b1lo[i]);
                   bn.setBig(i,b2lo[i]-1);
                   b_list.append(bn);
                   b1.setSmall(i,b2lo[i]);
               }
               if ((b1lo[i] <= b2hi[i]) && (b2hi[i] < b1hi[i]))
               {
                   Box bn(b1);
                   bn.setSmall(i,b2hi[i]+1);
                   bn.setBig(i,b1hi[i]);
                   b_list.append(bn);
                   b1.setBig(i,b2hi[i]);
               }
           }
       }
   }
   return b_list;
}

int
BoxList::simplify ()
{
    //
    // Try to merge adjacent boxes.
    //
    int count = 0;
    int lo[BL_SPACEDIM];
    int hi[BL_SPACEDIM];

    for (Iterator bla = begin(); bla != end(); )
    {
        const int* alo   = bla->loVect();
        const int* ahi   = bla->hiVect();
        bool       match = false;
        Iterator blb = bla;
        ++blb;
        while ( blb != end() )
        {
            const int* blo = blb->loVect();
            const int* bhi = blb->hiVect();
            //
            // Determine of a and b can be coalasced.
            // They must have equal extents in all index direciton
            // except possibly one and must abutt in that direction.
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
    for (int n; (n=simplify()) > 0; )
        cnt += n;
    return cnt;
}

Box
BoxList::minimalBox () const
{
    Box minbox;
    if ( !isEmpty() )
    {
        ConstIterator bli = begin();
        minbox = *bli;
        while ( bli != end() )
	{
            minbox.minBox(*bli++);
	}
    }
    return minbox;
}

BoxList&
BoxList::maxSize (const IntVect& chunk)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        const int* len = bli->length().getVect();

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

                    append(bli->chop(i,pos));
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
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->surroundingNodes();
    }
    return *this;
}

BoxList&
BoxList::surroundingNodes (int dir)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->surroundingNodes(dir);
    }
    return *this;
}

BoxList&
BoxList::enclosedCells ()
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->enclosedCells();
    }
    return *this;
}

BoxList&
BoxList::enclosedCells (int dir)
{
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->enclosedCells(dir);
    }
    return *this;
}

BoxList&
BoxList::convert (IndexType typ)
{
    btype = typ;
    for (Iterator bli = begin(); bli != end(); ++bli)
    {
        bli->convert(typ);
    }
    return *this;
}

std::ostream&
operator<< (std::ostream&  os,
            const BoxList& blist)
{
    BoxList::ConstIterator bli = blist.begin();
    os << "(BoxList " << blist.length() << ' ' << blist.ixType() << '\n';
    for (int count = 1; bli != blist.end(); ++bli, ++count)
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
    bool rc = true;
    if (length() != rhs.length())
    {
        rc = false;
    }
    else
    {
        BoxList::ConstIterator liter = begin(), riter = rhs.begin();
        for (; liter != end() && rc; ++liter, ++riter)
	{
            if ( !( *liter == *riter) )
	    {
                rc = false;
	    }
	}
    }
    return rc;
}
