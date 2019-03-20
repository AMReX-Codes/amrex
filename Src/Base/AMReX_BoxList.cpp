
#include <algorithm>
#include <iostream>
#include <cmath>

#include <AMReX_Print.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_BLProfiler.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

namespace {

static void chop_boxes (Box* bxv, const Box& bx, int nboxes)
{
    if (nboxes == 1)
    {
        *bxv = bx;
    }
    else
    {
        int longdir;
        int longlen = bx.longside(longdir);
        int chop_pnt = bx.smallEnd(longdir) + longlen/2;
        Box bxleft = bx;
        Box bxright = bxleft.chop(longdir, chop_pnt);

        int nleft = nboxes / 2;
        chop_boxes(bxv, bxleft, nleft);

        int nright = nboxes - nleft;
        chop_boxes(bxv+nleft, bxright, nright);
    }
}

static void chop_boxes_dir (Box* bxv, const Box& bx, int nboxes, int idir)
{
    if (nboxes == 1)
    {
        *bxv = bx;
    }
    else
    {
        int chop_pnt = bx.smallEnd(idir) + bx.length(idir)/2;
        Box bxleft = bx;
        Box bxright = bxleft.chop(idir, chop_pnt);

        int nleft = nboxes / 2;
        chop_boxes_dir(bxv, bxleft, nleft, idir);

        int nright = nboxes - nleft;
        chop_boxes_dir(bxv+nleft, bxright, nright, idir);
    }
}

}

void
BoxList::clear ()
{
    m_lbox.clear();
}

void
BoxList::join (const BoxList& blist)
{
    BL_ASSERT(ixType() == blist.ixType());
    m_lbox.insert(std::end(m_lbox), std::begin(blist), std::end(blist));
}

void
BoxList::join (const Vector<Box>& barr)
{
    BL_ASSERT(barr.size() == 0 || ixType() == barr[0].ixType());
    m_lbox.insert(std::end(m_lbox), std::begin(barr), std::end(barr));
}

void
BoxList::catenate (BoxList& blist)
{
    BL_ASSERT(ixType() == blist.ixType());
    m_lbox.insert(std::end(m_lbox), std::begin(blist), std::end(blist));
    blist.m_lbox.clear();
}

BoxList&
BoxList::removeEmpty()
{
    m_lbox.erase(std::remove_if(m_lbox.begin(), m_lbox.end(),
                                [](const Box& x) { return x.isEmpty(); }),
                 m_lbox.end());
    return *this;
}

BoxList
intersect (const BoxList& bl,
           const Box&     b)
{
    BL_ASSERT(bl.ixType() == b.ixType());
    BoxList newbl(bl);
    newbl.intersect(b);
    return newbl;
}

BoxList
refine (const BoxList& bl,
        int            ratio)
{
    BoxList nbl(bl);
    nbl.refine(ratio);
    return nbl;
}

BoxList
coarsen (const BoxList& bl,
         int            ratio)
{
    BoxList nbl(bl);
    nbl.coarsen(ratio);
    return nbl;
}

BoxList
accrete (const BoxList& bl,
         int            sz)
{
    BoxList nbl(bl);
    nbl.accrete(sz);
    return nbl;
}

BoxList
removeOverlap (const BoxList& bl)
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
    m_lbox(),
    btype(IndexType::TheCellType())
{}

BoxList::BoxList (const Box& bx)
    : 
    m_lbox(1,bx),
    btype(bx.ixType())
{
}

BoxList::BoxList (IndexType _btype)
    :
    m_lbox(),
    btype(_btype)
{}

BoxList::BoxList (const BoxArray &ba)
    :
    m_lbox(std::move(ba.boxList().data())),
    btype(ba.ixType())
{
}

BoxList::BoxList (Vector<Box>&& bxs)
    : m_lbox(std::move(bxs)),
      btype(IndexType::TheCellType())
{
    if (m_lbox.size() > 0) {
        btype = m_lbox[0].ixType();
    }
}

BoxList::BoxList(const Box& bx, const IntVect& tilesize)
    : btype(bx.ixType())
{
    int ntiles = 1;
    IntVect nt;
    for (int d=0; d<AMREX_SPACEDIM; d++) {
	nt[d] = (bx.length(d)+tilesize[d]-1)/tilesize[d];
	ntiles *= nt[d];
    }

    IntVect small, big, ijk;  // note that the initial values are all zero.
    ijk[0] = -1;
    for (int t=0; t<ntiles; ++t) {
	for (int d=0; d<AMREX_SPACEDIM; d++) {
	    if (ijk[d]<nt[d]-1) {
		ijk[d]++;
		break;
	    } else {
		ijk[d] = 0;
	    }
	}

	for (int d=0; d<AMREX_SPACEDIM; d++) {
	    small[d] = ijk[d]*tilesize[d];
	    big[d] = std::min(small[d]+tilesize[d]-1, bx.length(d)-1);
	}

	Box tbx(small, big, btype);
	tbx.shift(bx.smallEnd());
	push_back(tbx);
    }
}

BoxList::BoxList (const Box& bx, int nboxes)
    : btype(bx.ixType())
{
    AMREX_ASSERT(nboxes > 0);
    AMREX_ASSERT(bx.numPts() >= nboxes);

    m_lbox.resize(nboxes);
    chop_boxes(m_lbox.data(), bx, nboxes);
}

BoxList::BoxList (const Box& bx, int nboxes, Direction dir)
    : btype(bx.ixType())
{
    int idir = static_cast<int>(dir);

    AMREX_ASSERT(nboxes > 0);
    AMREX_ASSERT(bx.length(idir) >= nboxes);

    m_lbox.resize(nboxes);
    chop_boxes_dir(m_lbox.data(), bx, nboxes, idir);
}

bool
BoxList::ok () const noexcept
{
    for (const auto& b : *this) {
        if (!b.ok()) return false;
    }
    return true;
}

bool
BoxList::isDisjoint () const
{
    return BoxArray(*this).isDisjoint();
}

bool
BoxList::contains (const BoxList&  bl) const
{
    if (isEmpty() || bl.isEmpty()) return false;

    BL_ASSERT(ixType() == bl.ixType());

    BoxArray ba(*this);

    for (const Box& bx : bl) {
        if (!ba.contains(bx)) {
            return false;
        }
    }

    return true;
}

BoxList&
BoxList::intersect (const Box& b)
{
    BL_ASSERT(ixType() == b.ixType());

    for (Box& bx : m_lbox)
    {
        const Box& isect = bx & b;
        if (isect.ok())
        {
            bx = isect;
        }
        else
        {
            bx = Box();
        }
    }

    removeEmpty();

    return *this;
}

BoxList&
BoxList::intersect (const BoxList& bl)
{
    BL_ASSERT(ixType() == bl.ixType());
    *this = amrex::intersect(BoxArray{*this}, bl);
    return *this;
}

BoxList
complementIn (const Box&     b,
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
    BoxArray ba(bl);
    return complementIn(b, ba);
}

BoxList&
BoxList::complementIn (const Box& b,
                       BoxList&&  bl)
{
    BoxArray ba(std::move(bl));
    return complementIn(b, ba);
}

BoxList&
BoxList::complementIn (const Box& b, const BoxArray& ba)
{
    BL_PROFILE("BoxList::complementIn");

    if (ba.size() == 0)
    {
	clear();
	push_back(b);
    }
    else if (ba.size() == 1)
    {
	*this = amrex::boxDiff(b, ba[0]);
    }
    else
    {
        long npts_avgbox;
	Box mbox = ba.minimalBox(npts_avgbox);
        *this = amrex::boxDiff(b, mbox);
        auto mytyp = ixType();

	BoxList bl_mesh(mbox & b);

#if (AMREX_SPACEDIM == 1)
        Real s_avgbox = npts_avgbox;
#elif (AMREX_SPACEDIM == 2)
        Real s_avgbox = std::sqrt(npts_avgbox);
#elif (AMREX_SPACEDIM == 3)
        Real s_avgbox = std::cbrt(npts_avgbox);
#endif

        const int block_size = 4 * (static_cast<int>(std::ceil(s_avgbox/4.))*4);
	bl_mesh.maxSize(block_size);
	const int N = bl_mesh.size();

#ifdef _OPENMP
        bool start_omp_parallel = !omp_in_parallel();
        const int nthreads = omp_get_max_threads();
#else
        bool start_omp_parallel = false;
#endif

        if (start_omp_parallel)
        {
#ifdef _OPENMP
            Vector<BoxList> bl_priv(nthreads, BoxList(mytyp));
#pragma omp parallel
            {
                BoxList bl_tmp(mytyp);
                auto& vbox = bl_priv[omp_get_thread_num()].m_lbox;
#pragma omp for
                for (int i = 0; i < N; ++i)
                {
                    ba.complementIn(bl_tmp, bl_mesh.m_lbox[i]);
                    vbox.insert(std::end(vbox), std::begin(bl_tmp), std::end(bl_tmp));
                }
            }
            for (auto& bl : bl_priv) {
                m_lbox.insert(std::end(m_lbox), std::begin(bl), std::end(bl));
            }
#else
            amrex::Abort("BoxList::complementIn: how did this happen");
#endif
        }
        else
        {
            BoxList bl_tmp(mytyp);
            for (int i = 0; i < N; ++i)
            {
                ba.complementIn(bl_tmp, bl_mesh.m_lbox[i]);
                m_lbox.insert(std::end(m_lbox), std::begin(bl_tmp), std::end(bl_tmp));
            }
        }
    }

    return *this;
}

BoxList&
BoxList::refine (int ratio)
{
    for (auto& bx : m_lbox)
    {
        bx.refine(ratio);
    }
    return *this;
}

BoxList&
BoxList::refine (const IntVect& ratio)
{
    for (auto& bx : m_lbox)
    {
        bx.refine(ratio);
    }
    return *this;
}

BoxList&
BoxList::coarsen (int ratio)
{
    for (auto& bx : m_lbox)
    {
        bx.coarsen(ratio);
    }
    return *this;
}

BoxList&
BoxList::coarsen (const IntVect& ratio)
{
    for (auto& bx : m_lbox)
    {
        bx.coarsen(ratio);
    }
    return *this;
}

BoxList&
BoxList::accrete (int sz)
{
    for (auto& bx : m_lbox)
    {
        bx.grow(sz);
    }
    return *this;
}

BoxList&
BoxList::accrete (const IntVect& sz)
{
    for (auto& bx : m_lbox)
    {
        bx.grow(sz);
    }
    return *this;
}

BoxList&
BoxList::shift (int dir,
                int nzones)
{
    for (auto& bx : m_lbox)
    {
        bx.shift(dir, nzones);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (int dir,
                    int num_halfs)
{
    for (auto& bx : m_lbox)
    {
        bx.shiftHalf(dir, num_halfs);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (const IntVect& iv)
{
    for (auto& bx : m_lbox)
    {
        bx.shiftHalf(iv);
    }
    return *this;
}

//
// Returns a list of boxes defining the compliment of b2 in b1in.
//

BoxList
boxDiff (const Box& b1in,
         const Box& b2)
{
   BL_ASSERT(b1in.sameType(b2));  
   BoxList bl_diff(b1in.ixType());
   boxDiff(bl_diff,b1in,b2);
   return bl_diff;
}

void
boxDiff (BoxList& bl_diff, const Box& b1in, const Box& b2)
{
   AMREX_ASSERT(b1in.sameType(b2));

   bl_diff.clear();
   bl_diff.set(b2.ixType());

   if ( !b2.contains(b1in) )
   {
       Box b1(b1in);
       if ( !b1.intersects(b2) )
       {
           bl_diff.push_back(b1);
       }
       else
       {
           const int* b2lo = b2.loVect();
           const int* b2hi = b2.hiVect();

           for (int i = AMREX_SPACEDIM-1; i >= 0; i--)
           {
               const int* b1lo = b1.loVect();
               const int* b1hi = b1.hiVect();

               if ((b1lo[i] < b2lo[i]) && (b2lo[i] <= b1hi[i]))
               {
                   Box bn(b1);
                   bn.setSmall(i,b1lo[i]);
                   bn.setBig(i,b2lo[i]-1);
                   bl_diff.push_back(bn);
                   b1.setSmall(i,b2lo[i]);
               }
               if ((b1lo[i] <= b2hi[i]) && (b2hi[i] < b1hi[i]))
               {
                   Box bn(b1);
                   bn.setSmall(i,b2hi[i]+1);
                   bn.setBig(i,b1hi[i]);
                   bl_diff.push_back(bn);
                   b1.setBig(i,b2hi[i]);
               }
           }
       }
   }
}

int
BoxList::simplify (bool best)
{
    std::sort(m_lbox.begin(), m_lbox.end(), [](const Box& l, const Box& r) {
            return l.smallEnd() < r.smallEnd(); });

    return simplify_doit(best);
}

int
BoxList::simplify_doit (bool best)
{
    //
    // Try to merge adjacent boxes.
    //
    int count = 0, lo[AMREX_SPACEDIM], hi[AMREX_SPACEDIM];

    for (iterator bla = begin(), End = end(); bla != End; ++bla)
    {
        const int* alo   = bla->loVect();
        const int* ahi   = bla->hiVect();
        //
        // If we're not looking for the "best" we can do in one pass, we
        // limit how far afield we look for abutting boxes.  This greatly
        // speeds up this routine for large numbers of boxes.  It does not
        // do quite as good a job though as full brute force.
        //
        const int MaxCnt = (best ? size() : 100);

        iterator blb = bla + 1;
        for (int cnt = 0; blb != End && cnt < MaxCnt; ++cnt, ++blb)
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
            for (int i = 0; i < AMREX_SPACEDIM; i++)
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
                // Modify b and set a to empty
                //
                blb->setSmall(IntVect(lo));
                blb->setBig(IntVect(hi));
                *bla = Box();
                ++count;
                break;
            }
        }
    }

    removeEmpty();

    return count;
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
    Vector<Box> new_boxes;

    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        new_boxes.clear();
        for (auto& bx : m_lbox)
        {
            const IntVect& boxlen = bx.size();
            const int* len = boxlen.getVect();

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
                const int sz     = nlen/numblk;
                const int extra  = nlen%numblk;
                //
                // Number of cuts = number of blocks - 1.
                //
                for (int k = 0; k < numblk-1; k++)
                {
                    //
                    // Compute size of this chunk, expand by power of 2.
                    //
                    const int ksize = (k < extra ? sz+1 : sz) * ratio;
                    //
                    // Chop from high end.
                    //
                    const int pos = bx.bigEnd(i) - ksize + 1;

                    new_boxes.push_back(bx.chop(i,pos));
                }
            }
        }
        join(new_boxes);
    }
    return *this;
}

BoxList&
BoxList::maxSize (int chunk)
{
    return maxSize(IntVect(AMREX_D_DECL(chunk,chunk,chunk)));
}

BoxList&
BoxList::surroundingNodes () noexcept
{
    for (auto& bx : m_lbox)
    {
        bx.surroundingNodes();
    }
    return *this;
}

BoxList&
BoxList::surroundingNodes (int dir) noexcept
{
    for (auto& bx : m_lbox)
    {
        bx.surroundingNodes(dir);
    }
    return *this;
}

BoxList&
BoxList::enclosedCells () noexcept
{
    for (auto& bx : m_lbox)
    {
        bx.enclosedCells();
    }
    return *this;
}

BoxList&
BoxList::enclosedCells (int dir) noexcept
{
    for (auto& bx : m_lbox)
    {
        bx.enclosedCells(dir);
    }
    return *this;
}

BoxList&
BoxList::convert (IndexType typ) noexcept
{
    btype = typ;
    for (auto& bx : m_lbox)
    {
        bx.convert(typ);
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
        amrex::Error("operator<<(ostream&,BoxList&) failed");

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

}
