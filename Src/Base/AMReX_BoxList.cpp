
#include <AMReX_Print.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <algorithm>
#include <iostream>
#include <cmath>
#include <map>
#include <utility>

namespace amrex {

namespace {

void chop_boxes (Box* bxv, const Box& bx, int nboxes)
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

void chop_boxes_dir (Box* bxv, const Box& bx, int nboxes, int idir)
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
    BL_ASSERT(blist.size() == 0 || ixType() == blist.ixType());
    m_lbox.insert(std::end(m_lbox), std::begin(blist), std::end(blist));
}

void
BoxList::join (const Vector<Box>& barr)
{
    BL_ASSERT(barr.empty() || ixType() == barr[0].ixType());
    m_lbox.insert(std::end(m_lbox), std::begin(barr), std::end(barr));
}

void
BoxList::catenate (BoxList& blist)
{
    BL_ASSERT(blist.size() == 0 || ixType() == blist.ixType());
    m_lbox.insert(std::end(m_lbox), std::begin(blist), std::end(blist));
    blist.m_lbox.clear();
}

BoxList&
BoxList::removeEmpty()
{
    std::erase_if(m_lbox, [](const Box& x) { return x.isEmpty(); });
    return *this;
}

BoxList
intersect (const BoxList& bl, const Box& b)
{
    BL_ASSERT(bl.ixType() == b.ixType());
    BoxList newbl(bl);
    newbl.intersect(b);
    return newbl;
}

BoxList
refine (const BoxList& bl, int ratio)
{
    BoxList nbl(bl);
    nbl.refine(ratio);
    return nbl;
}

BoxList
coarsen (const BoxList& bl, int ratio)
{
    BoxList nbl(bl);
    nbl.coarsen(ratio);
    return nbl;
}

BoxList
accrete (const BoxList& bl, int sz)
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
    if (!m_lbox.empty()) {
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

    IntVect sml, big, ijk;  // note that the initial values are all zero.
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
            sml[d] = ijk[d]*tilesize[d];
            big[d] = std::min(sml[d]+tilesize[d]-1, bx.length(d)-1);
        }

        Box tbx(sml, big, btype);
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
    return std::ranges::all_of(*this,
                               [] (Box const& b) { return b.ok(); });
}

bool
BoxList::isDisjoint () const
{
    if (this->size() <= 1) {
        return true;
    } else {
        return BoxArray(*this).isDisjoint();
    }
}

bool
BoxList::contains (const BoxList& bl) const
{
    if (isEmpty() || bl.isEmpty()) { return false; }

    BL_ASSERT(ixType() == bl.ixType());

    BoxArray ba(*this);

    return std::ranges::all_of(bl,
                               [&ba] (Box const& b) { return ba.contains(b); });
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
complementIn (const Box& b, const BoxList& bl)
{
    BL_ASSERT(bl.ixType() == b.ixType());
    BoxList newb(b.ixType());
    newb.complementIn(b,bl);
    return newb;
}

BoxList&
BoxList::complementIn (const Box& b, const BoxList& bl)
{
    BoxArray ba(bl);
    return complementIn(b, ba);
}

BoxList&
BoxList::complementIn (const Box& b, BoxList&& bl)
{
    BoxArray ba(std::move(bl));
    return complementIn(b, ba);
}

BoxList&
BoxList::complementIn (const Box& b, const BoxArray& ba)
{
    BL_PROFILE("BoxList::complementIn");

    if (ba.empty())
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
        Long npts_avgbox;
        Box mbox = ba.minimalBox(npts_avgbox);
        *this = amrex::boxDiff(b, mbox);
        auto mytyp = ixType();

        BoxList bl_mesh(mbox & b);

#if (AMREX_SPACEDIM == 1)
        Real s_avgbox = static_cast<Real>(npts_avgbox);
#elif (AMREX_SPACEDIM == 2)
        Real s_avgbox = static_cast<Real>(std::sqrt(npts_avgbox));
#elif (AMREX_SPACEDIM == 3)
        Real s_avgbox = static_cast<Real>(std::cbrt(npts_avgbox));
#endif

        const int block_size = 4 * std::max(1,static_cast<int>(std::ceil(s_avgbox/4.))*4);
        bl_mesh.maxSize(block_size);
        const int N = static_cast<int>(bl_mesh.size());

#ifdef AMREX_USE_OMP
        bool start_omp_parallel = !omp_in_parallel();
        const int nthreads = omp_get_max_threads();
#else
        bool start_omp_parallel = false;
#endif

        if (start_omp_parallel)
        {
#ifdef AMREX_USE_OMP
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
BoxList::parallelComplementIn (const Box& b, const BoxList& bl)
{
    return parallelComplementIn(b, BoxArray(bl));
}

BoxList&
BoxList::parallelComplementIn (const Box& b, BoxList&& bl)
{
    return parallelComplementIn(b, BoxArray(std::move(bl)));
}

BoxList&
BoxList::parallelComplementIn (const Box& b, BoxArray const& ba)
{
    BL_PROFILE("BoxList::parallelComplementIn()");
#ifndef AMREX_USE_MPI
    return complementIn(b,ba);
#else
    if (ba.size() <= 8)
    {
        return complementIn(b,ba);
    }
    else
    {
        BL_PROFILE_VAR("BoxList::pci", boxlistpci);

        Long npts_avgbox;
        Box mbox = ba.minimalBox(npts_avgbox);
        *this = amrex::boxDiff(b, mbox);
        auto mytyp = ixType();

        BoxList bl_mesh(mbox & b);

#if (AMREX_SPACEDIM == 1)
        Real s_avgbox = static_cast<Real>(npts_avgbox);
#elif (AMREX_SPACEDIM == 2)
        Real s_avgbox = static_cast<Real>(std::sqrt(npts_avgbox));
#elif (AMREX_SPACEDIM == 3)
        Real s_avgbox = static_cast<Real>(std::cbrt(npts_avgbox));
#endif

        const int block_size = 4 * std::max(1,static_cast<int>(std::ceil(s_avgbox/4.))*4);
        bl_mesh.maxSize(block_size);
        const int N = static_cast<int>(bl_mesh.size());

        const int nprocs = ParallelContext::NProcsSub();
        const int myproc = ParallelContext::MyProcSub();
        const int navg = N / nprocs;
        const int nextra = N - navg*nprocs;
        const int ilo = (myproc < nextra) ? myproc*(navg+1) : myproc*navg+nextra;
        const int ihi = (myproc < nextra) ? ilo+navg+1-1 : ilo+navg-1;

        Vector<Box> local_boxes;

#ifdef AMREX_USE_OMP
        bool start_omp_parallel = !omp_in_parallel();
        const int nthreads = omp_get_max_threads();
#else
        bool start_omp_parallel = false;
#endif

        if (start_omp_parallel)
        {
#ifdef AMREX_USE_OMP
            Vector<BoxList> bl_priv(nthreads, BoxList(mytyp));
            int ntot = 0;
#pragma omp parallel reduction(+:ntot)
            {
                BoxList bl_tmp(mytyp);
                auto& vbox = bl_priv[omp_get_thread_num()].m_lbox;
#pragma omp for
                for (int i = ilo; i <= ihi; ++i)
                {
                    ba.complementIn(bl_tmp, bl_mesh.m_lbox[i]);
                    vbox.insert(std::end(vbox), std::begin(bl_tmp), std::end(bl_tmp));
                }
                ntot += static_cast<int>(bl_tmp.size());
            }
            local_boxes.reserve(ntot);
            for (auto& bl : bl_priv) {
                local_boxes.insert(std::end(local_boxes), std::begin(bl), std::end(bl));
            }
#else
            amrex::Abort("BoxList::complementIn: how did this happen");
#endif
        }
        else
        {
            BoxList bl_tmp(mytyp);
            for (int i = ilo; i <= ihi; ++i)
            {
                ba.complementIn(bl_tmp, bl_mesh.m_lbox[i]);
                local_boxes.insert(std::end(local_boxes), std::begin(bl_tmp), std::end(bl_tmp));
            }
        }

        amrex::AllGatherBoxes(local_boxes, static_cast<int>(this->size()));
        local_boxes.insert(std::end(local_boxes), std::begin(m_lbox), std::end(m_lbox));
        std::swap(m_lbox, local_boxes);

        return *this;
    }
#endif
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
BoxList::shift (int dir, int nzones)
{
    for (auto& bx : m_lbox)
    {
        bx.shift(dir, nzones);
    }
    return *this;
}

BoxList&
BoxList::shiftHalf (int dir, int num_halfs)
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
boxDiff (const Box& b1in, const Box& b2)
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
    std::ranges::sort(m_lbox, [](const Box& l, const Box& r) {
            return l.smallEnd() < r.smallEnd(); });

    //
    // If we're not looking for the "best" we can do in one pass, we
    // limit how far afield we look for abutting boxes.  This greatly
    // speeds up this routine for large numbers of boxes.  It does not
    // do quite as good a job though as full brute force.
    //
    int depth = best ? static_cast<int>(size()) : 100;
    return simplify_doit(depth);
}

int
BoxList::ordered_simplify ()
{
    int count;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        count = simplify_doit(1);
    }
    return count;
}

int
BoxList::simplify_doit (int depth)
{
    //
    // Try to merge adjacent boxes.
    //
    int count = 0, lo[AMREX_SPACEDIM], hi[AMREX_SPACEDIM];

    for (auto bla = begin(), End = end(); bla != End; ++bla)
    {
        const int* alo   = bla->loVect();
        const int* ahi   = bla->hiVect();

        auto blb = bla + 1;
        for (int cnt = 0; blb != End && cnt < depth; ++cnt, ++blb)
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
        auto bli = cbegin(), End = cend();
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
    for (auto const& bx : m_lbox) {
        const IntVect boxlen = amrex::enclosedCells(bx).size();
        const IntVect boxlo = bx.smallEnd();
        IntVect ratio{1}, numblk{1}, extra{0};
        IntVect sz = boxlen;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (boxlen[idim] > chunk[idim]) {
                int bs    = chunk[idim];
                int nlen  = boxlen[idim];
                while ((bs%2 == 0) && (nlen%2 == 0)) {
                    ratio[idim] *= 2;
                    bs    /= 2;
                    nlen  /= 2;
                }
                numblk[idim] = (nlen+bs-1)/bs;
                sz[idim] = nlen/numblk[idim];
                extra[idim] = nlen - sz[idim]*numblk[idim];
            }
        }
        if (numblk == 1) {
            new_boxes.push_back(bx);
        } else {
#if (AMREX_SPACEDIM == 3)
            for (int k = 0; k < numblk[2]; ++k) {
                int klo = (k < extra[2]) ? k*(sz[2]+1)*ratio[2] : (k*sz[2]+extra[2])*ratio[2];
                int khi = (k < extra[2]) ? klo+(sz[2]+1)*ratio[2]-1 : klo+sz[2]*ratio[2]-1;
                klo += boxlo[2];
                khi += boxlo[2];
#endif
#if (AMREX_SPACEDIM >= 2)
                for (int j = 0; j < numblk[1]; ++j) {
                    int jlo = (j < extra[1]) ? j*(sz[1]+1)*ratio[1] : (j*sz[1]+extra[1])*ratio[1];
                    int jhi = (j < extra[1]) ? jlo+(sz[1]+1)*ratio[1]-1 : jlo+sz[1]*ratio[1]-1;
                    jlo += boxlo[1];
                    jhi += boxlo[1];
#endif
                    for (int i = 0; i < numblk[0]; ++i) {
                        int ilo = (i < extra[0]) ? i*(sz[0]+1)*ratio[0] : (i*sz[0]+extra[0])*ratio[0];
                        int ihi = (i < extra[0]) ? ilo+(sz[0]+1)*ratio[0]-1 : ilo+sz[0]*ratio[0]-1;
                        ilo += boxlo[0];
                        ihi += boxlo[0];
                        new_boxes.push_back(Box(IntVect(AMREX_D_DECL(ilo,jlo,klo)),
                                                IntVect(AMREX_D_DECL(ihi,jhi,khi))).
                                            convert(ixType()));
            AMREX_D_TERM(},},})
        }
    }
    std::swap(new_boxes, m_lbox);

    return *this;
}

BoxList&
BoxList::maxSize (int chunk)
{
    return maxSize(IntVect(AMREX_D_DECL(chunk,chunk,chunk)));
}

BoxList&
BoxList::mergeAlongDir (int dir)
{
    AMREX_ASSERT(dir >= 0 && dir < AMREX_SPACEDIM);
    AMREX_ASSERT(btype.cellCentered());
    AMREX_ASSERT(isDisjoint());

    removeEmpty();

    // Overlap of a and b in the directions other than dir, using a's extent in dir.
    auto transverse_overlap = [dir] (Box const& a, Box const& b) -> Box
    {
        Box bp = b;
        bp.setSmall(dir, a.smallEnd(dir));
        bp.setBig(dir, a.bigEnd(dir));
        return a & bp;
    };

    // Whether a and b overlap in the directions other than dir.  Cheaper
    // than transverse_overlap(a,b).ok() for the pair search.
    auto transverse_touch = [dir] (Box const& a, Box const& b) -> bool
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (d != dir && (a.smallEnd(d) > b.bigEnd(d) || b.smallEnd(d) > a.bigEnd(d))) {
                return false;
            }
        }
        return true;
    };

    std::ranges::sort(m_lbox, [dir] (Box const& a, Box const& b) {
            return a.smallEnd(dir) < b.smallEnd(dir); });

    // Boxes whose hi face in dir may still be extended, keyed by bigEnd(dir).
    std::map<int, Vector<Box>> pending;
    // Boxes whose hi face in dir is known to be free.
    Vector<Box> done;

    // Sweep over layers of boxes with the same smallEnd(dir).
    Long const N = static_cast<Long>(m_lbox.size());
    Long ibeg = 0;
    while (ibeg < N) {
        int const p = m_lbox[ibeg].smallEnd(dir);
        Long iend = ibeg+1;
        while (iend < N && m_lbox[iend].smallEnd(dir) == p) { ++iend; }

        // Pieces of the boxes ending at p-1 (front) and starting at p (newb).
        Vector<Vector<Box>> fpieces;
        if (auto it = pending.find(p-1); it != pending.end()) {
            fpieces.reserve(it->second.size());
            for (auto const& b : it->second) { fpieces.push_back(Vector<Box>{b}); }
            pending.erase(it);
        }
        Vector<Vector<Box>> npieces;
        npieces.reserve(iend-ibeg);
        for (Long i = ibeg; i < iend; ++i) { npieces.push_back(Vector<Box>{m_lbox[i]}); }

        if (!fpieces.empty()) {
            int const nf = static_cast<int>(fpieces.size());
            int const nn = static_cast<int>(npieces.size());

            // Find pairs of original boxes that touch across the p-1/p interface.
            Vector<std::pair<int,int>> pairs;
#ifdef AMREX_USE_OMP
            if (static_cast<Long>(nf)*static_cast<Long>(nn) >= 4096 && !omp_in_parallel()) {
                Vector<Vector<std::pair<int,int>>> pairs_priv(omp_get_max_threads());
#pragma omp parallel
                {
                    auto& pp = pairs_priv[omp_get_thread_num()];
#pragma omp for schedule(static)
                    for (int i = 0; i < nf; ++i) {
                        for (int j = 0; j < nn; ++j) {
                            if (transverse_touch(fpieces[i][0], npieces[j][0])) {
                                pp.emplace_back(i,j);
                            }
                        }
                    }
                }
                for (auto const& pp : pairs_priv) {
                    pairs.insert(pairs.end(), pp.begin(), pp.end());
                }
            } else
#endif
            {
                for (int i = 0; i < nf; ++i) {
                    for (int j = 0; j < nn; ++j) {
                        if (transverse_touch(fpieces[i][0], npieces[j][0])) {
                            pairs.emplace_back(i,j);
                        }
                    }
                }
            }

            // Join the touching parts; split off the rest in the other directions.
            BoxList tmp(btype);
            for (auto const& [i,j] : pairs) {
                auto& A = fpieces[i];
                auto& B = npieces[j];
                for (std::size_t ia = 0; ia < A.size(); ++ia) {
                    for (std::size_t ib = 0; ib < B.size(); ++ib) {
                        if (!B[ib].ok()) { continue; }
                        Box const aT = transverse_overlap(A[ia], B[ib]);
                        if (aT.ok()) {
                            Box const a = A[ia];
                            Box const b = B[ib];
                            Box const bT = transverse_overlap(b, aT);
                            Box merged = aT;
                            merged.setBig(dir, bT.bigEnd(dir));
                            pending[merged.bigEnd(dir)].push_back(merged);
                            boxDiff(tmp, a, aT);
                            A[ia] = Box();
                            A.insert(A.end(), tmp.begin(), tmp.end());
                            boxDiff(tmp, b, bT);
                            B[ib] = Box();
                            B.insert(B.end(), tmp.begin(), tmp.end());
                            break;
                        }
                    }
                }
                std::erase_if(A, [] (Box const& b) { return !b.ok(); });
                std::erase_if(B, [] (Box const& b) { return !b.ok(); });
            }

            for (auto const& A : fpieces) {
                done.insert(done.end(), A.begin(), A.end());
            }
        }

        for (auto const& B : npieces) {
            for (auto const& b : B) {
                pending[b.bigEnd(dir)].push_back(b);
            }
        }

        ibeg = iend;
    }

    m_lbox = std::move(done);
    for (auto const& kv : pending) {
        m_lbox.insert(m_lbox.end(), kv.second.begin(), kv.second.end());
    }

    while (simplify() > 0) {}

    return *this;
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
operator<< (std::ostream& os, const BoxList& blist)
{
    auto bli = blist.begin(), End = blist.end();
    os << "(BoxList " << blist.size() << ' ' << blist.ixType() << '\n';
    for (int count = 1; bli != End; ++bli, ++count)
    {
        os << count << " : " << *bli << '\n';
    }
    os << ')' << '\n';

    if (os.fail()) {
        amrex::Error("operator<<(ostream&,BoxList&) failed");
    }

    return os;
}

bool
BoxList::operator== (const BoxList& rhs) const
{
    if ( !(size() == rhs.size()) ) { return false; }

    auto liter = begin(), riter = rhs.begin(), End = end();
    for (; liter != End; ++liter, ++riter) {
        if ( !( *liter == *riter) ) {
            return false;
        }
    }
    return true;
}

void
BoxList::Bcast ()
{
    int nboxes = static_cast<int>(this->size());
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::Bcast(&nboxes, 1, IOProcNumber);
    if (ParallelDescriptor::MyProc() != IOProcNumber) {
        m_lbox.resize(nboxes);
    }
    ParallelDescriptor::Bcast(m_lbox.data(), nboxes, IOProcNumber);
}

}
