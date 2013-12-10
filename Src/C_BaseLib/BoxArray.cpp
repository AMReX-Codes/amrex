
#include <iostream>

#include <BLassert.H>
#include <BoxArray.H>
#include <ParallelDescriptor.H>

void
BoxArray::reserve (long _truesize)
{
    if (!m_ref.unique())
        uniqify();
    m_ref->m_abox.reserve(_truesize);
}

BoxArray::BoxArray ()
    :
    m_ref(new BoxArray::Ref)
{}

BoxArray::BoxArray (const Box& bx)
    :
    m_ref(new BoxArray::Ref(1))
{
    m_ref->m_abox[0] = bx;
}

BoxArray::Ref::Ref (const BoxList& bl)
{
    define(bl);
}

BoxArray::BoxArray (const BoxList& bl)
    :
    m_ref(new BoxArray::Ref(bl))
{}

BoxArray::Ref::Ref (std::istream& is)
{
    define(is);
}

BoxArray::Ref::Ref (size_t size)
    :
    m_abox(size)
{}

BoxArray::BoxArray (size_t size)
    :
    m_ref(new BoxArray::Ref(size))
{}

BoxArray::BoxArray (const Box* bxvec,
                    int        nbox)
    :
    m_ref(new BoxArray::Ref(nbox))
{
    for (int i = 0; i < nbox; i++)
        m_ref->m_abox[i] = *bxvec++;
}

BoxArray::Ref::Ref (const Ref& rhs)
    :
    m_abox(rhs.m_abox)
{}

void
BoxArray::uniqify ()
{
    m_ref = new Ref(*m_ref);
}

void
BoxArray::clear ()
{
    if (!m_ref.unique())
        uniqify();
    m_ref->m_abox.clear();
}

void
BoxArray::resize (int len)
{
    if (!m_ref.unique())
        uniqify();
    m_ref->m_abox.resize(len);
}

void
BoxArray::set (int        i,
               const Box& ibox)
{
    if (!m_ref.unique())
        uniqify();
    m_ref->m_abox.set(i, ibox);
}

//
// Moved out of Utility.H
//
#define BL_IGNORE_MAX 100000

void
BoxArray::readFrom (std::istream& is)
{
    BL_ASSERT(size() == 0);
    if (!m_ref.unique())
        uniqify();
    m_ref->define(is);
}

void
BoxArray::Ref::define (std::istream& is)
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    BL_ASSERT(m_abox.size() == 0);
    int           maxbox;
    unsigned long hash;
    is.ignore(BL_IGNORE_MAX, '(') >> maxbox >> hash;
    m_abox.resize(maxbox);
    for (Array<Box>::iterator it = m_abox.begin(), End = m_abox.end(); it != End; ++it)
        is >> *it;
    is.ignore(BL_IGNORE_MAX, ')');
    if (is.fail())
        BoxLib::Error("BoxArray::define(istream&) failed");
}

void
BoxArray::define (const BoxList& bl)
{
    BL_ASSERT(size() == 0);
    if (!m_ref.unique())
        uniqify();
    m_ref->define(bl);
}

void
BoxArray::Ref::define (const BoxList& bl)
{
    BL_ASSERT(m_abox.size() == 0);
    m_abox.resize(bl.size());
    int count = 0;
    for (BoxList::const_iterator bli = bl.begin(), End = bl.end(); bli != End; ++bli)
    {
        m_abox.get(count++) = *bli;
    }
}

void
BoxArray::define (const BoxArray& bs)
{
    BL_ASSERT(size() == 0);
    m_ref = bs.m_ref;
}

BoxArray::~BoxArray () {}

bool
BoxArray::operator== (const BoxArray& rhs) const
{
    return m_ref == rhs.m_ref || m_ref->m_abox == rhs.m_ref->m_abox;
}

bool
BoxArray::operator!= (const BoxArray& rhs) const
{
    return !operator==(rhs);
}

BoxArray&
BoxArray::refine (int refinement_ratio)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].refine(refinement_ratio);
    return *this;
}

BoxArray&
BoxArray::refine (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].refine(iv);
    return *this;
}

BoxArray&
BoxArray::shift (int dir,
                 int nzones)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].shift(dir, nzones);
    return *this;
}


BoxArray&
BoxArray::shift (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].shift(iv);
    return *this;
}

BoxArray&
BoxArray::shiftHalf (int dir,
                     int num_halfs)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].shiftHalf(dir, num_halfs);
    return *this;
}

BoxArray&
BoxArray::shiftHalf (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].shiftHalf(iv);
    return *this;
}

BoxArray&
BoxArray::coarsen (int refinement_ratio)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].coarsen(refinement_ratio);
    return *this;
}

BoxArray&
BoxArray::coarsen (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].coarsen(iv);
    return *this;
}

BoxArray&
BoxArray::grow (int n)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].grow(n);
    return *this;
}

BoxArray&
BoxArray::grow (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].grow(iv);
    return *this;
}

BoxArray&
BoxArray::grow (int dir,
                int n_cell)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].grow(dir, n_cell);
    return *this;
}

bool
BoxArray::contains (const IntVect& iv) const
{
    if (size() > 0)
    {
        std::vector< std::pair<int,Box> > isects;
        intersections(Box(iv,iv,get(0).ixType()),isects);
        for (int i = 0, N = isects.size(); i < N; i++)
            if (get(isects[i].first).contains(iv))
                return true;
    }
    return false;
}

bool
BoxArray::contains (const Box& b) const
{
    if (size() > 0)
    {
        BL_ASSERT(get(0).sameType(b));

        std::vector< std::pair<int,Box> > isects;

        intersections(b,isects);

        if (isects.size() > 0)
        {
            BoxList bl(b.ixType());
            for (int i = 0, N = isects.size(); i < N; i++)
                bl.push_back(isects[i].second);
            BoxList blnew = BoxLib::complementIn(b, bl);
            return blnew.size() == 0;
        }
    }

    return false;
}

bool
BoxArray::contains (const BoxArray& bl) const
{
    if (size() == 0) return false;

    if (!minimalBox().contains(bl.minimalBox())) return false;

    for (BoxArray::const_iterator it = bl.begin(), End = bl.end(); it != End; ++it)
        if (!contains(*it))
            return false;

    return true;
}

BoxArray&
BoxArray::surroundingNodes ()
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].surroundingNodes();
    return *this;
}

BoxArray&
BoxArray::surroundingNodes (int dir)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].surroundingNodes(dir);
    return *this;
}

BoxArray&
BoxArray::enclosedCells ()
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].enclosedCells();
    return *this;
}

BoxArray&
BoxArray::enclosedCells (int dir)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].enclosedCells(dir);
    return *this;
}

BoxArray&
BoxArray::convert (IndexType typ)
{
    if (!m_ref.unique())
        uniqify();
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].convert(typ);
    return *this;
}

BoxArray&
BoxArray::convert (Box (*fp)(const Box&))
{
    BL_ASSERT(!(fp == 0));

    if (!m_ref.unique())
        uniqify();
    for (Array<Box>::iterator it = m_ref->m_abox.begin(), End = m_ref->m_abox.end(); it != End; ++it)
        *it = (*fp)(*it);
    return *this;
}

std::ostream&
BoxArray::writeOn (std::ostream& os) const
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    os << '(' << size() << ' ' << 0 << '\n';

    for (BoxArray::const_iterator it = begin(), End = end(); it != End; ++it)
        os << *it << '\n';

    os << ')';

    if (os.fail())
        BoxLib::Error("BoxArray::writeOn(ostream&) failed");

    return os;
}

bool
BoxArray::isDisjoint () const
{
    for (BoxArray::const_iterator it = begin(), End = end(); it != End; ++it)
    {
        std::vector< std::pair<int,Box> > isects;

        intersections(*it,isects);

        if ( !(isects.size() == 1 && isects[0].second == *it) )
        {
            return false;
        }
    }

    return true;
}

bool
BoxArray::ok () const
{
    bool isok = true;

    if (size() > 0)
    {
        const Box& bx0 = m_ref->m_abox[0];

        if (size() == 1) isok = bx0.ok();

        for (int i = 1, N = size(); i < N && isok; i++)
        {
            const Box& bxi = m_ref->m_abox[i];
            isok = bxi.ok() && bxi.sameType(bx0);
        }
    }

    return isok;
}

long
BoxArray::numPts () const
{
    long result = 0;
    for (BoxArray::const_iterator it = begin(), End = end(); it != End; ++it)
    {
        result += it->numPts();
    }
    return result;
}

double
BoxArray::d_numPts () const
{
    double result = 0;
    for (BoxArray::const_iterator it = begin(), End = end(); it != End; ++it)
    {
        result += it->d_numPts();
    }
    return result;
}


BoxArray&
BoxArray::maxSize (const IntVect& block_size)
{
    BoxList blst(*this);
    blst.maxSize(block_size);
    clear();
    m_ref->m_abox.resize(blst.size());
    BoxList::iterator bli = blst.begin(), End = blst.end();
    for (int i = 0; bli != End; ++bli)
        set(i++, *bli);
    return *this;
}

BoxArray&
BoxArray::maxSize (int block_size)
{
    return maxSize(IntVect(D_DECL(block_size,block_size,block_size)));
}

std::ostream&
operator<< (std::ostream&   os,
            const BoxArray& ba)
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    os << "(BoxArray maxbox("
       << ba.size()
       << ")\n       m_ref->m_hash_sig("
       << 0
       << ")\n       ";

    for (BoxArray::const_iterator it = ba.begin(), End = ba.end(); it != End; ++it)
        os << *it << ' ';

    os << ")\n";

    if (os.fail())
        BoxLib::Error("operator<<(ostream& os,const BoxArray&) failed");

    return os;
}

BoxList
BoxArray::boxList () const
{
    if ( size() == 0 ) return BoxList();
    BoxList newb(get(0).ixType());
    for (BoxArray::const_iterator it = begin(), End = end(); it != End; ++it)
        newb.push_back(*it);
    return newb;
}

Box
BoxArray::minimalBox () const
{
    Box minbox;
    if (size() > 0)
    {
        minbox = m_ref->m_abox.get(0);
        for (BoxArray::const_iterator it = begin(), End = end(); it != End; ++it)
            minbox.minBox(*it);
    }
    return minbox;
}

BoxArray
BoxLib::boxComplement (const Box& b1in,
		       const Box& b2)
{
    return BoxArray(BoxLib::boxDiff(b1in, b2));
}

BoxArray
BoxLib::complementIn (const Box&      b,
		      const BoxArray& ba)
{
    return BoxArray(BoxLib::complementIn(b,ba.boxList()));
}

BoxArray
BoxLib::intersect (const BoxArray& ba,
		   const Box&      b)
{
    std::vector< std::pair<int,Box> > isects;

    ba.intersections(b,isects);

    BoxArray r(isects.size());

    for (int i = 0, N = isects.size(); i < N; i++)
    {
        r.set(i, isects[i].second);
    }

    return r;
}

BoxArray
BoxLib::intersect (const BoxArray& lhs,
		   const BoxArray& rhs)
{
    if (lhs.size() == 0 || rhs.size() == 0) return BoxArray();
    BoxList bl(lhs[0].ixType());
    for (BoxArray::const_iterator it = lhs.begin(), End = lhs.end(); it != End; ++it)
    {
        BoxArray ba  = BoxLib::intersect(rhs, *it);
        BoxList  tmp = ba.boxList();
        bl.catenate(tmp);
    }
    return BoxArray(bl);
}

BoxList
BoxLib::GetBndryCells (const BoxArray& ba,
                       int             ngrow)
{
    BL_ASSERT(ba.ok());
    BL_ASSERT(ba.size() > 0);
    //
    // First get list of all ghost cells.
    //
    const IndexType btype = ba[0].ixType();

    BoxList bcells = ba.boxList(), gcells(btype), leftover(btype);

    bcells.simplify();

    BoxArray tba(bcells);

    bcells.clear();

    for (BoxArray::const_iterator it = tba.begin(), End = tba.end(); it != End; ++it)
    {
        bcells = BoxLib::boxDiff(BoxLib::grow(*it,ngrow),*it);

	gcells.catenate(bcells);
    }
    //
    // Now strip out intersections with original BoxArray.
    //
    std::vector< std::pair<int,Box> > isects;

    for (BoxList::const_iterator it = gcells.begin(), End = gcells.end(); it != End; ++it)
    {
        tba.intersections(*it,isects);

        if (isects.empty())
        {
            bcells.push_back(*it);
        }
        else
        {
            //
            // Collect all the intersection pieces.
            //
            BoxList pieces(btype);
            for (int i = 0, N = isects.size(); i < N; i++)
                pieces.push_back(isects[i].second);
            leftover = BoxLib::complementIn(*it,pieces);
            bcells.catenate(leftover);
        }
    }
    //
    // Now strip out overlaps.
    //
    gcells.clear();

    gcells = BoxLib::removeOverlap(bcells);

    bcells.clear();

    gcells.simplify();

    return gcells;
}

std::vector< std::pair<int,Box> >
BoxArray::intersections (const Box& bx) const
{
    std::vector< std::pair<int,Box> > isects;
    intersections(bx,isects);
    return isects;
}

void
BoxArray::clear_hash_bin () const
{
    m_ref->hash.clear();
}

void
BoxArray::intersections (const Box&                         bx,
                         std::vector< std::pair<int,Box> >& isects) const
{

#ifdef _OPENMP
#pragma omp critical(intersections_lock)
#endif
    {
        if (!m_ref->hash.isAllocated() && size() > 0)
        {
            BL_ASSERT(bx.sameType(get(0)));
            //
            // Calculate the bounding box & maximum extent of the boxes.
            //
            IntVect maxext(D_DECL(0,0,0));

            Box boundingbox = get(0);

            for (BoxArray::const_iterator it = begin(), End = end(); it != End; ++it)
            {
                boundingbox.minBox(*it);
                maxext = BoxLib::max(maxext, it->size());

            }

            m_ref->crsn = maxext;

            boundingbox.coarsen(maxext);

            m_ref->hash.resize(boundingbox, 1);

            long total = boundingbox.numPts()*sizeof(std::vector<int>);

            for (int i = 0, N = size(); i < N; i++)
            {
                m_ref->hash(BoxLib::coarsen(get(i).smallEnd(),maxext)).push_back(i);
            }
            //
            // Now compress out excess capacity in vectors.
            //
            for (IntVect iv = boundingbox.smallEnd(), End = boundingbox.bigEnd(); iv <= End; boundingbox.next(iv))
            {
                std::vector<int>& v = m_ref->hash(iv);

                if (!v.empty())
                {
                    const int N = v.size();

                    std::vector<int> tmp(N);

                    for (int i = 0; i < N; i++) tmp[i] = v[i];

                    v.swap(tmp);

                    total += N * sizeof(int);
                }
            }

            if (false && ParallelDescriptor::IOProcessor())
                std::cout << "*** BoxArray::intersections(): bytes in box hash: " << total << '\n';
        }
    }

    isects.resize(0);

    isects.reserve(27);

    if (m_ref->hash.isAllocated())
    {
        BL_ASSERT(bx.sameType(get(0)));

        Box     cbx = BoxLib::coarsen(bx, m_ref->crsn);
        IntVect sm  = BoxLib::max(cbx.smallEnd()-1, m_ref->hash.box().smallEnd());
        IntVect bg  = BoxLib::min(cbx.bigEnd(),     m_ref->hash.box().bigEnd());

        cbx = Box(sm,bg,bx.ixType());

        for (IntVect iv = cbx.smallEnd(), End = cbx.bigEnd(); iv <= End; cbx.next(iv))
        {
            std::vector<int>& v = m_ref->hash(iv);

            for (int i = 0, N = v.size(); i < N; i++)
            {
                const Box isect = bx & get(v[i]);

                if (isect.ok())
                {
                    isects.push_back(std::pair<int,Box>(v[i],isect));
                }
            }
        }
    }
}

//
// Currently this assumes your Boxes are cell-centered.
//

void
BoxArray::removeOverlap ()
{
    if (!m_ref.unique()) uniqify();

    BoxList bl;

    const Box EmptyBox;

    std::vector< std::pair<int,Box> > isects;
    //
    // Note that "size()" can increase in this loop!!!
    //
    for (int i = 0; i < size(); i++)
    {
        if (m_ref->m_abox[i].ok())
        {
            intersections(m_ref->m_abox[i],isects);

            for (int j = 0, N = isects.size(); j < N; j++)
            {
                if (isects[j].first == i) continue;

                Box& bx = m_ref->m_abox[isects[j].first];

                bl = BoxLib::boxDiff(bx, isects[j].second);

                bx = EmptyBox;

                for (BoxList::const_iterator it = bl.begin(), End = bl.end(); it != End; ++it)
                {
                    m_ref->m_abox.push_back(*it);

                    m_ref->hash(BoxLib::coarsen(it->smallEnd(),m_ref->crsn)).push_back(size()-1);
                }
            }
        }
    }
    //
    // We now have "holes" in our BoxArray. Make us good.
    //
    bl.clear();

    const Box& bb = m_ref->hash.box();

    for (IntVect iv = bb.smallEnd(), End = bb.bigEnd(); iv <= End; bb.next(iv))
    {
        std::vector<int>& v = m_ref->hash(iv);

        for (int i = 0, N = v.size(); i < N; i++)
        {
            if (m_ref->m_abox[v[i]].ok())
            {
                bl.push_back(m_ref->m_abox[v[i]]);
            }
        }
    }

    bl.simplify();

    BoxArray nba(bl);

    bl.clear();

    *this = nba;

    BL_ASSERT(isDisjoint());
}
