
#include <iostream>

#include <BLassert.H>
#include <BoxArray.H>
#include <ParallelDescriptor.H>

typedef std::map< IntVect,std::vector<int>,IntVect::Compare > BoxHashMapType;
typedef std::map< IntVect,std::vector<int>,IntVect::Compare >::iterator BoxHashMapIter;
typedef std::map< IntVect,std::vector<int>,IntVect::Compare >::const_iterator ConstBoxHashMapIter;

void
BoxArray::decrementCounters () const
{
    if (m_ref.linkCount() == 1)
    {
        clear_hash_bin();
    }
}

//
// Most heavily used BoxArray constructor.
//
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

BoxArray::BoxArray (const BoxList& bl)
    :
    m_ref(new BoxArray::Ref(bl))
{}

BoxArray::BoxArray (size_t n)
    :
    m_ref(new BoxArray::Ref(n))
{}

BoxArray::BoxArray (const Box* bxvec,
                    int        nbox)
    :
    m_ref(new BoxArray::Ref(nbox))
{
    for (int i = 0; i < nbox; i++)
        m_ref->m_abox[i] = *bxvec++;
}

//
// The copy constructor.
//
BoxArray::BoxArray (const BoxArray& rhs)
    :
    m_ref(rhs.m_ref)
{}

//
// The assignment operator.
//
BoxArray&
BoxArray::operator= (const BoxArray& rhs)
{
    decrementCounters();
    m_ref = rhs.m_ref;
    return *this;
}

void
BoxArray::uniqify ()
{
    BL_ASSERT(!m_ref.unique());
    m_ref = new BoxArray::Ref(*m_ref);
}

int
BoxArray::size () const
{
    return m_ref->m_abox.size();
}

bool
BoxArray::empty () const
{
    return m_ref->m_abox.empty();
}

void
BoxArray::clear ()
{
    if (!m_ref.unique())
        uniqify();

    if (m_ref->m_abox.capacity() > 0)
    {
        //
        // I want to forcefully clear out any memory in m_abox.
        // Normally a clear() would just reset some internal pointers.
        // I want to really let go of the memory so my byte count
        // statistics are more accurate.
        //
        Array<Box>().swap(m_ref->m_abox);
    }
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
    const int N = bl.size();
    m_abox.resize(N);
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
    decrementCounters();
    m_ref = bs.m_ref;
}

BoxArray::~BoxArray ()
{
    decrementCounters();
}

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
BoxArray::intersects (const Box& b) const
{
    std::vector< std::pair<int,Box> > isects;

    bool first_only = true;
    intersections(b,isects,first_only);

    return (isects.size() > 0) ;
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
    std::vector< std::pair<int,Box> > isects;

    for (BoxArray::const_iterator it = begin(), End = end(); it != End; ++it)
    {
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

long long
BoxArray::ll_numPts () const
{
    long long result = 0;
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
    const int N = blst.size();
    m_ref->m_abox.resize(N);
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

    tba.clear_hash_bin();
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
BoxArray::intersections (const Box& bx, bool first_only) const
{
    std::vector< std::pair<int,Box> > isects;
    intersections(bx,isects,first_only);
    return isects;
}

void
BoxArray::clear_hash_bin () const
{
    if (!m_ref->hash.empty())
    {
        m_ref->hash.clear();
    }
}

void
BoxArray::intersections (const Box&                         bx,
                         std::vector< std::pair<int,Box> >& isects,
			 bool first_only) const
{
    BoxHashMapType& BoxHashMap = m_ref->hash;

#ifdef _OPENMP
    #pragma omp critical(intersections_lock)
#endif
    {
        if (BoxHashMap.empty() && size() > 0)
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

            boundingbox.coarsen(maxext);

            m_ref->crsn = maxext;
            m_ref->bbox = boundingbox;

            for (int i = 0, N = size(); i < N; i++)
            {
                BoxHashMap[BoxLib::coarsen(get(i).smallEnd(),maxext)].push_back(i);
            }
        }
    }

    isects.resize(0);

    if (!BoxHashMap.empty())
    {
        BL_ASSERT(bx.sameType(get(0)));

        Box           cbx = BoxLib::coarsen(bx, m_ref->crsn);
        const IntVect  sm = BoxLib::max(cbx.smallEnd()-1, m_ref->bbox.smallEnd());
        const IntVect  bg = BoxLib::min(cbx.bigEnd(),     m_ref->bbox.bigEnd());

        cbx = Box(sm,bg,bx.ixType());

        ConstBoxHashMapIter TheEnd = BoxHashMap.end();

        for (IntVect iv = cbx.smallEnd(), End = cbx.bigEnd(); iv <= End; cbx.next(iv))
        {
            ConstBoxHashMapIter it = BoxHashMap.find(iv);

            if (it != TheEnd)
            {
                for (std::vector<int>::const_iterator v_it = it->second.begin(), v_end = it->second.end();
                     v_it != v_end;
                     ++v_it)
                {
                    const int index = *v_it;
                    const Box isect = bx & get(index);

                    if (isect.ok())
                    {
                        isects.push_back(std::pair<int,Box>(index,isect));
			if (first_only) return;
                    }
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

    std::map< IntVect,std::vector<int>,IntVect::Compare >& BoxHashMap = m_ref->hash;

    typedef std::map< IntVect,std::vector<int>,IntVect::Compare >::const_iterator ConstBoxHashMapIter;

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

                    BoxHashMap[BoxLib::coarsen(it->smallEnd(),m_ref->crsn)].push_back(size()-1);
                }
            }
        }
    }
    //
    // We now have "holes" in our BoxArray. Make us good.
    //
    bl.clear();

    const Box& bb = m_ref->bbox;

    ConstBoxHashMapIter TheEnd = BoxHashMap.end();

    for (IntVect iv = bb.smallEnd(), End = bb.bigEnd(); iv <= End; bb.next(iv))
    {
        ConstBoxHashMapIter it = BoxHashMap.find(iv);

        if (it != TheEnd)
        {
            for (std::vector<int>::const_iterator v_it = it->second.begin(), v_end = it->second.end();
                 v_it != v_end;
                 ++v_it)
            {
                const int index = *v_it;

                if (m_ref->m_abox[index].ok())
                {
                    bl.push_back(m_ref->m_abox[index]);
                }
            }
        }
    }

    bl.simplify();

    BoxArray nba(bl);

    bl.clear();

    *this = nba;

    BL_ASSERT(isDisjoint());
}
