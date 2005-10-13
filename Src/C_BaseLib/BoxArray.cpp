//
// $Id: BoxArray.cpp,v 1.43 2005-10-13 23:04:50 lijewski Exp $
//
#include <iostream>

#include <BLassert.H>
#include <BoxArray.H>

void
BoxArray::reserve (long _truesize)
{
    if (!m_ref.unique())
        uniqify();
    m_ref->m_abox.reserve(_truesize);
}

BoxArray::Ref::Ref () {}

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

BoxArray::BoxArray (const BoxArray& rhs)
    :
    m_ref(rhs.m_ref)
{}

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

BoxArray&
BoxArray::operator= (const BoxArray& rhs)
{
    m_ref = rhs.m_ref;
    return *this;
}

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
    for (int i = 0; i < m_abox.size(); i++)
        is >> m_abox.get(i);
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
    for (BoxList::const_iterator bli = bl.begin(); bli != bl.end(); ++bli)
        m_abox.get(count++) = *bli;
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
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).refine(refinement_ratio);
    return *this;
}

BoxArray&
BoxArray::refine (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i <size(); i++)
        m_ref->m_abox.get(i).refine(iv);
    return *this;
}

BoxArray&
BoxArray::shift (int dir,
                 int nzones)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).shift(dir, nzones);
    return *this;
}

BoxArray&
BoxArray::shiftHalf (int dir,
                     int num_halfs)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).shiftHalf(dir, num_halfs);
    return *this;
}

BoxArray&
BoxArray::shiftHalf (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).shiftHalf(iv);
    return *this;
}

BoxArray&
BoxArray::coarsen (int refinement_ratio)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).coarsen(refinement_ratio);
    return *this;
}

BoxArray&
BoxArray::coarsen (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).coarsen(iv);
    return *this;
}

BoxArray&
BoxArray::grow (int n)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).grow(n);
    return *this;
}

BoxArray&
BoxArray::grow (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).grow(iv);
    return *this;
}

BoxArray&
BoxArray::grow (int dir,
                int n_cell)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).grow(dir, n_cell);
    return *this;
}

bool
BoxArray::contains (const IntVect& v) const
{
    for (int i = 0; i < size(); i++)
        if (get(i).contains(v))
            return true;
    return false;
}

bool
BoxArray::contains (const Box& b) const
{
    BoxArray bnew = BoxLib::complementIn(b, *this);
    return bnew.size() == 0;
}

bool
BoxArray::contains (const BoxArray& bl) const
{
    for (int i = 0; i < bl.size(); i++)
       if (!contains(bl.m_ref->m_abox.get(i)))
           return false;
    return true;
}

BoxArray&
BoxArray::surroundingNodes ()
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).surroundingNodes();
    return *this;
}

BoxArray&
BoxArray::surroundingNodes (int dir)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).surroundingNodes(dir);
    return *this;
}

BoxArray&
BoxArray::enclosedCells ()
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).enclosedCells();
    return *this;
}

BoxArray&
BoxArray::enclosedCells (int dir)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).enclosedCells(dir);
    return *this;
}

BoxArray&
BoxArray::convert (IndexType typ)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox.get(i).convert(typ);
    return *this;
}

BoxArray&
BoxArray::convert (Box (*fp)(const Box&))
{
    BL_ASSERT(!(fp == 0));

    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < size(); i++)
        m_ref->m_abox[i] = (*fp)(m_ref->m_abox[i]);
    return *this;
}

std::ostream&
BoxArray::writeOn (std::ostream& os) const
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    os << '(' << size() << ' ' << 0 << '\n';

    for (int i = 0; i < size(); i++)
        os << get(i) << '\n';

    os << ')';

    if (os.fail())
        BoxLib::Error("BoxArray::writeOn(ostream&) failed");

    return os;
}

bool
BoxArray::isDisjoint () const
{
    for (int i = 0; i < size(); i++)
    {
        std::vector< std::pair<int,Box> > isects = intersections(get(i));
        if ( !(isects.size() == 1 && isects[0].second == get(i)) )
            return false;
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

        if (size() == 1)
            isok = bx0.ok();

        for (int i = 1; i < size() && isok; i++)
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
  for ( int i = 0; i < size(); ++i )
    {
      result += m_ref->m_abox.get(i).numPts();
    }
  return result;
}

double
BoxArray::d_numPts () const
{
  double result = 0;
  for ( int i = 0; i < size(); ++i )
  {
      result += m_ref->m_abox.get(i).d_numPts();
  }
  return result;
}

BoxArray&
BoxArray::maxSize (int block_size)
{
    BoxList blst(*this);
    blst.maxSize(block_size);
    clear();
    m_ref->m_abox.resize(blst.size());
    BoxList::iterator bli = blst.begin();
    for (int i = 0; bli != blst.end(); ++bli)
        set(i++, *bli);
    return *this;
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

    for (int i = 0; i < ba.size(); ++i)
        os << ba[i] << ' ';

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
    for (int i = 0; i < size(); ++i)
        newb.push_back(get(i));
    return newb;
}

Box
BoxArray::minimalBox () const
{
    Box minbox;
    if (size() > 0)
    {
        minbox = m_ref->m_abox.get(0);
        for (int i = 0; i < size(); i++)
            minbox.minBox(m_ref->m_abox.get(i));
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
    std::vector< std::pair<int,Box> > isects = ba.intersections(b);
    BoxList bl(b.ixType());
    for (int i = 0; i < isects.size(); ++i)
    {
        bl.push_back(isects[i].second);
    }
    return BoxArray(BoxLib::complementIn(b, bl));
}

BoxArray
BoxLib::intersect (const BoxArray& ba,
		   const Box&      b)
{
    std::vector< std::pair<int,Box> > isects = ba.intersections(b);
    BoxArray r(isects.size());
    for (int i = 0; i < isects.size(); i++)
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
    for (int i = 0; i < lhs.size(); i++)
    {
        BoxArray ba  = BoxLib::intersect(rhs, lhs[i]);
        BoxList  tmp = ba.boxList();
        bl.catenate(tmp);
    }
    return BoxArray(bl);
}

std::vector< std::pair<int,Box> >
BoxArray::intersections (const Box& bx) const
{
    if (!m_ref->hash.isAllocated() && size() > 0)
    {
        BL_ASSERT(bx.sameType(get(0)));
        //
        // Calculate the bounding box & maximum extent of the boxes.
        //
        IntVect maxext(D_DECL(0,0,0));

        Box boundingbox = get(0);

        for (int i = 0; i < size(); i++)
        {
            boundingbox.minBox(get(i));
            maxext = BoxLib::max(maxext, get(i).length());
        }
        m_ref->crsn = maxext;
        boundingbox.coarsen(maxext);

        m_ref->hash.resize(boundingbox, 1);

        for (int i = 0; i < size(); i++)
        {
            m_ref->hash(BoxLib::coarsen(get(i).smallEnd(),maxext)).push_back(i);
        }
    }

    std::vector< std::pair<int,Box> > isects;

    if (m_ref->hash.isAllocated())
    {
        BL_ASSERT(bx.sameType(get(0)));

        Box     cbx = BoxLib::coarsen(bx, m_ref->crsn);
        IntVect sm  = BoxLib::max(cbx.smallEnd()-1, m_ref->hash.box().smallEnd());
        IntVect bg  = BoxLib::min(cbx.bigEnd(),     m_ref->hash.box().bigEnd());

        cbx = Box(sm,bg,bx.ixType());

        for (IntVect iv = cbx.smallEnd(); iv <= cbx.bigEnd(); cbx.next(iv))
        {
            std::vector<int>& v = m_ref->hash(iv);

            for (int i = 0; i < v.size(); i++)
            {
                const Box isect = bx & get(v[i]);
                if (isect.ok())
                    isects.push_back(std::pair<int,Box>(v[i],isect));
            }
        }
    }

    return isects;
}

//
// Currently this assumes your Boxes are cell-centered.
//
BoxList
BoxArray::removeOverlap ()
{
    if (!m_ref.unique()) uniqify();

    const Box EmptyBox;

    for (int i = 0; i < size(); i++)
    {
        const Box b = m_ref->m_abox[i];

        if (b.ok())
        {
            std::vector< std::pair<int,Box> > isects = intersections(b);

            for (int j = 0; j < isects.size(); j++)
            {
                if (isects[j].first == i) continue;

                BoxList diff = BoxLib::boxDiff(m_ref->m_abox[isects[j].first], isects[j].second);

                m_ref->m_abox[isects[j].first] = EmptyBox;

                for (BoxList::const_iterator it = diff.begin(); it != diff.end(); ++it)
                {
                    m_ref->m_abox.push_back(*it);

                    m_ref->hash(BoxLib::coarsen(it->smallEnd(),m_ref->crsn)).push_back(size()-1);
                }
            }
        }
    }
    //
    // Now attempt to simplify the list of Boxes by just simplify()ing within hash bins.
    //
    BoxList nbl;
    Box     bb = m_ref->hash.box();

    for (IntVect iv = bb.smallEnd(); iv <= bb.bigEnd(); bb.next(iv))
    {
        std::vector<int>& v = m_ref->hash(iv);

        BoxList pieces;

        for (int i = 0; i < v.size(); i++)
            if (m_ref->m_abox[v[i]].ok())
                pieces.push_back(m_ref->m_abox[v[i]]);

        pieces.minimize();

        nbl.catenate(pieces);
    }

    BL_ASSERT(nbl.isDisjoint());

    return nbl;
}
