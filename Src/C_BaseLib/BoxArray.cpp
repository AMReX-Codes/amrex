//
// $Id: BoxArray.cpp,v 1.34 2003-06-24 17:18:21 lijewski Exp $
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
        if (m_ref->m_abox.get(i).contains(v))
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
        for (int j = i + 1; j < size(); j++)
            if (get(i).intersects(get(j)))
                return false;
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
    return BoxArray(BoxLib::complementIn(b, ba.boxList()));
}

BoxArray
BoxLib::intersect (const BoxArray& ba,
		   const Box&      b)
{
    return BoxArray(BoxLib::intersect(ba.boxList(), b));
}

BoxArray
BoxLib::intersect (const BoxArray& lhs,
		   const BoxArray& rhs)
{
    return BoxArray(BoxLib::intersect(lhs.boxList(), rhs.boxList()));
}
