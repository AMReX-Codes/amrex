//
// $Id: BoxArray.cpp,v 1.18 2001-07-17 23:02:19 lijewski Exp $
//

#include <BLassert.H>
#include <BoxArray.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

int
BoxArray::length () const
{
    return m_ref->m_abox.length();
}

bool
BoxArray::ready () const
{
    return m_ref->m_abox.ready();
}

const Box&
BoxArray::operator[] (int index) const
{
    return m_ref->m_abox.get(index);
}

const Box&
BoxArray::get (int index) const
{
    return m_ref->m_abox.get(index);
}

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

BoxArray::BoxArray (std::istream& is)
    :
    m_ref(new BoxArray::Ref(is))
{}

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

BoxArray::Ref::Ref (const Box* bxvec,
                    int        nbox)
    :
    m_abox(bxvec,nbox)
{}

BoxArray::BoxArray (const Box* bxvec,
                    int        nbox)
    :
    m_ref(new BoxArray::Ref(bxvec, nbox))
{}

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
BoxArray::define (std::istream& is)
{
    BL_ASSERT(length() == 0);
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
    BL_ASSERT(m_abox.length() == 0);
    int           maxbox;
    unsigned long hash;
    is.ignore(BL_IGNORE_MAX, '(') >> maxbox >> hash;
    m_abox.resize(maxbox);
    for (int i = 0; i < m_abox.length(); i++)
        is >> m_abox.get(i);
    is.ignore(BL_IGNORE_MAX, ')');
    if (is.fail())
        BoxLib::Error("BoxArray::define(istream&) failed");
}

void
BoxArray::define (const BoxList& bl)
{
    BL_ASSERT(length() == 0);
    if (!m_ref.unique())
        uniqify();
    m_ref->define(bl);
}

void
BoxArray::Ref::define (const BoxList& bl)
{
    BL_ASSERT(m_abox.length() == 0);
    m_abox.resize(bl.length());
    int count = 0;
    for (BoxListIterator bli(bl); bli; ++bli)
        m_abox.get(count++) = bli();
}

void
BoxArray::define (const BoxArray& bs)
{
    BL_ASSERT(length() == 0);
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
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).refine(refinement_ratio);
    return *this;
}

BoxArray&
BoxArray::refine (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i <length(); i++)
        m_ref->m_abox.get(i).refine(iv);
    return *this;
}

BoxArray&
BoxArray::shift (int dir,
                 int nzones)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).shift(dir, nzones);
    return *this;
}

BoxArray&
BoxArray::shiftHalf (int dir,
                     int num_halfs)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).shiftHalf(dir, num_halfs);
    return *this;
}

BoxArray&
BoxArray::shiftHalf (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).shiftHalf(iv);
    return *this;
}

BoxArray&
BoxArray::coarsen (int refinement_ratio)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).coarsen(refinement_ratio);
    return *this;
}

BoxArray&
BoxArray::coarsen (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).coarsen(iv);
    return *this;
}

BoxArray&
BoxArray::grow (int n)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).grow(n);
    return *this;
}

BoxArray&
BoxArray::grow (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).grow(iv);
    return *this;
}

BoxArray&
BoxArray::grow (int dir,
                int n_cell)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).grow(dir, n_cell);
    return *this;
}

bool
BoxArray::contains (const IntVect& v) const
{
    for (int i = 0; i < length(); i++)
        if (m_ref->m_abox.get(i).contains(v))
            return true;
    return false;
}

bool
BoxArray::contains (const Box& b) const
{
#ifndef BL_NAMESPACE
    BoxArray bnew = ::complementIn(b, *this);
#else
    BoxArray bnew = BL_NAMESPACE::complementIn(b, *this);
#endif
    return bnew.length() == 0;
}

bool
BoxArray::contains (const BoxArray& bl) const
{
    for (int i = 0; i < bl.length(); i++)
       if (!contains(bl.m_ref->m_abox.get(i)))
           return false;
    return true;
}

BoxArray&
BoxArray::surroundingNodes ()
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).surroundingNodes();
    return *this;
}

BoxArray&
BoxArray::surroundingNodes (int dir)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).surroundingNodes(dir);
    return *this;
}

BoxArray&
BoxArray::enclosedCells ()
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).enclosedCells();
    return *this;
}

BoxArray&
BoxArray::enclosedCells (int dir)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).enclosedCells(dir);
    return *this;
}

BoxArray&
BoxArray::convert (IndexType typ)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).convert(typ);
    return *this;
}

BoxArray&
BoxArray::convert (Box (*fp)(const Box&))
{
    BL_ASSERT(!(fp == 0));

    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox[i] = (*fp)(m_ref->m_abox[i]);
    return *this;
}

std::ostream&
BoxArray::writeOn (std::ostream& os) const
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    os << '(' << length() << ' ' << 0 << '\n';

    for (int i = 0; i < length(); i++)
        os << get(i) << '\n';

    os << ')';

    if (os.fail())
        BoxLib::Error("BoxArray::writeOn(ostream&) failed");

    return os;
}

bool
BoxArray::isDisjoint () const
{
    for (int i = 0; i < length(); i++)
        for (int j = i + 1; j < length(); j++)
            if (get(i).intersects(get(j)))
                return false;
    return true;
}

bool
BoxArray::ok () const
{
    bool isok = true;

    if (length() > 0)
    {
        const Box& bx0 = m_ref->m_abox[0];

        if (length() == 1)
            isok = bx0.ok();

        for (int i = 1; i < length() && isok; i++)
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
  for ( int i = 0; i < length(); ++i )
    {
      result += m_ref->m_abox.get(i).numPts();
    }
  return result;
}

BoxArray&
BoxArray::maxSize (int block_size)
{
    BoxList blst(*this);
    blst.maxSize(block_size);
    clear();
    m_ref->m_abox.resize(blst.length());
    BoxListIterator bli(blst);
    for (int i = 0; bli; ++bli)
        set(i++, bli());
    return *this;
}

std::ostream&
operator<< (std::ostream&   os,
            const BoxArray& ba)
{
    return ba.print(os);
}

std::ostream&
BoxArray::print (std::ostream& os) const
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    os << "(BoxArray maxbox("
       << m_ref->m_abox.length()
       << ")\n       m_ref->m_hash_sig("
       << 0
       << ")\n       ";

    for (int i = 0; i < m_ref->m_abox.length(); ++i)
        os << m_ref->m_abox[i] << ' ';

    os << ")\n";

    if (os.fail())
        BoxLib::Error("BoxArray::print(ostream& os) failed");

    return os;
}

BoxList
BoxArray::boxList () const
{
    BL_ASSERT(length() > 0);
    BoxList newb(get(0).ixType());
    for (int i = 0; i < length(); ++i)
        newb.append(get(i));
    return newb;
}

Box
BoxArray::minimalBox () const
{
    Box minbox;
    if (length() > 0)
    {
        minbox = m_ref->m_abox.get(0);
        for (int i = 0; i < length(); i++)
            minbox.minBox(m_ref->m_abox.get(i));
    }
    return minbox;
}

BoxArray
boxComplement (const Box& b1in,
               const Box& b2)
{
#ifndef BL_NAMESPACE
    return BoxArray(::boxDiff(b1in, b2));
#else
    return BoxArray(BL_NAMESPACE::boxDiff(b1in, b2));
#endif
}

BoxArray
complementIn (const Box&      b,
              const BoxArray& ba)
{
#ifndef BL_NAMESPACE
    return BoxArray(::complementIn(b, ba.boxList()));
#else
    return BoxArray(BL_NAMESPACE::complementIn(b, ba.boxList()));
#endif
}

BoxArray
intersect (const BoxArray& ba,
           const Box&      b)
{
#ifndef BL_NAMESPACE
    return BoxArray(::intersect(ba.boxList(), b));
#else
    return BoxArray(BL_NAMESPACE::intersect(ba.boxList(), b));
#endif
}

BoxArray
intersect (const BoxArray& lhs,
           const BoxArray& rhs)
{
#ifndef BL_NAMESPACE
    return BoxArray(::intersect(lhs.boxList(), rhs.boxList()));
#else
    return BoxArray(BL_NAMESPACE::intersect(lhs.boxList(), rhs.boxList()));
#endif
}


#ifdef BL_NAMESPACE
}
#endif

