//BL_COPYRIGHT_NOTICE

//
// $Id: BoxArray.cpp,v 1.8 1998-03-26 17:32:17 lijewski Exp $
//

#include <Assert.H>
#include <BoxArray.H>

//
// Returns a 24 bit (nearly) unique number.
//
static
unsigned long
Hash (const Array<Box>& arr)
{
    unsigned long hash = 0;

    for (int i = 0; i < arr.length(); i++)
    {
        const int* hiv = arr.get(i).hiVect();
        const int* lov = arr.get(i).loVect();

        for (int k = 0; k < SpaceDim; k++)
        {
            hash = (hash << 3) + hiv[k];
            //
            // ANSI C guarantees that unsigned longs have at least 32 bits.
            //
            unsigned long g;
            if ((g = (hash & 0xf0000000)))
            {
                hash ^= g >> 24;
                hash ^= g;
            }
            hash = (hash << 3) + lov[k];
            if ((g = (hash & 0xf0000000)))
            {
                hash ^= g >> 24;
                hash ^= g;
            }
        }
    }

    return hash;
}

void
BoxArray::rehash ()
{
    m_ref->m_hash_sig = Hash(m_ref->m_abox);
}

BoxArray::Ref::Ref ()
    :
    m_hash_sig(0)
{}

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

BoxArray::Ref::Ref (istream& is)
{
    define(is);
}

BoxArray::BoxArray (istream& is)
    :
    m_ref(new BoxArray::Ref(is))
{}

BoxArray::BoxArray (const BoxArray& rhs)
    :
    m_ref(rhs.m_ref)
{}

BoxArray::Ref::Ref (size_t size)
    :
    m_abox(size),
    m_hash_sig(Hash(m_abox))
{}

BoxArray::BoxArray (size_t size)
    :
    m_ref(new BoxArray::Ref(size))
{}

BoxArray::Ref::Ref (const Box* bxvec,
                    int        nbox)
    :
    m_abox(bxvec, nbox),
    m_hash_sig(Hash(m_abox))
{}

BoxArray::BoxArray (const Box* bxvec,
                    int        nbox)
    :
    m_ref(new BoxArray::Ref(bxvec, nbox))
{}

BoxArray::Ref::Ref (const Ref& rhs)
    :
    m_abox(rhs.m_abox),
    m_hash_sig(rhs.m_hash_sig)
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
    m_ref->m_hash_sig = 0;
}

void
BoxArray::resize (int len)
{
    if (!m_ref.unique())
        uniqify();
    m_ref->m_abox.resize(len);
    rehash();
}

void
BoxArray::set (int        i,
               const Box& ibox)
{
    if (!m_ref.unique())
        uniqify();
    m_ref->m_abox.set(i, ibox);
    rehash();
}

bool
BoxArray::operator== (const BoxArray& rhs) const
{
    return m_ref == rhs.m_ref || *m_ref == *rhs.m_ref;
}

bool
BoxArray::Ref::operator== (const BoxArray::Ref& rhs) const
{
    return m_hash_sig == rhs.m_hash_sig && m_abox == rhs.m_abox;
}

//
// Moved out of Utility.H
//
#define BL_IGNORE_MAX 100000

void
BoxArray::define (istream& is)
{
    assert(length() == 0);
    if (!m_ref.unique())
        uniqify();
    m_ref->define(is);
}

void
BoxArray::Ref::define (istream& is)
{
    assert(m_abox.length() == 0);

    int           maxbox;
    unsigned long in_hash;
    is.ignore(BL_IGNORE_MAX, '(') >> maxbox >> in_hash;
    m_abox.resize(maxbox);
    for (int i = 0; i < m_abox.length(); i++)
        is >> m_abox.get(i);
    is.ignore(BL_IGNORE_MAX, ')');
    m_hash_sig = Hash(m_abox);
    assert(m_hash_sig == in_hash);

    if (is.fail())
        BoxLib::Error("BoxArray::define(istream&) failed");
}

void
BoxArray::define (const BoxList& bl)
{
    assert(length() == 0);
    if (!m_ref.unique())
        uniqify();
    m_ref->define(bl);
}

void
BoxArray::Ref::define (const BoxList& bl)
{
    assert(m_abox.length() == 0);
    //
    // Init box's and compute m_ref->m_hash_sig at the same time.
    //
    m_abox.resize(bl.length());
    int count = 0;
    for (BoxListIterator bli(bl); bli; ++bli)
        m_abox.get(count++) = bli();
    m_hash_sig = Hash(m_abox);
}

void
BoxArray::define (const BoxArray& bs)
{
    assert(length() == 0);
    if (!m_ref.unique())
        uniqify();
    m_ref->define(bs.m_ref->m_abox, bs.m_ref->m_hash_sig);
}

void
BoxArray::Ref::define (const Array<Box>& ba,
                       unsigned long     hash)
{
    assert(m_abox.length() == 0);
    //
    // Init box's and compute m_hash_sig at the same time.
    //
    m_abox.resize(ba.length());
    for (int i = 0; i < ba.length(); i++)
        m_abox.set(i,ba[i]);
    m_hash_sig = hash;
}

BoxArray::~BoxArray () {}

BoxArray&
BoxArray::refine (int refinement_ratio)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).refine(refinement_ratio);
    rehash();
    return *this;
}

BoxArray&
BoxArray::refine (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i <length(); i++)
        m_ref->m_abox.get(i).refine(iv);
    rehash();
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
    rehash();
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
    rehash();
    return *this;
}

BoxArray&
BoxArray::shiftHalf (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).shiftHalf(iv);
    rehash();
    return *this;
}

BoxArray&
BoxArray::coarsen (int refinement_ratio)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).coarsen(refinement_ratio);
    rehash();
    return *this;
}

BoxArray&
BoxArray::coarsen (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).coarsen(iv);
    rehash();
    return *this;
}

BoxArray&
BoxArray::grow (int n)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).grow(n);
    rehash();
    return *this;
}

BoxArray&
BoxArray::grow (const IntVect& iv)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).grow(iv);
    rehash();
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
    rehash();
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
    BoxArray bnew = ::complementIn(b, *this);
    return bnew.length() == 0;
}

bool
BoxArray::contains (const BoxArray& bl) const
{
    for (int i = 0; i < length(); i++)
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
    rehash();
    return *this;
}

BoxArray&
BoxArray::surroundingNodes (int dir)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).surroundingNodes(dir);
    rehash();
    return *this;
}

BoxArray&
BoxArray::enclosedCells ()
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).enclosedCells();
    rehash();
    return *this;
}

BoxArray&
BoxArray::enclosedCells (int dir)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).enclosedCells(dir);
    rehash();
    return *this;
}

BoxArray&
BoxArray::convert (IndexType typ)
{
    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox.get(i).convert(typ);
    rehash();
    return *this;
}

BoxArray&
BoxArray::convert (Box (*fp)(const Box&))
{
    assert(!(fp == 0));

    if (!m_ref.unique())
        uniqify();
    for (int i = 0; i < length(); i++)
        m_ref->m_abox[i] = (*fp)(m_ref->m_abox[i]);
    rehash();
    return *this;
}

ostream&
BoxArray::writeOn (ostream& os) const
{
    os << '(' << length() << ' ' << m_ref->m_hash_sig << '\n';
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
    {
        for (int j = i + 1; j < length(); j++)
        {
            if (get(i).intersects(get(j)))
                return false;
        }
    }
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

ostream&
operator<< (ostream&        os,
            const BoxArray& ba)
{
    return ba.print(os);
}

ostream&
BoxArray::print (ostream& os) const
{
    os << "(BoxArray maxbox("
       << m_ref->m_abox.length()
       << ")\n       m_ref->m_hash_sig("
       << m_ref->m_hash_sig
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
    assert(length() > 0);
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
    return BoxArray(boxDiff(b1in, b2));
}

BoxArray
complementIn (const Box&      b,
              const BoxArray& ba)
{
    return BoxArray(complementIn(b, ba.boxList()));
}

BoxArray
intersect (const BoxArray& ba,
           const Box&      b)
{
    return BoxArray(intersect(ba.boxList(), b));
}

