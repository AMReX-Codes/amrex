//BL_COPYRIGHT_NOTICE

//
// $Id: BoxArray.cpp,v 1.4 1997-11-18 00:04:47 lijewski Exp $
//

#include <Assert.H>
#include <BoxArray.H>
#include <Utility.H>

BoxArray::BoxArray ()
    :
    abox(),
    hash_sig(0)
{}

BoxArray::BoxArray (const BoxList& bl)
    :
    abox(),
    hash_sig(0)
{
    define(bl);
}

BoxArray::BoxArray (istream& is)
    :
    abox(),
    hash_sig(0)
{
    define(is);
}

BoxArray&
BoxArray::operator= (const BoxArray& rhs)
{
    abox     = rhs.abox;
    hash_sig = rhs.hash_sig;
    return *this;
}

BoxArray::BoxArray (const BoxArray& bs)
    :
    abox(bs.abox)
{
    hash_sig = bs.hash_sig;
}

BoxArray::BoxArray (size_t size)
    :
    abox(size)
{
    hash_sig = do_hash();
}

BoxArray::BoxArray (const Box* bxvec,
                    int        nbox)
    :
    abox(bxvec, nbox)
{
    hash_sig = do_hash();
}

void
BoxArray::set (int        i,
               const Box& ibox)
{
    abox.set(i, ibox);
    hash_sig = do_hash();
}

void
BoxArray::define (istream& is)
{
    assert(length() == 0);
    int           maxbox;
    unsigned long in_hash;
    is.ignore(BL_IGNORE_MAX, '(') >> maxbox >> in_hash;
    abox.resize(maxbox);
    for (int i = 0; i < length(); i++)
        is >> abox.get(i);
    is.ignore(BL_IGNORE_MAX, ')');
    hash_sig = do_hash();
    assert(hash_sig == in_hash);

    if (is.fail())
        BoxLib::Error("BoxArray::define(istream&) failed");
}

void
BoxArray::define (const BoxList& bl)
{
    assert(length() == 0);
    //
    // Init box's and compute hash_sig at the same time.
    //
    abox.resize(bl.length());
    int count = 0;
    for (BoxListIterator bli(bl); bli; ++bli)
        abox.get(count++) = bli();
    hash_sig = do_hash();
}

unsigned long
BoxArray::do_hash() const
{
    unsigned long hash = 0;

    for (int i = 0; i < length(); i++)
    {
        const int* hiv = get(i).hiVect();
        const int* lov = get(i).loVect();
        for (int k=0; k<SpaceDim; k++)
        {
            hash = (hash << 3) + hiv[k];
            //
            // This uses the fact that ANSI C guarantees that unsigned
            // longs have at least 32 bits.
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
BoxArray::define (const BoxArray& bs)
{
    //assert(length() == 0);
    //
    // Init box's and compute hash_sig at the same time.
    //
    abox.resize(bs.length());
    for (int i=0; i<length(); i++)
        abox.set(i,bs[i]);
    hash_sig = bs.hash_sig;
}

BoxArray::~BoxArray ()
{
    clear();
}

BoxArray&
BoxArray::refine (int refinement_ratio)
{
    for (int i=0; i<length(); i++)
        abox.get(i).refine(refinement_ratio);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::refine (const IntVect& iv)
{
    for (int i=0; i<length(); i++)
        abox.get(i).refine(iv);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::shift (int dir,
                 int nzones)
{
    for (int i = 0; i < length(); i++)
        abox.get(i).shift(dir, nzones);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::shiftHalf (int dir,
                     int num_halfs)
{
    for (int i = 0; i < length(); i++)
        abox.get(i).shiftHalf(dir, num_halfs);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::shiftHalf (const IntVect& iv)
{
    for (int i = 0; i < length(); i++)
        abox.get(i).shiftHalf(iv);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::coarsen (int refinement_ratio)
{
    for (int i=0; i<length(); i++)
        abox.get(i).coarsen(refinement_ratio);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::coarsen (const IntVect& iv)
{
    for (int i=0; i<length(); i++)
        abox.get(i).coarsen(iv);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::grow (int n)
{
    for (int i=0; i<length(); i++)
        abox.get(i).grow(n);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::grow (const IntVect& iv)
{
    for (int i=0; i<length(); i++)
        abox.get(i).grow(iv);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::grow (int dir,
                int n_cell)
{
    for (int i=0; i<length(); i++)
        abox.get(i).grow(dir, n_cell);
    hash_sig = do_hash();
    return *this;
}

bool
BoxArray::contains (const IntVect& v) const
{
    bool contained = false;
    for (int i = 0; i < length() && !contained; i++)
        if (abox.get(i).contains(v))
            contained = true;
    return contained;
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
    bool contained = true;
    for (int i = 0; i < length() && !contained; i++)
       if (!contains(bl.abox.get(i)))
           contained = false;
    return contained;
}

BoxArray&
BoxArray::surroundingNodes ()
{
    for (int i=0; i<length(); i++)
        abox.get(i).surroundingNodes();
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::surroundingNodes (int dir)
{
    for (int i=0; i<length(); i++)
        abox.get(i).surroundingNodes(dir);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::enclosedCells ()
{
    for (int i=0; i<length(); i++)
        abox.get(i).enclosedCells();
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::enclosedCells (int dir)
{
    for (int i=0; i<length(); i++)
        abox.get(i).enclosedCells(dir);
    hash_sig = do_hash();
    return *this;
}

BoxArray&
BoxArray::convert (IndexType typ)
{
    for (int i=0; i<length(); i++)
        abox.get(i).convert(typ);
    hash_sig = do_hash();
    return *this;
}

ostream&
BoxArray::writeOn (ostream& os) const
{
    os << '(' << length() << ' ' << hash_sig << '\n';
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
    bool rc = true;
    for (int i = 0; i < length() && rc; i++)
    {
        for (int j = i + 1; j < length() && rc; j++)
            if (get(i).intersects(get(j)))
                rc = 0;
    }
    return rc;
}

bool
BoxArray::ok () const
{
    bool isok = true;
    if (length() > 0)
    {
        const Box& bx0 = abox[0];
        if (length() == 1)
            isok = bx0.ok();
        for (int i = 1; i < length() && isok; i++)
        {
            const Box& bxi = abox[i];
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
    abox.resize(blst.length());
    BoxListIterator bli(blst);
    for (int i = 0; bli; ++bli)
        set(i++, bli());
    return *this;
}

ostream&
operator<< (ostream&        os,
            const BoxArray& bs)
{
    os << "(BoxArray maxbox(" << bs.length() << ")\n";
    os << "       hash_sig("  << bs.hash_sig << ")\n";
    os << "       ";
    for (int i = 0; i < bs.length(); ++i)
        os << bs[i] << ' ';
    os << ')' << '\n';

    if (os.fail())
        BoxLib::Error("operator<<(ostream&,BoxArray&) failed");

    return os;
}

BoxList
BoxArray::boxList () const
{
    assert(length() > 0);
    BoxList newb(get(0).ixType());
    for (int i =0; i < length(); ++i)
        newb.append(get(i));
    return newb;
}

Box
BoxArray::minimalBox () const
{
    Box minbox;
    if (length() > 0)
    {
        int i = 0;
        minbox = abox.get(i);
        for ( ; i < length(); i++)
            minbox.minBox(abox.get(i));
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

