//
// $Id: BoxDomain.cpp,v 1.12 2001-07-23 22:30:23 car Exp $
//

#include <BoxDomain.H>

bool
BoxDomain::contains (const IntVect& v) const
{
    return BoxList::contains(v);
}

bool
BoxDomain::contains (const Box& b) const
{
    return BoxList::contains(b);
}

bool
BoxDomain::contains (const BoxList& bl) const
{
    return BoxList::contains(bl);
}

BoxDomain&
BoxDomain::intersect (const Box& b)
{
    BoxList::intersect(b);
    BL_ASSERT(ok());
    return *this;
}

void
BoxLib::intersect (BoxDomain&       dest,
		   const BoxDomain& fin,
		   const Box&       b)
{
   dest = fin;
   dest.intersect(b);
}

BoxDomain&
BoxDomain::refine (int ratio)
{
    BoxList::refine(ratio);
    BL_ASSERT(ok());
    return *this;
}

void
BoxLib::refine (BoxDomain&       dest,
		const BoxDomain& fin,
		int              ratio)
{
    dest = fin;
    dest.refine(ratio);
}

void
BoxLib::accrete (BoxDomain&       dest,
		 const BoxDomain& fin,
		 int              sz)
{
    dest = fin;
    dest.accrete(sz);
}

void
BoxLib::coarsen (BoxDomain&       dest,
		 const BoxDomain& fin,
		 int              ratio)
{
    dest = fin;
    dest.coarsen(ratio);
}

BoxDomain&
BoxDomain::complementIn (const Box&       b,
                         const BoxDomain& bl)
{
    BoxList::complementIn(b,bl);
    BL_ASSERT(ok());
    return *this;
}

BoxDomain
BoxLib::complementIn (const Box&       b,
		      const BoxDomain& bl)
{
    BoxDomain result;
    result.complementIn(b,bl);
    return result;
}

BoxDomain&
BoxDomain::shift (int dir,
                  int nzones)
{
    BoxList::shift(dir, nzones);
    return *this;
}

BoxDomain&
BoxDomain::shiftHalf (int dir,
                      int num_halfs)
{
    BoxList::shiftHalf(dir, num_halfs);
    return *this;
}

BoxDomain&
BoxDomain::shiftHalf (const IntVect& iv)
{
    BoxList::shiftHalf(iv);
    return *this;
}

BoxList
BoxDomain::boxList () const
{
    return BoxList(*this);
}

bool
BoxDomain::operator== (const BoxDomain& rhs) const
{
    return BoxList::operator==(rhs);
}

bool
BoxDomain::operator!= (const BoxDomain& rhs) const
{
    return !BoxList::operator==(rhs);
}

BoxDomain::BoxDomain ()
    : BoxList(IndexType::TheCellType())
{}

BoxDomain::BoxDomain (IndexType _ctype)
    : BoxList(_ctype)
{}

void
BoxDomain::add (const Box& b)
{
    BL_ASSERT(b.ixType() == ixType());

    std::list<Box> check;
    check.push_back(b);
    for (Iterator bli = lbox.begin(); bli != lbox.end(); ++bli)
    {
        std::list<Box> tmp;
        for (Iterator ci = check.begin(); ci != check.end(); )
        {
            if (ci->intersects(*bli))
            {
                //
                // Remove c from the check list, compute the
                // part of it that is outside bln and collect
                // those boxes in the tmp list.
                //
                BoxList tmpbl(BoxLib::boxDiff(*ci, *bli));
                tmp.splice(tmp.end(), tmpbl.listBox());
                ci = check.erase(ci);
            }
            else
                ++ci;
        }
        check.splice(check.end(), tmp);
    }
    //
    // At this point, the only thing left in the check list
    // are boxes that nowhere intersect boxes in the domain.
    //
    lbox.splice(lbox.end(), check);
    BL_ASSERT(ok());
}

void
BoxDomain::add (const BoxList& bl)
{
    for (BoxList::ConstIterator bli = bl.begin(); bli != bl.end(); ++bli)
        add(*bli);
}

BoxDomain&
BoxDomain::rmBox (const Box& b)
{
    BL_ASSERT(b.ixType() == ixType());

    std::list<Box> tmp;

    for (std::list<Box>::iterator  bli = lbox.begin(); bli != lbox.end(); )
    {
        if (bli->intersects(b))
        {
            BoxList tmpbl(BoxLib::boxDiff(*bli,b));
            tmp.splice(tmp.end(), tmpbl.listBox());
            bli = lbox.erase(bli);
        }
        else
            ++bli;
    }
    lbox.splice(lbox.end(), tmp);
    return *this;
}

bool
BoxDomain::ok () const
{
    //
    // First check to see if boxes are valid.
    //
    bool status = BoxList::ok();
    if (status)
    {
        //
        // Now check to see that boxes are disjoint.
        //
        for (ConstIterator bli = begin(); bli != end(); ++bli)
        {
            ConstIterator blii = bli; ++blii;
            while ( blii != end() )
            {
                if ( bli->intersects(*blii) )
                {
                    std::cout << "Invalid DOMAIN, boxes overlap" << '\n'
                              << "b1 = " << *bli << '\n'
                              << "b2 = " << *blii << '\n';
                    status = false;
                }
                ++blii;
            }
        }
    }
    return status;
}

BoxDomain&
BoxDomain::accrete (int sz)
{
    BoxList bl(*this);
    bl.accrete(sz);
    clear();
    add(bl);
    return *this;
}

BoxDomain&
BoxDomain::coarsen (int ratio)
{
    BoxList bl(*this);
    bl.coarsen(ratio);
    clear();
    add(bl);
    return *this;
}

std::ostream&
operator<< (std::ostream&    os,
            const BoxDomain& bd)
{
    os << "(BoxDomain " << bd.boxList() << ")" << std::flush;
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,BoxDomain&) failed");
    return os;
}

