//BL_COPYRIGHT_NOTICE

//
// $Id: BoxDomain.cpp,v 1.5 2000-04-24 17:52:33 car Exp $
//

#include <BoxDomain.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

BoxDomain::BoxDomain ()
    : BoxList(IndexType::TheCellType())
{}

BoxDomain::BoxDomain (IndexType _ctype)
    : BoxList(_ctype)
{}

BoxDomain::BoxDomain (const BoxDomain& rhs)
    : BoxList(rhs)
{}

BoxDomain&
BoxDomain::operator= (const BoxDomain& rhs)
{
    BoxList::operator=(rhs);
    return *this;
}

BoxDomain::~BoxDomain ()
{}

void
BoxDomain::add (const Box& b)
{
    BL_ASSERT(b.ixType() == ixType());

    List<Box> check;
    check.append(b);
    for (ListIterator<Box> bli(lbox); bli; ++bli)
    {
        List<Box> tmp;
        for (ListIterator<Box> ci(check); ci; )
        {
            if (ci().intersects(bli()))
            {
                //
                // Remove c from the check list, compute the
                // part of it that is outside bln and collect
                // those boxes in the tmp list.
                //
                BoxList tmpbl(boxDiff(ci(), bli()));
                tmp.catenate(tmpbl.listBox());
                check.remove(ci);
            }
            else
                ++ci;
        }
        check.catenate(tmp);
    }
    //
    // At this point, the only thing left in the check list
    // are boxes that nowhere intersect boxes in the domain.
    //
    lbox.catenate(check);
    BL_ASSERT(ok());
}

void
BoxDomain::add (const BoxList& bl)
{
    for (BoxListIterator bli(bl); bli; ++bli)
        add(*bli);
}

BoxDomain&
BoxDomain::rmBox (const Box& b)
{
    BL_ASSERT(b.ixType() == ixType());

    List<Box> tmp;

    for (ListIterator<Box> bli(lbox); bli; )
    {
        if (bli().intersects(b))
        {
            BoxList tmpbl(boxDiff(bli(),b));
            tmp.catenate(tmpbl.listBox());
            lbox.remove(bli);
        }
        else
            ++bli;
    }
    lbox.catenate(tmp);
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
        for (BoxListIterator bli(*this); bli; ++bli)
        {
            BoxListIterator blii(bli); ++blii;
            while (blii)
            {
                if (bli().intersects(blii()))
                {
                    cout << "Invalid DOMAIN, boxes overlap" << '\n';
                    cout << "b1 = " << bli() << '\n';
                    cout << "b2 = " << blii() << '\n';
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

ostream&
operator<< (ostream&         os,
            const BoxDomain& bd)
{
    os << "(BoxDomain " << BoxList(bd) << ")" << flush;
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,BoxDomain&) failed");
    return os;
}

#ifdef BL_NAMESPACE
}
#endif

