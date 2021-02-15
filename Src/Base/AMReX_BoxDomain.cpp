
#include <iostream>

#include <AMReX_BoxDomain.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>

namespace amrex {

BoxDomain&
BoxDomain::intersect (const Box& b)
{
    BoxList::intersect(b);
    BL_ASSERT(ok());
    return *this;
}

void
intersect (BoxDomain&       dest,
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
refine (BoxDomain&       dest,
        const BoxDomain& fin,
        int              ratio)
{
    dest = fin;
    dest.refine(ratio);
}

void
accrete (BoxDomain&       dest,
         const BoxDomain& fin,
         int              sz)
{
    dest = fin;
    dest.accrete(sz);
}

void
coarsen (BoxDomain&       dest,
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
    BL_PROFILE("BoxDomain::complementIn()");
    BoxList::complementIn(b,bl);
    BL_ASSERT(ok());
    return *this;
}

BoxDomain
complementIn (const Box&       b,
              const BoxDomain& bl)
{
    BoxDomain result;
    result.complementIn(b,bl);
    return result;
}

const BoxList&
BoxDomain::boxList () const
{
    return *this;
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
    :
    BoxList(IndexType::TheCellType())
{}

BoxDomain::BoxDomain (IndexType _ctype)
    :
    BoxList(_ctype)
{}

BoxDomain::BoxDomain (const Box& bx)
    :
    BoxList(bx)
{
}

void
BoxDomain::add (const Box& b)
{
    BL_ASSERT(b.ixType() == ixType());

    Vector<Box> tmp, check;

    check.push_back(b);

    for (const auto& bx : *this)
    {
        tmp.clear();
        for (auto& cbx : check)
        {
            if (cbx.intersects(bx))
            {
                const BoxList& tmpbl = amrex::boxDiff(cbx, bx);
                tmp.insert(std::end(tmp), std::begin(tmpbl), std::end(tmpbl));
                cbx = Box();
            }
        }
        check.erase(std::remove_if(check.begin(), check.end(),
                                   [](const Box& x) { return x.isEmpty(); }),
                    check.end());
        check.insert(std::end(check), std::begin(tmp), std::end(tmp));
    }
    join(check);
    BL_ASSERT(ok());
}

void
BoxDomain::add (const BoxList& bl)
{
    BoxList bl2 = bl;
    bl2.catenate(*this);
    join(amrex::removeOverlap(bl2));
}

BoxDomain&
BoxDomain::rmBox (const Box& b)
{
    BL_ASSERT(b.ixType() == ixType());

    Vector<Box> tmp;

    for (auto& bx : *this)
    {
        if (bx.intersects(b))
        {
            const BoxList& tmpbl = amrex::boxDiff(bx,b);
            tmp.insert(std::end(tmp), std::begin(tmpbl), std::end(tmpbl));
            bx = Box();
        }
    }
    removeEmpty();
    join(tmp);
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
        for (const_iterator bli = begin(); bli != end(); ++bli)
        {
            const_iterator blii = bli; ++blii;
            for ( ; blii != end(); ++blii)
            {
                if (bli->intersects(*blii))
                {
//		    amrex::Print(Print::AllProcs) << "Invalid DOMAIN, boxes overlap" << '\n'
//						  << "b1 = " << *bli << '\n'
//						  << "b2 = " << *blii << '\n';
                    status = false;
                }
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
        amrex::Error("operator<<(ostream&,BoxDomain&) failed");
    return os;
}

}
