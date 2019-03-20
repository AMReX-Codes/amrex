
#ifndef BL_BOXDOMAIN_H
#define BL_BOXDOMAIN_H

#include <iosfwd>

#include <AMReX_IndexType.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_BoxList.H>

namespace amrex
{
    class BoxDomain;

    //! Returns the complement of BoxDomain bl in Box b.
    BoxDomain complementIn (const Box& b, const BoxDomain& bl);

    /**
    * \brief Coarsen all Boxes in the domain by the refinement ratio.
    *           The result is placed into a new BoxDomain.
    */
    void coarsen (BoxDomain& dest, const BoxDomain& fin, int ratio);

    /**
    * \brief Refine all Boxes in the domain by the refinement ratio
    *           and return the result in dest.
    */
    void refine (BoxDomain&       dest,
                 const BoxDomain& fin,
                 int              ratio);

    /**
    * \brief Compute the intersection of BoxDomain fin with Box b
    *           and place the result into BoxDomain dest.
    */
    void intersect (BoxDomain&       dest,
                    const BoxDomain& fin,
                    const Box&       b);

    /**
    * \brief Grow each Box in BoxDomain fin by size sz and place
    *           the result into BoxDomain dest.
    */
    void accrete (BoxDomain&       dest,
                  const BoxDomain& fin,
                  int              sz = 1);


//! Output a BoxDomain to an ostream is ASCII format.
std::ostream& operator<< (std::ostream& os, const BoxDomain& bd);


/**
* \brief A List of Disjoint Boxes.
* A BoxDomain is a BoxList with the restriction that Boxes in the list
* are disjoint.
*/

//Note that a BoxDomain is NOT a BoxList due to the protected inheritance.
//This is a concrete class, not a polymorphic one.


class BoxDomain
    :
    protected BoxList
{
public:

    typedef BoxList::const_iterator const_iterator;
    //! Construct an empty BoxDomain of IndexType::TheCellType().
    BoxDomain ();

    //! A Box Domain out of single Box
    BoxDomain (const Box& bx);

    //! Construct an empty BoxDomain of IndexType itype.
    explicit BoxDomain (IndexType itype);

    BoxDomain (const BoxDomain& rhs) = default;
    BoxDomain (BoxDomain&& rhs) = default;
    BoxDomain& operator= (const BoxDomain& rhs) = default;
    BoxDomain& operator= (BoxDomain&& rhs) = default;

    //! Add a Box to the domain.
    void add (const Box& b);

    //! Add all Boxes in the BoxList to the domain.
    void add (const BoxList& bl);

    //! Remove a box from the domain.
    BoxDomain& rmBox (const Box& b);

    /**
    * Returns True if this BoxDomain is valid.  Checks for
    * disjointness. Returns true if the BoxDomain is empty.
    */
    bool ok () const;

    //! Promoted from BoxList.
    using BoxList::clear;

    //! Promoted from BoxList.
    using BoxList::size;

    //! Promoted from BoxList.
    using BoxList::isEmpty;

    //! Promoted from BoxList.
    using BoxList::isNotEmpty;

    //! Promoted from BoxList.
    using BoxList::minimalBox;

    //! Promoted from BoxList.
    using BoxList::ixType;

    //! Promoted from BoxList.
    using BoxList::simplify;

    //! Promoted from BoxList.
    using BoxList::begin;

    //! Promoted from BoxList.
    using BoxList::end;

    //! Promoted from BoxList
    using BoxList::contains;

    //! Promoted from BoxList
    using BoxList::shift;

    //! Promoted from BoxList
    using BoxList::shiftHalf;

    //! Return a const reference to the underlying BoxList of this BoxDomain.
    const BoxList& boxList () const;

    //! Are the BoxDomains equal?
    bool operator== (const BoxDomain& rhs) const;

    //! Are the BoxDomains not equal?
    bool operator!= (const BoxDomain& rhs) const;

    //! Refine all Boxes in the domain by the refinement ratio.
    BoxDomain& refine (int ratio);

    //! Coarsen all Boxes in the domain by the refinement ratio.
    BoxDomain& coarsen (int ratio);

    //! Intersect this BoxDomain with Box b.
    BoxDomain& intersect (const Box& b);

    //! Grow each Box by size sz.
    BoxDomain& accrete (int sz);

    //! Creates the complement of BoxDomain bl in Box b.
    BoxDomain& complementIn (const Box&       b,
                             const BoxDomain& bl);
};

}

#endif /*BL_BOXDOMAIN_H*/
