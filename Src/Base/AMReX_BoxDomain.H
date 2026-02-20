
#ifndef BL_BOXDOMAIN_H
#define BL_BOXDOMAIN_H
#include <AMReX_Config.H>

#include <AMReX_IndexType.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_BoxList.H>

#include <iosfwd>

/**
 * \file AMReX_BoxDomain.H
 * \brief Disjoint-list representation of boxes and associated utilities.
 */

namespace amrex
{
    class BoxDomain;
    class BoxArray;

    //! Returns the complement of BoxDomain bl in Box b.
    BoxDomain complementIn (const Box& b, const BoxDomain& bl);

    /**
     * \brief Coarsen all Boxes in the domain by the refinement ratio.
     *
     * \param dest Destination domain to populate.
     * \param fin Source domain to coarsen.
     * \param ratio Uniform coarsening factor.
     */
    void coarsen (BoxDomain& dest, const BoxDomain& fin, int ratio);

    /**
     * \brief Refine all Boxes in the domain by the refinement ratio.
     *
     * \param dest Destination domain to populate.
     * \param fin Source domain to refine.
     * \param ratio Uniform refinement factor.
     */
    void refine (BoxDomain&       dest,
                 const BoxDomain& fin,
                 int              ratio);

    /**
     * \brief Compute the intersection of BoxDomain fin with Box b.
     *
     * \param dest Destination domain to populate.
     * \param fin Source domain.
     * \param b Clipping box.
     */
    void intersect (BoxDomain&       dest,
                    const BoxDomain& fin,
                    const Box&       b);

    /**
     * \brief Grow each Box in BoxDomain fin by size sz.
     *
     * \param dest Destination domain to populate.
     * \param fin Source domain.
     * \param sz Number of cells to grow on all sides.
     */
    void accrete (BoxDomain&       dest,
                  const BoxDomain& fin,
                  int              sz = 1);


//! Output a BoxDomain to an ostream is ASCII format.
std::ostream& operator<< (std::ostream& os, const BoxDomain& bd);


/**
 * \class amrex::BoxDomain
 * \ingroup amrex_domain_decomposition
 * \brief A list of disjoint boxes built on top of `BoxList`.
 *
 * The underlying `BoxList` is exposed selectively through promoted members.
 */
class BoxDomain
    :
    protected BoxList
{
public:
    friend class ClusterList;

    using const_iterator = BoxList::const_iterator;

    //! Construct an empty BoxDomain of IndexType::TheCellType().
    BoxDomain ();

    //! Make a domain out of a single Box.
    BoxDomain (const Box& bx);

    //! Construct an empty BoxDomain of IndexType itype.
    explicit BoxDomain (IndexType itype);

    //! Add a Box to the domain, clipping any overlap to keep boxes disjoint.
    void add (const Box& b);

    //! Add all Boxes in the BoxList to the domain, clipping overlaps as needed.
    void add (const BoxList& bl);

    //! Remove a box from the domain.
    BoxDomain& rmBox (const Box& b);

    /**
     * \brief Returns true if this BoxDomain is valid (disjoint or empty).
     */
    [[nodiscard]] bool ok () const;

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
    [[nodiscard]] const BoxList& boxList () const;

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

private:
    //! Construct from a BoxList that has no overlap.
    explicit BoxDomain (BoxList&& bl);
};

}

#endif /*BL_BOXDOMAIN_H*/
