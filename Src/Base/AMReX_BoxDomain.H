
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

    /**
     * \brief Returns the complement of BoxDomain \p bl in Box \p b.
     *
     * \param b Bounding box that defines the full region.
     * \param bl BoxDomain describing the portion to remove.
     * \return A BoxDomain covering the indices of \p b not present in \p bl.
     */
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


/**
 * \brief Output a BoxDomain to an ostream in ASCII format.
 *
 * \param os Stream receiving the textual representation.
 * \param bd Domain to print.
 * \return Reference to \p os to enable chaining.
 */
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

    //! Make a BoxDomain out of a single Box.
    BoxDomain (const Box& bx);

    //! Construct an empty BoxDomain of IndexType itype.
    explicit BoxDomain (IndexType itype);

    /**
     * \brief Insert a box into the domain, clipping it so the stored boxes remain disjoint.
     * \param b Box (matching the domain index type) to merge into the domain.
     */
    void add (const Box& b);

    /**
     * \brief Insert every box from \p bl, clipping to enforce disjointedness.
     * \param bl List of boxes (matching the domain index type) to merge.
     */
    void add (const BoxList& bl);

    /**
     * \brief Erase the intersection between this domain and \p b.
     * \param b Box to subtract.
     * \return Reference to this domain for chaining.
     */
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

    //! Return true if both domains cover identical regions.
    bool operator== (const BoxDomain& rhs) const;

    //! Logical negation of operator==.
    bool operator!= (const BoxDomain& rhs) const;

    /**
     * \brief Refine every box by the supplied ratio.
     * \param ratio Uniform refinement factor applied in each dimension.
     */
    BoxDomain& refine (int ratio);

    /**
     * \brief Coarsen every box by the supplied ratio.
     * \param ratio Uniform coarsening factor applied in each dimension.
     */
    BoxDomain& coarsen (int ratio);

    //! Restrict the domain to the overlap with box \p b.
    BoxDomain& intersect (const Box& b);

    //! Grow each stored box by \p sz cells on all sides.
    BoxDomain& accrete (int sz);

    /**
     * \brief Replace this domain with the indices in \p b that are not present in \p bl.
     * \param b Bounding box that defines the superset region.
     * \param bl Domain of indices to subtract.
     * \return Reference to this domain for chaining.
     */
    BoxDomain& complementIn (const Box&       b,
                             const BoxDomain& bl);

private:
    //! Construct from a BoxList that has no overlap.
    explicit BoxDomain (BoxList&& bl);
};

}

#endif /*BL_BOXDOMAIN_H*/
