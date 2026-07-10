
#ifndef BL_BOXLIST_H
#define BL_BOXLIST_H
#include <AMReX_Config.H>

/**
 * \file AMReX_BoxList.H
 * \brief A list of Boxes sharing a common IndexType.
 */

#include <AMReX_IntVect.H>
#include <AMReX_IndexType.H>
#include <AMReX_Box.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

#include <iosfwd>

namespace amrex
{
    class BoxArray;
    class BoxList;

    /**
     * \brief Returns a BoxList defining the complement of BoxList \p bl in Box \p b.
     *
     * \param b Bounding box that defines the full domain.
     * \param bl BoxList describing the region to exclude.
     * \return A BoxList covering the indices of \p b not present in \p bl.
     */
    BoxList complementIn (const Box& b, const BoxList& bl);

    /**
     * \brief Returns a BoxList defining the complement of \p b2 in \p b1in.
     *
     * \param b1in Base box.
     * \param b2 Box describing the region to subtract.
     */
    BoxList boxDiff (const Box& b1in, const Box& b2);

    /**
     * \brief Populate \p bl_diff with the complement of \p b2 in \p b1in.
     *
     * \param bl_diff Output BoxList receiving the difference.
     * \param b1in Base box.
     * \param b2 Box describing the region to subtract.
     */
    void boxDiff (BoxList& bl_diff, const Box& b1in, const Box& b2);

    /**
     * \brief Returns a new BoxList in which each Box is refined by the given ratio.
     *
     * \param bl Input BoxList.
     * \param ratio Uniform refinement ratio applied to every box.
     */
    BoxList refine (const BoxList&  bl, int ratio);

    /**
     * \brief Returns a new BoxList in which each Box is coarsened by the given ratio.
     *
     * \param bl Input BoxList.
     * \param ratio Uniform coarsening ratio applied to every box.
     */
    BoxList coarsen (const BoxList& bl, int ratio);

    /**
     * \brief Returns a BoxList defining the intersection of \p bl with \p b.
     *
     * \param bl Input BoxList.
     * \param b Clipping box.
     */
    BoxList intersect (const BoxList& bl, const Box& b);

    /**
     * \brief Returns a new BoxList in which each Box is grown by the given size.
     *
     * \param bl Input BoxList.
     * \param sz Number of cells to grow in every direction.
     */
    BoxList accrete (const BoxList& bl, int sz);

    /**
     * \brief Return BoxList which covers the same area but has no overlapping boxes.
     *
     * \param bl Input BoxList (possibly overlapping).
     */
    BoxList removeOverlap (const BoxList& bl);


    /**
     * \brief Output a BoxList to an ostream in ASCII format.
     *
     * \param os Stream receiving the text output.
     * \param blist BoxList to print.
     * \return Reference to \p os to enable chaining.
     */
    std::ostream& operator<< (std::ostream&  os, const BoxList& blist);


/**
 * \class amrex::BoxList
 * \ingroup amrex_index_space
 * \brief A list of Boxes sharing a common IndexType.
 *
 * BoxList implements complement, intersection, refinement, coarsening, and
 * other operations on collections of Boxes.
 */
class BoxList
{
public:

    friend class BoxDomain;

    using       iterator = Vector<Box>::iterator;
    using const_iterator = Vector<Box>::const_iterator;

    //! Construct an empty BoxList with IndexType::TheCellType().
    BoxList ();

    //! Make a boxList consisting of a single Box
    explicit BoxList (const Box& bx);

    //! Construct an empty BoxList with IndexType _btype
    explicit BoxList (IndexType _btype);

    //!  Construct a BoxList from a BoxArray.
    explicit BoxList (const BoxArray& ba);

    //! Construct a boxlist from a rvalue Vector<Box>.
    explicit BoxList (Vector<Box>&& bxs);

    //! Make a boxList consisting of non-overlapping tile box from a single Box
    BoxList (const Box& bx, const IntVect& tilesize);

    //! Make a BoxList consisting of nboxes Boxes from a single Box
    BoxList (const Box& bx, int nboxes);

    //! Make a BoxList consisting of nboxes Boxes chopped in dir-direction from a single Box
    BoxList (const Box& bx, int nboxes, Direction dir);

    ~BoxList () = default;
    BoxList (const BoxList& rhs) = default;
    BoxList (BoxList&& rhs) = default;
    BoxList& operator= (const BoxList& rhs) = default;
    BoxList& operator= (BoxList&& rhs) = default;

    //! Reserve storage for at least \p n Boxes.
    void reserve (std::size_t n) { m_lbox.reserve(n); }

    //! Append a Box to this BoxList.
    void push_back (const Box& bn) {
        if (m_lbox.empty()) { btype = bn.ixType(); }
        BL_ASSERT(ixType() == bn.ixType());
        m_lbox.push_back(bn);
    }

    //! Return a reference to the first Box in the list.
    [[nodiscard]] Box& front () noexcept { BL_ASSERT(!m_lbox.empty()); return m_lbox.front(); }

    //! Return a const reference to the first Box in the list.
    [[nodiscard]] const Box& front () const noexcept { BL_ASSERT(!m_lbox.empty()); return m_lbox.front(); }

    //! Join the BoxList to ourselves.
    void join (const BoxList& blist);
    //! Join the Array of Boxes to ourselves.
    void join (const Vector<Box>& barr);
    //! Catenate the BoxList to ourselves. Removes entries from blist.
    void catenate (BoxList& blist);

    //! Remove all Boxes from this BoxList.
    void clear ();
    //! The number of Boxes in this BoxList.
    [[nodiscard]] Long size () const noexcept { return m_lbox.size(); }

    //! Return the current allocated capacity.
    [[nodiscard]] std::size_t capacity () const noexcept { return m_lbox.capacity(); }

    //! Return an iterator to the first Box.
    [[nodiscard]] iterator begin () noexcept { return m_lbox.begin(); }
    //! Return a const iterator to the first Box.
    [[nodiscard]] const_iterator begin () const noexcept { return m_lbox.cbegin(); }
    //! Return a const iterator to the first Box.
    [[nodiscard]] const_iterator cbegin () const noexcept { return m_lbox.cbegin(); }

    //! Return an iterator past the last Box.
    [[nodiscard]] iterator end () noexcept { return m_lbox.end(); }
    //! Return a const iterator past the last Box.
    [[nodiscard]] const_iterator end () const noexcept { return m_lbox.cend(); }
    //! Return a const iterator past the last Box.
    [[nodiscard]] const_iterator cend () const noexcept { return m_lbox.cend(); }
    /**
    * \brief True if this BoxList is valid; i.e. all the Boxes are
    * valid and they all have the same IndexType.  Also returns
    * true if the BoxList is empty.
    */
    [[nodiscard]] bool ok () const noexcept;
    //! Is this BoxList equal to rhs?
    bool operator== (const BoxList& rhs) const;
    //! Is this BoxList notequal to rhs?
    bool operator!= (const BoxList& rhs) const;
    //! Is this BoxList empty?
    [[nodiscard]] bool isEmpty () const noexcept { return m_lbox.empty(); }
    //! Is this BoxList not empty?
    [[nodiscard]] bool isNotEmpty () const noexcept { return !m_lbox.empty(); }
    //! True if the set of intersecting Boxes is empty.
    [[nodiscard]] bool isDisjoint () const;
    //! True if all Boxes in bl are contained in this BoxList.
    [[nodiscard]] bool contains (const BoxList& bl) const;
    //! Modify this BoxList to contain only its intersection with Box b.
    BoxList& intersect (const Box& b);
    //! Modify this BoxList to contain only its intersection with BoxList bl.
    BoxList& intersect (const BoxList& bl);
    //! Remove empty Boxes from this BoxList.
    BoxList& removeEmpty();

    /**
     * \brief Set this BoxList to the complement of BoxList \p bl in Box \p b.
     *
     * \param b Bounding box that defines the full domain.
     * \param bl Region to exclude.
     * \return Reference to \c *this.
     */
    BoxList& complementIn (const Box& b, const BoxList& bl);
    /**
     * \brief Set this BoxList to the complement of BoxList \p bl (moved) in Box \p b.
     *
     * \param b Bounding box that defines the full domain.
     * \param bl Region to exclude (consumed by the call).
     * \return Reference to \c *this.
     */
    BoxList& complementIn (const Box& b, BoxList&& bl);
    /**
     * \brief Set this BoxList to the complement of BoxArray \p ba in Box \p b.
     *
     * \param b Bounding box that defines the full domain.
     * \param ba Region to exclude (provided as a BoxArray).
     * \return Reference to \c *this.
     */
    BoxList& complementIn (const Box& b, const BoxArray& ba);
    /**
     * \brief Parallel version of complementIn using BoxList.
     *
     * \param b Bounding box that defines the full domain.
     * \param bl Region to exclude.
     * \return Reference to \c *this.
     */
    BoxList& parallelComplementIn (const Box& b, const BoxList& bl);
    /**
     * \brief Parallel version of complementIn using BoxList (moved).
     *
     * \param b Bounding box that defines the full domain.
     * \param bl Region to exclude (consumed by the call).
     * \return Reference to \c *this.
     */
    BoxList& parallelComplementIn (const Box& b, BoxList&& bl);
    /**
     * \brief Parallel version of complementIn using BoxArray.
     *
     * \param b Bounding box that defines the full domain.
     * \param ba Region to exclude (provided as BoxArray).
     * \return Reference to \c *this.
     */
    BoxList& parallelComplementIn (const Box& b, const BoxArray& ba);

    /**
     * \brief Refine each Box in the BoxList by a uniform ratio.
     *
     * \param ratio Refinement factor applied in every direction.
     */
    BoxList& refine (int ratio);
    /**
     * \brief Refine each Box in the BoxList by direction-dependent ratios.
     *
     * \param ratio Component-wise refinement factors.
     */
    BoxList& refine (const IntVect& ratio);
    /**
     * \brief Coarsen each Box in the BoxList by a uniform ratio.
     *
     * \param ratio Coarsening factor applied in every direction.
     */
    BoxList& coarsen (int ratio);
    /**
     * \brief Coarsen each Box in the BoxList by direction-dependent ratios.
     *
     * \param ratio Component-wise coarsening factors.
     */
    BoxList& coarsen (const IntVect& ratio);
    /**
     * \brief Grow each Box in the BoxList by \p sz cells in every direction.
     *
     * \param sz Cells added on each side in every direction (negative values shrink).
     */
    BoxList& accrete  (int sz);
    /**
     * \brief Grow each Box in the BoxList by \p sz cells (per direction).
     *
     * \param sz Direction-wise growth vector.
     */
    BoxList& accrete  (const IntVect& sz);
    /**
     * \brief Apply `Box::shift(dir, nzones)` to every Box in the list.
     *
     * \param dir Coordinate direction to shift.
     * \param nzones Number of cell widths to shift (positive or negative).
     */
    BoxList& shift (int dir,
                    int nzones);
    /**
     * \brief Apply `Box::shiftHalf(dir, num_halfs)` to every Box in the list.
     *
     * \param dir Coordinate direction to shift.
     * \param num_halfs Number of half-cell shifts.
     */
    BoxList& shiftHalf (int dir,
                        int num_halfs);
    /**
     * \brief Apply `Box::shiftHalf(const IntVect&)` to every Box in the list.
     *
     * \param iv Vector of half-cell shifts per direction.
     */
    BoxList& shiftHalf (const IntVect& iv);
    /**
    * \brief Merge adjacent Boxes in this BoxList. Return the number
    * of Boxes merged.  If "best" is specified we do a single
    * bruteforce pass over the list checking each Box against
    * all Boxes after it in the list to see if they can be
    * merged.  If "best" is not specified we limit how fair
    * afield we look for possible matches.  The "best" algorithm
    * is O(N-squared) while the other algorithm is O(N).
    */
    int simplify (bool best = false);
    //! Merge adjacent Boxes assuming the list is already sorted/ordered.
    int ordered_simplify ();
    //! Forces each Box in the BoxList to have sides of length <= chunk.
    BoxList& maxSize (int chunk);
    //! Forces each Box in the BoxList to have dimth side of length <= chunk[dim].
    BoxList& maxSize (const IntVect& chunk);
    //! Returns smallest Box that contains all Boxes in this BoxList.
    [[nodiscard]] Box minimalBox () const;
    //! Returns the IndexType of Boxes in this BoxList.
    [[nodiscard]] IndexType ixType () const noexcept { return btype; }
    //! Set the type of the BoxList.  It's an error if the BoxList isn't empty.
    void set (IndexType ixtyp) noexcept { BL_ASSERT(m_lbox.empty()); btype = ixtyp; }
    /**
    * \brief Applies surroundingNodes(Box) to each Box in BoxArray.
    * See the documentation of Box for details.
    */
    BoxList& surroundingNodes () noexcept;
    /**
    * \brief Applies surroundingNodes(Box,int) to each Box in
    * BoxList.  See the documentation of Box for details.
    */
    BoxList& surroundingNodes (int dir) noexcept;
    //! Applies Box::enclosedCells() to each Box in the BoxList.
    BoxList& enclosedCells () noexcept;
    //! Applies Box::enclosedCells(int) to each Box in the BoxList.
    BoxList& enclosedCells  (int dir) noexcept;
    //! Applies Box::convert(IndexType) to each Box in the BoxList.
    BoxList& convert (IndexType typ) noexcept;

    //! Returns a reference to the Vector<Box>.
    [[nodiscard]] Vector<Box>& data () noexcept { return m_lbox; }
    //! Returns a constant reference to the Vector<Box>.
    [[nodiscard]] const Vector<Box>& data () const noexcept { return m_lbox; }

    //! Swap contents with another BoxList.
    void swap (BoxList& rhs) noexcept {
        std::swap(m_lbox, rhs.m_lbox);
        std::swap(btype, rhs.btype);
    }

    //! Broadcast this BoxList from the I/O processor to all other processors.
    void Bcast ();

private:
    //! Core simplify routine.
    int simplify_doit (int depth);

    //! The list of Boxes.
    Vector<Box> m_lbox;
    //! The IndexType of Boxes in the BoxList.
    IndexType btype;

};

}

#endif /*BL_BOXLIST_H*/
