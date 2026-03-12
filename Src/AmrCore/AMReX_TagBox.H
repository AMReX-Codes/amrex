
#ifndef AMREX_TagBox_H_
#define AMREX_TagBox_H_
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BaseFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

/**
 * \file AMReX_TagBox.H
 * \brief Tag storage used for AMR error estimation and refinement decisions.
 */

namespace amrex {


/**
* \brief Tagged cells in a Box.
*
* This class is used to tag cells in a Box that need addition refinement.
*/
class TagBox final
    :
    public BaseFab<char>
{
public:

    //! The type of each tag.
    using TagType = char;

    //! Possible values for each cell.
    enum TagVal { CLEAR=0, BUF, SET };

    TagBox () noexcept = default;

    explicit TagBox (Arena* ar) noexcept;

    TagBox (const Box& bx, int n, Arena* ar);

    explicit TagBox (const Box& bx, int n = 1, bool alloc = true,
                     bool shared = false, Arena* ar = nullptr);

    TagBox (const TagBox& rhs, MakeType make_type, int scomp, int ncomp);

    ~TagBox () noexcept override = default;

    TagBox (TagBox&& rhs) noexcept = default;

    TagBox (const TagBox& rhs) = delete;
    TagBox& operator= (const TagBox& rhs) = delete;
    TagBox& operator= (TagBox&& rhs) = delete;


    /**
    * \brief Coarsen the TagBox in-place and mark coarse cells as tagged if any fine cell was tagged.
    *
    * \param ratio Refinement ratio used to coarsen the data.
    * \param cbox  Coarse box describing the region that should remain after coarsening.
    */
    void coarsen (const IntVect& ratio, const Box& cbox) noexcept;

    /**
    * \brief Grow tags by \p a_nbuff cells, ignoring seed tags within \p nwid of the boundary.
    *
    * \param a_nbuff Number of buffer cells to add in each direction.
    * \param nwid    Width of the boundary region to examine for seed tags.
    */
    void buffer (const IntVect& a_nbuff, const IntVect& nwid) noexcept;

    /**
    * \brief Returns Vector\<int\> of size domain.numPts() suitable for calling
    * Fortran, with positions set to same value as in the TagBox
    * dataPtr().
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    [[nodiscard]] Vector<int> tags () const noexcept;

    /**
    * \brief Since a TagBox is a BaseFab\<char\>, we can use this utility
    * function to allocate an integer array to have the same number
    * of elements as cells in tilebx
    *
    * \param ar
    * \param tilebx
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    void get_itags(Vector<int>& ar, const Box& tilebx) const noexcept;

    /**
    * \brief Set values as specified by \p ar, tagging cells only.
    *
    * \param ar Flat array of tag values; must match the size of the TagBox.
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    void tags (const Vector<int>& ar) noexcept;

    /**
    * \brief Set values as specified by \p ar, allowing both tagging and untagging.
    *
    * \param ar Flat array of tag values; must match the size of the TagBox.
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    void tags_and_untags (const Vector<int>& ar) noexcept;

    /**
    * \brief Same as tags(Vector<int>) but only modifies the region \p tilebx.
    *
    * \param ar      Flat array of tag values.
    * \param tilebx  Box describing the subregion to update.
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    void tags (const Vector<int>& ar, const Box& tilebx) noexcept;

    /**
    * \brief Same as tags_and_untags(Vector<int>) but restricted to \p tilebx.
    *
    * \param ar      Flat array of tag values.
    * \param tilebx  Box describing the subregion to update.
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    void tags_and_untags (const Vector<int>& ar, const Box& tilebx) noexcept;
};


/**
* \brief An array of TagBoxes.
*
* A container class for TagBoxes.
*/
class TagBoxArray
    :
    public FabArray<TagBox>
{
public:

    //! The type of each tag.
    using TagType = TagBox::TagType;

    /**
    * \brief Construct a TagBoxArray on BoxArray \p ba / DistributionMapping \p dm with grow vector \p _ngrow.
    */
    TagBoxArray (const BoxArray& ba, const DistributionMapping& dm, int _ngrow=0);
    TagBoxArray (const BoxArray& ba, const DistributionMapping& dm, const IntVect& _ngrow);

    ~TagBoxArray () = default;

    TagBoxArray (TagBoxArray&& rhs) noexcept = default;
    TagBoxArray& operator= (TagBoxArray&& rhs) noexcept = default;

    TagBoxArray (const TagBoxArray& rhs) = delete;
    TagBoxArray& operator= (const TagBoxArray& rhs) = delete;

    /**
    * \brief Call TagBox::buffer() on every FAB in the array with grow vector \p nbuf.
    */
    void buffer (const IntVect& nbuf);

    /**
    * \brief Map tags across periodic boundaries and remove duplicates using \p geom.
    */
    void mapPeriodicRemoveDuplicates (const Geometry& geom);

    /**
    * \brief Set all tags inside BoxArray \p ba to value \p val.
    */
    void setVal (const BoxArray& ba, TagBox::TagVal val);
    using FabArray<TagBox>::setVal;

    /**
     * \brief Calls TagBox::coarsen() on every FAB using ratio \p ratio.
     *
     * \param ratio Refinement ratio to apply when coarsening the tag layout.
     */
    void coarsen (const IntVect& ratio);

    /**
     * \brief Calls TagBox::collate() across all FABs, storing indices in \p TheGlobalCollateSpace.
     *
     * \param TheGlobalCollateSpace Output pinned vector that receives the global list of tag locations.
     */
    void collate (Gpu::PinnedVector<IntVect>& TheGlobalCollateSpace) const;

    //! Return true if any cells inside \p bx are tagged (SET or BUF), false otherwise.
    bool hasTags (Box const& bx) const;

    void local_collate_cpu (Gpu::PinnedVector<IntVect>& v) const;
#ifdef AMREX_USE_GPU
    void local_collate_gpu (Gpu::PinnedVector<IntVect>& v) const;
#endif
};

}

#endif /*_TagBox_H_*/
