
#ifndef AMREX_TagBox_H_
#define AMREX_TagBox_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BaseFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>

namespace amrex {


/**
* \brief Tagged cells in a Box.
*
* This class is used to tag cells in a Box that need addition refinement.
*/

class BoxDomain;

class TagBox
    :
    public BaseFab<char>
{
public:

    //! The type of each tag.
    typedef char TagType;

    //! Possible values for each cell.
    enum TagVal { CLEAR=0, BUF, SET };

    TagBox () noexcept;

    explicit TagBox (Arena* ar) noexcept;

    TagBox (const Box& bx, int n, Arena* ar);

    explicit TagBox (const Box& bx, int n = 1, bool alloc = true,
                     bool shared = false, Arena* ar = nullptr);

    TagBox (const TagBox& rhs, MakeType make_type, int scomp, int ncomp);

    ~TagBox () = default;

    TagBox (TagBox&& rhs) noexcept = default;

    TagBox (const TagBox& rhs) = delete;
    TagBox& operator= (const TagBox& rhs) = delete;
    TagBox& operator= (TagBox&& rhs) = delete;


    /**
    * \brief Construct and return a new tagbox in which the coarsened cell
    * is tagged of any of the corresponding fine cells are tagged.
    */
    void coarsen (const IntVect& ratio, const Box& cbox) noexcept;

    /**
    * \brief Mark neighbors of every tagged cell a distance nbuff away
    *
    * \param nbuff
    */
    void buffer (const IntVect& nbuf) noexcept;

    /**
    * \brief Returns Vector\<int\> of size domain.numPts() suitable for calling
    * Fortran, with positions set to same value as in the TagBox
    * dataPtr().
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    Vector<int> tags () const noexcept;

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
    * \brief Set values as specified by the array -- this only tags.
    * It's an error if ar.length() != domain.numPts().
    *
    * \param ar
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    void tags (const Vector<int>& ar) noexcept;

    /**
    * \brief Set values as specified by the array -- this tags and untags.
    * It's an error if ar.length() != domain.numPts().
    *
    * \param ar
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    void tags_and_untags (const Vector<int>& ar) noexcept;

    /**
    * \brief Set values as specified by the array -- this only tags.
    * only changes values in the tilebx region
    *
    * \param ar
    * \param tilebx
    */
//#if (__cplusplus >= 201402L)
//    [[deprecated("No need to use this unless calling Fortran < 2003")]]
//#endif
    void tags (const Vector<int>& ar, const Box& tilebx) noexcept;

    /**
    * \brief Set values as specified by the array -- this tags and untags.
    * only changes values in the tilebx region
    *
    * \param ar
    * \param tilebx
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
    typedef TagBox::TagType TagType;

    /**
    * \brief The constructor.
    *
    * \param bs
    * \param dm
    * \param _ngrow
    */
    TagBoxArray (const BoxArray& bs, const DistributionMapping& dm, int _ngrow=0);
    TagBoxArray (const BoxArray& bs, const DistributionMapping& dm, const IntVect& _ngrow);

    ~TagBoxArray () override = default;

    TagBoxArray (TagBoxArray&& rhs) noexcept = default;
    TagBoxArray& operator= (TagBoxArray&& rhs) noexcept = default;

    TagBoxArray (const TagBoxArray& rhs) = delete;
    TagBoxArray& operator= (const TagBoxArray& rhs) = delete;

    /**
    * \brief Calls buffer() on all contained TagBoxes.
    *
    * \param nbuf
    */
    void buffer (const IntVect& nbuf);

    /**
    * \brief This funciton does two things.  Map tagged cells through a periodic boundary to other
    * grids in TagBoxArray cells, and remove duplicates.
    *
    * \param geom
    */
    void mapPeriodicRemoveDuplicates (const Geometry& geom);

    /**
    * \brief Set values in bl to val.
    *
    * \param bl
    * \param val
    */
    void setVal (const BoxList& bl, TagBox::TagVal val);

    /**
    * \brief Set values in bd to val.
    *
    * \param bd
    * \param val
    */
    void setVal (const BoxDomain& bd, TagBox::TagVal val);

    /**
    * \brief Set values in ba to val.
    *
    * \param ba
    * \param val
    */
    void setVal (const BoxArray& ba, TagBox::TagVal val);
    using FabArray<TagBox>::setVal;

    /**
    * \brief Calls coarsen() on all contained TagBoxes.
    *
    * \param ratio
    */
    void coarsen (const IntVect& ratio);

    /**
    * \brief Calls collate() on all contained TagBoxes.
    *
    * \param TheGlobalCollateSpace
    */
    void collate (Vector<IntVect>& TheGlobalCollateSpace) const;

    // \brief Are there tags in the region defined by bx?
    bool hasTags (Box const& bx) const;

    void local_collate_cpu (Vector<IntVect>& v) const;
#ifdef AMREX_USE_GPU
    void local_collate_gpu (Vector<IntVect>& v) const;
#endif
};

}

#endif /*_TagBox_H_*/
