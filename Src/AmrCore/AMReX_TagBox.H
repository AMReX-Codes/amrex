
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

    /**
    * \brief Construct an invalid TagBox with no memory.
    */
    TagBox () noexcept;

    /**
    * \brief Construct a TagBox on Box bx with number of components n.
    *
    * \param bx
    * \param n
    * \param alloc
    * \param shared
    */
    explicit TagBox (const Box& bx, int n = 1, bool alloc = true, bool shared = false);

    TagBox (const TagBox& rhs, MakeType make_type, int scomp, int ncomp);

    ~TagBox () = default;

    TagBox (TagBox&& rhs) noexcept = default;

    TagBox (const TagBox& rhs) = delete;
    TagBox& operator= (const TagBox& rhs) = delete;
    TagBox& operator= (TagBox&& rhs) = delete;


    /**
    * \brief Construct and return a new tagbox in which the coarsened cell
    * is tagged of any of the corresponding fine cells are tagged.
    *
    * \param ratio
    * \param owner
    */
    void coarsen (const IntVect& ratio, bool owner) noexcept;

    /**
    * \brief Mark neighbors of every tagged cell a distance nbuff away
    * only search interior for initial tagged points where nwid
    * is given as the width of the bndry region.
    *
    * \param nbuff
    * \param nwid
    */
    void buffer (const IntVect& nbuf, const IntVect& nwid) noexcept;

    /**
    * \brief Tag cells on intersect with src if corresponding src cell is tagged.
    *
    * \param src
    */
    void merge (const TagBox& src) noexcept;

    /**
    * \brief Add location of every tagged cell to IntVect array,
    * starting at given location.  Returns the number of
    * collated points.
    *
    * \param ar
    * \param start
    */
    long collate (Vector<IntVect>& ar, int start) const noexcept;

    /**
    * \brief Returns number of tagged cells in specified Box.
    *
    * \param bx
    */
    long numTags (const Box& bx) const noexcept;

    /**
    * \brief Returns total number of tagged cells in the TagBox.
    */
    long numTags () const noexcept;

    /**
    * \brief Returns Vector\<int\> of size domain.numPts() suitable for calling
    * Fortran, with positions set to same value as in the TagBox
    * dataPtr().
    */
    Vector<int> tags () const noexcept;

    /**
    * \brief Since a TagBox is a BaseFab\<char\>, we can use this utility
    * function to allocate an integer array to have the same number
    * of elements as cells in tilebx
    *
    * \param ar
    * \param tilebx
    */
    void get_itags(Vector<int>& ar, const Box& tilebx) const noexcept;

    /**
    * \brief Set values as specified by the array -- this only tags.
    * It's an error if ar.length() != domain.numPts().
    *
    * \param ar
    */
    void tags (const Vector<int>& ar) noexcept;

    /**
    * \brief Set values as specified by the array -- this tags and untags.
    * It's an error if ar.length() != domain.numPts().
    *
    * \param ar
    */
    void tags_and_untags (const Vector<int>& ar) noexcept;

    /**
    * \brief Set values as specified by the array -- this only tags.
    * only changes values in the tilebx region
    *
    * \param ar
    * \param tilebx
    */
    void tags (const Vector<int>& ar, const Box& tilebx) noexcept;

    /**
    * \brief Set values as specified by the array -- this tags and untags.
    * only changes values in the tilebx region
    *
    * \param ar
    * \param tilebx
    */
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

    TagBoxArray (const TagBoxArray& rhs) = delete;
    TagBoxArray& operator= (const TagBoxArray& rhs) = delete;
    TagBoxArray& operator= (TagBoxArray&& rhs) = delete;


    /**
    * \brief Returns the grow factor for the TagBoxArray.
    */
    IntVect borderSize () const noexcept;

    /**
    * \brief Calls buffer() on all contained TagBoxes.
    *
    * \param nbuf
    */
    void buffer (const IntVect& nbuf);

    /**
    * \brief Map tagged cells through a periodic boundary to other grids in
    * TagBoxArray cells which were outside domain are set to TagBox::CLEAR.
    *
    * \param geom
    */
    void mapPeriodic (const Geometry& geom);

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
    * \brief The total number of Tags in all the contained TagBoxes.
    */
    long numTags () const;

    /**
    * \brief Calls collate() on all contained TagBoxes.
    *
    * \param TheGlobalCollateSpace
    */
    void collate (Vector<IntVect>& TheGlobalCollateSpace) const;
};

}

#endif /*_TagBox_H_*/
