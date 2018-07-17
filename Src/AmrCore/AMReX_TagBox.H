
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

//
// Tagged cells in a Box.
//
// This class is used to tag cells in a Box that need addition refinement.
//

class BoxDomain;

class TagBox
    :
    public BaseFab<char>
{
public:
    //
    // The type of each tag.
    //
    typedef char TagType;
    //
    // Possible values for each cell.
    //
    enum TagVal { CLEAR=0, BUF, SET };
    //
    // Construct an invalid TagBox with no memory.
    //
    TagBox ();
    //
    // Construct a TagBox on Box bx with number of components n.
    //
    explicit TagBox (const Box& bx, int n = 1, bool alloc = true, bool shared = false);

    ~TagBox () = default;
    TagBox (TagBox&& rhs) noexcept = default;

    TagBox (const TagBox& rhs) = delete;
    TagBox& operator= (const TagBox& rhs) = delete;
    TagBox& operator= (TagBox&& rhs) = delete;

    //
    // Construct and return a new tagbox in which the coarsened cell
    // is tagged of any of the corresponding fine cells are tagged.
    //
    void coarsen (const IntVect& ratio, bool owner);
    //
    // Mark neighbors of every tagged cell a distance nbuff away
    // only search interior for initial tagged points where nwid
    // is given as the width of the bndry region.
    //
    void buffer (int nbuff, int nwid);
    //
    // Tag cells on intersect with src if corresponding src cell is tagged.
    //
    void merge (const TagBox& src);
    //
    // Add location of every tagged cell to IntVect array,
    // starting at given location.  Returns the number of
    // collated points.
    //
    long collate (Vector<IntVect>& ar, int start) const;
    //
    // Returns number of tagged cells in specified Box.
    //
    long numTags (const Box& bx) const;
    //
    // Returns total number of tagged cells in the TagBox.
    //
    long numTags () const;
    //
    // Returns Vector<int> of size domain.numPts() suitable for calling
    // Fortran, with positions set to same value as in the TagBox
    // dataPtr().
    //
    Vector<int> tags () const;
    //
    // Since a TagBox is a BaseFab<char>, we can use this utility
    // function to allocate an integer array to have the same number 
    // of elements as cells in tilebx
    //
    void get_itags(Vector<int>& ar, const Box& tilebx) const;
    //
    // Set values as specified by the array -- this only tags.
    // It's an error if ar.length() != domain.numPts().
    //
    void tags (const Vector<int>& ar);
    //
    // Set values as specified by the array -- this tags and untags.
    // It's an error if ar.length() != domain.numPts().
    //
    void tags_and_untags (const Vector<int>& ar);
    //
    // Set values as specified by the array -- this only tags.
    // only changes values in the tilebx region
    // 
    void tags (const Vector<int>& ar, const Box& tilebx);
    //
    // Set values as specified by the array -- this tags and untags.
    // only changes values in the tilebx region
    //
    void tags_and_untags (const Vector<int>& ar, const Box& tilebx);
};

//
// An array of TagBoxes.
//
// A container class for TagBoxes.
//

class TagBoxArray
    :
    public FabArray<TagBox>
{
public:
    //
    // The type of each tag.
    //
    typedef TagBox::TagType TagType;
    //
    // The constructor.
    //
    TagBoxArray (const BoxArray& bs, const DistributionMapping& dm, int _ngrow=0);

    ~TagBoxArray () override = default;

    TagBoxArray (TagBoxArray&& rhs) noexcept = default;

    TagBoxArray (const TagBoxArray& rhs) = delete;
    TagBoxArray& operator= (const TagBoxArray& rhs) = delete;
    TagBoxArray& operator= (TagBoxArray&& rhs) = delete;

    //
    // Returns the grow factor for the TagBoxArray.
    //
    int borderSize () const;
    //
    // Calls buffer() on all contained TagBoxes.
    //
    void buffer (int nbuf);
    //
    // Map tagged cells through a periodic boundary to other grids in
    // TagBoxArray cells which were outside domain are set to TagBox::CLEAR.
    //
    void mapPeriodic (const Geometry& geom);
    //
    // Set values in bl to val.
    //
    void setVal (const BoxList& bl, TagBox::TagVal val);
    //
    // Set values in bd to val.
    //
    void setVal (const BoxDomain& bd, TagBox::TagVal val);
    //
    // Set values in ba to val. 
    //
    void setVal (const BoxArray& ba, TagBox::TagVal val);
    using FabArray<TagBox>::setVal;
    //
    // Calls coarsen() on all contained TagBoxes.
    //
    void coarsen (const IntVect& ratio);
    //
    // The total number of Tags in all the contained TagBoxes.
    //
    long numTags () const;
    //
    // Calls collate() on all contained TagBoxes.
    //
    void collate (Vector<IntVect>& TheGlobalCollateSpace) const;
};

}

#endif /*_TagBox_H_*/
