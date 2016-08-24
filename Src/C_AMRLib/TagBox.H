
#ifndef _TagBox_H_
#define _TagBox_H_

#include <IntVect.H>
#include <Box.H>
#include <Array.H>
#include <BaseFab.H>
#include <FabArray.H>
#include <BoxArray.H>
#include <PArray.H>
#include <Geometry.H>

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
    long collate (std::vector<IntVect>& ar, int start) const;
    //
    // Returns number of tagged cells in specified Box.
    //
    long numTags (const Box& bx) const;
    //
    // Returns total number of tagged cells in the TagBox.
    //
    long numTags () const;
    //
    // Returns Array<int> of size domain.numPts() suitable for calling
    // Fortran, with positions set to same value as in the TagBox
    // dataPtr().
    //
    Array<int> tags () const;
    //
    // Set values as specified by the array -- this only tags.
    // It's an error if ar.length() != domain.numPts().
    //
    void tags (const Array<int>& ar);
    //
    // Set values as specified by the array -- this tags and untags.
    // It's an error if ar.length() != domain.numPts().
    //
    void tags_and_untags (const Array<int>& ar);

    void get_itags(Array<int>& ar, const Box& tilebx) const;
    void tags (const Array<int>& ar, const Box& tilebx);
    void tags_and_untags (const Array<int>& ar, const Box&tilebx);

private:
    //
    // Disallowed.
    //
    TagBox (const TagBox& rhs);
    TagBox& operator= (const TagBox& rhs);
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
    TagBoxArray (const BoxArray& bs, int _ngrow=0);
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
    void collate (std::vector<IntVect>& TheGlobalCollateSpace) const;

    virtual void AddProcsToComp (int ioProcNumSCS, int ioProcNumAll,
                                 int scsMyId, MPI_Comm scsComm);

private:
    //
    // Disallowed.
    //
    TagBoxArray ();
    TagBoxArray (const TagBoxArray& rhs);
    TagBoxArray& operator= (const TagBoxArray& rhs);
};

#endif /*_TagBox_H_*/
