//BL_COPYRIGHT_NOTICE

#ifndef _TagBox_H_
#define _TagBox_H_

//
// $Id: TagBox.H,v 1.13 1998-07-29 19:09:09 lijewski Exp $
//

#include <IntVect.H>
#include <Box.H>
#include <Array.H>
#include <BaseFab.H>
#include <FabArray.H>
#include <BoxArray.H>
#include <BoxDomain.H>
#include <PArray.H>
#include <Geometry.H>

#if !defined(BL_SPACEDIM) || (BL_SPACEDIM != 1 && BL_SPACEDIM !=2 && BL_SPACEDIM != 3)
#error BL_SPACEDIM must be defined
#endif

#if (BL_SPACEDIM==1)
#define TagBox        TagBox1d
#elif (BL_SPACEDIM==2)
#define TagBox        TagBox2d
#elif (BL_SPACEDIM==3)
#define TagBox        TagBox3d
#endif

//
//@Man:
//@Memo: Tagged cells in a Box.
/*@Doc:

  This class is used to tag cells in a Box that need addition refinement.
*/

class TagBox
    :
    public BaseFab<char>
{
public:
    //
    //@ManDoc: The type of each tag.
    //
    typedef char TagType;
    //
    //@ManDoc: Possible values for each cell.
    //
    enum TagVal { CLEAR=0, BUF, SET };
    //
    //@ManDoc: Construct an invalid TagBox with no memory.
    //
    TagBox ();
    //
    //@ManDoc: Construct a TagBox on Box `bx' with number of components `n'.
    //
    TagBox (const Box& bx, int n = 1);
    //
    //@ManDoc: The destructor.
    //
    ~TagBox();
    //
    //@ManDoc: Resize the box on which the TagBox is defined.
    //
    void resize (const Box& b, int ncomp=1)
    {
        BaseFab<TagType>::resize(b,ncomp);
    }

    /*@ManDoc: Construct and return a new tagbox in which the coarsened cell
               is tagged of any of the corresponding fine cells are tagged.
    */
    TagBox* coarsen (const IntVect & ratio);

    /*@ManDoc: Mark neighbors of every tagged cell a distance nbuff away
               only search interior for initial tagged points where nwid
               is given as the width of the bndry region.
    */
    void buffer (int nbuff, int nwid);
    //
    //@ManDoc: Tag cells on intersect with src if corresponding src cell is tagged.
    //
    void merge (const TagBox& src);

    /*@ManDoc: Add location of every tagged cell to IntVect array,
               starting at given location.  Returns # of collated points.
    */
    int collate (IntVect* ar, int start) const;
    //
    //@ManDoc: Returns number of tagged cells in specified Box.
    //
    int numTags (const Box& bx) const;
    //
    //@ManDoc: Returns total number of tagged cells in the TagBox.
    //
    int numTags () const;

    /*@ManDoc: Returns Array<int> of size domain.numPts() suitable for calling
               Fortran, with positions set to same value as in the TagBox
               dataPtr().
    */
    Array<int> tags () const;

    /*@ManDoc: Set values as specified by the array.
               It's an error if ar.length() != domain.numPts().
    */
    void tags (const Array<int>& ar);

private:
    //
    // Disallowed.
    //
    TagBox (const TagBox& rhs);
    TagBox& operator= (const TagBox& rhs);
};

//
//@Man:
//@Memo: An array of TagBoxes.
/*@Doc:

  A container class for TagBoxes.
*/

class TagBoxArray
    :
    public FabArray<TagBox::TagType,TagBox>
{
public:
    //
    //@ManDoc: The type of each tag.
    //
    typedef TagBox::TagType TagType;
    //
    //@ManDoc: The constructor.
    //
    TagBoxArray (const BoxArray& bs, int _ngrow=0);
    //
    //@ManDoc: The destructor.
    //
    ~TagBoxArray ();
    //
    //@ManDoc: Returns the grow factor for the TagBoxArray.
    //
    int borderSize () { return m_border; }
    //
    //@ManDoc: Calls buffer() on all contained TagBoxes.
    //
    void buffer (int nbuf);

    /*@ManDoc: Map tagged cells through a periodic boundary to other grids in
               TagBoxArray cells which were outside domain are set to TagBox::CLEAR.
    */
    void mapPeriodic (const Geometry& geom);
    //
    //@ManDoc: Set values in `bd' to `val'.
    //
    void setVal (BoxDomain& bd, TagBox::TagVal val);
    //
    //@ManDoc: Set values in `ba' to `val'. 
    //
    void setVal (BoxArray& ba, TagBox::TagVal val);
    //
    //@ManDoc: Calls coarsen() on all contained TagBoxes.
    //
    void coarsen (const IntVect& ratio);
    //
    //@ManDoc: The total number of Tags in all the contained TagBoxes.
    //
    long numTags () const;
    //
    //@ManDoc: Calls collate() on all contained TagBoxes.
    //
    IntVect* collate (long& numtags) const;

private:
    //
    // Disallowed.
    //
    TagBoxArray ();
    TagBoxArray (const TagBoxArray& rhs);
    TagBoxArray& operator= (const TagBoxArray& rhs);
    //
    // Static space used by collate().
    //
    static IntVect* m_CollateSpace;
    //
    // Number of IntVects that can fit into m_CollateSpace
    //
    static long m_CollateCount;
    //
    // The initial size of `TagBoxArray::m_CollateCount' and the
    // amount we bump it each time it is not large enough.
    //
    enum { ChunkSize = 50000 };
    //
    // Helper functions.
    //
    static void BumpCollateSpace (long numtags);
    //
    // Not to be confused with FabArray::n_grow which defines
    // the invalid regions of the fabs in the FabArray.
    //
    int m_border;
};

#endif /*_TagBox_H_*/
