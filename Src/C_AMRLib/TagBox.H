//BL_COPYRIGHT_NOTICE

#ifndef _TagBox_H_
#define _TagBox_H_

//
// $Id: TagBox.H,v 1.5 1997-12-11 23:27:54 lijewski Exp $
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

class TagBox
    :
    public BaseFab<int>
{
public:
    //
    // The type of each tag.
    //
    typedef int TagType;
    //
    // Possible values for each cell.
    //
    enum TagVal { CLEAR=0, BUF, SET };
    //
    // Construct an invalid TagBox with no memory.
    //
    TagBox ();

    TagBox (const Box& bx, int n = 1);

    virtual ~TagBox();
    
    void resize (const Box& b) { BaseFab<TagType>::resize(b,1); }

    TagType& operator() (const IntVect& p)
    {
        return dptr[domain.index(p)];
    }
    const TagType& operator() (const IntVect& p) const
    {
        return dptr[domain.index(p)];
    }
    //
    // Construct and return a new tagbox in which the coarsened cell is
    // tagged of any of the corresponding fine cells are tagged.
    //
    TagBox* coarsen (const IntVect & ratio);
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
    // starting at given location.
    //
    int colate (Array<IntVect> &ar, int start) const;

    int numTags (const Box& bx) const;
    int numTags () const;

private:
    //
    // Disallowed.
    //
    TagBox (const TagBox& rhs);
    TagBox& operator= (const TagBox& rhs);
};

class TagBoxArray
    :
    public FabArray<TagBox::TagType,TagBox>
{
public:
    //
    // The type of each tag.
    //
    typedef TagBox::TagType TagType;

    TagBoxArray (const BoxArray& bs, int _ngrow=0);

    virtual ~TagBoxArray ();

    int borderSize () { return m_border; }
    void buffer (int nbuf);
    void mergeUnique ();
    //
    // Map tagged cells through a periodic boundary to other grids in
    // TagBoxArray cells which were outside domain are set to TagBox::CLEAR.
    //
    void mapPeriodic (const Geometry& geom);

    void setVal (BoxDomain& bd, TagBox::TagVal val);
    void setVal (BoxArray& ba, TagBox::TagVal val);
    void coarsen (const IntVect& ratio);

    long numTags() const;

    Array<IntVect>* colate() const;

private:
    //
    // Disallowed.
    //
    TagBoxArray ();
    TagBoxArray (const TagBoxArray& rhs);
    TagBoxArray& operator= (const TagBoxArray& rhs);
    //
    // Not to be confused with FabArray::n_grow which defines
    // the invalid regions of the fabs in the FabArray.
    //
    int m_border;
};

#endif /*_TagBox_H_*/
