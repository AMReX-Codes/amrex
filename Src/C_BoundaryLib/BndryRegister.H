//BL_COPYRIGHT_NOTICE

#ifndef _BNDRYREGISTER_H_
#define _BNDRYREGISTER_H_

//
// $Id: BndryRegister.H,v 1.2 1997-12-10 21:56:01 lijewski Exp $
//

#include <BoxArray.H>
#include <FArrayBox.H>
#include <FabSet.H>
#include <Orientation.H>

//@Man:
/*@Memo:
        A BndryRegister organizes FabSets bounding each grid in a BoxArray.
	A FabSet is maintained for each boundary orientation, as well as
	the BoxArray domain of definition.
*/	
/*@Doc:
        A BndryRegister object contains a list of FabSets bounding the grids
	in a BoxArray.  The FabSet FABs are at the same refinement level
	as the grids they bound, and are accessed and modified via a variety
	of member functions.

	Non-default instantiation allocates a set of FABs, grown into and
	out of the bounding surfaces of each box in the BoxArray.  The width of
	the layer (in/out), as well as the "extent" of a bounding FABs (the
	number of nodes beyond the box boundaries, parallel to the grid
	surface) are determined by constructor argument.  All boxes and
	FABs in this context are assumed to be cell-centered.

	A small number of linear mathematical operations are provided for
	BndryRegisters, as well as overloaded [] operators for access based
	on grid boundary orientation.  The BoxArray domain of definition is
	accessible, but not resettable,
*/
class BndryRegister
{
protected:
    FabSet    bndry[2*BL_SPACEDIM];
    BoxArray  grids;

public:
//@ManMemo: administrative functions
//@ManDoc: default constructor
    BndryRegister();
//@ManDoc: constructor, given number of cells in/out, extent and number of components (assumes cell-centered boxes, and allocates cell-centered FABs)
    BndryRegister(const BoxArray& _grids, int _in_rad,
		  int _out_rad, int _extent_rad, int _ncomp);
//@ManDoc: construct from input stream
    BndryRegister(istream& is);
//@ManDoc: copy constructor
    BndryRegister(const BndryRegister &src);
//@ManDoc: copy operator
    BndryRegister& operator = (const BndryRegister& src);
//@ManDoc: destructor
    virtual ~BndryRegister();

//@ManMemo: I/O functions
//@ManDoc: write to stream, FAB-style
    ostream& writeOn(ostream &os) const;
//@ManDoc: write to stream
    friend ostream& operator << (ostream &os, const BndryRegister &br);
//@ManDoc: read from stream
    istream& readFrom(istream &is);

//@ManMemo: access functions
//@ManDoc: get box domain (as an array of boxes)
    const BoxArray& boxes() const
    {
	return grids;
    }
//@ManDoc: return the number of grids in this domain
    int length() const
    {
	return grids.length();
    }
//@ManDoc: return const set of FABs bounding the domain grid boxes on a given orientation
    const FabSet& operator[](const Orientation& _face) const
    {
	return bndry[_face];
    }
//@ManDoc: return set of FABs bounding the domain grid boxes on a given orientation
    FabSet& operator[](const Orientation& _face)
    {
	return bndry[_face];
    }


//@ManMemo: linear algebra functions on member data (over all orientations)
//@ManDoc: set all boundary FABs to given value
    void setVal(Real v);
//@ManDoc: fill the boundary FABs on intersect with given MultiFab
    BndryRegister& copyFrom(const MultiFab& src, int nghost, int src_comp,
			    int dest_comp, int num_comp);
//@ManDoc: increment the boundary FABs on intersect with given MultiFab
    BndryRegister& plusFrom(const MultiFab& src, int nghost, int src_comp,
			    int dest_comp, int num_comp);
//@ManDoc: fill the given MultiFab on intersect with the boundary FABs
    const BndryRegister& copyTo(MultiFab& dest, int nghost, int src_comp,
				int dest_comp, int num_comp) const;
//@ManDoc: Linear combination: this := a*mfa + b*mfb on intersection of MultiFabs with the boundary FABs
    BndryRegister& linComb(Real a, const MultiFab& mfa, int a_comp,
			   Real b, const MultiFab& mfb, int b_comp,
			   int dest_comp, int num_comp, int n_ghost=0);

//@ManMemo: internal administrative functions
//@ManDoc: set box domain, if not set previously
    void setBoxes(const BoxArray& _grids);
//@ManDoc: build FABs along given face
    void define(const Orientation& _face, const IndexType &_typ,
		int _in_rad,  int _out_rad, int _extent_rad,
		int _ncomp);

};

#endif

