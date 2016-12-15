
#ifndef _MCINTERPBNDRYDATA_H_
#define _MCINTERPBNDRYDATA_H_

#include <BLassert.H>
#include <BoxArray.H>
#include <FArrayBox.H>
#include <Orientation.H>
#include <Mask.H>
#include <BCRec.H>
#include <Geometry.H>
#include <BoundCond.H>
#include <BndryData.H>


/*
        An MCInterpBndryData object adds to a BndryData object the ability to
	manipulate and set the data stored in the boundary cells.

        The "Interpbndrydata" class is a virtual base class derived from
	BndryData.  It is intended to provide a more physical method for
	filling boundary-related data.  Boundary values in a BndryData object
	are stored in FabSets around each grid in the domain, and the
	MCInterpBndryData class provides a mechanism for filling these FabSets,
	consistent with AMR-like numerical discretizations.  When asked to
	set it's boundary values, an MCInterpBndryData object:

	Fills with physical boundary values if the FAB is on the
	domain boundary - the corresponding values are presumed to be
	stored in the ghost cells of a MultiFab given to the boundary filling
	routine

	Fills on intersection with data from the VALID region of the
	input MultiFab, ensuring that adjacent FABs use consistent data at
	their intersection, and otherwise,

	Fills with values interpolated from a coarser FAB that
	bounds the cells that do not meet the above two criteria
*/
class MCInterpBndryData
    :
    public BndryData
{
public:
    //
    // Default constructor.
    //
    MCInterpBndryData () : BndryData() {}
    //
    // Constructor for given BoxArray, etc.
    //
    MCInterpBndryData (const BoxArray& _grids,
		       int             _ncomp,
		       const Geometry& geom);
    //
    // Interpret user BCs to those expected by LinOp.
    //
    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio,
                                int          comp=0) = 0;
    //
    // Set bndry values at coarse level (non interpolation performed).
    //
    void setBndryValues (const MultiFab&     mf,
			 int                 mf_start,
			 int                 bnd_start,
			 int                 num_comp,
			 const Array<BCRec>& phys_bc);
    //
    // Set bndry values at fine level, performing necessary interpolations.
    //
    void setBndryValues (const ::BndryRegister& crse,
			 int                  c_start,
			 const MultiFab&      fine,
			 int                  f_start,
			 int                  bnd_start,
			 int                  num_comp,
			 int                  ratio,
			 const Array<BCRec>&  phys_bc);
};

#endif /*_MCINTERPBNDRYDATA_H_*/
