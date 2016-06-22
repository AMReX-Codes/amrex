
#ifndef _INTERPBNDRYDATA_H_
#define _INTERPBNDRYDATA_H_

#include <BoxArray.H>
#include <FArrayBox.H>
#include <Orientation.H>
#include <Mask.H>
#include <BCRec.H>
#include <Geometry.H>
#include <BoundCond.H>
#include <BndryData.H>
#include <BL_CXX11.H>

/*
        An InterpBndryData object adds to a BndryData object the ability to
        manipulate and set the data stored in the boundary cells.

        The Interpbndrydata class is a virtual base class derived from
        BndryData.  It is intended to provide a more physical method for
        filling boundary-related data.  Boundary values in a BndryData object
        are stored in FabSets around each grid in the domain, and the
        InterpBndryData class provides a mechanism for filling these FabSets,
        consistent with AMR-like numerical discretizations.  When asked to
        set its boundary values, an InterpBndryData object:

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

class InterpBndryData
    :
    public BndryData
{
public:
    //
    // default constructor
    //
    InterpBndryData ();
    //
    // constructor for given BoxArray, etc
    //
    InterpBndryData (const BoxArray& _grids,
                     int             _ncomp,
                     const Geometry& geom,
		     ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //
    // Copy constructor.
    //
    InterpBndryData (const InterpBndryData& rhs);

    virtual ~InterpBndryData () BL_OVERRIDE;
    //
    // Copy assignment operator.
    //
    InterpBndryData& operator= (const InterpBndryData& rhs);
    //
    // interpret user BCs to those expected by LinOp: takes integer ratio and
    // is a wrapper to the IntVect version
    //
    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio);
    //
    // interpret user BCs to those expected by LinOp: takes IntVect ratio
    //
    virtual void setBndryConds (const BCRec&   phys_bc,
                                const IntVect& ratio,
				int            comp = 0) = 0;
    //
    // set bndry values at coarse level (non interpolation performed)
    //
    virtual void setBndryValues (const MultiFab& mf,
                                 int             mf_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const BCRec&    phys_bc);
    //
    // set bndry values at fine level, performing necessary interpolations
    //
    virtual void setBndryValues (::BndryRegister&  crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 int             ratio, 
                                 const BCRec&    phys_bc,
                                 int             max_order = IBD_max_order_DEF);    
    //
    // set bndry values at fine level, performing necessary interpolations
    //
    virtual void setBndryValues (::BndryRegister&  crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const IntVect&  ratio, 
                                 const BCRec&    phys_bc,
                                 int             max_order = IBD_max_order_DEF);
    
    static int IBD_max_order_DEF;

    static int maxOrderDEF() {return IBD_max_order_DEF;}
};

#endif /*_INTERPBNDRYDATA_H_*/

