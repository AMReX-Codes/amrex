
#ifndef AMREX_MCINTERPBNDRYDATA_H_
#define AMREX_MCINTERPBNDRYDATA_H_

#include <AMReX_BLassert.H>
#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Orientation.H>
#include <AMReX_Mask.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoundCond.H>
#include <AMReX_BndryData.H>

namespace amrex {

/**
        \brief An MCInterpBndryData object adds to a BndryData object the ability to
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

    /**
    * \brief Default constructor.
    */
    MCInterpBndryData () : BndryData() {}

    /**
    * \brief Constructor for given BoxArray, etc.
    *
    * \param _grids
    * \param _dmap
    * \param _ncomp
    * \param geom
    */
    MCInterpBndryData (const BoxArray& _grids,
		       const DistributionMapping& _dmap,
		       int             _ncomp,
		       const Geometry& geom);

    virtual ~MCInterpBndryData () = default;

    MCInterpBndryData (MCInterpBndryData&& rhs) = delete;
    MCInterpBndryData (const MCInterpBndryData& rhs) = delete;
    MCInterpBndryData& operator= (const MCInterpBndryData& rhs) = delete;
    MCInterpBndryData& operator= (MCInterpBndryData&& rhs) = delete;


    /**
    * \brief Interpret user BCs to those expected by LinOp.
    *
    * \param phys_bc
    * \param ratio
    * \param comp
    */
    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio,
                                int          comp=0) = 0;

    /**
    * \brief Set bndry values at coarse level (non interpolation performed).
    *
    * \param mf
    * \param mf_start
    * \param bnd_start
    * \param num_comp
    * \param phys_bc
    */
    void setBndryValues (const MultiFab&     mf,
			 int                 mf_start,
			 int                 bnd_start,
			 int                 num_comp,
			 const Vector<BCRec>& phys_bc);

    /**
    * \brief Set bndry values at fine level, performing necessary interpolations.
    *
    * \param crse
    * \param c_start
    * \param fine
    * \param f_start
    * \param bnd_start
    * \param num_comp
    * \param ratio
    * \param phys_bc
    */
    void setBndryValues (const BndryRegister& crse,
			 int                  c_start,
			 const MultiFab&      fine,
			 int                  f_start,
			 int                  bnd_start,
			 int                  num_comp,
			 int                  ratio,
			 const Vector<BCRec>&  phys_bc);
};

}

#endif /*_MCINTERPBNDRYDATA_H_*/
