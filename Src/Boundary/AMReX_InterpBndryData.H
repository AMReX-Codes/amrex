
#ifndef AMREX_INTERPBNDRYDATA_H_
#define AMREX_INTERPBNDRYDATA_H_

#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Orientation.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoundCond.H>
#include <AMReX_BndryData.H>

namespace amrex {

/**
        An InterpBndryData object adds to a BndryData object the ability to
        manipulate and set the data stored in the boundary cells.

        The InterpBndryData class is a virtual base class derived from
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

    /**
    * \brief default constructor
    */
    InterpBndryData () noexcept;

    /**
    * \brief constructor for given BoxArray, etc
    *
    * \param _grids
    * \param _dmap
    * \param _ncomp
    * \param geom
    */
    InterpBndryData (const BoxArray& _grids,
		     const DistributionMapping& _dmap,
                     int             _ncomp,
                     const Geometry& geom);

    /**
    * \brief Copy constructor.
    *
    * \param rhs
    */
    InterpBndryData (const InterpBndryData& rhs);

    virtual ~InterpBndryData () override;
    //
    //! Copy assignment operator.
    InterpBndryData& operator= (const InterpBndryData& rhs);

    InterpBndryData (InterpBndryData&& rhs) = delete;
    InterpBndryData& operator= (InterpBndryData&& rhs) = delete;


    /**
    * \brief interpret user BCs to those expected by LinOp: takes integer ratio and
    * is a wrapper to the IntVect version
    *
    * \param phys_bc
    * \param ratio
    */
    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio);

    /**
    * \brief interpret user BCs to those expected by LinOp: takes IntVect ratio
    *
    * \param phys_bc
    * \param ratio
    * \param comp
    */
    virtual void setBndryConds (const BCRec&   phys_bc,
                                const IntVect& ratio,
				int            comp = 0) {}

    /**
    * \brief set bndry values at coarse level (non interpolation performed)
    *
    * \param mf
    * \param mf_start
    * \param bnd_start
    * \param num_comp
    * \param phys_bc
    */
    virtual void setBndryValues (const MultiFab& mf,
                                 int             mf_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const BCRec&    phys_bc);

    virtual void setBndryValues (const MultiFab& mf,
                                 int             mf_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 int             ratio,
                                 const BCRec&    phys_bc);

    virtual void setBndryValues (const MultiFab& mf,
                                 int             mf_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const IntVect&  ratio,
                                 const BCRec&    phys_bc);

    /**
    * \brief set bndry values at fine level, performing necessary interpolations
    *
    * \param crse
    * \param c_start
    * \param fine
    * \param f_start
    * \param bnd_start
    * \param num_comp
    * \param ratio
    * \param phys_bc
    * \param max_order
    */
    virtual void setBndryValues (BndryRegister&  crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 int             ratio,
                                 const BCRec&    phys_bc,
                                 int             max_order = IBD_max_order_DEF);

    /**
    * \brief set bndry values at fine level, performing necessary interpolations
    *
    * \param crse
    * \param c_start
    * \param fine
    * \param f_start
    * \param bnd_start
    * \param num_comp
    * \param ratio
    * \param phys_bc
    * \param max_order
    */
    virtual void setBndryValues (BndryRegister&  crse,
                                 int             c_start,
                                 const MultiFab& fine,
                                 int             f_start,
                                 int             bnd_start,
                                 int             num_comp,
                                 const IntVect&  ratio,
                                 const BCRec&    phys_bc,
                                 int             max_order = IBD_max_order_DEF);

    void updateBndryValues (BndryRegister& crse, int c_start, int bnd_start, int num_comp,
                            const IntVect& ratio, int max_order = IBD_max_order_DEF);

    void updateBndryValues (BndryRegister& crse, int c_start, int bnd_start, int num_comp,
                            int ratio, int max_order = IBD_max_order_DEF);

    void setHomogValues ();

    static int IBD_max_order_DEF;

    static int maxOrderDEF() {return IBD_max_order_DEF;}

// private:  // public for cuda

    void BndryValuesDoIt (BndryRegister&  crse,
                          int             c_start,
                          const MultiFab* fine,
                          int             f_start,
                          int             bnd_start,
                          int             num_comp,
                          const IntVect&  ratio,
                          const BCRec*    phys_bc,
                          int             max_order);
};

}

#endif /*_INTERPBNDRYDATA_H_*/
