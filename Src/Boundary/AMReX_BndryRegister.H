
#ifndef AMREX_BNDRYREGISTER_H_
#define AMREX_BNDRYREGISTER_H_

#include <AMReX_BoxArray.H>
#include <AMReX_FabSet.H>

namespace amrex {

class Orientation;


// \cond GENCODE
struct BndryBATransformer
    : public BATBase<BndryBATransformer>
{
    BndryBATransformer (Orientation face, IndexType typ,
			int in_rad, int out_rad, int extent_rad);

    virtual ~BndryBATransformer () override {}

    virtual IntVect doiLo () const override { return m_doilo; }
    virtual IntVect doiHi () const override { return m_doihi; }

    virtual Box operator() (const Box& bx) const override;

    bool operator== (const BndryBATransformer& rhs) const noexcept;

private:
    int     m_dir;
    bool    m_lo_face;
    int     m_in_rad;
    int     m_out_rad;
    int     m_extent_rad;
    // The variables below are derived from the variables above.
    IntVect m_nodal;
    IntVect m_loshft;
    IntVect m_hishft;
    IntVect m_doilo;
    IntVect m_doihi;
};
// \endcond

/**
* \brief A BndryRegister organizes FabSets bounding each grid in a BoxArray.
*        A FabSet is maintained for each boundary orientation, as well as
*        the BoxArray domain of definition.
*
*        A BndryRegister object contains a list of FabSets bounding the grids
*        in a BoxArray.  The FabSet FABs are at the same refinement level
*        as the grids they bound, and are accessed and modified via a variety
*        of member functions.
*
*        Non-default instantiation allocates a set of FABs, grown into and
*        out of the bounding surfaces of each box in the BoxArray.  The width of
*        the layer (in/out), as well as the "extent" of a bounding FABs (the
*        number of nodes beyond the box boundaries, parallel to the grid
*        surface) are determined by constructor argument.  All boxes and
*        FABs in this context are assumed to be cell-centered.
*
*        A small number of linear mathematical operations are provided for
*        BndryRegisters, as well as overloaded [] operators for access based
*        on grid boundary orientation.  The BoxArray domain of definition is
*        accessible, but not resettable,
*/
class BndryRegister
{

public:
    //! The default constructor.
    BndryRegister () noexcept;

    //! The constructor, given number of cells in/out, extent and number of components (assumes cell-centered boxes, and allocates cell-centered FABs)
    BndryRegister (const BoxArray& grids_,
		   const DistributionMapping& dmap,
                   int             in_rad,
                   int             out_rad,
                   int             extent_rad,
                   int             ncomp);
    //
    //! The copy constructor.
    BndryRegister (const BndryRegister& src);

    //! The copy assignment operator.
    BndryRegister& operator= (const BndryRegister& src);

    //! The destructor.
    virtual ~BndryRegister();

    BndryRegister (BndryRegister&& rhs) noexcept = default;

    BndryRegister& operator= (BndryRegister&& rhs) = delete;

    void define (const BoxArray& grids_,
		 const DistributionMapping& dmap,
		 int             in_rad,
		 int             out_rad,
		 int             extent_rad,
		 int             ncomp);
    //! Build FABs along given face, specifying the DistributionMapping.
    void define (Orientation                face,
		 IndexType                  typ,
		 int                        in_rad,
		 int                        out_rad,
		 int                        extent_rad,
		 int                        ncomp,
		 const DistributionMapping& dm);

    void clear ();

    //! Get box domain (as an array of boxes).
    const BoxArray& boxes () const noexcept { return grids; }

    //! Return the number of grids in this domain.
    int size () const noexcept { return grids.size(); }

    //! Return const set of FABs bounding the domain grid boxes on a given orientation
    const FabSet& operator[] (Orientation face) const noexcept { return bndry[face]; }

    //! Return set of FABs bounding the domain grid boxes on a given orientation
    FabSet& operator[] (Orientation face) noexcept { return bndry[face]; }

    //! Set all boundary FABs to given value.
    void setVal (Real v);

    //! register += rhs
    BndryRegister& operator+= (const BndryRegister& rhs);
    BndryRegister& plus (const BndryRegister& rhs);

    //! Fill the boundary FABs on intersection with given MultiFab.
    BndryRegister& copyFrom (const MultiFab& src,
                             int             nghost,
                             int             src_comp,
                             int             dest_comp,
                             int             num_comp,
                             const Periodicity& period = Periodicity::NonPeriodic());

    //! Increment the boundary FABs on intersect with given MultiFab.
    BndryRegister& plusFrom (const MultiFab& src,
                             int             nghost,
                             int             src_comp,
                             int             dest_comp,
                             int             num_comp,
                             const Periodicity& period = Periodicity::NonPeriodic());

    //! Linear combination: this := a*mfa + b*mfb on intersection of MultiFabs with the boundary FABs
    BndryRegister& linComb (Real            a,
                            const MultiFab& mfa,
                            int             a_comp,
                            Real            b,
                            const MultiFab& mfb,
                            int             b_comp,
                            int             dest_comp,
                            int             num_comp,
                            int             n_ghost = 0);

    //! Set box domain, if not set previously.
    void setBoxes (const BoxArray& grids);

    //! Returns constant reference to associated DistributionMapping.
    const DistributionMapping& DistributionMap () const noexcept { return bndry[0].DistributionMap(); }

    //! Write (used for writing to checkpoint)
    void write (const std::string& name, std::ostream& os) const;

    //! Read (used for reading from checkpoint)
    void read (const std::string& name, std::istream& is);

    //! Local copy function
    static void Copy (BndryRegister& dst, const BndryRegister& src);

    //! Used by the copy constructor and assignment operator.
    void init (const BndryRegister& src);

protected:

    //! The data.
    FabSet    bndry[2*AMREX_SPACEDIM];
    BoxArray  grids;
};

}

#endif /*_BNDRYREGISTER_H_*/
