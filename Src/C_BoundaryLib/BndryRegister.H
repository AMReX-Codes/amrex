
#ifndef _BNDRYREGISTER_H_
#define _BNDRYREGISTER_H_

#include <BoxArray.H>
#include <FabSet.H>

class Orientation;

/*
        A BndryRegister organizes FabSets bounding each grid in a BoxArray.
        A FabSet is maintained for each boundary orientation, as well as
        the BoxArray domain of definition.

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

struct BndryBATransformer
    : public BATBase<BndryBATransformer>
{
    BndryBATransformer (Orientation face, IndexType typ,
			int in_rad, int out_rad, int extent_rad);

    virtual ~BndryBATransformer () BL_OVERRIDE {}

    virtual IntVect doiLo () const BL_OVERRIDE { return m_doilo; }
    virtual IntVect doiHi () const BL_OVERRIDE { return m_doihi; }

    virtual Box operator() (const Box& bx) const BL_OVERRIDE;
    
    bool operator== (const BndryBATransformer& rhs) const;

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

class BndryRegister
{

public:
    //
    // The default constructor.
    //
    BndryRegister ();
    //
    // The constructor, given number of cells in/out, extent and number of components (assumes cell-centered boxes, and allocates cell-centered FABs)
    //
    BndryRegister (const BoxArray& grids,
                   int             in_rad,
                   int             out_rad,
                   int             extent_rad,
                   int             ncomp,
		   ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //
    // The copy constructor.
    //
    BndryRegister (const BndryRegister& src);
    //
    // The copy assignment operator.
    //
    BndryRegister& operator= (const BndryRegister& src);
    //
    // The destructor.
    //
    virtual ~BndryRegister();
    //
    // Get box domain (as an array of boxes).
    //
    const BoxArray& boxes () const { return grids; }
    //
    // Return the number of grids in this domain.
    //
    int size () const { return grids.size(); }
    //
    // Return color
    //
    ParallelDescriptor::Color color() const { return bndry[0].color(); }
    //
    // Return const set of FABs bounding the domain grid boxes on a given orientation
    //
    const FabSet& operator[] (Orientation face) const { return bndry[face]; }
    //
    // Return set of FABs bounding the domain grid boxes on a given orientation
    //
    FabSet& operator[] (Orientation face) { return bndry[face]; }
    //
    // Set all boundary FABs to given value.
    //
    void setVal (Real v);
    //
    // register += rhs
    //
    BndryRegister& operator+= (const BndryRegister& rhs);
    BndryRegister& plus (const BndryRegister& rhs);
    //
    // Fill the boundary FABs on intersection with given MultiFab.
    //
    BndryRegister& copyFrom (const MultiFab& src,
                             int             nghost,
                             int             src_comp,
                             int             dest_comp,
                             int             num_comp);
    //
    // Increment the boundary FABs on intersect with given MultiFab.
    //
    BndryRegister& plusFrom (const MultiFab& src,
                             int             nghost,
                             int             src_comp,
                             int             dest_comp,
                             int             num_comp);
    //
    // Linear combination: this := a*mfa + b*mfb on intersection of MultiFabs with the boundary FABs
    //
    BndryRegister& linComb (Real            a,
                            const MultiFab& mfa,
                            int             a_comp,
                            Real            b,
                            const MultiFab& mfb,
                            int             b_comp,
                            int             dest_comp,
                            int             num_comp,
                            int             n_ghost = 0);
    //
    // Set box domain, if not set previously.
    //
    void setBoxes (const BoxArray& grids);
    //
    // Returns constant reference to associated DistributionMapping.
    //
    const DistributionMapping& DistributionMap () const { return bndry[0].DistributionMap(); }
    //
    // Build FABs along given face.
    //
    void define (Orientation face,
                 IndexType   typ,
                 int         in_rad,
                 int         out_rad,
                 int         extent_rad,
                 int         ncomp,
		 ParallelDescriptor::Color color = ParallelDescriptor::DefaultColor());
    //
    // Build FABs along given face, specifying the DistributionMapping.
    //
    void define (Orientation               face,
                 IndexType                 typ,
                 int                        in_rad,
                 int                        out_rad,
                 int                        extent_rad,
                 int                        ncomp,
                 const DistributionMapping& dm);

    //
    // Write (used for writing to checkpoint)
    //
    void write (const std::string& name, std::ostream& os) const;
    //
    // Read (used for reading from checkpoint)
    //
    void read (const std::string& name, std::istream& is);

    //
    // AddProcsToComp
    //
    virtual void AddProcsToComp(int ioProcNumSCS, int ioProcNumAll,
                                int scsMyId, MPI_Comm scsComm);

protected:
    //
    // Used by the copy constructor and assignment operator.
    //
    void init (const BndryRegister& src);
    //
    // The data.
    //
    FabSet    bndry[2*BL_SPACEDIM];
    BoxArray  grids;
};

#endif /*_BNDRYREGISTER_H_*/
