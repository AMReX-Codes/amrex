
#ifndef _LINOP_H_
#define _LINOP_H_

#include <memory>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Tuple.H>
#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BndryData.H>

namespace amrex {

/*
        A LinOp is a virtual base class for general linear operators capable
        of acting on MultiFabs.  All implementation and access functions are
        designed to make a LinOp object useful for representing and solving
        a set of linear equations on a union of rectangles in 2D or 3D.

        A LinOp contains all the information necessary to construct and
        maintain a linear operator which acts on a cell-centered MultiFab.  It
        is a virtual base class, and therefore cannot itself be
        instantiated.  The virtual functions implement "apply" 
         and "smooth" operators (for the internal nodes), and therefore
        determine the mathematical character of the operator.  LinOp,
        on the other hand, determines the implementation,
        filling ghost cells and creating coarsened versions
        of the domain and operator, etc. necessary for solving linear systems.

        LinOp is designed specifically for representing discrete approximations
        to partial differential operators.  It requires a BndryData object,
        which maintains ghost cell data surrounding each rectangle in the
        domain, as well as position and boundary condition type
        specifiers.  It also requires a (vector/scalar) grid spacing.  On
        levels above the base level, LinOp internally recognizes adjacent
        domain rectangles (using methods similar to those in the BndryData
        class) and when applying boundary conditions (filling ghost cell
        regions) will give preference to "valid" data (i.e. ghost cells
        covered by adjacent grids are filled with the corresponding data from
        the interior of that adjacent grid rather than via the auxiliary
        boundary instructions).

        A LinOp constructs a set of "levels", which are useful for linear
        solution methods such as multigrid.  On each grid, a new level is
        created by uniformly coarsening the grid structure by a factor of
        two in each coordinate direction (and then allocating and initializing
        any internal data necessary--new level grid spacing, for example).
        A LinOp can fill boundary ghost cells, compute a "norm" and coordinate
        the "apply" and "smooth"  operations at each level.
        Note that there are the same number of levels on each grid in the
        LinOp's domain.

        If a boundary type specifier indicates the type "BndryData::
        LO_DIRICHLET", the ghost cells not covered by adjacent grids are
        filled on demand by a polynomial interpolant (of settable order) to
        compute the value specified in the BndryData FabSets along the edge
        of the grid boxes.  A value is put into the ghost cell of the input
        multifab that can be assumed to exist outside the domain, at equal
        spacing with the rest of the nodes.  
        "BndryData::LO_NEUMANN" conditions result in ghost cells filled by
        second order extrapolation, and results in again filling the
        ghost cell with a value appropriate for application in a constant
        spacing cell-centered discretization.
        The order of the interpolant is set
        by default in this class to 2, which is compatible with most
        simple algebraic relaxation methods for linear systems based on
        numerical PDE's.  The interpolant can be queried/set via member
        functions, and is applied uniformly over the domain.  The boundary
        location is specified as a distance from the wall of the grid block,
        and is given in the same units that the grid spacing, h, is specified.

        All the member functions of LinOp which interact directly with the
        apply() member take a BC_Mode flag.  The BC_mode can be either
        Homogeneous_BC, or Inhomogeneous_BC.  It is a strict requirement of
        the linear operator that LinOp::apply(out,in,level,bc_mode=Homogeneous_BC)
        acting on in=0 returns out=0.
        
        This class does NOT provide a copy constructor or assignment operator.
*/

class LinOp
{
public:

    enum BC_Mode { Homogeneous_BC = 0, Inhomogeneous_BC };
    //
    // Allocate a LinOp for this box array, boundary and (uniform) spacing info.
    //
    LinOp (const BndryData& mgb,
           const Real       _h);
    //
    // Allocate a LinOp for this box array, boundary data and spacing array info.
    //
    LinOp (const BndryData& mgb,
           const Real*      _h);
    //
    // Allocate a LinOp for this box array, boundary data and spacing array info.
    //
    // LinOp assumes ownership of the pointer.
    //
    LinOp (BndryData*  mgb,
           const Real* _h);
    //
    // The destructor.
    //
    virtual ~LinOp ();
    //
    // Applies level LinOp to "in", returns "out", uses BC_mode flag.
    //
    virtual void apply (MultiFab&      out,
                        MultiFab&      in,
                        int            level   = 0,
                        LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC,
                        bool           local   = false,
			int            src_comp = 0,
                        int            dst_comp = 0,
                        int            num_comp = 1,
                        int            bndry_comp = 0);
    //
    // Fills level boundary cells using BC_mode flag, int. BC data if reqd.
    //
    virtual void applyBC (MultiFab&      inout,
                  int            src_comp,
                  int            num_comp,
                  int            level   = 0,
                  LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC,
                  bool           local   = false,
		  int            bndryComp = 0);
    //
    // Compute the level residual = rhsL - L(solnL).
    //
    virtual void residual (MultiFab&       residL,
                           const MultiFab& rhsL,
                           MultiFab&       solnL,
                           int             level   = 0,
                           LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC,
                           bool            local   = false);
    //
    // Smooth the level system L(solnL)=rhsL.
    //
    virtual void smooth (MultiFab&       solnL,
                         const MultiFab& rhsL,
                         int             level   = 0,
                         LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);

    virtual void jacobi_smooth (MultiFab&       solnL,
                                const MultiFab& rhsL,
                                int             level   = 0,
                                LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);
    //
    // Estimate the norm of the operator.
    //
    virtual Real norm (int nm = 0, int level = 0, const bool local = false);
    //
    // Compute flux associated with the op
    //
    virtual void compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=Inhomogeneous_BC,
			   int sComp=0, int dComp=0, int nComp=1, int bndComp=0) = 0;
    //
    // Return the boundary data object.
    //
    const BndryData& bndryData () const { return *bgb; }
    //
    ParallelDescriptor::Color color () const { return bgb->color(); }
    //
    const DistributionMapping& DistributionMap () const { return bgb->DistributionMap(); }
    //
    // Set the boundary data object.
    //
    void bndryData (const BndryData& bd);
    //
    // Return the box array.
    //
    virtual const BoxArray& boxArray (int level = 0) const
    {
        BL_ASSERT(level < numLevels()); return gbox[level];
    }
    //
    // Return the number of grids.
    //
    int numGrids () const { return gbox[0].size(); }
    //
    // Return the number of levels.
    //
    virtual int numLevels () const { return h.size(); }
    //
    // Return the order of the boundary condition interpolant.
    //
    virtual int maxOrder () const { return maxorder; }
    //
    // Set the order of the boundary condition interpolant.
    //
    virtual int maxOrder (int maxorder_);
    //
    // Return the number of grow cells this operator expects in the input state to compute "apply"
    //
    virtual int NumGrow (int level = 0) const {return LinOp_grow;}
    //
    // Construct/allocate internal data necessary for adding a new level.
    //
    virtual void prepareForLevel (int level);
    //
    // Output operator internal to an ASCII stream.
    //
    friend std::ostream& operator<< (std::ostream& os, const LinOp& lp);

    const Geometry& getGeom(int level);
    const Real * getDx(int level);
    //
    // Get scalar alpha coefficient.
    //
    virtual Real get_alpha () const;
    //
    // Get scalar beta coefficient.
    //
    virtual Real get_beta () const;
    //
    // Return reference to "a" coefficients for base level.
    //
    virtual const MultiFab& aCoefficients (int level = 0);
    //
    // Return reference to "b" coefficients for base level.
    //
    virtual const MultiFab& bCoefficients (int dir, int level=0);
    
protected:
    //
    // Remove internal data necessary for a level and all higher.
    //
    virtual void clearToLevel (int level);
    //
    // Virtual to apply the level operator to the internal nodes of
    // "in", return result in "out"
    //
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level) = 0;
    virtual void Fapply (MultiFab&       out,
			 int             dst_comp,
                         const MultiFab& in,
			 int             src_comp,
			 int             num_comp,
                         int             level) = 0;
    //
    // Virtual to carry out the level smoothing operation for
    //  L(solnL)=rhsL on internal nodes.  Modify solnL in place.
    //
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag) = 0;

    virtual void Fsmooth_jacobi (MultiFab&       solnL,
                                 const MultiFab& rhsL,
                                 int             level) = 0;
    //
    // Build coefficients at coarser level by interpolating "fine"
    //  (builds in appropriate node/cell centering)
    //
    void makeCoefficients (MultiFab&       crs,
                           const MultiFab& fine,
                           int             level);
    //
    // Initialize LinOp internal data.
    //
    static void Initialize ();

    static void Finalize ();
    //
    // Helper function for object construction.
    //
    void initConstruct (const Real* _h);
    //
    // Array (on level) of Tuples (on dimension) of grid spacings
    //
    Array< Tuple<Real,BL_SPACEDIM> > h;
    //
    // Array (on level) of BoxArray's of LinOp's domain
    //
    Array< BoxArray > gbox;
    //
    // Array (on level) of pointers to BndryRegisters along each grid
    //  for scratch data required to modify internal stencil on boundary
    //
    Array<BndryRegister> undrrelxr;
    //
    // Array (on level) of Arrays (on orientation) of MultiMasks for 
    // whether boundary Fabs are covered, not_covered, outside_domain
    //
    Array<Array<MultiMask> > maskvals;
    //
    // Similar to maskvals, but only for level 0.
    // It's only made if it's needed.  Currently used only in applyBC().
    //
    Array<Array<MultiMask> > lmaskvals;
    //
    // boundary data class
    //
    BndryData* bgb;
    //
    // Array (on level) of geometry objects.  Needed for determining
    // whether stuff intersects on periodic domains
    //
    Array< Geometry > geomarray;
    //
    // flag (=1 if use harmonic averaged interpolation for coefficients,
    // =0 is arithmetic averaging)
    //
    int harmavg;
    //
    // flag (>0 is verbose execution)
    //
    int verbose;
    //
    // maximum interpolation order used for constructing Dirichlet
    // ghost node values
    //
    int maxorder;
    //
    // default value for harm_avg
    //
    static int def_harmavg;
    //
    // default value for verbose
    //
    static int def_verbose;
    //
    // default maximum BC interpolant order
    //
    static int def_maxorder;
    //
    // Number of grow cells required for this operator
    //
   static int LinOp_grow;

private:
    //
    // Not implemented.
    //
    LinOp (const LinOp& rhs);
    void operator = (const LinOp& rhs);
};

}

#endif /*_LINOP_H_*/
