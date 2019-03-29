
#ifndef AMREX_LINOP_H_
#define AMREX_LINOP_H_

#include <memory>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BndryData.H>

namespace amrex {

/**
        \brief A LinOp is a virtual base class for general linear operators capable
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

    /**
    * \brief Allocate a LinOp for this box array, boundary and (uniform) spacing info.
    *
    * \param mgb
    * \param _h
    */
    LinOp (const BndryData& mgb,
           const Real       _h);

    /**
    * \brief Allocate a LinOp for this box array, boundary data and spacing array info.
    *
    * \param mgb
    * \param _h
    */
    LinOp (const BndryData& mgb,
           const Real*      _h);

    /**
    * \brief Allocate a LinOp for this box array, boundary data and spacing array info.
    *
    * LinOp assumes ownership of the pointer.
    *
    * \param mgb
    * \param _h
    */
    LinOp (BndryData*  mgb,
           const Real* _h);

    /**
    * \brief The destructor.
    */
    virtual ~LinOp ();

    /**
    * \brief Applies level LinOp to "in", returns "out", uses BC_mode flag.
    *
    * \param out
    * \param in
    * \param level
    * \param bc_mode
    * \param local
    * \param src_comp
    * \param dst_comp
    * \param num_comp
    * \param bndry_comp
    */
    virtual void apply (MultiFab&      out,
                        MultiFab&      in,
                        int            level   = 0,
                        LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC,
                        bool           local   = false,
			            int            src_comp = 0,
                        int            dst_comp = 0,
                        int            num_comp = 1,
                        int            bndry_comp = 0);

    /**
    * \brief Fills level boundary cells using BC_mode flag, int. BC data if reqd.
    *
    * \param inout
    * \param src_comp
    * \param num_comp
    * \param level
    * \param bc_mode
    * \param local
    * \param bndryComp
    */
    virtual void applyBC (MultiFab&      inout,
                  int            src_comp,
                  int            num_comp,
                  int            level   = 0,
                  LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC,
                  bool           local   = false,
		          int            bndryComp = 0);

    /**
    * \brief Compute the level residual = rhsL - L(solnL).
    *
    * \param residL
    * \param rhsL
    * \param solnL
    * \param level
    * \param bc_mode
    * \param local
    */
    virtual void residual (MultiFab&       residL,
                           const MultiFab& rhsL,
                           MultiFab&       solnL,
                           int             level   = 0,
                           LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC,
                           bool            local   = false);

    /**
    * \brief Smooth the level system L(solnL)=rhsL.
    *
    * \param solnL
    * \param rhsL
    * \param level
    * \param bc_mode
    */
    virtual void smooth (MultiFab&       solnL,
                         const MultiFab& rhsL,
                         int             level   = 0,
                         LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);

    virtual void jacobi_smooth (MultiFab&       solnL,
                                const MultiFab& rhsL,
                                int             level   = 0,
                                LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);

    /**
    * \brief Estimate the norm of the operator.
    *
    * \param nm
    * \param level
    * \param local
    */
    virtual Real norm (int nm = 0, int level = 0, const bool local = false);

    /**
    * \brief Compute flux associated with the op
    *
    * \paramxflux
    * \paramyflux
    * \paramzflux)
    * \param in
    * \param bc_mode
    * \param sComp
    * \param dComp
    * \param nComp
    * \param bndComp
    */
    virtual void compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=Inhomogeneous_BC,
			   int sComp=0, int dComp=0, int nComp=1, int bndComp=0) = 0;

    /**
    * \brief Return the boundary data object.
    */
    const BndryData& bndryData () const { return *bgb; }

    const DistributionMapping& DistributionMap () const { return bgb->DistributionMap(); }

    /**
    * \brief Set the boundary data object.
    *
    * \param bd
    */
    void bndryData (const BndryData& bd);

    /**
    * \brief Return the box array.
    *
    * \param level
    */
    virtual const BoxArray& boxArray (int level = 0) const
    {
        BL_ASSERT(level < numLevels()); return gbox[level];
    }

    /**
    * \brief Return the number of grids.
    */
    int numGrids () const { return gbox[0].size(); }

    /**
    * \brief Return the number of levels.
    */
    virtual int numLevels () const { return h.size(); }

    /**
    * \brief Return the order of the boundary condition interpolant.
    */
    virtual int maxOrder () const { return maxorder; }

    /**
    * \brief Set the order of the boundary condition interpolant.
    *
    * \param maxorder_
    */
    virtual int maxOrder (int maxorder_);

    /**
    * \brief Return the number of grow cells this operator expects in the input state to compute "apply"
    *
    * \param level
    */
    virtual int NumGrow (int level = 0) const {return LinOp_grow;}

    /**
    * \brief Construct/allocate internal data necessary for adding a new level.
    *
    * \param level
    */
    virtual void prepareForLevel (int level);

    /**
    * \brief Output operator internal to an ASCII stream.
    *
    * \param os
    * \param lp
    */
    friend std::ostream& operator<< (std::ostream& os, const LinOp& lp);

    const Geometry& getGeom(int level);
    const Real * getDx(int level);

    /**
    * \brief Get scalar alpha coefficient.
    */
    virtual Real get_alpha () const;

    /**
    * \brief Get scalar beta coefficient.
    */
    virtual Real get_beta () const;

    /**
    * \brief Return reference to "a" coefficients for base level.
    *
    * \param level
    */
    virtual const MultiFab& aCoefficients (int level = 0);

    /**
    * \brief Return reference to "b" coefficients for base level.
    *
    * \param dir
    * \param level
    */
    virtual const MultiFab& bCoefficients (int dir, int level=0);

protected:

    /**
    * \brief Remove internal data necessary for a level and all higher.
    *
    * \param level
    */
    virtual void clearToLevel (int level);

    /**
    * \brief Virtual to apply the level operator to the internal nodes of
    * "in", return result in "out"
    *
    * \param out
    * \param in
    * \param level
    */
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level) = 0;
    virtual void Fapply (MultiFab&       out,
			 int             dst_comp,
                         const MultiFab& in,
			 int             src_comp,
			 int             num_comp,
                         int             level) = 0;

    /**
    * \brief Virtual to carry out the level smoothing operation for
    * L(solnL)=rhsL on internal nodes.  Modify solnL in place.
    *
    * \param solnL
    * \param rhsL
    * \param level
    * \param rgbflag
    */
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag) = 0;

    virtual void Fsmooth_jacobi (MultiFab&       solnL,
                                 const MultiFab& rhsL,
                                 int             level) = 0;

    /**
    * \brief Build coefficients at coarser level by interpolating "fine"
    * (builds in appropriate node/cell centering)
    *
    * \param crs
    * \param fine
    * \param level
    */
    void makeCoefficients (MultiFab&       crs,
                           const MultiFab& fine,
                           int             level);

    /**
    * \brief Initialize LinOp internal data.
    */
    static void Initialize ();

    static void Finalize ();

    /**
    * \brief Helper function for object construction.
    *
    * \param _h
    */
    void initConstruct (const Real* _h);

    //! Vector (on level) of Arrays (on dimension) of grid spacings
    Vector< Array<Real,BL_SPACEDIM> > h;

    //! Vector (on level) of BoxArray's of LinOp's domain
    Vector< BoxArray > gbox;

    /**
    * \brief Vector (on level) of pointers to BndryRegisters along each grid
    * for scratch data required to modify internal stencil on boundary
    */
    Vector<BndryRegister> undrrelxr;

    /**
    * \brief Vector (on level) of Vectors (on orientation) of MultiMasks for
    * whether boundary Fabs are covered, not_covered, outside_domain
    */
    Vector<Vector<MultiMask> > maskvals;

    /**
    * \brief Similar to maskvals, but only for level 0.
    * It's only made if it's needed.  Currently used only in applyBC().
    */
    Vector<Vector<MultiMask> > lmaskvals;

    //! boundary data class
    BndryData* bgb;

    /**
    * \brief Vector (on level) of geometry objects.  Needed for determining
    * whether stuff intersects on periodic domains
    */
    Vector< Geometry > geomarray;

    /**
    * \brief flag (=1 if use harmonic averaged interpolation for coefficients,
    * =0 is arithmetic averaging)
    */
    int harmavg;

    //! flag (\>0 is verbose execution)
    int verbose;

    /**
    * \brief maximum interpolation order used for constructing Dirichlet
    * ghost node values
    */
    int maxorder;

    //! default value for harm_avg
    static int def_harmavg;

    //! default value for verbose
    static int def_verbose;

    //! default maximum BC interpolant order
    static int def_maxorder;

   //! Number of grow cells required for this operator
   static int LinOp_grow;

private:

    /**
    * \brief Not implemented.
    *
    * \param rhs
    */
    LinOp (const LinOp& rhs);
    void operator = (const LinOp& rhs);
};

}

#endif /*_LINOP_H_*/
