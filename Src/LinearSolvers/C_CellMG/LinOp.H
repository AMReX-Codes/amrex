//BL_COPYRIGHT_NOTICE

#ifndef _LINOP_H_
#define _LINOP_H_

//
// $Id: LinOp.H,v 1.12 2000-08-24 16:02:46 car Exp $
//

#include <Array.H>
#include <Tuple.H>
#include <Pointers.H>
#include <REAL.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <BndryData.H>

//@Man:
/*@Memo:
  A LinOp is a virtual base class for general linear operators capable
  of acting on MultiFabs.  All implementation and access functions are
  designed to make a LinOp object useful for representing and solving
  a set of linear equations on a union of rectangles in 2D or 3D.
*/        
/*@Doc:
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
        LO\_DIRICHLET", the ghost cells not covered by adjacent grids are
        filled on demand by a polynomial interpolant (of settable order) to
        compute the value specified in the BndryData FabSets along the edge
        of the grid boxes.  A value is put into the ghost cell of the input
        multifab that can be assumed to exist outside the domain, at equal
        spacing with the rest of the nodes.  
        "BndryData::LO\_NEUMANN" conditions result in ghost cells filled by
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
        apply() member take a BC\_Mode flag.  The BC\_mode can be either
        Homogeneous\_BC, or Inhomogeneous\_BC.  It is a strict requirement of
        the linear operator that LinOp::apply(out,in,level,bc\_mode=Homogeneous\_BC)
        acting on in=0 returns out=0.
        
        This class does NOT provide a copy constructor or assignment operator.
*/

class LinOp
{
public:

    enum BC_Mode { Homogeneous_BC = 0, Inhomogeneous_BC };
    //
    //@ManDoc: allocate a LinOp for this box array, boundary and (uniform) spacing info
    //
    LinOp (const BndryData& mgb,
           const Real       _h);
    //
    //@ManDoc: allocate a LinOp for this box array, boundary data and spacing array info
    //
    LinOp (const BndryData& mgb,
           const Real*      _h);
    //
    //@ManDoc: alocate a LinOp w/level structure beginning at the
    //  specified level of this LinOp
    //
    LinOp (const LinOp& _Lp,
           int          level);
    //
    //@ManDoc: destructor
    //
    virtual ~LinOp ();
    //
    //@ManDoc: applies level LinOp to "in", returns "out", uses BC\_mode flag
    //
    virtual void apply (MultiFab&      out,
                        MultiFab&      in,
                        int            level = 0,
                        LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC);
    //
    //@ManDoc: fills level boundary cells using BC\_mode flag, int. BC data if reqd
    //
    void applyBC (MultiFab&      inout,
                  int            src_comp,
                  int            num_comp,
                  int            level = 0,
                  LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC);
    //
    //@ManDoc: compute the level residual = rhsL - L(solnL)
    //
    virtual void residual (MultiFab&       residL,
                           const MultiFab& rhsL,
                           MultiFab&       solnL,
                           int             level = 0,
                           LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);
    //
    //@ManDoc: smooth the level system L(solnL)=rhsL
    //
    void smooth (MultiFab&       solnL,
                 const MultiFab& rhsL,
                 int             level = 0,
                 LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);
#if 0
    //
    //@ManDoc: Compute the norm of "in"
    //
    Real norm (const MultiFab& in,
               int             level = 0) const;
#endif
    //
    //@ManDoc: Compute flux associated with the op
    //
    virtual void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=Inhomogeneous_BC) = 0;
    //
    //@ManDoc: return the boundary data object
    //
    const BndryData& bndryData () const;
    //
    //@ManDoc: set the boundary data object
    //
    void bndryData (const BndryData& bd);
    //
    //@ManDoc: return the box array
    //
    const BoxArray& boxArray (int level = 0) const;
    //
    //@ManDoc: return the number of grids
    //
    int numGrids () const;
    //
    //@ManDoc: return the number of levels
    //
    int numLevels () const;
    //
    //@ManDoc: return the order of the boundary condition interpolant
    //
    int maxOrder () const;
    //
    //@ManDoc: set the order of the boundary condition interpolant
    //
    int maxOrder (int maxorder_);
    //
    //@ManDoc: construct/allocate internal data necessary for adding a new level
    //
    virtual void prepareForLevel (int level);
    //
    //@ManDoc: Output operator internal to an ASCII stream
    //
    friend ostream& operator<< (ostream&     os,
                                const LinOp& lp);
    
protected:
    //
    //@ManDoc: Remove internal data necessary for a level and all higher.
    //
    virtual void clearToLevel (int level) {}
    //
    //@ManDoc: Virtual to apply the level operator to the internal nodes of
    // "in", return result in "out"
    //
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level) = 0;
    //
    //@ManDoc: Virtual to carry out the level smoothing operation for
    //  L(solnL)=rhsL on internal nodes.  Modify solnL in place.
    //
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag) = 0;
    //
    //@ManDoc: build coefficients at coarser level by interpolating "fine"
    //  (builds in appropriate node/cell centering)
    //
    void makeCoefficients (MultiFab&       crs,
                           const MultiFab& fine,
                           int             level);
    //
    //@ManDoc: Initialize LinOp internal data.
    //
    static void initialize ();
    //
    //@ManDoc: Helper function for object construction.
    //
    void initConstruct (const Real* _h);
    //
    //@ManDoc: Array (on level) of Tuples (on dimension) of grid spacings
    //
    vector< Tuple< Real, BL_SPACEDIM > > h;
    //
    //@ManDoc: Array (on level) of BoxArray's of LinOp's domain
    //
    vector< BoxArray > gbox;
    //
    //@ManDoc: Array (on level) of pointers to BndryRegisters along each grid
    //  for scratch data required to modify internal stencil on boundary
    //
    vector< CpClassPtr<BndryRegister> > undrrelxr;
    //
    //@ManDoc: Array (on level) of Arrays (on grid) of Arrays (on orientation)
    // of pointers to Masks for whether boundary Fabs are covered,
    // not\_covered, outside\_domain
    //
    //  FIXME: This should really be an array of integer BndryRegisters
    //
    Array< Array< Array< Mask*> > > maskvals;
    //
    //@ManDoc: boundary data class
    //
    BndryData bgb;
    //
    //@ManDoc: Array (on level) of geometry objects.  Needed for determining
    // whether stuff intersects on periodic domains
    //
    vector< Geometry > geomarray;
    //
    //@ManDoc: flag (=1 if use harmonic averaged interpolation for coefficients,
    // =0 is arithmetic averaging)
    //
    int harmavg;
    //
    //@ManDoc: flag (>0 is verbose execution)
    //
    int verbose;
    //
    //@ManDoc: maximum interpolation order used for constructing Dirichlet
    // ghost node values
    //
    int maxorder;
    //
    //@ManDoc: flag (=1 if internal data initialized)
    //
    static bool initialized;
    //
    //@ManDoc: default value for harm\_avg
    //
    static int def_harmavg;
    //
    //@ManDoc: default value for verbose
    //
    static int def_verbose;
    //
    //@ManDoc: default maximum BC interpolant order
    //
    static int def_maxorder;
    
private:
    //
    // Not implemented.
    //
    LinOp (const LinOp& rhs);
    void operator = (const LinOp& rhs);
};

inline
const BndryData&
LinOp::bndryData () const
{
    return bgb;
}

inline
void
LinOp::bndryData (const BndryData& bd)
{
    BL_ASSERT(gbox[0] == bd.boxes());
    bgb = bd;
}

inline
int
LinOp::numLevels () const
{
    return h.size();
}

inline
const BoxArray&
LinOp::boxArray (int level) const
{
    BL_ASSERT(level < numLevels());
    return gbox[level];
}

inline
int
LinOp::numGrids () const
{
    return gbox[0].length();
}

inline
int
LinOp::maxOrder () const
{
    return maxorder;
}

inline
int
LinOp::maxOrder (int maxorder_)
{
    BL_ASSERT(maxorder_ >= 2);
    maxorder_ = (maxorder_ < 2 ? 2 : maxorder_ );
    int omaxorder = maxorder;
    maxorder = maxorder_;
    return omaxorder;
}

#endif /*_LINOP_H_*/
