#ifndef _MCLINOP_H_
#define _MCLINOP_H_ "$Id: MCLinOp.H,v 1.1 1998-03-26 19:15:35 marc Exp $"

#include <Assert.H>
#include <Array.H>
#include <Tuple.H>
#include <Pointers.H>
#include <REAL.H>
#include <BoxArray.H>
#include <MultiFab.H>

#include <MCBndryData.H>

enum MCBC_Mode { MCHomogeneous_BC = 0, MCInhomogeneous_BC };

//@Man:
/*@Memo:
        A MCLinOp is a virtual base class for general linear operators capable
	of acting on multicomponent MultiFabs.  Eventually, this should
	be merged back into LinOp. All implementation and access functions are
	designed to make a MCLinOp object useful for representing and solving
	a set of linear equations on a union of rectangles in 2D or 3D.
*/	
/*@Doc:
        A MCLinOp contains all the information necessary to construct and
	maintain a linear operator which acts on a cell-centered MultiFab.  It
	is a virtual base class, and therefore cannot itself be
	instantiated.  The virtual functions implement "apply" 
 	and "smooth" operators (for the internal nodes), and therefore
	determine the mathematical character of the operator.  MCLinOp,
	on the other hand, determines the implementation,
	filling ghost cells and creating coarsened versions
	of the domain and operator, etc. necessary for solving linear systems.

	MCLinOp is designed specifically for representing discrete 
	approximations
	to partial differential operators.  It requires a MCBndryData object,
	which maintains ghost cell data surrounding each rectangle in the
	domain, as well as position and boundary condition type
	specifiers.  It also requires a (vector/scalar) grid spacing.  On
	levels above the base level, MCLinOp internally recognizes adjacent
	domain rectangles (using methods similar to those in the MCBndryData
	class) and when applying boundary conditions (filling ghost cell
	regions) will give preference to "valid" data (i.e. ghost cells
	covered by adjacent grids are filled with the corresponding data from
	the interior of that adjacent grid rather than via the auxiliary
	boundary instructions).

	A MCLinOp constructs a set of "levels", which are useful for linear
	solution methods such as multigrid.  On each grid, a new level is
	created by uniformly coarsening the grid structure by a factor of
	two in each coordinate direction (and then allocating and initializing
	any internal data necessary--new level grid spacing, for example).
	A MCLinOp can fill boundary ghost cells, 
	compute a "norm" and coordinate
	the "apply" and "smooth"  operations at each level.
	Note that there are the same number of levels on each grid in the
	MCLinOp's domain.

	If a boundary type specifier indicates the type "MCBndryData::
	LO_DIRICHLET", the ghost cells not covered by adjacent grids are
	filled on demand by a polynomial interpolant (of settable order) to
	compute the value specified in the MCBndryData FabSets along the edge
	of the grid boxes.  A value is put into the ghost cell of the input
	multifab that can be assumed to exist outside the domain, at equal
	spacing with the rest of the nodes.  
	"MCBndryData::LO_NEUMANN" conditions result in ghost cells filled by
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

	All the member functions of MCLinOp which interact directly with the
	apply() member take a MCBC_Mode flag.  The BC_mode can be either
	MCHomogeneous_BC, or MCInhomogeneous_BC.  It is a strict requirement of
	the linear operator that 
	MCLinOp::apply(out,in,level,bc_mode=MCHomogeneous_BC)
	acting on in=0 returns out=0.
*/
class MCLinOp
{
public:
//@ManMemo: administrative functions
//@ManDoc: allocate a MCLinOp for this box array, boundary and (uniform) spacing info
    MCLinOp(const BoxArray &ba, const MCBndryData& mgb, const Real _h);
//@ManDoc: allocate a MCLinOp for this box array, boundary data and spacing array info
    MCLinOp(const BoxArray &ba, const MCBndryData& mgb, const Real* _h);
//@ManDoc: alocate a MCLinOp w/level structure beginning at the specified level of this MCLinOp
    MCLinOp(const MCLinOp& _Lp, int level);
//@ManDoc: destructor
    virtual ~MCLinOp();

//@ManMemo: apply functions
//@ManDoc: applies level MCLinOp to "in", returns solution in "out"
    void operator()(MultiFab &out, MultiFab& in, int level)
    {
	apply(out, in, level);
    }
//@ManDoc: applies level MCLinOp to "in", returns solution in "out" at level = 0
    void operator()(MultiFab &out, MultiFab& in)
    {
        int level=0;
	apply(out, in, level);
    }
//@ManDoc: applies level MCLinOp to "in", returns "out", uses BC_mode flag
    void apply(MultiFab &out, MultiFab& in, int level = 0,
	       MCBC_Mode bc_mode=MCInhomogeneous_BC);
//@ManDoc: fills level boundary cells using BC_mode flag, int. BC data if reqd
    void applyBC(MultiFab &inout, int level = 0,
		 MCBC_Mode bc_mode=MCInhomogeneous_BC);
//@ManDoc: compute the level residual = rhsL - L(solnL)
    void residual(MultiFab &residL, const MultiFab &rhsL, MultiFab &solnL,
		  int level = 0, MCBC_Mode bc_mode=MCInhomogeneous_BC);
//@ManDoc: smooth the level system L(solnL)=rhsL
    void smooth(MultiFab &solnL, const MultiFab &rhsL,
		int level = 0, MCBC_Mode bc_mode=MCInhomogeneous_BC);

//@ManMemo: utility functions
//@ManDoc: Compute the norm of "in"
    Real norm(const MultiFab &in, int level = 0) const;

//@ManMemo: access functions
//@ManDoc: return the boundary data object
    const MCBndryData& bndryData() const
    {
	return bgb;
    }
//@ManDoc: set the boundary data object
    void bndryData(const MCBndryData& bd)
    {
	assert ( *gbox[0] == bd.boxes() );
	bgb = bd;
    }
//@ManDoc: return the box array
    const BoxArray& boxArray(int level = 0) const
    {
	assert(level < numLevels());
	return *gbox[level];
    }
//@ManDoc: return the number of grids
    int numGrids() const
    {
	return gbox[0]->length();
    }
//@ManDoc: return the number of levels
    int numLevels() const
    {
	return h.length();
    }
//@ManDoc: return the order of the boundary condition interpolant
    int maxOrder() const
    {
	return maxorder;
    }
//@ManDoc: set the order of the boundary condition interpolant
    int maxOrder(int maxorder_)
    {
	assert(maxorder_ >= 2);
	maxorder_ = (maxorder_ < 2 ? 2 : maxorder_ );
	int omaxorder = maxorder;
	maxorder = maxorder_;
	return omaxorder;
    }
//@ManDoc: construct/allocate internal data necessary for adding a new level
    virtual void prepareForLevel(int level);
//@ManDoc: remove internal data necessary for a level and all higher
    virtual void clearToLevel(int level);
  //@ManDoc: return number of components.  This is virtual since only the derived knows
  virtual int numberComponents() = 0;
  //@ManDoc: number of relaxation phases.  On derived class knows.
  virtual int numberPhases() = 0;
protected:
  //@ManDoc: number of components
  int numcomp;
  //@ManDoc: number of relaxation phases
  int numphase;
  struct PeriodicIntersectionRec {
    int srcno;
    int destno;
    Box srcbox;
    Box destbox;
#ifdef __GNUG__
    bool operator == (const struct PeriodicIntersectionRec& other) const {
      return (srcno == other.srcno ) &&
	     (destno == other.destno )&&
	     (srcbox == other.srcbox ) &&
	     (destbox == other.destbox ) ;
    }
#endif
  };
//@ManDoc: virtual to apply the level operator to the internal nodes of "in", return result in "out"
    virtual void Fapply(MultiFab &out, const MultiFab &in, int level) = 0;
//@ManDoc: virtual to carry out the level smoothing operation for L(solnL)=rhsL on internal nodes.  Modify solnL in place.
    virtual void Fsmooth(MultiFab &solnL, const MultiFab &rhsL, int level,
			 int rgbflag) = 0;
protected:
//@ManDoc: build coefficients at coarser level by interpolating "fine" (builds in appropriate node/cell centering)
    void makeCoefficients(MultiFab& crs, const MultiFab &fine, int level);
//@ManDoc: initialize MCLinOp internal data
    static void initialize();
//@ManDoc: Array (on level) of Tuples (on dimension) of grid spacings
    Array< Tuple<Real, SPACEDIM> > h;
//@ManDoc: Array (on level) of pointers to BoxArray's of MCLinOp's domain
    Array< LnClassPtr<BoxArray> > gbox;
//@ManDoc: Array (on level) of pointers to BndryRegisters along each grid for scratch data required to modify internal stencil on boundary
    Array< CpClassPtr<BndryRegister> > undrrelxr;
  //@ManDoc: Array (on level) of pointers to BndryRegisters along each grid for tangential derivative data
  Array< CpClassPtr<BndryRegister> > tangderiv;
//@ManDoc: Array (on level) of Arrays (on grid) of Arrays (on orientation) of pointers to Masks for whether boundary Fabs are covered, not_covered, outside_domain
    Array< Array< Array< Mask*> > > maskvals;
//@ManDoc: boundary data class
    MCBndryData bgb;

  //@ManDoc: Array (on level) of geometry objects.  Needed for determining whether stuff intersects on periodic domains
    Array< Geometry > geomarray;
  //@ManDoc: Array (on level) of Array of PeriodicIntersectionRec
  Array< Array< struct PeriodicIntersectionRec > > pirarray;

  //@ManDoc: internal function which constructs pirarray
  void computePeriodicIntersections(const Geometry&geo, const BoxArray& boxarr,
				   Array<struct PeriodicIntersectionRec>&pir);
    
//@ManDoc: flag (=1 if use harmonic averaged interpolation for coefficients, =0 is arithmetic averaging)
    int harmavg;
//@ManDoc: flag (>0 is verbose execution)
    int verbose;
//@ManDoc: maximum interpolation order used for constructing Dirichlet ghost node values
    int maxorder;
//@ManDoc: flag (=1 if internal data initialized)
    static bool initialized;
//@ManDoc: default value for harm_avg
    static int def_harmavg;
//@ManDoc: default value for verbose
    static int def_verbose;
//@ManDoc: default maximum BC interpolant order
    static int def_maxorder;
};


#endif
