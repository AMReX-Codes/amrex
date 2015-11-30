
#ifndef _ABecLAPLACIAN_H_
#define _ABecLAPLACIAN_H_

#include <Tuple.H>
#include <LinOp.H>
#include <BL_CXX11.H>

/*
        An ABecLaplacian tailors the description of a linear operator to apply
        the second-order central difference approximation to the differential
        operator, alpha*a(x).phi - beta*div[b(x).grad(phi)], on a cell-centered
        MultiFab, phi.  Here, alpha and beta are constants, a(x) is a
        cell-centered MultiFab, and b(x) is a vector of wall-centered
        MultiFabs, 1 per dimension.  phi can have multiple components, but
        all coefficient info is scalar, and norm reports on ALL components
        together.

        An ABecLaplacian tailors a LinOp (a virtual base class for general
        linear operators) to compute the second-order central difference
        approximation to the differential operator,

              alpha*a(x).phi - beta*div[b(x).grad(phi)]

        on a cell-centered
        MultiFab, phi.  Here, alpha and beta are constants, a(x) is a
        cell-centered MultiFab, and b(x) is a vector of wall-centered
        MultiFabs, 1 per dimension.  phi can have multiple components, but
        all coefficient info is scalar, and norm reports on ALL components
        together.

        This class provides the necessary
        implementations for applying the operator, and for 
        smoothing a guessed solution to the linear equation system,
        L(phi) = rhs (where rhs is another cell-centered MultiFab).  This
        class also provides access functions to set/control the coefficient
        MulitFabs a and b, as well as the scalars, alpha and beta.  These
        customizations are designed for 2D and 3D, with constant (but
        not necessarily equal) mesh spacing in each direction.

        Application of the operator and smoothing are
        "level" operations, and therefore must access "hidden" level data
        (such as boundary conditions, etc) as required, in concert with the
        semantics of the underlying LinOp defintions.  Accordingly, these
        implementations are "protected" and called only by the publically
        accessible members of the virtual base class, LinOp.  Note that these
        implementations may therefore assume that boundary (ghost) nodes and
        other pertinent information are filled properly by the base class
        functions prior to call.

        Defaults are as follows for the coefficients:

         alpha = 1.0
         beta = 1.0
         a (MultiFab) = 0.0
         b (MultiFab) = 1.0
        
        This class does NOT provide a copy constructor or assignment operator.
*/

class ABecLaplacian
    :
    public LinOp
{
public:
    //
    // constructor for box array, boundary data and scalar mesh spacing
    //
    ABecLaplacian (const BndryData& bd,
                   Real             _h);
    //
    // constructor for box array, boundary data and vector mesh spacing
    //
    ABecLaplacian (const BndryData& bd,
                   const Real*      _h);
    //
    // constructor for box array, boundary data and vector mesh spacing
    //
    // The ABecLaplacian assumes ownership of the pointer.
    //
    ABecLaplacian (BndryData*   bd,
                   const Real* _h);
    //
    // destructor
    //
    virtual ~ABecLaplacian () BL_OVERRIDE;
    //
    // Compute extensive (area-weighted) flux associated with the op
    //
    virtual void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
			   int src_comp=0, int dst_comp=0, int num_comp=1, int bnd_comp=0) BL_OVERRIDE;

    void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                   MultiFab& in, bool do_ApplyBC, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
		   int src_comp=0, int dst_comp=0, int num_comp=1, int bnd_comp=0);
    //
    // Set scalar coefficients.
    //
    void setScalars (Real _alpha, Real _beta) { alpha = _alpha; beta = _beta;}
    //
    // get scalar alpha coefficient
    //
    virtual Real get_alpha () const BL_OVERRIDE { return alpha; }
    //
    // get scalar beta coefficient
    //
    virtual Real get_beta () const BL_OVERRIDE { return beta; }
    //
    // return reference to "a" coefficients for base level
    //
    virtual const MultiFab& aCoefficients (int level = 0) BL_OVERRIDE;
    //
    // return reference to "b" coefficients for base level
    //
    virtual const MultiFab& bCoefficients (int dir,
					   int level=0) BL_OVERRIDE;
    //
    // copy _a into "a" coeffs for base level
    //
    void aCoefficients (const MultiFab& _a);
    //
    // Zero out the "a" coeffs for base level.
    //
    void ZeroACoefficients ();
    //
    // copy _b into "b" coeffs in "dir" coordinate direction for base level
    //
    void bCoefficients (const MultiFab& _b,
                        int             dir);
    //
    // copy _b into "b" coeffs in "dir" coordinate direction for base level
    // at supplied gridno
    //
    void bCoefficients (const FArrayBox& _b,
                        int              dir,
                        int              gridno);
    //
    // alternative (older) version of above members
    //
    void setCoefficients (const MultiFab& _a,
                          const MultiFab& _bX,
                          const MultiFab& _bY);
    //
    // collective version of above members (taking an array of MultiFabs for "b")
    //
    void setCoefficients (const MultiFab& _a,
                          const MultiFab* _b);
    void setCoefficients (const MultiFab& _a,
                          const PArray<MultiFab>& _b);
    //
    // allocate/fill internal data for new level
    //
    virtual void prepareForLevel (int level) BL_OVERRIDE;
    //
    // remove internal data for this level and all levels above
    //
    virtual void clearToLevel (int level) BL_OVERRIDE;
    //
    // set flags so that a coeffs at lev and coarser require recalculation
    //
    void invalidate_a_to_level (int lev);
    //
    // set flags so that b coeffs at lev and coarser require recalculation
    //
    void invalidate_b_to_level (int lev);

    virtual Real norm (int nm = 0, int level = 0, const bool local = false) BL_OVERRIDE;
  
protected:
    //
    // initialize a full set (a,b) of coefficients on the box array
    //
    void initCoefficients (const BoxArray &_ba);
    //
    // compute out=L(in) at level=level
    //
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level) BL_OVERRIDE;
    virtual void Fapply (MultiFab&       out,
			 int             dst_comp,
                         const MultiFab& in,
			 int             src_comp,
			 int             num_comp,
                         int             level) BL_OVERRIDE;
    //
    // apply GSRB smoother to improve residual to L(solnL)=rhsL
    //
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag) BL_OVERRIDE;
    //
    // apply Jacobi smoother to improve residual to L(solnL)=rhsL
    //
    virtual void Fsmooth_jacobi (MultiFab&       solnL,
                                 const MultiFab& rhsL,
                                 int             level) BL_OVERRIDE;
private:
    //
    //
    // Array (on level) of "a" coefficients
    //
    Array< MultiFab* > acoefs;
    //
    // Array (on level) of Tuple (on dimension) of "b" coefficients
    //
    Array< Tuple< MultiFab*, BL_SPACEDIM> > bcoefs;
    //
    // Scalar "alpha" coefficient
    //
    Real alpha;
    //
    // Scalar "beta" coefficient
    //
    Real beta;
    //
    // Flag, can a coeffs be trusted at a level.
    //
    Array<int> a_valid;
    //
    // Flag, can b coeffs be trusted at a level.
    //
    Array<int> b_valid;
    //
    // Default value for a (MultiFab) coefficient.
    //
    static Real a_def;
    //
    // Default value for b (MultiFab) coefficient.
    //
    static Real b_def;
    //
    // Default value for alpha (scalar) coefficient.
    //
    static Real alpha_def;
    //
    // Default value for beta (scalar) coefficient.
    //
    static Real beta_def;
    //
    // Disallow copy constructors (for now...to be fixed)
    //
    ABecLaplacian (const ABecLaplacian&);
    ABecLaplacian& operator= (const ABecLaplacian&);
};

#endif /*_ABecLAPLACIAN_H_*/
