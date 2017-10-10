
#ifndef _ABec4_H_
#define _ABec4_H_

#include <AMReX_Tuple.H>
#include <AMReX_LinOp.H>
#include <AMReX_ABec2.H>

namespace amrex {

/*
        An ABec4 tailors the description of a linear operator to apply
        the fourth-order central difference approximation to the differential
        operator, alpha*a(x).phi - beta*div[b(x).grad(phi)], on a cell-centered
        MultiFab, phi.  Here, alpha and beta are constants, a(x) is a
        cell-centered MultiFab, and b(x) is a vector of wall-centered
        MultiFabs, 1 per dimension.  phi can have multiple components, but
        all coefficient info is scalar, and norm reports on ALL components
        together.

        An ABec4 tailors a LinOp (a virtual base class for general
        linear operators) to compute the fourth-order central difference
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

class ABec4
    : public LinOp
{
public:
    //
    // constructor for box array, boundary data and vector mesh spacing
    //
    ABec4 (const BndryData& bd,
	    const Real*      _h);
    //
    // destructor
    //
    virtual ~ABec4 ();
    //
    // Compute extensive (area-weighted) flux associated with the op
    //
    virtual void compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
			   int src_comp=0, int dst_comp=0, int num_comp=1, int bnd_comp=0);

    void compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                   MultiFab& in, bool do_ApplyBC, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
		   int src_comp=0, int dst_comp=0, int num_comp=1, int bnd_comp=0);
    //
    // Compute the residual. If level==0, use 4th order otherwise use lower order
    //
    virtual void residual (MultiFab&       residL,
                           const MultiFab& rhsL,
                           MultiFab&       solnL,
                           int             level   = 0,
                           LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC,
                           bool            local   = false);

    virtual void smooth (MultiFab&       solnL,
                         const MultiFab& rhsL,
                         int             level   = 0,
                         LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);

    virtual void jacobi_smooth (MultiFab&       solnL,
                                const MultiFab& rhsL,
                                int             level   = 0,
                                LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);

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
    // Set scalar coefficients.
    //
    void setScalars (Real _alpha, Real _beta) { alpha = _alpha; beta = _beta;}
    //
    // get scalar alpha coefficient
    //
    Real get_alpha () const { return alpha; }
    //
    // get scalar beta coefficient
    //
    Real get_beta () const { return beta; }
    //
    // return reference to "a" coefficient
    //
    const MultiFab& aCoefficients (int level = 0);
    //
    // return reference to "b" coefficients
    //
    const MultiFab& bCoefficients (int level=0);
    //
    // copy _a into "a" coeffs for base level
    //
    void aCoefficients (const MultiFab& _a);
    //
    // Zero out the "a" coeffs for base level.
    //
    void ZeroACoefficients ();
    //
    // copy _b into "b" coeffs for base level
    //
    void bCoefficients (const MultiFab& _b);
    //
    // copy _b into "b" coeffs for base level at supplied gridno
    //
    void bCoefficients (const FArrayBox& _b,
                        int              gridno);
    //
    // alternative (older) version of above members
    //
    void setCoefficients (const MultiFab& _a,
                          const MultiFab& _b);
    // allocate/fill internal data for new level
    //
    virtual void prepareForLevel (int level);
    //
    // remove internal data for this level and all levels above
    //
    virtual void clearToLevel (int level);
    //
    // set flags so that a coeffs at lev and coarser require recalculation
    //
    void invalidate_a_to_level (int lev);
    //
    // set flags so that b coeffs at lev and coarser require recalculation
    //
    void invalidate_b_to_level (int lev);

    virtual Real norm (int nm = 0, int level = 0, const bool local = false);

  virtual int NumGrow(int level = 0) const;

  virtual void applyBC (MultiFab&      inout,
			int            src_comp,
			int            num_comp,
			int            level   = 0,
			LinOp::BC_Mode bc_mode = LinOp::Inhomogeneous_BC,
			bool           local   = false,
			int            bndryComp = 0);

  virtual int numLevels () const;
  virtual int numLevelsHO () const;

  virtual const BoxArray& boxArray (int level = 0) const;

  static void ca2cc(const MultiFab& ca, MultiFab& cc,
                     int sComp, int dComp, int nComp);

  static void cc2ca(const MultiFab& cc, MultiFab& ca,
                    int sComp, int dComp, int nComp);

  static void lo_cc2ec(const MultiFab& cc, MultiFab& ec,
                       int sComp, int dComp, int nComp, int dir, bool do_harm);

protected:
  void buildWorkSpace();

  MultiFab resL; // A temporary needed for relaxations at mg level = 0

    //
    // initialize a full set (a,b) of coefficients on the box array
    //
    void initCoefficients (const BoxArray &_ba, const DistributionMapping& _dm);
    //
    // compute out=L(in) at level=level
    //
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level);
    virtual void Fapply (MultiFab&       out,
			 int             dst_comp,
                         const MultiFab& in,
			 int             src_comp,
			 int             num_comp,
                         int             level);
    //
    // apply GSRB smoother to improve residual to L(solnL)=rhsL
    //
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag);
    //
    // apply Jacobi smoother to improve residual to L(solnL)=rhsL
    //
    virtual void Fsmooth_jacobi (MultiFab&       solnL,
                                 const MultiFab& rhsL,
                                 int             level);

    ABec2* LO_Op;

private:
    //
    //
    // Array (on level) of "a" coefficients
    //
    Vector< MultiFab* > acoefs;
    //
    // Array (on level) of "b" coefficients
    //
    Vector< MultiFab* > bcoefs;
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
    Vector<int> a_valid;
    //
    // Flag, can b coeffs be trusted at a level.
    //
    Vector<int> b_valid;
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
    ABec4 (const ABec4&);
    ABec4& operator= (const ABec4&);
};

}

#endif /*_ABec4_H_*/
