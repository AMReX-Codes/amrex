
#ifndef _ABecLAPLACIAN_H_
#define _ABecLAPLACIAN_H_

//
// $Id: ABecLaplacian.H,v 1.11 2002-11-20 16:55:55 lijewski Exp $
//

#include <Tuple.H>
#include <LinOp.H>

//@Man:
/*@Memo:
        An ABecLaplacian tailors the description of a linear operator to apply
        the second-order central difference approximation to the differential
        operator, alpha*a(x).phi - beta*div[b(x).grad(phi)], on a cell-centered
        MultiFab, phi.  Here, alpha and beta are constants, a(x) is a
        cell-centered MultiFab, and b(x) is a vector of wall-centered
        MultiFabs, 1 per dimension.  phi can have multiple components, but
        all coefficient info is scalar, and norm reports on ALL components
        together.
*/        
/*@Doc:
        An ABecLaplacian tailors a LinOp (a virtual base class for general
        linear operators) to compute the second-order central difference
        approximation to the differential operator,
        \begin{center}
        alpha*a(x).phi - beta*div[b(x).grad(phi)],
        \end{center}
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
        \begin{itemize}
        \item alpha = 1.0
        \item beta = 1.0
        \item a (MultiFab) = 0.0
        \item b (MultiFab) = 1.0
        \end{itemize}
        
        This class does NOT provide a copy constructor or assignment operator.
*/

class ABecLaplacian
    :
    public LinOp
{
public:
    //
    //@ManDoc: constructor for box array, boundary data and scalar mesh spacing
    //
    ABecLaplacian (const BndryData& bd,
                   Real             _h);
    //
    //@ManDoc: constructor for box array, boundary data and vector mesh spacing
    //
    ABecLaplacian (const BndryData& bd,
                   const Real*      _h);
    //
    //@ManDoc: destructor
    //
    virtual ~ABecLaplacian ();
    //
    //@ManDoc: Compute extensive (area-weighted) flux associated with the op
    //
    virtual void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=Inhomogeneous_BC);

    void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                   MultiFab& in, const BC_Mode& bc_mode, bool do_ApplyBC);
    //
    //@ManDoc: Set scalar coefficients.
    //
    void setScalars (Real _alpha,
                     Real _beta);
    //
    //@ManDoc: get scalar alpha coefficient
    //
    Real get_alpha () const;
    //
    //@ManDoc: get scalar beta coefficient
    //
    Real get_beta () const;
    //
    //@ManDoc: return reference to "a" coefficients for base level
    //
    const MultiFab& aCoefficients (int level = 0);
    //
    //@ManDoc: return reference to "b" coefficients for base level
    //
    const MultiFab& bCoefficients (int dir,
                                   int level=0);
    //
    //@ManDoc: copy \_a into "a" coeffs for base level
    //
    void aCoefficients (const MultiFab& _a);
    //
    //@ManDoc: copy \_b into "b" coeffs in "dir" coordinate direction for base level
    //
    void bCoefficients (const MultiFab& _b,
                        int             dir);
    //
    //@ManDoc: alternative (older) version of above members
    //
    void setCoefficients (const MultiFab& _a,
                          const MultiFab& _bX,
                          const MultiFab& _bY);
    //
    //@ManDoc: collective version of above members (taking an array of MultiFabs for "b")
    //
    void setCoefficients (const MultiFab& _a,
                          const MultiFab* _b);
    //
    //@ManDoc: allocate/fill internal data for new level
    //
    virtual void prepareForLevel (int level);
    //
    //@ManDoc: remove internal data for this level and all levels above
    //
    virtual void clearToLevel (int level);
    //
    //@ManDoc: set flags so that a coeffs at lev and coarser require recalculation
    //
    void invalidate_a_to_level (int lev);
    //
    //@ManDoc: set flags so that b coeffs at lev and coarser require recalculation
    //
    void invalidate_b_to_level (int lev);

    virtual Real norm (int nm = 0, int level = 0);
  
protected:
    //
    //@ManDoc: initialize a full set (a,b) of coefficients on the box array
    //
    void initCoefficients (const BoxArray &_ba);
    //
    //@ManDoc: compute out=L(in) at level=level
    //
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level);
    //
    //@ManDoc: apply smoother to improve residual to L(solnL)=rhsL
    //
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag);
    //
    //@ManDoc: Array (on level) of "a" coefficients
    //
    Array< MultiFab* > acoefs;
    //
    //@ManDoc: Array (on level) of Tuple (on dimension) of "b" coefficients
    //
    Array< Tuple< MultiFab*, BL_SPACEDIM> > bcoefs;
    //
    //@ManDoc: Scalar "alpha" coefficient
    //
    Real alpha;
    //
    //@ManDoc: Scalar "beta" coefficient
    //
    Real beta;
    
private:
    //
    // Flag, can a coeffs be trusted at a level.
    //
    std::vector<bool> a_valid;
    //
    // Flag, can b coeffs be trusted at a level.
    //
    std::vector<bool> b_valid;
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

inline
Real
ABecLaplacian::get_alpha () const
{
    return alpha;
}

inline
Real
ABecLaplacian::get_beta () const
{
    return beta;
}

inline
void
ABecLaplacian::aCoefficients (const MultiFab& _a)
{
    BL_ASSERT(_a.ok());
    BL_ASSERT(_a.boxArray() == (acoefs[0])->boxArray());
    invalidate_a_to_level(0);
    (*acoefs[0]).copy(_a,0,0,1);
}

inline
void
ABecLaplacian::bCoefficients (const MultiFab& _b,
                              int             dir)
{
    BL_ASSERT(_b.ok());
    BL_ASSERT(_b.boxArray() == (bcoefs[0][dir])->boxArray());
    invalidate_b_to_level(0);
    (*bcoefs[0][dir]).copy(_b,0,0,1);
}

inline
const MultiFab&
ABecLaplacian::aCoefficients (int level)
{
    prepareForLevel(level);
    return *acoefs[level];
}

inline
const MultiFab&
ABecLaplacian::bCoefficients (int dir,int level)
{
    prepareForLevel(level);
    return *bcoefs[level][dir];
}

inline
void
ABecLaplacian::setCoefficients (const MultiFab &_a,
                                const MultiFab &_bX,
                                const MultiFab &_bY)
{
    aCoefficients(_a);
    bCoefficients(_bX, 0);
    bCoefficients(_bY, 1);
}

inline
void
ABecLaplacian::setCoefficients (const MultiFab& _a,
                                const MultiFab* _b)
{
    aCoefficients(_a);
    for (int n = 0; n < BL_SPACEDIM; ++n)
        bCoefficients(_b[n], n);
}

inline
void
ABecLaplacian::setScalars (Real _alpha,
                           Real _beta)
{
    alpha = _alpha;
    beta  = _beta;
}

#endif /*_ABecLAPLACIAN_H_*/
