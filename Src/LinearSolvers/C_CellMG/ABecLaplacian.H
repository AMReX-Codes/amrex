#ifndef _AABBecLAPLACIAN_H_
#define _AABBecLAPLACIAN_H_

//
// $Id: ABecLaplacian.H,v 1.1 1998-03-24 07:05:47 almgren Exp $
//

#include <Tuple.H>
#include <Boolean.H>
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
    : public LinOp
{
public:
//@ManMemo: administrative functions
//@ManDoc: constructor for box array, boundary data and scalar mesh spacing
    ABecLaplacian(const BoxArray &ba, const BndryData& bd, Real _h);
//@ManDoc: constructor for box array, boundary data and vector mesh spacing
    ABecLaplacian(const BoxArray &ba, const BndryData& bd, const Real * _h);

//@ManDoc: destructor
    virtual ~ABecLaplacian();

//@ManMemo: access functions
//@ManDoc: set scalar coefficients
    void setScalars(Real _alpha, Real _beta);
//@ManDoc: get scalar alpha coefficient
    Real get_alpha() const
    {
        return alpha;
    }
//@ManDoc: get scalar beta coefficient
    Real get_beta() const
    {
        return beta;
    }
//@ManDoc: return reference to "a" coefficients for base level
    const MultiFab& aCoefficients(int level=0);
//@ManDoc: return reference to "b" coefficients for base level
    const MultiFab& bCoefficients(int dir, int level=0);
//@ManDoc: copy \_a into "a" coeffs for base level
    void aCoefficients(const MultiFab& _a);
//@ManDoc: copy \_b into "b" coeffs in "dir" coordinate direction for base level
    void bCoefficients(const MultiFab& _b, int dir);

//@ManDoc: alternative (older) version of above members
    void setCoefficients(const MultiFab &_a,
                         const MultiFab &_bX,
                         const MultiFab &_bY)
    {
        aCoefficients(_a);
        bCoefficients(_bX, 0);
        bCoefficients(_bY, 1);
    }
//@ManDoc: collective version of above members (taking an array of MultiFabs for "b")
    void setCoefficients(const MultiFab& _a,
                         const MultiFab* _b)
    {
        aCoefficients(_a);
        for (int n=0; n<BL_SPACEDIM; ++n) bCoefficients(_b[n], n);
    }
//@ManDoc: allocate/fill internal data for new level
    virtual void prepareForLevel(int level);
//@ManDoc: remove internal data for this level and all levels above
    virtual void clearToLevel(int level);
//@ManDoc: set flags so that a coeffs at lev and coarser require recalculation
    void invalidate_a_to_level(int lev);
//@ManDoc: set flags so that b coeffs at lev and coarser require recalculation
    void invalidate_b_to_level(int lev);
protected:
//@ManDoc: initialize a full set (a,b) of coefficients on the box array
    void initCoefficients(const BoxArray &_ba);
//@ManMemoDoc: apply the operator
//@ManDoc: compute out=L(in) at level=level
    virtual void Fapply(MultiFab &out, const MultiFab &in, int level);
//@ManMemoDoc: smooth the residual 
//@ManDoc: apply smoother to improve residual to L(solnL)=rhsL
    virtual void Fsmooth(MultiFab &solnL, const MultiFab &rhsL, int level,
                         int rgbflag);
protected:
//@ManDoc: Array (on level) of "a" coefficients
    Array< MultiFab* > acoefs;
//@ManDoc: Array (on level) of Tuple (on dimension) of "b" coefficients
    Array< Tuple< MultiFab*, BL_SPACEDIM> > bcoefs;
//@ManDoc: Scalar "alpha" coefficient
    Real alpha;
//@ManDoc: Scalar "beta" coefficient
    Real beta;
private:
// Flag, can a coeffs be trusted at a level
    Array< bool > a_valid;
// Flag, can b coeffs be trusted at a level
    Array< bool > b_valid;
// Default value for a (MultiFab) coeficient
    static Real a_def;
// Default value for b (MultiFab) coeficient
    static Real b_def;
// Default value for alpha (scalar) coeficient
    static Real alpha_def;
// Default value for beta (scalar) coeficient
    static Real beta_def;
// Disallow copy constructors (for now...to be fixed)
    ABecLaplacian(const ABecLaplacian&);
    ABecLaplacian& operator =(const ABecLaplacian);
};

#endif
