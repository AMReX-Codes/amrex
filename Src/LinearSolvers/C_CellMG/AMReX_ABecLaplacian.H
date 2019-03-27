
#ifndef AMREX_ABecLAPLACIAN_H_
#define AMREX_ABecLAPLACIAN_H_

#include <AMReX_Array.H>
#include <AMReX_LinOp.H>

namespace amrex {

/**
        An ABecLaplacian tailors the description of a linear operator to apply
        the second-order central difference approximation to the differential
        operator, ``alpha*a(x).phi - beta*div[b(x).grad(phi)]``, on a cell-centered
        MultiFab, phi.  Here, alpha and beta are constants, a(x) is a
        cell-centered MultiFab, and b(x) is a vector of wall-centered
        MultiFabs, 1 per dimension.  phi can have multiple components, but
        all coefficient info is scalar, and norm reports on ALL components
        together.

        An ABecLaplacian tailors a LinOp (a virtual base class for general
        linear operators) to compute the second-order central difference
        approximation to the differential operator,

              ``alpha*a(x).phi - beta*div[b(x).grad(phi)]``

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

         ``alpha = 1.0
         beta = 1.0
         a (MultiFab) = 0.0
         b (MultiFab) = 1.0``

        This class does NOT provide a copy constructor or assignment operator.
*/

class ABecLaplacian
    :
    public LinOp
{
public:

    /**
    * \brief constructor for box array, boundary data and scalar mesh spacing
    *
    * \param bd
    * \param _h
    */
    ABecLaplacian (const BndryData& bd,
                   Real             _h);

    /**
    * \brief constructor for box array, boundary data and vector mesh spacing
    *
    * \param bd
    * \param _h
    */
    ABecLaplacian (const BndryData& bd,
                   const Real*      _h);

    /**
    * \brief constructor for box array, boundary data and vector mesh spacing
    *
    * The ABecLaplacian assumes ownership of the pointer.
    *
    * \param bd
    * \param _h
    */
    ABecLaplacian (BndryData*   bd,
                   const Real* _h);

    /**
    * \brief destructor
    */
    virtual ~ABecLaplacian () override;

    /**
    * \brief Compute extensive (area-weighted) flux associated with the op
    *
    * \param xflux
    * \param yflux
    * \param zflux
    * \param in
    * \param bc_mode
    * \param src_comp
    * \param dst_comp
    * \param num_comp
    * \param bnd_comp
    */
    virtual void compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
			   int src_comp=0, int dst_comp=0, int num_comp=1, int bnd_comp=0) override;

    void compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
                   MultiFab& in, bool do_ApplyBC, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
		   int src_comp=0, int dst_comp=0, int num_comp=1, int bnd_comp=0);

    /**
    * \brief Set scalar coefficients.
    *
    * \param _alpha
    * \param _beta
    */
    void setScalars (Real _alpha, Real _beta) { alpha = _alpha; beta = _beta;}

    /**
    * \brief get scalar alpha coefficient
    */
    virtual Real get_alpha () const override { return alpha; }

    /**
    * \brief get scalar beta coefficient
    */
    virtual Real get_beta () const override { return beta; }

    /**
    * \brief return reference to "a" coefficients for base level
    *
    * \param level
    */
    virtual const MultiFab& aCoefficients (int level = 0) override;

    /**
    * \brief return reference to "b" coefficients for base level
    *
    * \param dir
    * \param level
    */
    virtual const MultiFab& bCoefficients (int dir,
					   int level=0) override;

    /**
    * \brief copy _a into "a" coeffs for base level
    *
    * \param _a
    */
    void aCoefficients (const MultiFab& _a);

    /**
    * \brief Zero out the "a" coeffs for base level.
    */
    void ZeroACoefficients ();

    /**
    * \brief copy _b into "b" coeffs in "dir" coordinate direction for base level
    *
    * \param _b
    * \param dir
    */
    void bCoefficients (const MultiFab& _b,
                        int             dir);

    /**
    * \brief copy _b into "b" coeffs in "dir" coordinate direction for base level
    * at supplied gridno
    *
    * \param _b
    * \param dir
    * \param gridno
    */
    void bCoefficients (const FArrayBox& _b,
                        int              dir,
                        int              gridno);

    /**
    * \brief alternative (older) version of above members
    *
    * \param _a
    * \param _bX
    * \param _bY
    */
    void setCoefficients (const MultiFab& _a,
                          const MultiFab& _bX,
                          const MultiFab& _bY);

    /**
    * \brief collective version of above members (taking an array of MultiFabs for "b")
    *
    * \param _a
    * \param _b
    */
    void setCoefficients (const MultiFab& _a,
                          const MultiFab* _b);
    void setCoefficients (const MultiFab& _a,
                          const Vector<MultiFab*>& _b);

    /**
    * \brief allocate/fill internal data for new level
    *
    * \param level
    */
    virtual void prepareForLevel (int level) override;

    /**
    * \brief remove internal data for this level and all levels above
    *
    * \param level
    */
    virtual void clearToLevel (int level) override;

    /**
    * \brief set flags so that a coeffs at lev and coarser require recalculation
    *
    * \param lev
    */
    void invalidate_a_to_level (int lev);

    /**
    * \brief set flags so that b coeffs at lev and coarser require recalculation
    *
    * \param lev
    */
    void invalidate_b_to_level (int lev);

    virtual Real norm (int nm = 0, int level = 0, const bool local = false) override;

protected:

    /**
    * \brief initialize a full set (a,b) of coefficients on the box array
    *
    * \param _ba
    * \param _dm
    */
    void initCoefficients (const BoxArray &_ba, const DistributionMapping& _dm);

    /**
    * \brief compute out=L(in) at level=level
    *
    * \param out
    * \param in
    * \param level
    */
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level) override;
    virtual void Fapply (MultiFab&       out,
			 int             dst_comp,
                         const MultiFab& in,
			 int             src_comp,
			 int             num_comp,
                         int             level) override;

    /**
    * \brief apply GSRB smoother to improve residual to L(solnL)=rhsL
    *
    * \param solnL
    * \param rhsL
    * \param level
    * \param rgbflag
    */
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag) override;

    /**
    * \brief apply Jacobi smoother to improve residual to L(solnL)=rhsL
    *
    * \param solnL
    * \param rhsL
    * \param level
    */
    virtual void Fsmooth_jacobi (MultiFab&       solnL,
                                 const MultiFab& rhsL,
                                 int             level) override;
private:

    //! Vector (on level) of "a" coefficients
    Vector< MultiFab* > acoefs;

    //! Vector (on level) of Array (on dimension) of "b" coefficients
    Vector< Array< MultiFab*, BL_SPACEDIM> > bcoefs;

    //! Scalar "alpha" coefficient
    Real alpha;

    //! Scalar "beta" coefficient
    Real beta;

    //! Flag, can a coeffs be trusted at a level.
    Vector<int> a_valid;

    //! Flag, can b coeffs be trusted at a level.
    Vector<int> b_valid;

    //! Default value for a (MultiFab) coefficient.
    static Real a_def;

    //! Default value for b (MultiFab) coefficient.
    static Real b_def;

    //! Default value for alpha (scalar) coefficient.
    static Real alpha_def;

    //! Default value for beta (scalar) coefficient.
    static Real beta_def;

    /**
    * \brief Disallow copy constructors (for now...to be fixed)
    *
    */
    ABecLaplacian (const ABecLaplacian&);
    ABecLaplacian& operator= (const ABecLaplacian&);
};

}

#endif /*_ABecLAPLACIAN_H_*/
