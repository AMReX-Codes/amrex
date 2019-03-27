#ifndef AMREX_DivVis_H_
#define AMREX_DivVis_H_

#include <AMReX_Array.H>
#include <AMReX_MCLinOp.H>

namespace amrex {

/**

  \brief A DivVis tailors the description of a linear operator to apply
  the second-order central difference approximation to the differential
  operator, ``alpha*a(x).u(x) - beta*div[tau(u(x))]``, on a cell-centered
  MultiFab, u(x).  Here, alpha and beta are constants, a(x) is a
  cell-centered MultiFab, and b(x) is a vector of wall-centered
  MultiFabs, 1 per dimension.  u(x) can have multiple components, but
  all coefficient info is scalar, and norm reports on ALL components
  together.
  tau is a three by three tensor

  ``b(x)*(u_{i,j}+u_{j,i})``

  A DivVis tailors a LinOp (a virtual base class for general
  linear operators) to compute the second-order central difference
  approximation to the differential operator,

  ``alpha*a(x).u(x) - beta*div[tau(u(x))]``

  on a cell-centered
  MultiFab, u(x).  Here, alpha and beta are constants, a(x) is a
  cell-centered MultiFab, and b(x) is a vector of wall-centered
  MultiFabs, 1 per dimension.  u(x) can have multiple components, but
  all coefficient info is scalar, and norm reports on ALL components
  together.
  Tau is the 3 by 3 tensor

  ``b(x)*(u_{i,j}+u_{j,i})``


  This class provides the necessary
  implementations for applying the operator, and for
  smoothing a guessed solution to the linear equation system,
  L(u(x)) = rhs (where rhs is another cell-centered MultiFab).  This
  class also provides access functions to set/control the coefficient
  MultiFabs a and b, as well as the scalars, alpha and beta.  These
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
  */

class DivVis
    :
    public MCLinOp
{
public:

    /**
    * \brief Constructor for box array, boundary data and scalar mesh spacing.
    *
    * \param bd
    * \param _h
    */
    DivVis (const BndryData& bd,
	    Real             _h);

    /**
    * \brief Constructor for box array, boundary data and vector mesh spacing.
    *
    * \param bd
    * \param _h
    */
    DivVis (const BndryData& bd,
	    const Real*      _h);

    /**
    * \brief The destructor.
    */
    virtual ~DivVis () override;

    /**
    * \brief Set scalar coefficients.
    *
    * \param _alpha
    * \param _beta
    */
    void setScalars (Real _alpha,
		     Real _beta);

    /**
    * \brief Get scalar alpha coefficient.
    */
    Real get_alpha () const;

    /**
    * \brief Get scalar beta coefficient.
    */
    Real get_beta () const;

    /**
    * \brief Return reference to "a" coefficients for base level.
    *
    * \param level
    */
    const MultiFab& aCoefficients (int level = 0);

    /**
    * \brief Return reference to "b" coefficients for base level.
    *
    * \param dir
    * \param level
    */
    const MultiFab& bCoefficients (int dir,
				   int level = 0);

    /**
    * \brief Copy _a into "a" coeffs for base level.
    *
    * \param _a
    */
    void aCoefficients (const MultiFab& _a);

    /**
    * \brief Zero out the "a" coeffs for base level.
    */
    void ZeroACoefficients ();

    /**
    * \brief Copy _b into "b" coeffs in "dir" coordinate direction for base level.
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

    //! Alternative (older) version of above members.
    void setCoefficients (const MultiFab &_a,
			  const MultiFab &_bX,
			  const MultiFab &_bY
#if BL_SPACEDIM==3
			  ,const MultiFab &_bZ
#endif
	);

    /**
    * \brief Collective version of above members (taking an array of MultiFabs for "b").
    *
    * \param _a
    * \param _b
    */
    void setCoefficients (const MultiFab& _a,
			  const MultiFab* _b);

    /**
    * \brief Allocate/fill internal data for new level.
    *
    * \param level
    */
    virtual void prepareForLevel (int level) override;

    /**
    * \brief Remove internal data for this level and all levels above.
    *
    * \param level
    */
    virtual void clearToLevel (int level) override;

    /**
    * \brief Set flags so that a coeffs at lev and coarser require recalculation.
    *
    * \param lev
    */
    void invalidate_a_to_level (int lev);

    /**
    * \brief Set flags so that b coeffs at lev and coarser require recalculation.
    *
    * \param lev
    */
    void invalidate_b_to_level (int lev);

    //! Computes flux associated with operator.
    void compFlux (MultiFab &xflux,
		   MultiFab &yflux,
#if BL_SPACEDIM>2
		   MultiFab &zflux,
#endif
		   MultiFab& in);

protected:

    /**
    * \brief Initialize a full set (a,b) of coefficients on the box array.
    *
    * \param &_ba
    * \param _dm
    */
    void initCoefficients (const BoxArray &_ba, const DistributionMapping& _dm);

    /**
    * \brief Compute out=L (in) at level=level.
    *
    * \param out
    * \param in
    * \param level
    */
    virtual void Fapply (MultiFab&       out,
			 const MultiFab& in,
			 int             level) override;

    /**
    * \brief Apply smoother to improve residual to L(solnL)=rhsL.
    *
    * \param solnL
    * \param rhsL
    * \param level
    * \param phaseflag
    */
    virtual void Fsmooth (MultiFab&       solnL,
			  const MultiFab& rhsL,
			  int             level,
			  int             phaseflag) override;

    /**
    * \brief Return number of components.  This is virtual since only the derived knows.
    */
    virtual int numberComponents () override;

    /**
    * \brief Number of relaxation phases.  On derived class knows.
    */
    virtual int numberPhases () override;

    /**
    * \brief Helper function for object construction.
    *
    * \param _h
    */
    void initConstruct (const Real* _h);

protected:

    //! Vector (on level) of "a" coefficients
    Vector< MultiFab* > acoefs;

    //! Vector (on level) of Array (on dimension) of "b" coefficients
    Vector< Array< MultiFab*, BL_SPACEDIM> > bcoefs;

    //! Scalar "alpha" coefficient
    Real alpha;

    //! Scalar "beta" coefficient
    Real beta;

private:

    //! Flag, can a coeffs be trusted at a level.
    Vector<int> a_valid;

    //! Flag, can b coeffs be trusted at a level.
    Vector<int> b_valid;

    //! Default value for a (MultiFab) coeficient.
    static Real a_def;

    //! Default value for b (MultiFab) coeficient.
    static Real b_def;

    //! Default value for alpha (scalar) coeficient.
    static Real alpha_def;

    //! Default value for beta (scalar) coeficient.
    static Real beta_def;

    /**
    * \brief Disallow copy constructors (for now...to be fixed).
    *
    */
    DivVis (const DivVis&);
    DivVis& operator = (const DivVis&);
};

inline
Real
DivVis::get_alpha () const
{
    return alpha;
}

inline
Real
DivVis::get_beta () const
{
    return beta;
}

inline
void
DivVis::setCoefficients (const MultiFab &_a,
			 const MultiFab &_bX,
			 const MultiFab &_bY
#if BL_SPACEDIM==3
			 ,const MultiFab &_bZ
#endif
    )
{
    aCoefficients(_a);
    bCoefficients(_bX, 0);
    bCoefficients(_bY, 1);
#if BL_SPACEDIM==3
    bCoefficients(_bZ, 2);
#endif
}

inline
void
DivVis::setCoefficients (const MultiFab& _a,
			 const MultiFab* _b)
{
    aCoefficients(_a);
    for (int n = 0; n < BL_SPACEDIM; ++n)
        bCoefficients(_b[n], n);
}

}

#endif /*_DivVis_H_*/
