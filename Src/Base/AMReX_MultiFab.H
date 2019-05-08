
#ifndef BL_MULTIFAB_H
#define BL_MULTIFAB_H


#include <cstdint>

#include <AMReX_BLassert.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_Periodicity.H>

namespace amrex
{

class iMultiFab;

/**
 * \brief
 * A collection (stored as an array) of FArrayBox objects.
 * This class is useful for storing floating point data on a domain defined by
 * a union of rectangular regions embedded in a uniform index space.
 * MultiFab class extends the function of the underlying FabArray class just
 * as the FArrayBox class extends the funtion of BaseFab<Real>.
 *Additional member functions are defined for I/O and simple arithmetic operations on these aggregate objects.
 * This class does NOT provide a copy constructor or assignment operator.
 */

class MultiFab
    :
    public FabArray<FArrayBox>
{
public:
    /**
    * \brief Constructs an empty MultiFab.  Data can be defined at a later
    * time using the define member functions inherited
    * from FabArray.
    */
    MultiFab () noexcept;
    /**
    * \brief
    * Constructs a MultiFab
    * \param bs a valid region
    * \param dm a DistribuionMapping
    * \param ncomp number of components
    * \param ngrow number of cells the region grows
    * \param info MFInfo

    * The size of the FArrayBox is given by the Box grown by ngrow and
    * the number of components is given by ncomp.  If info is set to
    * not allocateing memroy, then no FArrayBoxes are allocated at
    * this time but can be defined later.
    */
    MultiFab (const BoxArray&            bs,
              const DistributionMapping& dm,
              int                        ncomp,
              int                        ngrow,
#ifdef AMREX_STRICT_MODE
	      const MFInfo&              info,
              const FabFactory<FArrayBox>& factory);
#else
	      const MFInfo&              info = MFInfo(),
              const FabFactory<FArrayBox>& factory = FArrayBoxFactory());
#endif

    MultiFab (const BoxArray&            bs,
              const DistributionMapping& dm,
              int                        ncomp,
              const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
	      const MFInfo&              info,
              const FabFactory<FArrayBox>& factory);
#else
	      const MFInfo&              info = MFInfo(),
              const FabFactory<FArrayBox>& factory = FArrayBoxFactory());
#endif

    /**
     * \brief Make an alias MultiFab. maketype must be
     * amrex::make_alias.  scomp is the starting component of the
     * alias and ncomp is the number of components in the new aliasing
     * MultiFab.
     */
    MultiFab (const MultiFab& rhs, MakeType maketype, int scomp, int ncomp);

    virtual ~MultiFab () override;

    MultiFab (MultiFab&& rhs) noexcept;
    MultiFab& operator= (MultiFab&& rhs) noexcept = default;

    MultiFab (const MultiFab& rhs) = delete;
    MultiFab& operator= (const MultiFab& rhs) = delete;

    virtual void define (const BoxArray&            bxs,
			 const DistributionMapping& dm,
			 int                        nvar,
			 int                        ngrow,
#ifdef AMREX_STRICT_MODE
			 const MFInfo&              info,
                         const FabFactory<FArrayBox>& factory) override;
#else
			 const MFInfo&              info = MFInfo(),
                         const FabFactory<FArrayBox>& factory = FArrayBoxFactory()) override;
#endif

    virtual void define (const BoxArray&            bxs,
			 const DistributionMapping& dm,
			 int                        nvar,
			 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
			 const MFInfo&              info,
                         const FabFactory<FArrayBox>& factory) override;
#else
			 const MFInfo&              info = MFInfo(),
                         const FabFactory<FArrayBox>& factory = FArrayBoxFactory()) override;
#endif

    void operator= (Real r);
    //
    //! Initialize FArrayBoxes
    void initVal ();
    //
    /**
    * \brief Returns the minimum value contained in component comp of the
    * MultiFab.  The parameter nghost determines the number of
    * boundary cells to search for the minimum.  The default is to
    * search only the valid regions of the FArrayBoxes.
    */
    Real min (int comp,
              int nghost = 0,
	      bool local = false) const;
    /**
    * \brief Identical to min() function, but confines its
    * search to intersection of Box b and the MultiFab.
    */
    Real min (const Box& b,
              int        comp,
              int        nghost = 0,
	      bool       local = false) const;
    /**
    * \brief Returns the maximum value contained in component comp of the
    * MultiFab.  The parameter nghost determines the number of
    * boundary cells to search for the maximum.  The default is to
    * search only the valid regions of the FArrayBoxes.
    */
    Real max (int comp,
              int nghost = 0,
	      bool local = false) const;
    /**
    * \brief Identical to the previous max() function, but confines its
    * search to intersection of Box b and the MultiFab.
    */
    Real max (const Box& b,
              int        comp,
              int        nghost = 0,
	      bool       local = false) const;
    /**
    * \brief Returns the maximum *absolute* value contained in
    * component comp of the MultiFab.
    */
    Real norm0 (int comp = 0, int nghost = 0, bool local = false) const;

    Real norm0 (const iMultiFab& mask, int comp = 0, int nghost = 0, bool local = false) const;

    /**
    * \brief Returns the maximum *absolute* values contained in
    * each component of "comps" of the MultiFab.  No ghost cells are used.
    */

    Vector<Real> norm0 (const Vector<int>& comps, int nghost = 0, bool local = false) const;
    /**
    * \brief Returns the L1 norm of component "comp" over the MultiFab.
    * No ghost cells are used.  This version has no double counting for nodal data.
    */
    Real norm1 (int comp, const Periodicity& period) const;
    /**
    * \brief Returns the L1 norm of component "comp" over the MultiFab.
    * ngrow ghost cells are used.
    */
    Real norm1 (int comp = 0, int ngrow = 0, bool local = false) const;
    /**
    * \brief Returns the L1 norm of each component of "comps" over the MultiFab.
    * ngrow ghost cells are used.
    */
    Vector<Real> norm1 (const Vector<int>& comps, int ngrow = 0, bool local = false) const;
    /**
    * \brief Returns the L2 norm of component "comp" over the MultiFab.
    * No ghost cells are used.
    */
    Real norm2 (int comp = 0) const;
    /**
    * \brief Returns the L2 norm of component "comp" over the MultiFab.
    * No ghost cells are used. This version has no double counting for nodal data.
    */
    Real norm2 (int comp, const Periodicity& period) const;
    /**
    * \brief Returns the L2 norm of each component of "comps" over the MultiFab.
    * No ghost cells are used.
    */
    Vector<Real> norm2 (const Vector<int>& comps) const;
    /**
    * \brief Returns the sum of component "comp" over the MultiFab -- no ghost cells are included.
    */
    Real sum (int comp = 0, bool local = false) const;
    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * specified subregion of the MultiFab.  The subregion consists
    * of the num_comp components starting at component comp.
    * The value of nghost specifies the number of cells in the
    * boundary region of each FArrayBox in the subregion that should
    * be modified.
    */
    void plus (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    /**
    * \brief Identical to the previous version of plus(), with the
    * restriction that the subregion is further constrained to
    * the intersection with Box region.
    */
    void plus (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * valid region of each component of the MultiFab.  The value
    * of nghost specifies the number of cells in the boundary
    * region that should be modified.
    */
    void plus (Real val,
               int  nghost);
    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * valid region of each component of the MultiFab, that also
    * intersects the Box region.  The value of nghost specifies the
    * number of cells in the boundary region of each FArrayBox in
    * the subregion that should be modified.
    */
    void plus (Real       val,
               const Box& region,
               int        nghost);
    /**
    * \brief Scales the value of each cell in the specified subregion of the
    * MultiFab by the scalar val (a[i] <- a[i]*val). The subregion
    * consists of the num_comp components starting at component comp.
    * The value of nghost specifies the number of cells in the
    * boundary region of each FArrayBox in the subregion that should
    * be modified.
    */
    void mult (Real val,
               int  comp,
               int  num_comp,
               int  nghost = 0);
    /**
    * \brief Identical to the previous version of mult(), with the
    * restriction that the subregion is further constrained to the
    * intersection with Box region.  The value of nghost specifies the
    * number of cells in the boundary region of each FArrayBox in
    * the subregion that should be modified.
    */
    void mult (Real       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);
    /**
    * \brief Scales the value of each cell in the valid region of each
    * component of the MultiFab by the scalar val (a[i] <- a[i]*val).
    * The value of nghost specifies the number of cells in the
    * boundary region that should be modified.
    */
    void mult (Real val,
               int  nghost = 0);
    /**
    * \brief Scales the value of each cell in the valid region of each
    * component of the MultiFab by the scalar val (a[i] <- a[i]*val),
    * that also intersects the Box region.  The value of nghost
    * specifies the number of cells in the boundary region of each
    * FArrayBox in the subregion that should be modified.
    */
    void mult (Real       val,
               const Box& region,
               int        nghost = 0);
    /**
    * \brief Replaces the value of each cell in the specified subregion of
    * the MultiFab with its reciprocal multiplied by the value of
    * numerator. The subregion consists of the num_comp components
    * starting at component comp.  The value of nghost specifies the
    * number of cells in the boundary region of each FArrayBox in the
    * subregion that should be modified.
    */
    void invert (Real numerator,
                 int  comp,
                 int  num_comp,
                 int  nghost = 0);
    /**
    * \brief Identical to the previous version of invert(), with the
    * restriction that the subregion is further constrained to the
    * intersection with Box region.  The value of nghost specifies the
    * number of cells in the boundary region of each FArrayBox in the
    * subregion that should be modified.
    */
    void invert (Real       numerator,
                 const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    /**
    * \brief Replaces the value of each cell in the specified subregion of
    * the MultiFab with its reciprocal multiplied by the value of
    * numerator.  The value of nghost specifies the number of cells
    * in the boundary region that should be modified.
    */
    void invert (Real numerator,
                 int  nghost);
    /**
    * \brief Replaces the value of each cell in the specified subregion of
    * the MultiFab, that also intersects the Box region, with its
    * reciprocal multiplied by the value of numerator.  The value
    * of nghost specifies the number of cells in the boundary region
    * of each FArrayBox in the subregion that should be modified.
    */
    void invert (Real       numerator,
                 const Box& region,
                 int        nghost);
    /**
    * \brief Negates the value of each cell in the specified subregion of
    * the MultiFab.  The subregion consists of the num_comp
    * components starting at component comp.  The value of nghost
    * specifies the number of cells in the boundary region of each
    * FArrayBox in the subregion that should be modified.
    */
    void negate (int comp,
                 int num_comp,
                 int nghost = 0);
    /**
    * \brief Identical to the previous version of negate(), with the
    * restriction that the subregion is further constrained to
    * the intersection with Box region.
    */
    void negate (const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);
    /**
    * \brief Negates the value of each cell in the valid region of
    * the MultiFab.  The value of nghost specifies the number of
    * cells in the boundary region that should be modified.
    */
    void negate (int nghost = 0);
    /**
    * \brief Negates the value of each cell in the valid region of
    * the MultiFab that also intersects the Box region.  The value
    * of nghost specifies the number of cells in the boundary region
    * that should be modified.
    */
    void negate (const Box& region,
                 int        nghost = 0);

    IntVect minIndex (int comp,
                      int nghost = 0) const;

    IntVect maxIndex (int comp,
                      int nghost = 0) const;
    /**
    * \brief This function adds the values of the cells in mf to the corresponding
    * cells of this MultiFab.  mf is required to have the same BoxArray or
    * "valid region" as this MultiFab.  The addition is done only to num_comp
    * components, starting with component number strt_comp.  The parameter
    * nghost specifies the number of boundary cells that will be modified.
    * If nghost == 0, only the valid region of each FArrayBox will be
    * modified.
    */
    void plus (const MultiFab& mf,
               int             strt_comp,
               int             num_comp,
               int             nghost);
    /**
    * \brief This function subtracts the values of the cells in mf from the
    * corresponding cells of this MultiFab.  mf is required to have the
    * same BoxArray or "valid region" as this MultiFab.  The subtraction is
    * done only to num_comp components, starting with component number
    * strt_comp.  The parameter nghost specifies the number of boundary
    * cells that will be modified.  If nghost == 0, only the valid region of
    * each FArrayBox will be modified.
    */
    void minus (const MultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost);
    /**
    * \brief This function divides the values of the cells in mf from the
    * corresponding cells of this MultiFab.  mf is required to have the
    * same BoxArray or "valid region" as this MultiFab.  The division is
    * done only to num_comp components, starting with component number
    * strt_comp.  The parameter nghost specifies the number of boundary
    * cells that will be modified.  If nghost == 0, only the valid region of
    * each FArrayBox will be modified.  Note, nothing is done to protect
    * against divide by zero.
    */
    void divide (const MultiFab& mf,
		 int             strt_comp,
		 int             num_comp,
		 int             nghost);
    /**
    * \brief Returns the dot product of two MultiFabs.
    */
    static Real Dot (const MultiFab& x, int xcomp,
		     const MultiFab& y, int ycomp,
		     int num_comp, int nghost, bool local = false);

    static Real Dot (const iMultiFab& mask,
                     const MultiFab& x, int xcomp,
		     const MultiFab& y, int ycomp,
		     int num_comp, int nghost, bool local = false);
    /**
    * \brief Add src to dst including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Add (MultiFab&       dst,
		     const MultiFab& src,
		     int             srccomp,
		     int             dstcomp,
		     int             numcomp,
		     int             nghost);

    static void Add (MultiFab&       dst,
                     const MultiFab& src,
                     int             srccomp,
                     int             dstcomp,
                     int             numcomp,
                     const IntVect&  nghost);
    /**
    * \brief Copy from src to dst including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    * The copy is local.  The parallel copy function is in FabArray.H
    */
    static void Copy (MultiFab&       dst,
                      const MultiFab& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);

    static void Copy (MultiFab&       dst,
                      const MultiFab& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      const IntVect&  nghost);

#ifdef USE_PERILLA
    static void Copy (MultiFab&       dst,
                const MultiFab& src,
                int             f,
                int             srccomp,
                int             dstcomp,
                int             numcomp,
                const Box&      bx);
#endif

    /**
    * \brief Subtract src from dst including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Subtract (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);

    static void Subtract (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  const IntVect&  nghost);
    /**
    * \brief Multiply dst by src including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Multiply (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);

    static void Multiply (MultiFab&       dst,
			  const MultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  const IntVect&  nghost);
    /**
    * \brief Divide dst by src including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    */
    static void Divide (MultiFab&       dst,
			const MultiFab& src,
			int             srccomp,
			int             dstcomp,
			int             numcomp,
			int             nghost);

    static void Divide (MultiFab&       dst,
			const MultiFab& src,
			int             srccomp,
			int             dstcomp,
			int             numcomp,
			const IntVect&  nghost);
    /**
    * \brief dst += a*src
    */
    static void Saxpy (MultiFab&       dst,
		       Real            a,
		       const MultiFab& src,
		       int             srccomp,
		       int             dstcomp,
		       int             numcomp,
		       int             nghost);
    /**
    * \brief dst = src + a*dst
    */
    static void Xpay (MultiFab&       dst,
		      Real            a,
		      const MultiFab& src,
		      int             srccomp,
		      int             dstcomp,
		      int             numcomp,
		      int             nghost);
    /**
    * \brief dst = a*x + b*y
    */
    static void LinComb (MultiFab&       dst,
			 Real            a,
			 const MultiFab& x,
			 int             xcomp,
			 Real            b,
			 const MultiFab& y,
			 int             ycomp,
			 int             dstcomp,
			 int             numcomp,
			 int             nghost);
    /**
    * \brief dst += src1*src2
    */
    static void AddProduct (MultiFab&       dst,
			    const MultiFab& src1,
			    int             comp1,
			    const MultiFab& src2,
			    int             comp2,
			    int             dstcomp,
			    int             numcomp,
			    int             nghost);

    /**
    * \brief Are there any NaNs in the MF?
    * This may return false, even if the MF contains NaNs, if the machine
    * doesn't support the appropriate NaN testing functions.
    *
    * This version checks all components and grow cells.
    */
    bool contains_nan (bool local=false) const;
    bool contains_nan (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    bool contains_nan (int scomp, int ncomp, const IntVect& ngrow, bool local=false) const;
    /**
    * \brief Are there any Infs in the MF?
    * This may return false, even if the MF contains Infs, if the machine
    * doesn't support the appropriate Inf testing functions.
    * This version checks all components and grow cells.
    */
    bool contains_inf (bool local=false) const;
    bool contains_inf (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    bool contains_inf (int scomp, int ncomp, const IntVect& ngrow, bool local=false) const;

    /**
    * \brief This tests on whether the MultiFab is fully nodal.
    */
    bool is_nodal () const noexcept;
    /**
    * \brief This tests on whether the MultiFab is nodal in direction dir.
    */
    bool is_nodal (int dir) const noexcept;
    /**
    * \brief This tests on whether the MultiFab is cell-centered.
    */
    bool is_cell_centered () const noexcept;

    /**
    * \brief Sum values in overlapped cells.  The destination is limited to valid cells.
    */
    void SumBoundary (const Periodicity& period = Periodicity::NonPeriodic());

    void SumBoundary (int scomp, int ncomp, const Periodicity& period = Periodicity::NonPeriodic());

    /**
    * \brief Sum values in overlapped cells.  The destination is limited to valid + ngrow cells.
    */
    void SumBoundary (int scomp, int ncomp, IntVect const& ngrow,
                      const Periodicity& period = Periodicity::NonPeriodic());

    /**
     * \brief Return a mask indicating how many duplicates are in each point.
     */
    std::unique_ptr<MultiFab> OverlapMask (const Periodicity& period = Periodicity::NonPeriodic()) const;
    //! Owner is the grid with the lowest grid number containing the data.
    std::unique_ptr<iMultiFab> OwnerMask (const Periodicity& period = Periodicity::NonPeriodic()) const;

    //! Sync up nodal data via averaging
    void AverageSync (const Periodicity& period = Periodicity::NonPeriodic());
    //! Sync up nodal data with weights
    void WeightedSync (const MultiFab& wgt, const Periodicity& period = Periodicity::NonPeriodic());
    //! Sync up nodal data with owners overriding non-owners
    void OverrideSync (const Periodicity& period = Periodicity::NonPeriodic());
    void OverrideSync (const iMultiFab& msk, const Periodicity& period = Periodicity::NonPeriodic());

    static void Initialize ();
    static void Finalize ();

private:
    //
    //! Some useful typedefs.
    typedef FabArrayBase::CopyComTagsContainer CopyComTagsContainer;
    typedef FabArrayBase::MapOfCopyComTagContainers MapOfCopyComTagContainers;
};

inline void GccPlacaterMF ()
{
    std::allocator<MultiFab*> a1;
    std::allocator<MultiFab const*> a2;
    std::allocator<FabArray<FArrayBox>*> a3;
    std::allocator<FabArray<FArrayBox> const*> a4;
}

//!  This is a special version of FillBoundary for warpx
void FillBoundary (Vector<MultiFab*> const& mf, const Periodicity& period);

}

#endif /*BL_MULTIFAB_H*/
