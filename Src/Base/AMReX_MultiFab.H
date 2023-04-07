
#ifndef BL_MULTIFAB_H
#define BL_MULTIFAB_H
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_Periodicity.H>
#include <AMReX_NonLocalBC.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#include <cstdint>

namespace amrex
{

typedef FabArray<BaseFab<float> > fMultiFab;

class iMultiFab;

/**
 * \brief
 * A collection (stored as an array) of FArrayBox objects.
 * This class is useful for storing floating point data on a domain defined by
 * a union of rectangular regions embedded in a uniform index space.
 * MultiFab class extends the function of the underlying FabArray class just
 * as the FArrayBox class extends the function of BaseFab<Real>.
 * Additional member functions are defined for I/O and simple arithmetic operations on these aggregate objects.
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
    * \brief Constructs an empty MultiFab.  Data can be defined at a later
    * time using the define member functions inherited from FabArray.  If
    * `define` is called later with a nulltpr as MFInfo's arena, the default
    * Arena `a` will be used.  If the arena in MFInfo is not a nullptr, the
    * MFInfo's arena will be used.
    */
    explicit MultiFab (Arena* a) noexcept;

    /**
    * \brief
    * Constructs a MultiFab
    * \param bs a valid region
    * \param dm a DistribuionMapping
    * \param ncomp number of components
    * \param ngrow number of cells the region grows
    * \param info MFInfo

    * The size of the FArrayBox is given by the Box grown by ngrow, and
    * the number of components is given by ncomp. If info is set to
    * not allocating memory, then no FArrayBoxes are allocated at
    * this time but can be defined later.
    */
    MultiFab (const BoxArray&            bxs,
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

    MultiFab (const BoxArray&            bxs,
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

    ~MultiFab ();

    MultiFab (MultiFab&& rhs) noexcept;
    MultiFab& operator= (MultiFab&& rhs) noexcept = default;

    MultiFab (const MultiFab& rhs) = delete;
    MultiFab& operator= (const MultiFab& rhs) = delete;

    void define (const BoxArray&            bxs,
                 const DistributionMapping& dm,
                 int                        nvar,
                 int                        ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FArrayBox>& factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const FabFactory<FArrayBox>& factory = FArrayBoxFactory());
#endif

    void define (const BoxArray&            bxs,
                 const DistributionMapping& dm,
                 int                        nvar,
                 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FArrayBox>& factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const FabFactory<FArrayBox>& factory = FArrayBoxFactory());
#endif

    MultiFab& operator= (Real r);
    //
    /**
    * \brief Returns the minimum value contained in component comp of the
    * MultiFab.  The parameter nghost determines the number of
    * boundary cells to search for the minimum.  The default is to
    * search only the valid regions of the FArrayBoxes.
    */
    [[nodiscard]] Real min (int comp,
              int nghost = 0,
              bool local = false) const;
    /**
    * \brief Identical to min() function, but confines its
    * search to intersection of Box b and the MultiFab.
    */
    [[nodiscard]] Real min (const Box& region,
              int        comp,
              int        nghost = 0,
              bool       local = false) const;
    /**
    * \brief Returns the maximum value contained in component comp of the
    * MultiFab.  The parameter nghost determines the number of
    * boundary cells to search for the maximum.  The default is to
    * search only the valid regions of the FArrayBoxes.
    */
    [[nodiscard]] Real max (int comp,
              int nghost = 0,
              bool local = false) const;
    /**
    * \brief Identical to the previous max() function, but confines its
    * search to intersection of Box b and the MultiFab.
    */
    [[nodiscard]] Real max (const Box& region,
              int        comp,
              int        nghost = 0,
              bool       local = false) const;
    /**
    * \brief Returns the maximum *absolute* value contained in
    * component comp of the MultiFab.
    */
    [[nodiscard]] Real norm0 (int comp = 0, int nghost = 0, bool local = false, bool ignore_covered = false ) const;
    [[nodiscard]] Real norminf (int comp = 0, int nghost = 0, bool local = false, bool ignore_covered = false ) const {
        return norm0(comp,nghost,local,ignore_covered);
    }

    [[nodiscard]] Real norm0 (const iMultiFab& mask, int comp = 0, int nghost = 0, bool local = false) const;
    [[nodiscard]] Real norminf (const iMultiFab& mask, int comp = 0, int nghost = 0, bool local = false) const {
        return norm0(mask,comp,nghost,local);
    }

    [[nodiscard]] Real norm0 (int comp, int ncomp, IntVect const& nghost, bool local = false,
                bool ignore_covered = false) const;

    using FabArray<FArrayBox>::norminf;

    /**
    * \brief Returns the maximum *absolute* values contained in
    * each component of "comps" of the MultiFab.  "nghost" ghost cells are used.
    */
    [[nodiscard]] Vector<Real> norm0 (const Vector<int>& comps, int nghost = 0, bool local = false, bool ignore_covered = false ) const;
    [[nodiscard]] Vector<Real> norminf (const Vector<int>& comps, int nghost = 0, bool local = false, bool ignore_covered = false) const {
        return norm0(comps,nghost,local,ignore_covered);
    }

    /**
    * \brief Returns the L1 norm of component "comp" over the MultiFab.
    * No ghost cells are used.  This version has no double counting for nodal data.
    */
    [[nodiscard]] Real norm1 (int comp, const Periodicity& period, bool ignore_covered = false) const;
    /**
    * \brief Returns the L1 norm of component "comp" over the MultiFab.
    * ngrow ghost cells are used.
    */
    [[nodiscard]] Real norm1 (int comp = 0, int ngrow = 0, bool local = false) const;
    /**
    * \brief Returns the L1 norm of each component of "comps" over the MultiFab.
    * ngrow ghost cells are used.
    */
    [[nodiscard]] Vector<Real> norm1 (const Vector<int>& comps, int ngrow = 0, bool local = false) const;
    /**
    * \brief Returns the L2 norm of component "comp" over the MultiFab.
    * No ghost cells are used.
    */
    [[nodiscard]] Real norm2 (int comp = 0) const;
    /**
    * \brief Returns the L2 norm of component "comp" over the MultiFab.
    * No ghost cells are used. This version has no double counting for nodal data.
    */
    [[nodiscard]] Real norm2 (int comp, const Periodicity& period) const;
    /**
    * \brief Returns the L2 norm of each component of "comps" over the MultiFab.
    * No ghost cells are used.
    */
    [[nodiscard]] Vector<Real> norm2 (const Vector<int>& comps) const;
    /**
    * \brief Returns the sum of component "comp" over the MultiFab -- no ghost cells are included.
    */
    [[nodiscard]] Real sum (int comp = 0, bool local = false) const;

    using FabArray<FArrayBox>::sum;

    /**
    * \brief Same as sum with local=false, but for non-cell-centered data, this
    *        skips non-unique points that are owned by multiple boxes.
    */
    [[nodiscard]] Real sum_unique (int comp = 0,
                     bool local = false,
                     const Periodicity& period = Periodicity::NonPeriodic()) const;
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

    [[nodiscard]] IntVect minIndex (int comp,
                      int nghost = 0) const;

    [[nodiscard]] IntVect maxIndex (int comp,
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
                     int numcomp, int nghost, bool local = false);

    /**
    * \brief Returns the dot product of a MultiFab with itself.
    */
    static Real Dot (const MultiFab& x, int xcomp,
                     int numcomp, int nghost, bool local = false);

    static Real Dot (const iMultiFab& mask,
                     const MultiFab& x, int xcomp,
                     const MultiFab& y, int ycomp,
                     int numcomp, int nghost, bool local = false);
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

    /**
    * \brief Swap from src to dst including nghost ghost cells.
    * The two MultiFabs MUST have the same underlying BoxArray.
    * The swap is local.
    */
    static void Swap (MultiFab&       dst,
                      MultiFab&       src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);

    static void Swap (MultiFab&       dst,
                      MultiFab&       src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      const IntVect&  nghost);

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

    using FabArray<FArrayBox>::Saxpy;

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

    using FabArray<FArrayBox>::Xpay;

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

    using FabArray<FArrayBox>::LinComb;

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

    static void AddProduct (MultiFab&       dst,
                            const MultiFab& src1,
                            int             comp1,
                            const MultiFab& src2,
                            int             comp2,
                            int             dstcomp,
                            int             numcomp,
                            const IntVect&  nghost);

    /**
    * \brief Are there any NaNs in the MF?
    * This may return false, even if the MF contains NaNs, if the machine
    * doesn't support the appropriate NaN testing functions.
    *
    * This version checks all components and grow cells.
    */
    [[nodiscard]] bool contains_nan (bool local=false) const;
    [[nodiscard]] bool contains_nan (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    [[nodiscard]] bool contains_nan (int scomp, int ncomp, const IntVect& ngrow, bool local=false) const;
    /**
    * \brief Are there any Infs in the MF?
    * This may return false, even if the MF contains Infs, if the machine
    * doesn't support the appropriate Inf testing functions.
    * This version checks all components and grow cells.
    */
    [[nodiscard]] bool contains_inf (bool local=false) const;
    [[nodiscard]] bool contains_inf (int scomp, int ncomp, int ngrow = 0, bool local=false) const;
    [[nodiscard]] bool contains_inf (int scomp, int ncomp, const IntVect& ngrow, bool local=false) const;

    /**
     * \brief Return a mask indicating how many duplicates are in each point.
     */
    [[nodiscard]] std::unique_ptr<MultiFab> OverlapMask (const Periodicity& period = Periodicity::NonPeriodic()) const;
    //! Owner is the grid with the lowest grid number containing the data.
    [[nodiscard]] std::unique_ptr<iMultiFab> OwnerMask (const Periodicity& period = Periodicity::NonPeriodic()) const;

    //! Sync up nodal data via averaging
    void AverageSync (const Periodicity& period = Periodicity::NonPeriodic());
    //! Sync up nodal data with weights
    void WeightedSync (const MultiFab& wgt, const Periodicity& period = Periodicity::NonPeriodic());
    //! Sync up nodal data with owners overriding non-owners
    void OverrideSync (const iMultiFab& msk, const Periodicity& period = Periodicity::NonPeriodic());

    using FabArray<FArrayBox>::OverrideSync;

    static void Initialize ();
    static void Finalize ();

private:
    //
    //! Some useful typedefs.
    typedef FabArrayBase::CopyComTagsContainer CopyComTagsContainer;
    typedef FabArrayBase::MapOfCopyComTagContainers MapOfCopyComTagContainers;

    void initVal ();
};

#ifndef _MSC_VER
inline void GccPlacaterMF ()
{
    std::allocator<MultiFab*> a1;
    std::allocator<MultiFab const*> a2;
    std::allocator<FabArray<FArrayBox>*> a3;
    std::allocator<FabArray<FArrayBox> const*> a4;

    amrex::ignore_unused(a1);
    amrex::ignore_unused(a2);
    amrex::ignore_unused(a3);
    amrex::ignore_unused(a4);
}
#endif

}

#endif /*BL_MULTIFAB_H*/
