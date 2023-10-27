#ifndef AMREX_FABARRAYSET_H_
#define AMREX_FABARRAYSET_H_

#include <AMReX_FabArray.H>
#include <AMReX_IndexType.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

namespace amrex {

/**
 * \brief A set of FabArray's with distinct IndexType's living on the same mesh.
 *
 * Intended for use as a single object to use as input/output in a linear solver
 * class that inherits from MLLinOpT<FabArraySet< \p MF >>.
 *
 * Not to be confused with FabSet.
 *
 * \tparam MF FabArray or derived class
 */
template <class MF>
class FabArraySet
{
public:

    using value_type = typename MF::value_type;
    using fab_type = typename MF::fab_type;
    using FAB = typename MF::fab_type;

    //! debug variable
    // using MF = typename amrex::MultiFab;
    // using value_type = typename amrex::Real;
    // using fab_type = typename amrex::FArrayBox;
    // using FAB = typename amrex::FArrayBox;

    //! Constructs an empty FabArraySet.
    FabArraySet () noexcept;

    /**
     * \brief Constructs a new FabArraySet.
     *
     * The size of the each FabArray ( \p MF ) element is given by the BoxArray
     * \p bxs modified by the IndexType's \p ixtypes and grown by \p ngrow.
     * Each \p MF has \p ncomp components, and the size of \p ixtypes determines
     * the number of \p MFs in the set. If \p info is set to not allocating memory,
     * then no FABs are allocated at this time but can be defined later.
     *
     * \param bxs a valid region
     * \param ixtypes IndexType for each \p MF. Its size determines the number of \p MFs in the set.
     * \param dm a DistributionMapping
     * \param ncomp number of components for each \p MF (not to be confused with the number of \p MFs in the set)
     * \param ngrow number of cells the region grows
     * \param info MFInfo
     * \param factory FabFactory
     */
    FabArraySet (const BoxArray&            bxs,
                 const Vector<IndexType>    ixtypes,
                 const DistributionMapping& dm,
                 int                        nvar,
                 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FAB>&     factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    /**
     * \brief Make an alias FabArraySet<MF>.
     *
     * \param rhs Make an alias of rhs
     * \param maketype must be amrex::make_alias
     * \param scomp starting component in the MFs of the alias
     * \param ncomp number of components in the MFs in the new aliasing FabArraySet
     */
    FabArraySet (const FabArraySet& rhs, MakeType maketype, int scomp, int ncomp);

    //! The destructor -- deletes all FabArrays in the set.
    virtual ~FabArraySet () = default;

    FabArraySet (FabArraySet&& rhs) noexcept;
    FabArraySet& operator= (FabArraySet&& rhs) noexcept;
    FabArraySet (const FabArraySet& rhs) = delete;
    FabArraySet& operator= (const FabArraySet& rhs) = delete;

    /**
     * \brief Define this FabArraySet identically to that performed by
     * the constructor having an analogous function signature.
     * This is only valid if this FabArraySet was defined using
     * the default constructor.
     */
    void define (const BoxArray&            bxs,
                 const Vector<IndexType>    ixtypes,
                 const DistributionMapping& dm,
                 int                        nvar,
                 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FAB>&     factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif
    //! get a reference to set element of index \p iSet
    MF& operator[] (int iSet) noexcept { return *m_mf_set[iSet]; }
    //! get a const reference to set element of index \p iSet
    const MF& operator[] (int iSet) const noexcept { return *m_mf_set[iSet]; }

    //! get a reference to element of index \p iSet
    MF& getElem (int iSet) noexcept { return *m_mf_set[iSet]; }
    //! get a const reference to element of index iSet
    const MF& getElem (int iSet) const noexcept { return *m_mf_set[iSet]; }

    [[nodiscard]] IntVect nGrowVect (int iSet = 0) const noexcept { return getElem(iSet).nGrowVect(); }

    //! Return the grow factor that defines the region of definition for element of index `iSet`
    [[nodiscard]] int nGrow (int direction = 0, int iSet = 0) const noexcept { return getElem(iSet).nGrow(direction); }

    //! Return number of variables (aka components) associated with each point for `iSet` element
    [[nodiscard]] int nComp (int iSet = 0) const noexcept { return getElem(iSet).nComp(); }

    //Return whether every element of this FabArraySet is empty
    [[nodiscard]] bool empty () const noexcept;

    /**
    * \brief Return a constant reference to the nodal BoxArray that defines the
    * valid region associated with this FabArraySet.
    */
    [[nodiscard]] const BoxArray& boxArray (int iSet = 0) const noexcept { return getElem(iSet).boxArray(); }

    //! Return constant reference to associated DistributionMapping.
    [[nodiscard]] const DistributionMapping& DistributionMap (int iSet = 0) const noexcept { return getElem(iSet).DistributionMap(); }

    //! number of MF elements in the set
    [[nodiscard]] std::size_t nSet () const noexcept { return m_mf_set.size(); }

    const FabFactory<FAB>& Factory (int iSet = 0) const noexcept { return getElem(iSet).Factory(); }

    //! Used in MLLinOpT<MF>::make
    void clear ();

    /**
     * \brief Perform local copy of FabArraySet data.
     *
     * \param src source FabArraySet
     * \param scomp starting component of source
     * \param dcomp starting component of this FabArraySet
     * \param ncomp number of components
     * \param nghost number of ghost cells
     */
    template <typename SFA, typename DFA = MF,
              std::enable_if_t<std::conjunction_v<
                  IsFabArray<DFA>, IsFabArray<SFA>,
                  std::is_convertible<typename SFA::value_type,
                                      typename DFA::value_type>>, int> = 0>
    void LocalCopy (const FabArraySet<SFA>& src, int scomp, int dcomp, int ncomp,
                    IntVect const& nghost);

    /**
     * \brief Perform local addition of FabArraySet data.
     *
     * The two FabArraySets must have the same BoxArray and
     * DistributionMapping.
     *
     * \param src    source FabArraySet
     * \param scomp  starting component of source
     * \param dcomp  starting component of this FabArraySet
     * \param ncomp  number of components
     * \param nghost number of ghost cells
     */
    // template <class F=FA, std::enable_if_t<IsFabArray<F>::value,int> = 0>
    void LocalAdd (FabArraySet<MF> const& src, int scomp, int dcomp, int ncomp,
                   IntVect const& nghost);

    bool isAllRegular () const noexcept;

    //! Set all components in the entire region of each FAB to val.
    void setVal (value_type val);

    //! Set ncomp values in the boundary region, starting at start_comp to val.
    // template <class F=FA, typename std::enable_if<IsFabArray<F>::value,int>::type = 0>
    void setBndry (value_type val, int strt_comp, int ncomp);

    /**
     * \brief Return the largest infinity norm of the MF's in the FabArraySet
     *
     * \param comp           starting component
     * \param ncomp          number of components
     * \param nghost         number of ghost cells
     * \param local          If true, MPI communciation is skipped.
     * \param ignore_covered ignore covered cells. Only relevant for cell-centered EB data.
     */
    // template <typename F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    // typename F::value_type
    value_type norminf (int comp, int ncomp, IntVect const& nghost, bool local = false,
             [[maybe_unused]] bool ignore_covered = false) const;

    void ParallelCopy (const FabArraySet& src,
                       const Periodicity& period = Periodicity::NonPeriodic(),
                       FabArrayBase::CpOp op = FabArrayBase::COPY);

    void ParallelCopy (const FabArraySet&   src,
                       int                  scomp,
                       int                  dcomp,
                       int                  ncomp,
                       const IntVect&       snghost,
                       const IntVect&       dnghost,
                       const Periodicity&   period = Periodicity::NonPeriodic(),
                       FabArrayBase::CpOp   op = FabArrayBase::COPY,
                       const FabArrayBase::CPC* a_cpc = nullptr);

    /**
     * \brief y += a*x
     *
     * \param y      FabArraySet y
     * \param a      scalar a
     * \param x      FabArraySet x
     * \param xcomp  starting component of x
     * \param ycomp  starting component of y
     * \param ncomp  number of components
     * \param nghost number of ghost cells
     */
    // template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    static void Saxpy (FabArraySet& y, value_type a, FabArraySet<MF> const& x,
                       int xcomp, int ycomp, int ncomp, IntVect const& nghost);

    /**
     * \brief y = x + a*y
     *
     * \param y      FabArraySet y
     * \param a      scalar a
     * \param x      FabArraySet x
     * \param xcomp  starting component of x
     * \param ycomp  starting component of y
     * \param ncomp  number of components
     * \param nghost number of ghost cells
     */
    // template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    static void Xpay (FabArraySet& y, value_type a, FabArraySet<MF> const& x,
                      int xcomp, int ycomp, int ncomp, IntVect const& nghost);

    /**
     * \brief dst = a*x + b*y
     *
     * \param dst     destination FabArray
     * \param a       scalar a
     * \param x       FabArraySet x
     * \param xcomp   starting component of x
     * \param b       scalar b
     * \param y       FabArraySet y
     * \param ycomp   starting component of y
     * \param dstcomp starting component of destination
     * \param numcomp number of components
     * \param nghost  number of ghost cells
     */
    // template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    static void LinComb (FabArraySet& dst,
                         value_type a, const FabArraySet& x, int xcomp,
                         value_type b, const FabArraySet& y, int ycomp,
                         int dstcomp, int numcomp, const IntVect& nghost);

    /**
     * \brief Negates the value of each cell in the valid region of
     * the MultiFab.  The value of nghost specifies the number of
     * cells in the boundary region that should be modified.
     */
    template <typename A=MF, std::enable_if_t<IsMultiFab<A>::value||IsiMultiFab<A>::value,int> = 0>
    void negate (int nghost = 0);

protected:

    Vector<std::unique_ptr<MF>> m_mf_set;
    // std::unique_ptr<MF>* m_mf_set_2;
    // MF* m_mf_set_3;
    // Array<MF*,10> m_mf_set_4;
    // Array<MF,10> m_mf_set_4;


};

template <class MF>
FabArraySet<MF>::FabArraySet () noexcept
{
    // for (int iSet = 0; iSet < nSet(); ++iSet) {
    //     m_mf_set[iSet] = std::make_unique<MF>();
    // }
}

template <class MF>
FabArraySet<MF>::FabArraySet (const BoxArray& bxs, const Vector<IndexType> ixTypes, const DistributionMapping& dm,
                              int nvar, const IntVect& ngrow, const MFInfo& info, const FabFactory<FAB>& factory)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(ixTypes.size());
    // m_mf_set_2 = new std::unique_ptr<MF>();
    // m_mf_set_3 = new MF[ixTypes.size()];
    // m_mf_set_4 = (MF*)std::malloc(sizeof(MF)*ixTypes.size());
    for (int iSet = 0; iSet < ixTypes.size(); ++iSet) {
        m_mf_set[iSet] = std::make_unique<MF>(amrex::convert(bxs,ixTypes[iSet]),dm,nvar,ngrow,info,factory);
    }
}

template <class MF>
FabArraySet<MF>::FabArraySet (const FabArraySet<MF>& rhs, MakeType maketype, int scomp, int ncomp)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int iSet = 0; iSet < rhs.nSet(); ++iSet) {
        m_mf_set[iSet] = std::make_unique<MF>(rhs[iSet], maketype, scomp, ncomp);
    }
}

template <class MF>
FabArraySet<MF>::FabArraySet (FabArraySet<MF>&& rhs) noexcept
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int iSet = 0; iSet < rhs.nSet(); ++iSet) {
        m_mf_set[iSet] = std::make_unique<MF>(std::move(rhs[iSet]));
    }
}

template <class MF>
FabArraySet<MF>&
FabArraySet<MF>::operator= (FabArraySet<MF>&& rhs) noexcept
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int iSet = 0; iSet < rhs.nSet(); ++iSet) {
        *m_mf_set[iSet] = std::move(rhs[iSet]);
        // getElem(iSet).operator=(std::move(rhs[iSet]));
        // getElem(iSet) = std::move(rhs[iSet]);
    }
    return *this;
}

template <class MF>
void
FabArraySet<MF>::define (const BoxArray& bxs, const Vector<IndexType> ixTypes, const DistributionMapping& dm,
                         int nvar, const IntVect& ngrow, const MFInfo& info, const FabFactory<FAB>& factory)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(ixTypes.size());
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).define(amrex::convert(bxs,ixTypes[iSet]),dm,nvar,ngrow,info,factory);
    }
}

template <class MF>
bool
FabArraySet<MF>::empty () const noexcept
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        if (!getElem(iSet).empty()) {
            return false;
        }
    }
    return true;
}

template <class MF>
void
FabArraySet<MF>::setVal (value_type val)
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).setVal(val);
    }
}

template <class MF>
void
FabArraySet<MF>::setBndry (value_type val, int strt_comp, int ncomp)
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).setBndry(val,strt_comp,ncomp);
    }
}

template <class MF>
FabArraySet<MF>::value_type
FabArraySet<MF>::norminf (int comp, int ncomp, IntVect const& nghost, bool local,
                          [[maybe_unused]] bool ignore_covered) const
{
    value_type result = value_type(0);
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        result = std::max(result, getElem(iSet).norminf(comp,ncomp,nghost,local,ignore_covered));
    }
    return result;
}

template <class MF>
void
FabArraySet<MF>::clear()
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).clear();
    }
}

template <class MF>
template <typename SFA, typename DFA,
              std::enable_if_t<std::conjunction_v<
                  IsFabArray<DFA>, IsFabArray<SFA>,
                  std::is_convertible<typename SFA::value_type,
                                      typename DFA::value_type>>, int>>
void
FabArraySet<MF>::LocalCopy (const FabArraySet<SFA>& src, int scomp, int dcomp, int ncomp,
                            IntVect const& nghost)
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).LocalCopy(src[iSet],scomp,dcomp,ncomp,nghost);
    }
}

template <class MF>
void
FabArraySet<MF>::LocalAdd (FabArraySet<MF> const& src, int scomp, int dcomp, int ncomp,
                       IntVect const& nghost)
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).LocalAdd(src[iSet],scomp,dcomp,ncomp,nghost);
    }
}

template <class MF>
bool
FabArraySet<MF>::isAllRegular() const noexcept
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        if (!getElem(iSet).isAllRegular()) {
            return false;
        }
    }
    return true;
}

template <class MF>
void
FabArraySet<MF>::ParallelCopy (const FabArraySet& src, const Periodicity& period, FabArrayBase::CpOp op)
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).ParallelCopy(src[iSet],period,op);
    }
}

template <class MF>
void
FabArraySet<MF>::ParallelCopy (const FabArraySet<MF>& src, int scomp, int dcomp, int ncomp,
                               const IntVect& snghost, const IntVect& dnghost,
                               const Periodicity& period, FabArrayBase::CpOp op,
                               const FabArrayBase::CPC* a_cpc)
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).ParallelCopy(src[iSet],scomp,dcomp,ncomp,snghost,dnghost,period,op,a_cpc);
    }
}

template <class MF>
void
FabArraySet<MF>::Saxpy (FabArraySet<MF>& y, value_type a, const FabArraySet<MF>& x,
                        int xcomp, int ycomp, int ncomp, const IntVect& nghost)
{
    for (int iSet = 0; iSet < y.nSet(); ++iSet) {
        MF::Saxpy(y[iSet],a,x[iSet],xcomp,ycomp,ncomp,nghost);
    }
}

template <class MF>
void
FabArraySet<MF>::Xpay (FabArraySet<MF>& y, value_type a, const FabArraySet<MF>& x,
                       int xcomp, int ycomp, int ncomp, const IntVect& nghost)
{
    for (int iSet = 0; iSet < y.nSet(); ++iSet) {
        MF::Xpay(y[iSet],a,x[iSet],xcomp,ycomp,ncomp,nghost);
    }
}

/**
 * \brief dst = a*x + b*y
 *
 * \param dst     destination FabArray
 * \param a       scalar a
 * \param x       FabArraySet x
 * \param xcomp   starting component of x
 * \param b       scalar b
 * \param y       FabArraySet y
 * \param ycomp   starting component of y
 * \param dstcomp starting component of destination
 * \param numcomp number of components
 * \param nghost  number of ghost cells
 */
template <class MF>
void
FabArraySet<MF>::LinComb (FabArraySet<MF>& dst,
                          value_type a, const FabArraySet<MF>& x, int xcomp,
                          value_type b, const FabArraySet<MF>& y, int ycomp,
                          int dstcomp, int numcomp, const IntVect& nghost)
{
    for (int iSet = 0; iSet < y.nSet(); ++iSet) {
        MF::LinComb(dst[iSet],a,x[iSet],xcomp,b,y[iSet],ycomp,dstcomp,numcomp,nghost);
    }
}

template <class MF>
template <typename A, std::enable_if_t<IsMultiFab<A>::value||IsiMultiFab<A>::value,int>>
void
FabArraySet<MF>::negate (int nghost)
{
    for (int iSet = 0; iSet < nSet(); ++iSet) {
        getElem(iSet).negate(nghost);
    }
}

extern template class FabArraySet<MultiFab>;

using MultiFabSet = FabArraySet<MultiFab>;

};

#endif /* AMREX_FABARRAYSET_H_ */