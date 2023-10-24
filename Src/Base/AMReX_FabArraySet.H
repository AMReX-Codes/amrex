#include <AMReX_FabArray.H>
#include <AMReX_IndexType.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

namespace amrex {

/**
 * \brief A set of FabArray objects with distinct IndexType's
 *
 * Intended for use as a single object to use as input/output in a linear solver
 * class that inherits from MLLinOpT<FabArraySet<MF>>.
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

    // /**
    //  * \brief Construct an empty FabArraySet that has a default Arena.  If
    //  * `define` is called later with a nullptr as MFInfo's arena, the
    //  * default Arena `a` will be used.  If the arena in MFInfo is not a
    //  * nullptr, the MFInfo's arena will be used.
    //  */
    // explicit FabArraySet (Arena* a) noexcept;

    /**
     * \brief Construct a FabArraySet with a valid region defined by \p bxs
     * and a region of definition defined by IndexType's \p ixtypes,
     * grow factor \p ngrow, and the number of components per MF \p nvar.
     */

    /**
     * \brief Constructs a new FabArraySet<MF>
     *
     * The size of the MFs is given by the Box Array \p bxs modified by the IndexType's
     * \p ixtypes and grown by \p ngrow . The number of components for each MF is given
     * by \p ncomp , and the size of \p ixtypes determines number of MFs in the set.
     * If \p info is set to not allocating memory, then no FABs are allocated at this
     * time but can be defined later.
     *
     * \param bxs a valid region
     * \param ixtypes IndexType for each MF. The number of MFs in the set will be the size of \p ixTypes
     * \param dm a DistribuionMapping
     * \param ncomp number of components for each MF (different from number of MFs in the set)
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
    //! get a const reference to element of index i
    const MF& getElem (int iSet) const noexcept { return *m_mf_set[iSet]; }

    [[nodiscard]] IntVect nGrowVect (int i = 0) const noexcept { return getElem(i).nGrowVect(); }

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
    [[nodiscard]] const BoxArray& boxArray (int i = 0) const noexcept { return getElem(i).boxArray(); }

    //! Return constant reference to associated DistributionMapping.
    [[nodiscard]] const DistributionMapping& DistributionMap (int i = 0) const noexcept { return getElem(i).DistributionMap(); }

    //! number of MF elements in the set
    [[nodiscard]] std::size_t nSet () const noexcept { return m_mf_set.size(); }

    // //! Return the set of index types
    // [[nodiscard]] const Vector<IndexType>& ixTypeSet () const noexcept
    // {
    //     Vector<IndexType> ixTypes(nSet());
    //     for (int i = 0; i < nSet(); ++i) {
    //         ixTypes[i] = getElem(i).ixType();
    //     }
    //     return ixTypes;
    // }

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
    // template <typename SFA, typename DFA = FA,
    //           std::enable_if_t<std::conjunction_v<
    //               IsFabArray<DFA>, IsFabArray<SFA>,
    //               std::is_convertible<typename SFA::value_type,
    //                                   typename DFA::value_type>>, int> = 0>
    void LocalCopy (const FabArraySet& src, int scomp, int dcomp, int ncomp,
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
    void LocalAdd (FabArraySet const& src, int scomp, int dcomp, int ncomp,
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
    static void Saxpy (FabArraySet& y, value_type a, FabArraySet const& x,
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
    static void Xpay (FabArraySet& y, value_type a, FabArraySet const& x,
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
    // for (int i = 0; i < nSet(); ++i) {
    //     m_mf_set[i] = std::make_unique<MF>();
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
    for (int i = 0; i < ixTypes.size(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(amrex::convert(bxs,ixTypes[i]),dm,nvar,ngrow,info,factory);
    }
}

template <class MF>
FabArraySet<MF>::FabArraySet (const FabArraySet<MF>& rhs, MakeType maketype, int scomp, int ncomp)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int i = 0; i < rhs.nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(rhs[i], maketype, scomp, ncomp);
    }
}

template <class MF>
FabArraySet<MF>::FabArraySet (FabArraySet<MF>&& rhs) noexcept
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int i = 0; i < rhs.nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(std::move(rhs[i]));
    }
}

template <class MF>
FabArraySet<MF>&
FabArraySet<MF>::operator= (FabArraySet<MF>&& rhs) noexcept
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int i = 0; i < rhs.nSet(); ++i) {
        *m_mf_set[i] = std::move(rhs[i]);
        // getElem(i).operator=(std::move(rhs[i]));
        // getElem(i) = std::move(rhs[i]);
    }
    return *this;
}

template <class MF>
void
FabArraySet<MF>::define (const BoxArray& bxs, const Vector<IndexType> ixTypes, const DistributionMapping& dm,
                         int nvar, const IntVect& ngrow, const MFInfo& info, const FabFactory<FAB>& factory)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(ixTypes.size());
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).define(amrex::convert(bxs,ixTypes[i]),dm,nvar,ngrow,info,factory);
    }
}

template <class MF>
bool
FabArraySet<MF>::empty () const noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        if (!getElem(i).empty()) {
            return false;
        }
    }
    return true;
}

template <class MF>
void
FabArraySet<MF>::setVal (value_type val)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).setVal(val);
    }
}

template <class MF>
void
FabArraySet<MF>::setBndry (value_type val, int strt_comp, int ncomp)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).setBndry(val,strt_comp,ncomp);
    }
}

template <class MF>
FabArraySet<MF>::value_type
FabArraySet<MF>::norminf (int comp, int ncomp, IntVect const& nghost, bool local,
                          [[maybe_unused]] bool ignore_covered) const
{
    value_type result = value_type(0);
    for (int i = 0; i < nSet(); ++i) {
        result = std::max(result, getElem(i).norminf(comp,ncomp,nghost,local,ignore_covered));
    }
    return result;
}

template <class MF>
void
FabArraySet<MF>::clear()
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).clear();
    }
}

template <class MF>
void
FabArraySet<MF>::LocalCopy (const FabArraySet& src, int scomp, int dcomp, int ncomp,
                            IntVect const& nghost)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).LocalCopy(src[i],scomp,dcomp,ncomp,nghost);
    }
}

template <class MF>
void
FabArraySet<MF>::LocalAdd (FabArraySet const& src, int scomp, int dcomp, int ncomp,
                       IntVect const& nghost)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).LocalAdd(src[i],scomp,dcomp,ncomp,nghost);
    }
}

template <class MF>
bool
FabArraySet<MF>::isAllRegular() const noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        if (!getElem(i).isAllRegular()) {
            return false;
        }
    }
    return true;
}

template <class MF>
void
FabArraySet<MF>::ParallelCopy (const FabArraySet& src, const Periodicity& period, FabArrayBase::CpOp op)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).ParallelCopy(src[i],period,op);
    }
}

template <class MF>
void
FabArraySet<MF>::ParallelCopy (const FabArraySet<MF>& src, int scomp, int dcomp, int ncomp,
                               const IntVect& snghost, const IntVect& dnghost,
                               const Periodicity& period, FabArrayBase::CpOp op,
                               const FabArrayBase::CPC* a_cpc)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).ParallelCopy(src[i],scomp,dcomp,ncomp,snghost,dnghost,period,op,a_cpc);
    }
}

template <class MF>
void
FabArraySet<MF>::Saxpy (FabArraySet<MF>& y, value_type a, const FabArraySet<MF>& x,
                        int xcomp, int ycomp, int ncomp, const IntVect& nghost)
{
    for (int i = 0; i < y.nSet(); ++i) {
        MF::Saxpy(y[i],a,x[i],xcomp,ycomp,ncomp,nghost);
    }
}

template <class MF>
void
FabArraySet<MF>::Xpay (FabArraySet<MF>& y, value_type a, const FabArraySet<MF>& x,
                       int xcomp, int ycomp, int ncomp, const IntVect& nghost)
{
    for (int i = 0; i < y.nSet(); ++i) {
        MF::Xpay(y[i],a,x[i],xcomp,ycomp,ncomp,nghost);
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
    for (int i = 0; i < y.nSet(); ++i) {
        MF::LinComb(dst[i],a,x[i],xcomp,b,y[i],ycomp,dstcomp,numcomp,nghost);
    }
}

template <class MF>
template <typename A, std::enable_if_t<IsMultiFab<A>::value||IsiMultiFab<A>::value,int>>
void
FabArraySet<MF>::negate (int nghost)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).negate(nghost);
    }
}

};