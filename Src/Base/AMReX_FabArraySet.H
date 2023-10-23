#include <AMReX_FabArray.H>
#include <AMReX_IndexType.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

namespace amrex {

/**
 * \brief FabArraySet
 *
 * \tparam MF FabArray or derived class
 */
// template <class MF>
class FabArraySet
{
public:
    //! debug variable
    using MF = typename amrex::MultiFab;
    // using MF = typename amrex::FabArray<FArrayBox>;

    // using value_type = typename MF::value_type;
    // using fab_type = typename MF::fab_type;
    // using FAB = typename MF::fab_type;

    using value_type = typename amrex::Real;
    using fab_type = typename amrex::FArrayBox;
    using FAB = typename amrex::FArrayBox;

    //! Constructs an empty FabArraySet.
    FabArraySet () noexcept;

    // /**
    //  * \brief Construct an empty FabArraySet that has a default Arena.  If
    //  * `define` is called later with a nulltpr as MFInfo's arena, the
    //  * default Arena `a` will be used.  If the arena in MFInfo is not a
    //  * nullptr, the MFInfo's arena will be used.
    //  */
    // explicit FabArraySet (Arena* a) noexcept;

    /**
     * \brief Construct a FabArraySet with a valid region defined by bxs
     * and a region of definition defined by IndexType's \p ixtypes,
     * grow factor \p ngrow, and the number of components per MF \p nvar.
     */

    /**
     * \brief Constructs a new FabArraySet<MF>
     *
     * The size of the FArrayBox is given by the Box modified by \p ixtypes
     * and grown by \p `ngrow`. The number of components for each MF is given
     * by \p `ncomp` . If info is set to
     * not allocating memory, then no FArrayBoxes are allocated at
     * this time but can be defined later.
     *
     * \param bxs a valid region
     * \param ixTypes IndexType for each MF in the set. ixTypes.size() specifies the number of MFs in the set
     * \param dm a DistribuionMapping
     * \param ncomp number of components for each MF (different from number of MFs in the set)
     * \param ngrow number of cells the region grows
     * \param info MFInfo
     *
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
<<<<<<< HEAD
=======
                 const DistributionMapping& dm,
                 int                        nvar,
                 int                        ngrow,
#ifdef AMREX_STRICT_MODE
                 const MFInfo&              info,
                 const FabFactory<FAB>&     factory);
#else
                 const MFInfo&              info = MFInfo(),
                 const FabFactory<FAB>&     factory = DefaultFabFactory<FAB>());
#endif

    void define (const BoxArray&            bxs,
                 const Vector<IndexType>    ixtypes,
>>>>>>> 6c1f74744aee54f2bdd4b38491c83ee3ff7d39c8
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
    //! get a reference to element of index i
    MF& operator[] (int i) noexcept { return *m_mf_set[i]; }
    //! get a const reference to element of index i
    const MF& operator[] (int i) const noexcept { return *m_mf_set[i]; }

    //! get a reference to element of index i
    MF& getElem (int i) noexcept { return *m_mf_set[i]; }
    //! get a const reference to element of index i
    const MF& getElem (int i) const noexcept { return *m_mf_set[i]; }

    [[nodiscard]] IntVect nGrowVect (int i = 0) const noexcept { return getElem(i).nGrowVect(); }

    //! Return the grow factor that defines the region of definition for i-th element
    [[nodiscard]] int nGrow (int direction = 0, int i = 0) const noexcept { return getElem(i).nGrow(direction); }

    //! Return number of variables (aka components) associated with each point for i-th element
    [[nodiscard]] int nComp (int i = 0) const noexcept { return getElem(i).nComp(); }

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
     * \brief Return the largest infinity norm of the MultiFabs
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

FabArraySet::FabArraySet () noexcept
{
    // for (int i = 0; i < nSet(); ++i) {
    //     m_mf_set[i] = std::make_unique<MF>();
    // }
}

/**
* \brief Construct a FabArray<FAB> with a valid region defined by bxs
* and a region of definition defined by the grow factor ngrow
* and the number of components nvar.
*/
FabArraySet::FabArraySet (const BoxArray& bxs, const Vector<IndexType> ixTypes, const DistributionMapping& dm,
                          int nvar, int ngrow, const MFInfo& info, const FabFactory<FAB>& factory)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(ixTypes.size());
    // m_mf_set_2 = new std::unique_ptr<MF>();
    // m_mf_set_3 = new MF[ixTypes.size()];
    // m_mf_set_4 = (MF*)std::malloc(sizeof(MF)*ixTypes.size());

    for (int i = 0; i < ixTypes.size(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(amrex::convert(bxs,ixTypes[i]),dm,nvar,ngrow,info,factory);
    }
}


FabArraySet::FabArraySet (const BoxArray& bxs, const Vector<IndexType> ixTypes, const DistributionMapping& dm,
                          int nvar, const IntVect& ngrow, const MFInfo& info, const FabFactory<FAB>& factory)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(ixTypes.size());
    for (int i = 0; i < ixTypes.size(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(amrex::convert(bxs,ixTypes[i]),dm,nvar,ngrow,info,factory);
    }
}

FabArraySet::FabArraySet (const FabArraySet& rhs, MakeType maketype, int scomp, int ncomp)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int i = 0; i < rhs.nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(rhs[i], maketype, scomp, ncomp);
    }
}

FabArraySet::FabArraySet (FabArraySet&& rhs) noexcept
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int i = 0; i < rhs.nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(std::move(rhs[i]));
    }
}

FabArraySet&
FabArraySet::operator= (FabArraySet&& rhs) noexcept
{
    m_mf_set = Vector<std::unique_ptr<MF>>(rhs.nSet());
    for (int i = 0; i < rhs.nSet(); ++i) {
        *m_mf_set[i] = std::move(rhs[i]);
        // getElem(i).operator=(std::move(rhs[i]));
        getElem(i) = std::move(rhs[i]);
    }
    return *this;
}

void
FabArraySet::define (const BoxArray& bxs, const Vector<IndexType> ixTypes, const DistributionMapping& dm,
                     int nvar, int ngrow, const MFInfo& info, const FabFactory<FAB>& factory)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(ixTypes.size());
    for (int i = 0; i < ixTypes.size(); ++i) {
        getElem(i).define(amrex::convert(bxs,ixTypes[i]),dm,nvar,ngrow,info,factory);
    }
}

void
FabArraySet::define (const BoxArray& bxs, const Vector<IndexType> ixTypes, const DistributionMapping& dm,
                     int nvar, const IntVect& ngrow, const MFInfo& info, const FabFactory<FAB>& factory)
{
    m_mf_set = Vector<std::unique_ptr<MF>>(ixTypes.size());
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).define(amrex::convert(bxs,ixTypes[i]),dm,nvar,ngrow,info,factory);
    }
}

bool
FabArraySet::empty () const noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        if (!getElem(i).empty()) {
            return false;
        }
    }
    return true;
}

void
FabArraySet::setVal (value_type val)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).setVal(val);
    }
}

void
FabArraySet::setBndry (value_type val, int strt_comp, int ncomp)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).setBndry(val,strt_comp,ncomp);
    }
}

FabArraySet::value_type
FabArraySet::norminf (int comp, int ncomp, IntVect const& nghost, bool local,
                      [[maybe_unused]] bool ignore_covered) const
{
    value_type result = value_type(0);
    for (int i = 0; i < nSet(); ++i) {
        result = std::max(result, getElem(i).norminf(comp,ncomp,nghost,local,ignore_covered));
    }
    return result;
}

void
FabArraySet::clear()
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).clear();
    }
}

void
FabArraySet::LocalCopy (const FabArraySet& src, int scomp, int dcomp, int ncomp,
                        IntVect const& nghost)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).LocalCopy(src[i],scomp,dcomp,ncomp,nghost);
    }
}

void
FabArraySet::LocalAdd (FabArraySet const& src, int scomp, int dcomp, int ncomp,
                       IntVect const& nghost)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).LocalAdd(src[i],scomp,dcomp,ncomp,nghost);
    }
}

bool
FabArraySet::isAllRegular() const noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        if (!getElem(i).isAllRegular()) {
            return false;
        }
    }
    return true;
}

void
FabArraySet::ParallelCopy (const FabArraySet& src, const Periodicity& period, FabArrayBase::CpOp op)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).ParallelCopy(src[i],period,op);
    }
}

void
FabArraySet::ParallelCopy (const FabArraySet& src, int scomp, int dcomp, int ncomp,
                           const IntVect& snghost, const IntVect& dnghost,
                           const Periodicity& period, FabArrayBase::CpOp op,
                           const FabArrayBase::CPC* a_cpc)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).ParallelCopy(src[i],scomp,dcomp,ncomp,snghost,dnghost,period,op,a_cpc);
    }
}

// template <class FA>
// template <class F=FA, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
void
FabArraySet::Saxpy (FabArraySet& y, value_type a, FabArraySet const& x,
                    int xcomp, int ycomp, int ncomp, IntVect const& nghost)
{
    for (int i = 0; i < y.nSet(); ++i) {
        MF::Saxpy(y[i],a,x[i],xcomp,ycomp,ncomp,nghost);
    }
}

// template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
void
FabArraySet::Xpay (FabArraySet& y, value_type a, FabArraySet const& x,
                   int xcomp, int ycomp, int ncomp, IntVect const& nghost)
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
// template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
void
FabArraySet::LinComb (FabArraySet& dst,
                     value_type a, const FabArraySet& x, int xcomp,
                     value_type b, const FabArraySet& y, int ycomp,
                     int dstcomp, int numcomp, const IntVect& nghost)
{
    for (int i = 0; i < y.nSet(); ++i) {
        MF::LinComb(dst[i],a,x[i],xcomp,b,y[i],ycomp,dstcomp,numcomp,nghost);
    }
}

template <typename A, std::enable_if_t<IsMultiFab<A>::value||IsiMultiFab<A>::value,int>>
void
FabArraySet::negate (int nghost = 0)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).negate(nghost);
    }
}

};