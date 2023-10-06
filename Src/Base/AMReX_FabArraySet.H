#include <AMReX_FabArray.H>
#include <AMReX_BaseFabArraySet.H>
#include <AMReX_IndexType.H>

namespace amrex {

static constexpr std::size_t
amrex::GetSetSizeFromEncodedLong (unsigned long a_EncodedLong)
{
    return 1UL + (a_EncodedLong & (IndexTypeSet::MAX_SIZE - 1UL));
}

/**
 * \brief FabArraySet
 * 
 * \tparam MF FabArray or derived class
 * \tparam encodedIxTypes 
 */
// template <class MF, unsigned long encodedIxTypes>
// template <class MF, >
// template <unsigned long encodedIxTypes>
class FabArraySet
{
public:
    //! debug variable
    static constexpr unsigned long encodedIxTypes = 40UL;
    //! debug variable
    // using MF = typename amrex::MultiFab;
    using MF = typename amrex::FabArray<FArrayBox>;

protected:

    static constexpr std::size_t m_nSet = amrex::GetSetSizeFromEncodedLong(encodedIxTypes);
    static constexpr Array<IndexType, m_nSet> m_ixtype_set;

public:

    using value_type = typename MF::value_type;
    using fab_type = typename MF::fab_type;
    using FAB = typename MF::fab_type;

    //! Constructs an empty FabArraySet.
    FabArraySet () noexcept;
    
    /**
     * \brief Construct an empty FabArray<FAB> that has a default Arena.  If
     * `define` is called later with a nulltpr as MFInfo's arena, the
     * default Arena `a` will be used.  If the arena in MFInfo is not a
     * nullptr, the MFInfo's arena will be used.
     */
    explicit FabArraySet (Arena* a) noexcept;

    /**
    * \brief Construct a FabArray<FAB> with a valid region defined by bxs
    * and a region of definition defined by the grow factor ngrow
    * and the number of components nvar.
    */
    FabArraySet (const BoxArray&            bxs,
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

    FabArraySet (const BoxArray&            bxs,
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
    ~FabArraySet () = default;

    FabArraySet (FabArraySet&& rhs) noexcept;
    FabArraySet& operator= (FabArraySet&& rhs) noexcept;
    
    FabArraySet (const FabArraySet& rhs) = delete;
    FabArraySet& operator= (const FabArray<FAB>& rhs) = delete;

    /**
    * \brief Define this FabArraySet identically to that performed by
    * the constructor having an analogous function signature.
    * This is only valid if this FabArraySet was defined using
    * the default constructor.
    */
    void define (const BoxArray&            bxs,
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

    void setElem (int i, std::unique_ptr<MF>&& a_mf_ptr) { m_mf_set[i] = std::move(a_mf_ptr); };

    static IndexType getElemType (int i) { return m_ixtype_set[i]; }

    //! make a BoxArray with the appropriate index type for the i-th element
    static BoxArray MakeElemBoxArray (const BoxArray bxs, int i) { return amrex::convert(bxs, getElemType(i)); }
    
    [[nodiscard]] IntVect nGrowVect (int i = 0) const noexcept { return getElem(i).nGrowVect(); }

    //! Return the grow factor that defines the region of definition.
    [[nodiscard]] int nGrow (int direction = 0, int i = 0) const noexcept { return getElem(i).nGrow(direction); }

    //! Return number of variables (aka components) associated with each point.
    //! Need to decide if nComp reutrns the number of elements in the set or the number of components in each MF
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
    
    // MF* getElemPtr (int i) noexcept { return m_mf_set[i]; }
    // const MF* getElemPtr (int i) const noexcept { return &m_mf_set[i]; }

    // void setElem (int i, MF& mf) { m_mf_set[i] = mf; }
    // void setElem (int i, MF* mfptr) { m_mf_set[i] = *mfptr; }

    //! number of MF elements in the set
    [[nodiscard]] std::size_t nSet () const noexcept { return m_nSet; }

    //! Return the set of index type set
    // [[nodiscard]] const Vector<IndexType>& ixTypeSet () const noexcept { return m_ixtype_set; }
    [[nodiscard]] const Array<IndexType,m_nSet>& ixTypeSet () const noexcept { return m_ixtype_set; }

    const FabFactory<FAB>& Factory () const noexcept { return *m_factory; }

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
    // template <typename SFAB, typename DFAB = FAB,
    //           std::enable_if_t<std::conjunction_v<
    //               IsBaseFab<DFAB>, IsBaseFab<SFAB>,
    //               std::is_convertible<typename SFAB::value_type,
    //                                   typename DFAB::value_type>>, int> = 0>
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
    // template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
    void LocalAdd (FabArraySet const& src, int scomp, int dcomp, int ncomp,
                   IntVect const& nghost);

    bool isAllRegular () const noexcept;

    //! Set all components in the entire region of each FAB to val.
    void setVal (value_type val);

    //! Set ncomp values in the boundary region, starting at start_comp to val.
    // template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
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

protected:
    // MF* m_mf_set;
    // Array<MF,m_nSet> m_mf_set;
    Array<std::unique_ptr<MF>,m_nSet> m_mf_set;
    // bool m_NeedToDelete = false;
    std::unique_ptr<FabFactory<FAB>> m_factory;

};

// template <class MF, unsigned long encodedIxTypes>
// constexpr Vector<IndexType> FabArraySet::m_ixtype_set = IndexTypeSet::ConvertToIxTypes(FabArraySet::encodedIxTypes);

// template <class MF, unsigned long encodedIxTypes>
// constexpr int FabArraySet::m_nSet = IndexTypeSet::ConvertToIxTypes(FabArraySet::encodedIxTypes).size();


FabArraySet::FabArraySet () noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>();
    }
}

/**
* \brief Construct a FabArray<FAB> with a valid region defined by bxs
* and a region of definition defined by the grow factor ngrow
* and the number of components nvar.
*/
FabArraySet::FabArraySet (const BoxArray& bxs, const DistributionMapping& dm,
                          int nvar, int ngrow, const MFInfo& info,
                          const FabFactory<FAB>& factory)
{
    for (int i = 0; i < nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(MakeElemBoxArray(bxs,i),dm,nvar,ngrow,info,factory);
    }
}


FabArraySet::FabArraySet (const BoxArray& bxs, const DistributionMapping& dm,
                          int nvar, const IntVect& ngrow, const MFInfo& info,
                          const FabFactory<FAB>& factory)
{
    for (int i = 0; i < nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(MakeElemBoxArray(bxs,i),dm,nvar,ngrow,info,factory);
    }
}

FabArraySet::FabArraySet (const FabArraySet& rhs, MakeType maketype, int scomp, int ncomp)
{
    for (int i = 0; i < nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(rhs[i], maketype, scomp, ncomp);
    }
}

FabArraySet::FabArraySet (FabArraySet&& rhs) noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        m_mf_set[i] = std::make_unique<MF>(std::move(rhs[i]));
    }
}

FabArraySet&
FabArraySet::operator= (FabArraySet&& rhs) noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).operator=(std::move(rhs[i]));
    }
    return *this;
}

void 
FabArraySet::define (const BoxArray& bxs, const DistributionMapping& dm, 
                     int nvar, int ngrow, const MFInfo& info, 
                     const FabFactory<FAB>& factory)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).define(MakeElemBoxArray(bxs,i),dm,nvar,ngrow,info,factory);
    }
}

void
FabArraySet::define (const BoxArray& bxs, const DistributionMapping& dm,
                     int nvar, const IntVect& ngrow, const MFInfo& info,
                     const FabFactory<FAB>& factory)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).define(MakeElemBoxArray(bxs,i),dm,nvar,ngrow,info,factory);
    }
}

bool
FabArraySet::empty () const noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        if (!getElem(i).empty()) {
            return false;
        }
    };
    return true;
}

void
FabArraySet::setVal (value_type val)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).setVal(val);
    };
}

void
FabArraySet::setBndry (value_type val, int strt_comp, int ncomp)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).setBndry(val,strt_comp,ncomp);
    };
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
FabArraySet::ParallelCopy (const FabArraySet& src, int scomp, int dcomp, int ncomp,
                           const IntVect& snghost, const IntVect& dnghost,
                           const Periodicity& period, FabArrayBase::CpOp op,
                           const FabArrayBase::CPC* a_cpc)
{
    for (int i = 0; i < nSet(); ++i) {
        getElem(i).ParallelCopy(src[i],scomp,dcomp,ncomp,snghost,dnghost,period,op,a_cpc);
    }
}

// template <class F=FAB, std::enable_if_t<IsBaseFab<F>::value,int> = 0>
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

};