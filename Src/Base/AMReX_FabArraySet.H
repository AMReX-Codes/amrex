#include <AMReX_FabArray.H>
#include <AMReX_BaseFabArraySet.H>

namespace amrex {

/**
 * \brief FabArraySet
 * 
 * \tparam MF FabArray or derived class
 * \tparam encodedIxTypes 
 */
// template <class MF, unsigned long encodedIxTypes>
// template <class MF, >
class FabArraySet : public BaseFabArraySet
{

public:
    //! debug variable
    static const unsigned long encodedIxTypes = 40UL;
    //! debug variable
    using MF = typename amrex::MultiFab;
    
    using value_type = typename MF::value_type;
    using fab_type = typename MF::fab_type;
    
    using FAB = typename MF::fab_type;

    //! Constructs an empty FabArraySet.
    FabArraySet () noexcept;
    //! The destructor -- deletes all FabArrays in the set.
    ~FabArraySet ();

    FabArraySet (FabArraySet&& rhs) noexcept = default;
    FabArraySet (const FabArraySet& rhs) = default;
    FabArraySet& operator= (const FabArraySet& rhs) = default;
    FabArraySet& operator= (FabArraySet&& rhs) = default;

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

    MF& operator[] (int i) noexcept { return m_mfptr_set[i]; }
    const MF& operator[] (int i) const noexcept { return m_mfptr_set[i]; }

    MF& getElem (int i) noexcept { return m_mfptr_set[i]; }
    const MF& getElem (int i) const noexcept { return m_mfptr_set[i]; }
    
    MF* getElemPtr (int i) noexcept { return &m_mfptr_set[i]; }
    const MF* getElemPtr (int i) const noexcept { return &m_mfptr_set[i]; }

    void setElem (int i, MF& mf) { m_mfptr_set[i] = &mf; }
    void setElem (int i, MF* mfptr) { m_mfptr_set[i] = mfptr; }
    
    //! number of MF elements in the set
    [[nodiscard]] int nSet () const noexcept { return m_nSet; }

    //! Return the set of index type set
    [[nodiscard]] const Vector<IndexType>& ixTypeSet () const noexcept { return m_ixtype_set; }

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
    template <class F=FAB, typename std::enable_if<IsBaseFab<F>::value,int>::type = 0>
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
    value_type
    norminf (int comp, int ncomp, IntVect const& nghost, bool local = false,
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

    static const Vector<IndexType>   m_ixtype_set;
    static const int                 m_nSet;
    MF*                              m_mfptr_set;
    bool                             m_NeedToDelete = false;
    std::unique_ptr<FabFactory<FAB>> m_factory;
};

// template <class MF, unsigned long encodedIxTypes>
const Vector<IndexType> FabArraySet::m_ixtype_set = IndexTypeSet::ConvertToIxTypes(FabArraySet::encodedIxTypes);

// template <class MF, unsigned long encodedIxTypes>
const int FabArraySet::m_nSet = IndexTypeSet::ConvertToIxTypes(FabArraySet::encodedIxTypes).size();

// template <class FAB>
FabArraySet::FabArraySet () noexcept
{
    for (int i = 0; i < nSet(); ++i) {
        m_mfptr_set.push_back(new MF());
    }
    m_NeedToDelete = true;
}

// template <class FAB>
FabArraySet::~FabArraySet ()
{
    if (m_NeedToDelete) {
        for (MF* mfptr : m_mfptr_set) {
            delete mfptr;
            mfptr = nullptr;
        }
        m_NeedToDelete = false;
    }
}

void
FabArraySet::setVal (value_type val)
{
    for (MF* mfptr : m_mfptr_set) {
        mfptr->setVal(val);
    };
}

FabArraySet::value_type
FabArraySet::norminf (int comp, int ncomp, IntVect const& nghost, bool local,
                      [[maybe_unused]] bool ignore_covered) const
{
    value_type result = value_type(0);
    for (MF* mfptr : m_mfptr_set) {
        result = std::max(result, mfptr->norminf(comp, ncomp, nghost, local, ignore_covered));
    }
    return result;
}

void
FabArraySet::clear()
{
    for (MF* mfptr : m_mfptr_set) {
        mfptr->clear();
    }
}

void
FabArraySet::LocalCopy (const FabArraySet& src, int scomp, int dcomp, int ncomp,
                        IntVect const& nghost)
{
    for (int i = 0; i < nSet(); ++i) {
        (*this)[i].LocalCopy(src[i], scomp, dcomp, ncomp, nghost);
    }
}

bool
FabArraySet::isAllRegular() const noexcept
{
    for (MF* mfptr : m_mfptr_set) {
        if (!(mfptr->isAllRegular())) {
            return false;
        }
    }
    return true;
}

};