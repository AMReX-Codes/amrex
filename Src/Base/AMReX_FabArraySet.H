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
class FabArraySet : BaseFabArraySet
{

public:
    //! debug variable
    static const bool m_ncompRefersToMFComps = true;
    //! debug variable
    static const unsigned long encodedIxTypes = 40UL;
    //! debug variable
    using MF = typename amrex::MultiFab;
    
    using value_type = typename MF::value_type;
    using fab_type = typename MF::fab_type;
    
    using FAB = typename MF::fab_type;

    //! Constructs an empty FabArray<FAB>.
    FabArraySet () noexcept;
    //! The destructor -- deletes all FABs in the array.
    ~FabArraySet ();

    FabArraySet (FabArraySet&& rhs) noexcept = default;
    FabArraySet (const FabArraySet& rhs) = default;
    FabArraySet& operator= (const FabArraySet& rhs) = default;
    FabArraySet& operator= (FabArraySet&& rhs) = default;

    MF& operator[] (int i) noexcept { return *m_mfptr_set[i]; }
    const MF& operator[] (int i) const noexcept { return *m_mfptr_set[i]; }

    MF& getElem (int i) noexcept { return *m_mfptr_set[i]; }
    const MF& getElem (int i) const noexcept { return *m_mfptr_set[i]; }
    
    MF* getElemPtr (int i) noexcept { return m_mfptr_set[i]; }
    const MF* getElemPtr (int i) const noexcept { return m_mfptr_set[i]; }

    void setElem (int i, MF& mf) { this->m_mfptr_set[i] = &mf; }
    void setElem (int i, MF* mfptr) { this->m_mfptr_set[i] = mfptr; }
    
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
     * \param use_nGrowVect Use nGrowVect from each FabArray instead of nghost
     */
    void LocalCopy (const FabArraySet& src, int scomp, int dcomp, int ncomp,
                    IntVect const& nghost, bool use_nGrowVect = false);

    bool isAllRegular () const noexcept;

    //! Set all components in the entire region of each FAB to val.
    void setVal (value_type val);

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

protected:

    Vector<MF*>                      m_mfptr_set;
    static const Vector<IndexType>   m_ixtype_set;
    static const int                 m_nSet;
    std::unique_ptr<FabFactory<FAB>> m_factory;
};

// template <class MF, unsigned long encodedIxTypes>
const Vector<IndexType> FabArraySet::m_ixtype_set = IndexTypeSet::ConvertToIxTypes(FabArraySet::encodedIxTypes);

// template <class MF, unsigned long encodedIxTypes>
const int FabArraySet::m_nSet = IndexTypeSet::ConvertToIxTypes(FabArraySet::encodedIxTypes).size();

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
    if (m_ncompRefersToMFComps) {
        for (MF* mfptr : m_mfptr_set) {
            result = std::max(result, mfptr->norminf(comp, ncomp, nghost, local, ignore_covered));
        }
    } else {
        for (int i = comp; i < comp + ncomp; ++i) {
            const auto& mf = getElem(i);
            result = std::max(result, mf.norminf(0, mf.nComp(), nghost, local, ignore_covered));
        }
    }
    return result;
}

void
FabArraySet::clear()
{
    for (MF* mfptr : m_mfptr_set) {
        mfptr->clear();
    }
    m_mfptr_set.clear();
}

void
FabArraySet::LocalCopy (const FabArraySet& src, int scomp, int dcomp, int ncomp,
                        IntVect const& nghost, bool use_nGrowVect)
{
    for (int i = 0; i < nSet(); ++i) {
        (*this)[i].LocalCopy(src[i], scomp, dcomp, ncomp, nghost);
    }

    // Define common variables:
    // int i, n;
    // IntVect iv;

    // Two possible ways to write this:
    // 1) Make scomp, dcomp, and ncomp refer to each component of the MultiFabSet:
    // BL_ASSERT(scomp == dcomp);
    // n = ncomp < 1 ? src.nSet() : ncomp;
    // for (i = scomp; i < n + scomp; ++i) {
    //     MF& dstmf = (*this)[i];
    //     const MF& srcmf = src[i];
    //     iv = use_nGrowVect ? srcmf.nGrowVect() : nghost;
    //     dstmf.LocalCopy(srcmf, 0, 0, srcmf.nComp(), iv);
    // }
    // 2) Make scomp, dcomp, and ncomp refer to each component of the underlying MultiFabs:
    // for (i = 0; i < src.nSet(); ++i) {
    //     const MF& srcmf = src[i];
    //     n = ncomp < 1 ? srcmf.nComp() : ncomp;
    //     iv = use_nGrowVect ? srcmf.nGrowVect() : nghost;
    //     getElem(i).LocalCopy(srcmf, scomp, dcomp, n, iv);
    // }
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