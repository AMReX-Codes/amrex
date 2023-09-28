#include <AMReX_FabArray.H>

namespace amrex {

/**
 * \brief BaseMultiFabSet
 * 
 * \tparam MF FabArray or derived class
 * \tparam N Number of FabArray elements in the set
 */
template <class MF>
class BaseMultiFabSet {

public:
    using fab_type = typename MF::fab_type;
    using value_type = typename MF::value_type;

    void initializeSet (int setSize);

    int nSet () const { return m_nSet; }

    

    // BaseMultiFabSet (BaseMultiFabSet<MF>&& rhs) noexcept = default;
    // BaseMultiFabSet (const BaseMultiFabSet<MF>& rhs) = default;
    // BaseMultiFabSet& operator= (const BaseMultiFabSet<MF>& rhs) = default;
    // BaseMultiFabSet& operator= (BaseMultiFabSet<MF>&& rhs) = default;

    // BaseMultiFabSet (BaseMultiFabSet<MF>&& rhs) noexcept;
    // BaseMultiFabSet (const BaseMultiFabSet<MF>& rhs) = delete;
    // BaseMultiFabSet& operator= (BaseMultiFabSet<MF>&& rhs) noexcept;
    // BaseMultiFabSet& operator= (const BaseMultiFabSet<MF>& rhs) = delete;

protected:
    int m_nSet;
    bool m_ncompRefersToSetComps = true;

    // MF** m_mfptr_set = nullptr;
    // IndexType* m_ixtype_set = nullptr;
    // IntVect* m_ngrow_set = nullptr;

    // std::vector<std::shared_ptr<MF>> m_mfptr_set;
    
    std::vector<MF*> m_mfptr_set;
    std::vector<IndexType> m_ixtype_set;
    std::vector<IntVect> m_ngrow_set;

    
    

};



    // template <std::size_t N> BaseMultiFabSet (const Array<const MF*, N>& mfarray);
    // template <std::size_t N> BaseMultiFabSet (const Array<      MF , N>& mfarray);
    // template <std::size_t N> BaseMultiFabSet (const Array<      MF*, N>& mfarray);
    // template <std::size_t N> BaseMultiFabSet (      Array<      MF*, N>  mfarray);
    // template <std::size_t N> BaseMultiFabSet (      Array<      MF , N>& mfarray);
    // template <std::size_t N> BaseMultiFabSet (      Array<const MF*, N>  mfarray);
    // template <std::size_t N> BaseMultiFabSet (      Array<const MF*, N>& mfarray);


};

template <std::size_t N>
struct MultiFabSetInfo {

    MultiFabSetInfo () {
        for (int i = 0; i < N; i++) {
            ix_type_set[i] = IndexType::TheNodeType();
            ngrow_set[i] = IntVect::TheZeroVector();
        }
    };

    MultiFabSetInfo (const Array<IndexType,N>& a_ix_type_set) {
        for (int i = 0; i < N; i++) {
            ix_type_set[i] = IndexType::TheNodeType();
            ngrow_set[i] = IntVect::TheZeroVector();
        }
    };

    Array<IndexType,N> ix_type_set;
    Array<IntVect,N> ngrow_set;
    IndexType[] ix_types_test;

};

