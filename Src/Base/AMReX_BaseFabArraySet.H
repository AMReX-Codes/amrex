#include <AMReX_FabArray.H>

namespace amrex {

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


/**
 * \brief 
 * 
 * \tparam MF FabArray or derived class
 * \tparam N Number of FabArray elements in the set
 */
template <class MF>
class BaseMultiFabSet {

public:

    template <std::size_t N> BaseMultiFabSet (const Array<const MF*, N>& mfarray);
    template <std::size_t N> BaseMultiFabSet (const Array<      MF , N>& mfarray);
    template <std::size_t N> BaseMultiFabSet (const Array<      MF*, N>& mfarray);
    template <std::size_t N> BaseMultiFabSet (      Array<      MF*, N>  mfarray);
    template <std::size_t N> BaseMultiFabSet (      Array<      MF , N>& mfarray);
    template <std::size_t N> BaseMultiFabSet (      Array<const MF*, N>  mfarray);
    template <std::size_t N> BaseMultiFabSet (      Array<const MF*, N>& mfarray);
    
    // BaseMultiFabSet (BaseMultiFabSet<MF>&& rhs) noexcept = default;
    // BaseMultiFabSet (const BaseMultiFabSet<MF>& rhs) = default;
    // BaseMultiFabSet& operator= (const BaseMultiFabSet<MF>& rhs) = default;
    // BaseMultiFabSet& operator= (BaseMultiFabSet<MF>&& rhs) = default;

    // BaseMultiFabSet (BaseMultiFabSet<MF>&& rhs) noexcept;
    // BaseMultiFabSet (const BaseMultiFabSet<MF>& rhs) = delete;
    // BaseMultiFabSet& operator= (BaseMultiFabSet<MF>&& rhs) noexcept;
    // BaseMultiFabSet& operator= (const BaseMultiFabSet<MF>& rhs) = delete;

    unsigned int getSize () { return nSetComp; }

protected:

    unsigned int nSetComp;
    IndexType[] ix_type_array_test;


};




// template<class T, std::size_t N>
// class BaseMultiFabSet {};

// template <class MF, std::size_t N>
// class BaseMultiFabSet<MF, typename std::enable_if<std::is_base_of<FabArrayBase>::value>::type, N> : public std::array<MF, N> 
// {
//     FabArrayBase x();
// };

// // the partial specialization of A is enabled via a template parameter
// template<class T, class Enable = void, std::size_t N>
// class A {}; // primary template
 
// template<class T, std::size_t N>
// class A<T, typename std::enable_if<std::is_floating_point<T>::value>::type, N>
// {
// }; // specialization for floating point types


};