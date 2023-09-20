#include <AMReX_MultiFabSet.H>
// #include <AMReX_Config.H>
// #include <AMReX_MultiFab.H>
// #include <AMReX_BoxArray.H>
// #include <AMReX_DistributionMapping.H>
// #include <AMReX_REAL.H>
// #include <AMReX_Vector.H>

namespace amrex {

using RT = Real;

// template <int N_test>
// void
// amrex::TestFunction () {
//     MultiFabSet<N_test> test_mfs;
//     auto x = test_mfs[1];
//     auto y = test_mfs.m_mf_array[1];
//     auto z = test_mfs.entry(1);
//     int i;
//     auto w = test_mfs.entry(i);

// }

template <std::size_t N>
MultiFabSet<N>::MultiFabSet (const MultiFabSet<N>& src, MakeType maketype, int scomp, int ncomp) {
    for (int i = 0; i < N; ++i) {
        this->setPtr(*MultiFab(src.getElem(i), maketype, scomp, ncomp))
        // (*this).setElem(i, MultiFab(src.getElem(i), maketype, scomp, ncomp));
    }
}

template <std::size_t N>
MultiFabSet<N>::MultiFabSet (const BoxArray& bxs, const DM& dm, int ncomp, const IntVect& ngrow, 
                             const MFInfo& info, const FABFactory& factory, const Array<IntVect,N>* ix_type_array) {
    // Vector<const IntVect&> test_iv_vec;
    // Vector<IntVect> test_iv_vec;
    // Array<IntVect,3> test_iv_vec;
    // Array<IntVect,N> test_iv_vec;
    // const IntVect& test_iv;
    const Array<IntVect,3>* ix_type_array_2;
    Array<const IntVect&,3>* ix_type_array_3;
    if (ix_type_array) {
        for (int i = 0; i < N; ++i) {
            BoxArray bxs_tmp = bxs;
            bxs_tmp.convert((*ix_type_array_2)[i]);
            this->setPtr(*MultiFab(bxs_tmp, dm, ncomp, ngrow, info, factory));
        }
    } else {
        for (int i = 0; i < N; ++i) {
            this->setPtr(*MultiFab(bxs, dm, ncomp, ngrow, info, factory));
        }
    }
}

template <std::size_t N>
MultiFabSet<N>::MultiFabSet (Array<MF,N> mfarray){

}

template <std::size_t N>
void
MultiFabSet<N>::setVal (RT val, int start_index, int num_of_MF) {
    for (int i = start_index; i < start_index + num_of_MF; ++i) {
        (*this).getElem(i).setVal(val);
    };
};

template <std::size_t N>
RT
MultiFabSet<N>::norminf (int comp, int ncomp, IntVect const& nghost, bool local,
                        [[maybe_unused]] bool ignore_covered) const {
    RT result = RT(0);
    for (int i = comp; i < comp + ncomp; ++i) {
        // const MF& mf = (*this)[i];
        const auto& mf = (*this).getElem(i);
        result = std::max(result, mf.norminf(0, mf.nComp(), nghost, local, nghost));
    }
    return result;
}

template<std::size_t N>
void 
MultiFabSet<N>::LocalCopy (const MultiFabSet<N>& src, int scomp, int dcomp, int ncomp,
                           IntVect const& nghost, bool use_nGrowVect = false) {
    // Define common variables:
    int i, n;
    MF& dstmf;
    const MF& srcmf;
    IntVect iv;

    // Two possible ways to write this:
    // 1) Make scomp, dcomp, and ncomp refer to each component of the MultiFabSet:
    BL_ASSERT(scomp == dcomp);
    n = ncomp < 1 ? N : ncomp;
    for (i = scomp; i < n + scomp; ++i) {
        dstmf = (*this)[i];
        srcmf = src[i];
        iv = use_nGrowVect ? srcmf.nGrowVect() : nghost;
        dstmf.LocalCopy(srcmf, 0, 0, srcmf.nComp(), iv);
    }
    // 2) Make scomp, dcomp, and ncomp refer to each component of the underlying MultiFabs:
    for (i = 0; i < N; ++i) {
        dstmf = (*this)[i];
        srcmf = src[i];
        n = ncomp < 1 ? srcmf.nComp() : ncomp;
        iv = use_nGrowVect ? srcmf.nGrowVect() : nghost;
        dstmf.LocalCopy(srcmf, scomp, dcomp, n, iv);
    }
}




}