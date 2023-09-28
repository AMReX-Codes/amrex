#include <AMReX_MultiFabSet.H>
// #include <AMReX_Config.H>
// #include <AMReX_MultiFab.H>
// #include <AMReX_BoxArray.H>
// #include <AMReX_DistributionMapping.H>
// #include <AMReX_REAL.H>
// #include <AMReX_Vector.H>

namespace amrex {

using RT = Real;

MultiFabSet::MultiFabSet (const MultiFabSet& src, MakeType maketype, int scomp, int ncomp)
{
    if (m_ncompRefersToSetComps) {
        initializeSet(ncomp);
        for (int i = 0; i < ncomp; ++i) {
            const MF& srcmf = src.getElem(scomp+i);
            setElem(i, new MF(srcmf, maketype, 0, srcmf.nComp()));
        }
    } else {
        initializeSet(src.nSet());
        for (int i = 0; i < nSet(); ++i) {
            setElem(i, new MF(src.getElem(i), maketype, scomp, ncomp));
        }
    }
}

MultiFabSet::MultiFabSet (const BoxArray& bxs, const DM& dm, int ncomp, const IntVect& ngrow, 
                             const MFInfo& info, const FABFactory& factory, const Vector<IndexType>* ix_type_array)
{
    int ncompMF;
    if (m_ncompRefersToSetComps) {
        initializeSet(ncomp);
        ncompMF = 1;
    } else {
        // initializeSet(???); ???
        ncompMF = ncomp;
    }
    if (ix_type_array) {
        for (int i = 0; i < nSet(); ++i) {
            BoxArray bxs_tmp = bxs;
            bxs_tmp.convert((*ix_type_array)[i]);
            setElem(i, new MF(bxs_tmp, dm, ncompMF, ngrow, info, factory));
        }
    } else {
        for (int i = 0; i < nSet(); ++i) {
            setElem(i, new MF(bxs, dm, ncompMF, ngrow, info, factory));
        }
    }
}

template <std::size_t N>
MultiFabSet::MultiFabSet (Array<MF*,N>& mfarray) 
{
    initializeSet(mfarray.size());
    for (int i = 0; i < ; ++i) {
        setElem(i, mfarray[i]);
    }
}

template <std::size_t N>
MultiFabSet::MultiFabSet (const Array<MF*,N>& mfarray)
{
    initializeSet(mfarray.size());
    for (int i = 0; i < ; ++i) {
        setElem(i, mfarray[i]);
    }
}

void
MultiFabSet::setVal (RT val, int comp, int ncomp)
{
    if (m_ncompRefersToSetComps) {
        for (int i = comp; i < comp + ncomp; ++i) {
            MF& mf = getElem(i);
            mf.setVal(val, 0, mf.nComp());
        };
    } else {
        for (int i = 0; i < nSet(); ++i) {
            getElem(i).setVal(val, comp, ncomp);
        };
    }
};

RT
MultiFabSet::norminf (int comp, int ncomp, IntVect const& nghost, bool local,
                      [[maybe_unused]] bool ignore_covered) const
{
    RT result = RT(0);
    if (m_ncompRefersToSetComps) {
        for (int i = comp; i < comp + ncomp; ++i) {
            const auto& mf = getElem(i);
            result = std::max(result, mf.norminf(0, mf.nComp(), nghost, local, ignore_covered));
        }
    } else {
        for (int i = 0; i < nSet(); ++i) {
            const auto& mf = getElem(i);
            result = std::max(result, mf.norminf(comp, ncomp, nghost, local, ignore_covered));
        }
    }
    return result;
}

void
MultiFabSet::LocalCopy (const MultiFabSet& src, int scomp, int dcomp, int ncomp,
                        IntVect const& nghost, bool use_nGrowVect)
{
    // Define common variables:
    int i, n;
    IntVect iv;

    // Two possible ways to write this:
    // 1) Make scomp, dcomp, and ncomp refer to each component of the MultiFabSet:
    BL_ASSERT(scomp == dcomp);
    n = ncomp < 1 ? src.nSet() : ncomp;
    for (i = scomp; i < n + scomp; ++i) {
        MF& dstmf = (*this)[i];
        const MF& srcmf = src[i];
        iv = use_nGrowVect ? srcmf.nGrowVect() : nghost;
        dstmf.LocalCopy(srcmf, 0, 0, srcmf.nComp(), iv);
    }
    // 2) Make scomp, dcomp, and ncomp refer to each component of the underlying MultiFabs:
    for (i = 0; i < src.nSet(); ++i) {
        MF& dstmf = (*this)[i];
        const MF& srcmf = src[i];
        n = ncomp < 1 ? srcmf.nComp() : ncomp;
        iv = use_nGrowVect ? srcmf.nGrowVect() : nghost;
        dstmf.LocalCopy(srcmf, scomp, dcomp, n, iv);
    }
}

void
MultiFabSet::clear()
{
    for (int i; i < m_nSet; ++i) {
        getElem(i).clear();
    }
    m_mfptr_set.clear();
}

bool
MultiFabSet::isAllRegular() const noexcept
{
    for (int i = 0; i < m_nSet; ++i) {
        if (!(getElem(i).isAllRegular())) {
            return false;
        }
    }
    return true;
}

}