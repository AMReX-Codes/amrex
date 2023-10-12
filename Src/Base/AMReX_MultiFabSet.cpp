#include <AMReX_MultiFabSet.H>
#include <AMReX_BoxArray.H>

namespace amrex {

MultiFabSet::MultiFabSet () noexcept
{}

MultiFabSet::MultiFabSet () noexcept
    : MultiFab()
{}

MultiFabSet::MultiFabSet (Arena* a) noexcept
    : MultiFab(a)
{}

MultiFabSet::MultiFabSet (const BoxArray&              bxs,
                          const DistributionMapping&   dm,
                          const Vector<IndexType>      ixTypes,
                          const IntVect&               ngrow,
                          const MFInfo&                info,
                          const FabFactory<FArrayBox>& factory)
    : MultiFab(amrex::convert(bxs,IndexType::TheNodeType()),dm,ixTypes.size(),ngrow,info,factory)
    , m_ixTypes(ixTypes)
{}

MultiFabSet::MultiFabSet (const MultiFabSet& rhs, MakeType maketype, int scomp, int ncomp)
    : MultiFab(rhs, maketype, scomp, ncomp)
    , m_ixTypes(rhs.m_ixTypes)
{}

MultiFabSet::~MultiFabSet ()
{}

MultiFabSet::MultiFabSet (MultiFabSet&& rhs) noexcept
    : MultiFab(std::move(rhs))
    , m_ixTypes(std::move(rhs.m_ixTypes))
{}

void
MultiFabSet::define (const BoxArray&              bxs,
                     const DistributionMapping&   dm,
                     const Vector<IndexType>      ixTypes,
                     const IntVect&               ngrow,
                     const MFInfo&                info,
                     const FabFactory<FArrayBox>& factory)
{
    MultiFab::define(amrex::convert(bxs,IndexType::TheNodeType()),dm,ixTypes.size(),ngrow,info,factory);
    m_ixTypes = ixTypes;
}

void
MultiFabSet::LocalCopy (const MultiFabSet& src, int scomp, int dcomp, int ncomp, const IntVect& nghost)
{
    if (m_ixTypes.size() < 1) {
        m_ixTypes = Vector<IndexType>(dcomp+ncomp);
    }
    for (int n = 0; n < ncomp; n++) {
        m_ixTypes[dcomp+n] = src.m_ixTypes[scomp+n];
    }
    MultiFab::LocalCopy(src, scomp, dcomp, ncomp, nghost);
}


}