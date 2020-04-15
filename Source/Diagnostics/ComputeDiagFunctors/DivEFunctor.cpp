#include "DivEFunctor.H"
#include "Utils/Average.H"

using namespace amrex;

DivEFunctor::DivEFunctor(const std::array<const amrex::MultiFab* const, 3> arr_mf_src, const int lev, const amrex::IntVect crse_ratio, const int ncomp)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_arr_mf_src(arr_mf_src), m_lev(lev)
{}

void
DivEFunctor::operator()(amrex::MultiFab& mf_dst, const int dcomp) const
{
    auto& warpx = WarpX::GetInstance();
    // Guard cell is set to 1 for generality. However, for a cell-centered
    // output Multifab, mf_avg, the guard-cell data is not needed especially considering
    // the operations performend in the CoarsenAndInterpolate function.
    constexpr int ng = 1;
    // For staggered and nodal calculations, divE is computed on the nodes.
    // The temporary divE MultiFab is generated to comply with the location of divE.
    const BoxArray& ba = amrex::convert(warpx.boxArray(m_lev),IntVect::TheUnitVector());
    MultiFab divE(ba, warpx.DistributionMap(m_lev), 1, ng );
    warpx.ComputeDivE(divE, m_lev);
    // Coarsen and interpolate from divE on the entire domain to the cell-centered mf_dst.
    Average::CoarsenAndInterpolate(mf_dst, divE, dcomp, 0, nComp(), 0, m_crse_ratio);

}
