
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MultiFabUtil.H>

#include <AMReX_MLABecLap_F.H>
#include <AMReX_ABec_F.H>

namespace amrex {

MLEBABecLap::MLEBABecLap (const Vector<Geometry>& a_geom,
                          const Vector<BoxArray>& a_grids,
                          const Vector<DistributionMapping>& a_dmap,
                          const LPInfo& a_info,
                          const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

void
MLEBABecLap::define (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info,
                     const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    BL_PROFILE("MLEBABecLap::define()");

    Vector<FabFactory<FArrayBox> const*> _factory;
    for (auto x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }

    MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, _factory);
}

MLEBABecLap::~MLEBABecLap ()
{}

void
MLEBABecLap::prepareForSolve ()
{

}

void
MLEBABecLap::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{

}

void
MLEBABecLap::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const
{
}

void
MLEBABecLap::FFlux (int amrlev, const MFIter& mfi, const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                    const FArrayBox& sol, const int face_only) const
{
}

void
MLEBABecLap::normalize (int amrlev, int mglev, MultiFab& mf) const
{
}

}
