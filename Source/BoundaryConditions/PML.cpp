
#include <PML.H>
#include <WarpX.H>
#include <WarpXConst.H>

#include <AMReX_Print.H>
#include <AMReX_VisMF.H>

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

namespace
{
    static void FillLo (int idim, Sigma& sigma, Sigma& sigma_star,
                        const Box& overlap, const Box& grid, Real fac)
    {
        int glo = grid.smallEnd(idim);
        int olo = overlap.smallEnd(idim);
        int ohi = overlap.bigEnd(idim);
        int slo = sigma.m_lo;
        int sslo = sigma_star.m_lo;
        for (int i = olo; i <= ohi+1; ++i)
        {
            Real offset = static_cast<Real>(glo-i);
            sigma[i-slo] = fac*(offset*offset);
        }
        for (int i = olo; i <= ohi; ++i)
        {
            Real offset = static_cast<Real>(glo-i) - 0.5;
            sigma_star[i-sslo] = fac*(offset*offset);
        }
    }

    static void FillHi (int idim, Sigma& sigma, Sigma& sigma_star,
                        const Box& overlap, const Box& grid, Real fac)
    {
        int ghi = grid.bigEnd(idim);
        int olo = overlap.smallEnd(idim);
        int ohi = overlap.bigEnd(idim);
        int slo = sigma.m_lo;
        int sslo = sigma_star.m_lo;
        for (int i = olo; i <= ohi+1; ++i)
        {
            Real offset = static_cast<Real>(i-ghi-1);
            sigma[i-slo] = fac*(offset*offset);
        }
        for (int i = olo; i <= ohi; ++i)
        {
            Real offset = static_cast<Real>(i-ghi) - 0.5;
            sigma_star[i-sslo] = fac*(offset*offset);
        }
    }

    static void FillZero (int idim, Sigma& sigma, Sigma& sigma_star, const Box& overlap)
    {
        int olo = overlap.smallEnd(idim);
        int ohi = overlap.bigEnd(idim);
        int slo = sigma.m_lo;
        int sslo = sigma_star.m_lo;
        std::fill(sigma.begin()+(olo-slo), sigma.begin()+(ohi+2-slo), 0.0);
        std::fill(sigma_star.begin()+(olo-sslo), sigma_star.begin()+(ohi+1-sslo), 0.0);
    }
}

SigmaBox::SigmaBox (const Box& box, const BoxArray& grids, const Real* dx, int ncell, int delta)
{
    BL_ASSERT(box.cellCentered());

    const IntVect& sz = box.size();
    const int*     lo = box.loVect();
    const int*     hi = box.hiVect();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        sigma         [idim].resize(sz[idim]+1);
        sigma_star    [idim].resize(sz[idim]  );
        sigma_fac     [idim].resize(sz[idim]+1);
        sigma_star_fac[idim].resize(sz[idim]  );

        sigma         [idim].m_lo = lo[idim];
        sigma         [idim].m_hi = hi[idim]+1;
        sigma_star    [idim].m_lo = lo[idim];
        sigma_star    [idim].m_hi = hi[idim];
        sigma_fac     [idim].m_lo = lo[idim];
        sigma_fac     [idim].m_hi = hi[idim]+1;
        sigma_star_fac[idim].m_lo = lo[idim];
        sigma_star_fac[idim].m_hi = hi[idim];
    }

    Array<Real,AMREX_SPACEDIM> fac;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fac[idim] = 4.0*PhysConst::c/(dx[idim]*static_cast<Real>(delta*delta));
    }

    const std::vector<std::pair<int,Box> >& isects = grids.intersections(box, false, ncell);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        int jdim = (idim+1) % AMREX_SPACEDIM;
#if (AMREX_SPACEDIM == 3)
        int kdim = (idim+2) % AMREX_SPACEDIM;
#endif

        Vector<int> direct_faces, side_faces, direct_side_edges, side_side_edges, corners;
        for (const auto& kv : isects)
        {
            const Box& grid_box = grids[kv.first];

            if (amrex::grow(grid_box, idim, ncell).intersects(box))
            {
                direct_faces.push_back(kv.first);
            }
            else if (amrex::grow(grid_box, jdim, ncell).intersects(box))
            {
                side_faces.push_back(kv.first);
            }
#if (AMREX_SPACEDIM == 3)
            else if (amrex::grow(grid_box, kdim, ncell).intersects(box))
            {
                side_faces.push_back(kv.first);
            }
            else if (amrex::grow(amrex::grow(grid_box,idim,ncell),
                                 jdim,ncell).intersects(box))
            {
                direct_side_edges.push_back(kv.first);
            }
            else if (amrex::grow(amrex::grow(grid_box,idim,ncell),
                                 kdim,ncell).intersects(box))
            {
                direct_side_edges.push_back(kv.first);
            }
            else if (amrex::grow(amrex::grow(grid_box,jdim,ncell),
                                 kdim,ncell).intersects(box))
            {
                side_side_edges.push_back(kv.first);
            }
#endif
            else
            {
                corners.push_back(kv.first);
            }
        }

        for (auto gid : corners)
        {
            const Box& grid_box = grids[gid];

            Box lobox = amrex::adjCellLo(grid_box, idim, ncell);
            lobox.grow(jdim,ncell);
#if (AMREX_SPACEDIM == 3)
            lobox.grow(kdim,ncell);
#endif
            Box looverlap = lobox & box;
            if (looverlap.ok()) {
                FillLo(idim, sigma[idim], sigma_star[idim], looverlap, grid_box, fac[idim]);
            }

            Box hibox = amrex::adjCellHi(grid_box, idim, ncell);
            hibox.grow(jdim,ncell);
#if (AMREX_SPACEDIM == 3)
            hibox.grow(kdim,ncell);
#endif
            Box hioverlap = hibox & box;
            if (hioverlap.ok()) {
                FillHi(idim, sigma[idim], sigma_star[idim], hioverlap, grid_box, fac[idim]);
            }

            if (!looverlap.ok() && !hioverlap.ok()) {
                amrex::Abort("SigmaBox::SigmaBox(): corners, how did this happen?\n");
            }
        }

#if (AMREX_SPACEDIM == 3)
        for (auto gid : side_side_edges)
        {
            const Box& grid_box = grids[gid];
            const Box& overlap = amrex::grow(amrex::grow(grid_box,jdim,ncell),kdim,ncell) & box;
            if (overlap.ok()) {
                FillZero(idim, sigma[idim], sigma_star[idim], overlap);
            } else {
                amrex::Abort("SigmaBox::SigmaBox(): side_side_edges, how did this happen?\n");
            }
        }

        for (auto gid : direct_side_edges)
        {
            const Box& grid_box = grids[gid];

            Box lobox = amrex::adjCellLo(grid_box, idim, ncell);
            Box looverlap = lobox.grow(jdim,ncell).grow(kdim,ncell) & box;
            if (looverlap.ok()) {
                FillLo(idim, sigma[idim], sigma_star[idim], looverlap, grid_box, fac[idim]);
            }

            Box hibox = amrex::adjCellHi(grid_box, idim, ncell);
            Box hioverlap = hibox.grow(jdim,ncell).grow(kdim,ncell) & box;
            if (hioverlap.ok()) {
                FillHi(idim, sigma[idim], sigma_star[idim], hioverlap, grid_box, fac[idim]);
            }

            if (!looverlap.ok() && !hioverlap.ok()) {
                amrex::Abort("SigmaBox::SigmaBox(): direct_side_edges, how did this happen?\n");
            }
        }
#endif

        for (auto gid : side_faces)
        {
            const Box& grid_box = grids[gid];
#if (AMREX_SPACEDIM == 2)
            const Box& overlap = amrex::grow(grid_box,jdim,ncell) & box;
#else
            const Box& overlap = amrex::grow(amrex::grow(grid_box,jdim,ncell),kdim,ncell) & box;
#endif
            if (overlap.ok()) {
                FillZero(idim, sigma[idim], sigma_star[idim], overlap);
            } else {
                amrex::Abort("SigmaBox::SigmaBox(): side_faces, how did this happen?\n");
            }
        }

        for (auto gid : direct_faces)
        {
            const Box& grid_box = grids[gid];

            const Box& lobox = amrex::adjCellLo(grid_box, idim, ncell);
            Box looverlap = lobox & box;
            if (looverlap.ok()) {
                FillLo(idim, sigma[idim], sigma_star[idim], looverlap, grid_box, fac[idim]);
            }

            const Box& hibox = amrex::adjCellHi(grid_box, idim, ncell);
            Box hioverlap = hibox & box;
            if (hioverlap.ok()) {
                FillHi(idim, sigma[idim], sigma_star[idim], hioverlap, grid_box, fac[idim]);
            }

            if (!looverlap.ok() && !hioverlap.ok()) {
                amrex::Abort("SigmaBox::SigmaBox(): direct faces, how did this happen?\n");
            }
        }

        if (direct_faces.size() > 1) {
            amrex::Abort("SigmaBox::SigmaBox(): direct_faces.size() > 1, Box gaps not wide enough?\n");
        }
    }
}

void
SigmaBox::ComputePMLFactorsB (const Real* dx, Real dt)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        for (int i = 0, N = sigma_star[idim].size(); i < N; ++i)
        {
            if (sigma_star[idim][i] == 0.0)
            {
                sigma_star_fac[idim][i] = 1.0;
            }
            else
            {
                sigma_star_fac[idim][i] = std::exp(-sigma_star[idim][i]*dt);
            }
        }
    }
}

void
SigmaBox::ComputePMLFactorsE (const Real* dx, Real dt)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        for (int i = 0, N = sigma[idim].size(); i < N; ++i)
        {
            if (sigma[idim][i] == 0.0)
            {
                sigma_fac[idim][i] = 1.0;
            }
            else
            {
                sigma_fac[idim][i] = std::exp(-sigma[idim][i]*dt);
            }
        }
    }
}

MultiSigmaBox::MultiSigmaBox (const BoxArray& ba, const DistributionMapping& dm,
                              const BoxArray& grid_ba, const Real* dx, int ncell, int delta)
    : FabArray<SigmaBox>(ba,dm,1,0,MFInfo(),
                         FabFactory<SigmaBox>(grid_ba,dx,ncell,delta))
{}

void
MultiSigmaBox::ComputePMLFactorsB (const Real* dx, Real dt)
{
    if (dt == dt_B) return;

    dt_B = dt;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        (*this)[mfi].ComputePMLFactorsB(dx, dt);
    }
}

void
MultiSigmaBox::ComputePMLFactorsE (const Real* dx, Real dt)
{
    if (dt == dt_E) return;

    dt_E = dt;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        (*this)[mfi].ComputePMLFactorsE(dx, dt);
    }
}

PML::PML (const BoxArray& grid_ba, const DistributionMapping& grid_dm,
          const Geometry* geom, const Geometry* cgeom,
          int ncell, int delta, int ref_ratio, int do_dive_cleaning, int do_moving_window)
    : m_geom(geom),
      m_cgeom(cgeom)
{
    const BoxArray& ba = MakeBoxArray(*geom, grid_ba, ncell);
    if (ba.size() == 0) {
        m_ok = false;
        return;
    } else {
        m_ok = true;
    }

    DistributionMapping dm{ba};

    int nge = 2;
    int ngb = 2;
    int ngf = (do_moving_window) ? 2 : 0;
    if (WarpX::maxwell_fdtd_solver_id == 1) ngf = std::max( ngf, 1 );

    pml_E_fp[0].reset(new MultiFab(amrex::convert(ba,WarpX::Ex_nodal_flag), dm, 3, nge));
    pml_E_fp[1].reset(new MultiFab(amrex::convert(ba,WarpX::Ey_nodal_flag), dm, 3, nge));
    pml_E_fp[2].reset(new MultiFab(amrex::convert(ba,WarpX::Ez_nodal_flag), dm, 3, nge));
    pml_B_fp[0].reset(new MultiFab(amrex::convert(ba,WarpX::Bx_nodal_flag), dm, 2, ngb));
    pml_B_fp[1].reset(new MultiFab(amrex::convert(ba,WarpX::By_nodal_flag), dm, 2, ngb));
    pml_B_fp[2].reset(new MultiFab(amrex::convert(ba,WarpX::Bz_nodal_flag), dm, 2, ngb));

    pml_E_fp[0]->setVal(0.0);
    pml_E_fp[1]->setVal(0.0);
    pml_E_fp[2]->setVal(0.0);
    pml_B_fp[0]->setVal(0.0);
    pml_B_fp[1]->setVal(0.0);
    pml_B_fp[2]->setVal(0.0);

    if (do_dive_cleaning)
    {
        pml_F_fp.reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()), dm, 3, ngf));
        pml_F_fp->setVal(0.0);
    }

    sigba_fp.reset(new MultiSigmaBox(ba, dm, grid_ba, geom->CellSize(), ncell, delta));

    if (cgeom)
    {

        nge = 1;
        ngb = 1;

        BoxArray grid_cba = grid_ba;
        grid_cba.coarsen(ref_ratio);
        const BoxArray& cba = MakeBoxArray(*cgeom, grid_cba, ncell);

        DistributionMapping cdm{cba};

        pml_E_cp[0].reset(new MultiFab(amrex::convert(cba,WarpX::Ex_nodal_flag), cdm, 3, nge));
        pml_E_cp[1].reset(new MultiFab(amrex::convert(cba,WarpX::Ey_nodal_flag), cdm, 3, nge));
        pml_E_cp[2].reset(new MultiFab(amrex::convert(cba,WarpX::Ez_nodal_flag), cdm, 3, nge));
        pml_B_cp[0].reset(new MultiFab(amrex::convert(cba,WarpX::Bx_nodal_flag), cdm, 2, ngb));
        pml_B_cp[1].reset(new MultiFab(amrex::convert(cba,WarpX::By_nodal_flag), cdm, 2, ngb));
        pml_B_cp[2].reset(new MultiFab(amrex::convert(cba,WarpX::Bz_nodal_flag), cdm, 2, ngb));

        pml_E_cp[0]->setVal(0.0);
        pml_E_cp[1]->setVal(0.0);
        pml_E_cp[2]->setVal(0.0);
        pml_B_cp[0]->setVal(0.0);
        pml_B_cp[1]->setVal(0.0);
        pml_B_cp[2]->setVal(0.0);

        if (do_dive_cleaning)
        {
            pml_F_cp.reset(new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()), cdm, 3, ngf));
            pml_F_cp->setVal(0.0);
        }

        sigba_cp.reset(new MultiSigmaBox(cba, cdm, grid_cba, cgeom->CellSize(), ncell, delta));
    }

}

BoxArray
PML::MakeBoxArray (const amrex::Geometry& geom, const amrex::BoxArray& grid_ba, int ncell)
{
    Box domain = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if ( ! Geometry::isPeriodic(idim) ) {
            domain.grow(idim, ncell);
        }
    }

    BoxList bl;
    for (int i = 0, N = grid_ba.size(); i < N; ++i)
    {
        const Box& grid_bx = grid_ba[i];
        const IntVect& grid_bx_sz = grid_bx.size();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(grid_bx.shortside() > ncell,
                                         "Consider using larger amr.blocking_factor");

        Box bx = grid_bx;
        bx.grow(ncell);
        bx &= domain;

        Vector<Box> bndryboxes;
#if (AMREX_SPACEDIM == 3)
        int kbegin = -1, kend = 1;
#else
        int kbegin =  0, kend = 0;
#endif
        for (int kk = kbegin; kk <= kend; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
                for (int ii = -1; ii <= 1; ++ii) {
                    if (ii != 0 || jj != 0 || kk != 0) {
                        Box b = grid_bx;
                        b.shift(grid_bx_sz * IntVect{AMREX_D_DECL(ii,jj,kk)});
                        b &= bx;
                        if (b.ok()) {
                            bndryboxes.push_back(b);
                        }
                    }
                }
            }
        }

        const BoxList& noncovered = grid_ba.complementIn(bx);
        for (const Box& b : noncovered) {
            for (const auto& bb : bndryboxes) {
                Box ib = b & bb;
                if (ib.ok()) {
                    bl.push_back(ib);
                }
            }
        }
    }

    BoxArray ba(bl);
    ba.removeOverlap(false);

    return ba;
}

void
PML::ComputePMLFactors (amrex::Real dt)
{
    if (sigba_fp) {
        sigba_fp->ComputePMLFactorsB(m_geom->CellSize(), dt);
        sigba_fp->ComputePMLFactorsE(m_geom->CellSize(), dt);
    }
    if (sigba_cp) {
        sigba_cp->ComputePMLFactorsB(m_cgeom->CellSize(), dt);
        sigba_cp->ComputePMLFactorsE(m_cgeom->CellSize(), dt);
    }
}

std::array<MultiFab*,3>
PML::GetE_fp ()
{
    return {pml_E_fp[0].get(), pml_E_fp[1].get(), pml_E_fp[2].get()};
}

std::array<MultiFab*,3>
PML::GetB_fp ()
{
    return {pml_B_fp[0].get(), pml_B_fp[1].get(), pml_B_fp[2].get()};
}

std::array<MultiFab*,3>
PML::GetE_cp ()
{
    return {pml_E_cp[0].get(), pml_E_cp[1].get(), pml_E_cp[2].get()};
}

std::array<MultiFab*,3>
PML::GetB_cp ()
{
    return {pml_B_cp[0].get(), pml_B_cp[1].get(), pml_B_cp[2].get()};
}

MultiFab*
PML::GetF_fp ()
{
    return pml_F_fp.get();
}

MultiFab*
PML::GetF_cp ()
{
    return pml_F_cp.get();
}

void
PML::ExchangeB (const std::array<amrex::MultiFab*,3>& B_fp,
                const std::array<amrex::MultiFab*,3>& B_cp)
{
  ExchangeB(PatchType::fine, B_fp);
  ExchangeB(PatchType::coarse, B_cp);
}

void
PML::ExchangeB (PatchType patch_type,
                const std::array<amrex::MultiFab*,3>& Bp)
{
    if (patch_type == PatchType::fine && pml_B_fp[0] && Bp[0])
    {
        Exchange(*pml_B_fp[0], *Bp[0], *m_geom);
        Exchange(*pml_B_fp[1], *Bp[1], *m_geom);
        Exchange(*pml_B_fp[2], *Bp[2], *m_geom);
    }
    else if (patch_type == PatchType::coarse && pml_B_cp[0] && Bp[0])
    {
        Exchange(*pml_B_cp[0], *Bp[0], *m_cgeom);
        Exchange(*pml_B_cp[1], *Bp[1], *m_cgeom);
        Exchange(*pml_B_cp[2], *Bp[2], *m_cgeom);
    }
}

void
PML::ExchangeE (const std::array<amrex::MultiFab*,3>& E_fp,
                const std::array<amrex::MultiFab*,3>& E_cp)
{
    ExchangeE(PatchType::fine, E_fp);
    ExchangeE(PatchType::coarse, E_cp);
}

void
PML::ExchangeE (PatchType patch_type,
                const std::array<amrex::MultiFab*,3>& Ep)
{
    if (patch_type == PatchType::fine && pml_E_fp[0] && Ep[0])
    {
        Exchange(*pml_E_fp[0], *Ep[0], *m_geom);
        Exchange(*pml_E_fp[1], *Ep[1], *m_geom);
        Exchange(*pml_E_fp[2], *Ep[2], *m_geom);
    }
    else if (patch_type == PatchType::coarse && pml_E_cp[0] && Ep[0])
    {
        Exchange(*pml_E_cp[0], *Ep[0], *m_cgeom);
        Exchange(*pml_E_cp[1], *Ep[1], *m_cgeom);
        Exchange(*pml_E_cp[2], *Ep[2], *m_cgeom);
    }
}

void
PML::ExchangeF (MultiFab* F_fp, MultiFab* F_cp)
{
    ExchangeF(PatchType::fine, F_fp);
    ExchangeF(PatchType::coarse, F_cp);
}

void
PML::ExchangeF (PatchType patch_type, MultiFab* Fp)
{
    if (patch_type == PatchType::fine && pml_F_fp && Fp) {
        Exchange(*pml_F_fp, *Fp, *m_geom);
    } else if (patch_type == PatchType::coarse && pml_F_cp && Fp) {
        Exchange(*pml_F_cp, *Fp, *m_cgeom);
    }
}

void
PML::Exchange (MultiFab& pml, MultiFab& reg, const Geometry& geom)
{
    const IntVect& ngr = reg.nGrowVect();
    const IntVect& ngp = pml.nGrowVect();
    const int ncp = pml.nComp();
    const auto& period = geom.periodicity();

    MultiFab tmpregmf(reg.boxArray(), reg.DistributionMap(), ncp, ngr);

    if (ngp.max() > 0)  // Copy from pml to the ghost cells of regular data
    {
        MultiFab totpmlmf(pml.boxArray(), pml.DistributionMap(), 1, 0);
        MultiFab::LinComb(totpmlmf, 1.0, pml, 0, 1.0, pml, 1, 0, 1, 0);
        if (ncp == 3) {
            MultiFab::Add(totpmlmf,pml,2,0,1,0);
        }

        MultiFab::Copy(tmpregmf, reg, 0, 0, 1, ngr);
        tmpregmf.ParallelCopy(totpmlmf, 0, 0, 1, IntVect(0), ngr, period);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(reg); mfi.isValid(); ++mfi)
        {
            const FArrayBox& src = tmpregmf[mfi];
            FArrayBox& dst = reg[mfi];
            const BoxList& bl = amrex::boxDiff(dst.box(), mfi.validbox());
            for (const Box& bx : bl)
            {
                dst.copy(src, bx, 0, bx, 0, 1);
            }
        }
    }

    // Copy from regular data to PML's first component
    // Zero out the second (and third) component
    MultiFab::Copy(tmpregmf,reg,0,0,1,0);
    tmpregmf.setVal(0.0, 1, ncp-1, 0);
    pml.ParallelCopy(tmpregmf, 0, 0, ncp, IntVect(0), ngp, period);
}

void
PML::FillBoundary ()
{
    FillBoundaryE();
    FillBoundaryB();
    FillBoundaryF();
}

void
PML::FillBoundaryE ()
{
    FillBoundaryE(PatchType::fine);
    FillBoundaryE(PatchType::coarse);
}

void
PML::FillBoundaryE (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_E_fp[0] && pml_E_fp[0]->nGrowVect().max() > 0)
    {
        const auto& period = m_geom->periodicity();
        Vector<MultiFab*> mf{pml_E_fp[0].get(),pml_E_fp[1].get(),pml_E_fp[2].get()};
        amrex::FillBoundary(mf, period);
    }
    else if (patch_type == PatchType::coarse && pml_E_cp[0] && pml_E_cp[0]->nGrowVect().max() > 0)
    {
        const auto& period = m_cgeom->periodicity();
        Vector<MultiFab*> mf{pml_E_cp[0].get(),pml_E_cp[1].get(),pml_E_cp[2].get()};
        amrex::FillBoundary(mf, period);
    }
}

void
PML::FillBoundaryB ()
{
    FillBoundaryB(PatchType::fine);
    FillBoundaryB(PatchType::coarse);
}

void
PML::FillBoundaryB (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_B_fp[0])
    {
        const auto& period = m_geom->periodicity();
        Vector<MultiFab*> mf{pml_B_fp[0].get(),pml_B_fp[1].get(),pml_B_fp[2].get()};
        amrex::FillBoundary(mf, period);
    }
    else if (patch_type == PatchType::coarse && pml_B_cp[0])
    {
        const auto& period = m_cgeom->periodicity();
        Vector<MultiFab*> mf{pml_B_cp[0].get(),pml_B_cp[1].get(),pml_B_cp[2].get()};
        amrex::FillBoundary(mf, period);
    }
}

void
PML::FillBoundaryF ()
{
    FillBoundaryF(PatchType::fine);
    FillBoundaryF(PatchType::coarse);
}

void
PML::FillBoundaryF (PatchType patch_type)
{
    if (patch_type == PatchType::fine && pml_F_fp && pml_F_fp->nGrowVect().max() > 0)
    {
        const auto& period = m_geom->periodicity();
        pml_F_fp->FillBoundary(period);
    }
    else if (patch_type == PatchType::coarse && pml_F_cp && pml_F_cp->nGrowVect().max() > 0)
    {
        const auto& period = m_cgeom->periodicity();
        pml_F_cp->FillBoundary(period);
    }
}

void
PML::CheckPoint (const std::string& dir) const
{
    if (pml_E_fp[0])
    {
        VisMF::Write(*pml_E_fp[0], dir+"_Ex_fp");
        VisMF::Write(*pml_E_fp[1], dir+"_Ey_fp");
        VisMF::Write(*pml_E_fp[2], dir+"_Ez_fp");
        VisMF::Write(*pml_B_fp[0], dir+"_Bx_fp");
        VisMF::Write(*pml_B_fp[1], dir+"_By_fp");
        VisMF::Write(*pml_B_fp[2], dir+"_Bz_fp");
    }

    if (pml_E_cp[0])
    {
        VisMF::Write(*pml_E_cp[0], dir+"_Ex_cp");
        VisMF::Write(*pml_E_cp[1], dir+"_Ey_cp");
        VisMF::Write(*pml_E_cp[2], dir+"_Ez_cp");
        VisMF::Write(*pml_B_cp[0], dir+"_Bx_cp");
        VisMF::Write(*pml_B_cp[1], dir+"_By_cp");
        VisMF::Write(*pml_B_cp[2], dir+"_Bz_cp");
    }
}

void
PML::Restart (const std::string& dir)
{
    if (pml_E_fp[0])
    {
        VisMF::Read(*pml_E_fp[0], dir+"_Ex_fp");
        VisMF::Read(*pml_E_fp[1], dir+"_Ey_fp");
        VisMF::Read(*pml_E_fp[2], dir+"_Ez_fp");
        VisMF::Read(*pml_B_fp[0], dir+"_Bx_fp");
        VisMF::Read(*pml_B_fp[1], dir+"_By_fp");
        VisMF::Read(*pml_B_fp[2], dir+"_Bz_fp");
    }

    if (pml_E_cp[0])
    {
        VisMF::Read(*pml_E_cp[0], dir+"_Ex_cp");
        VisMF::Read(*pml_E_cp[1], dir+"_Ey_cp");
        VisMF::Read(*pml_E_cp[2], dir+"_Ez_cp");
        VisMF::Read(*pml_B_cp[0], dir+"_Bx_cp");
        VisMF::Read(*pml_B_cp[1], dir+"_By_cp");
        VisMF::Read(*pml_B_cp[2], dir+"_Bz_cp");
    }
}
