
#include <WarpXPML.H>
#include <WarpX.H>
#include <WarpXConst.H>

#include <AMReX_Print.H>

#include <algorithm>

#include <omp.h>

using namespace amrex;

SigmaBox::SigmaBox (const Box& box, const BoxArray& grids, const Real* dx, int ncell)
{
    BL_ASSERT(box.cellCentered());

    const IntVect& sz = box.size();
    const int*     lo = box.loVect();
    const int*     hi = box.hiVect();

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        sigma     [idim].resize(sz[idim]+1, 0.0);
        sigma_star[idim].resize(sz[idim]  , 0.0);

        sigma_fac1     [idim].resize(sz[idim]+1);
        sigma_fac2     [idim].resize(sz[idim]+1);
        sigma_star_fac1[idim].resize(sz[idim]  );
        sigma_star_fac2[idim].resize(sz[idim]  );

        sigma     [idim].m_lo = lo[idim];
        sigma     [idim].m_hi = hi[idim]+1;
        sigma_star[idim].m_lo = lo[idim];
        sigma_star[idim].m_hi = hi[idim];

        sigma_fac1     [idim].m_lo = lo[idim];
        sigma_fac1     [idim].m_hi = hi[idim]+1;
        sigma_fac2     [idim].m_lo = lo[idim];
        sigma_fac2     [idim].m_hi = hi[idim]+1;
        sigma_star_fac1[idim].m_lo = lo[idim];
        sigma_star_fac1[idim].m_hi = hi[idim];
        sigma_star_fac2[idim].m_lo = lo[idim];
        sigma_star_fac2[idim].m_hi = hi[idim];
    }

    Array<Real> fac(BL_SPACEDIM);
    for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
        fac[idim] = 4.0*PhysConst::c/(dx[idim]*static_cast<Real>(ncell*ncell));
    }

    const std::vector<std::pair<int,Box> >& isects = grids.intersections(box, false, ncell);

    for (const auto& kv : isects)
    {
        const Box& grid_box = grids[kv.first];
        const int* glo = grid_box.loVect();
        const int* ghi = grid_box.hiVect();

        for (int idim = 0; idim < BL_SPACEDIM; ++idim)
        {
            // lo side
            for (int i = lo[idim], end=std::min(hi[idim]+1,glo[idim]); i <= end; ++i)
            {
                Real offset = static_cast<Real>(glo[idim]-i);
                sigma[idim][i-lo[idim]] = fac[idim]*(offset*offset);
            }
            for (int i = lo[idim], end=std::min(hi[idim],glo[idim]-1); i <= end; ++i)
            {
                Real offset = static_cast<Real>(glo[idim]-i) - 0.5;
                sigma_star[idim][i-lo[idim]] = fac[idim]*(offset*offset);
            }

            // hi side
            for (int i = std::max(ghi[idim]+1,lo[idim]); i <= hi[idim]+1; ++i)
            {
                Real offset = static_cast<Real>(i-ghi[idim]-1);
                sigma[idim][i-lo[idim]] = fac[idim]*(offset*offset);
            }
            for (int i = std::max(ghi[idim]+1,lo[idim]); i <= hi[idim]; ++i)
            {                
                Real offset = static_cast<Real>(i-ghi[idim]) - 0.5;
                sigma_star[idim][i-lo[idim]] = fac[idim]*(offset*offset);
            }
        }
    }

    // second pass
    for (const auto& kv : isects)
    {
        const Box& grid_box = grids[kv.first];
        for (int idim = 0; idim < BL_SPACEDIM; ++idim)
        {
            for (int orientation = 0; orientation < 2; ++orientation)
            {
                const Box& adjbx = (orientation == 0)
                    ? amrex::adjCellLo(grid_box, idim, ncell)
                    : amrex::adjCellHi(grid_box, idim, ncell);

                const Box& overlap = adjbx & box;
                if (overlap.ok())
                {
                    const int* olo = overlap.loVect();
                    const int* ohi = overlap.hiVect();
                    for (int jdim = 0; jdim < BL_SPACEDIM; ++jdim)
                    {
                        if (jdim != idim)
                        {
                            for (int i = olo[jdim]; i <= ohi[jdim]; ++i) {
                                sigma[jdim][i-lo[jdim]] = 0.0;
                                sigma_star[jdim][i-lo[jdim]] = 0.0;
                            }
                            sigma[jdim][ohi[jdim]+1-lo[jdim]] = 0.0;
                        }
                    }
                }
            }
        }
    }

    // third pass
    for (const auto& kv : isects)
    {
        const Box& grid_box = grids[kv.first];
        const int* glo = grid_box.loVect();
        const int* ghi = grid_box.hiVect();

        for (int idim = 0; idim < BL_SPACEDIM; ++idim)
        {
            const Box& adjlobx = amrex::adjCellLo(grid_box, idim, ncell);
            const Box& ovllobx = adjlobx & box;
            if (ovllobx.ok())
            {
                const int* olo = ovllobx.loVect();
                const int* ohi = ovllobx.hiVect();
                for (int i = olo[idim]; i <= ohi[idim]+1; ++i)
                {
                    Real offset = static_cast<Real>(glo[idim]-i);
                    sigma[idim][i-lo[idim]] = fac[idim]*(offset*offset);
                }
                for (int i = olo[idim]; i <= ohi[idim]; ++i)
                {
                    Real offset = static_cast<Real>(glo[idim]-i) - 0.5;
                    sigma_star[idim][i-lo[idim]] = fac[idim]*(offset*offset);
                }
            }

            const Box& adjhibx = amrex::adjCellHi(grid_box, idim, ncell);
            const Box& ovlhibx = adjhibx & box;
            if (ovlhibx.ok())
            {
                const int* olo = ovlhibx.loVect();
                const int* ohi = ovlhibx.hiVect();
                for (int i = olo[idim]; i <= ohi[idim]+1; ++i)
                {
                    Real offset = static_cast<Real>(i-ghi[idim]-1);
                    sigma[idim][i-lo[idim]] = fac[idim]*(offset*offset);
                }
                for (int i = olo[idim]; i <= ohi[idim]; ++i)
                {
                    Real offset = static_cast<Real>(i-ghi[idim]) - 0.5;
                    sigma_star[idim][i-lo[idim]] = fac[idim]*(offset*offset);
                }
            }
        }
    }
}

void
SigmaBox::ComputePMLFactorsB (const Real* dx, Real dt)
{
    const std::array<Real,BL_SPACEDIM> dtsdx {D_DECL(dt/dx[0], dt/dx[1], dt/dx[2])};

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        for (int i = 0, N = sigma_star[idim].size(); i < N; ++i)
        {
            if (sigma_star[idim][i] == 0.0)
            {
                sigma_star_fac1[idim][i] = 1.0;
                sigma_star_fac2[idim][i] = dtsdx[idim];
            }
            else
            {
                sigma_star_fac1[idim][i] = std::exp(-sigma_star[idim][i]*dt);
                sigma_star_fac2[idim][i] = (1.0-sigma_star_fac1[idim][i])
                    / (sigma_star[idim][i]*dt) * dtsdx[idim];
            }
        }
    }
}

void
SigmaBox::ComputePMLFactorsE (const Real* dx, Real dt)
{
    const std::array<Real,BL_SPACEDIM> dtsdx {D_DECL(dt/dx[0], dt/dx[1], dt/dx[2])};
    
    const Real c2 = PhysConst::c*PhysConst::c;
    const std::array<Real,BL_SPACEDIM> dtsdx_c2 {D_DECL(dtsdx[0]*c2, dtsdx[1]*c2, dtsdx[2]*c2)};

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        for (int i = 0, N = sigma[idim].size(); i < N; ++i)
        {
            if (sigma[idim][i] == 0.0)
            {
                sigma_fac1[idim][i] = 1.0;
                sigma_fac2[idim][i] = dtsdx_c2[idim];
            }
            else
            {
                sigma_fac1[idim][i] = std::exp(-sigma[idim][i]*dt);
                sigma_fac2[idim][i] = (1.0-sigma_fac1[idim][i])
                    / (sigma[idim][i]*dt) * dtsdx_c2[idim];
            }
        }
    }
}

MultiSigmaBox::MultiSigmaBox (const BoxArray& ba, const DistributionMapping& dm,
                              const BoxArray& grid_ba, const Real* dx, int ncell)
    : FabArray<SigmaBox>(ba,dm,1,0,MFInfo(),
                         FabFactory<SigmaBox>(grid_ba,dx,ncell))
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
          int ncell, int ref_ratio)
    : m_geom(geom),
      m_cgeom(cgeom)
{
    const BoxArray& ba = MakeBoxArray(*geom, grid_ba, ncell);
    DistributionMapping dm{ba};

    pml_E_fp[0].reset(new MultiFab(amrex::convert(ba,WarpX::Ex_nodal_flag), dm, 2, 0));
    pml_E_fp[1].reset(new MultiFab(amrex::convert(ba,WarpX::Ey_nodal_flag), dm, 2, 0));
    pml_E_fp[2].reset(new MultiFab(amrex::convert(ba,WarpX::Ez_nodal_flag), dm, 2, 0));
    pml_B_fp[0].reset(new MultiFab(amrex::convert(ba,WarpX::Bx_nodal_flag), dm, 2, 1));
    pml_B_fp[1].reset(new MultiFab(amrex::convert(ba,WarpX::By_nodal_flag), dm, 2, 1));
    pml_B_fp[2].reset(new MultiFab(amrex::convert(ba,WarpX::Bz_nodal_flag), dm, 2, 1));

    pml_E_fp[0]->setVal(0.0);
    pml_E_fp[1]->setVal(0.0);
    pml_E_fp[2]->setVal(0.0);
    pml_B_fp[0]->setVal(0.0);
    pml_B_fp[1]->setVal(0.0);
    pml_B_fp[2]->setVal(0.0);

    sigba_fp.reset(new MultiSigmaBox(ba, dm, grid_ba, geom->CellSize(), ncell));

    if (cgeom)
    {
        BoxArray grid_cba = grid_ba;
        grid_cba.coarsen(ref_ratio);
        const BoxArray& cba = MakeBoxArray(*cgeom, grid_cba, ncell);

        DistributionMapping cdm{cba};

        pml_E_cp[0].reset(new MultiFab(amrex::convert(cba,WarpX::Ex_nodal_flag), cdm, 2, 0));
        pml_E_cp[1].reset(new MultiFab(amrex::convert(cba,WarpX::Ey_nodal_flag), cdm, 2, 0));
        pml_E_cp[2].reset(new MultiFab(amrex::convert(cba,WarpX::Ez_nodal_flag), cdm, 2, 0));
        pml_B_cp[0].reset(new MultiFab(amrex::convert(cba,WarpX::Bx_nodal_flag), cdm, 2, 1));
        pml_B_cp[1].reset(new MultiFab(amrex::convert(cba,WarpX::By_nodal_flag), cdm, 2, 1));
        pml_B_cp[2].reset(new MultiFab(amrex::convert(cba,WarpX::Bz_nodal_flag), cdm, 2, 1));
        
        pml_E_cp[0]->setVal(0.0);
        pml_E_cp[1]->setVal(0.0);
        pml_E_cp[2]->setVal(0.0);
        pml_B_cp[0]->setVal(0.0);
        pml_B_cp[1]->setVal(0.0);
        pml_B_cp[2]->setVal(0.0);

        sigba_cp.reset(new MultiSigmaBox(cba, cdm, grid_cba, cgeom->CellSize(), ncell));
    }

}

BoxArray
PML::MakeBoxArray (const amrex::Geometry& geom, const amrex::BoxArray& grid_ba, int ncell)
{
    Box domain = geom.Domain();
    for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
        if ( ! Geometry::isPeriodic(idim) ) {
            domain.grow(idim, ncell);
        }
    }
    
    BoxList bl;
    for (int i = 0, N = grid_ba.size(); i < N; ++i)
        {
            Box bx = grid_ba[i];
            bx.grow(ncell);
            bx &= domain;
            
            const BoxList& noncovered = grid_ba.complementIn(bx);
            for (const Box& b : noncovered) {
                bl.push_back(b);
            }
        }
    
    BoxArray ba(bl);
    ba.removeOverlap();
    ba.maxSize(64);

    return ba;
}

void
PML::ComputePMLFactorsB (amrex::Real dt)
{
    if (sigba_fp) sigba_fp->ComputePMLFactorsB(m_geom->CellSize(), dt);
    if (sigba_cp) sigba_cp->ComputePMLFactorsB(m_cgeom->CellSize(), dt);
}

void
PML::ComputePMLFactorsE (amrex::Real dt)
{
    if (sigba_fp) sigba_fp->ComputePMLFactorsE(m_geom->CellSize(), dt);
    if (sigba_cp) sigba_cp->ComputePMLFactorsE(m_cgeom->CellSize(), dt);
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

void
PML::Exchange (const std::array<amrex::MultiFab*,3>& E_fp,
               const std::array<amrex::MultiFab*,3>& B_fp,
               const std::array<amrex::MultiFab*,3>& E_cp,
               const std::array<amrex::MultiFab*,3>& B_cp)
{
    Exchange(*pml_E_fp[0], *E_fp[0], *m_geom);
    Exchange(*pml_E_fp[1], *E_fp[1], *m_geom);
    Exchange(*pml_E_fp[2], *E_fp[2], *m_geom);
    Exchange(*pml_B_fp[0], *B_fp[0], *m_geom);
    Exchange(*pml_B_fp[1], *B_fp[1], *m_geom);
    Exchange(*pml_B_fp[2], *B_fp[2], *m_geom);
    if (E_cp[0])
    {
        Exchange(*pml_E_cp[0], *E_cp[0], *m_cgeom);
        Exchange(*pml_E_cp[1], *E_cp[1], *m_cgeom);
        Exchange(*pml_E_cp[2], *E_cp[2], *m_cgeom);
        Exchange(*pml_B_cp[0], *B_cp[0], *m_cgeom);
        Exchange(*pml_B_cp[1], *B_cp[1], *m_cgeom);
        Exchange(*pml_B_cp[2], *B_cp[2], *m_cgeom);
    }
}

void
PML::Exchange (MultiFab& pml, MultiFab& reg, const Geometry& geom)
{
    const int ngr = reg.nGrow();
    const int ngp = pml.nGrow();
    const auto& period = geom.periodicity();

    MultiFab totpmlmf(pml.boxArray(), pml.DistributionMap(), 1, 0);
    MultiFab::LinComb(totpmlmf, 1.0, pml, 0, 1.0, pml, 1, 0, 1, 0);

    MultiFab tmpregmf(reg.boxArray(), reg.DistributionMap(), 2, ngr);
    MultiFab::Copy(tmpregmf, reg, 0, 0, 1, ngr);
    tmpregmf.copy(totpmlmf, 0, 0, 1, 0, ngr, period);

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

    // Copy from regular data to PML's first component
    // Zero out the second component
    MultiFab::Copy(tmpregmf,reg,0,0,1,0);
    tmpregmf.setVal(0.0, 1, 1, 0);
    pml.copy (tmpregmf, 0, 0, 2, 0, ngp, period);
}

void
PML::FillBoundary ()
{
    const auto& period = m_geom->periodicity();
    pml_B_fp[0]->FillBoundary(period);
    pml_B_fp[1]->FillBoundary(period);
    pml_B_fp[2]->FillBoundary(period);
    
    if (pml_B_cp[0])
    {
        const auto& period = m_cgeom->periodicity();
        pml_B_cp[0]->FillBoundary(period);
        pml_B_cp[1]->FillBoundary(period);
        pml_B_cp[2]->FillBoundary(period);
    }
}
