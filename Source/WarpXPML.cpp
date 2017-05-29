
#include <WarpXPML.H>
#include <WarpX.H>
#include <WarpXConst.H>

#include <AMReX_Print.H>

#include <algorithm>

#include <omp.h>

using namespace amrex;

namespace
{
    static void FillLo (int glo, int lo, int hi, Real fac,
                        Array<Real>& sigma, Array<Real>& sigma_star)
    {
        for (int i = lo; i <= hi+1; ++i)
        {
            Real offset = static_cast<Real>(glo-i);
            sigma[i-lo] = fac*(offset*offset);
        }
        for (int i = lo; i <= hi; ++i)
        {
            Real offset = static_cast<Real>(glo-i) - 0.5;
            sigma_star[i-lo] = fac*(offset*offset);
        }
    }

    static void FillHi (int ghi, int lo, int hi, Real fac,
                        Array<Real>& sigma, Array<Real>& sigma_star)
    {
        for (int i = lo; i <= hi+1; ++i)
        {
            Real offset = static_cast<Real>(i-ghi-1);
            sigma[i-lo] = fac*(offset*offset);
        }
        for (int i = lo; i <= hi; ++i)
        {
            Real offset = static_cast<Real>(i-ghi) - 0.5;
            sigma_star[i-lo] = fac*(offset*offset);
        }
    }
}

SigmaBox::SigmaBox (const Box& box, const BoxArray& grids, const Real* dx, int ncell)
{
    BL_ASSERT(box.cellCentered());

    const IntVect& sz = box.size();
    const int*     lo = box.loVect();
    const int*     hi = box.hiVect();

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        sigma     [idim].resize(sz[idim]+1);
        sigma_star[idim].resize(sz[idim]  );

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

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        Array<int> direct_faces, side_faces, direct_side_edges, side_side_edges, corners;
        for (const auto& kv : isects)
        {
            const Box& grid_box = grids[kv.first];
            int jdim = (idim+1) % BL_SPACEDIM;
            int kdim = (idim+2) % BL_SPACEDIM;

            if (amrex::adjCellLo(grid_box, idim, ncell).contains(box))
            {
                direct_faces.push_back(kv.first);
            }
            else if (amrex::adjCellHi(grid_box, idim, ncell).contains(box))
            {
                direct_faces.push_back(kv.first);
            }
            else if (amrex::adjCellLo(grid_box, jdim, ncell).contains(box))
            {
                side_faces.push_back(kv.first);
            }
            else if (amrex::adjCellHi(grid_box, jdim, ncell).contains(box))
            {
                side_faces.push_back(kv.first);
            }
            else if (amrex::adjCellLo(grid_box, kdim, ncell).contains(box))
            {
                side_faces.push_back(kv.first);
            }
            else if (amrex::adjCellHi(grid_box, kdim, ncell).contains(box))
            {
                side_faces.push_back(kv.first);
            }
            else if (Box(grid_box).grow(idim, ncell).grow(jdim, ncell).contains(box))
            {
                direct_side_edges.push_back(kv.first);
            }
            else if (Box(grid_box).grow(idim, ncell).grow(kdim, ncell).contains(box))
            {
                direct_side_edges.push_back(kv.first);
            }
            else if (Box(grid_box).grow(jdim, ncell).grow(kdim, ncell).contains(box))
            {
                side_side_edges.push_back(kv.first);
            }
            else if (Box(grid_box).grow(ncell).contains(box))
            {
                corners.push_back(kv.first);
            }
#if 0
            else
            {
                amrex::Abort("SigmaBox::SigmaBox(): box not contained");
            }
#endif
        }

        if (direct_faces.size() > 1)
        {
            amrex::AllPrint() << "SigmaBox::SigmaBox(), box is " << box << "\n";
            for (auto i : direct_faces) {
                amrex::AllPrint() << "Box " << i << ": " << grids[i] << "\n";
            }
            amrex::Abort("SigmaBox::SigmBox(): box gaps not wide enough");
        }

        if (direct_faces.size() == 1)
        {
            const Box& grid_box = grids[direct_faces[0]];
            const int* glo = grid_box.loVect();
            const int* ghi = grid_box.hiVect();

            if (hi[idim] < glo[idim])
            {
                FillLo(glo[idim], lo[idim], hi[idim], fac[idim],
                       sigma[idim], sigma_star[idim]);
            }
            else if (lo[idim] > ghi[idim])
            {
                FillHi(ghi[idim], lo[idim], hi[idim], fac[idim],
                       sigma[idim], sigma_star[idim]);
            }
            else
            {
                amrex::Print() << "SigmaBox::Sigmabox(): direct_faces, how did this happen?";
            }
        }
        else if (side_faces.size() > 0)
        {
            std::fill(     sigma[idim].begin(),      sigma[idim].end(), 0.0);
            std::fill(sigma_star[idim].begin(), sigma_star[idim].end(), 0.0);
        }
        else if (direct_side_edges.size() > 0)
        {
            const Box& grid_box = grids[direct_side_edges[0]];
            const int* glo = grid_box.loVect();
            const int* ghi = grid_box.hiVect();

            if (hi[idim] < glo[idim])
            {
                FillLo(glo[idim], lo[idim], hi[idim], fac[idim],
                       sigma[idim], sigma_star[idim]);
            }
            else if (lo[idim] > ghi[idim])
            {

                FillHi(ghi[idim], lo[idim], hi[idim], fac[idim],
                       sigma[idim], sigma_star[idim]);
            }
            else
            {
                amrex::Print() << "SigmaBox::Sigmabox(): direc_side_edges, how did this happend?";
            }          
        }
        else if (side_side_edges.size() > 0)
        {
            std::fill(     sigma[idim].begin(),      sigma[idim].end(), 0.0);
            std::fill(sigma_star[idim].begin(), sigma_star[idim].end(), 0.0);
        }
        else
        {
            const Box& grid_box = grids[corners[0]];
            const int* glo = grid_box.loVect();
            const int* ghi = grid_box.hiVect();

            if (hi[idim] < glo[idim])
            {
                FillLo(glo[idim], lo[idim], hi[idim], fac[idim],
                       sigma[idim], sigma_star[idim]);
            }
            else if (lo[idim] > ghi[idim])
            {
                FillHi(ghi[idim], lo[idim], hi[idim], fac[idim],
                       sigma[idim], sigma_star[idim]);
            }
            else
            {
                amrex::Print() << "SigmaBox::Sigmabox(): corners, how did this happen?";
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
        const Box& grid_bx = grid_ba[i];
        const IntVect& grid_bx_sz = grid_bx.size();
        BL_ASSERT(grid_bx.shortside() > ncell);

        Box bx = grid_bx;
        bx.grow(ncell);
        bx &= domain;
        
        Array<Box> bndryboxes;
        for (int kk = -1; kk <= 1; ++kk) {
            for (int jj = -1; jj <= 1; ++jj) {
                for (int ii = -1; ii <= 1; ++ii) {
                    if (ii != 0 || jj != 0 || kk != 0) {
                        Box b = grid_bx;
                        b.shift(grid_bx_sz * IntVect{ii,jj,kk});
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
