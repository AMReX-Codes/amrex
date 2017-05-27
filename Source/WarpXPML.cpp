
#include <WarpXPML.H>
#include <WarpX.H>
#include <WarpXConst.H>

#include <AMReX_Print.H>

#include <algorithm>

#include <omp.h>

using namespace amrex;

SigmaBox::SigmaBox (const Box& a_box, const Box& a_grid_box, const Real* dx, int ncell)
    : box(a_box),
      grid_box(a_grid_box)
{
    BL_ASSERT(box.cellCentered());
    BL_ASSERT(grid_box.cellCentered());

    const IntVect& sz = box.size();

    const int* lo  =      box.loVect();
    const int* hi  =      box.hiVect();
    const int* glo = grid_box.loVect();
    const int* ghi = grid_box.hiVect();

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        sigma     [idim].resize(sz[idim]+1, 0.0);
        sigma_star[idim].resize(sz[idim]  , 0.0);

        sigma_fac1     [idim].resize(sz[idim]+1, 1.0);
        sigma_fac2     [idim].resize(sz[idim]+1     );
        sigma_star_fac1[idim].resize(sz[idim]  , 1.0);
        sigma_star_fac2[idim].resize(sz[idim]       );

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

        const Real fac = 4.0*PhysConst::c/(dx[idim]*static_cast<Real>(ncell*ncell));

        // lo side
        for (int i = lo[idim], end=std::min(hi[idim]+1,glo[idim]); i < end; ++i)
        {
            Real offset = static_cast<Real>(glo[idim]-i);
            sigma[idim][i-lo[idim]] = fac*(offset*offset);

            offset -= 0.5;
            sigma_star[idim][i-lo[idim]] = fac*(offset*offset);
        }

        // hi side
        for (int i = std::max(ghi[idim]+1,lo[idim]); i <= hi[idim]; ++i)
        {
            Real offset = static_cast<Real>(i-ghi[idim]);
            sigma[idim][i-lo[idim]+1] = fac*(offset*offset);
            
            offset -= 0.5;
            sigma_star[idim][i-lo[idim]] = fac*(offset*offset);
        }
    }
}

void
SigmaBox::ComputePMLFactorsB (const Real* dx, Real dt)
{
    const std::array<Real,BL_SPACEDIM> dtsdx {D_DECL(dt/dx[0], dt/dx[1], dt/dx[2])};

    const int* lo  =      box.loVect();
    const int* hi  =      box.hiVect();
    const int* glo = grid_box.loVect();
    const int* ghi = grid_box.hiVect();

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        std::fill(sigma_star_fac2[idim].begin(), sigma_star_fac2[idim].end(), dtsdx[idim]);

        // lo side
        for (int i = lo[idim], end=std::min(hi[idim]+1,glo[idim]); i < end; ++i)
        {
            int ii = i-lo[idim];
            sigma_star_fac1[idim][ii] = std::exp(-sigma_star[idim][ii]*dt);
            sigma_star_fac2[idim][ii] = (1.0-sigma_star_fac1[idim][ii])
                                      / (sigma_star[idim][ii]*dt) * dtsdx[idim];
        }

        // hi side
        for (int i = std::max(ghi[idim]+1,lo[idim]); i <= hi[idim]; ++i)
        {
            int ii = i-lo[idim];
            sigma_star_fac1[idim][ii] = std::exp(-sigma_star[idim][ii]*dt);
            sigma_star_fac2[idim][ii] = (1.0-sigma_star_fac1[idim][ii])
                                      / (sigma_star[idim][ii]*dt) * dtsdx[idim];
        }
    }
}

void
SigmaBox::ComputePMLFactorsE (const Real* dx, Real dt)
{
    const std::array<Real,BL_SPACEDIM> dtsdx {D_DECL(dt/dx[0], dt/dx[1], dt/dx[2])};
    
    const Real c2 = PhysConst::c*PhysConst::c;
    const std::array<Real,BL_SPACEDIM> dtsdx_c2 {D_DECL(dtsdx[0]*c2, dtsdx[1]*c2, dtsdx[2]*c2)};

    const int* lo  =      box.loVect();
    const int* hi  =      box.hiVect();
    const int* glo = grid_box.loVect();
    const int* ghi = grid_box.hiVect();

    for (int idim = 0; idim < BL_SPACEDIM; ++idim)
    {
        std::fill(sigma_fac2[idim].begin(), sigma_fac2[idim].end(), dtsdx_c2[idim]);

        // lo side
        for (int i = lo[idim], end=std::min(hi[idim]+1,glo[idim]); i < end; ++i)
        {
            int ii = i-lo[idim];
            sigma_fac1[idim][ii] = std::exp(-sigma[idim][ii]*dt);
            sigma_fac2[idim][ii] = (1.0-sigma_fac1[idim][ii])
                / (sigma[idim][ii]*dt) * dtsdx_c2[idim];
        }

        // hi side
        for (int i = std::max(ghi[idim]+1,lo[idim]); i <= hi[idim]; ++i)
        {
            int ii = i-lo[idim]+1;
            sigma_fac1[idim][ii] = std::exp(-sigma[idim][ii]*dt);
            sigma_fac2[idim][ii] = (1.0-sigma_fac1[idim][ii])
                / (sigma[idim][ii]*dt) * dtsdx_c2[idim];            
        }
    }
}

SigmaBoxArray::SigmaBoxArray (const BoxArray& ba, const DistributionMapping& dm,
                              const std::map<int,Box>& grid_bxs, const Real* dx, int ncell)
    : FabArray<SigmaBox>(ba,dm,1,0,MFInfo(),
                         FabFactory<SigmaBox>(grid_bxs,dx,ncell))
{}

void
SigmaBoxArray::ComputePMLFactorsB (const Real* dx, Real dt)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        (*this)[mfi].ComputePMLFactorsB(dx, dt);
    }
}

void
SigmaBoxArray::ComputePMLFactorsE (const Real* dx, Real dt)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        (*this)[mfi].ComputePMLFactorsE(dx, dt);
    }
}

PML::PML (const BoxArray& a_ba, const DistributionMapping& a_dm, 
          const Geometry* geom, const Geometry* cgeom,
          int ncell, int ref_ratio)
    : m_mf(a_ba, a_dm, 1, 0, MFInfo().SetAlloc(false)),
      m_cfinfo(FabArrayBase::TheCFinfo(m_mf,*geom,ncell,false,true)),
      m_geom(geom),
      m_cgeom(cgeom)
{
    const BoxArray& ba = m_cfinfo.ba_cfb;
    const DistributionMapping& dm = m_cfinfo.dm_cfb;

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

    std::map<int,Box> bxs;
    for (MFIter mfi(*pml_E_fp[0]); mfi.isValid(); ++mfi)
    {
        int i = mfi.index();
        int li = mfi.LocalIndex();
        bxs[i] = a_ba[m_cfinfo.fine_grid_idx[li]];
    }

    sigba_fp.reset(new SigmaBoxArray(ba, dm, bxs, geom->CellSize(), ncell));

    if (cgeom)
    {
        BoxArray cba = ba;
        cba.coarsen(ref_ratio);

        pml_E_cp[0].reset(new MultiFab(amrex::convert(cba,WarpX::Ex_nodal_flag), dm, 2, 0));
        pml_E_cp[1].reset(new MultiFab(amrex::convert(cba,WarpX::Ey_nodal_flag), dm, 2, 0));
        pml_E_cp[2].reset(new MultiFab(amrex::convert(cba,WarpX::Ez_nodal_flag), dm, 2, 0));
        pml_B_cp[0].reset(new MultiFab(amrex::convert(cba,WarpX::Bx_nodal_flag), dm, 2, 1));
        pml_B_cp[1].reset(new MultiFab(amrex::convert(cba,WarpX::By_nodal_flag), dm, 2, 1));
        pml_B_cp[2].reset(new MultiFab(amrex::convert(cba,WarpX::Bz_nodal_flag), dm, 2, 1));
        
        pml_E_cp[0]->setVal(0.0);
        pml_E_cp[1]->setVal(0.0);
        pml_E_cp[2]->setVal(0.0);
        pml_B_cp[0]->setVal(0.0);
        pml_B_cp[1]->setVal(0.0);
        pml_B_cp[2]->setVal(0.0);

        std::map<int,Box> cbxs;
        for (const auto& kv : bxs)
        {
            cbxs[kv.first] = amrex::coarsen(kv.second,ref_ratio);
        }

        sigba_cp.reset(new SigmaBoxArray(cba, dm, cbxs, cgeom->CellSize(), ncell));
    }

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
    Exchange(*pml_E_fp[0], *E_fp[0]);
    Exchange(*pml_E_fp[1], *E_fp[1]);
    Exchange(*pml_E_fp[2], *E_fp[2]);
    Exchange(*pml_B_fp[0], *B_fp[0]);
    Exchange(*pml_B_fp[1], *B_fp[1]);
    Exchange(*pml_B_fp[2], *B_fp[2]);
    if (E_cp[0])
    {
        Exchange(*pml_E_fp[0], *E_fp[0]);
        Exchange(*pml_E_fp[1], *E_fp[1]);
        Exchange(*pml_E_fp[2], *E_fp[2]);
        Exchange(*pml_B_fp[0], *B_fp[0]);
        Exchange(*pml_B_fp[1], *B_fp[1]);
        Exchange(*pml_B_fp[2], *B_fp[2]);
    }
}

void
PML::Exchange (MultiFab& pml, MultiFab& reg)
{
    const Array<int>& ridx = m_cfinfo.fine_grid_idx;
    const int ngr = reg.nGrow();
    const int ngp = pml.nGrow();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(pml); mfi.isValid(); ++mfi)
    {
        FArrayBox& pfab = pml[mfi];
        FArrayBox& rfab = reg[ridx[mfi.LocalIndex()]];
        const Box& pbx = pfab.box();
        const Box& rbx = rfab.box();
        const Box& vrbx = amrex::grow(rbx, -ngr);
        const Box& vpbx = amrex::grow(pbx, -ngp);
        const BoxList& bl = amrex::boxDiff(vpbx & rbx, vrbx);
        for (const Box& bx: bl)
        {
            rfab.linComb(pfab, bx, 0, pfab, bx, 1, 1.0, 1.0, bx, 0, 1);
        }

        const Box& bx = pbx & vrbx;
        if (bx.ok()) {
            pfab.copy(rfab, bx, 0, bx, 0, 1);
            pfab.setVal(0.0, bx, 1, 1);
        }
    }
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
