/* Copyright 2019 Andrew Myers, Aurore Blelly, Axel Huebl
 * Maxence Thevenet, Remi Lehe, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PML_H_
#define WARPX_PML_H_

#include "Utils/WarpXProfilerWrapper.H"

#ifdef WARPX_USE_PSATD
#   include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#include <array>


struct Sigma : amrex::Gpu::ManagedVector<amrex::Real>
{
    int lo() const { return m_lo; }
    int hi() const { return m_hi; }
    int m_lo, m_hi;
};

struct SigmaBox
{
    SigmaBox (const amrex::Box& box, const amrex::BoxArray& grids,
              const amrex::Real* dx, int ncell, int delta);

    void ComputePMLFactorsB (const amrex::Real* dx, amrex::Real dt);
    void ComputePMLFactorsE (const amrex::Real* dx, amrex::Real dt);

    using SigmaVect = std::array<Sigma,AMREX_SPACEDIM>;

    SigmaVect sigma;
    SigmaVect sigma_cumsum;
    SigmaVect sigma_star;
    SigmaVect sigma_star_cumsum;
    SigmaVect sigma_fac;
    SigmaVect sigma_cumsum_fac;
    SigmaVect sigma_star_fac;
    SigmaVect sigma_star_cumsum_fac;

};

namespace amrex {
    template<>
    class FabFactory<SigmaBox>
    {
    public:
        FabFactory<SigmaBox> (const BoxArray& grid_ba, const Real* dx, int ncell, int delta)
            : m_grids(grid_ba), m_dx(dx), m_ncell(ncell), m_delta(delta) {}
        virtual SigmaBox* create (const Box& box, int ncomps,
                                  const FabInfo& info, int box_index) const final
            { return new SigmaBox(box, m_grids, m_dx, m_ncell, m_delta); }
        virtual void destroy (SigmaBox* fab) const final {
            delete fab;
        }
        virtual FabFactory<SigmaBox>* clone () const {
            return new FabFactory<SigmaBox>(*this);
        }
    private:
        const BoxArray& m_grids;
        const Real* m_dx;
        int m_ncell;
        int m_delta;
    };
}

class MultiSigmaBox
    : public amrex::FabArray<SigmaBox>
{
public:
    MultiSigmaBox(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                  const amrex::BoxArray& grid_ba, const amrex::Real* dx, int ncell, int delta);
    void ComputePMLFactorsB (const amrex::Real* dx, amrex::Real dt);
    void ComputePMLFactorsE (const amrex::Real* dx, amrex::Real dt);
private:
    amrex::Real dt_B = -1.e10;
    amrex::Real dt_E = -1.e10;
};

enum struct PatchType : int;

class PML
{
public:
    PML (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
         const amrex::Geometry* geom, const amrex::Geometry* cgeom,
         int ncell, int delta, int ref_ratio,
#ifdef WARPX_USE_PSATD
         amrex::Real dt, int nox_fft, int noy_fft, int noz_fft, bool do_nodal,
#endif
         int do_dive_cleaning, int do_moving_window,
         int pml_has_particles, int do_pml_in_domain,
         const amrex::IntVect do_pml_Lo = amrex::IntVect::TheUnitVector(),
         const amrex::IntVect do_pml_Hi = amrex::IntVect::TheUnitVector());

    void ComputePMLFactors (amrex::Real dt);

    std::array<amrex::MultiFab*,3> GetE_fp ();
    std::array<amrex::MultiFab*,3> GetB_fp ();
    std::array<amrex::MultiFab*,3> Getj_fp ();
    std::array<amrex::MultiFab*,3> GetE_cp ();
    std::array<amrex::MultiFab*,3> GetB_cp ();
    std::array<amrex::MultiFab*,3> Getj_cp ();

    amrex::MultiFab* GetF_fp ();
    amrex::MultiFab* GetF_cp ();

    const MultiSigmaBox& GetMultiSigmaBox_fp () const
        { return *sigba_fp; }

    const MultiSigmaBox& GetMultiSigmaBox_cp () const
        { return *sigba_cp; }

#ifdef WARPX_USE_PSATD
    void PushPSATD ();
#endif

    void ExchangeB (const std::array<amrex::MultiFab*,3>& B_fp,
                    const std::array<amrex::MultiFab*,3>& B_cp, int do_pml_in_domain);
    void ExchangeE (const std::array<amrex::MultiFab*,3>& E_fp,
                    const std::array<amrex::MultiFab*,3>& E_cp, int do_pml_in_domain);
    void CopyJtoPMLs (const std::array<amrex::MultiFab*,3>& j_fp,
                    const std::array<amrex::MultiFab*,3>& j_cp);

    void ExchangeB (PatchType patch_type,
                    const std::array<amrex::MultiFab*,3>& Bp, int do_pml_in_domain);
    void ExchangeE (PatchType patch_type,
                    const std::array<amrex::MultiFab*,3>& Ep, int do_pml_in_domain);

    void CopyJtoPMLs (PatchType patch_type,
                    const std::array<amrex::MultiFab*,3>& jp);

    void ExchangeF (amrex::MultiFab* F_fp, amrex::MultiFab* F_cp, int do_pml_in_domain);
    void ExchangeF (PatchType patch_type, amrex::MultiFab* Fp, int do_pml_in_domain);

    void FillBoundary ();
    void FillBoundaryE ();
    void FillBoundaryB ();
    void FillBoundaryF ();
    void FillBoundaryE (PatchType patch_type);
    void FillBoundaryB (PatchType patch_type);
    void FillBoundaryF (PatchType patch_type);

    bool ok () const { return m_ok; }

    void CheckPoint (const std::string& dir) const;
    void Restart (const std::string& dir);

    static void Exchange (amrex::MultiFab& pml, amrex::MultiFab& reg, const amrex::Geometry& geom, int do_pml_in_domain);

private:
    bool m_ok;

    const amrex::Geometry* m_geom;
    const amrex::Geometry* m_cgeom;

    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E_fp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B_fp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_j_fp;

    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E_cp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B_cp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_j_cp;

    std::unique_ptr<amrex::MultiFab> pml_F_fp;
    std::unique_ptr<amrex::MultiFab> pml_F_cp;

    std::unique_ptr<MultiSigmaBox> sigba_fp;
    std::unique_ptr<MultiSigmaBox> sigba_cp;

#ifdef WARPX_USE_PSATD
    std::unique_ptr<SpectralSolver> spectral_solver_fp;
    std::unique_ptr<SpectralSolver> spectral_solver_cp;
#endif

    static amrex::BoxArray MakeBoxArray (const amrex::Geometry& geom,
                                         const amrex::BoxArray& grid_ba,
                                         int ncell, int do_pml_in_domain,
                                         const amrex::IntVect do_pml_Lo = amrex::IntVect::TheUnitVector(),
                                         const amrex::IntVect do_pml_Hi = amrex::IntVect::TheUnitVector());

    static void CopyToPML (amrex::MultiFab& pml, amrex::MultiFab& reg, const amrex::Geometry& geom);
};

#ifdef WARPX_USE_PSATD
void PushPMLPSATDSinglePatch( SpectralSolver& solver,
    std::array<std::unique_ptr<amrex::MultiFab>,3>& pml_E,
    std::array<std::unique_ptr<amrex::MultiFab>,3>& pml_B );
#endif

#endif
