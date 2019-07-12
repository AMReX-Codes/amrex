#include <array>

#ifndef WARPX_PML_H_
#define WARPX_PML_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#if (AMREX_SPACEDIM == 3)

#define WRPX_PML_TO_FORTRAN(x)                              \
    (x).sigma_fac[0].data(), (x).sigma_fac[0].m_lo, (x).sigma_fac[0].m_hi, \
    (x).sigma_fac[1].data(), (x).sigma_fac[1].m_lo, (x).sigma_fac[1].m_hi, \
    (x).sigma_fac[2].data(), (x).sigma_fac[2].m_lo, (x).sigma_fac[2].m_hi, \
    (x).sigma_star_fac[0].data(), (x).sigma_star_fac[0].m_lo, (x).sigma_star_fac[0].m_hi, \
    (x).sigma_star_fac[1].data(), (x).sigma_star_fac[1].m_lo, (x).sigma_star_fac[1].m_hi, \
    (x).sigma_star_fac[2].data(), (x).sigma_star_fac[2].m_lo, (x).sigma_star_fac[2].m_hi

#else

#define WRPX_PML_TO_FORTRAN(x)                              \
    (x).sigma_fac[0].data(), (x).sigma_fac[0].m_lo, (x).sigma_fac[0].m_hi, \
    (x).sigma_fac[1].data(), (x).sigma_fac[1].m_lo, (x).sigma_fac[1].m_hi, \
    (x).sigma_star_fac[0].data(), (x).sigma_star_fac[0].m_lo, (x).sigma_star_fac[0].m_hi, \
    (x).sigma_star_fac[1].data(), (x).sigma_star_fac[1].m_lo, (x).sigma_star_fac[1].m_hi

#endif

struct Sigma : amrex::Vector<amrex::Real>
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

    SigmaVect sigma;      // sigma/epsilon
    SigmaVect sigma_star; // sigma_star/mu
    SigmaVect sigma_fac;
    SigmaVect sigma_star_fac;
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
         int ncell, int delta, int ref_ratio, int do_dive_cleaning, int do_moving_window);

    void ComputePMLFactors (amrex::Real dt);

    std::array<amrex::MultiFab*,3> GetE_fp ();
    std::array<amrex::MultiFab*,3> GetB_fp ();
    std::array<amrex::MultiFab*,3> GetE_cp ();
    std::array<amrex::MultiFab*,3> GetB_cp ();

    amrex::MultiFab* GetF_fp ();
    amrex::MultiFab* GetF_cp ();

    const MultiSigmaBox& GetMultiSigmaBox_fp () const
        { return *sigba_fp; }

    const MultiSigmaBox& GetMultiSigmaBox_cp () const
        { return *sigba_cp; }

    void ExchangeB (const std::array<amrex::MultiFab*,3>& B_fp,
                    const std::array<amrex::MultiFab*,3>& B_cp);
    void ExchangeE (const std::array<amrex::MultiFab*,3>& E_fp,
                    const std::array<amrex::MultiFab*,3>& E_cp);
    void ExchangeB (PatchType patch_type,
		    const std::array<amrex::MultiFab*,3>& Bp);
    void ExchangeE (PatchType patch_type,
		    const std::array<amrex::MultiFab*,3>& Ep);

    void ExchangeF (amrex::MultiFab* F_fp, amrex::MultiFab* F_cp);
    void ExchangeF (PatchType patch_type, amrex::MultiFab* Fp);

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

private:
    bool m_ok;

    const amrex::Geometry* m_geom;
    const amrex::Geometry* m_cgeom;

    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E_fp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B_fp;

    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E_cp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B_cp;

    std::unique_ptr<amrex::MultiFab> pml_F_fp;
    std::unique_ptr<amrex::MultiFab> pml_F_cp;

    std::unique_ptr<MultiSigmaBox> sigba_fp;
    std::unique_ptr<MultiSigmaBox> sigba_cp;

    static amrex::BoxArray MakeBoxArray (const amrex::Geometry& geom,
                                         const amrex::BoxArray& grid_ba, int ncell);

    static void Exchange (amrex::MultiFab& pml, amrex::MultiFab& reg, const amrex::Geometry& geom);
};

#endif
