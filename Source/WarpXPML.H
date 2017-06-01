#include <array>

#ifndef WARPX_PML_H_
#define WARPX_PML_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

#if (BL_SPACEDIM == 3)

#define WRPX_PML_SIGMA_STAR_TO_FORTRAN(x)  \
    (x).sigma_star_fac1[0].data(), (x).sigma_star_fac1[0].m_lo, (x).sigma_star_fac1[0].m_hi, \
    (x).sigma_star_fac2[0].data(), (x).sigma_star_fac2[0].m_lo, (x).sigma_star_fac2[0].m_hi, \
    (x).sigma_star_fac1[1].data(), (x).sigma_star_fac1[1].m_lo, (x).sigma_star_fac1[1].m_hi, \
    (x).sigma_star_fac2[1].data(), (x).sigma_star_fac2[1].m_lo, (x).sigma_star_fac2[1].m_hi, \
    (x).sigma_star_fac1[2].data(), (x).sigma_star_fac1[2].m_lo, (x).sigma_star_fac1[2].m_hi, \
    (x).sigma_star_fac2[2].data(), (x).sigma_star_fac2[2].m_lo, (x).sigma_star_fac2[2].m_hi

#define WRPX_PML_SIGMA_TO_FORTRAN(x)  \
    (x).sigma_fac1[0].data(), (x).sigma_fac1[0].m_lo, (x).sigma_fac1[0].m_hi, \
    (x).sigma_fac2[0].data(), (x).sigma_fac2[0].m_lo, (x).sigma_fac2[0].m_hi, \
    (x).sigma_fac1[1].data(), (x).sigma_fac1[1].m_lo, (x).sigma_fac1[1].m_hi, \
    (x).sigma_fac2[1].data(), (x).sigma_fac2[1].m_lo, (x).sigma_fac2[1].m_hi, \
    (x).sigma_fac1[2].data(), (x).sigma_fac1[2].m_lo, (x).sigma_fac1[2].m_hi, \
    (x).sigma_fac2[2].data(), (x).sigma_fac2[2].m_lo, (x).sigma_fac2[2].m_hi

#else

#define WRPX_PML_SIGMA_STAR_TO_FORTRAN(x)  \
    (x).sigma_star_fac1[0].data(), (x).sigma_star_fac1[0].m_lo, (x).sigma_star_fac1[0].m_hi, \
    (x).sigma_star_fac2[0].data(), (x).sigma_star_fac2[0].m_lo, (x).sigma_star_fac2[0].m_hi, \
    (x).sigma_star_fac1[1].data(), (x).sigma_star_fac1[1].m_lo, (x).sigma_star_fac1[1].m_hi, \
    (x).sigma_star_fac2[1].data(), (x).sigma_star_fac2[1].m_lo, (x).sigma_star_fac2[1].m_hi

#define WRPX_PML_SIGMA_TO_FORTRAN(x)  \
    (x).sigma_fac1[0].data(), (x).sigma_fac1[0].m_lo, (x).sigma_fac1[0].m_hi, \
    (x).sigma_fac2[0].data(), (x).sigma_fac2[0].m_lo, (x).sigma_fac2[0].m_hi, \
    (x).sigma_fac1[1].data(), (x).sigma_fac1[1].m_lo, (x).sigma_fac1[1].m_hi, \
    (x).sigma_fac2[1].data(), (x).sigma_fac2[1].m_lo, (x).sigma_fac2[1].m_hi

#endif

struct Sigma : amrex::Array<amrex::Real>
{
    int lo() const { return m_lo; }
    int hi() const { return m_hi; }
    int m_lo, m_hi;
};

struct SigmaBox
{
    SigmaBox (const amrex::Box& box, const amrex::BoxArray& grids,
              const amrex::Real* dx, int ncell);

    void ComputePMLFactorsB (const amrex::Real* dx, amrex::Real dt);
    void ComputePMLFactorsE (const amrex::Real* dx, amrex::Real dt);

    std::array<Sigma,BL_SPACEDIM> sigma;      // sigma/epsilon
    std::array<Sigma,BL_SPACEDIM> sigma_star; // sigma_star/mu
    std::array<Sigma,BL_SPACEDIM> sigma_fac1, sigma_fac2;
    std::array<Sigma,BL_SPACEDIM> sigma_star_fac1, sigma_star_fac2;
};

namespace amrex {
    template<>
    class FabFactory<SigmaBox>
    {
    public:
        FabFactory<SigmaBox> (const BoxArray& grid_ba, const Real* dx, int ncell) 
            : m_grids(grid_ba), m_dx(dx), m_ncell(ncell) {}
        virtual SigmaBox* create (const Box& box, int ncomps,
                                  const FabInfo& info, int box_index) const final
            { return new SigmaBox(box, m_grids, m_dx, m_ncell); }
    private:
        const BoxArray& m_grids;
        const Real* m_dx;
        int m_ncell;
    };
}

class MultiSigmaBox
    : public amrex::FabArray<SigmaBox>
{
public:
    MultiSigmaBox(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                  const amrex::BoxArray& grid_ba, const amrex::Real* dx, int ncell);
    void ComputePMLFactorsB (const amrex::Real* dx, amrex::Real dt);
    void ComputePMLFactorsE (const amrex::Real* dx, amrex::Real dt);
private:
    amrex::Real dt_B = -1.e10;
    amrex::Real dt_E = -1.e10;
};

class PML
{
public:
    PML (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, 
         const amrex::Geometry* geom, const amrex::Geometry* cgeom,
         int ncell, int ref_ratio);

    void ComputePMLFactorsB (amrex::Real dt);
    void ComputePMLFactorsE (amrex::Real dt);

    std::array<amrex::MultiFab*,3> GetE_fp ();
    std::array<amrex::MultiFab*,3> GetB_fp ();
    std::array<amrex::MultiFab*,3> GetE_cp ();
    std::array<amrex::MultiFab*,3> GetB_cp ();

    const MultiSigmaBox& GetMultiSigmaBox_fp () const 
        { return *sigba_fp; }

    const MultiSigmaBox& GetMultiSigmaBox_cp () const 
        { return *sigba_cp; }

    void Exchange (const std::array<amrex::MultiFab*,3>& E_fp,
                   const std::array<amrex::MultiFab*,3>& B_fp,
                   const std::array<amrex::MultiFab*,3>& E_cp,
                   const std::array<amrex::MultiFab*,3>& B_cp);

    void FillBoundary ();
    void FillBoundaryE ();
    void FillBoundaryB ();

private:
    const amrex::Geometry* m_geom;
    const amrex::Geometry* m_cgeom;

    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E_fp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B_fp;

    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E_cp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B_cp;

    std::unique_ptr<MultiSigmaBox> sigba_fp;
    std::unique_ptr<MultiSigmaBox> sigba_cp;

    static amrex::BoxArray MakeBoxArray (const amrex::Geometry& geom,
                                         const amrex::BoxArray& grid_ba, int ncell);

    static void Exchange (amrex::MultiFab& pml, amrex::MultiFab& reg, const amrex::Geometry& geom);
};

#endif
