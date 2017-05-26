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

#endif

struct Sigma : amrex::Array<amrex::Real>
{
    int lo() const { return m_lo; }
    int hi() const { return m_hi; }
    int m_lo, m_hi;
};

struct SigmaBox
{
    SigmaBox (const amrex::Box& a_box, const amrex::Box& grid_box, const amrex::Real* dx, int ncell);

    void ComputePMLFactorsB (const amrex::Real* dx, amrex::Real dt);
    void ComputePMLFactorsE (const amrex::Real* dx, amrex::Real dt);

    std::array<Sigma,BL_SPACEDIM> sigma;      // sigma/epsilon
    std::array<Sigma,BL_SPACEDIM> sigma_star; // sigma_star/mu
    std::array<Sigma,BL_SPACEDIM> sigma_fac1, sigma_fac2;
    std::array<Sigma,BL_SPACEDIM> sigma_star_fac1, sigma_star_fac2;

    amrex::Box box;
    amrex::Box grid_box;
};

namespace amrex {
    template<>
    class FabFactory<SigmaBox>
    {
    public:
        FabFactory<SigmaBox> (const BoxArray& ba, const Real* dx, int ncell) 
            : m_ba(ba), m_dx(dx), m_ncell(ncell) {}
        virtual SigmaBox* create (const Box& box, int ncomps,
                                  const FabInfo& info, int box_index) const final
            { 
                return new SigmaBox(box, m_ba[box_index], m_dx, m_ncell);
            }
    private:
        const BoxArray& m_ba;
        const Real* m_dx;
        int m_ncell;
    };
}

class SigmaBoxArray
    : public amrex::FabArray<SigmaBox>
{
public:
    SigmaBoxArray(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                  const amrex::Real* dx, int ncell);
    void ComputePMLFactorsB (const amrex::Real* dx, amrex::Real dt);
    void ComputePMLFactorsE (const amrex::Real* dx, amrex::Real dt);
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

    const SigmaBoxArray& GetSigmaBoxArray_fp () const 
        { return *sigba_fp; }

    const SigmaBoxArray& GetSigmaBoxArray_cp () const 
        { return *sigba_cp; }

private:
    amrex::MultiFab m_mf;
    const amrex::FabArrayBase::CFinfo& m_cfinfo;
    const amrex::Geometry* m_geom;
    const amrex::Geometry* m_cgeom;

    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E_fp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B_fp;

    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E_cp;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B_cp;

    std::unique_ptr<SigmaBoxArray> sigba_fp;
    std::unique_ptr<SigmaBoxArray> sigba_cp;
};


#if 0
    //
    std::array<PMLSigma,BL_SPACEDIM> pml_sigma, pml_sigma_star; // sigma/epsilon, sigma_star/mu
    std::array<PMLSigma,BL_SPACEDIM> pml_sigma_fac1, pml_sigma_fac2, 
                                     pml_sigma_star_fac1, pml_sigma_star_fac2;

    int pml_ncell = 10;

    struct PMLSigma : amrex::Array<amrex::Real>
    {
        int lo() const { return m_lo; }
        int hi() const { return m_hi; }
        int m_lo, m_hi;
    };


    const amrex::Array<amrex::Real>& getPMLSigma(int dir) {return pml_sigma[dir];};
    const amrex::Array<amrex::Real>& getPMLSigmaStar(int dir) {return pml_sigma_star[dir];};
    void ComputePMLFactors (int lev, amrex::Real dt);
    void ComputePMLFactorsE (int lev, amrex::Real dt);
    void ComputePMLFactorsB (int lev, amrex::Real dt);
    void ExchangeWithPML (amrex::MultiFab& regmf, amrex::MultiFab& pmlmf,
                          const amrex::Geometry& gm);

    void InitPML ();
#endif

#endif
