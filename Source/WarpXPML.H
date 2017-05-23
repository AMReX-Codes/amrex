#ifndef WARPX_PML_H_
#define WARPX_PML_H_

class WarpXPML
{

};


#if 0
        // PML
        if (Geometry::isAllPeriodic()) {
            do_pml = 0;  // no PML for all periodic boundaries
        } else {
            pp.query("do_pml", do_pml);
            pp.query("pml_ncell", pml_ncell);
        }



    amrex::BoxArray pml_ba;
    amrex::DistributionMapping pml_dm;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_E;
    std::array<std::unique_ptr<amrex::MultiFab>,3> pml_B;
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
