
#include <AMReX_MLABecLaplacian.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_abeclaplacian (MLLinOp*& linop, int nlevels, const Geometry* geom[],
                                     const BoxArray* ba[], const DistributionMapping* dm[],
                                     int metric_term, int agglomeration, int consolidation,
                                     int max_coarsening_level)
    {
        LPInfo info;
        if (metric_term >= 0) info.setMetricTerm(metric_term);
        if (agglomeration >= 0) info.setAgglomeration(agglomeration);
        if (consolidation >= 0) info.setConsolidation(consolidation);
        info.setMaxCoarseningLevel(max_coarsening_level);
        Vector<Geometry> g;
        Vector<BoxArray> b;
        Vector<DistributionMapping> d;
        for (int i = 0; i < nlevels; ++i) {
            g.push_back(*geom[i]);
            b.push_back(*ba[i]);
            d.push_back(*dm[i]);
        }

        MLABecLaplacian* abeclap = new MLABecLaplacian(g,b,d,info);
        linop = static_cast<MLLinOp*>(abeclap);
    }

    void amrex_fi_abeclap_set_scalars (MLLinOp* linop, Real a, Real b)
    {
        MLABecLaplacian& abeclap = dynamic_cast<MLABecLaplacian&>(*linop);
        abeclap.setScalars(a,b);
    }

    void amrex_fi_abeclap_set_acoeffs (MLLinOp* linop, int amrlev, const MultiFab* alpha)
    {
        MLABecLaplacian& abeclap = dynamic_cast<MLABecLaplacian&>(*linop);
        abeclap.setACoeffs(amrlev, *alpha);
    }

    void amrex_fi_abeclap_set_bcoeffs (MLLinOp* linop, int amrlev, const MultiFab* beta[])
    {
        MLABecLaplacian& abeclap = dynamic_cast<MLABecLaplacian&>(*linop);
        std::array<MultiFab const*, AMREX_SPACEDIM> b{{AMREX_D_DECL(beta[0],beta[1],beta[2])}};
        abeclap.setBCoeffs(amrlev, b);
    }
}
