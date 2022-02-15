
#include <AMReX_MLPoisson.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_poisson (MLLinOp*& linop, int nlevels, const Geometry* geom[],
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

        MLPoisson* poisson = new MLPoisson(g,b,d,info);
        linop = static_cast<MLLinOp*>(poisson);
    }

}
