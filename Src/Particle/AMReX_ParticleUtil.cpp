#include <AMReX_ParticleUtil.H>

namespace amrex
{

IntVect computeRefFac (const ParGDBBase* a_gdb, int src_lev, int lev)
{
    IntVect ref_fac = IntVect(AMREX_D_DECL(1,1,1));
    if (src_lev < lev) {
        for (int l = src_lev; l < lev; ++l) {
            ref_fac *= a_gdb->refRatio(l);
        }
    } else if (src_lev > lev) {
        for (int l = src_lev; l > lev; --l) {
            ref_fac *= a_gdb->refRatio(l-1);
        }
        ref_fac *= -1;
    }
    return ref_fac;
}

Vector<int> computeNeighborProcs (const ParGDBBase* a_gdb, int ngrow)
{
    BL_PROFILE("amrex::computeNeighborProcs");

    Vector<int> neighbor_procs;
    for (int src_lev = 0; src_lev < a_gdb->finestLevel()+1; ++src_lev)
    {
        const auto& src_ba = a_gdb->ParticleBoxArray(src_lev);
        const auto& src_dm = a_gdb->ParticleDistributionMap(src_lev);
        for (MFIter mfi(src_ba, src_dm); mfi.isValid(); ++mfi)
        {
            const Box& src_box = mfi.validbox();
            std::vector< std::pair<int, Box> > isects;
            for (int lev = 0; lev < a_gdb->finestLevel()+1; ++lev)
            {
                Box box = src_box;
                const IntVect& ref_fac = computeRefFac(a_gdb, src_lev, lev);
                if (ref_fac < IntVect::TheZeroVector()) box.coarsen(-1*ref_fac);
                else if (ref_fac > IntVect::TheZeroVector()) box.refine(ref_fac);
                box.grow(computeRefFac(a_gdb, 0, src_lev)*ngrow);

                const Periodicity& periodicity = a_gdb->Geom(lev).periodicity();
                const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();
                const BoxArray& ba = a_gdb->ParticleBoxArray(lev);

                for (auto pit=pshifts.cbegin(); pit!=pshifts.cend(); ++pit)
                {
                    const Box& pbox = box + (*pit);
                    bool first_only = false;
                    ba.intersections(pbox, isects, first_only, 0);
                    for (const auto& isec : isects)
                    {
                        const int grid = isec.first;
                        const int global_proc = a_gdb->ParticleDistributionMap(lev)[grid];
                        const int proc = ParallelContext::global_to_local_rank(global_proc);
                        neighbor_procs.push_back(proc);
                    }
                }
            }
        }
    }

    RemoveDuplicates(neighbor_procs);
    return neighbor_procs;
}

}
