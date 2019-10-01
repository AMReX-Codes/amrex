#include <AMReX_ParticleUtil.H>

namespace amrex
{

AMREX_GPU_HOST_DEVICE
int getTileIndex (const IntVect& iv, const Box& box, const bool a_do_tiling, 
		  const IntVect& a_tile_size, Box& tbx)
{
    if (a_do_tiling == false) {
        tbx = box;
        return 0;
    } else {
        //
        // This function must be consistent with FabArrayBase::buildTileArray function!!!
        //
        auto tiling_1d = [](int i, int lo, int hi, int tilesize,
                            int& ntile, int& tileidx, int& tlo, int& thi) {
            int ncells = hi-lo+1;
            ntile = amrex::max(ncells/tilesize, 1);
            int ts_right = ncells/ntile;
            int ts_left  = ts_right+1;
            int nleft = ncells - ntile*ts_right;
	    int ii = i - lo;
            int nbndry = nleft*ts_left;
            if (ii < nbndry) {
                tileidx = ii / ts_left; // tiles on the left of nbndry have size of ts_left
                tlo = lo + tileidx * ts_left;
                thi = tlo + ts_left - 1;
            } else {
                tileidx = nleft + (ii-nbndry) / ts_right;  // tiles on the right: ts_right
                tlo = lo + tileidx * ts_right + nleft;
                thi = tlo + ts_right - 1;
            }
        };
        const IntVect& small = box.smallEnd();
        const IntVect& big   = box.bigEnd();
        IntVect ntiles, ivIndex, tilelo, tilehi;

        AMREX_D_TERM(int iv0 = amrex::min(amrex::max(iv[0], small[0]), big[0]);,
		     int iv1 = amrex::min(amrex::max(iv[1], small[1]), big[1]);,
		     int iv2 = amrex::min(amrex::max(iv[2], small[2]), big[2]););

        AMREX_D_TERM(tiling_1d(iv0, small[0], big[0], a_tile_size[0], ntiles[0], ivIndex[0], tilelo[0], tilehi[0]);,
		     tiling_1d(iv1, small[1], big[1], a_tile_size[1], ntiles[1], ivIndex[1], tilelo[1], tilehi[1]);,
		     tiling_1d(iv2, small[2], big[2], a_tile_size[2], ntiles[2], ivIndex[2], tilelo[2], tilehi[2]););

        tbx = Box(tilelo, tilehi);

        return AMREX_D_TERM(ivIndex[0], + ntiles[0]*ivIndex[1], + ntiles[0]*ntiles[1]*ivIndex[2]);
    }
}

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
                        const int proc = a_gdb->ParticleDistributionMap(lev)[grid];
                        if (proc != ParallelDescriptor::MyProc())
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
