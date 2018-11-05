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

}
