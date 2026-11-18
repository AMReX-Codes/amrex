// FillBoundary in spherical (theta, phi, r) index space with SphThetaPhiRIndexMapping.
// Every filled ghost cell is checked against its source cell, with Identity and with
// SphThetaPhiRComponentMapping (coordinate-basis signs).
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_NonLocalBC.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParReduce.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace {

// Value of valid cell (i,j,k), component c. Never zero, so unfilled ghost cells (0) are detected.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real cell_value (int i, int j, int k, int c)
{
    return Real(1 + i + 10*j + 100*k + 1000*c);
}

// Number of cell values that differ from the expected ones. vcomp < 0: no vector components.
int count_errors (MultiFab const& mf, IntVect const& n_cell, int vcomp)
{
    int const nx = n_cell[0];
    int const ny = n_cell[1];
    int const nz = n_cell[2];
    auto const& ma = mf.const_arrays();
    int nerr = ParReduce(TypeList<ReduceOpSum>{}, TypeList<int>{}, mf, mf.nGrowVect(), mf.nComp(),
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int c) -> GpuTuple<int>
        {
            bool const ilo = (i < 0);
            bool const ihi = (i >= nx);
            bool const klo = (k < 0);
            bool const khi = (k >= nz);
            // Not filled: the outer r boundary and the edges where a pole meets the origin.
            if (khi || (klo && (ilo || ihi))) { return {0}; }
            // Source cell at the same physical point.
            int si = i;
            int sk = k;
            int sj = ((j % ny) + ny) % ny;      // periodic in phi
            if (klo) {                          // through the origin
                si = nx-1-i;
                sk = -1-k;
                sj = (sj + ny/2) % ny;
            } else if (ilo) {                   // beyond theta = 0
                si = -1-i;
                sj = (sj + ny/2) % ny;
            } else if (ihi) {                   // beyond theta = pi
                si = 2*nx-1-i;
                sj = (sj + ny/2) % ny;
            }
            Real sign = Real(1.0);
            if (vcomp >= 0) {
                if ((c == vcomp) && (klo || ilo || ihi)) { sign = Real(-1.0); }
                if ((c == vcomp+2) && klo) { sign = Real(-1.0); }
            }
            return {(ma[b](i,j,k,c) == sign*cell_value(si,sj,sk,c)) ? 0 : 1};
        });
    ParallelDescriptor::ReduceIntSum(nerr);
    return nerr;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        IntVect n_cell(8);     // the number of cells in phi must be even
        int ng = 2;
        int max_grid_size = 4;
        int vcomp = 1;         // the vector (theta, phi, r) in vcomp..vcomp+2, scalars elsewhere
        {
            ParmParse pp;
            Vector<int> nc_in;
            if (pp.queryarr("n_cell", nc_in)) { n_cell = IntVect(nc_in); }
            pp.query("ng", ng);
            pp.query("max_grid_size", max_grid_size);
            pp.query("vcomp", vcomp);
        }
        int const nc = 4;
        AMREX_ALWAYS_ASSERT((n_cell[1] % 2 == 0) && (2*ng <= n_cell[1]) &&
                            (ng <= n_cell[0]) && (ng <= n_cell[2]));

        Box const domain(IntVect(0), n_cell-1);
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        Geometry geom(domain, RealBox({0.,0.,0.}, {1.,1.,1.}), 0, {0,1,0});  // periodic in phi

        MultiFab mf(ba, dm, nc, ng);
        NonLocalBC::SphThetaPhiRIndexMapping dtos(domain);
        auto const cmd = NonLocalBC::makeFillBoundaryMetaData(mf, mf.nGrowVect(), geom, dtos);

        auto init = [&] () {
            mf.setVal(Real(0.0));
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                auto const& a = mf.array(mfi);
                ParallelFor(mfi.validbox(), nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int c)
                {
                    a(i,j,k,c) = cell_value(i,j,k,c);
                });
            }
        };

        init();
        NonLocalBC::FillBoundary(mf, cmd, 0, nc, dtos, NonLocalBC::Identity{});
        int const nerr_identity = count_errors(mf, n_cell, -1);

        init();
        NonLocalBC::FillBoundary(mf, cmd, 0, nc, dtos,
                                 NonLocalBC::SphThetaPhiRComponentMapping(domain, vcomp));
        int const nerr_vector = count_errors(mf, n_cell, vcomp);

        amrex::Print() << "Identity: " << nerr_identity << " wrong values; "
                       << "SphThetaPhiRComponentMapping: " << nerr_vector << " wrong values\n";
        if ((nerr_identity != 0) || (nerr_vector != 0)) {
            amrex::Abort("SphFillBoundary test failed");
        }
    }
    amrex::Finalize();
}
