#include <AMReX_FFT_Poisson.H> // Put this at the top for testing

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    static_assert(AMREX_SPACEDIM == 3);

    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

        int n_cell_x = 128;
        int n_cell_y = 128;
        int n_cell_z = 128;

        int max_grid_size_x = 32;
        int max_grid_size_y = 32;
        int max_grid_size_z = 32;

        ParmParse pp;
        pp.query("n_cell_x", n_cell_x);
        pp.query("n_cell_y", n_cell_y);
        pp.query("n_cell_z", n_cell_z);
        pp.query("max_grid_size_x", max_grid_size_x);
        pp.query("max_grid_size_y", max_grid_size_y);
        pp.query("max_grid_size_z", max_grid_size_z);

        Box domain(IntVect(0), IntVect(n_cell_x-1,n_cell_y-1,n_cell_z-1));
        BoxArray ba(domain);
        ba.maxSize(IntVect(max_grid_size_x, max_grid_size_y, max_grid_size_z));
        DistributionMapping dm(ba);

        Geometry geom(domain, RealBox(-1._rt, -1._rt, -1._rt, 1._rt, 1._rt, 1._rt),
                      CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});

        auto const& dx = geom.CellSizeArray();
        auto const& problo = geom.ProbLoArray();

        std::array<IndexType,2> ixtypes{IndexType::TheCellType(),
                                        IndexType::TheNodeType()};
        for (auto const ixtype : ixtypes)
        {
            amrex::Print() << "\nTesting " << ixtype << "\n";

            BoxArray const& iba = amrex::convert(ba, ixtype);
            int ng = ixtype.cellCentered() ? 1 : 0;
            MultiFab rho(iba,dm,1,0);
            MultiFab phi(iba,dm,1,ng);
            phi.setVal(std::numeric_limits<Real>::max());

            auto const& rhoma = rho.arrays();

            constexpr int nsub = 4;
            Real dxsub = dx[0]/nsub;
            Real dysub = dx[1]/nsub;
            Real dzsub = dx[2]/nsub;

            ParallelFor(rho, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real x = (i+0.5_rt/nsub)*dx[0] + problo[0];
                Real y = (j+0.5_rt/nsub)*dx[1] + problo[1];
                Real z = (k+0.5_rt/nsub)*dx[2] + problo[2];
                if (ixtype.nodeCentered()) {
                    x -= 0.5_rt*dx[0];
                    y -= 0.5_rt*dx[1];
                    z -= 0.5_rt*dx[2];
                }
                int n = 0;
                for (int isub = 0; isub < nsub; ++isub) {
                for (int jsub = 0; jsub < nsub; ++jsub) {
                for (int ksub = 0; ksub < nsub; ++ksub) {
                    auto xs = x + isub*dxsub;
                    auto ys = y + jsub*dysub;
                    auto zs = z + ksub*dzsub;
                    if ((xs*xs+ys*ys+zs*zs) < 0.25_rt) { ++n; }
                }}}
                rhoma[b](i,j,k) = Real(n) / Real(nsub*nsub*nsub);
            });

            FFT::PoissonOpenBC solver(geom, ixtype, IntVect(ng));
            solver.solve(phi, rho);

            Real mass = rho.sum_unique(0) * dx[0]*dx[1]*dx[2];
            Real offset = ixtype.cellCentered() ? 0.5_rt : 0.0_rt;
            auto x0 = -1._rt + offset*dx[0];
            auto y0 = -1._rt + offset*dx[1];
            auto z0 = -1._rt + offset*dx[2];
            auto r0 = std::sqrt(x0*x0+y0*y0+z0*z0); // radius of the corner cell
            auto expected = -mass/(4._rt*Math::pi<Real>()*r0);
            amrex::Print() << "  Expected phi: " << expected << "\n";

            int iextra = ixtype.cellCentered() ? 1 : 0;

            for (int k = 0; k < 2; ++k) {
            for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                int ii = (i == 0) ? 0 : n_cell_x-iextra;
                int jj = (j == 0) ? 0 : n_cell_y-iextra;
                int kk = (k == 0) ? 0 : n_cell_z-iextra;
                IntVect corner(ii,jj,kk);
                auto v = amrex::get_cell_data(phi, corner);
                if (!v.empty()) {
                    amrex::AllPrint() << "  phi at " << corner << " is " << v[0] << "\n";
                    auto error = std::abs(expected-v[0])/std::max(std::abs(expected),std::abs(v[0]));
                    amrex::AllPrint() << "  error " << error << "\n";
#ifdef AMREX_USE_FLOAT
                    constexpr Real eps = 1.e-5;
#else
                    constexpr Real eps = 1.e-6;
#endif
                    AMREX_ALWAYS_ASSERT(error < eps);
                }
            }}}
        }
    }
    amrex::Finalize();
}
