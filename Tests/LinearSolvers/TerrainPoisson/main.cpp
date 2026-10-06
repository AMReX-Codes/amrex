// Solve the terrain-following Poisson equation with MLTerrainPoisson, using
// MLMG and GMRESMLMG, and compare with a known discrete solution.

#include <AMReX.H>
#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLTerrainPoisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>

using namespace amrex;

namespace {

LinOpBCType parse_bc (std::string const& s)
{
    if (s == "periodic")  { return LinOpBCType::Periodic; }
    if (s == "neumann")   { return LinOpBCType::Neumann; }
    if (s == "dirichlet") { return LinOpBCType::Dirichlet; }
    amrex::Abort("Unknown bc " + s);
    return LinOpBCType::bogus;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        IntVect n_cell(32,32,32);
        pp.queryarr("n_cell", n_cell);
        int max_grid_size = 16;
        pp.query("max_grid_size", max_grid_size);
        IntVect max_grid_size_v(max_grid_size);
        pp.queryarr("max_grid_size_v", max_grid_size_v);

        std::vector<Real> prob_hi{Real(1000.), Real(1000.), Real(500.)};
        pp.queryarr("prob_hi", prob_hi);

        // Hill height, phase shift (nonzero gives slopes at the domain faces)
        // and vertical stretching (0 for none).
        Real hill_height = Real(100.);
        Real hill_shift = Real(0.2);
        Real stretch = Real(2.0);
        pp.query("hill_height", hill_height);
        pp.query("hill_shift", hill_shift);
        pp.query("stretch", stretch);

        // Perturbation of the face area factors, as from map factors.
        Real area_noise = Real(0.1);
        pp.query("area_noise", area_noise);
        int seed = 42;
        pp.query("seed", seed);

        std::vector<std::string> bc_lo{"periodic","periodic","neumann"};
        std::vector<std::string> bc_hi{"periodic","periodic","dirichlet"};
        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);

        int hidden_direction = -1;
        pp.query("hidden_direction", hidden_direction);
        int max_coarsening_level = 30;
        pp.query("max_coarsening_level", max_coarsening_level);
        bool agglomeration = true;
        bool consolidation = true;
        pp.query("agglomeration", agglomeration);
        pp.query("consolidation", consolidation);

        std::string solver = "both";
        pp.query("solver", solver);
        int verbose = 1;
        pp.query("verbose", verbose);
        int bottom_verbose = 0;
        pp.query("bottom_verbose", bottom_verbose);
        int max_iter = 200;
        pp.query("max_iter", max_iter);
        std::string bottom_solver = "bicgstab";
        pp.query("bottom_solver", bottom_solver);
        std::string zsplit_solver = "spike";
        pp.query("zsplit_solver", zsplit_solver);
        Real reltol = std::is_same_v<Real,double> ? Real(1.e-10) : Real(1.e-4);
        pp.query("reltol", reltol);
        Real check_tol = std::is_same_v<Real,double> ? Real(1.e-7) : Real(5.e-2);
        pp.query("check_tol", check_tol);

        amrex::InitRandom(seed);

        Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
        Array<int,AMREX_SPACEDIM> is_periodic{};
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            lobc[idim] = parse_bc(bc_lo[idim]);
            hibc[idim] = parse_bc(bc_hi[idim]);
            AMREX_ALWAYS_ASSERT((lobc[idim] == LinOpBCType::Periodic) ==
                                (hibc[idim] == LinOpBCType::Periodic));
            is_periodic[idim] = (lobc[idim] == LinOpBCType::Periodic);
        }

        Box domain(IntVect(0), n_cell-1);
        RealBox rb({Real(0.),Real(0.),Real(0.)}, {prob_hi[0],prob_hi[1],prob_hi[2]});
        Geometry geom(domain, rb, CoordSys::cartesian, is_periodic);
        BoxArray ba(domain);
        ba.maxSize(max_grid_size_v);
        DistributionMapping dm(ba);

        auto const dx = geom.CellSizeArray();
        auto const dxinv = geom.InvCellSizeArray();
        Real const Lx = prob_hi[0];
        Real const Ly = prob_hi[1];
        Real const Lz = prob_hi[2];
        Real const nz = Real(n_cell[2]);

        // Nodal heights, with ghost nodes from the analytic mapping.
        MultiFab zp(amrex::convert(ba,IntVect(1)), dm, 1, 1);
        for (MFIter mfi(zp); mfi.isValid(); ++mfi) {
            auto const& z = zp.array(mfi);
            ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real x = i*dx[0];
                Real y = j*dx[1];
                Real s = Real(k)/nz;
                Real hx = Real(0.5) - Real(0.5)*std::cos(Real(2.)*Math::pi<Real>()*(x/Lx-hill_shift));
                Real hy = Real(0.5) - Real(0.5)*std::cos(Real(2.)*Math::pi<Real>()*(y/Ly-hill_shift));
                Real h = hill_height*hx*hy;
                Real zeta = (stretch == Real(0.)) ? s*Lz
                    : Lz*(std::exp(stretch*s)-Real(1.))/(std::exp(stretch)-Real(1.));
                z(i,j,k) = h + zeta*(Lz-h)/Lz;
            });
        }

        // Area and volume factors as in ERF.
        Array<MultiFab,AMREX_SPACEDIM> area;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            area[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 0);
        }
        area[2].setVal(Real(1.0));
        MultiFab detJ(ba, dm, 1, 0);
        for (MFIter mfi(detJ); mfi.isValid(); ++mfi) {
            auto const& z = zp.const_array(mfi);
            auto const& axa = area[0].array(mfi);
            auto const& aya = area[1].array(mfi);
            auto const& ja = detJ.array(mfi);
            ParallelFor(mfi.nodaltilebox(0), mfi.nodaltilebox(1), mfi.validbox(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                axa(i,j,k) = Real(0.5)*dxinv[2]*(z(i,j,k+1)+z(i,j+1,k+1)-z(i,j,k)-z(i,j+1,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                aya(i,j,k) = Real(0.5)*dxinv[2]*(z(i,j,k+1)+z(i+1,j,k+1)-z(i,j,k)-z(i+1,j,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                ja(i,j,k) = Real(0.25)*dxinv[2]*
                    ( z(i,j,k+1) + z(i+1,j,k+1) + z(i,j+1,k+1) + z(i+1,j+1,k+1)
                     -z(i,j,k  ) - z(i+1,j,k  ) - z(i,j+1,k  ) - z(i+1,j+1,k  ) );
            });
        }
        // Smooth periodic perturbation, so that faces shared by boxes agree.
        if (area_noise > Real(0.)) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                for (MFIter mfi(area[idim]); mfi.isValid(); ++mfi) {
                    auto const& a = area[idim].array(mfi);
                    ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        Real x = (i+Real(0.5)*(idim!=0))*dx[0]/Lx;
                        Real y = (j+Real(0.5)*(idim!=1))*dx[1]/Ly;
                        Real z = (k+Real(0.5)*(idim!=2))/nz;
                        a(i,j,k) *= Real(1.0) + area_noise
                            * std::sin(Real(2.)*Math::pi<Real>()*(x+Real(0.3)*idim))
                            * std::cos(Real(2.)*Math::pi<Real>()*(Real(2.)*y-Real(0.1)))
                            * std::cos(Real(3.)*z);
                    });
                }
            }
        }

        LPInfo info;
        info.setAgglomeration(agglomeration).setConsolidation(consolidation)
            .setMaxCoarseningLevel(max_coarsening_level);
        if (hidden_direction >= 0) { info.setHiddenDirection(hidden_direction); }

        MLTerrainPoisson linop({geom}, {ba}, {dm}, info);
        if (zsplit_solver == "column") {
            linop.setZSplitSolver(MLTerrainPoisson::ZSplitSolver::Column);
        } else if (zsplit_solver != "spike") {
            amrex::Abort("Unsupported zsplit_solver " + zsplit_solver);
        }
        linop.setDomainBC(lobc, hibc);
        linop.setLevelBC(0, nullptr);
        linop.setZPhys(0, zp);
        linop.setAreas(0, GetArrOfConstPtrs(area));
        linop.setDetJ(0, detJ);

        // Exact solution and its right-hand side, L(phi) per unit volume.
        MultiFab phi_exact(ba, dm, 1, 1);
        FillRandom(phi_exact, 0, 1);
        for (MFIter mfi(phi_exact); mfi.isValid(); ++mfi) {
            auto const& p = phi_exact.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real x = (i+Real(0.5))*dx[0]/Lx;
                Real y = (j+Real(0.5))*dx[1]/Ly;
                Real z = (k+Real(0.5))/nz;
                p(i,j,k) = std::cos(Real(2.)*Math::pi<Real>()*x)
                    *      std::cos(Real(2.)*Math::pi<Real>()*y)
                    *      std::cos(Math::pi<Real>()*z) + Real(0.1)*p(i,j,k);
            });
        }

        MLMG mlmg(linop);
        mlmg.setMaxIter(max_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
        if (bottom_solver == "smoother") {
            mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
        } else if (bottom_solver == "bicgstab") {
            mlmg.setBottomSolver(MLMG::BottomSolver::bicgstab);
        } else {
            amrex::Abort("Unsupported bottom_solver " + bottom_solver);
        }

        MultiFab rhs(ba, dm, 1, 0);
        mlmg.apply({&rhs}, {&phi_exact});
        MultiFab::Divide(rhs, detJ, 0, 0, 1, 0);

        bool const singular = linop.isSingular(0);

        auto check = [&] (MultiFab& phi, std::string const& name)
        {
            if (singular) {
                Real const shift = (phi.sum(0) - phi_exact.sum(0)) / Real(domain.numPts());
                phi.plus(-shift, 0, 1);
            }
            MultiFab::Subtract(phi, phi_exact, 0, 0, 1, 0);
            Real const err = phi.norminf(0) / phi_exact.norminf(0);
            amrex::Print() << name << ": relative max error " << err << "\n";
            if (!(err < check_tol)) {
                amrex::Abort(name + " failed");
            }
        };

        if (solver == "mlmg" || solver == "both") {
            MultiFab phi(ba, dm, 1, 1);
            phi.setVal(Real(0.0));
            mlmg.solve({&phi}, {&rhs}, reltol, Real(0.0));
            amrex::Print() << "MLMG iterations: " << mlmg.getNumIters() << "\n";
            check(phi, "MLMG");
        }

        if (solver == "gmres" || solver == "both") {
            MultiFab phi(ba, dm, 1, 1);
            phi.setVal(Real(0.0));
            GMRESMLMG gmsolver(mlmg);
            gmsolver.setVerbose(verbose);
            gmsolver.setMaxIters(max_iter);
            gmsolver.solve(phi, rhs, reltol, Real(0.0));
            amrex::Print() << "GMRES iterations: " << gmsolver.getNumIters() << "\n";
            check(phi, "GMRESMLMG");
        }
    }
    amrex::Finalize();
}
