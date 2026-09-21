#include <AMReX_BiCGStab_MV.H>
#include <AMReX_PCG_MV.H>
#include <AMReX_Smoother_MV.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <cmath>
#include <limits>

using namespace amrex;

// Periodic a*phi - lap(phi) with phi = prod sin^5, solved with BiCGStab and
// PCG, each with the Jacobi and the Chebyshev preconditioner.
int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        int n_cell = 16;
        ParmParse pp;
        pp.query("n_cell", n_cell);
        Box domain(IntVect(0),IntVect(n_cell-1));
        Long n = domain.numPts();
        AlgVector<Real> xvec(n);
        AlgVector<Real> bvec(xvec.partition());
        AlgVector<Real> exact(xvec.partition());

        Real a = Real(1);
        Real dx = Real(2)*amrex::Math::pi<Real>()/Real(domain.length(0));

        BoxIndexer box_indexer(domain);

        {
            auto* rhs = bvec.data();
            auto* phi = exact.data();
            auto nrows = bvec.numLocalRows();
            auto ib = bvec.globalBegin();
            ParallelFor(nrows, [=] AMREX_GPU_DEVICE (Long lrow)
            {
                auto row = lrow + ib; // global row index
                IntVect cell = box_indexer.intVect(row);
#if (AMREX_SPACEDIM == 1)
                auto x = (cell[0]+Real(0.5))*dx;
                auto phi0 = Math::powi<5>(std::sin(x));
                auto phixm = Math::powi<5>(std::sin(x-dx));
                auto phixp = Math::powi<5>(std::sin(x+dx));
                rhs[lrow] = a*phi0 + (Real(2)*phi0-phixm-phixp) / (dx*dx);
#elif (AMREX_SPACEDIM == 2)
                auto x = (Real(cell[0])+Real(0.5))*dx;
                auto y = (Real(cell[1])+Real(0.5))*dx;
                auto phi0 = Math::powi<5>(std::sin(x)*std::sin(y));
                auto phixm = Math::powi<5>(std::sin(x-dx)*std::sin(y));
                auto phixp = Math::powi<5>(std::sin(x+dx)*std::sin(y));
                auto phiym = Math::powi<5>(std::sin(x)*std::sin(y-dx));
                auto phiyp = Math::powi<5>(std::sin(x)*std::sin(y+dx));
                rhs[lrow] = a*phi0 + (Real(4)*phi0-phixm-phixp-phiym-phiyp) / (dx*dx);
#else
                auto x = (Real(cell[0])+Real(0.5))*dx;
                auto y = (Real(cell[1])+Real(0.5))*dx;
                auto z = (Real(cell[2])+Real(0.5))*dx;
                auto phi0 = Math::powi<5>(std::sin(x)*std::sin(y)*std::sin(z));
                auto phixm = Math::powi<5>(std::sin(x-dx)*std::sin(y)*std::sin(z));
                auto phixp = Math::powi<5>(std::sin(x+dx)*std::sin(y)*std::sin(z));
                auto phiym = Math::powi<5>(std::sin(x)*std::sin(y-dx)*std::sin(z));
                auto phiyp = Math::powi<5>(std::sin(x)*std::sin(y+dx)*std::sin(z));
                auto phizm = Math::powi<5>(std::sin(x)*std::sin(y)*std::sin(z-dx));
                auto phizp = Math::powi<5>(std::sin(x)*std::sin(y)*std::sin(z+dx));
                rhs[lrow] = a*phi0 + (Real(6)*phi0-phixm-phixp-phiym-phiyp-phizm-phizp) / (dx*dx);
#endif
                phi[lrow] = phi0;
            });
        }

        // cross stencil w/ periodic boundaries
        auto set_stencil = [=] AMREX_GPU_DEVICE (Long row, Long* col, Real* val)
        {
            IntVect cell = box_indexer.intVect(row);
            int i = 0;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                IntVect cell2 = cell;
                if (cell[idim] == domain.smallEnd(idim)) {
                    cell2[idim] = domain.bigEnd(idim);
                } else {
                    cell2[idim] = cell[idim] - 1;
                }
                Long row2 = domain.index(cell2);
                col[i] = row2;
                val[i] = Real(-1.0)/(dx*dx);
                ++i;

                if (cell[idim] == domain.bigEnd(idim)) {
                    cell2[idim] = domain.smallEnd(idim);
                } else {
                    cell2[idim] = cell[idim] + 1;
                }
                row2 = domain.index(cell2);
                col[i] = row2;
                val[i] = Real(-1.0)/(dx*dx);
                ++i;
            }
            col[i] = row;
            val[i] = Real(2*AMREX_SPACEDIM)/(dx*dx) + a;
        };

        int num_non_zeros = 2*AMREX_SPACEDIM+1;
        SpMatrix<Real> mat(xvec.partition(), num_non_zeros);
        mat.setVal(set_stencil, CsrSorted{false});

        auto eps = (sizeof(Real) == 4) ? Real(1.e-5) : Real(1.e-10);
        auto check = [&] (char const* name)
        {
            amrex::Axpy(xvec, Real(-1.0), exact);
            auto error = xvec.norminf();
            amrex::Print() << name << ": max norm error " << error << "\n";
            AMREX_ALWAYS_ASSERT(error < Real(1.e3)*eps);
        };

        for (int pc = 0; pc < 2; ++pc) {
            using PC = std::function<void(AlgVector<Real>&, AlgVector<Real> const&)>;
            PC precond = (pc == 0) ? PC(JacobiSmoother<Real>(&mat, true))     // l1-Jacobi
                                   : PC(ChebyshevSmoother<Real>(&mat));
            char const* pcname = (pc == 0) ? "l1-Jacobi" : "Chebyshev";

            // l1-Jacobi: default zero guess, so a NaN x must be ignored.
            // Chebyshev: nonzero initial guess.
            xvec.setVal((pc == 0) ? std::numeric_limits<Real>::quiet_NaN() : Real(1));
            BiCGStab_MV<Real> bicgstab(&mat);
            bicgstab.getSolver().setInitialGuessNonzero(pc == 1);
            bicgstab.setPrecond(precond);
            bicgstab.setVerbose(1);
            bicgstab.solve(xvec, bvec, eps, Real(0.0));
            AMREX_ALWAYS_ASSERT(bicgstab.getSolver().getStatus() == 0);
            check((std::string("BiCGStab/") + pcname).c_str());

            xvec.setVal((pc == 0) ? std::numeric_limits<Real>::quiet_NaN() : Real(1));
            PCG_MV<Real> pcg(&mat);
            pcg.getSolver().setInitialGuessNonzero(pc == 1);
            pcg.setPrecond(precond);
            pcg.setVerbose(1);
            pcg.solve(xvec, bvec, eps, Real(0.0));
            AMREX_ALWAYS_ASSERT(pcg.getSolver().getStatus() == 0);
            check((std::string("PCG/") + pcname).c_str());
        }

        // Starting from the solution must take no iterations.
        auto const atol = std::sqrt(eps) * bvec.norm2();
        {
            xvec.copyAsync(exact);
            BiCGStab_MV<Real> bicgstab(&mat);
            bicgstab.getSolver().setInitialGuessNonzero(true);
            bicgstab.solve(xvec, bvec, Real(0.0), atol);
            AMREX_ALWAYS_ASSERT(bicgstab.getSolver().getStatus() == 0 &&
                                bicgstab.getSolver().getNumIters() == 0);
        }
        {
            xvec.copyAsync(exact);
            PCG_MV<Real> pcg(&mat);
            pcg.getSolver().setInitialGuessNonzero(true);
            pcg.solve(xvec, bvec, Real(0.0), atol);
            AMREX_ALWAYS_ASSERT(pcg.getSolver().getStatus() == 0 &&
                                pcg.getSolver().getNumIters() == 0);
        }
    }
    amrex::Finalize();
}
