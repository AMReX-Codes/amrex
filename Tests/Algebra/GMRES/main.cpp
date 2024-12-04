#include <AMReX_GMRES_MV.H>

#include <AMReX.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        Box domain(IntVect(0),IntVect(15));
        Long n = domain.numPts();
        AlgVector<Real> xvec(n);
        AlgVector<Real> bvec(xvec.partition());
        AlgVector<Real> exact(xvec.partition());

        Real a = Real(1.e-6);
        Real dx = Real(2)*amrex::Math::pi<Real>()/Real(domain.length(0));

        // The system is a * phi - del dot grad phi.
        // Where phi = sin^5(x)*sin^5(y)*sin^5(z)

        BoxIndexer box_indexer(domain);

        // Initialzie bvec
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
                auto x = (cell[0]+Real(0.5))*dx;
                auto y = (cell[1]+Real(0.5))*dx;
                auto phi0 = Math::powi<5>(std::sin(x)*std::sin(y));
                auto phixm = Math::powi<5>(std::sin(x-dx)*std::sin(y));
                auto phixp = Math::powi<5>(std::sin(x+dx)*std::sin(y));
                auto phiym = Math::powi<5>(std::sin(x)*std::sin(y-dx));
                auto phiyp = Math::powi<5>(std::sin(x)*std::sin(y+dx));
                rhs[lrow] = a*phi0 + (Real(4)*phi0-phixm-phixp-phiym-phiyp) / (dx*dx);
#else
                auto x = (cell[0]+Real(0.5))*dx;
                auto y = (cell[1]+Real(0.5))*dx;
                auto z = (cell[2]+Real(0.5))*dx;
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

        // Initial guess
        xvec.setVal(0);

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
        mat.setVal(set_stencil);

        GMRES_MV<Real> gmres(&mat);
        gmres.setPrecond(JacobiSmoother<Real>(&mat));
        gmres.setVerbose(2);

        auto eps = (sizeof(Real) == 4) ? Real(1.e-5) : Real (1.e-12);
        gmres.solve(xvec, bvec, eps, Real(0.0));

        // Check the solution
        amrex::Axpy(xvec, Real(-1.0), exact);
        auto error = xvec.norminf();
        amrex::Print() << " Max norm error: " << error << "\n";
        AMREX_ALWAYS_ASSERT(error*10 < eps);
    }
    amrex::Finalize();
}
