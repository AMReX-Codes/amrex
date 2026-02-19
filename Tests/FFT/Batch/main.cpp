#include <AMReX_FFT.H> // Put this at the top for testing

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");

        AMREX_D_TERM(int n_cell_x = 64;,
                     int n_cell_y = 16;,
                     int n_cell_z = 32);

        AMREX_D_TERM(int max_grid_size_x = 32;,
                     int max_grid_size_y = 16;,
                     int max_grid_size_z = 16);

        AMREX_D_TERM(Real prob_lo_x = 0.;,
                     Real prob_lo_y = 0.;,
                     Real prob_lo_z = 0.);
        AMREX_D_TERM(Real prob_hi_x = 1.;,
                     Real prob_hi_y = 1.;,
                     Real prob_hi_z = 1.);

        int batch_size = 4;

        {
            ParmParse pp;
            AMREX_D_TERM(pp.query("n_cell_x", n_cell_x);,
                         pp.query("n_cell_y", n_cell_y);,
                         pp.query("n_cell_z", n_cell_z));
            AMREX_D_TERM(pp.query("max_grid_size_x", max_grid_size_x);,
                         pp.query("max_grid_size_y", max_grid_size_y);,
                         pp.query("max_grid_size_z", max_grid_size_z));
            pp.query("batch_size", batch_size);
        }

        Box domain(IntVect(0),IntVect(AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)));
        BoxArray ba(domain);
        ba.maxSize(IntVect(AMREX_D_DECL(max_grid_size_x,
                                        max_grid_size_y,
                                        max_grid_size_z)));
        DistributionMapping dm(ba);

        Geometry geom;
        {
            geom.define(domain,
                        RealBox(AMREX_D_DECL(prob_lo_x,prob_lo_y,prob_lo_z),
                                AMREX_D_DECL(prob_hi_x,prob_hi_y,prob_hi_z)),
                        CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        }
        auto const& dx = geom.CellSizeArray();

        MultiFab mf(ba,dm,batch_size,0);
        auto const& ma = mf.arrays();
        ParallelFor(mf, IntVect(0), batch_size,
                    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
        {
            AMREX_D_TERM(Real x = (i+0.5_rt) * dx[0] - 0.5_rt;,
                         Real y = (j+0.5_rt) * dx[1] - 0.5_rt;,
                         Real z = (k+0.5_rt) * dx[2] - 0.5_rt);
            ma[b](i,j,k,n) = std::exp(-10._rt*
                (AMREX_D_TERM(x*x*1.05_rt, + y*y*0.90_rt, + z*z))) + Real(n);
        });

        MultiFab mf2(ba,dm,batch_size,0);

        auto scaling = Real(1) / Real(geom.Domain().d_numPts());

        cMultiFab cmf;

        // forward
        {
            FFT::Info info{};
            info.setDomainStrategy(FFT::DomainStrategy::pencil);
            info.setBatchSize(batch_size);
            FFT::R2C<Real,FFT::Direction::forward> r2c(geom.Domain(), info);
            auto const& [cba, cdm] = r2c.getSpectralDataLayout();
            cmf.define(cba, cdm, batch_size, 0);
            r2c.forward(mf,cmf);
        }

        // backward
        {
            FFT::Info info{};
            info.setDomainStrategy(FFT::DomainStrategy::slab);
            info.setBatchSize(batch_size);
            FFT::R2C<Real,FFT::Direction::backward> r2c(geom.Domain(), info);
            r2c.backward(cmf,mf2);
        }

        {
            auto const& ma2 = mf2.arrays();
            ParallelFor(mf2, IntVect(0), batch_size,
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
            {
                ma2[b](i,j,k,n) = ma[b](i,j,k,n) - ma2[b](i,j,k,n)*scaling;
            });

            auto error = mf2.norminf(0, batch_size, IntVect(0));
            amrex::Print() << "  Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 3.e-6F;
#else
            auto eps = 1.e-13;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }

        {
            FFT::R2C<Real,FFT::Direction::forward> r2c(geom.Domain());
            cMultiFab cmf2(cmf.boxArray(), cmf.DistributionMap(), 2, 0);
            MultiFab errmf(cmf.boxArray(), cmf.DistributionMap(), cmf.nComp(), 0);
            for (int icomp = 0; icomp < batch_size; ++icomp) {
                r2c.forward(mf, cmf2, icomp, 1);
                auto const& cma = cmf.const_arrays();
                auto const& cma2 = cmf2.const_arrays();
                auto const& ema = errmf.arrays();
                ParallelFor(errmf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                {
                    auto c = cma[b](i,j,k,icomp) - cma2[b](i,j,k,1);
                    ema[b](i,j,k,icomp) = amrex::norm(c);
                });
                Gpu::streamSynchronize();
            }

            auto error = errmf.norminf(0, batch_size, IntVect(0));
            amrex::Print() << "  Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 3.e-6F;
#else
            auto eps = 1.e-15;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }

        {
            FFT::R2C<Real,FFT::Direction::backward> r2c(geom.Domain());
            for (int icomp = 0; icomp < batch_size; ++icomp) {
                r2c.backward(cmf, mf2, icomp, icomp);
            }

            auto const& ma2 = mf2.arrays();
            ParallelFor(mf2, IntVect(0), batch_size,
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
            {
                ma2[b](i,j,k,n) = ma[b](i,j,k,n) - ma2[b](i,j,k,n)*scaling;
            });

            auto error = mf2.norminf(0, batch_size, IntVect(0));
            amrex::Print() << "  Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 3.e-6F;
#else
            auto eps = 1.e-13;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }
    }
    amrex::Finalize();
}
