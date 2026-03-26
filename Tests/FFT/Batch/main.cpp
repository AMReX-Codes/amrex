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

#if (AMREX_SPACEDIM >= 2)
#if (AMREX_SPACEDIM == 2)
        constexpr const char* oned_mode_dim_tag = "2D";
#else
        constexpr const char* oned_mode_dim_tag = "3D";
#endif
        {
            MultiFab oned_mf(ba, dm, 1, 0);
            auto const& oa = oned_mf.arrays();
            int nx = geom.Domain().length(0);
            Real two_pi = 2._rt * Math::pi<Real>();
            ParallelFor(oned_mf, IntVect(0), 1,
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
            {
                Real base = 0._rt;
                Real cos_amp = 0._rt;
#if (AMREX_SPACEDIM == 2)
                amrex::ignore_unused(k);
                base = (1._rt + Real(j)) * 0.25_rt;
                cos_amp = 0.05_rt * Real(j+1);
#else
                base = 0.2_rt + 0.05_rt*Real(j+1) + 0.02_rt*Real(k+1);
                cos_amp = 0.01_rt * Real(j+1) * Real(k+1);
#endif
                Real theta = two_pi * Real(i) / Real(nx);
                oa[b](i,j,k,n) = base + cos_amp * std::cos(theta);
            });

            FFT::Info info{};
            info.setOneDMode(true);
            FFT::R2C<Real,FFT::Direction::both> r2c(geom.Domain(), info);
            auto const& [cba, cdm] = r2c.getSpectralDataLayout();
            cMultiFab spec(cba, cdm, 1, 0);
            r2c.forward(oned_mf, spec);

            MultiFab err(spec.boxArray(), spec.DistributionMap(), 1, 0);
            auto const& ca = spec.const_arrays();
            auto const& ea = err.arrays();
            ParallelFor(err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real base = 0._rt;
                Real cos_amp = 0._rt;
#if (AMREX_SPACEDIM == 2)
                amrex::ignore_unused(k);
                base = (1._rt + Real(j)) * 0.25_rt;
                cos_amp = 0.05_rt * Real(j+1);
#else
                base = 0.2_rt + 0.05_rt*Real(j+1) + 0.02_rt*Real(k+1);
                cos_amp = 0.01_rt * Real(j+1) * Real(k+1);
#endif
                Real expected_real = 0._rt;
                if (i == 0) {
                    expected_real = base * Real(nx);
                } else if (i == 1) {
                    expected_real = cos_amp * Real(nx) * 0.5_rt;
                }
                auto c = ca[b](i,j,k,0);
                ea[b](i,j,k) = amrex::norm(c - GpuComplex<Real>(expected_real, 0._rt));
            });

#ifdef AMREX_USE_FLOAT
            Real tol = 5.e-6_rt;
#else
            Real tol = 1.e-14_rt;
#endif
            Real err_norm = err.norminf();
            amrex::Print() << "  Expected to be close to zero (" << oned_mode_dim_tag
                           << " R2C one-d-mode): " << err_norm << "\n";
            AMREX_ALWAYS_ASSERT(err_norm < tol);

            MultiFab recon(ba, dm, 1, 0);
            r2c.backward(spec, recon);
            MultiFab back_err(ba, dm, 1, 0);
            auto const back_scaling = r2c.scalingFactor();
            auto const& orig_a = oned_mf.const_arrays();
            auto const& recon_a = recon.const_arrays();
            auto const& back_a = back_err.arrays();
            ParallelFor(back_err, IntVect(0), 1,
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
            {
                Real diff = orig_a[b](i,j,k,n) - recon_a[b](i,j,k,n) * back_scaling;
                back_a[b](i,j,k,n) = amrex::Math::abs(diff);
            });
            Real back_norm = back_err.norminf();
            amrex::Print() << "  Expected to be close to zero (" << oned_mode_dim_tag
                           << " R2C one-d-mode backward): " << back_norm << "\n";
            AMREX_ALWAYS_ASSERT(back_norm < tol);
        }

        {
            cMultiFab cin(ba, dm, 1, 0);
            auto const& cia = cin.arrays();
            int nx = geom.Domain().length(0);
            Real two_pi = 2._rt * Math::pi<Real>();
            ParallelFor(cin, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real theta = two_pi * Real(i) / Real(nx);
                GpuComplex<Real> base(0._rt, 0._rt);
                GpuComplex<Real> cos_amp(0._rt, 0._rt);
#if (AMREX_SPACEDIM == 2)
                amrex::ignore_unused(k);
                Real mag = 0.01_rt * Real(j+1);
                base = GpuComplex<Real>(mag, -0.125_rt*mag);
                cos_amp = GpuComplex<Real>(0.2_rt*mag, 0.05_rt*mag);
#else
                Real mag = 0.01_rt * Real(j+1) + 0.005_rt * Real(k+1);
                base = GpuComplex<Real>(mag, -0.1_rt*mag);
                cos_amp = GpuComplex<Real>(0.15_rt*Real(j+1), 0.05_rt*Real(k+1));
#endif
                cia[b](i,j,k) = base + cos_amp * std::cos(theta);
            });

            FFT::Info info{};
            info.setOneDMode(true);
            FFT::C2C<Real,FFT::Direction::both> c2c(geom.Domain(), info);
            auto const& [cba, cdm] = c2c.getSpectralDataLayout();
            cMultiFab spec(cba, cdm, 1, 0);
            c2c.forward(cin, spec);

            MultiFab err(spec.boxArray(), spec.DistributionMap(), 1, 0);
            auto const& ca = spec.const_arrays();
            auto const& ea = err.arrays();
            ParallelFor(err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                GpuComplex<Real> base(0._rt, 0._rt);
                GpuComplex<Real> cos_amp(0._rt, 0._rt);
#if (AMREX_SPACEDIM == 2)
                amrex::ignore_unused(k);
                Real mag = 0.01_rt * Real(j+1);
                base = GpuComplex<Real>(mag, -0.125_rt*mag);
                cos_amp = GpuComplex<Real>(0.2_rt*mag, 0.05_rt*mag);
#else
                Real mag = 0.01_rt * Real(j+1) + 0.005_rt * Real(k+1);
                base = GpuComplex<Real>(mag, -0.1_rt*mag);
                cos_amp = GpuComplex<Real>(0.15_rt*Real(j+1), 0.05_rt*Real(k+1));
#endif
                GpuComplex<Real> expected(0._rt, 0._rt);
                if (i == 0) {
                    expected = base * Real(nx);
                } else if (i == 1 || i == nx-1) {
                    expected = cos_amp * (Real(nx) * 0.5_rt);
                }
                auto c = ca[b](i,j,k,0);
                ea[b](i,j,k) = amrex::norm(c - expected);
            });

#ifdef AMREX_USE_FLOAT
            Real tol = 5.e-6_rt;
#else
            Real tol = 1.e-14_rt;
#endif
            Real err_norm = err.norminf();
            amrex::Print() << "  Expected to be close to zero (" << oned_mode_dim_tag
                           << " C2C one-d-mode): " << err_norm << "\n";
            AMREX_ALWAYS_ASSERT(err_norm < tol);

            cMultiFab recon(cin.boxArray(), cin.DistributionMap(), 1, 0);
            c2c.backward(spec, recon);

            MultiFab back_err(cin.boxArray(), cin.DistributionMap(), 1, 0);
            auto const back_scaling = c2c.scalingFactor();
            auto const& cin_a = cin.const_arrays();
            auto const& recon_a = recon.const_arrays();
            auto const& back_a = back_err.arrays();
            ParallelFor(back_err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto diff = cin_a[b](i,j,k);
                auto tmp = recon_a[b](i,j,k);
                tmp *= back_scaling;
                diff -= tmp;
                back_a[b](i,j,k) = amrex::norm(diff);
            });
            Real back_norm = back_err.norminf();
            amrex::Print() << "  Expected to be close to zero (" << oned_mode_dim_tag
                           << " C2C one-d-mode backward): " << back_norm << "\n";
            AMREX_ALWAYS_ASSERT(back_norm < tol);
        }
#endif

#if (AMREX_SPACEDIM == 3)
        {
            MultiFab twod_mf(ba, dm, 1, 0);
            auto const& ta = twod_mf.arrays();
            int nx = geom.Domain().length(0);
            int ny = geom.Domain().length(1);
            Real two_pi = 2._rt * Math::pi<Real>();
            ParallelFor(twod_mf, IntVect(0), 1,
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
            {
                Real base = 1._rt + 0.05_rt*Real(k+1);
                Real cos_amp = 0.02_rt * Real(k+1);
                Real theta = two_pi * Real(i) / Real(nx);
                ta[b](i,j,k,n) = base + cos_amp * std::cos(theta);
            });

            FFT::Info info{};
            info.setTwoDMode(true);
            FFT::R2C<Real,FFT::Direction::both> r2c(geom.Domain(), info);
            auto const& [cba, cdm] = r2c.getSpectralDataLayout();
            cMultiFab spec(cba, cdm, 1, 0);
            r2c.forward(twod_mf, spec);

            MultiFab err(spec.boxArray(), spec.DistributionMap(), 1, 0);
            auto const& ca = spec.const_arrays();
            auto const& ea = err.arrays();
            ParallelFor(err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real base = 1._rt + 0.05_rt*Real(k+1);
                Real cos_amp = 0.02_rt * Real(k+1);
                Real expected_real = 0._rt;
                bool x_mode = (i == 1 || i == nx-1) && (j == 0);
                if (i == 0 && j == 0) {
                    expected_real = base * Real(nx*ny);
                } else if (x_mode) {
                    expected_real = cos_amp * Real(nx) * 0.5_rt * Real(ny);
                }
                auto c = ca[b](i,j,k,0);
                ea[b](i,j,k) = amrex::norm(c - GpuComplex<Real>(expected_real, 0._rt));
            });

#ifdef AMREX_USE_FLOAT
            Real tol = 5.e-6_rt;
#else
            Real tol = 1.e-14_rt;
#endif
            Real err_norm = err.norminf();
            amrex::Print() << "  Expected to be close to zero (3D R2C two-d-mode): "
                           << err_norm << "\n";
            AMREX_ALWAYS_ASSERT(err_norm < tol);

            MultiFab recon(ba, dm, 1, 0);
            r2c.backward(spec, recon);
            MultiFab back_err(ba, dm, 1, 0);
            auto const back_scaling = r2c.scalingFactor();
            auto const& orig_a = twod_mf.const_arrays();
            auto const& recon_a = recon.const_arrays();
            auto const& back_a = back_err.arrays();
            ParallelFor(back_err, IntVect(0), 1,
                        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
            {
                Real diff = orig_a[b](i,j,k,n) - recon_a[b](i,j,k,n) * back_scaling;
                back_a[b](i,j,k,n) = amrex::Math::abs(diff);
            });
            Real back_norm = back_err.norminf();
            amrex::Print() << "  Expected to be close to zero (3D R2C two-d-mode backward): "
                           << back_norm << "\n";
            AMREX_ALWAYS_ASSERT(back_norm < tol);
        }

        {
            cMultiFab cin(ba, dm, 1, 0);
            auto const& cia = cin.arrays();
            int nx = geom.Domain().length(0);
            int ny = geom.Domain().length(1);
            Real two_pi = 2._rt * Math::pi<Real>();
            ParallelFor(cin, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real mag = 0.02_rt * Real(k+1);
                GpuComplex<Real> base(mag, 0.1_rt*mag);
                GpuComplex<Real> cos_amp(0.3_rt*mag, -0.15_rt*mag);
                Real theta = two_pi * Real(i) / Real(nx);
                cia[b](i,j,k) = base + cos_amp * std::cos(theta);
            });

            FFT::Info info{};
            info.setTwoDMode(true);
            FFT::C2C<Real,FFT::Direction::both> c2c(geom.Domain(), info);
            auto const& [cba, cdm] = c2c.getSpectralDataLayout();
            cMultiFab spec(cba, cdm, 1, 0);
            c2c.forward(cin, spec);

            MultiFab err(spec.boxArray(), spec.DistributionMap(), 1, 0);
            auto const& ca = spec.const_arrays();
            auto const& ea = err.arrays();
            ParallelFor(err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real mag = 0.02_rt * Real(k+1);
                GpuComplex<Real> base(mag, 0.1_rt*mag);
                GpuComplex<Real> cos_amp(0.3_rt*mag, -0.15_rt*mag);
                GpuComplex<Real> expected(0._rt, 0._rt);
                bool x_mode = (i == 1 || i == nx-1) && (j == 0);
                if (i == 0 && j == 0) {
                    expected = base * Real(nx*ny);
                } else if (x_mode) {
                    expected = cos_amp * (Real(nx) * 0.5_rt * Real(ny));
                }
                auto c = ca[b](i,j,k,0);
                ea[b](i,j,k) = amrex::norm(c - expected);
            });

#ifdef AMREX_USE_FLOAT
            Real tol = 5.e-6_rt;
#else
            Real tol = 1.e-14_rt;
#endif
            Real err_norm = err.norminf();
            amrex::Print() << "  Expected to be close to zero (3D C2C two-d-mode): "
                           << err_norm << "\n";
            AMREX_ALWAYS_ASSERT(err_norm < tol);

            cMultiFab recon(ba, dm, 1, 0);
            c2c.backward(spec, recon);
            MultiFab back_err(ba, dm, 1, 0);
            auto const back_scaling = c2c.scalingFactor();
            auto const& cin_a = cin.const_arrays();
            auto const& recon_a = recon.const_arrays();
            auto const& back_a = back_err.arrays();
            ParallelFor(back_err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto diff = cin_a[b](i,j,k);
                auto tmp = recon_a[b](i,j,k);
                tmp *= back_scaling;
                diff -= tmp;
                back_a[b](i,j,k) = amrex::norm(diff);
            });
            Real back_norm = back_err.norminf();
            amrex::Print() << "  Expected to be close to zero (3D C2C two-d-mode backward): "
                           << back_norm << "\n";
            AMREX_ALWAYS_ASSERT(back_norm < tol);
        }
#endif
    }
    amrex::Finalize();
}
