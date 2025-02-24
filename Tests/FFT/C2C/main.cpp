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

        AMREX_D_TERM(int n_cell_x = 128;,
                     int n_cell_y = 32;,
                     int n_cell_z = 64);

        AMREX_D_TERM(int max_grid_size_x = 16;,
                     int max_grid_size_y = 32;,
                     int max_grid_size_z = 32);

        AMREX_D_TERM(Real prob_lo_x = 0.;,
                     Real prob_lo_y = 0.;,
                     Real prob_lo_z = 0.);
        AMREX_D_TERM(Real prob_hi_x = 1.;,
                     Real prob_hi_y = 1.;,
                     Real prob_hi_z = 1.);

        {
            ParmParse pp;
            AMREX_D_TERM(pp.query("n_cell_x", n_cell_x);,
                         pp.query("n_cell_y", n_cell_y);,
                         pp.query("n_cell_z", n_cell_z));
            AMREX_D_TERM(pp.query("max_grid_size_x", max_grid_size_x);,
                         pp.query("max_grid_size_y", max_grid_size_y);,
                         pp.query("max_grid_size_z", max_grid_size_z));
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

        cMultiFab mf(ba,dm,1,0);
        auto const& ma = mf.arrays();
        ParallelFor(mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            AMREX_D_TERM(Real x = (i+0.5_rt) * dx[0] - 0.5_rt;,
                         Real y = (j+0.5_rt) * dx[1] - 0.5_rt;,
                         Real z = (k+0.5_rt) * dx[2] - 0.5_rt);
            auto tmp = std::exp(-10._rt*
                (AMREX_D_TERM(x*x*1.05_rt, + y*y*0.90_rt, + z*z)));
            ma[b](i,j,k) = GpuComplex<Real>{tmp, 1._rt-tmp};
        });

        cMultiFab mf2(ba,dm,1,0);

        auto scaling = Real(1) / Real(geom.Domain().d_numPts());

        {
            FFT::Info info{};
            info.setDomainStrategy(FFT::DomainStrategy::slab);
            FFT::C2C<Real,FFT::Direction::both> c2c(geom.Domain(), info);
            c2c.forward(mf);
            c2c.backward(mf2);
        }

        {
            MultiFab errmf(ba,dm,1,0);
            auto const& errma = errmf.arrays();

            auto const& ma2 = mf2.const_arrays();
            ParallelFor(mf2, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto err = ma[b](i,j,k) - ma2[b](i,j,k)*scaling;
                errma[b](i,j,k) = amrex::norm(err);
            });

            auto error = errmf.norminf();
            amrex::Print() << "  Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 1.e-6f;
#else
            auto eps = 1.e-13;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }

        mf2.setVal(0);

        {
            FFT::Info info{};
            info.setDomainStrategy(FFT::DomainStrategy::pencil);
            FFT::C2C<Real,FFT::Direction::both> c2c(geom.Domain(), info);
            c2c.forward(mf);
            c2c.backward(mf2);
        }

        {
            MultiFab errmf(ba,dm,1,0);
            auto const& errma = errmf.arrays();

            auto const& ma2 = mf2.const_arrays();
            ParallelFor(mf2, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto err = ma[b](i,j,k) - ma2[b](i,j,k)*scaling;
                errma[b](i,j,k) = amrex::norm(err);
            });

            auto error = errmf.norminf();
            amrex::Print() << "  Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 1.e-6f;
#else
            auto eps = 1.e-13;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }

        {
            auto sba = amrex::decompose(domain, ParallelDescriptor::NProcs());
            DistributionMapping sdm{sba};
            cMultiFab smf(sba,sdm,1,0);
            smf.ParallelCopy(mf);

            auto domain_size = domain.length().toArray();
            FFT::C2C<Real,FFT::Direction::both> c2c(domain_size);

            GpuComplex<Real>* pio = nullptr;
            std::array<int,AMREX_SPACEDIM> local_start{AMREX_D_DECL(0,0,0)};
            std::array<int,AMREX_SPACEDIM> local_size{AMREX_D_DECL(0,0,0)};
            auto const& imap = smf.IndexArray();
            AMREX_ALWAYS_ASSERT(imap.size() <= 1);
            if (imap.size() == 1) {
                pio = smf[imap[0]].dataPtr();
                auto const& box = sba[imap[0]];
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    local_start[idim] = box.smallEnd(idim);
                    local_size[idim] = box.length(idim);
                }
            }
            c2c.setLocalDomain(local_start, local_size);

            auto const& [sp_start, sp_size] = c2c.getLocalSpectralDomain();
            auto npts = AMREX_D_TERM(std::size_t(sp_size[0]),
                                    *std::size_t(sp_size[1]),
                                    *std::size_t(sp_size[2]));
            Gpu::DeviceVector<GpuComplex<Real>> spv(npts);

            c2c.forward(pio, spv.data());
            c2c.backward(spv.data(), pio);

            amrex::Scale(smf, -scaling, 0, 1, 0);
            smf.ParallelAdd(mf);

            MultiFab errmf(sba,sdm,1,0);
            auto const& errma = errmf.arrays();

            auto const& sma = smf.const_arrays();
            ParallelFor(smf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                errma[b](i,j,k) = amrex::norm(sma[b](i,j,k));
            });

            auto error = errmf.norminf();
            amrex::Print() << "  Expected to be close to zero: " << error << "\n";
#ifdef AMREX_USE_FLOAT
            auto eps = 1.e-6f;
#else
            auto eps = 1.e-13;
#endif
            AMREX_ALWAYS_ASSERT(error < eps);
        }
    }
    amrex::Finalize();
}
