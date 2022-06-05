#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Random.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <cmath>

using namespace amrex;

void RandomNumGen();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    RandomNumGen();
    amrex::Finalize();
}

void RandomNumGen ()
{
    BL_PROFILE("main");

    ParmParse pp;

    int Ndraw = 1000000;

    pp.query("num_draw", Ndraw);

    Gpu::DeviceVector<Real> x_d(Ndraw);
    Gpu::DeviceVector<Real> y_d(Ndraw);
    Gpu::DeviceVector<Real> z_d(Ndraw);

    {
        BL_PROFILE("Draw");

        auto x_d_ptr = x_d.dataPtr();
        auto y_d_ptr = y_d.dataPtr();
        auto z_d_ptr = z_d.dataPtr();

        amrex::ParallelForRNG(Ndraw,
        [=] AMREX_GPU_DEVICE (int i, RandomEngine const& engine) noexcept
        {
            x_d_ptr[i] = amrex::Random(engine);
            y_d_ptr[i] = amrex::Random(engine);
            z_d_ptr[i] = amrex::Random(engine);
        });

        Gpu::streamSynchronize();
    }

    std::vector<Real> x_h(Ndraw);
    std::vector<Real> y_h(Ndraw);
    std::vector<Real> z_h(Ndraw);
    Gpu::copyAsync(Gpu::deviceToHost, x_d.begin(), x_d.end(), x_h.begin());
    Gpu::copyAsync(Gpu::deviceToHost, y_d.begin(), y_d.end(), y_h.begin());
    Gpu::copyAsync(Gpu::deviceToHost, z_d.begin(), z_d.end(), z_h.begin());
    Gpu::streamSynchronize();

    Real xmean=0., ymean=0., zmean=0., xvar=0., yvar=0., zvar=0.;
    for (int i = 0; i < Ndraw; ++i) {
        xmean += x_h[i];
        ymean += y_h[i];
        zmean += z_h[i];
        xvar += std::pow(x_h[i]-0.5,2);
        yvar += std::pow(y_h[i]-0.5,2);
        zvar += std::pow(z_h[i]-0.5,2);
    }
    xmean /= Ndraw;
    ymean /= Ndraw;
    zmean /= Ndraw;
    xvar /= Ndraw;
    yvar /= Ndraw;
    zvar /= Ndraw;
    amrex::Print() << "\n  Means: " << xmean << ", " << ymean << ", " << zmean
                   << "  Variances: " << xvar << ", " << yvar << ", " << zvar
                   << std::endl;
}
