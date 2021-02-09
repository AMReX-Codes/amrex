#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Tuple.H>
#include <AMReX_PODVector.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParticleDataLayout.H>

using namespace amrex;

/**
   /brief Particle type is defined as a GpuTuple.
   x, y, z, id, and cpu methods are required,
   but can be implemented freely, i.e. using
   a cell + offset representation of position.
 */
template <typename... T>
struct Particle : public amrex::GpuTuple<T...>
{
    using amrex::GpuTuple<T...>::GpuTuple;
    auto& x () { return amrex::get<0>(*this); }
    auto& y () { return amrex::get<1>(*this); }
    auto& z () { return amrex::get<2>(*this); }
    auto& id () { return amrex::get<3>(*this); }
    auto& cpu () { return amrex::get<4>(*this); }
};

template <DataLayout DataLayoutTag>
void testLayout ()
{
    using ParticleType = Particle<double, double, double, int, int>;
    using ParticleTileType = ParticleTile<amrex::Gpu::DeviceVector, ParticleType, DataLayoutTag>;

    ParticleTileType ptile(1000);
    auto particles = ptile.get_particle_data();

    amrex::ParallelFor( particles.size(), [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        auto&& p = particles[i];
	p.x() = 7.0;
        p.y() = 8.0;
	p.z() = 9.0;
    });

    amrex::ParallelFor( particles.size(), [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        auto&& p = particles[i];
        p.x() -= 6.0;
        p.y() -= 7.0;
        p.z() -= 8.0;
    });

    amrex::ParallelFor( particles.size(), [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        auto&& p = particles[i];
        AMREX_ALWAYS_ASSERT(p.x() == p.y() == p.z() == 1.0);
    });
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Print() << "Running data layout test \n";
    testLayout<DataLayout::AoS>();
    testLayout<DataLayout::SoA>();

    amrex::Finalize();
}
