#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Algorithm.H>

void testCLZ ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Print() << "Running CLZ test. \n";
    testCLZ();

    amrex::Finalize();
}

void testCLZ ()
{
    AMREX_ALWAYS_ASSERT(amrex::clz(std::uint8_t(10) ) == 4 );
    AMREX_ALWAYS_ASSERT(amrex::clz(std::uint16_t(10)) == 12);
    AMREX_ALWAYS_ASSERT(amrex::clz(std::uint32_t(10)) == 28);
    AMREX_ALWAYS_ASSERT(amrex::clz(std::uint64_t(10)) == 60);

    AMREX_ALWAYS_ASSERT(amrex::clz(std::uint8_t (1 << 7  )) == 0);
    AMREX_ALWAYS_ASSERT(amrex::clz(std::uint16_t(1 << 15 )) == 0);
    AMREX_ALWAYS_ASSERT(amrex::clz(std::uint32_t(1 << 31 )) == 0);
    AMREX_ALWAYS_ASSERT(amrex::clz(std::uint64_t(1L << 63)) == 0);

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int /*i*/) noexcept {
            AMREX_ALWAYS_ASSERT(amrex::clz(std::uint8_t(10) ) == 4 );
            AMREX_ALWAYS_ASSERT(amrex::clz(std::uint16_t(10)) == 12);
            AMREX_ALWAYS_ASSERT(amrex::clz(std::uint32_t(10)) == 28);
            AMREX_ALWAYS_ASSERT(amrex::clz(std::uint64_t(10)) == 60);

            AMREX_ALWAYS_ASSERT(amrex::clz(std::uint8_t (1 << 7  )) == 0);
            AMREX_ALWAYS_ASSERT(amrex::clz(std::uint16_t(1 << 15 )) == 0);
            AMREX_ALWAYS_ASSERT(amrex::clz(std::uint32_t(1 << 31 )) == 0);
            AMREX_ALWAYS_ASSERT(amrex::clz(std::uint64_t(1L << 63)) == 0);
        });
}
