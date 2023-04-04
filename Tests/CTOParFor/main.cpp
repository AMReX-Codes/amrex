#include <AMReX.H>
#include <AMReX_IArrayBox.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
#if (__cplusplus >= 201703L)
    {
        enum A_options: int {
            A0 = 0, A1
        };

        enum B_options: int {
            B0 = 0, B1, B2
        };

        Box box(IntVect(0),IntVect(7));
        IArrayBox fab(box,2);
        fab.setVal<RunOn::Device>(-10);

        auto const& arr = fab.array();

        for (int ia = 0; ia < 2; ++ia) {
            for (int ib = 0; ib < 3; ++ib) {
                ParallelFor(TypeList<CompileTimeOptions<A0,A1>,
                                     CompileTimeOptions<B0,B1,B2>>{},
                            {ia, ib},
                            box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                       auto A_control,
                                                       auto B_control)
                {
                    auto const& larr = arr;
                    int a, b;
                    if constexpr (A_control.value == 0) {
                        a = 0;
                    } else if constexpr (A_control.value == 1) {
                        a = 1;
                    } else {
                        a = -1;
                    }
                    if constexpr (B_control.value == 0) {
                        b = 0;
                    } else if constexpr (B_control.value == 1) {
                        b = 1;
                    } else if constexpr (B_control.value == 2) {
                        b = 2;
                    } else if constexpr (B_control.value == 3) {
                        b = 3;
                    }
                    larr(i,j,k) = a*10 + b;
                });

                auto s = fab.sum<RunOn::Device>(0);
                AMREX_ALWAYS_ASSERT(s == box.numPts()*(ia*10+ib));
            }
        }
    }
#else
    amrex::Print() << "This test requires C++17." << std::endl;
#endif
    amrex::Finalize();
}
