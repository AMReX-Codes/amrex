#include <AMReX.H>
#include <AMReX_IArrayBox.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
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

        // Test the For overloads with compile time options. Unlike
        // ParallelFor, For makes no promise that iterations are independent,
        // so exercise it with the pattern it exists for: all iterations
        // accumulate into the same memory location. Gpu::Atomic::AddNoRet
        // is an atomic update on GPU and a plain += on CPU, which is only
        // safe because For carries no SIMD pragma.
        Box accbox(IntVect(0),IntVect(0));
        IArrayBox accfab(accbox,1);
        auto const& accarr = accfab.array();

        for (int ia = 0; ia < 2; ++ia) {
            for (int ib = 0; ib < 3; ++ib) {
                accfab.setVal<RunOn::Device>(0);
                For(TypeList<CompileTimeOptions<A0,A1>,
                             CompileTimeOptions<B0,B1,B2>>{},
                    {ia, ib},
                    box, [=] AMREX_GPU_DEVICE (int, int, int,
                                               auto A_control,
                                               auto B_control)
                {
                    auto const& lacc = accarr;
                    int a, b;
                    if constexpr (A_control.value == 0) {
                        a = 0;
                    } else {
                        a = 1;
                    }
                    if constexpr (B_control.value == 0) {
                        b = 0;
                    } else if constexpr (B_control.value == 1) {
                        b = 1;
                    } else {
                        b = 2;
                    }
                    Gpu::Atomic::AddNoRet(&lacc(0,0,0), a*10 + b);
                });

                auto s = accfab.sum<RunOn::Device>(0);
                AMREX_ALWAYS_ASSERT(s == box.numPts()*(ia*10+ib));

                accfab.setVal<RunOn::Device>(0);
                For(TypeList<CompileTimeOptions<A0,A1>,
                             CompileTimeOptions<B0,B1,B2>>{},
                    {ia, ib},
                    box.numPts(), [=] AMREX_GPU_DEVICE (Long,
                                                        auto A_control,
                                                        auto B_control)
                {
                    auto const& lacc = accarr;
                    int a, b;
                    if constexpr (A_control.value == 0) {
                        a = 0;
                    } else {
                        a = 1;
                    }
                    if constexpr (B_control.value == 0) {
                        b = 0;
                    } else if constexpr (B_control.value == 1) {
                        b = 1;
                    } else {
                        b = 2;
                    }
                    Gpu::Atomic::AddNoRet(&lacc(0,0,0), a*10 + b);
                });

                s = accfab.sum<RunOn::Device>(0);
                AMREX_ALWAYS_ASSERT(s == box.numPts()*(ia*10+ib));
            }
        }
    }
    amrex::Finalize();
}
