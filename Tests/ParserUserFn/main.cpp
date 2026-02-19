#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Parser.H>

#include "my_fn.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    double x = 3.0, y = 2.0;
    double expected = x + y +  f1_h(x) +  f2_h(x+y,x-2*y) +  f3_h(x,y,x-y) -  f4_h(x,y,x,7);
    std::string expr("x + y + uf1  (x) + uf2  (x+y,x-2*y) + uf3  (x,y,x-y) - uf4  (x,y,x,7)");

    { // host only
        Parser parser(expr);
        parser.registerVariables({"x","y"});
        {
            std::map<std::string,int> const& ufs = parser.userFunctions();
            for (auto const& [fname, nargs] : ufs) {
                std::cout << "User function: " << fname << "(";
                for (int iarg = 0; iarg < nargs; ++iarg) {
                    std::cout << "double";
                    if (iarg != nargs-1) {
                        std::cout << ",";
                    }
                }
                std::cout << ")\n";
            }
        }
        parser.registerUserFn1("uf1", f1_h, nullptr);
        parser.registerUserFn2("uf2", f2_h, nullptr);
        parser.registerUserFn3("uf3", f3_h, nullptr);
        parser.registerUserFn4("uf4", f4_h, nullptr);

        auto const exe = parser.compile<2>();
        auto const result = exe(x,y);
        AMREX_ALWAYS_ASSERT(result == expected);
        amrex::Print() << "SUCCESS on host\n";
    }

#if !defined(AMREX_USE_GPU) || defined(AMREX_USE_GPU_RDC)

    { // device only
        Parser parser(expr);
        parser.registerVariables({"x","y"});

        auto* fp1 = AMREX_GET_DEVICE_FUNC_PTR(ParserUserFn1, f1_d);
        auto* fp2 = AMREX_GET_DEVICE_FUNC_PTR(ParserUserFn2, f2_d);
        auto* fp3 = AMREX_GET_DEVICE_FUNC_PTR(ParserUserFn3, f3_d);
        auto* fp4 = AMREX_GET_DEVICE_FUNC_PTR(ParserUserFn4, f4_d);

        parser.registerUserFn1("uf1", nullptr, fp1);
        parser.registerUserFn2("uf2", nullptr, fp2);
        parser.registerUserFn3("uf3", nullptr, fp3);
        parser.registerUserFn4("uf4", nullptr, fp4);

        auto const exe = parser.compile<2>();

        Gpu::PinnedVector<double> result(1);
        auto* pr = result.data();
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
        {
            *pr = exe(x,y);
        });
        Gpu::streamSynchronize();
        AMREX_ALWAYS_ASSERT(*pr == expected);
        amrex::Print() << "SUCCESS on device\n";
    }

    { // host and device
        Parser parser(expr);
        parser.registerVariables({"x","y"});

        auto* fp1_d = AMREX_GET_DEVICE_FUNC_PTR(ParserUserFn1, f1_hd);
        auto* fp2_d = AMREX_GET_DEVICE_FUNC_PTR(ParserUserFn2, f2_hd);
        auto* fp3_d = AMREX_GET_DEVICE_FUNC_PTR(ParserUserFn3, f3_hd);
        auto* fp4_d = AMREX_GET_DEVICE_FUNC_PTR(ParserUserFn4, f4_hd);

        parser.registerUserFn1("uf1", f1_hd, fp1_d);
        parser.registerUserFn2("uf2", f2_hd, fp2_d);
        parser.registerUserFn3("uf3", f3_hd, fp3_d);
        parser.registerUserFn4("uf4", f4_hd, fp4_d);

        auto const exe = parser.compile<2>();

        AMREX_ALWAYS_ASSERT(exe(x,y) == expected); // run on host

        Gpu::PinnedVector<double> result(1);
        auto* pr = result.data();
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) // run one device
        {
            *pr = exe(x,y);
        });
        Gpu::streamSynchronize();
        AMREX_ALWAYS_ASSERT(*pr == expected);
        amrex::Print() << "SUCCESS on host and device\n";
    }

#endif

    amrex::Finalize();
}
