#include <Thread.H>
#include <Profiler.H>
#include <ParallelDescriptor.H>

namespace BoxLib3
{
namespace testing
{
void profiler_main(int& argc, char**& argv);
}
}

extern "C"
void* test_profiler_routine(void*)
{
    BL_PROFILE("tp_routine");
    // BoxLib3::Thread::sleep(1.0);
    return 0;
}

namespace
{
void
thread_timing()
{
    BL_PROFILE("a_thread_timing()");
    FunctionThread ft(test_profiler_routine);
    FunctionThread ft1(test_profiler_routine);
    test_profiler_routine(0);
    ft.join();
    ft1.join();
}
}

int
main(int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);

    BL_PROFILE("BoxLib3::testing::profiler_main()");
    BoxLib::WallTimer wt;
    wt.start();
    wt.stop();
    //  std::cout << "Wall timer reports = " << wt << std::endl;
    thread_timing();

    BoxLib::Finalize();
}
