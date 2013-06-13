//
// A test program for FillBoundary().
//

#include <Utility.H>
#include <MultiFab.H>

int
main (int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);

//    Box bx(IntVect(0,0,0),IntVect(1023,1023,1023));
//    Box bx(IntVect(0,0,0),IntVect(127,127,127));
    Box bx(IntVect(0,0,0),IntVect(255,255,255));

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Domain: " << bx << '\n';

    BoxArray ba(bx);

    ba.maxSize(64);

    const int N = 100;  // This should be divisible by 4 !!!

    if (ParallelDescriptor::IOProcessor())
        std::cout << "# boxes in BoxArray: " << ba.size() << '\n';

    {
        //
        // A test of FillBoundary() on 1 grow cell with cross stencil.
        //
        MultiFab mf(ba,1,1); mf.setVal(1.23);

        double beg = ParallelDescriptor::second();
        for (int i = 0; i < N; i++)
            mf.FillBoundary(false,true);
        double end = (ParallelDescriptor::second() - beg);

        ParallelDescriptor::ReduceRealMax(end,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
            std::cout << "cross x 1: " << end << std::endl;
    }

    {
        //
        // A test of FillBoundary() on 1 grow cell with dense stencil.
        //
        MultiFab mf(ba,1,1); mf.setVal(1.23);

        double beg = ParallelDescriptor::second();
        for (int i = 0; i < N; i++)
            mf.FillBoundary(false,false);
        double end = (ParallelDescriptor::second() - beg);

        ParallelDescriptor::ReduceRealMax(end,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
            std::cout << "dense x 1: " << end << std::endl;
    }

    {
        //
        // First a test of FillBoundary() on 2 grow cells with dense stencil.
        //
        MultiFab mf(ba,1,2); mf.setVal(1.23);

        double beg = ParallelDescriptor::second();
        for (int i = 0; i < N/2; i++)
            mf.FillBoundary(false,false);
        double end = (ParallelDescriptor::second() - beg);

        ParallelDescriptor::ReduceRealMax(end,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
            std::cout << "dense x 2: " << end << std::endl;
    }

    {
        //
        // First a test of FillBoundary() on 4 grow cells with dense stencil.
        //
        MultiFab mf(ba,1,4); mf.setVal(1.23);

        double beg = ParallelDescriptor::second();
        for (int i = 0; i < N/4; i++)
            mf.FillBoundary(false,false);
        double end = (ParallelDescriptor::second() - beg);

        ParallelDescriptor::ReduceRealMax(end,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
            std::cout << "dense x 4: " << end << std::endl;
    }

    BoxLib::Finalize();

    return 0;
}
