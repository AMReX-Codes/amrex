#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Utility.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    long iters = 0;
    int boxsize = 0;

    {
        ParmParse pp;
        pp.get("iters", iters);
        pp.get("boxsize",boxsize);
    }

#ifdef AMREX_USE_GPU_PRAGMA
        amrex::Print() << "Fortran version of BaseFab testing suite." << std::endl;
#else
        amrex::Print() << "C++ version of BaseFab testing suite." << std::endl;
#endif
        amrex::Print() << "Cubic boxes of length: " << boxsize << std::endl;
        amrex::Print() << "Number of iterations of each test: " << iters << std::endl;
        amrex::Print() << "=========================================" << iters << std::endl << std::endl;

    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,1);
        fab1.setVal(1.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.setVal(5.5);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::setVal() test." << std::endl;
        amrex::Print() << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl;
        amrex::Print() << "Completed in: "                <<  timer << " seconds." << std::endl;
        amrex::Print() << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl;
        amrex::Print() << "                         or: " << double(iters)/timer << " iters/second." << std::endl; 
    }

    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,1);
        fab1.setVal(1.0);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,1);
        fab2.setVal(2.0);

        Box bx3(IntVect(10), IntVect(10+boxsize-1));
        BaseFab<Real> fab3(bx3,1);
        fab3.setVal(3.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.linComb(fab2, bx2, 0, fab3, bx3, 0, 0.5, 1.5, bx1, 0, 1);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::linComb() test." << std::endl;
        amrex::Print() << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl;
        amrex::Print() << "Completed in: "                <<  timer << " seconds." << std::endl;
        amrex::Print() << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl;
        amrex::Print() << "                         or: " << double(iters)/timer << " iters/second." << std::endl; 
    }

    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,1);
        fab1.setVal(5.5);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,1);
        fab2.setVal(0.0);

        Vector<Real> buffer(bx1.numPts(), 0.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.copyToMem(bx1, 0, 1, buffer.data());
           fab2.copyFromMem(bx2, 0, 1, buffer.data());
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::copyToMem() and copyFromMem() test." << std::endl;
        amrex::Print() << "Result: " << fab2(IntVect::TheZeroVector())  << std::endl;
        amrex::Print() << "Completed in: "                <<  timer << " seconds." << std::endl;
        amrex::Print() << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl;
        amrex::Print() << "                         or: " << double(iters)/timer << " iters/second." << std::endl; 
    }

    amrex::Finalize();
}
