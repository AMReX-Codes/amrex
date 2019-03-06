#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BaseFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    long iters = 0;
    int boxsize = 0;
    int ncomps = 0;

    {
        ParmParse pp;
        pp.get("iters", iters);
        pp.get("boxsize",boxsize);
        pp.get("ncomps", ncomps);
    }

        amrex::Print() << std::endl;

#ifdef AMREX_USE_GPU_PRAGMA
        amrex::Print() << "Fortran version of BaseFab testing suite." << std::endl;
#else
        amrex::Print() << "C++ version of BaseFab testing suite." << std::endl;
#endif
        amrex::Print() << "A lot of results should equal 5.5." << std::endl  
                       << "Cubic boxes of length: " << boxsize << std::endl
                       << "Number of components: " << ncomps << std::endl
                       << "Number of iterations of each test: " << iters << std::endl
                       << "=========================================" << std::endl << std::endl;

    // ====================================================================
    // SetVal
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);

        fab1.setVal(1.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.setVal(5.5);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::setVal() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ====================================================================

    // SetValIfNot 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        BaseFab<int> fab2(bx1,ncomps);

        fab1.setVal(1.0);
        fab2.setVal(0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.setValIfNot(5.5, bx1, fab2, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::setValIfNot() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ====================================================================
    // Invert 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);

        fab1.setVal( (iters%2 == 0) ? 5.5 : (1/5.5) );

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.invert(1.0, bx1, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::invert() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ====================================================================
    // Abs 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(-5.5);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.abs(bx1, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::abs() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ====================================================================
    // negate 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal( (iters%2 ==0) ? 5.5 : -5.5 );

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.negate(bx1, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::negate() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ====================================================================
    // Norm 
    {
        double timer;
        double result = 0.0;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5/(bx1.numPts()*ncomps));

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           result = fab1.norm(bx1, 1, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::norm() test." << std::endl
                       << "Result: " << result << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ====================================================================

    // Norminfmask 
    {
        double timer;
        double result = 0.0;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5);

        Box bx2(IntVect(0), IntVect(boxsize-1));
        BaseFab<int> fab2(bx2,ncomps);
        fab2.setVal(1);

        for (int i = 0; i<ncomps; ++i)
        {
           int count = 0;
           IntVect currIndx = bx2.smallEnd();
           while (bx2.contains(currIndx))
           {
              int random = Random_int(2);
              fab2(currIndx, i) = random;
              count += random;
              bx2.next(currIndx);
           }
           amrex::Print() << " *** norm using " << count << " 'true' masks out of " 
                          << bx2.numPts() << " points in component " << i << "." << std::endl;
        }

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           result = fab1.norminfmask(bx1, fab2, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::norminfmask() test." << std::endl
                       << "Result: " << result << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ====================================================================
    // Min
    {
        double timer;
        double result = 0.0;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           result = fab1.min(bx1, 0);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::min() test." << std::endl
                       << "Result: " << result << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ====================================================================
    // Max 
    {
        double timer;
        double result = 0.0;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           result = fab1.max(bx1, 0);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::max() test." << std::endl
                       << "Result: " << result << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ====================================================================
    // minIndex
    {
        double timer;
        IntVect result(-100);

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           result = fab1.minIndex(bx1, 0);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::minIndex() test." << std::endl
                       << "Result: " << result << " = " << fab1(result) << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ====================================================================
    // maxIndex
    {
        double timer;
        IntVect result(-100);

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           result = fab1.maxIndex(bx1, 0);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::maxIndex() test." << std::endl
                       << "Result: " << result << " = " << fab1(result) << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ====================================================================
    // Sum 
    {
        double timer;
        double total = 0;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5/(bx1.numPts()*ncomps));

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           total = fab1.sum(bx1, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::sum() test." << std::endl
                       << "Result: " << total << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ===================================================================
    // plus 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(2.75);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.75/iters);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.plus(fab2, bx2, bx1, 0, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::plus() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ===================================================================
    // minus 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(8.25);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.75/iters);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.minus(fab2, bx2, bx1, 0, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::minus() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // mult
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(2);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.75);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.mult(fab2, bx2, bx1, 0, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::mult() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/ops." << std::endl
                       << "                         or: " << double(iters)/timer <<  " ops/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // mult & divide 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(11.0);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.divide(fab2, bx2, bx1, 0, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::divide() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/ops." << std::endl
                       << "                         or: " << double(iters)/timer <<  " ops/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // protected_divide 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(11.0);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.protected_divide(fab2, bx2, bx1, 0, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::protected_divide() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/ops." << std::endl
                       << "                         or: " << double(iters)/timer   <<  " ops/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // protected_divide w/ zeros
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(0.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.protected_divide(fab2, bx2, bx1, 0, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::protected_divide() test with zeroes." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/ops." << std::endl
                       << "                         or: " << double(iters)/timer   <<  " ops/second." << std::endl << std::endl; 
    }




    // ===================================================================
    // dot 
    {
        double timer;
        double total = 0;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(2.75);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(double(2.0)/(bx2.numPts()*ncomps));

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           total = fab1.dot(bx1, 0, fab2, bx2, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::dot() test." << std::endl
                       << "Result: " << total  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ===================================================================
    // dotmask
    {
        double timer;
        double total = 0;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(2.75/boxsize);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.0);

        Box bx3(IntVect(0), IntVect(boxsize-1));
        BaseFab<int> fab3(bx3,ncomps);
        fab3.setVal(1);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           total = fab1.dotmask(fab3, bx1, 0, fab2, bx2, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::dotmask() test." << std::endl
                       << "Result: " << total << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ===================================================================
    // masks
    {
        double timer;
        int totalLT = 0;
        int totalLE = 0;
        int totalEQ = 0;
        int totalGE = 0;
        int totalGT = 0;

        Box bx1(IntVect(10), IntVect(10+boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5);

        Box bx2(IntVect(10), IntVect(10+boxsize-1));
        BaseFab<int> fab2(bx2,ncomps);
        fab2.setVal(0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           totalLT = fab1.maskLT(fab2, 19, 0);
           totalLE = fab1.maskLE(fab2, 5.5, 0);
           totalEQ = fab1.maskEQ(fab2, 5.5, 0);
           totalGE = fab1.maskGE(fab2, 5.5, 0);
           totalGT = fab1.maskGT(fab2, 1, 0);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::maskLT() test." << std::endl
                       << "Results (LT/LE/EQ/GE/GT): " << totalLT << "/" << totalLE << "/" << totalEQ << "/"
                                                       << totalGE << "/" << totalGT << " of " << fab1.numPts() << " points." << std::endl 
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/(double(iters)/5) << " seconds/ops." << std::endl
                       << "                         or: " << double(iters)/(5*timer) << " ops/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // saxpy 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(1);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.25/iters);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.saxpy(2.0, fab2, bx2, bx1, 0, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::saxpy() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector()) << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // xpay 
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(2.25/(iters));

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(1.0/iters);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.xpay(2/iters, fab2, bx2, bx1, 0, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::xpay() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // addproduct
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(-0.75);

        Box bx2(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.5);

        Box bx3(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab3(bx3,ncomps);
        fab3.setVal(2.5/iters);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.addproduct(bx1, 0, ncomps, fab2, 0, fab3, 0);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::addproduct() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // LinComb
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(1.0);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(2.0);

        Box bx3(IntVect(10), IntVect(10+boxsize-1));
        BaseFab<Real> fab3(bx3,ncomps);
        fab3.setVal(3.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.linComb(fab2, bx2, 0, fab3, bx3, 0, 0.5, 1.5, bx1, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::linComb() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }
    // ===================================================================
    // LinInterp
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(1.0);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(3.0);

        Box bx3(IntVect(10), IntVect(10+boxsize-1));
        BaseFab<Real> fab3(bx3,ncomps);
        fab3.setVal(4.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.linInterp(fab2, bx2, 0, fab3, bx3, 0, 0.5, 1.5, 3.0, bx1, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::linInterp() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // Copy
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(0.0);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(5.5);

        Vector<Real> buffer(bx1.numPts()*ncomps, 0.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.copy(fab2, bx2, 0, bx1, 0, ncomps);
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::copy() test." << std::endl
                       << "Result: " << fab1(IntVect::TheZeroVector())  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/iter." << std::endl
                       << "                         or: " << double(iters)/timer << " iters/second." << std::endl << std::endl; 
    }

    // ===================================================================
    // CopyToMem & CopyFromMem
    {
        double timer;

        Box bx1(IntVect(0), IntVect(boxsize-1));
        BaseFab<Real> fab1(bx1,ncomps);
        fab1.setVal(5.5);

        Box bx2(IntVect(1000), IntVect(1000+boxsize-1));
        BaseFab<Real> fab2(bx2,ncomps);
        fab2.setVal(0.0);

        Vector<Real> buffer(bx1.numPts()*ncomps, 0.0);

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab1.copyToMem(bx1, 0, ncomps, buffer.data());
        }
        timer = second() - timer;

        timer = second();
        for (int i=0; i<iters; ++i)
        {
           fab2.copyFromMem(bx2, 0, ncomps, buffer.data());
        }
        timer = second() - timer;

        amrex::Print() << "BaseFab<Real>::copyToMem() test." << std::endl
                       << "Result: " << fab2(IntVect(1000))  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/op." << std::endl
                       << "                         or: " << double(iters)/timer << " ops/second." << std::endl << std::endl; 

        amrex::Print() << "BaseFab<Real>::copyFromMem() test." << std::endl
                       << "Result: " << fab2(IntVect(1000))  << std::endl
                       << "Completed in: "                <<  timer << " seconds." << std::endl
                       << " or, completed at a rate of: " <<         timer/iters << " seconds/op." << std::endl
                       << "                         or: " << double(iters)/timer << " ops/second." << std::endl << std::endl; 
    }
    // ===================================================================

    amrex::Finalize();
}
