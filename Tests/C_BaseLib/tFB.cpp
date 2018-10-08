//
// A test program for FillBoundary().
//

#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

const int nTimes(5);
const int nStrategies(4);


int
main (int argc, char** argv)
{
  amrex::Initialize(argc, argv);

  BL_PROFILE_VAR("main()", pmain);

  Vector<DistributionMapping::Strategy> dmStrategies(nStrategies);
  dmStrategies[0] = DistributionMapping::ROUNDROBIN;
  dmStrategies[1] = DistributionMapping::KNAPSACK;
  dmStrategies[2] = DistributionMapping::SFC;
  dmStrategies[3] = DistributionMapping::PFC;

  Vector<std::string> dmSNames(nStrategies);
  dmSNames[0] = "ROUNDROBIN";
  dmSNames[1] = "KNAPSACK";
  dmSNames[2] = "SFC";
  dmSNames[3] = "PFC";

  Vector<double> dmSTimes(nStrategies, 0.0);

  for(int iS(0); iS < nStrategies * nTimes; ++iS) {

    int whichStrategy(iS % nStrategies);

    DistributionMapping::strategy(dmStrategies[whichStrategy]);

//    Box bx(IntVect(0,0,0),IntVect(511,511,255));
//    Box bx(IntVect(0,0,0),IntVect(1023,1023,255));
    Box bx(IntVect(0,0,0),IntVect(1023,1023,1023));
//    Box bx(IntVect(0,0,0),IntVect(2047,2047,1023));
//    Box bx(IntVect(0,0,0),IntVect(127,127,127));
//    Box bx(IntVect(0,0,0),IntVect(255,255,255));

    BoxArray ba(bx);
    ba.maxSize(64);

    DistributionMapping dm{ba};

    const int N = 2000;  // This should be divisible by 4 !!!

    if (ParallelDescriptor::IOProcessor() && iS == 0) {
        std::cout << "Domain: " << bx << "  # boxes in BoxArray:  " << ba.size() << '\n';
    }
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Strategy: " << dmSNames[DistributionMapping::strategy()] << '\n';


    ParallelDescriptor::Barrier();

    {
        //
        // A test of FillBoundary() on 1 grow cell with cross stencil.
        //
        MultiFab mf(ba,dm,1,1); mf.setVal(1.23);

        ParallelDescriptor::Barrier();
        double beg = ParallelDescriptor::second();
        for (int i = 0; i < N; i++)
            mf.FillBoundary(true);
        double end = (ParallelDescriptor::second() - beg);

        ParallelDescriptor::ReduceRealMax(end,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor()) {
          std::cout << N << " cross x 1: " << end << std::endl;
	  dmSTimes[whichStrategy] += end;
	}
    }


    {
        //
        // A test of FillBoundary() on 1 grow cell with dense stencil.
        //
        MultiFab mf(ba,dm,1,1); mf.setVal(1.23);

        ParallelDescriptor::Barrier();
        double beg = ParallelDescriptor::second();
        for (int i = 0; i < N; i++)
            mf.FillBoundary();
        double end = (ParallelDescriptor::second() - beg);

        ParallelDescriptor::ReduceRealMax(end,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor()) {
          std::cout << N << " dense x 1: " << end << std::endl;
	  dmSTimes[whichStrategy] += end;
	}
    }

    {
        //
        // First a test of FillBoundary() on 2 grow cells with dense stencil.
        //
        MultiFab mf(ba,dm,1,2); mf.setVal(1.23);

        ParallelDescriptor::Barrier();
        double beg = ParallelDescriptor::second();
        for (int i = 0; i < N/2; i++)
            mf.FillBoundary();
        double end = (ParallelDescriptor::second() - beg);

        ParallelDescriptor::ReduceRealMax(end,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor()) {
          std::cout << (N/2) << " dense x 2: " << end << std::endl;
	  dmSTimes[whichStrategy] += end;
	}
    }

    {
        //
        // First a test of FillBoundary() on 4 grow cells with dense stencil.
        //
        MultiFab mf(ba,dm,1,4); mf.setVal(1.23);

        ParallelDescriptor::Barrier();
        double beg = ParallelDescriptor::second();
        for (int i = 0; i < N/4; i++)
            mf.FillBoundary();
        double end = (ParallelDescriptor::second() - beg);

        ParallelDescriptor::ReduceRealMax(end,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor()) {
          std::cout << (N/4) << " dense x 4: " << end << std::endl;
	  dmSTimes[whichStrategy] += end;
	}
    }
    if (ParallelDescriptor::IOProcessor())
        std::cout << std::endl;

  }  // end for iS


    if(ParallelDescriptor::IOProcessor()) {
      for(int i(0); i < nStrategies; ++i) {
        std::cout << std::endl << "Total times:" << std::endl;
	std::cout << dmSNames[i] << " time = " << dmSTimes[i] << std::endl;
      }
      std::cout << std::endl << std::endl;
    }

    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();

    return 0;
}
