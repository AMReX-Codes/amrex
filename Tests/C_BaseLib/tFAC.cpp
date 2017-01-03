
//
// A simple program to test FabArray<T>::copy() in parallel.
//

#if !(BL_SPACEDIM==2)
#error "This code assumes BL_SPACEDIM==2"
#endif

#include <AMReX_MultiFab.H>

using namespace amrex;

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);

    BL_ASSERT(ParallelDescriptor::NProcs() == 2);

    BoxArray ba_1(1);
    BoxArray ba_2(5);

    ba_1.set(0, Box(IntVect(0,0), IntVect(4,4)));

    ba_2.set(0, Box(IntVect( 1,0), IntVect(3,0)));
    ba_2.set(1, Box(IntVect( 0,1), IntVect(4,1)));
    ba_2.set(2, Box(IntVect(-1,2), IntVect(3,2)));
    ba_2.set(3, Box(IntVect( 1,3), IntVect(5,3)));
    ba_2.set(4, Box(IntVect(-1,4), IntVect(5,4)));

    DistributionMapping dm_1 {ba_1};
    DistributionMapping dm_2 {ba_2};

    MultiFab mf_1(ba_1,dm_1,1,0);

    MultiFab mf_2(ba_2,dm_2,2,0);
    //
    // Set all on mf_1 to zero.
    //
    mf_1.setVal(0);
    //
    // Set the first component of mf_2 to zero.
    //
    mf_2.setVal(0,0,1,0);
    //
    // Set second component to relevent index.
    //
    for (int i = 0; i < mf_2.size(); i++)
    {
        mf_2.setVal(i+1,ba_2[i],1,1,0);
    }

    mf_1.copy(mf_2, 1, 0, 1);

    if (ParallelDescriptor::IOProcessor())
    {
	std::cout << mf_1[0] << std::endl;
    }

    amrex::Finalize();
}
