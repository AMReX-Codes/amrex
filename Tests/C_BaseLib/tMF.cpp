//
// A test program for MultiFab.
//

#include <winstd.H>
#include <map>
#include <vector>

#include <Utility.H>
#include <MultiFab.H>
 
int
main (int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);

    std::map<int,Real> bins;

    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //
    // Each CPU will put a different # of things in map.
    //
    for (int i = 0; i < 10 + (10*MyProc); i++)
    {
        bins[1000*MyProc+i] = Real(MyProc);
    }

#if BL_USE_MPI
    //
    // The number of keys and their respective offsets.
    // Used by IOProc but the "1" is so others can call dataPtr()
    // without error.
    //
    Array<int> nmkeys(1);
    Array<int> offset(1);

    if (ParallelDescriptor::IOProcessor())
    {
         nmkeys.resize(NProcs,0);
         offset.resize(NProcs,0);
    }
    //
    // Tell root CPU how many tags each CPU will be sending.
    // The I/O processor doesn't send anything so doesn't add to this count.
    //
    int lcount = bins.size();
    //
    // Tell root CPU how many tags each CPU will be sending.
    // Fills in "nmkeys" array on IOProc.
    //
    MPI_Gather(&lcount,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nmkeys.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());
    //
    // The "global" (keys,vals) pairs into which all processors send their data.
    //
    Array<int>  keys(1);
    Array<Real> vals(1);

    if (ParallelDescriptor::IOProcessor())
    {
        int tcount = 0;

        for (int i = 0; i < NProcs; i++)
            tcount += nmkeys[i];
        //
        // Where IOProc stores sent data.
        //
        keys.resize(tcount);
        vals.resize(tcount);

        for (int i = 1, N = offset.size(); i < N; i++)
            offset[i] = offset[i-1] + nmkeys[i-1];
    }
    //
    // Each CPU must pack its data into simple arrays.
    //
    Array<int>  lkeys(bins.size() + 1);
    Array<Real> lvals(bins.size() + 1);

    int idx = 0;

    for (std::map<int,Real>::const_iterator it = bins.begin(), End = bins.end();
         it != End;
         ++it, ++idx)
    {
        lkeys[idx] = it->first;
        lvals[idx] = it->second;
    }
    //
    // Clear all the "bins" including IOProc's.
    // For IOProc we do this so that when we add back
    // things into "bin" at the end we don't double
    // count anything.  For all others it's just the cut
    // out unnecessary memory since MPI will need to use
    // some memory for the Gatherv().
    //
    bins.clear();
    //
    // First pass the keys.
    //
    MPI_Gatherv(lkeys.dataPtr(),
                lcount,
                ParallelDescriptor::Mpi_typemap<int>::type(),
                keys.dataPtr(),
                nmkeys.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<int>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    ParallelDescriptor::Barrier();
    //
    // Then the values.
    //
    MPI_Gatherv(lvals.dataPtr(),
                lcount,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                vals.dataPtr(),
                nmkeys.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Put the (keys,vals) into our map and then print out map.
        //
        for (int i = 0; i < keys.size(); i++)
        {
            bins[keys[i]] += vals[i];
        }

        std::cout << "Got " << bins.size() << " (key,val) pairs:\n";

        for (std::map<int,Real>::const_iterator it = bins.begin(), End = bins.end();
             it != End;
             ++it)
        {
            std::cout << it->first << ' ' << it->second << '\n';
        }
    }
#endif

    BoxLib::Finalize();

    return 0;
}
