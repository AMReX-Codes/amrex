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

    typedef std::map<IntVect,Real,IntVect::Compare> OurBinMap;

    OurBinMap bins;

    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //
    // Each CPU will put a different # of things in map.
    //
    for (int i = 0; i < 10 + (10*MyProc); i++)
    {
        IntVect iv(D_DECL(1000*MyProc+i, 1000*MyProc+i, 1000*MyProc+i));

        bins[iv] = Real(MyProc);
    }
    //
    // The number of "ints" in the key to our map.
    //
    const int KeySize = BL_SPACEDIM;

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
    // Here lcount is the number of keys.
    //
    int lcount = bins.size();

    MPI_Gather(&lcount,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nmkeys.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());
    //
    // Each CPU must pack its data into simple arrays.
    //
    Array<int>  lkeys;
    Array<Real> lvals;

    for (OurBinMap::const_iterator it = bins.begin(), End = bins.end();
         it != End;
         ++it)
    {
        //
        // Our "KeySize" ints per key.
        //
        const IntVect& iv = it->first;

        for (int i = 0; i < BL_SPACEDIM; i++)
            lkeys.push_back(iv[i]);
        //
        // And our value.
        //
        lvals.push_back(it->second);
    }
    //
    // Clear all the "bins" including IOProc's.
    // For IOProc we do this so that when we add back
    // things into "bin" at the end we don't double
    // count anything.  For all others it's just to cut
    // out unnecessary memory since MPI will need to use
    // some memory for the Gatherv().
    //
    bins.clear();
    //
    // The "global" (keys,vals) pairs into which all processors send their data.
    //
    Array<int>  keys(1);
    Array<Real> vals(1);
    //
    // First pass the vals.
    //
    if (ParallelDescriptor::IOProcessor())
    {
        int tcount = 0;
        for (int i = 0; i < NProcs; i++)
            tcount += nmkeys[i];

        vals.resize(tcount);

        for (int i = 1, N = offset.size(); i < N; i++)
            offset[i] = offset[i-1] + nmkeys[i-1];
    }

    MPI_Gatherv(lcount == 0 ? 0 : lvals.dataPtr(),
                lcount,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                vals.dataPtr(),
                nmkeys.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());
    //
    // Then the values.
    // Don't forget to update offset and nmkeys appropriately.
    //
    if (ParallelDescriptor::IOProcessor())
    {
        int tcount = 0;
        for (int i = 0; i < NProcs; i++)
        {
            nmkeys[i] *= KeySize;
            tcount    += nmkeys[i];
        }
        keys.resize(tcount);

        for (int i = 1, N = offset.size(); i < N; i++)
            offset[i] = offset[i-1] + nmkeys[i-1];
    }
    //
    // There are KeySize ints in each Key.
    //
    lcount *= KeySize;

    MPI_Gatherv(lcount == 0 ? 0 : lkeys.dataPtr(),
                lcount,
                ParallelDescriptor::Mpi_typemap<int>::type(),
                keys.dataPtr(),
                nmkeys.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<int>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Put the (keys,vals) into our map and then print out map.
        //
        BL_ASSERT(keys.size() % KeySize == 0);
        BL_ASSERT(vals.size() * KeySize == keys.size());

        for (int i = 0; i < vals.size(); i++)
        {
            int ik = i * KeySize;

            IntVect iv(D_DECL(keys[ik],keys[ik+1],keys[ik+2]));

            bins[iv] += vals[i];
        }

        std::cout << "Got " << bins.size() << " (key,val) pairs:\n";

        for (OurBinMap::const_iterator it = bins.begin(), End = bins.end();
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
