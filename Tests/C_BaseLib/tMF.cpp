//
// A test program for MultiFab.
//

#include <map>
#include <vector>

#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);

    typedef std::map<Vector<int>,Vector<Real> > OurBinMap;

    OurBinMap bins;

    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    //
    // The number of "ints" in the key to our map.
    //
    const int nKeys = BL_SPACEDIM;
    const int nVals = 1;
    //
    // Each CPU will put a different # of things in map.
    //
    for (int i = 0; i < 10 + (10*MyProc); i++)
    {
      Vector<int> iv(nKeys);
      for (int j=0; j<nKeys; ++j) {
        iv[j] = 1000*MyProc+i;
      }
      bins[iv].resize(1, Real(MyProc) );
    }

#ifdef BL_USE_MPI
    //
    // The number of keys and their respective offsets.
    // Used by IOProc but the "1" is so others can call dataPtr()
    // without error.
    //
    Vector<int> nmkeys(1);
    Vector<int> nmvals(1);
    Vector<int> nmentries(1);
    Vector<int> offset(1);

    if (ParallelDescriptor::IOProcessor())
    {
         nmentries.resize(NProcs,0);
         nmkeys.resize(NProcs,0);
         nmvals.resize(NProcs,0);
         offset.resize(NProcs,0);
    }
    //
    // Here lcount is the number of entries.
    //
    int lcount = bins.size();

    MPI_Gather(&lcount,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nmentries.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());

    //
    // Each CPU must pack its data into simple arrays.
    //
    Vector<int>  lkeys;
    Vector<Real> lvals;

    for (OurBinMap::const_iterator it = bins.begin(), End = bins.end();
         it != End;
         ++it)
    {
        const Vector<int>& iv = it->first;
        BL_ASSERT(iv.size() == nKeys);

        for (int i = 0; i < nKeys; i++)
            lkeys.push_back(iv[i]);

        const Vector<Real>& vals = it->second;
        BL_ASSERT(vals.size() == nVals);

        for (int i = 0; i < nVals; i++)
            lvals.push_back(vals[i]);
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
    Vector<int>  keys(1);
    Vector<Real> vals(1);
    //
    // First pass the vals.
    //
    if (ParallelDescriptor::IOProcessor())
    {
        int tcount = 0;
        for (int i = 0; i < NProcs; i++)
        {
            nmvals[i] = nmentries[i] * nVals;
            tcount    += nmvals[i];
        }

        vals.resize(tcount);

        for (int i = 1, N = offset.size(); i < N; i++) {
            offset[i] = offset[i-1] + nmvals[i-1];
        }
    }

    ParallelDescriptor::Barrier();

    int vcount = lcount * nVals;

    MPI_Gatherv(lcount == 0 ? 0 : lvals.dataPtr(),
                vcount,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                vals.dataPtr(),
                nmvals.dataPtr(),
                offset.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    //
    // Then the keys
    //
    if (ParallelDescriptor::IOProcessor())
    {

        int tcount = 0;
        for (int i = 0; i < NProcs; i++)
        {
            nmkeys[i] = nmentries[i] * nKeys;
            tcount    += nmkeys[i];
        }
        keys.resize(tcount);

        for (int i = 1, N = offset.size(); i < N; i++) {
            offset[i] = offset[i-1] + nmkeys[i-1];
        }
    }

    //
    // There are nKeys ints in each Key.
    //
    int kcount = lcount * nKeys;

    MPI_Gatherv(kcount == 0 ? 0 : lkeys.dataPtr(),
                kcount,
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
        int nentries = vals.size() / nVals;
        BL_ASSERT(nentries * nVals == vals.size());

        for (int i = 0; i < nentries; i++)
        {
            int ik = i * nKeys;
            Vector<int> key(nKeys);
            for (int k=0; k<nKeys; ++k) {
                key[k] = keys[ik+k];
            }

            int iv = i * nVals;
            Vector<Real> val(nVals);
            for (int k=0; k<nVals; ++k) {
                val[k] = vals[iv+k];
            }

            bins[key] = val;
        }

        std::cout << "Got " << bins.size() << " (key,val) pairs:\n";
    }
#endif

    amrex::Finalize();

    return 0;
}
