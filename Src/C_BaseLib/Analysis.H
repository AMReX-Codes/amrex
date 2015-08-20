#ifndef BL_ANALYSIS_H
#define BL_ANALYSIS_H

/* An abstract base class for doing analysis on BoxLib data sets on the fly.
 * This can be used either in-situ or in-transit (on the "sidecar" MPI
 * processes). A derived class of this type should define both the
 * functionality to perform the analysis itself, as well as the storage for
 * holding the results, in whatever form it may take (a MultiFab, an array of
 * ints, etc.). */

#include <MultiFab.H>

class Analysis
{
    public:
        virtual void DoAnalysis () = 0;
};

#endif /* BL_ANALYSIS_H */
