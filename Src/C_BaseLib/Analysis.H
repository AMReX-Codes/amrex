#ifndef BL_ANALYSIS_H
#define BL_ANALYSIS_H

#include <MultiFab.H>
#include <Geometry.H>

/* An abstract base class for doing analysis on BoxLib data sets on the fly.
 * This can be used either in-situ or in-transit (on the "sidecar" MPI
 * processes). A derived class of this type should define both the
 * functionality to perform the analysis itself, as well as the storage for
 * holding the results, in whatever form it may take (a MultiFab, an array of
 * ints, etc.).
 *
 * Different types of analyses will require different types of input data. We
 * can't have virtual constructors, so the way we handle instantiation of any
 * type of analysis class is by making Initialize() virtual but *not* pure.
 * Furthermore, we must include here every type of analysis class
 * initialization we know of.
 *
 * */

class AnalysisCallbackInterface
{
    public:
        virtual void Initialize () { BoxLib::Error("I AM A VIRTUAL FUNCTION. WHY ARE YOU CALLING ME? :-("); };
        virtual void Initialize (MultiFab&, Geometry&, int) { BoxLib::Error("I AM A VIRTUAL FUNCTION. WHY ARE YOU CALLING ME? :-("); };
        virtual void DoAnalysis () = 0;
        virtual void Finalize () = 0;
};

/* A generic class for emulating callbacks to user-defined analysis functions. */

class AnalysisContainer
{
    public:
        void connectCallback(AnalysisCallbackInterface *aci) { m_aci = aci; }

        void DoAnalysis () { m_aci->DoAnalysis (); }
        void Finalize   () { m_aci->Finalize   (); }

    private:
        AnalysisCallbackInterface *m_aci;
};

namespace Analysis
{
    // Various signals to send to the sidecar group.
    const extern int NyxHaloFinderSignal;
    const extern int QuitSignal;

    extern AnalysisContainer *analysis;
}

#endif /* BL_ANALYSIS_H */
