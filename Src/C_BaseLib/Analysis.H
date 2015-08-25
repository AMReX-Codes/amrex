#ifndef BL_ANALYSIS_H
#define BL_ANALYSIS_H

/* An abstract base class for doing analysis on BoxLib data sets on the fly.
 * This can be used either in-situ or in-transit (on the "sidecar" MPI
 * processes). A derived class of this type should define both the
 * functionality to perform the analysis itself, as well as the storage for
 * holding the results, in whatever form it may take (a MultiFab, an array of
 * ints, etc.). */

class AnalysisCallbackInterface
{
    public:
        virtual void Initialize () = 0;
        virtual void DoAnalysis () = 0;
        virtual void Finalize () = 0;
};

class Analysis
{
    public:
        void connectCallback(AnalysisCallbackInterface *aci) { m_aci = aci; }

        void Initialize () { m_aci->Initialize (); }
        void DoAnalysis () { m_aci->DoAnalysis (); }
        void Finalize   () { m_aci->Finalize   (); }

    private:
        AnalysisCallbackInterface *m_aci;
};

#endif /* BL_ANALYSIS_H */
