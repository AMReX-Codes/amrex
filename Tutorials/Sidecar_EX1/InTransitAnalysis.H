#ifndef INTRANSITANALYSIS_H
#define INTRANSITANALYSIS_H

#include <Analysis.H>
#include <Geometry.H>
#include <MultiFab.H>

class InTransitAnalysis: public AnalysisCallbackInterface
{
    public:
        InTransitAnalysis (MultiFab&, Geometry&, int);
        void Initialize ();
        void DoAnalysis ();
        void Finalize ();
        void PrintResults () const;

    private:
        MultiFab* mf;
        Geometry* geom;
        Real norm0;
        Real probsize;
        int time_step;
};

#endif /* INTRANSITANALYSIS_H */
