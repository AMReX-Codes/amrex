#ifndef INSITUANALYSIS_H
#define INSITUANALYSIS_H

#include <Analysis.H>
#include <Geometry.H>
#include <MultiFab.H>

class InSituAnalysis: public Analysis
{
    public:
        InSituAnalysis(MultiFab&, Geometry&);
        void DoAnalysis ();
        void PrintResults () const;

    private:
        MultiFab* mf;
        Geometry* geom;
        Real norm0;
        Real probsize;
};

#endif /* INSITUANALYSIS_H */
