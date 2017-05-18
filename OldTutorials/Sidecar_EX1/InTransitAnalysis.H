#ifndef INTRANSITANALYSIS_H
#define INTRANSITANALYSIS_H

class InTransitAnalysis
{
    public:
        void Initialize (amrex::MultiFab&, amrex::Geometry&, int);
        void DoAnalysis ();
        void Finalize ();
        void PrintResults () const;

    private:
        amrex::MultiFab* mf;
        amrex::Geometry* geom;
        Real norm0;
        Real probsize;
        int time_step;
};

#endif /* INTRANSITANALYSIS_H */
