#ifndef _PERIODICITY_H_
#define _PERIODICITY_H_

class Periodicity
{
public:
    Periodicity () : period(D_DECL(0,0,0)) {}
    explicit Periodicity (const IntVect& v) : period(v) {}

    bool isAnyPeriodic () 
	{ return D_TERM(period[0]>0, || period[1]>0, || period[2]>0); }
    bool isAllPeriodic () 
	{ return D_TERM(period[0]>0, && period[1]>0, && period[2]>0); }
    bool isPeriodic (int dir)
	{ return period[dir]>0; }
private:
    IntVect period;
};

#endif
