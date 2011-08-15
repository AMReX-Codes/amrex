#ifndef _SlabStat_H_
#define _SlabStat_H_

#include <list>
#include <ArrayLim.H>
#include <REAL.H>
#include <Box.H>
#include <MultiFab.H>
//
// Forward declaration.
//
class AmrLevel;

extern "C"
{
    //
    // Type of extern "C" function called by SlabStat to compute a statistic.
    //
    // dst contains the statistics to be updated in place.
    // src contains the variables need to calculate the statistics.
    //
    typedef void (*SlabStatFunc)(const Real* src,
                                 ARLIM_P(src_lo), ARLIM_P(src_hi),
                                 const int* nsrc,
                                 Real* dst, ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                 const int* ndst,
                                 const Real* dt,
                                 const Real* dx);
}

//
// Slab Statistics Record
//
//
//  SlabStatRecs are designed to control the output of running averages as 
//  a simulation proceeds.  A required list of state and derived variables
//  is specified, as well as a list of boxes on which the statistics are
//  to be calculated.  Every time step on the level where these boxes are 
//  defined, the specified fortran function is called and the running statistics
//  are accumulated.  The statistics are written out every check_int level 0
//  timesteps.
//
class SlabStatRec
{
   friend class SlabStatList;

public:
    //
    // The name of the statistic.
    //
    const std::string& name () const;
    //
    // Number of components in the statistic.
    //
    int nComp () const;
    //
    // The names of State/Derived variables needed for statistic.
    //
    const Array<std::string>& vars () const;
    //
    // Number of variables used to calculate statistic.
    //
    long nVariables () const;
    //
    // Ghost cells needed in the variables when calculating statistic.
    //
    int nGrow () const;
    //
    // Level on which to calculate the statistic
    //
    int level () const;
    //
    // Boxes, at specified level, on which to calculate the statistic.
    //
    const BoxArray& boxes () const;
    //
    // The SlabStatFunc used to calculate the statistic.
    //
    SlabStatFunc func () const;
    //
    // The MultiFab in which statistics are accumulated.
    //
    MultiFab& mf ();
    //
    // Temporary (grown) MultiFab used in calculating the statistic.
    //
    MultiFab& tmp_mf ();
    //
    // Interval over which the stats have been calculated.
    //
    Real interval () const;

protected:
    //
    // A Constructor.
    //
    SlabStatRec (const std::string&  name,
                 int             ncomp,
                 Array<std::string>& vars,
                 int             ngrow,
                 int             level,
                 const BoxArray& boxes,
                 SlabStatFunc    func);
    //
    // Another Constructor.
    //
    SlabStatRec (const std::string&  name,
                 int             ncomp,
                 Array<std::string>& vars,
                 int             ngrow,
                 SlabStatFunc    func);
    //
    // Destructor.
    //
    ~SlabStatRec ();
    //
    // Name of statistic
    //
    std::string m_name;
    //
    // Number of components in the statistic.
    //
    int m_ncomp;
    //    
    // Names of State/Derived variables needed to calculate statistic.
    //
    Array<std::string> m_vars;
    //
    // Ghost cells needed in variables when calculating the statistic.
    //
    int m_ngrow;
    //
    // Level on which to calculate statistic.
    //
    int m_level;
    //
    // Boxes, at specified level, on which to calculate the statistic.
    //
    BoxArray m_boxes;
    //
    // Function that computes the statistic.
    //
    SlabStatFunc m_func;
    //
    // MultiFab in which to accumulate statistics across checkpoints.
    //
    MultiFab m_mf;
    //
    // Temporary (grown) MultiFab  used in calculating the statistic.
    //
    MultiFab m_tmp_mf;
    //
    // Time interval over which the Stats have been summed.
    //
    Real m_interval;

private:
    //
    // Disallowed.
    //
    SlabStatRec (const SlabStatRec&);
    SlabStatRec& operator= (const SlabStatRec&);
};

//
// A list of SlabStatRecs.
//
// SlabStatList manages and provides access to the list of SlabStatRecs.
//

class SlabStatList
{
public:
    //
    // The default constructor.
    //
    SlabStatList ();
    //
    // The destructor.
    //
    ~SlabStatList ();
    //
    // Adds another entry to the list of SlabStats.
    //
    void add (const std::string&  name,
              int                 ncomp,
              Array<std::string>& vars,
              int                 ngrow,
              int                 level,
              const BoxArray&     boxes,
              SlabStatFunc        func);
    //
    // Another way to add an entry to the list of SlabStats.
    //
    void add (const std::string&  name,
              int                 ncomp,
              Array<std::string>& vars,
              int                 ngrow,
              SlabStatFunc        func);
    //
    // Update the statistics in the list.
    //
    void update (AmrLevel& amrlevel, Real time, Real dt);
    //
    // Write out any SlabStats in the list.
    //
    void checkPoint (PArray<AmrLevel>& amrLevels, int level0_step);

    std::list<SlabStatRec*>& list ();

private:
    //
    // Disallowed.
    //
    SlabStatList (const SlabStatList&);
    SlabStatList& operator= (const SlabStatList&);

    std::list<SlabStatRec*> m_list;
};

#endif /*_SlabStat_H_*/
