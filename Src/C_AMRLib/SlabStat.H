//BL_COPYRIGHT_NOTICE

#ifndef _SlabStat_H_
#define _SlabStat_H_

//
// $Id: SlabStat.H,v 1.4 2000-04-21 22:28:16 sstanley Exp $
//

#include <ArrayLim.H>
#include <REAL.H>
#include <List.H>
#include <aString.H>
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
// `dst' contains the statistics to be updated in place.
// `src' contains the variables need to calculate the statistics.
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
//@Man:
//@Memo: Slab Statistics Record
/*@Doc:

  SlabStatRecs are designed to control the output of running averages as 
  a simulation proceeds.  A required list of state and derived variables
  is specified, as well as a list of boxes on which the statistics are
  to be calculated.  Every time step on the level where these boxes are 
  defined, the specified fortran function is called and the running statistics
  are accumulated.  The statistics are written out every check_int level 0
  timesteps.
*/
class SlabStatRec
{
   friend class SlabStatList;

public:
    //
    //@ManDoc: The name of the statistic.
    //
    const aString& name () const { return m_name; }
    //
    //@ManDoc: Number of components in the statistic.
    //
    int nComp () const { return m_ncomp; }
    //
    //@ManDoc: The names of State/Derived variables needed for statistic.
    //
    const Array<aString>& vars () const { return m_vars; }
    //
    //@ManDoc: Number of variables used to calculate statistic.
    //
    long nVariables () const { return m_vars.length(); }
    //
    //@ManDoc: Ghost cells needed in the variables when calculating statistic.
    //
    int nGrow () const { return m_ngrow; }
    //
    //@ManDoc: Level on which to calculate the statistic
    //
    int level () const { return m_level; }
    //
    //@ManDoc: Boxes, at specified level, on which to calculate the statistic.
    //
    const BoxArray& boxes () const { return m_boxes; }
    //
    //@ManDoc: The SlabStatFunc used to calculate the statistic.
    //
    SlabStatFunc func () const { return m_func; }
    //
    //@ManDoc: The MultiFab in which statistics are accumulated.
    //
    MultiFab& mf () { return m_mf; }
    //
    //@ManDoc: Temporary (grown) MultiFab used in calculating the statistic.
    //
    MultiFab& tmp_mf () { return m_tmp_mf; }
    //
    //@ManDoc: Interval over which the stats have been calculated.
    //
    Real interval () const { return m_interval; }

protected:
    //
    // A Constructor.
    //
    SlabStatRec (const aString&  name,
                 int             ncomp,
                 Array<aString>& vars,
                 int             ngrow,
                 int             level,
                 const BoxArray& boxes,
                 SlabStatFunc    func);
    //
    // Another Constructor.
    //
    SlabStatRec (const aString&  name,
                 int             ncomp,
                 Array<aString>& vars,
                 int             ngrow,
                 int             level,
                 SlabStatFunc    func);
    //
    // Destructor.
    //
    ~SlabStatRec ();
    //
    // Name of statistic
    //
    aString m_name;
    //
    // Number of components in the statistic.
    //
    int m_ncomp;
    //    
    // Names of State/Derived variables needed to calculate statistic.
    //
    Array<aString> m_vars;
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
//@Man:
//@Memo: A list of SlabStatRecs.
/*@Doc:

  SlabStatList manages and provides access to the list of SlabStatRecs.
*/

class SlabStatList
{
public:
    //
    //@ManDoc: The default constructor.
    //
    SlabStatList () {}
    //
    //@ManDoc: The destructor.
    //
    ~SlabStatList ();
    //
    //@ManDoc: Adds another entry to the list of SlabStats.
    //
    void add (const aString&  name,
              int             ncomp,
              Array<aString>& vars,
              int             ngrow,
              int             level,
              const BoxArray& boxes,
              SlabStatFunc    func)
    {
        m_list.add(new SlabStatRec(name,ncomp,vars,ngrow,level,boxes,func));
    }
    //
    //@ManDoc: Another way to add an entry to the list of SlabStats.
    //
    void add (const aString&  name,
              int             ncomp,
              Array<aString>& vars,
              int             ngrow,
              int             level,
              SlabStatFunc    func)
    {
        m_list.add(new SlabStatRec(name,ncomp,vars,ngrow,level,func));
    }
    //
    //@ManDoc: Update the statistics in the list.
    //
    void update (AmrLevel& amrlevel, Real time, Real dt);
    //
    //@ManDoc: Write out any SlabStats in the list.
    //
    void checkPoint (const aString& ckfile);

    List<SlabStatRec*>& list () { return m_list; }

protected:

    List<SlabStatRec*> m_list;

private:
    //
    // Disallowed.
    //
    SlabStatList (const SlabStatList&);
    SlabStatList& operator= (const SlabStatList&);
};


#endif /*_SlabStat_H_*/
