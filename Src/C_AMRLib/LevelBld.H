
#ifndef _LEVELBLD_H_
#define _LEVELBLD_H_

//
// $Id: LevelBld.H,v 1.6 2000-10-02 20:48:42 lijewski Exp $
//

#include <REAL.H>
#include <Box.H>
#include <Geometry.H>
#include <Amr.H>
#include <AmrLevel.H>

extern LevelBld* getLevelBld ();

//
//@Man:
//@Memo: Builds problem-specific AmrLevels.
/*@Doc:

  Abstract base class specifying an interface for building problem-specific
  AmrLevels.  
*/

class LevelBld
{
public:

    /*@ManDoc: Perform any problem-dependent setup such as physical
               boundary condition and derived quantities.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual void variableSetUp () = 0;

    /*@ManDoc: Perform any problem-dependent cleanup.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual void variableCleanUp () = 0;

    /*@ManDoc: This is a virtual constructor for types derived
               from AmrLevel.  The derived type is initialized
               with the default constructor.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual AmrLevel* operator() () = 0;

    /*@ManDoc: This is a virtual constructor for types derived
               from AmrLevel.  The derived type is initialized
               with the five specified variables.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual AmrLevel* operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& geom_lev,
                                  const BoxArray& ba,
                                  Real            time) = 0;
};

#endif /*_LEVELBLD_H_*/
