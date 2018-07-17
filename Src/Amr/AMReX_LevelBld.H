
#ifndef AMREX_LEVELBLD_H_
#define AMREX_LEVELBLD_H_

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>

namespace amrex {

//
// Builds problem-specific AmrLevels.
//
// Abstract base class specifying an interface for building problem-specific
// AmrLevels.  
//

class LevelBld
{
public:

    virtual ~LevelBld () {}
    //
    // Perform any problem-dependent setup such as physical
    // boundary condition and derived quantities.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual void variableSetUp () = 0;
    virtual void variableSetUpForNewCompProcs () { } 
    //
    // Perform any problem-dependent cleanup.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual void variableCleanUp () = 0;
    //
    // This is a virtual constructor for types derived
    // from AmrLevel.  The derived type is initialized
    // with the default constructor.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual AmrLevel* operator() () = 0;
    //
    // This is a virtual constructor for types derived
    // from AmrLevel.  The derived type is initialized
    // with the five specified variables.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual AmrLevel* operator() (Amr&                       papa,
                                  int                        lev,
                                  const Geometry&            geom_lev,
                                  const BoxArray&            ba,
                                  const DistributionMapping& dm,
                                  Real                       time) = 0;
};

}

extern amrex::LevelBld* getLevelBld ();

#endif /*_LEVELBLD_H_*/
