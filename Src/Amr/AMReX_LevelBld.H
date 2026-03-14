
#ifndef AMREX_LEVELBLD_H_
#define AMREX_LEVELBLD_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>

/**
 * \file AMReX_LevelBld.H
 * \brief Factory interface for creating problem-specific AmrLevel objects.
 */

namespace amrex {


/**
* \brief Builds problem-specific AmrLevels.
*
* Abstract base class specifying an interface for building problem-specific
* AmrLevels.
*/
class LevelBld
{
public:

    //! Polymorphic base classes need virtual destructors.
    virtual ~LevelBld () = default;
    //! Defaulted special members keep the factory trivially constructible/copyable.
    LevelBld () noexcept = default;
    LevelBld (LevelBld const&) = default;
    LevelBld (LevelBld &&) noexcept = default;
    LevelBld& operator= (LevelBld const&) = default;
    LevelBld& operator= (LevelBld &&) noexcept = default;

    /**
    * \brief Perform any problem-dependent setup such as physical
    * boundary condition and derived quantities.
    * This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    */
    virtual void variableSetUp () = 0;
    //! Optional hook invoked when new compute processes join (multi-job setups).
    virtual void variableSetUpForNewCompProcs () { }

    /**
    * \brief Perform any problem-dependent cleanup.
    * This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    */
    virtual void variableCleanUp () = 0;

    /**
    * \brief This is a virtual constructor for types derived
    * from AmrLevel.  The derived type is initialized
    * with the default constructor.
    * This is a pure virtual function and hence MUST
    * be implemented by derived classes.
    */
    virtual AmrLevel* operator() () = 0;

    /**
     * \brief Create a new AmrLevel instance for level \p lev using the supplied metadata.
     *
     * \param papa     Owning Amr object.
     * \param lev      Level index being constructed.
     * \param geom_lev Geometry describing the level.
     * \param ba       Grid layout for the level.
     * \param dm       Distribution mapping for the level.
     * \param time     Simulation time at which the level is initialized.
     */
    virtual AmrLevel* operator() (Amr&                       papa,
                                  int                        lev,
                                  const Geometry&            geom_lev,
                                  const BoxArray&            ba,
                                  const DistributionMapping& dm,
                                  Real                       time) = 0;
};

}

#endif /*_LEVELBLD_H_*/
