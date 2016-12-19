
#ifndef _MACOPERATOR_H_
#define _MACOPERATOR_H_

//==========================================================================
// This class defines how mass continuity is defined for the projection operators
// in class NavierStokes.  For the Mac project, this is done through setting a
// linop and a few driver functions, while for the mac project this is done
// through a few driver functions which scale variables.
//==========================================================================

#include <AMReX_ABecLaplacian.H>
#include <AMReX_BndryData.H>
#include <AMReX_AmrCore.H>

class MacOperator
    :
    public amrex::ABecLaplacian
{
protected:

    amrex::AmrCore*      parent;
public:

    MacOperator (amrex::AmrCore*             Parent,
                 const amrex::BndryData& mgb,
                 const Real*      _h);

    virtual ~MacOperator ();
    //
    // Define the gradient coeficients used by the multigrid object.
    //
    void setCoefficients (amrex::MultiFab*   area,
                          amrex::MultiFab&   rho,
                          int         rho_comp,
                          const Real* dx);
    //
    // Initialize the rhs for the level solves.
    //
    void defRHS (amrex::MultiFab* area,
                 amrex::MultiFab& vol,
                 amrex::MultiFab& rhs,
                 amrex::MultiFab* vel,
                 Real      scale);
    //
    // Apply the pressure gradient to the mac velocities.
    //
    void velUpdate (amrex::MultiFab*       vel,
                    amrex::MultiFab&       phi,
                    const amrex::MultiFab& Rho,
                    int             rho_comp,
                    const Real*     dx,
                    Real            scale);
    //
    // Scale the rhs for the sync solves.
    //
    void syncRhs (const amrex::MultiFab& Volume,
                  amrex::MultiFab&       Rhs,
                  Real            rhs_scale,
                  const Real*     dx);

    static void Initialize ();
    static void Finalize ();
};

#endif
