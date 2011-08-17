
#ifndef _MACOPERATOR_H_
#define _MACOPERATOR_H_

//==========================================================================
// This class defines how mass continuity is defined for the projection operators
// in class NavierStokes.  For the Mac project, this is done through setting a
// linop and a few driver functions, while for the mac project this is done
// through a few driver functions which scale variables.
//==========================================================================

#include <AmrLevel.H>
#include <ABecLaplacian.H>
#include <BndryData.H>

class MacOperator
    :
    public ABecLaplacian
{
protected:

    Amr*        parent;
public:

    MacOperator (Amr*             Parent,
                 const BndryData& mgb,
                 const Real*      _h);

    virtual ~MacOperator ();
    //
    // Define the gradient coeficients used by the multigrid object.
    //
    void setCoefficients (MultiFab*   area,
                          MultiFab&   rho,
                          int         rho_comp,
                          const Real* dx);
    //
    // Initialize the rhs for the level solves.
    //
    void defRHS (MultiFab* area,
                 MultiFab& vol,
                 MultiFab& rhs,
                 MultiFab* vel,
                 Real      scale);
    //
    // Apply the pressure gradient to the mac velocities.
    //
    void velUpdate (MultiFab*       vel,
                    MultiFab&       phi,
                    const MultiFab& Rho,
                    int             rho_comp,
                    const Real*     dx,
                    Real            scale);
    //
    // Scale the rhs for the sync solves.
    //
    void syncRhs (const MultiFab& Volume,
                  MultiFab&       Rhs,
                  Real            rhs_scale,
                  const Real*     dx);

    static void Initialize ();
    static void Finalize ();
};

#endif
