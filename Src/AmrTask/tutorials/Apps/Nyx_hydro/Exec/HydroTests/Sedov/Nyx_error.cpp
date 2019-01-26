
#include "Nyx.H"
#include "Nyx_error_F.H"

using namespace amrex;

using std::string;

void
Nyx::error_setup()
{
    // The lines below define routines to be called to tag cells for error
    // estimation -- the arguments of each "add" call are:
    //   1. Name of variable (state variable or derived quantity) which will be
    //      passed into the Fortran subroutine.
    //   2. Number of ghost cells each array needs in each call to the Fortran
    //      subroutine.
    //   3. Type of Fortran subroutine -- this determines the argument list of
    //      the Fortran subroutine.  These types are pre-defined and are
    //      currently restricted to ErrorRec::Standard and ErrorRec::UseAverage.
    //   4. Name of Fortran subroutine.

    // The routine LAPLAC_ERROR uses the special evaluation of the second
    // derivative and can be called with any variable. Note that two ghost cells
    // are needed.
    err_list.add("density",1,ErrorRec::Standard,
                 BL_FORT_PROC_CALL(TAG_DENERROR,tag_denerror));
    err_list.add("pressure",1,ErrorRec::Standard,
                 BL_FORT_PROC_CALL(TAG_PRESSERROR,tag_presserror));
    err_list.add("x_velocity",1,ErrorRec::Standard,
                 BL_FORT_PROC_CALL(TAG_VELERROR,tag_velerror));
}

void
Nyx::manual_tags_placement (TagBoxArray&    tags,
                            const Vector<IntVect>& bf_lev)
{
}
