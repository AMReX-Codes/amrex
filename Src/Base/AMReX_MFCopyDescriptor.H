
#ifndef BL_MFCOPYDESCRIPTOR_H_
#define BL_MFCOPYDESCRIPTOR_H_
#include <AMReX_Config.H>

#include <AMReX_FACopyDescriptor.H>
#include <AMReX_MultiFab.H>

/**
 * \file AMReX_MFCopyDescriptor.H
 * \brief MultiFab-specific wrappers around FabArrayCopyDescriptor.
 */

namespace amrex {

class MultiFabCopyDescriptor;

using MultiFabId = FabArrayId;

/**
 * \brief Add an interpolation request spanning two source MultiFabs.
 *
 * Populates \p returnedFillBoxIds describing the fill operations required to
 * interpolate data from MultiFabs \p faid1 and \p faid2 at time \p t.
 */
void InterpAddBox (MultiFabCopyDescriptor& fabCopyDesc,
                   BoxList*                returnUnfilledBoxes,
                   Vector<FillBoxId>&       returnedFillBoxIds,
                   const Box&              subbox,
                   MultiFabId              faid1,
                   MultiFabId              faid2,
                   Real                    t1,
                   Real                    t2,
                   Real                    t,
                   int                     src_comp,
                   int                     dest_comp,
                   int                     num_comp,
                   bool                    extrap);

/**
 * \brief Execute interpolation fills previously staged via InterpAddBox.
 */
void InterpFillFab (MultiFabCopyDescriptor& fabCopyDesc,
                    const Vector<FillBoxId>& fillBoxIds,
                    MultiFabId              faid1,
                    MultiFabId              faid2,
                    FArrayBox&              dest,
                    Real                    t1,
                    Real                    t2,
                    Real                    t,
                    int                     src_comp,
                    int                     dest_comp,
                    int                     num_comp,
                    bool                    extrap);


/**
 * \class amrex::MultiFabCopyDescriptor
 * \brief Thin wrapper that specializes FabArrayCopyDescriptor for FArrayBox.
 */
class MultiFabCopyDescriptor
    :
    public FabArrayCopyDescriptor<FArrayBox>
{
  public:

    MultiFabCopyDescriptor () = default;
    ~MultiFabCopyDescriptor () = default;

    MultiFabCopyDescriptor (const MultiFabCopyDescriptor&) = delete;
    MultiFabCopyDescriptor (MultiFabCopyDescriptor&&) = delete;
    MultiFabCopyDescriptor& operator= (const MultiFabCopyDescriptor&) = delete;
    MultiFabCopyDescriptor& operator= (MultiFabCopyDescriptor&&) = delete;

    MultiFabId RegisterMultiFab (MultiFab* mf) { return RegisterFabArray(mf); }
};

}

#endif
