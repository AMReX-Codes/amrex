
#ifndef BL_MFCOPYDESCRIPTOR_H_
#define BL_MFCOPYDESCRIPTOR_H_

#include <AMReX_FACopyDescriptor.H>
#include <AMReX_MultiFab.H>

namespace amrex {

class MultiFabCopyDescriptor;

typedef FabArrayId MultiFabId;

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


class MultiFabCopyDescriptor
    :
    public FabArrayCopyDescriptor<FArrayBox>
{
  public:

    MultiFabCopyDescriptor ();

    ~MultiFabCopyDescriptor ();

    MultiFabCopyDescriptor (const MultiFabCopyDescriptor&) = delete;
    MultiFabCopyDescriptor& operator= (const MultiFabCopyDescriptor&) = delete;

    MultiFabId RegisterMultiFab (MultiFab* mf) { return RegisterFabArray(mf); }
};

}

#endif
