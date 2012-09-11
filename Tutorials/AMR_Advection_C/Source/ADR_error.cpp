#include <winstd.H>

#include "ADR.H"
#include "ADR_F.H"

using std::string;

static Box the_same_box (const Box& b) { return b; }
static Box grow_box_by_one (const Box& b) { return BoxLib::grow(b,1); }
static Box grow_box_by_two (const Box& b) { return BoxLib::grow(b,2); }

typedef StateDescriptor::BndryFunc BndryFunc;

void
ADR::ErrorSetUp ()
{
    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    err_list.add("StateErr",1,ErrorRec::Special,
                 BL_FORT_PROC_CALL(STATE_ERROR,state_error));

}
