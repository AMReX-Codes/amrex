#ifndef _MULTIMASK_H_
#define _MULTIMASK_H_

#include <Mask.H>
#include <FabArray.H>

class MultiMask
    :
    public FabArray<Mask>
{
public:
    MultiMask ();

private:
    //
    // These are disabled.
    //
    MultiMask (const MultiMask& rhs);
    MultiMask& operator= (const MultiMask& rhs);
};

#endif

