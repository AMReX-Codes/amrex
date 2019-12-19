#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_

class ReducedDiags
{
public:

    /** constructor */
    ReducedDiags ();

    /** destructor */
    ~ReducedDiags ();

    /** output frequency of the reduced diagnostics */
    int m_freq;

    void ComputeDiags();
};

#endif //WARPX_DIAGNOSTICS_REDUCEDDIAGS_REDUCEDDIAGS_H_
