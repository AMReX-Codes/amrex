#ifndef Prob_H_
#define Prob_H_

#include "Adv_F.H"

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        // Nothing needs to be done here,
        // since there are no extra inputs to be read from probin file
    }
}
#endif
