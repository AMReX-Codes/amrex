#ifndef AMR_Adv_C_V2_Adv_H_
#define AMR_Adv_C_V2_Adv_H_

#include <AmrCore.H>

class AmrAdv
    : public AmrCore
{
public:
    AmrAdv ();
    virtual ~AmrAdv();

    void InitData ();
    void Evolve ();

protected:
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override;
};

#endif
