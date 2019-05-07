#include <Filter.H>
#include <NCIGodfreyTables.H>

#ifndef WARPX_GODFREY_FILTER_H_
#define WARPX_GODFREY_FILTER_H_

enum class godfrey_coeff_set { Ex_Ey_Bz=0, Bx_By_Ez=1 };

class NCIGodfreyFilter : public Filter
{
public:

    NCIGodfreyFilter(godfrey_coeff_set coeff_set_in)
        : Filter(),
        coeff_set(coeff_set_in)
    {
        if        (coeff_set == godfrey_coeff_set::Ex_Ey_Bz){
            coeff_table = (amrex::Real*) table_nci_godfrey_Ex_Ey_Bz;
        } else if (coeff_set == godfrey_coeff_set::Bx_By_Ez){
            coeff_table = (amrex::Real*) table_nci_godfrey_Bx_By_Ez;
        }
    };

    virtual void ComputeStencils() override final;

    void getGodfreyCoeffs(godfrey_coeff_set coeff_set_in);

private:
    
    godfrey_coeff_set coeff_set;
    amrex::Real* coeff_table;

};

#endif // #ifndef WARPX_GODFREY_FILTER_H_
