#include <Filter.H>

#ifndef WARPX_GODFREY_FILTER_H_
#define WARPX_GODFREY_FILTER_H_

enum class godfrey_coeff_set { Ex_Ey_Bz=0, Bx_By_Ez=1 };

class NCIGodfreyFilter : public Filter
{
public:

    NCIGodfreyFilter () = default;

    NCIGodfreyFilter(godfrey_coeff_set coeff_set_, amrex::Real cdtodz_, amrex::Real l_lower_order_in_v_);

    void ComputeStencils();

    void getGodfreyCoeffs(godfrey_coeff_set coeff_set_in);

    static constexpr int stencil_width = 4;

private:
    
    godfrey_coeff_set coeff_set;
    amrex::Real cdtodz;
    int l_lower_order_in_v;

};

#endif // #ifndef WARPX_GODFREY_FILTER_H_
