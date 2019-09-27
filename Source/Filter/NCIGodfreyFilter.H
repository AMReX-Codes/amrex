#include <Filter.H>

#ifndef WARPX_GODFREY_FILTER_H_
#define WARPX_GODFREY_FILTER_H_

enum class godfrey_coeff_set { Ex_Ey_Bz=0, Bx_By_Ez=1 };

/**
 * \brief Class for Godrey's filter to suppress Numerical Cherenkov Instability
 *
 * It derives from the base class Filter.
 * The filter stencil is initialized in method ComputeStencils. Computing the
 * stencil requires to read parameters from a table, which is done in
 * getGodfreyCoeffs.
 * The filter is applied using the base class' method ApplyStencil.
 */
class NCIGodfreyFilter : public Filter
{
public:

    NCIGodfreyFilter () = default;

    NCIGodfreyFilter(godfrey_coeff_set coeff_set_, amrex::Real cdtodz_,
                     amrex::Real l_lower_order_in_v_, bool nodal_gather_);

    void ComputeStencils();

    void getGodfreyCoeffs(godfrey_coeff_set coeff_set_in);

    static constexpr int stencil_width = 4;

private:

    // Set of coefficients (different fields require to read
    // different coefficients from the table)
    godfrey_coeff_set coeff_set;
    // The stencil depends on parameter c*dt/dz
    amrex::Real cdtodz;
    // If 1, use Galerkin method for gather. If 0, uniform gather
    int l_lower_order_in_v;
    // Whether the gather is from nodal fields or staggered fields
    bool nodal_gather;
};

#endif // #ifndef WARPX_GODFREY_FILTER_H_
