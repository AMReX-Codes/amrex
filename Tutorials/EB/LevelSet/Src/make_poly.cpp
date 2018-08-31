
#include <AMReX_ParmParse.H>

#include <make_shapes.H>

using namespace amrex;



std::unique_ptr<CylinderIF>
make_poly_eb2_geom(int lev, int max_order, std::string field_prefix)
{
    // Construct the ParamParse database field names based on the
    // `field_prefix` string:
    ParmParse pp("poly");

    // Coefficients vector is stored in the inputs database with the field name:
    //      <field_prefix>_[x,y,z]_coeffs
    const std::array<const std::string, 3> var_names{"x", "y", "z"};
    std::array<std::string, 3> field_names;
    for(int i = 0; i < 3; i++) {
        std::stringstream field_name;
        field_name << field_prefix;
        field_name << "_" << var_names[i] << "_coeffs";
        field_names[i] = field_name.str();
    }

    // There are two more fields associated with the PolynomialIF:
    //      <field_prefix>_mirror    (true if fluid is inside the PolynomialIF)
    //      <field_prefix>_translate (vector representing center-axis position)
    std::stringstream mirror_field, translate_field;
    mirror_field << field_prefix << "_mirror";
    translate_field << field_prefix << "_translate";


    // Generate vector representing polynomial
    Vector<EB2::PolyTerm> poly;
    for(int idir = 0; idir < 3; idir++) {
        Vector<Real> coefvec(SpaceDim);

        if(idir == 0)      pp.getarr(field_names[idir].c_str(), coefvec, 0, max_order);
        else if(idir == 1) pp.getarr(field_names[idir].c_str(), coefvec, 0, max_order);
        else if(idir == 2) pp.getarr(field_names[idir].c_str(), coefvec, 0, max_order);

        for(int lc = 0; lc < max_order; lc++) {
            // x^(lc) term
            Real coef = coefvec[lc];
            IntVect powers = IntVect::Zero;
            powers[idir] = lc;

            EB2::PolyTerm mono = {.coef = coef, .powers = powers};
            poly.push_back(mono);
        }
    }


    /****************************************************************************
     * Construct PolynomialIF (called CylinderIF here when used within a        *
     * TranslationIF, see make_shapes.H )                                       *
     ****************************************************************************/

    bool flip = true;
    pp.query(mirror_field.str().c_str(), flip);

    Vector<Real> transvec(SpaceDim);
    pp.getarr(translate_field.str().c_str(), transvec, 0, SpaceDim);

    RealArray offset = {AMREX_D_DECL(transvec[0], transvec[1], transvec[2])};
    CylinderIF poly2(EB2::PolynomialIF(poly, flip), offset);

    std::unique_ptr<CylinderIF> ret = std::unique_ptr<CylinderIF>(new CylinderIF(poly2));

    return ret;

}
