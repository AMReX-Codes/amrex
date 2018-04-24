#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>

#include <AMReX_ParmParse.H>
#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <ParticleContainer.H>

std::string UserConstants::replaceStringValue(std::string math_expr){
    std::string pattern, value_str;
    std::string patternexp = "e+";
    std::size_t found;
    amrex::Real value;
    for (int i=0; i<nb_constants; i++){
        pattern    = constant_names[i];
        value      = constant_values[i];


        // Convert value to string, with scientific notation
        std::ostringstream streamObj;
        streamObj << value;
        value_str = streamObj.str();

        // Replace "e+" by "e" in scientific notation for Fortran compatibility
        found = value_str.find(patternexp);
        if (found != std::string::npos){
            value_str.replace(found,patternexp.length(),"e");
        }

        // Replace all occurences of variable pattern by a string with its value value_str
        found = math_expr.find(pattern);
        while (found != std::string::npos){
            if ((found==0               
                        && !isalnum(math_expr[pattern.length()      ])) ||
               (!isalnum(math_expr[found-1])  
                        && !isalnum(math_expr[found+pattern.length()]))){ 
                math_expr.replace(found,pattern.length(),value_str);
             }
            found = math_expr.find(pattern, found + pattern.length());
        }
    }
    return math_expr;
}

void UserConstants::ReadParameters()
{
    if (!initialized){
	    amrex::ParmParse pp("constants");
        pp.query("use_my_constants", use_my_constants);
        if (use_my_constants){
            pp.getarr("constant_names", constant_names);
            nb_constants = constant_names.size();
            pp.getarr("constant_values", constant_values);
            BL_ASSERT(constant_values.size() == nb_constants);
        }
	initialized = true;
    }
}
