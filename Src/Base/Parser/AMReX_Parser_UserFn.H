// IWYU pragma: private, include "AMReX_Parser.H"
#ifndef AMREX_PARSER_USERFN_H_
#define AMREX_PARSER_USERFN_H_

namespace amrex {

using ParserUserFn1 = double (*)(double);
using ParserUserFn2 = double (*)(double,double);
using ParserUserFn3 = double (*)(double,double,double);
using ParserUserFn4 = double (*)(double,double,double,double);

}

#endif
