%option noyywrap nodefault
%option never-interactive
%{
#include "AMReX_IParser_Y.H"
#include "amrex_iparser.tab.h"
#ifdef _WIN32
#define YY_NO_UNISTD_H
#endif
%}

%option nounput
%option noinput
%option prefix="amrex_iparser"

 /* Tokens NUMBER, SYMBOL, F1, POW, F2, etc. are defined in amrex_iparser.y. */
 /* Types IPARSER_ABS, IPARSER_MIN etc. are defined in AMReX_IParser_y.H. */

%%

"+" |
"-" |
"*" |
"/" |
"=" |
"," |
"<" |
">" |
";" |
"(" |
")"     { return amrex_iparsertext[0]; } /* simply pass through */

 /* amrex_iparserlval is union type defined in amrex_iparser.tab.h that is generated
  * by bison with amrex_iparser.y */

"abs"         { amrex_iparserlval.f1 = amrex::IPARSER_ABS;       return F1; }
"iabs"        { amrex_iparserlval.f1 = amrex::IPARSER_ABS;       return F1; }
"//"          { amrex_iparserlval.f2 = amrex::IPARSER_FLRDIV;    return FLRDIV; }
"**"          { amrex_iparserlval.f2 = amrex::IPARSER_POW;       return POW;}
"^"           { amrex_iparserlval.f2 = amrex::IPARSER_POW;       return POW;}
">="          { amrex_iparserlval.f2 = amrex::IPARSER_GEQ;       return GEQ;}
"<="          { amrex_iparserlval.f2 = amrex::IPARSER_LEQ;       return LEQ;}
"=="          { amrex_iparserlval.f2 = amrex::IPARSER_EQ;        return EQ;}
"!="          { amrex_iparserlval.f2 = amrex::IPARSER_NEQ;       return NEQ;}
"and"         { amrex_iparserlval.f2 = amrex::IPARSER_AND;       return AND;}
"or"          { amrex_iparserlval.f2 = amrex::IPARSER_OR;        return OR;}
"min"         { amrex_iparserlval.f2 = amrex::IPARSER_MIN;       return F2; }
"max"         { amrex_iparserlval.f2 = amrex::IPARSER_MAX;       return F2; }
"if"          { amrex_iparserlval.f3 = amrex::IPARSER_IF;        return F3; }

 /* We use SYMBOL to hold variables and constants */
[a-zA-Z_][a-zA-Z0-9_]*  { amrex_iparserlval.s = amrex::iparser_makesymbol(amrex_iparsertext); return SYMBOL; }

 /* Number */
[0-9]+        { amrex_iparserlval.d = std::atoi(amrex_iparsertext);   return NUMBER; }

 /* Special characters */
[ \t]   /* ignore white space */
\\\n    /* ignore line continuation */
"\n"    { return EOL; }

 /* everything else */
.	{ amrex_iparsererror("Unknown character %c", *amrex_iparsertext); }

%%
