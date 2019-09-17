%option noyywrap nodefault
%{
#include "wp_parser_y.h"
#include "wp_parser.tab.h"
%}

 /* Tokens NUMBER, SYMBOL, F1, POW, F2, etc. are defined in wp_parser.y. */
 /* Types WP_SQRT, WP_SQRT, etc. are defined in wp_parser_y.h. */

 /* Used leater to define NUMBER */
EXP	([Ee][-+]?[0-9]+)

%%

"+" |
"-" |
"*" |
"/" |
"=" |
"|" |
"," |
"<" |
">" |
"(" |
")"     { return yytext[0]; } /* simply pass through */

 /* yylval is union type defined in wp_parser.tab.h that is generated
  * by bison with wp_parser.y */

"sqrt"        { yylval.f1 = WP_SQRT;      return F1; }
"exp"         { yylval.f1 = WP_EXP;       return F1; }
"log"         { yylval.f1 = WP_LOG;       return F1; }
"log10"       { yylval.f1 = WP_LOG10;     return F1; }
"sin"         { yylval.f1 = WP_SIN;       return F1; }
"cos"         { yylval.f1 = WP_COS;       return F1; }
"tan"         { yylval.f1 = WP_TAN;       return F1; }
"asin"        { yylval.f1 = WP_ASIN;      return F1; }
"acos"        { yylval.f1 = WP_ACOS;      return F1; }
"atan"        { yylval.f1 = WP_ATAN;      return F1; }
"sinh"        { yylval.f1 = WP_SINH;      return F1; }
"cosh"        { yylval.f1 = WP_COSH;      return F1; }
"tanh"        { yylval.f1 = WP_TANH;      return F1; }
"abs"         { yylval.f1 = WP_ABS;       return F1; }
"fabs"        { yylval.f1 = WP_ABS;       return F1; }
"**"          { yylval.f2 = WP_POW;       return POW;}
"^"           { yylval.f2 = WP_POW;       return POW;}
">="          { yylval.f2 = WP_GEQ;       return GEQ;}
"<="          { yylval.f2 = WP_LEQ;       return LEQ;}
"=="          { yylval.f2 = WP_EQ;        return EQ;}
"!="          { yylval.f2 = WP_NEQ;       return NEQ;}
"and"         { yylval.f2 = WP_AND;       return AND;}
"or"          { yylval.f2 = WP_OR;        return OR;}
"pow"         { yylval.f2 = WP_POW;       return F2; }
"heaviside"   { yylval.f2 = WP_HEAVISIDE; return F2; }
"min"         { yylval.f2 = WP_MIN;       return F2; }
"max"         { yylval.f2 = WP_MAX;       return F2; }

 /* We use SYMBOL to hold variables and constants */
[a-zA-Z_][a-zA-Z0-9_]*  { yylval.s = wp_makesymbol(yytext); return SYMBOL; }

 /* Number */
[0-9]+"."[0-9]*{EXP}? |
"."?[0-9]+{EXP}?  { yylval.d = atof(yytext);   return NUMBER; }

 /* Special characters */
"//".*
[ \t]   /* ignore white space */
\\\n    /* ignore line continuation */
"\n"    { return EOL; }

 /* everything else */
.	{ yyerror("Unknow character %c\n", *yytext); }

%%
