
%{
#include <AMReX_IParser_Y.H>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int amrex_iparserlex (void);
/* Bison seems to have a bug. yyalloc etc. do not have the api.prefix. */
#ifndef yyalloc
#  define yyalloc amrex_iparseralloc
#endif
#ifndef yysymbol_kind_t
#  define yysymbol_kind_t amrex_iparsersymbol_kind_t
#endif
%}

/* We do not need to make this reentrant safe, because we use flex and
   bison for generating AST only and this part doesn't need to be
   thread safe.
*/
/*%define api.pure full */
%define api.prefix {amrex_iparser}

/* This is the type returned by functions iparser_new* declared in
   AMReX_IParser_Y.H.  See also bison rules at the end of this file.
*/
%union {
    struct amrex::iparser_node* n;
    long long d;
    struct amrex::iparser_symbol* s;
    enum amrex::iparser_f1_t f1;
    enum amrex::iparser_f2_t f2;
    enum amrex::iparser_f3_t f3;
}

/* Define tokens.  They are used by flex too. */
%token <n>  NODE
%token <d>  NUMBER
%token <s>  SYMBOL
%token <f1> F1
%token <f2> F2
%token <f3> F3
%token EOL
%token FLRDIV "//"
%token POW "**" '^'
%token GEQ ">="
%token LEQ "<="
%token EQ "=="
%token NEQ "!="
%token AND "and"
%token OR "or"

%left ';'
%nonassoc F1 F2 F3
%right '='
%left OR
%left AND
%left EQ NEQ
%left '<' '>' GEQ LEQ
%left '+' '-'
%left '*' '/' FLRDIV
%nonassoc NEG UPLUS
%right POW

/* This specifies the type of expressions */
%type <n> exp stmt or_exp and_exp cmp_exp add_exp mul_exp pow_exp unary_exp primary_exp

%start input

%%

/* Given `\n` terminated input, a tree is generated and passed to
 * function iparser_defexpr defined in AMReX_IParser_y.cpp.
 */
input:
  %empty
| input exp EOL {
    amrex::iparser_defexpr($2);
  }
;

/* Top level - handles lists and assignments */
exp:
  stmt                       { $$ = $1; }
| exp ';' stmt               { $$ = amrex::iparser_newlist($1, $3); }
| exp ';'                    { $$ = amrex::iparser_newlist($1, nullptr); }
;

/* Statements - handles assignments and expressions */
stmt:
  or_exp                     { $$ = $1; }
| SYMBOL '=' or_exp          { $$ = amrex::iparser_newassign($1, $3); }

/* OR expressions */
or_exp:
  and_exp                    { $$ = $1; }
| or_exp OR and_exp          { $$ = amrex::iparser_newf2(amrex::IPARSER_OR, $1, $3); }
;

/* AND expressions */
and_exp:
  cmp_exp                    { $$ = $1; }
| and_exp AND cmp_exp        { $$ = amrex::iparser_newf2(amrex::IPARSER_AND, $1, $3); }
;

/* Comparison expressions - handles all comparison operators and chaining */
cmp_exp:
  add_exp                    { $$ = $1; }
| cmp_exp '<' add_exp        { $$ = amrex::iparser_newcmpchain($1, amrex::IPARSER_LT, $3); }
| cmp_exp '>' add_exp        { $$ = amrex::iparser_newcmpchain($1, amrex::IPARSER_GT, $3); }
| cmp_exp LEQ add_exp        { $$ = amrex::iparser_newcmpchain($1, amrex::IPARSER_LEQ,$3); }
| cmp_exp GEQ add_exp        { $$ = amrex::iparser_newcmpchain($1, amrex::IPARSER_GEQ,$3); }
| cmp_exp EQ add_exp         { $$ = amrex::iparser_newcmpchain($1, amrex::IPARSER_EQ ,$3); }
| cmp_exp NEQ add_exp        { $$ = amrex::iparser_newcmpchain($1, amrex::IPARSER_NEQ,$3); }
;

/* Addition and subtraction */
add_exp:
  mul_exp                    { $$ = $1; }
| add_exp '+' mul_exp        { $$ = amrex::iparser_newnode(amrex::IPARSER_ADD, $1, $3); }
| add_exp '-' mul_exp        { $$ = amrex::iparser_newnode(amrex::IPARSER_SUB, $1, $3); }
;

/* Multiplication and division */
mul_exp:
  unary_exp                  { $$ = $1; }
| mul_exp '*' unary_exp      { $$ = amrex::iparser_newnode(amrex::IPARSER_MUL, $1, $3); }
| mul_exp '/' unary_exp      { $$ = amrex::iparser_newnode(amrex::IPARSER_DIV, $1, $3); }
| mul_exp FLRDIV unary_exp   { $$ = amrex::iparser_newf2(amrex::IPARSER_FLRDIV, $1, $3); };

/* Unary expressions */
unary_exp:
  pow_exp                    { $$ = $1; }
| '-' unary_exp              { $$ = amrex::iparser_newnode(amrex::IPARSER_NEG, $2, nullptr); }
| '+' unary_exp              { $$ = $2; }
;

/* Power (right associative) */
pow_exp:
  primary_exp                { $$ = $1; }
| primary_exp POW unary_exp  { $$ = amrex::iparser_newf2(amrex::IPARSER_POW, $1, $3); }
;

/* Primary expressions */
primary_exp:
  NUMBER                     { $$ = amrex::iparser_newnumber($1); }
| SYMBOL                     { $$ = amrex::iparser_newsymbol($1); }
| '(' or_exp ')'                { $$ = $2; }
| F1 '(' or_exp ')'             { $$ = amrex::iparser_newf1($1, $3); }
| F2 '(' or_exp ',' or_exp ')'     { $$ = amrex::iparser_newf2($1, $3, $5); }
| F3 '(' or_exp ',' or_exp ',' or_exp ')' { $$ = amrex::iparser_newf3($1, $3, $5, $7); }
;

%%
