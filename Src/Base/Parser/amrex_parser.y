
%{
#include <AMReX_Parser_Y.H>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int amrex_parserlex (void);
/* Bison seems to have a bug. yyalloc etc. do not have the api.prefix. */
#ifndef yyalloc
#  define yyalloc amrex_parseralloc
#endif
#ifndef yysymbol_kind_t
#  define yysymbol_kind_t amrex_parsersymbol_kind_t
#endif
%}

/* We do not need to make this reentrant safe, because we use flex and
   bison for generating AST only and this part doesn't need to be
   thread safe.
*/
/*%define api.pure full */
%define api.prefix {amrex_parser}

/* This is the type returned by functions parser_new* declared in
   AMReX_Parser_y.H.  See also bison rules at the end of this file.
*/
%union {
    struct amrex::parser_node* n;
    double d;
    struct amrex::parser_symbol* s;
    enum amrex::parser_f1_t f1;
    enum amrex::parser_f2_t f2;
    enum amrex::parser_f3_t f3;
}

/* Define tokens.  They are used by flex too. */
%token <n>  NODE
%token <d>  NUMBER
%token <s>  SYMBOL
%token <f1> F1
%token <f2> F2
%token <f3> F3
%token EOL
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
%left '*' '/'
%nonassoc NEG UPLUS
%right POW

/* This specifies the type of expressions */
%type <n> exp stmt or_exp and_exp cmp_exp add_exp mul_exp pow_exp unary_exp primary_exp

%start input

%%

/* Given `\n` terminated input, a tree is generated and passed to
 * function parser_defexpr defined in AMReX_Parser_Y.cpp.
 */
input:
  %empty
| input exp EOL {
    amrex::parser_defexpr($2);
  }
;

/* Top level - handles lists and assignments */
exp:
  stmt                       { $$ = $1; }
| exp ';' stmt               { $$ = amrex::parser_newlist($1, $3); }
| exp ';'                    { $$ = amrex::parser_newlist($1, nullptr); }
;

/* Statements - handles assignments and expressions */
stmt:
  or_exp                     { $$ = $1; }
| SYMBOL '=' or_exp          { $$ = amrex::parser_newassign($1, $3); }

/* OR expressions */
or_exp:
  and_exp                    { $$ = $1; }
| or_exp OR and_exp          { $$ = amrex::parser_newf2(amrex::PARSER_OR, $1, $3); }
;

/* AND expressions */
and_exp:
  cmp_exp                    { $$ = $1; }
| and_exp AND cmp_exp        { $$ = amrex::parser_newf2(amrex::PARSER_AND, $1, $3); }
;

/* Comparison expressions - handles all comparison operators and chaining */
cmp_exp:
  add_exp                    { $$ = $1; }
| cmp_exp '<' add_exp        { $$ = amrex::parser_newcmpchain($1, amrex::PARSER_LT, $3); }
| cmp_exp '>' add_exp        { $$ = amrex::parser_newcmpchain($1, amrex::PARSER_GT, $3); }
| cmp_exp LEQ add_exp        { $$ = amrex::parser_newcmpchain($1, amrex::PARSER_LEQ,$3); }
| cmp_exp GEQ add_exp        { $$ = amrex::parser_newcmpchain($1, amrex::PARSER_GEQ,$3); }
| cmp_exp EQ add_exp         { $$ = amrex::parser_newcmpchain($1, amrex::PARSER_EQ ,$3); }
| cmp_exp NEQ add_exp        { $$ = amrex::parser_newcmpchain($1, amrex::PARSER_NEQ,$3); }
;

/* Addition and subtraction */
add_exp:
  mul_exp                    { $$ = $1; }
| add_exp '+' mul_exp        { $$ = amrex::parser_newnode(amrex::PARSER_ADD, $1, $3); }
| add_exp '-' mul_exp        { $$ = amrex::parser_newnode(amrex::PARSER_SUB, $1, $3); }
;

/* Multiplication and division */
mul_exp:
  unary_exp                  { $$ = $1; }
| mul_exp '*' unary_exp      { $$ = amrex::parser_newnode(amrex::PARSER_MUL, $1, $3); }
| mul_exp '/' unary_exp      { $$ = amrex::parser_newnode(amrex::PARSER_DIV, $1, $3); }
;

/* Unary expressions */
unary_exp:
  pow_exp                    { $$ = $1; }
| '-' unary_exp              { $$ = amrex::parser_newneg($2); }
| '+' unary_exp              { $$ = $2; }
;

/* Power (right associative) */
pow_exp:
  primary_exp                { $$ = $1; }
| primary_exp POW unary_exp  { $$ = amrex::parser_newf2(amrex::PARSER_POW, $1, $3); }
;

/* Primary expressions */
primary_exp:
  NUMBER                     { $$ = amrex::parser_newnumber($1); }
| SYMBOL                     { $$ = amrex::parser_newsymbol($1); }
| '(' or_exp ')'                { $$ = $2; }
| F1 '(' or_exp ')'             { $$ = amrex::parser_newf1($1, $3); }
| F2 '(' or_exp ',' or_exp ')'     { $$ = amrex::parser_newf2($1, $3, $5); }
| F3 '(' or_exp ',' or_exp ',' or_exp ')' { $$ = amrex::parser_newf3($1, $3, $5, $7); }
| SYMBOL '(' or_exp ')'                 { $$ = amrex::parser_newusrf1($1, $3); }
| SYMBOL '(' or_exp ',' or_exp ')'         { $$ = amrex::parser_newusrf2($1, $3, $5); }
| SYMBOL '(' or_exp ',' or_exp ',' or_exp ')' { $$ = amrex::parser_newusrfn($1, {$3, $5, $7}); }
| SYMBOL '(' or_exp ',' or_exp ',' or_exp ',' or_exp ')' { $$ = amrex::parser_newusrfn($1, {$3, $5, $7, $9}); }
;

%%
