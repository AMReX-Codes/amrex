
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
    int d;
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

/* This specifies the type of `exp` (i.e., struct iparser_node*).  Rules
   specified later pass `exp` to iparser_new* functions declared in
   AMReX_IParser_y.H.
*/
%type <n> exp

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

/* Enum types IPARSER_ADD, IPARSER_SUB, etc. are defined in AMReX_IParser_Y.H
 * Functions iparser_new* are also declared in that file.
 */
exp:
  NUMBER                     { $$ = amrex::iparser_newnumber($1); }
| SYMBOL                     { $$ = amrex::iparser_newsymbol($1); }
| exp '+' exp                { $$ = amrex::iparser_newnode(amrex::IPARSER_ADD, $1, $3); }
| exp '-' exp                { $$ = amrex::iparser_newnode(amrex::IPARSER_SUB, $1, $3); }
| exp '*' exp                { $$ = amrex::iparser_newnode(amrex::IPARSER_MUL, $1, $3); }
| exp '/' exp                { $$ = amrex::iparser_newnode(amrex::IPARSER_DIV, $1, $3); }
| exp FLRDIV exp             { $$ = amrex::iparser_newf2(amrex::IPARSER_FLRDIV, $1, $3); }
| '(' exp ')'                { $$ = $2; }
| exp '<' exp                { $$ = amrex::iparser_newf2(amrex::IPARSER_LT, $1, $3); }
| exp '>' exp                { $$ = amrex::iparser_newf2(amrex::IPARSER_GT, $1, $3); }
| exp LEQ exp                { $$ = amrex::iparser_newf2(amrex::IPARSER_LEQ, $1, $3); }
| exp GEQ exp                { $$ = amrex::iparser_newf2(amrex::IPARSER_GEQ, $1, $3); }
| exp EQ exp                 { $$ = amrex::iparser_newf2(amrex::IPARSER_EQ, $1, $3); }
| exp NEQ exp                { $$ = amrex::iparser_newf2(amrex::IPARSER_NEQ, $1, $3); }
| exp AND exp                { $$ = amrex::iparser_newf2(amrex::IPARSER_AND, $1, $3); }
| exp OR exp                 { $$ = amrex::iparser_newf2(amrex::IPARSER_OR, $1, $3); }
| '-'exp %prec NEG           { $$ = amrex::iparser_newnode(amrex::IPARSER_NEG, $2, nullptr); }
| '+'exp %prec UPLUS         { $$ = $2; }
| exp POW exp                { $$ = amrex::iparser_newf2(amrex::IPARSER_POW, $1, $3); }
| F1 '(' exp ')'             { $$ = amrex::iparser_newf1($1, $3); }
| F2 '(' exp ',' exp ')'     { $$ = amrex::iparser_newf2($1, $3, $5); }
| F3 '(' exp ',' exp ',' exp ')' { $$ = amrex::iparser_newf3($1, $3, $5, $7); }
| SYMBOL '=' exp             { $$ = amrex::iparser_newassign($1, $3); }
| exp ';' exp                { $$ = amrex::iparser_newlist($1, $3); }
| exp ';'                    { $$ = amrex::iparser_newlist($1, nullptr); }
;

%%
