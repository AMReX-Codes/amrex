
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

/* This specifies the type of `exp` (i.e., struct parser_node*).  Rules
   specified later pass `exp` to parser_new* functions declared in
   AMReX_Parser_Y.H.
*/
%type <n> exp

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

/* Enum types PARSER_ADD, PARSER_SUB, etc. are defined in AMReX_Parser_Y.H
 * Functions parser_new* are also declared in that file.
 */
exp:
  NUMBER                     { $$ = amrex::parser_newnumber($1); }
| SYMBOL                     { $$ = amrex::parser_newsymbol($1); }
| exp '+' exp                { $$ = amrex::parser_newnode(amrex::PARSER_ADD, $1, $3); }
| exp '-' exp                { $$ = amrex::parser_newnode(amrex::PARSER_SUB, $1, $3); }
| exp '*' exp                { $$ = amrex::parser_newnode(amrex::PARSER_MUL, $1, $3); }
| exp '/' exp                { $$ = amrex::parser_newnode(amrex::PARSER_DIV, $1, $3); }
| '(' exp ')'                { $$ = $2; }
| exp '<' exp                { $$ = amrex::parser_newf2(amrex::PARSER_LT, $1, $3); }
| exp '>' exp                { $$ = amrex::parser_newf2(amrex::PARSER_GT, $1, $3); }
| exp LEQ exp                { $$ = amrex::parser_newf2(amrex::PARSER_LEQ, $1, $3); }
| exp GEQ exp                { $$ = amrex::parser_newf2(amrex::PARSER_GEQ, $1, $3); }
| exp EQ exp                 { $$ = amrex::parser_newf2(amrex::PARSER_EQ, $1, $3); }
| exp NEQ exp                { $$ = amrex::parser_newf2(amrex::PARSER_NEQ, $1, $3); }
| exp AND exp                { $$ = amrex::parser_newf2(amrex::PARSER_AND, $1, $3); }
| exp OR exp                 { $$ = amrex::parser_newf2(amrex::PARSER_OR, $1, $3); }
| '-'exp %prec NEG           { $$ = amrex::parser_newnode(amrex::PARSER_NEG, $2, nullptr); }
| '+'exp %prec UPLUS         { $$ = $2; }
| exp POW exp                { $$ = amrex::parser_newf2(amrex::PARSER_POW, $1, $3); }
| F1 '(' exp ')'             { $$ = amrex::parser_newf1($1, $3); }
| F2 '(' exp ',' exp ')'     { $$ = amrex::parser_newf2($1, $3, $5); }
| F3 '(' exp ',' exp ',' exp ')' { $$ = amrex::parser_newf3($1, $3, $5, $7); }
| SYMBOL '=' exp             { $$ = amrex::parser_newassign($1, $3); }
| exp ';' exp                { $$ = amrex::parser_newlist($1, $3); }
| exp ';'                    { $$ = amrex::parser_newlist($1, nullptr); }
;

%%
