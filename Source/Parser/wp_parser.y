
%{
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include "wp_parser_y.h"
    int yylex (void);
%}

/* We do not need to make this reentrant safe, because we use flex and
   bison for generating AST only and this part doesn't need to be
   thread safe.
*/
/*%define api.pure full */

/* This is the type returned by functions wp_new* declared in
   wp_parser_y.h.  See also bison rules at the end of this file.
*/
%union {
    struct wp_node* n;
    amrex_real d;
    struct wp_symbol* s;
    enum wp_f1_t f1;
    enum wp_f2_t f2;
}

/* Define tokens.  They are used by flex too. */
%token <n>  NODE
%token <d>  NUMBER
%token <s>  SYMBOL
%token <f1> F1
%token <f2> F2
%token EOL
%token POW "**" '^'
%token GEQ ">="
%token LEQ "<="
%token EQ "=="
%token NEQ "!="
%token AND "and"
%token OR "or"

%nonassoc F1 F2
%right '='
%left OR
%left AND
%left EQ NEQ
%left '<' '>' GEQ LEQ
%left '+' '-'
%left '*' '/'
%nonassoc NEG UPLUS
%right POW

/* This specifies the type of `exp` (i.e., struct wp_node*).  Rules
   specified later pass `exp` to wp_new* functions declared in
   wp_parser_y.h.
*/
%type <n> exp

%start input

%%

/* Given `\n` terminated input, a tree is generated and passed to
 * function wp_defexpr defined in wp_parser_y.c.
 */
input:
  %empty
| input exp EOL {
    wp_defexpr($2);
  }
;

/* Enum types WP_ADD, WP_SUB, etc. are defined in wp_parser_y.h.
 * Functions wp_new* are also declared in that file.
 */
exp:
  NUMBER                     { $$ = wp_newnumber($1); }
| SYMBOL                     { $$ = wp_newsymbol($1); }
| exp '+' exp                { $$ = wp_newnode(WP_ADD, $1, $3); }
| exp '-' exp                { $$ = wp_newnode(WP_SUB, $1, $3); }
| exp '*' exp                { $$ = wp_newnode(WP_MUL, $1, $3); }
| exp '/' exp                { $$ = wp_newnode(WP_DIV, $1, $3); }
| '(' exp ')'                { $$ = $2; }
| exp '<' exp                { $$ = wp_newf2(WP_LT, $1, $3); }
| exp '>' exp                { $$ = wp_newf2(WP_GT, $1, $3); }
| exp LEQ exp                { $$ = wp_newf2(WP_LEQ, $1, $3); }
| exp GEQ exp                { $$ = wp_newf2(WP_GEQ, $1, $3); }
| exp EQ exp                 { $$ = wp_newf2(WP_EQ, $1, $3); }
| exp NEQ exp                { $$ = wp_newf2(WP_NEQ, $1, $3); }
| exp AND exp                { $$ = wp_newf2(WP_AND, $1, $3); }
| exp OR exp                 { $$ = wp_newf2(WP_OR, $1, $3); }
| '-'exp %prec NEG           { $$ = wp_newnode(WP_NEG, $2, NULL); }
| '+'exp %prec UPLUS         { $$ = $2; }
| exp POW exp                { $$ = wp_newf2(WP_POW, $1, $3); }
| F1 '(' exp ')'             { $$ = wp_newf1($1, $3); }
| F2 '(' exp ',' exp ')'     { $$ = wp_newf2($1, $3, $5); }
;

%%
