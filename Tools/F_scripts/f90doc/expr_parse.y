%{
package expr_parse;

# On failure, print out this as the line we were working on.
$expr_parse::line = "";

# Portion of line left to parse
$expr_parse::left = "";
%}

%token COMMA LPAREN RPAREN NOT OR AND EQV NEQV COMPARISON DBLSLASH PERCENT
%token PLUS MINUS UPLUS UMINUS ASTERIK SLASH DBLASTERIK CONST NAME COLON
%token LARRAY RARRAY EQUALS

%left EQV NEQV
%left OR
%left AND
%nonassoc NOT
%nonassoc COMPARISON
%left DBLSLASH
%left PLUS MINUS
%nonassoc UPLUS UMINUS
%left ASTERIK SLASH
%right DBLASTERIK
%left PERCENT

%%

expr_with_abort: expr           { $$ = $1; return 1; }
    | expr COMMA                { $$ = $1; return "s,"; }

expr:
    CONST                         { $$ = [ "%const", @{$1} ]; }
  | expr_without_const            { $$ = $1; }

expr_without_const:
    chain                         { $$ = $1; }
  | LARRAY array RARRAY           { $$ = [ "%array", @{$2} ]; }
  | PLUS expr %prec UPLUS         { $$ = [ "u+", $2 ]; }
  | MINUS expr %prec UMINUS       { $$ = [ "u-", $2 ]; }
  | NOT expr                      { $$ = [ $1, $2 ]; }
  | LPAREN potential_complex_or_implied_do RPAREN
                                  { $$ = $2; }
  | expr DBLASTERIK expr          { $$ = [ $2, $1, $3 ]; }
  | expr ASTERIK    expr          { $$ = [ $2, $1, $3 ]; }
  | expr SLASH      expr          { $$ = [ $2, $1, $3 ]; }
  | expr PLUS       expr          { $$ = [ $2, $1, $3 ]; }
  | expr MINUS      expr          { $$ = [ $2, $1, $3 ]; }
  | expr DBLSLASH   expr          { $$ = [ $2, $1, $3 ]; }
  | expr COMPARISON expr          { $$ = [ $2, $1, $3 ]; }
  | expr AND        expr          { $$ = [ $2, $1, $3 ]; }
  | expr OR         expr          { $$ = [ $2, $1, $3 ]; }
  | expr EQV        expr          { $$ = [ $2, $1, $3 ]; }
  | expr NEQV       expr          { $$ = [ $2, $1, $3 ]; }

potential_complex_or_implied_do:
    CONST                         { $$ = [ "%const", @{$1} ]; }
  | CONST COMMA CONST
      { my ($type1, $val1) = @{$1};
        my ($type2, $val2) = @{$3};
        $$ = ["%const", typing::make_complex_type ($type1, $type2),
              [$val1, $val2]];
      }
  | expr_without_const            { $$ = $1; }
  | expr_without_const COMMA do_args
                                  { $$ = [ "%do", $1, @{$3} ]; }
  | CONST COMMA do_args
                                  { $$ = [ "%do", [ "%const", @{$1} ], @{$3} ];
                                  }

array:
    array COMMA array_piece       { $$ = [ @{$1}, $3 ]; }
  | array_piece                   { $$ = [ $1 ]; }

array_piece:
    expr                          { $$ = $1; }
# | implied_do is handled within expr

do_args:
    NAME EQUALS expr COMMA expr   { $$ = [ $1, $3, $5 ]; }
  | NAME EQUALS expr COMMA expr COMMA expr
                                  { $$ = [ $1, $3, $5, $7 ]; }

chain:
    NAME                          { $$ = [ "%var", $1 ]; }
  | chain PERCENT NAME            { $$ = [ $2, $1, $3 ]; }
  | chain LPAREN exprlist RPAREN  { $$ = [ "%call", $1, @{$3} ]; }

exprlist:
                                  { $$ = []; }
  | exprlist_ne                   { $$ = $1; }

exprlist_ne:
    exprlist_ne COMMA argument    { $$ = [ @{$1}, $3 ]; }
  | argument                      { $$ = [ $1 ]; }

argument:
    expr                          { $$ = $1; }
  | colonexpr                     { $$ = $1; }
  | namedargument                 { $$ = $1; }

namedargument:
    NAME EQUALS expr              { $$ = [ "%namedarg", $1, $3 ]; }

colonexpr:
    COLON                         { $$ = [ "%colon", "", "" ]; }
  | expr COLON                    { $$ = [ "%colon", $1, "" ]; }
  | COLON expr                    { $$ = [ "%colon", "", $2 ]; }
  | expr COLON expr               { $$ = [ "%colon", $1, $2 ]; }

%%

sub yylex {
   $expr_parse::left =~ s/^\s*//;
   return 0 if $expr_parse::left eq "";
   my ($ncharsread, $token, $value) = expr_parse::good_yylex ($expr_parse::left);
   # print "yylex: token eof\n" unless $ncharsread;
   return 0 unless $ncharsread;
   # print "yylex: token $token (" . substr ($expr_parse::left, 0, $ncharsread) . ") with value $value\n";
   # print join (";", @$value) . "\n";
   $expr_parse::left = substr ($expr_parse::left, $ncharsread);
   $yylval = $value;
   return $token;
}

# returns (ncharsread, token, value)
sub good_yylex {
   my ($s) = @_;
   my ($c) = substr ($s, 0, 1);

   if ($c eq "") {
      return 0;
   } elsif ($s =~ /^(\d+(?:\.\d*)?|\.\d+)D[+-]?\d+/i) {
      return (length ($&), $CONST, [$typing::double_precision, $&]);
   } elsif ($s =~ /^(\d+E[+-]?\d+|(?:\d+\.\d*|\.\d+)(?:E[+-]?\d+)?)(_\w+)?/i) {
      if (defined $2) {
         return (length ($&), $CONST, [typing::make_type ('real', substr ($2, 1)), $1]);
      } else {
         return (length ($&), $CONST, [$typing::default_type{'real'}, $1]);
      }
   } elsif ($s =~ /^(\d+)(_\w+)?/) {
      if ($2) {
         return (length ($&), $CONST, [typing::make_type ('integer', substr ($2, 1)), $1]);
      } else {
         return (length ($&), $CONST, [$typing::default_type{'integer'}, $1]);
      }
   } elsif ($s =~ /^(\.true\.|\.false\.)(_\w+)?/i) {
      if (defined $2) {
         return (length ($&), $CONST, [typing::make_type ('logical', substr ($2, 1)), $1]);
      } else {
         return (length ($&), $CONST, [$typing::default_type{'logical'}, $1]);
      }
   } elsif ($s =~ /^'(\d+)'(_\w+)?/) {
      # Interior of string is digits because it has been grabbed already.
      my ($str) = stmts::get_string ($1);
      if (defined $2) {
         return (length ($&), $CONST, [typing::make_character_type (substr ($2, 1), length ($str)), $str]);
      } else {
         return (length ($&), $CONST, [typing::make_character_type ($typing::default_character_kind, length ($str)), $str]);
      }
   } elsif ($s =~ /^\w+/) {
      return (length ($&), $NAME, $&);
   } else {
      switch: {
         $s =~ /^==/      && return (2, $COMPARISON, "==");
         $s =~ /^<=/      && return (2, $COMPARISON, "<=");
         $s =~ /^>=/      && return (2, $COMPARISON, ">=");
         $s =~ /^</       && return (1, $COMPARISON, "<");
         $s =~ /^>/       && return (1, $COMPARISON, ">");
         $s =~ /^\/=/     && return (2, $COMPARISON, "/=");
         $s =~ /^=/       && return (1, $EQUALS, "=");
         $s =~ /^\.eq\./i && return (4, $COMPARISON, "==");
         $s =~ /^\.le\./i && return (4, $COMPARISON, "<=");
         $s =~ /^\.ge\./i && return (4, $COMPARISON, ">=");
         $s =~ /^\.lt\./i && return (4, $COMPARISON, "<");
         $s =~ /^\.gt\./i && return (4, $COMPARISON, ">");
         $s =~ /^\.ne\./i && return (4, $COMPARISON, "/=");
         $s =~ /^\.neqv\./i && return (6, $NEQV, ".neqv.");
         $s =~ /^\.eqv\./i && return (5, $EQV, ".eqv.");
         $s =~ /^\.and\./i && return (5, $AND, ".and.");
         $s =~ /^\.or\./i && return (4, $OR, ".or.");
         $s =~ /^\.not\./i && return (5, $NOT, ".not.");
         $s =~ /^\*\*/    && return (2, $DBLASTERIK, "**");
         $s =~ /^\/\//    && return (2, $DBLSLASH, "//");
         $s =~ /^\(\//    && return (2, $LARRAY, "(/");
         $s =~ /^\/\)/    && return (2, $RARRAY, "/)");
         $c eq ","        && return (1, $COMMA, ",");
         $c eq "+"        && return (1, $PLUS, "+");
         $c eq "-"        && return (1, $MINUS, "-");
         $c eq "*"        && return (1, $ASTERIK, "*");
         $c eq "/"        && return (1, $SLASH, "/");
         $c eq "("        && return (1, $LPAREN, "(");
         $c eq ")"        && return (1, $RPAREN, ")");
         $c eq "%"        && return (1, $PERCENT, "%");
         $c eq ":"        && return (1, $COLON, ":");
      }
      die "Lexer failed on `$s'";
   }
}

#####
# Takes a string that consists entirely of an expression, and returns a
# reference to the parse tree it defines.
#####
sub parse_expr {
  my ($s) = @_;
  # print "parsing string: $s.\n";
  $expr_parse::left = $expr_parse::line = $s;
  die "Expression `$expr_parse::line' has trailing garbage `$1$expr_parse::left'"
    if yyparse () =~ /^s(.*)$/;
  return $yyval;
}

#####
# Takes a string that consists partly of an expression.  (The first part
# is an expression.)  Returns (parse tree ref, rest string, separator string).
#####
sub parse_part_as_expr {
  my ($s) = @_;
  # print "parsing part of string: $s.\n";
  $expr_parse::left = $expr_parse::line = $s;
  if (yyparse () =~ /^s(.*)$/) {
    return ($yyval, $expr_parse::left, $1);
  } else {
    return ($yyval);
  }
}

sub yyerror {
  my ($s) = @_;
  die "yyerror: $s during parsing of F90 code `$expr_parse::line'";
}

1;
