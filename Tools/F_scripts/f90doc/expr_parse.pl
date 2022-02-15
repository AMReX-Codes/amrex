$yysccsid = "@(#)yaccpar 1.8 (Berkeley) 01/20/91 (Perl 2.0 12/31/92)";
#define YYBYACC 1
#line 2 "expr_parse.y"
package expr_parse;

;# On failure, print out this as the line we were working on.
$expr_parse::line = "";

;# Portion of line left to parse
$expr_parse::left = "";
#line 12 "y.tab.pl"
$COMMA=257;
$LPAREN=258;
$RPAREN=259;
$NOT=260;
$OR=261;
$AND=262;
$EQV=263;
$NEQV=264;
$COMPARISON=265;
$DBLSLASH=266;
$PERCENT=267;
$PLUS=268;
$MINUS=269;
$UPLUS=270;
$UMINUS=271;
$ASTERIK=272;
$SLASH=273;
$DBLASTERIK=274;
$CONST=275;
$NAME=276;
$COLON=277;
$LARRAY=278;
$RARRAY=279;
$EQUALS=280;
$YYERRCODE=256;
@yylhs = (                                               -1,
    0,    0,    1,    1,    2,    2,    2,    2,    2,    2,
    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
    2,    5,    5,    5,    5,    5,    4,    4,    7,    6,
    6,    3,    3,    3,    8,    8,    9,    9,   10,   10,
   10,   12,   11,   11,   11,   11,
);
@yylen = (                                                2,
    1,    2,    1,    1,    1,    3,    2,    2,    2,    3,
    3,    3,    3,    3,    3,    3,    3,    3,    3,    3,
    3,    1,    3,    1,    3,    3,    3,    1,    1,    5,
    7,    1,    3,    4,    0,    1,    3,    1,    1,    1,
    1,    3,    1,    2,    2,    3,
);
@yydefred = (                                             0,
    0,    0,    0,    0,    3,   32,    0,    0,    0,    4,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
   28,    2,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,   10,    0,    6,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,   38,   40,   41,   33,
   23,    0,   26,   25,   27,    0,    0,    0,   34,    0,
    0,    0,    0,   37,    0,    0,    0,    0,    0,
);
@yydgoto = (                                              8,
   19,   10,   11,   20,   15,   63,   21,   55,   56,   57,
   58,   59,
);
@yysindex = (                                          -212,
 -157, -212, -212, -212,    0,    0, -212,    0, -137,    0,
 -246, -241,  -29, -234, -235,  -19, -223, -223,  -29, -257,
    0,    0, -212, -212, -212, -212, -212, -212, -212, -212,
 -212, -212, -212, -216, -229, -267, -222,    0, -212,    0,
 -255,  -19,  227,  227,  236, -164, -223, -223, -233, -233,
 -233, -205, -212,  -76, -174, -162,    0,    0,    0,    0,
    0, -180,    0,    0,    0, -212,  -29, -212,    0, -216,
 -212,  -29,  -29,    0, -118, -212,  -95, -212,  -29,
);
@yyrindex = (                                             0,
    0,    0,    0,    0,    0,    0,    0,    0,  106,    0,
    1,  -59,    0,  -43,    0,  163,   77,   96, -242,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0, -152,    0,    0,    0,    0,    0,    0,
  191,  172,  199,  208,  182,  153,  115,  134,   20,   39,
   58, -175, -219, -214,    0, -146,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0, -192, -188,    0,    0,
    0, -183, -178,    0,    0,    0, -145,    0, -143,
);
@yygindex = (                                             0,
    2,  116,    0,    0,    0,   85,   84,    0,    0,   60,
    0,    0,
);
$YYTABLESIZE=510;
@yytable = (                                             39,
    5,    9,   13,   16,   17,   18,   24,   61,   62,   27,
   28,   34,   29,   30,   29,   36,   31,   32,   33,   12,
   35,   40,   37,   38,   41,   42,   43,   44,   45,   46,
   47,   48,   49,   50,   51,   54,   29,   43,   13,   43,
   33,    1,   39,    2,   39,    1,   60,    2,   31,   32,
   33,    3,    4,   62,   67,    3,    4,   11,    5,   52,
   53,    7,    5,    6,   45,    7,   45,   72,   44,   73,
   44,   54,   75,   42,   66,   42,    7,   77,   46,   79,
   46,   32,   32,   32,   69,   32,   32,   32,   32,   32,
   32,   32,   32,   32,   70,    8,   32,   32,   32,   71,
    1,   32,    2,   29,   30,    1,   35,   31,   32,   33,
    3,    4,   36,   30,   14,   31,   14,   12,    6,   22,
    7,   64,   65,   23,   24,   25,   26,   27,   28,   74,
   29,   30,    0,   15,   31,   32,   33,    0,   76,    0,
    0,    0,   23,   24,   25,   26,   27,   28,    0,   29,
   30,    0,   16,   31,   32,   33,    0,    0,    0,    0,
    0,   78,    9,    0,    0,   23,   24,   25,   26,   27,
   28,   18,   29,   30,    0,    0,   31,   32,   33,    0,
    0,   17,    0,    0,   23,   24,   25,   26,   27,   28,
   19,   29,   30,    0,    0,   31,   32,   33,   20,   22,
   68,    3,    3,    3,    3,    3,    3,   21,    3,    3,
    0,    0,    3,    3,    3,   24,    0,    4,    4,    4,
    4,    4,    4,    0,    4,    4,    0,    0,    4,    4,
    4,   23,   24,   25,   26,   27,   28,    0,   29,   30,
    0,    0,   31,   32,   33,   27,   28,    0,   29,   30,
    0,    0,   31,   32,   33,    0,    0,    5,    0,    5,
    0,    5,    5,    5,    5,    5,    5,    0,    5,    5,
    0,    0,    5,    5,    5,    0,   12,    5,   12,    5,
   12,   12,   12,   12,   12,   12,    0,   12,   12,    0,
    0,   12,   12,    0,    0,   13,   12,   13,   12,   13,
   13,   13,   13,   13,   13,    0,   13,   13,    0,    0,
   13,   13,    0,    0,   11,   13,   11,   13,   11,   11,
   11,   11,   11,   11,    0,   11,   11,    0,    0,   11,
   11,    0,    0,    7,   11,    7,   11,    7,    7,    7,
    7,    7,    7,    0,    7,    7,    0,    0,    0,    0,
    0,    0,    8,    7,    8,    7,    8,    8,    8,    8,
    8,    8,    0,    8,    8,    0,    0,    0,    0,    0,
    0,   14,    8,   14,    8,   14,   14,   14,   14,   14,
   14,    0,   14,   14,    0,    0,    0,    0,    0,    0,
   15,   14,   15,   14,   15,   15,   15,   15,   15,   15,
    0,   15,   15,    0,    0,    0,    0,    0,    0,   16,
   15,   16,   15,   16,   16,   16,   16,   16,   16,    9,
    0,    9,    0,    9,    9,    9,    9,    0,   18,   16,
   18,   16,   18,   18,   18,   18,    0,    0,   17,    9,
   17,    9,   17,   17,   17,   17,    0,   19,   18,   19,
   18,   19,    0,   19,   19,   20,    0,   20,   17,    0,
   17,   20,   20,    0,   21,    0,   21,   19,    0,   19,
   21,   21,    0,    0,    0,   20,    0,   20,    0,    0,
    0,    0,    0,    0,   21,    0,   21,   23,   24,    0,
    0,   27,   28,    0,   29,   30,    0,    0,   31,   32,
   33,   28,    0,   29,   30,    0,    0,   31,   32,   33,
);
@yycheck = (                                            257,
    0,    0,    1,    2,    3,    4,  262,  275,  276,  265,
  266,  258,  268,  269,  257,  257,  272,  273,  274,    0,
  267,  279,  257,  259,   23,   24,   25,   26,   27,   28,
   29,   30,   31,   32,   33,   34,  279,  257,    0,  259,
  274,  258,  257,  260,  259,  258,  276,  260,  272,  273,
  274,  268,  269,  276,   53,  268,  269,    0,  275,  276,
  277,  278,  275,  276,  257,  278,  259,   66,  257,   68,
  259,   70,   71,  257,  280,  259,    0,   76,  257,   78,
  259,  257,  258,  259,  259,  261,  262,  263,  264,  265,
  266,  267,  268,  269,  257,    0,  272,  273,  274,  280,
  258,  277,  260,  268,  269,    0,  259,  272,  273,  274,
  268,  269,  259,  259,    0,  259,    1,  275,  276,  257,
  278,   37,   39,  261,  262,  263,  264,  265,  266,   70,
  268,  269,   -1,    0,  272,  273,  274,   -1,  257,   -1,
   -1,   -1,  261,  262,  263,  264,  265,  266,   -1,  268,
  269,   -1,    0,  272,  273,  274,   -1,   -1,   -1,   -1,
   -1,  257,    0,   -1,   -1,  261,  262,  263,  264,  265,
  266,    0,  268,  269,   -1,   -1,  272,  273,  274,   -1,
   -1,    0,   -1,   -1,  261,  262,  263,  264,  265,  266,
    0,  268,  269,   -1,   -1,  272,  273,  274,    0,  259,
  277,  261,  262,  263,  264,  265,  266,    0,  268,  269,
   -1,   -1,  272,  273,  274,  259,   -1,  261,  262,  263,
  264,  265,  266,   -1,  268,  269,   -1,   -1,  272,  273,
  274,  261,  262,  263,  264,  265,  266,   -1,  268,  269,
   -1,   -1,  272,  273,  274,  265,  266,   -1,  268,  269,
   -1,   -1,  272,  273,  274,   -1,   -1,  257,   -1,  259,
   -1,  261,  262,  263,  264,  265,  266,   -1,  268,  269,
   -1,   -1,  272,  273,  274,   -1,  257,  277,  259,  279,
  261,  262,  263,  264,  265,  266,   -1,  268,  269,   -1,
   -1,  272,  273,   -1,   -1,  257,  277,  259,  279,  261,
  262,  263,  264,  265,  266,   -1,  268,  269,   -1,   -1,
  272,  273,   -1,   -1,  257,  277,  259,  279,  261,  262,
  263,  264,  265,  266,   -1,  268,  269,   -1,   -1,  272,
  273,   -1,   -1,  257,  277,  259,  279,  261,  262,  263,
  264,  265,  266,   -1,  268,  269,   -1,   -1,   -1,   -1,
   -1,   -1,  257,  277,  259,  279,  261,  262,  263,  264,
  265,  266,   -1,  268,  269,   -1,   -1,   -1,   -1,   -1,
   -1,  257,  277,  259,  279,  261,  262,  263,  264,  265,
  266,   -1,  268,  269,   -1,   -1,   -1,   -1,   -1,   -1,
  257,  277,  259,  279,  261,  262,  263,  264,  265,  266,
   -1,  268,  269,   -1,   -1,   -1,   -1,   -1,   -1,  257,
  277,  259,  279,  261,  262,  263,  264,  265,  266,  257,
   -1,  259,   -1,  261,  262,  263,  264,   -1,  257,  277,
  259,  279,  261,  262,  263,  264,   -1,   -1,  257,  277,
  259,  279,  261,  262,  263,  264,   -1,  257,  277,  259,
  279,  261,   -1,  263,  264,  257,   -1,  259,  277,   -1,
  279,  263,  264,   -1,  257,   -1,  259,  277,   -1,  279,
  263,  264,   -1,   -1,   -1,  277,   -1,  279,   -1,   -1,
   -1,   -1,   -1,   -1,  277,   -1,  279,  261,  262,   -1,
   -1,  265,  266,   -1,  268,  269,   -1,   -1,  272,  273,
  274,  266,   -1,  268,  269,   -1,   -1,  272,  273,  274,
);
$YYFINAL=8;
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
$YYMAXTOKEN=280;
#if YYDEBUG
@yyname = (
"end-of-file",'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','',
'','','','','','','','','','','','','','','','','','','','','','','',"COMMA","LPAREN","RPAREN","NOT",
"OR","AND","EQV","NEQV","COMPARISON","DBLSLASH","PERCENT","PLUS","MINUS",
"UPLUS","UMINUS","ASTERIK","SLASH","DBLASTERIK","CONST","NAME","COLON","LARRAY",
"RARRAY","EQUALS",
);
@yyrule = (
"\$accept : expr_with_abort",
"expr_with_abort : expr",
"expr_with_abort : expr COMMA",
"expr : CONST",
"expr : expr_without_const",
"expr_without_const : chain",
"expr_without_const : LARRAY array RARRAY",
"expr_without_const : PLUS expr",
"expr_without_const : MINUS expr",
"expr_without_const : NOT expr",
"expr_without_const : LPAREN potential_complex_or_implied_do RPAREN",
"expr_without_const : expr DBLASTERIK expr",
"expr_without_const : expr ASTERIK expr",
"expr_without_const : expr SLASH expr",
"expr_without_const : expr PLUS expr",
"expr_without_const : expr MINUS expr",
"expr_without_const : expr DBLSLASH expr",
"expr_without_const : expr COMPARISON expr",
"expr_without_const : expr AND expr",
"expr_without_const : expr OR expr",
"expr_without_const : expr EQV expr",
"expr_without_const : expr NEQV expr",
"potential_complex_or_implied_do : CONST",
"potential_complex_or_implied_do : CONST COMMA CONST",
"potential_complex_or_implied_do : expr_without_const",
"potential_complex_or_implied_do : expr_without_const COMMA do_args",
"potential_complex_or_implied_do : CONST COMMA do_args",
"array : array COMMA array_piece",
"array : array_piece",
"array_piece : expr",
"do_args : NAME EQUALS expr COMMA expr",
"do_args : NAME EQUALS expr COMMA expr COMMA expr",
"chain : NAME",
"chain : chain PERCENT NAME",
"chain : chain LPAREN exprlist RPAREN",
"exprlist :",
"exprlist : exprlist_ne",
"exprlist_ne : exprlist_ne COMMA argument",
"exprlist_ne : argument",
"argument : expr",
"argument : colonexpr",
"argument : namedargument",
"namedargument : NAME EQUALS expr",
"colonexpr : COLON",
"colonexpr : expr COLON",
"colonexpr : COLON expr",
"colonexpr : expr COLON expr",
);
#endif
sub yyclearin { $yychar = -1; }
sub yyerrok { $yyerrflag = 0; }
$YYSTACKSIZE = $YYSTACKSIZE || $YYMAXDEPTH || 500;
$YYMAXDEPTH = $YYMAXDEPTH || $YYSTACKSIZE || 500;
$yyss[$YYSTACKSIZE] = 0;
$yyvs[$YYSTACKSIZE] = 0;
sub YYERROR { ++$yynerrs; &yy_err_recover; }
sub yy_err_recover
{
  if ($yyerrflag < 3)
  {
    $yyerrflag = 3;
    while (1)
    {
      if (($yyn = $yysindex[$yyss[$yyssp]]) && 
          ($yyn += $YYERRCODE) >= 0 && 
          $yycheck[$yyn] == $YYERRCODE)
      {
#if YYDEBUG
       print "yydebug: state $yyss[$yyssp], error recovery shifting",
             " to state $yytable[$yyn]\n" if $yydebug;
#endif
        $yyss[++$yyssp] = $yystate = $yytable[$yyn];
        $yyvs[++$yyvsp] = $yylval;
        next yyloop;
      }
      else
      {
#if YYDEBUG
        print "yydebug: error recovery discarding state ",
              $yyss[$yyssp], "\n"  if $yydebug;
#endif
        return(1) if $yyssp <= 0;
        --$yyssp;
        --$yyvsp;
      }
    }
  }
  else
  {
    return (1) if $yychar == 0;
#if YYDEBUG
    if ($yydebug)
    {
      $yys = '';
      if ($yychar <= $YYMAXTOKEN) { $yys = $yyname[$yychar]; }
      if (!$yys) { $yys = 'illegal-symbol'; }
      print "yydebug: state $yystate, error recovery discards ",
            "token $yychar ($yys)\n";
    }
#endif
    $yychar = -1;
    next yyloop;
  }
0;
} # yy_err_recover

sub yyparse
{
#ifdef YYDEBUG
  if ($yys = $ENV{'YYDEBUG'})
  {
    $yydebug = int($1) if $yys =~ /^(\d)/;
  }
#endif

  $yynerrs = 0;
  $yyerrflag = 0;
  $yychar = (-1);

  $yyssp = 0;
  $yyvsp = 0;
  $yyss[$yyssp] = $yystate = 0;

yyloop: while(1)
  {
    yyreduce: {
      last yyreduce if ($yyn = $yydefred[$yystate]);
      if ($yychar < 0)
      {
        if (($yychar = &yylex) < 0) { $yychar = 0; }
#if YYDEBUG
        if ($yydebug)
        {
          $yys = '';
          if ($yychar <= $#yyname) { $yys = $yyname[$yychar]; }
          if (!$yys) { $yys = 'illegal-symbol'; };
          print "yydebug: state $yystate, reading $yychar ($yys)\n";
        }
#endif
      }
      if (($yyn = $yysindex[$yystate]) && ($yyn += $yychar) >= 0 &&
              $yycheck[$yyn] == $yychar)
      {
#if YYDEBUG
        print "yydebug: state $yystate, shifting to state ",
              $yytable[$yyn], "\n"  if $yydebug;
#endif
        $yyss[++$yyssp] = $yystate = $yytable[$yyn];
        $yyvs[++$yyvsp] = $yylval;
        $yychar = (-1);
        --$yyerrflag if $yyerrflag > 0;
        next yyloop;
      }
      if (($yyn = $yyrindex[$yystate]) && ($yyn += $yychar) >= 0 &&
            $yycheck[$yyn] == $yychar)
      {
        $yyn = $yytable[$yyn];
        last yyreduce;
      }
      if (! $yyerrflag) {
        &yyerror('syntax error');
        ++$yynerrs;
      }
      return(1) if &yy_err_recover;
    } # yyreduce
#if YYDEBUG
    print "yydebug: state $yystate, reducing by rule ",
          "$yyn ($yyrule[$yyn])\n"  if $yydebug;
#endif
    $yym = $yylen[$yyn];
    $yyval = $yyvs[$yyvsp+1-$yym];
    switch:
    {
if ($yyn == 1) {
#line 29 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; return 1; 
last switch;
} }
if ($yyn == 2) {
#line 30 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-1]; return "s,"; 
last switch;
} }
if ($yyn == 3) {
#line 33 "expr_parse.y"
{ $yyval = [ "%const", @{$yyvs[$yyvsp-0]} ]; 
last switch;
} }
if ($yyn == 4) {
#line 34 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; 
last switch;
} }
if ($yyn == 5) {
#line 37 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; 
last switch;
} }
if ($yyn == 6) {
#line 38 "expr_parse.y"
{ $yyval = [ "%array", @{$yyvs[$yyvsp-1]} ]; 
last switch;
} }
if ($yyn == 7) {
#line 39 "expr_parse.y"
{ $yyval = [ "u+", $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 8) {
#line 40 "expr_parse.y"
{ $yyval = [ "u-", $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 9) {
#line 41 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 10) {
#line 43 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-1]; 
last switch;
} }
if ($yyn == 11) {
#line 44 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 12) {
#line 45 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 13) {
#line 46 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 14) {
#line 47 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 15) {
#line 48 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 16) {
#line 49 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 17) {
#line 50 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 18) {
#line 51 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 19) {
#line 52 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 20) {
#line 53 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 21) {
#line 54 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 22) {
#line 57 "expr_parse.y"
{ $yyval = [ "%const", @{$yyvs[$yyvsp-0]} ]; 
last switch;
} }
if ($yyn == 23) {
#line 59 "expr_parse.y"
{ my ($type1, $val1) = @{$yyvs[$yyvsp-2]};
        my ($type2, $val2) = @{$yyvs[$yyvsp-0]};
        $yyval = ["%const", typing::make_complex_type ($type1, $type2),
              [$val1, $val2]];
      
last switch;
} }
if ($yyn == 24) {
#line 64 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; 
last switch;
} }
if ($yyn == 25) {
#line 66 "expr_parse.y"
{ $yyval = [ "%do", $yyvs[$yyvsp-2], @{$yyvs[$yyvsp-0]} ]; 
last switch;
} }
if ($yyn == 26) {
#line 68 "expr_parse.y"
{ $yyval = [ "%do", [ "%const", @{$yyvs[$yyvsp-2]} ], @{$yyvs[$yyvsp-0]} ];
                                  
last switch;
} }
if ($yyn == 27) {
#line 72 "expr_parse.y"
{ $yyval = [ @{$yyvs[$yyvsp-2]}, $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 28) {
#line 73 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 29) {
#line 76 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; 
last switch;
} }
if ($yyn == 30) {
#line 80 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-4], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 31) {
#line 82 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-6], $yyvs[$yyvsp-4], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 32) {
#line 85 "expr_parse.y"
{ $yyval = [ "%var", $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 33) {
#line 86 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-1], $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 34) {
#line 87 "expr_parse.y"
{ $yyval = [ "%call", $yyvs[$yyvsp-3], @{$yyvs[$yyvsp-1]} ]; 
last switch;
} }
if ($yyn == 35) {
#line 90 "expr_parse.y"
{ $yyval = []; 
last switch;
} }
if ($yyn == 36) {
#line 91 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; 
last switch;
} }
if ($yyn == 37) {
#line 94 "expr_parse.y"
{ $yyval = [ @{$yyvs[$yyvsp-2]}, $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 38) {
#line 95 "expr_parse.y"
{ $yyval = [ $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 39) {
#line 98 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; 
last switch;
} }
if ($yyn == 40) {
#line 99 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; 
last switch;
} }
if ($yyn == 41) {
#line 100 "expr_parse.y"
{ $yyval = $yyvs[$yyvsp-0]; 
last switch;
} }
if ($yyn == 42) {
#line 103 "expr_parse.y"
{ $yyval = [ "%namedarg", $yyvs[$yyvsp-2], $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 43) {
#line 106 "expr_parse.y"
{ $yyval = [ "%colon", "", "" ]; 
last switch;
} }
if ($yyn == 44) {
#line 107 "expr_parse.y"
{ $yyval = [ "%colon", $yyvs[$yyvsp-1], "" ]; 
last switch;
} }
if ($yyn == 45) {
#line 108 "expr_parse.y"
{ $yyval = [ "%colon", "", $yyvs[$yyvsp-0] ]; 
last switch;
} }
if ($yyn == 46) {
#line 109 "expr_parse.y"
{ $yyval = [ "%colon", $yyvs[$yyvsp-2], $yyvs[$yyvsp-1] ]; 
last switch;
} }
#line 624 "y.tab.pl"
    } # switch
    $yyssp -= $yym;
    $yystate = $yyss[$yyssp];
    $yyvsp -= $yym;
    $yym = $yylhs[$yyn];
    if ($yystate == 0 && $yym == 0)
    {
#if YYDEBUG
      print "yydebug: after reduction, shifting from state 0 ",
            "to state $YYFINAL\n" if $yydebug;
#endif
      $yystate = $YYFINAL;
      $yyss[++$yyssp] = $YYFINAL;
      $yyvs[++$yyvsp] = $yyval;
      if ($yychar < 0)
      {
        if (($yychar = &yylex) < 0) { $yychar = 0; }
#if YYDEBUG
        if ($yydebug)
        {
          $yys = '';
          if ($yychar <= $#yyname) { $yys = $yyname[$yychar]; }
          if (!$yys) { $yys = 'illegal-symbol'; }
          print "yydebug: state $YYFINAL, reading $yychar ($yys)\n";
        }
#endif
      }
      return(0) if $yychar == 0;
      next yyloop;
    }
    if (($yyn = $yygindex[$yym]) && ($yyn += $yystate) >= 0 &&
        $yyn <= $#yycheck && $yycheck[$yyn] == $yystate)
    {
        $yystate = $yytable[$yyn];
    } else {
        $yystate = $yydgoto[$yym];
    }
#if YYDEBUG
    print "yydebug: after reduction, shifting from state ",
        "$yyss[$yyssp] to state $yystate\n" if $yydebug;
#endif
    $yyss[++$yyssp] = $yystate;
    $yyvs[++$yyvsp] = $yyval;
  } # yyloop
} # yyparse
#line 112 "expr_parse.y"

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
#line 794 "y.tab.pl"
