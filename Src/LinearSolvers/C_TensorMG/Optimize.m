(* :Name: Optimize` *)

(* :Title: Expression Optimization. *)

(* :Author: Mark Sofroniou *)

(* :Summary:
 This package adds the procedure Optimize for performing common
 sub-expression optimization on Mathematica expressions.
 An optimized Module and an optimized compiled Module is also possible.
 The function Horner factors uni-variate and multi-variate polynomials
 in efficient computational form.
 The function Cost measures the computational expense of numerical
 evaluation of an expression.
 The additional package Format.m (MathSource) enables the production of
 efficient C and Fortran code. *)

(* :Context: Optimize` *)

(* :Package Version: 1.2 *)

(* :Copyright: Copyright 1993-4,  Mark Sofroniou.
 Permission is hereby granted to modify and/or make copies of
 this file for any purpose other than direct profit, or as part
 of a commercial product, provided this copyright notice is left
 intact. Sale, other than for the cost of media, is prohibited.

 Permission is hereby granted to reproduce part or all of
 this file, provided that the source is acknowledged. *)

(* :History:
 Modified Optimize syntax and improved performance, August 1994.
 Significantly revised and publicly released May, 1994.
 Original Version by Mark Sofroniou, September, 1993. *)

(* :Keywords:
 Assign, C, CAssign, Common, Compile, Cost, CForm, FORTRAN,
 FortranAssign, FortranForm, Horner, Optimize, Optimization,
 Polynomial, Sub-Expression, Syntactic. *)

(* :Source:
 Mark Sofroniou, Ph.D. Thesis, Loughborough University, Loughborough,
 Leicestershire LE11 3TU, England. *)

(* :Mathematica Version: 2.2 *)

(* :Limitations:
 This package enables syntactic optimization in linear time.
 Optimization is used in the sense of improving the arithmetic
 operation count and compactness of the code, rather than the
 best possible.
 Expressions containing non-binary arithmetic operators (Plus
 and Times) are optimized by matching only entire sub-expressions
 possibly after extracting numeric coefficients.
 These operations rarely dominate the computational cost.
 Options enable control over the optimization process.
 The issue of numerical stability has not been addressed. *)

BeginPackage["Optimize`"]

Cost::usage = "Cost[expr,options] returns a list of Mathematica operators
present in expr, providing a count of the basic arithmetic operations and
function calls. Cost has attribute HoldAll so that arguments are maintained
in unevaluated form. Some examples of the behaviour (using options):\n
1+2*3 is counted as an addition and a multiplication.
a+b+c is counted as two additions and a*b*c as two multiplications.\n
Power[E,x] is counted as the exponential function Exp[x].
Power[x,-1] is counted as a division.\n
Power[x,y] is calculated as y-1 multiplications for integer y."

Horner::usage = "Horner[poly] puts the polynomial poly in Horner or nested form
with respect to the default variables Variables[poly]. This is an efficient form
for numerical evaluation. Horner[poly,vars] specifies the variable ordering
explicitly as the List vars. Multinomial conversion is applied recursively.
Horner factorisation of rational polynomials is possible as Horner[p1/p2] or
Horner[p1/p2,v1,v2]. Pre-conditioning of the coefficients is sometimes necessary
to attain numerical stability and more efficient methods exist in this case.
Horner's rule is optimal if no pre-conditioning is assumed, requiring n
multiplications and n additions to evaluate a polynomial of degree n."

Optimize::usage = "Optimize[expr,opts]\n
Optimize transforms expr into a sequence of optimization statements and an
optimized expression. Optimization is performed in linear time (O(n) operations)
providing an efficient means of reducing the arithmetic operation count - not the
best possible. The optimization performed is mainly syntactic (only exact common
sub-expressions are matched) with a few heuristics. Options opts control the
type of optimizations performed and the format of the output."


(* Options: *)

CostDivide::usage = "CostDivide is an option of Cost specifying whether
to consider negative exponentiations as divisions. The setting All
corresponds to the standard evaluation procedure and the compiler.
The setting Share counts only one division in an expression such as
(x^-2)*(y^-2) in correspondence with OuputForm and C and FORTRAN code
translation. CostDivide may evaluate to All, Share or False. 
See also CostPower"

CostExp::usage = "CostExp is an option of Cost specifying whether to
consider Power[E,x] as a call to the exponential function Exp[x].
CostExp may evaluate to True or False."

CostNull::usage = "CostNull is an option of Cost specifying a
list of symbols to be omitted for the purposes of costing. This is
useful, for example, for removing named arrays from consideration
(which have the same syntax as functions).
CostNull may evaluate to a (possibly empty) list of symbols."

CostPower::usage="CostPower is an option of Cost specifying whether to
consider Power[x,y] as calculated using two function calls, namely
Exp[y*Log[x]]. The setting Integer assumes in addition that integer
powers are calculated using y-1 multiplications. This is particularly
useful when costing expressions optimized using the option
OptimizePower->Binary. CostPower may evaluate to Integer, True or False.
See also CostDivide"

Binary::usage = "Binary is an argument of the option OptimizePower
specifying whether to perform repeated binary factoring of exponents."

OptimizeCoefficients::usage = "OptimizeCoefficients is an option of
Optimize specifying whether to extract numerical coefficients in
expressions with head Plus and Times. OptimizeCoefficients may
evaluate to True or False."

OptimizeFunction::usage = "OptimizeFunction is an option of Optimize
specifying whether to optimize common sub-expressions which are
neither atoms (?AtomQ) nor expressions with head Plus, Power,
or Times. OptimizeFunction may evaluate to True or False."

OptimizeNull::usage = "OptimizeNull is an option of Optimize specifying
a list of symbols to be disregarded during the optimization process.
This is useful, for example, for specifying a named array (which has the
same syntax as a function). OptimizeNull may evaluate to a (possibly empty)
list of symbols."

OptimizePlus::usage = "OptimizePlus is an option of Optimize specifying
whether to extract common sub-expressions with head Plus. Only exact
sub-expressions are matched. Numerical coefficients may be extracted
using OptimizeCoefficients. OptimizePlus may evaluate to True or False."

OptimizePower::usage = "OptimizePower is an option of Optimize
specifying whether to optimize sub-expressions with head Power.
Integer and rational powers can be factored efficiently by repeated
binary exponentiation. OptimizePower may evaluate to True, False,
or Binary."

OptimizeProcedure::usage = "OptimizeProcedure is an option of Optimize
specifying the form of the output returned. The setting True returns a
Module wrapped in Hold. Function, Compile and SetDelayed all know
about Optimize and make use of the option OptimizeProcedure (automatically
removing Hold). Applying ReleaseHold recovers the original (unoptimized)
expression. With the setting False a list of transformation
rules is returned together with the optimized expression. A definition
{optrules,optexpr} = Optimize[expr] enables the original
(unoptimized) expression to be recovered as optexpr //. Dispatch[optrules]."

OptimizeTimes::usage = "OptimizeTimes is an option of Optimize specifying
whether to extract common sub-expressions with head Times. Only exact
sub-expressions are matched. Numerical coefficients may be extracted
using OptimizeCoefficients. OptimizeTimes may evaluate to True or False."

OptimizeVariable::usage = "OptimizeVariable is an option of Optimize
specifying a pair {symb,form} where symb is the optimization variable
to introduce and form is the type of format to use.
The format can be either Sequence (o1,o2,...) or Array (o[1],o[2],...).
Consecutively numbered variables are returned."

(* Default optimization symbol. *)

SetAttributes[o,NHoldAll];

o::usage="o is the default optimization symbol introduced by
Optimize when performing common sub-expression optimization."

Unprotect[Binary, Cost, CostDivide, CostExp, CostNull, CostPower, Horner, o,
Optimize, OptimizeCoefficients, OptimizeFunction, OptimizeNull, OptimizePlus,
OptimizePower, OptimizeProcedure, OptimizeTimes, OptimizeVariable];

Begin["`Private`"]

(* Function to check the data types of Cost options with error
 messages. *)

Options[Cost] = {CostDivide->All,CostExp->True,
CostNull->{Block,Module,With,CompoundExpression,Hold},
CostPower->Integer};

Cost::args = "The `1` did not evaluate to `2`.";

costerrmssgs = {{"option CostDivide","All, Share or False"},
{"option CostExp","True or False"},
{"option CostNull","a (possibly empty) list of symbols"},
{"option CostAsAtom","Integer, True or False"}};

CostOptTest[opts___]:=
  Module[{costdiv,costexp,costnull,costpower,types,defaults=Options[Cost]},

    optlist = {costdiv,costexp,costnull,costpower} =
      Map[First,defaults] /. {opts} /. defaults;

    types = {MatchQ[costdiv,All|Share|False],
      MatchQ[costexp,True|False],
      MatchQ[costnull,{___Symbol}],
      MatchQ[costpower,Integer|True|False]};

    Check[
      MapThread[
        If[#1,#1,Message[Cost::args,Apply[Sequence,#2]]]&,
        {types,costerrmssgs}
      ]; optlist,      (* Return list of option values. *)
      $Failed,         (* Option of wrong type. *)
      Cost::args       (* Check only for these messages. *)
    ]
  ];


SetAttributes[{Cost,MainCost},HoldAll];

Cost[expr_,opts___?OptionQ]:=
  Module[{optvals},
    optvals /;
      And[
        (optvals = CostOptTest[opts])=!=$Failed,
        optvals = MainCost[Unevaluated[expr],Evaluate[optvals]]; True
      ]
  ];



(* Function for counting cost of basic operations. *)

MainCost[expr_,{costdiv_,costexp_,costnull_,costpower_}]:=
  Block[{CostFunction},

(* Prevent evaluation during costing. *)

    SetAttributes[{CostFunction},HoldAll];

(* Rules to cost expressions. *)

    If[costexp, CostFunction[Power[E,_]]:= Exp ];

(* Integer powers as multiplications and divisions. *)

    If[costpower===Integer,
      If[MatchQ[costdiv,All|Share],
        CostFunction[Power[_,y_Integer?Positive]]:= (y-1) Times;
        CostFunction[Power[_,y_Integer?Negative]]:= Divide + (-y-1) Times,
        CostFunction[Power[_,y_Integer]]:= (Abs[y]-1) Times  (* No divisions. *)
      ]
    ];

(* x^y calculated as Exp[y Log[x]]. *)

    If[MatchQ[costpower,Integer|True],
      If[MatchQ[costdiv,All|Share],
        CostFunction[Power[_,_?Negative]]:= Divide+Exp+Log ];
      CostFunction[_Power]:= Exp+Log
    ];

(* Add back in occurrances of shared divisions. *)

    If[costdiv===Share,

(* Pure reciprocal case. *)

      Literal[CostFunction[e:Times[Power[_,_?Negative]..]]]:=
        (1-Length[Unevaluated[e]]) Divide + (Length[Unevaluated[e]]-1) Times;

(* Mixed case: (a*b..)/(c*d..). Similar to above, but don't count first
 division as a multiplication. Hence: Length[num]-1+Length[denom]-1 = Length[e]-2. *)

      Literal[CostFunction[e:(Times[___,Power[_,_?Negative],___]..)]]:=
        (1-Count[Unevaluated[e],Power[_,_?Negative]]) Divide +
          (Length[Unevaluated[e]]-2) Times;
    ];

(* Arithmetic operators. *)

    CostFunction[x_Plus]:= (Length[Unevaluated[x]]-1) Plus;
    CostFunction[x_Times]:= (Length[Unevaluated[x]]-1) Times;

(* Ignore Lists. *)

    CostFunction[x_List]:= 0;

(* Cost remaining functions. *)

    CostFunction[x_]:= Head[Unevaluated[x]];

(* Cost required operators. *)

    DeleteCases[
      If[Head[#]===Plus,Apply[List,#],{#}]&[
        Apply[Plus,
          Map[CostFunction, Level[Unevaluated[{expr}],-2,Unevaluated] ]
        ]
      ],
      Apply[Alternatives,Times[costnull,_.]] (* Remove specified operands. *)
    ]

  ]; (* End of MainCost *)



(* Horner's rule. *)

(* Default variables. *)

Horner[p1_/p2_]:=
  Block[{$RecursionLimit=Infinity},
    HornerRule[ Expand[p1], Variables[p1] ]/
      HornerRule[ Expand[p2], Variables[p2] ]
  ];

Horner[ser_SeriesData]:=
  Block[{$RecursionLimit=Infinity},
    HornerRule[Expand[#],Variables[#]]& @ Normal[ser]
  ];

Horner[poly_]:=
  Block[{$RecursionLimit=Infinity},
    HornerRule[Expand[poly],Variables[poly]]
  ];


(* Specified variables. *)

Horner[p1_/p2_,varsp1_?VectorQ,varsp2_?VectorQ]:=
  Block[{$RecursionLimit=Infinity},
    Times[
      HornerRule[ Expand[p1], varsp1 ],
      Power[ HornerRule[ Expand[p2], varsp2 ], -1]
    ]
  ];

Horner[ser_SeriesData,vars_?VectorQ]:=
  Block[{$RecursionLimit=Infinity},
    HornerRule[ Expand[Normal[ser]], vars]
  ];

Horner[poly_,vars_?VectorQ]:=
  Block[{$RecursionLimit=Infinity},
    HornerRule[ Expand[poly], vars ]
  ];

Horner[poly_,var_]:=
  Block[{$RecursionLimit=Infinity},
    HornerRule[ Expand[poly], var ]
  ];

(* No variables found by Variables. *)

HornerRule[poly_,{}]:= poly;

HornerRule[0,_]:= 0;

(* Horner's rule for multi-variate polynomials as recursive univariate
 decomposition. *)

HornerRule[poly_,{v_,rem__}]:=
  Fold[
    SumCoeffs,
    0,

(* Pair off variable with coefficients as {var^exp,hcoeff} where hcoeff
 is the coefficient in Horner form. *)

    Thread[{
      v^Offsets[#],
      Map[ HornerRule[#,{rem}]&, GetCoeffs[poly,v,#] ]
    }]& @ Reverse[ Exponent[poly,v,Union[{##}]&] ] (* Powers sorted in descending order. *)
  ];

(* Horner's rule for uni-variate polynomials. *)

HornerRule[poly_,{v_}]:=
  Fold[
    SumCoeffs,
    0,

(* Pair off variable with coefficients as {var^exp,coeff}. *)

    Thread[{
      v^Offsets[#], GetCoeffs[poly,v,#]
    }]& @ Reverse[ Exponent[poly,v,Union[{##}]&] ] (* Powers sorted in descending order. *)
  ];

(* Calculate offset powers as successive differences (not necessarily
 incremental). *)

Offsets[e:{_}]:= e;
Offsets[e_]:= Join[-Drop[e,1],{0}] + e;

(* Accumulate Horner form of polynomial. *)

SumCoeffs[sum_,{v_,c_}]:= v (c + sum);

(* Can eventually be replaced by CoefficientList when bugs are fixed. *)

SetAttributes[GetCoeffs,Listable];
GetCoeffs[c_,v_,0]:= Coefficient[c,v,0];
GetCoeffs[c_,v_,p_]:= Coefficient[c,v^p];


(* Optimization of expressions. *)

(* Function to check the data types of Optimize options with
 error messages. *)

Options[Optimize] = {OptimizeCoefficients->False,OptimizeFunction->True,
OptimizeNull->{List},OptimizePlus->True,OptimizePower->True,
OptimizeProcedure->False,OptimizeTimes->True,OptimizeVariable->{o,Sequence}};

Optimize::args = "The `1` did not evaluate to `2`.";

opterrmsgs = {
      {"option OptimizeCoefficients","True or False"},
      {"option OptimizeFunction","True or False"},
      {"option OptimizeNull","a (possibly empty) list of symbols"},
      {"option OptimizePlus","True or False"},
      {"option OptimizePower","True, False or Binary"},
      {"option OptimizeProcedure","True or False"},
      {"option OptimizeTimes","True or False"},
      {"option OptimizeVariable","a pair of the form {Symbol,Sequence|Array}"}};

OptimizeOptionTest[opts___]:= 
  Module[{defaults,types,optcoeff,optfunc,optlist,optnull,optplus,
    optpower,optproc,opttimes,optoutv},

    defaults = Options[Optimize];

    optlist = {optcoeff,optfunc,optnull,optplus,optpower,optproc,opttimes,
      optoutv} = Map[First,defaults] /. {opts} /. defaults;

    types = {
      MatchQ[optcoeff,True|False],
      MatchQ[optfunc,True|False],
      MatchQ[optnull,{___Symbol}],
      MatchQ[optplus,True|False],
      MatchQ[optpower,Binary|True|False],
      MatchQ[optproc,True|False],
      MatchQ[opttimes,True|False],
      MatchQ[optoutv,{_Symbol,Sequence|Array}]};

    Check[
      MapThread[
        If[#1,#1,Message[Optimize::args,Apply[Sequence,#2]]]&,
        {types,opterrmsgs}
      ]; optlist,      (* Return list of option values. *)
      $Failed,         (* Option of wrong type. *)
      Optimize::args   (* Check only for these messages. *)
    ]
  ]; (* End of OptimizeOptionTest. *)


(* Default optimization symbol is o. *)

Optimize[expr_,opts___?OptionQ]:=
  Module[{optvals},
    optvals /;
      And[
        (optvals = OptimizeOptionTest[opts])=!=$Failed,
	    optvals = MainOptimize[expr,optvals]; True
	  ]
  ];



(* Define optimization function. All expression heads are wrapped
 in Hold to prevent re-evaluation during the optimization process. *)

MainOptimize[expr_,{optcoeff_,optfunc_,optnull_,optplus_,optpower_,optproc_,
  opttimes_,{outvar_,outform_}}]:=
Block[{$RecursionLimit=Infinity,BinPowDecomp,BinPowRule,downvals,HoldHead,
    index=0,keeprules,MakeBins,MakeRule,outvar,OptRule,optexpr,optvar},

(* Store sub-expression by making a new DownValue for OptRule and replace
 expression by a new optimization variable. *)

  SetAttributes[MakeRule,HoldAll];

  MakeRule[e_]:= (OptRule[e]:=OptRule[e]=#; #)& @ optvar[++index];

(* Rules for reciprocals. *)

  If[optpower=!=False,
    OptRule[e:Hold[Power][_,-1]]:= MakeRule[e];

    OptRule[Hold[Power][x_,y_?(NumberQ[#]&&Negative[#]&)]]:=
      OptRule[Hold[Power][OptRule[Hold[Power][x,-y]],-1]];
  ];

(* Rules to binary factor integer and rational powers. *)

  Switch[optpower,

(* Rules for binary decomposition of positive powers. *)

    Binary,

(* Save all binary power optimizations (even if they only occur once). *)

    MakeBins[e_]:= OptRule[e]=optvar[++index];

(* Calculate and store binary power decomposition. *)

    BinPowDecomp[p_]:= BinPowDecomp[p]=
      2^(-1+Flatten[Position[Reverse[IntegerDigits[p,2]],1]]);

(* Decompose as products of binary powers. *)

    binprod[{p_}]:= p;
      binprod[{p__}]:= OptRule[ Hold[Times][p] ];

    OptRule[Hold[Power][x_,p_Integer]]:=
      binprod[ BinPowRule[x,BinPowDecomp[p]] ];

    SetAttributes[BinPowRule,{Listable}];

(* Binary power stopping criterion. *)

    BinPowRule[e_,1]:= e;

(* Store all binary powers. *)

    BinPowRule[x_,p_]:=
      BinPowRule[x,p]=
        MakeBins[ Hold[Power][BinPowRule[x,Quotient[p,2]],2] ];

(* Rule for positive rational power with unit numerator. *)

    OptRule[e:Hold[Power][_,Rational[1,_]]]:= MakeBins[e];

(* Rules to binary decompose fractional exponents. *)

    OptRule[Hold[Power][x_,Rational[y_,z_]]]:=
      If[y<z,

(* Rule to decompose rational power as: x^(n/m) -> (x^(1/m))^n. *)

        OptRule[ Hold[Power][ OptRule[Hold[Power][x,Rational[1,z]]] ,y] ],

(* Rule to decompose rational power as: x^(c/d) -> x^(q+r/d) -> x^(q)*x^(r/d). *)

        OptRule[
          Hold[Times][
            OptRule[ Hold[Power][ x, Quotient[y,z] ] ],
            OptRule[ Hold[Power][ OptRule[Hold[Power][x,Rational[1,z]]], Mod[y,z]] ]
          ]
        ]
      ];

(* Rule for general (non-numeric) and non-binary powers. *)

    OptRule[e:Blank[Hold[Power]]]:= MakeRule[e],

(* Rule for literal powers. *)

    True,

    OptRule[e:Blank[Hold[Power]]]:= MakeRule[e],

    False,

    OptRule[e:Blank[Hold[Power]]]:= e (* Else disregard head Power. *)
  ];


(* Plus rules. *)

  If[optplus,

(* Extract numerical coefficient. *)

    If[optcoeff,
      OptRule[ Hold[Plus][n_?NumberQ,x_,y__] ]:=
        OptRule[ Hold[Plus][ n, OptRule[Hold[Plus][x,y]] ] ]
    ];

(* Store literal sub-expression. *)

    OptRule[e:Blank[Hold[Plus]]]:= MakeRule[e],

    OptRule[e:Blank[Hold[Plus]]]:= e (* Else disregard head Plus. *)

  ]; (* End of rules for expressions with head Plus. *)


(* Times rules. *)

  If[opttimes,

(* Extract numerical coefficient. *)

    If[optcoeff,
      OptRule[ Hold[Times][n_?NumberQ,x_,y__] ]:=
        OptRule[ Hold[Times][ n, OptRule[Hold[Times][x,y]] ] ]
    ];

(* Store literal sub-expression. *)

    OptRule[e:Blank[Hold[Times]]]:= MakeRule[e],

    OptRule[e:Blank[Hold[Times]]]:= e (* Else disregard head Times. *)

  ]; (* End of rules for expressions with head Times. *)


(* Rule for functions. *)

  If[optfunc,
    OptRule[e_]:= MakeRule[e],

    OptRule[e_]:= e (* Else disregard remaining expressions. *)
  ];


(* Don't make rules for specified symbols, but enable argument evaluation. *)

  Map[(HoldHead[e:Blank[#]]:= Operate[Hold,e])&, optnull ];

(* Wrap up function heads to prevent re-evaluation. Enables held
 functions arguments to be evaluated. *)

  HoldHead[e_]:= OptRule[Operate[Hold,e]];

(* Map function for creating optimization rules at non-atomic levels.
 Include special case of binary power factorisation. *)

  optexpr = If[ Head[expr]===Power&&optpower===Binary, HoldHead[#], # ]& @
              Map[ HoldHead, expr, -2 ];


(* Find occurances of temporary variables in optimized expression and
 in optimization rules. Decide which temporary variables to keep and
 which to discard by pattern look-up. *)

  downvals = DownValues[OptRule];

(* Discard rules by assigning optimization variable to subexpression.
 Head replacement makes this efficient by avoiding re-evaluation. *)

  Cases[downvals,Literal[_:>(OptRule[rhs_]=lhs_)]:>(lhs=rhs)];

(* Substitute discarded rules into optimized expression. *)

  optexpr = optexpr;

(* Save ordered list of repeated rules. *)

  keeprules =
    Sort[
      Cases[downvals,Literal[_[OptRule[rhs:_]]:>lhs:_optvar]:>lhs->rhs]
    ];


(* Format optimization variable. Reset index for consecutively numbered
 output variables. *)

  index = 0;

  If[outform===Sequence,

(* Write optimization variable as consecutive symbols o1, o2, ...
 Store indexing string function (not localised). *)

    istr[i_]:= istr[i]=ToString[i];
    (optvar[i_]:= optvar[i]=ToExpression[#<>istr[++index]])& @ ToString[outvar],

(* Write optimization variable as a consecutive array o[1], o[2], ... *)

    optvar[i_]:= optvar[i]=outvar[++index]
  ];


(* Create a Module or leave as a list of replacement rules and
 apply formatting rules for optimization variable. *)

  If[keeprules==={},
    If[optproc, expr, {{},expr} ], (* No optimization performed *)
    If[optproc,
      MakeProc[outvar,outform,ReleaseHold[keeprules],ReleaseHold[optexpr]],
      ReleaseHold[ {keeprules, optexpr} ]
    ]
  ]

]; (* End of MainOptimize.*)



(* Convert output to a held module. *)

MakeProc[optvar_,optform_,{optseq__Rule},optexpr_]:=
  (Hold[ Module[#,optseq; optexpr] ] /. Rule->Set)& @
    If[optform===Sequence, First[Thread[{optseq},Rule]], {optvar}];

(* Create an optimized procedure. *)

RemoveHold[expr_Hold]:= Apply[Unevaluated,expr];
RemoveHold[expr_]:= Unevaluated[expr];

(* SetDelayed definition. *)

Optimize /:
  (lhs_:=Optimize[expr_,opts___?OptionQ]):=
    (Unevaluated[lhs]:=#)& @
      RemoveHold[ Optimize[expr,OptimizeProcedure->True,opts] ];

(* Compile definition. *)

Optimize /:
  Compile[args_,Optimize[expr_,opts___?OptionQ],info___]:=
    Compile[Unevaluated[args],#,info]& @
      RemoveHold[ Optimize[expr,OptimizeProcedure->True,opts] ];

(* Function definition. *)

Optimize /:
  Function[args_,Optimize[expr_,opts___?OptionQ],attr___]:=
    ReleaseHold[
      Function[
        args,
        Evaluate[Optimize[expr,OptimizeProcedure->True,opts]],
        attr
      ]
    ];


End[];  (* End `Private` Context. *)

(* Protect exported symbols. *)

SetAttributes[{Cost,Horner,Optimize},ReadProtected];

Protect[Binary, Cost, CostDivide, CostExp, CostNull, CostPower, Horner, o,
Optimize, OptimizeCoefficients, OptimizeFunction, OptimizeNull, OptimizePlus,
OptimizePower, OptimizeProcedure, OptimizeTimes, OptimizeVariable];

EndPackage[];    (* End package Context. *)
