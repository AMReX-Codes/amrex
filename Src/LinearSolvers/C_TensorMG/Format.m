(* :Name: Format` *)

(* :Title: Extensions to the built-in Format rules and more... *)

(* :Author: Mark Sofroniou *)

(* :Summary:
 This package extends Mathematica's built-in format rules.
 Assignments to expressions and lists are now possible.
 The package adds definitions Assign, CAssign and FortranAssign
 and MapleAssign. Many shortcomings of the built-in formatting code
 have also been addressed, such as the limit on continuation lines
 in FORTRAN77 and assignments to Expressions.
 Code optimization is possible via the auxiliary package Optimize.m
 and the option AssignOptimize.
 The functions are primarily intended for use with the Splice command.
 When using Splice, the option FormatType->OutputForm should be
 specified.
 Interactive output within a Mathematica session is also possible
 (see also the AssignToFile option).
 All expressions are written as Strings. This enable more precise
 formatting of expressions, removing the need for text editing.
 Any Mathematica print form (e.g. TeXForm) can be specified as an
 argument of the Assign command. *)

(* :Context: Format` *)

(* :Package Version: 1.5 *)

(* :Copyright: Copyright 1992-4,  Mark Sofroniou.
 Permission is hereby granted to modify and/or make copies of
 this file for any purpose other than direct profit, or as part
 of a commercial product, provided this copyright notice is left
 intact. Sale, other than for the cost of media, is prohibited.

 Permission is hereby granted to reproduce part or all of
 this file, provided that the source is acknowledged. *)

(* :History:
 December 1994 - Modifications to MapleAssign.

 August 1994 - Version 1.5. Removed option AssignExp.
 Added array and sequence formatting for temporaries.
 Better handling of Arrays in C.

 April 1994 - Version 1.4. Removed options AssignComplexRules,
 AssignLevel, AssignMinSize, AssignRecursive. Renamed AssignNProtect
 as AssignToArray. Improved linebreaking and option/argument
 testing and scoping. Improved MapleAssign functionality.

 September 1993 - Version 1.3. Added option AssignOptimize
 and modified evaluation process accordingly. 1.3.1 minor
 modifications to FORTRAN numbering and optimization.

 August 1993 - Version 1.2. New MapleAssign accepts standard
 Mathematica input and returns analogous maple code. Options
 AssignIndex and AssignZero added.

 June 1993 - Version 1.1. Changed evaluation process.
 Removed support for log10 function and changed representation
 of rational powers in CAssign and FortranAssign.

 April 1993 - added lists of assignments.

 October 1992 - automatic breaking of long FORTRAN expressions
 and Maple format added.

 Original Version by Mark Sofroniou, July, 1992.
 
 Those who have contributed suggestions (historical order):
 Dave Withoff, Rolf Mertig, Emily Martin, Troels Petersen,
 Alain Kessi, Richard Fateman, Christophe Pichon. *)

(* :Keywords:
 Assign, CAssign, CForm, InputForm, FortranAssign, FortranForm,
 OutputForm, TeXForm. *)

(* :Source:
 Mark Sofroniou, Ph.D. Thesis, Loughborough University, Loughborough,
 Leicestershire LE11 3TU, England. *)

(* :Mathematica Version: 2.1 *)

(* :Limitations:
 This package has been developed to address limitations and problems
 encountered when using Mathematica's format rules. Commands are
 encapsulated to avoid interference with Mathematica's built-in
 definitions. Necessary rules are therefore added and removed upon
 each execution.
 Recursive code breaking may be slow for large FORTRAN
 expressions especially if too strict a tolerance is specified.
 This is a particular weakness of FortranForm.
 Suggestions for improvements and enhancements are welcome. *)

BeginPackage["Format`","Utilities`FilterOptions`"]

Assign::usage = "Assign[lhs,rhs,outputformat,options]\n
Assign converts the assignment of lhs to rhs into specified
outputformat strings (such as TeXForm). If assignments to
expressions are not required, the lhs argument may be omitted.\n
When used with Splice, the option FormatType->OutputForm
should be specified."

CAssign::usage = "CAssign[lhs,rhs,options]\n
CAssign converts the assignment of lhs to rhs into C compatible
strings. Options enable control over the conversion process, such as
the precision of real numbers. If assignments to expressions
are not required, the lhs argument may be omitted.\n When used
with Splice, the option FormatType->OutputForm should be specified."

FortranAssign::usage = "FortranAssign[lhs,rhs,options]\n
FortranAssign converts the assignment of lhs to rhs into
FORTRAN compatible strings. Options enable control over the conversion
process. Expressions are broken up and continuation characters are
added by default. The precision of real numbers in expressions may
be specified together with single or double precision exponents.
Generic FORTRAN functions are used since compilers can infer the
function type from the precision of the argument. If assignments
to expressions are not required, the lhs argument may be omitted.\n
When used with Splice, the option FormatType->OutputForm should be
specified."

MapleAssign::usage = "MapleAssign[lhs,rhs,options]\n
converts Mathematica expressions into Maple expressions and assignments.
MapleAssign converts the assignment of lhs to rhs into strings
suitable as input to Maple. If assignments to expressions are not
required, the lhs argument may be omitted.\n When used with Splice,
the option FormatType->OutputForm should be specified."

(* Options: *)

AssignBreak::usage = "AssignBreak specifies how long lines of
code should be broken up. AssignBreak may evaluate to a List of
{linewidth,string} or False. One of the string characters is
assumed to be \\n."

AssignCase::usage = "AssignCase specifies whether case conversion
of characters should be performed. AssignCase may evaluate to
Default, Lower or Upper."

AssignEnd::usage = "AssignEnd is a string appended to expressions.
It can be used to add C-style statement delimiters (AssignEnd->\";\")
or separate multiple expressions (AssignEnd->\"\\n\")."

AssignFortranNumbers::usage= "AssignFortranNumbers specifies whether
real numbers should be formatted in standard single or double precision
FORTRAN notation (e.g. 7.2d0, 10.3e1). Its value may be True or False
(uses default exponentiation)."

AssignFunction::usage = "The message AssignFunction::undef is generated
whenever a non ANSI C or FORTRAN function is encountered for the first time.
The message may be suppressed using Off."

AssignHyperbolic::usage = "AssignHyperbolic is an option of
CAssign and FortranAssign specifying whether to transform
reciprocal and inverse hyperbolic functions (these are not
supported by some compilers). Unique principal values are
assumed by writing in terms of hyperbolic functions and/or
log and sqrt. E.g. atanh(x) = log((1+x)/(1-x))/2.
AssignHyperbolic may evaluate to True or False."

AssignIndent::usage = "AssignIndent specifies a string to prepend
to an expression."

AssignIndex::usage = "AssignIndex specifies the starting index of an
assignment array. Its value may be any positive integer or zero."

AssignLabel::usage = "AssignLabel specifies a string or positive integer
to attach to the first in a list of expressions. Less than 6 characters
or digits must be involved."

AssignMaxSize::usage = "AssignMaxSize specifies an upper bound for the
maximum number of bytes of a single expression. In FORTRAN77, there is
compiler dependent limit on the allowable number of continuation lines
an expression may occupy (typically 19). Some editors such as vi in UNIX
impose limitations on the permissible length of a single line.
AssignMaxSize specifies a limit on the number of bytes using ByteCount.
This is an approximate heuristic which is roughly proportional to the
number of characters in an expression. The setting AssignMaxSize->Infinity
ensures that no expressions are broken up. AssignMaxSize
may evaluate to a positive integer (>= 200) or Infinity. See also the option
AssignTemporary."

AssignOptimize::usage = "AssignOptimize is an option of CAssign and
FortranAssign specifying whether to generate an optimized computational
sequence. This option requires the auxiliary package Optimize.m.
The degreee of optimization performed can be set as options to the
function Optimize. AssignOptimize may evaluate to True or False."

AssignPrecision::usage = "AssignPrecision specifies the precision
of real numbers in expressions. Its value may be any positive integer
or infinity."

AssignRange::usage = "AssignRange is an option of CAssign and FortranAssign.
AssignRange is used to check the range of real and integer numbers in an
expression. Numbers are checked against IEEE standards for single and double
precision according to the setting of AssignPrecision. AssignRange may evaluate
to True or False."

AssignReplace::usage = "AssignReplace specifies a list of String
replacement rules. AssignReplace can be used to compact expressions
by AssignReplace->{\" \"->\"\"}."

AssignTemporary::usage = "AssignTemporary specifies the name and format
of the temporary assignment variable used. Temporary variables are introduced
in order to break up large expressions when specified bounds are exceeded.
For example, specifying {t,Sequence} introduces the variables t1, t2,... etc.
User definitions for a symbol may interfere with the assigment process.
AssignTemporary may evaluate to an empty list or a pair {var,form} where var
is a symbol or a string and form is either Array or Sequence."

AssignTemporaryIndex::usage = "AssignTemporaryIndex stores the maximum number
of temporary variables introduced during each assignment. This is useful for
array dimensioning."

AssignToArray::usage= "AssignToArray is used to convert
Mathematica arrays and functions into arrays in C and FORTRAN.
Arguments are protected from N and maintained in exact form.
AssignToArray may evaluate to any list of symbols."

AssignToFile::usage = "AssignToFile specifies the name of an output
file to write results. Any previous contents of the file will be
overwritten. AssignToFile may evaluate to any string. See also the
Mathematica function Splice."

AssignTrig::usage = "AssignTrig is an option of CAssign and
FortranAssign specifying whether to transform reciprocal and
inverse trigonometric functions (these are not supported by
some compilers). Unique principal values in terms of trigonometric
functions are assumed. E.g. cot(x) = tan(1/x).
AssignTrig may evaluate to True or False."

AssignZero::usage = "AssignZero specifies whether zero-valued
elements in an array should be assigned or removed. This is
useful when assigning large arrays in ANSI C (default values
are zero). AssignZero may evaluate to True or False."


(* Set general error message for arguments. *)

Assign::args = "The `1` did not evaluate to `2`."
CAssign::args = "The `1` did not evaluate to `2`."
FortranAssign::args = "The `1` did not evaluate to `2`."
MapleAssign::args = "The `1` did not evaluate to `2`."

AssignFunction::undef = "Expression contains the function
`1` which is not part of the ANSI `2` standard."

AssignOptimize::fail = "Unable to optimize expression - check
default options for Optimize. Continuing with unoptimized
expression."

(* Symbols used to format real numbers in FORTRAN. *)

d::usage="d is the exponent used to format double precision
numbers in FortranAssign.";
e::usage="e is the exponent used to format single precision
numbers in FortranAssign.";

(* Symbols used to format CAssign, FortranAssign and MapleAssign functions. *)

abs::usage="The symbol abs is used to format Abs in CAssign, FortranAssign
and MapleAssign.";
acos::usage="The symbol acos is used to format ArcCos in CAssign and
FortranAssign.";
aimag::usage="The symbol aimag is used in the test for ANSI compatible
functions in FortranAssign.";
aint::usage="The symbol aint is used in the test for ANSI compatible
functions in FortranAssign.";
alog::usage="The symbol alog is used in the test for ANSI compatible
functions in FortranAssign.";
alog10::usage="The symbol alog10 is used in the test for ANSI compatible
functions in FortranAssign.";
amax0::usage="The symbol amax0 is used in the test for ANSI compatible
functions in FortranAssign.";
amax1::usage="The symbol amax1 is used in the test for ANSI compatible
functions in FortranAssign.";
amin0::usage="The symbol amin0 is used in the test for ANSI compatible
functions in FortranAssign.";
amin1::usage="The symbol amin1 is used in the test for ANSI compatible
functions in FortranAssign.";
amod::usage="The symbol amod is used in the test for ANSI compatible
functions in FortranAssign.";
and::usage="The symbol and is used to format And in MapleAssign.";
anint::usage="The symbol anint is used in the test for ANSI compatible
functions in FortranAssign.";
arccos::usage="The symbol arccos is used to format ArcCos in MapleAssign.";
acosh::usage="The symbol acosh is used to format ArcCosh in CAssign and
FortranAssign. If the compiler does not support an acosh function,
then the option AssignHyperbolic may be used.";
Ai::usage="The symbol Ai is used to format AiryAi in MapleAssign.";
arccosh::usage="The symbol arccosh is used to format ArcCosh in MapleAssign.";
arccot::usage="The symbol arccot is used to format ArcCot in MapleAssign.";
arccoth::usage="The symbol arccoth is used to format ArcCoth in MapleAssign.";
arccsc::usage="The symbol arccsc is used to format ArcCsc in MapleAssign.";
arccsch::usage="The symbol arccsch is used to format ArcCsch in MapleAssign.";
arcsec::usage="The symbol arcsec is used to format ArcSec in MapleAssign.";
arcsech::usage="The symbol arcsech is used to format ArcSech in MapleAssign.";
asin::usage="The symbol asin is used to format ArcSin in CAssign and
FortranAssign.";
arcsin::usage="The symbol arcsin is used to format ArcSin in MapleAssign.";
asinh::usage="The symbol asinh is used to format ArcSinh in CAssign and 
FortranAssign. If the compiler does not support an asinh function,
then the option AssignHyperbolic may be used.";
arcsinh::usage="The symbol arcsinh is used to format ArcSinh in MapleAssign.";
atan::usage="The symbol atan is used to format ArcTan in CAssign and
FortranAssign. If the compiler does not support an atanh function,
then the option AssignHyperbolic may be used.";
arctan::usage="The symbol arctan is used to format ArcTan in MapleAssign.";
atan2::usage="The symbol atan2 is used in the test for ANSI compatible
functions in CAssign and FortranAssign.";
atanh::usage="The symbol atanh is used to format ArcTanh in CAssign and
FortranAssign.";
arctanh::usage="The symbol arctanh is used to format ArcTanh in MapleAssign.";
bernoulli::usage="The symbol bernoulli is used to format BernoulliB in
MapleAssign.";
Bi::usage="The symbol Bi is used to format AiryBi in MapleAssign.";
binomial::usage="The symbol binomial is used to format Binomial in
MapleAssign.";
cabs::usage="The symbol cabs is used in the test for ANSI compatible
functions in FortranAssign.";
ccos::usage="The symbol ccos is used in the test for ANSI compatible
functions in FortranAssign.";
ceil::usage="The symbol ceil is used to format Round in CAssign.";
cexp::usage="The symbol cexp is used in the test for ANSI compatible
functions in FortranAssign.";
char::usage="The symbol char is used in the test for ANSI compatible functions
in FortranAssign.";
Ci::usage="The symbol Ci is used to format CosIntegral in MapleAssign.";
clog::usage="The symbol clog is used in the test for ANSI compatible
functions in FortranAssign.";
cmplx::usage="The symbol cmplx is used in the test for ANSI compatible
functions in FortranAssign.";
collect::usage="The symbol collect is used to format Collect in MapleAssign.";
conjg::usage="The symbol cmplx is used to format Conjugate in FortranAssign.";
cos::usage="The symbol cos is used to format Cos in CAssign, FortranAssign
and MapleAssign.";
cosh::usage="The symbol cosh is used to format Cosh in CAssign, FortranAssign
and MapleAssign.";
cot::usage="The symbol cot is used to format Cot in MapleAssign.";
coth::usage="The symbol coth is used to format Coth in MapleAssign.";
csc::usage="The symbol csc is used to format Csc in MapleAssign.";
csch::usage="The symbol csch is used to format Csch in MapleAssign.";
csin::usage="The symbol csin is used in the test for ANSI compatible
functions in FortranAssign.";
csqrt::usage="The symbol csqrt is used in the test for ANSI compatible
functions in FortranAssign.";
dabs::usage="The symbol dabs is used in the test for ANSI compatible
functions in FortranAssign.";
dacos::usage="The symbol dacos is used in the test for ANSI compatible
functions in FortranAssign.";
dasin::usage="The symbol dasin is used in the test for ANSI compatible
functions in FortranAssign.";
datan::usage="The symbol datan is used in the test for ANSI compatible
functions in FortranAssign.";
datan2::usage="The symbol datan is used in the test for ANSI compatible
functions in FortranAssign.";
dble::usage="The symbol dble is used in the test for ANSI compatible functions
in FortranAssign.";
dcos::usage="The symbol dcos is used in the test for ANSI compatible
functions in FortranAssign.";
dcosh::usage="The symbol dcosh is used in the test for ANSI compatible
functions in FortranAssign.";
ddim::usage="The symbol ddim is used in the test for ANSI compatible
functions in FortranAssign.";
denom::usage="The symbol denom is used to format Denominator in MapleAssign.";
dexp::usage="The symbol dexp is used in the test for ANSI compatible
functions in FortranAssign.";
diff::usage="The symbol diff is used to format D in MapleAssign.";
dilog::usage="The symbol dilog is used to format PolyLog in MapleAssign.";
dim::usage="The symbol dim is used in the test for ANSI compatible functions in
FortranAssign.";
dint::usage="The symbol dint is used in the test for ANSI compatible
functions in FortranAssign.";
div::usage="The symbol div is used in the test for ANSI compatible functions in
CAssign.";
dlog::usage="The symbol dlog is used in the test for ANSI compatible
functions in FortranAssign.";
dlog10::usage="The symbol dlog10 is used in the test for ANSI compatible
functions in FortranAssign.";
dmax1::usage="The symbol dmax1 is used in the test for ANSI compatible
functions in FortranAssign.";
dmin1::usage="The symbol dmin1 is used in the test for ANSI compatible
functions in FortranAssign.";
dmod::usage="The symbol dmod is used in the test for ANSI compatible
functions in FortranAssign.";
dnint::usage="The symbol dnint is used in the test for ANSI compatible
functions in FortranAssign.";
dprod::usage="The symbol dprod is used in the test for ANSI compatible
functions in FortranAssign.";
dsign::usage="The symbol dsign is used in the test for ANSI compatible
functions in FortranAssign.";
dsin::usage="The symbol dsin is used in the test for ANSI compatible
functions in FortranAssign.";
dsinh::usage="The symbol dsinh is used in the test for ANSI compatible
functions in FortranAssign.";
dsqrt::usage="The symbol dsqrt is used in the test for ANSI compatible
functions in FortranAssign.";
dtan::usage="The symbol dtan is used in the test for ANSI compatible
functions in FortranAssign.";
dtanh::usage="The symbol dtanh is used in the test for ANSI compatible
functions in FortranAssign.";
Ei::usage="The symbol Ei is used to format ExpIntegral in MapleAssign.";
erf::usage="The symbol erf is used to format Erf in MapleAssign.";
erfc::usage="The symbol erfc is used to format Erfc in MapleAssign.";
euler::usage="The symbol euler is used to format EulerE in MapleAssign.";
evalf::usage="The symbol evalf is used to format N in MapleAssign.";
exp::usage="The symbol exp is used to format Exp in CAssign, FortranAssign
and MapleAssign.";
expand::usage="The symbol expand is used to format Expand in MapleAssign.";
fabs::usage="The symbol fabs is used in the test for ANSI compatible functions
in CAssign.";
factor::usage="The symbol factor is used to format Factor in MapleAssign.";
false::usage="The symbol false is used to format False in MapleAssign.";
float::usage="The symbol float is used in the test for ANSI compatible
functions in FortranAssign.";
floor::usage="The symbol floor is used to format Floor in CAssign.";
fmod::usage="The symbol fmod is used in the test for ANSI compatible functions
in CAssign.";
frexp::usage="The symbol frexp is used in the test for ANSI compatible
functions in CAssign.";
fsolve::usage="The symbol fsolve is used to format NSolve in MapleAssign.";
GAMMA::usage="The symbol GAMMA is used to format Gamma in MapleAssign.";
iabs::usage="The symbol iabs is used in the test for ANSI compatible
functions in FortranAssign.";
ichar::usage="The symbol ichar is used in the test for ANSI compatible
functions in FortranAssign.";
idim::usage="The symbol idim is used in the test for ANSI compatible
functions in FortranAssign.";
idint::usage="The symbol idint is used in the test for ANSI compatible
functions in FortranAssign.";
idnint::usage="The symbol idnint is used in the test for ANSI compatible
functions in FortranAssign.";
ifix::usage="The symbol ifix is used in the test for ANSI compatible
functions in FortranAssign.";
index::usage="The symbol index is used in the test for ANSI compatible
functions in FortranAssign.";
infinity::usage="The symbol infinity is used to format ComplexInfinity,
and Infinity in MapleAssign.";
int::usage="The symbol int is used to format Floor in FortranAssign and
Integrate and NIntegrate in MapleAssign.";
isign::usage="The symbol isign is used in the test for ANSI compatible
functions in FortranAssign.";
labs::usage="The symbol labs is used in the test for ANSI compatible
functions in CAssign.";
ldexp::usage="The symbol ldexp is used in the test for ANSI compatible
functions in CAssign.";
ldiv::usage="The symbol ldiv is used in the test for ANSI compatible
functions in CAssign.";
len::usage="The symbol len is used in the test for ANSI compatible
functions in FortranAssign.";
lge::usage="The symbol lge is used in the test for ANSI compatible
functions in FortranAssign.";
lgt::usage="The symbol lgt is used in the test for ANSI compatible
functions in FortranAssign.";
lle::usage="The symbol lle is used in the test for ANSI compatible
functions in FortranAssign.";
llt::usage="The symbol llt is used in the test for ANSI compatible
functions in FortranAssign.";
lnGAMMA::usage="The symbol lnGAMMA is used to format LogGamma in MapleAssign.";
log::usage="The symbol log is used to format Log in CAssign, FortranAssign
and MapleAssign.";
log10::usage="The symbol log10 is used in the test for ANSI compatible
functions in CAssign and FortranAssign.";
map::usage="The symbol map is used to format Map in MapleAssign.";
max::usage="The symbol max is used to format Max in FortranAssign and
MapleAssign.";
max0::usage="The symbol max0 is used in the test for ANSI compatible
functions in FortranAssign.";
max1::usage="The symbol max1 is used in the test for ANSI compatible
functions in FortranAssign.";
min::usage="The symbol min is used to format Min in FortranAssign and
MapleAssign.";
min0::usage="The symbol min0 is used in the test for ANSI compatible
functions in FortranAssign.";
min1::usage="The symbol min1 is used in the test for ANSI compatible
functions in FortranAssign.";
mod::usage="The symbol mod is used to format Mod in CAssign, FortranAssign
and MapleAssign.";
mtaylor::usage="The symbol mtaylor is used to format Series in MapleAssign.";
modf::usage="The symbol modf is used in the test for ANSI compatible
functions in CAssign.";
nint::usage="The symbol nint is used to format Round in FortranAssign.";
nops::usage="The symbol nops is used to format Length in MapleAssign.";
normal::usage="The symbol normal is used to format Together in MapleAssign.";
not::usage="The symbol not is used to format Not in MapleAssign.";
NULL::usage="The symbol NULL is used to format Null in MapleAssign.";
num::usage="The symbol num is used to format Numerator in MapleAssign.";
op::usage="The symbol op is used to format Part in MapleAssign. op(0,expr)
is analogous to Head[expr]";
or::usage="The symbol or is used to format Or in MapleAssign.";
pow::usage="The symbol pow is used to format Power in CAssign.";
product::usage="The symbol product is used to format Product in MapleAssign.";
Psi::usage="The symbol Psi is used to format PolyGamma in MapleAssign.";
rand::usage="The symbol rand is used to format Random in CAssign.";
real::usage="The symbol real is used to format Re in FortranAssign.";
RootOf::usage="The symbol RootOf is used to format Roots in MapleAssign.";
round::usage="The symbol round is used to format Round in MapleAssign.";
sec::usage="The symbol sec is used to format Sec in MapleAssign.";
sech::usage="The symbol sech is used to format Sech in MapleAssign.";
series::usage="The symbol series is used to format Series in MapleAssign.";
Si::usage="The symbol Si is used to format SinIntegral in MapleAssign.";
sign::usage="The symbol sign is used to format Sign in MapleAssign and
in the test for ANSI compatible functions in FortranAssign.";
simplify::usage="The symbol simplify is used to format Simplify in
MapleAssign.";
sin::usage="The symbol sin is used to format Sin in CAssign, FortranAssign
and MapleAssign.";
sinh::usage="The symbol sinh is used to format Sinh in CAssign, FortranAssign
and MapleAssign.";
sngl::usage="The symbol sngl is used in the test for ANSI compatible
functions in FortranAssign.";
solve::usage="The symbol solve is used to format Solve in MapleAssign.";
sqrt::usage="The symbol sqrt is used to format Sqrt in CAssign, FortranAssign
and MapleAssign.";
srand::usage="The symbol srand is used in the test for ANSI compatible
functions in CAssign and FortranAssign.";
subs::usage="The symbol subs is used to format ReplaceAll in MapleAssign.";
sum::usage="The symbol sum is used to format Sum in MapleAssign.";
tan::usage="The symbol tan is used to format Tan in CAssign, FortranAssign
and MapleAssign.";
tanh::usage="The symbol tanh is used to format Tanh in CAssign, FortranAssign
and MapleAssign.";
true::usage="The symbol true is used to format True in MapleAssign.";
trig::usage="The symbol trig is used to format the option Trig->True
(used in Simplify and related functions) in MapleAssign.";

(* Symbols used to format ASCII strings. Default already exists
 in the global context. *)

Lower::usage="Lower is used to format ASCII strings via
the option AssignCase.";
Upper::usage="Upper is used to format ASCII strings via
the option AssignCase.";

Unprotect[d,e,abs,acos,acosh,Ai,aimag,aint,alog,alog10,amax0,amax1,amin0,amin1,
amod,and,anint,arccos,arccosh,arccot,arccoth,arccsc,arccsch,arcsec,arcsech,arcsin,
arcsinh,arctan,arctanh,asin,asinh,atan,atan2,atanh,bernoulli,Bi,binomial,cabs,
ccos,ceil,cexp,char,Ci,clog,cmplx,collect,conjg,cos,cosh,cot,coth,csc,csch,csin,
csqrt,dabs,dacos,dasin,datan,datan2,dble,dcos,dcosh,ddim,denom,dexp,dilog,dim,
dint,dlog,dlog10,dmax1,dmin1,dmod,dnint,dprod,dsign,dsin,dsinh,dsqrt,dtan,dtanh,
Ei,erf,erfc,euler,evalf,exp,expand,factor,factorial,false,float,fsolve,GAMMA,iabs,
ichar,idim,idint,idnint,ifix,index,infinity,int,isign,len,lge,lgt,lle,llt,log,log10,lnGAMMA,
map,max,max0,max1,min,min0,min1,mod,mtaylor,nint,not,NULL,num,op,or,pow,Psi,real,RootOf,
round,sec,sech,series,Si,sign,sin,sinh,sngl,solve,sqrt,subs,tan,tanh,true,
Lower,Upper,Assign,AssignBreak,AssignCase,AssignEnd,AssignFortranNumbers,
AssignIndent,AssignHyperbolic,AssignLabel,AssignMaxSize,AssignPrecision,
AssignRange,AssignReplace,AssignTemporary,AssignToArray,AssignToFile,AssignTrig,
CAssign,FortranAssign,MapleAssign];


Begin["`Private`"]

errmsgs = {
  {"argument lhs","a (flat) list of the same length as rhs"},
  {"option AssignBreak","False or a List of a positive integer and a string"},
  {"option AssignCase","Default, Lower, or Upper"},
  {"option AssignEnd","a string"},
  {"option AssignFortranNumbers","True or False"},
  {"option AssignHyperbolic","True or False"},
  {"option AssignIndent","a string or a positive integer"},
  {"option AssignIndex","a positive integer or zero"},
  {"option AssignLabel","a string or positive integer"},
  {"option AssignMaxSize","a positive integer (>= 200) or infinity"},
  {"option AssignOptimize","True or False"},
  {"option AssignPrecision","a positive integer or infinity"},
  {"option AssignRange","True or False"},
  {"option AssignReplace","a (possibly empty) list of string replacement rules"},
  {"option AssignTemporary","a list of the form {_Symbol|_String,Sequence|Array}"},
  {"option AssignToArray","a (possibly empty) list of symbols"},
  {"option AssignToFile","a string"},
  {"option AssignTrig","True or False"},
  {"option AssignZero","True or False"}};

(* Function to check the data types of options with error messages. *)

OptionTest[expr_,var_,assignfn_,opts___?OptionQ]:= 
  Module[{defaults=Options[assignfn],optlist,types,linbrk,acase,
    aend,fnumsQ,hypQ,indent,index,albl,amxsz,optQ,
    prec,rangeQ,arep,tvar,atoarry,atofile,trigQ,zeroQ},

    optlist =
      {linbrk,acase,aend,fnumsQ,hypQ,indent,index,albl,amxsz,optQ,
       prec,rangeQ,arep,tvar,atoarry,atofile,trigQ,zeroQ} =
         Map[First,defaults] /. {opts} /. defaults;

    types = {
If[VectorQ[var],
  MatchQ[expr,_List]&&Length[var]===Length[expr],
  If[ListQ[var],False,True]
],
MatchQ[linbrk,False|{_Integer?Positive,_String}],
MatchQ[acase,Default|Lower|Upper],
StringQ[aend],
MatchQ[fnumsQ,True|False],
MatchQ[hypQ,True|False],
StringQ[indent],
MatchQ[index,_Integer?Positive|0],
MatchQ[albl,_Integer?(0<#<100000&)|_String?(StringLength[#]<6&)],
MatchQ[amxsz,_Integer?(#>=200&)|Infinity],
MatchQ[optQ,True|False],
MatchQ[prec,_Integer?Positive|Infinity],
MatchQ[rangeQ,True|False],
MatchQ[arep,{}|{(_String->_String)...}],
MatchQ[tvar,{}|{_Symbol|_String,Sequence|Array}],
MatchQ[atoarry,{___Symbol}],
StringQ[atofile],
MatchQ[trigQ,True|False],
MatchQ[zeroQ,True|False] };

(* Add optimization variable to list of arrays and avoid duplicates. *)

If[optQ&&MatchQ[#,{_Symbol,Array}],
  optlist[[-4]] = Union[ Join[ atoarry, {First[#]} ] ]
]& @ (Optimize`OptimizeVariable /. {opts} /. Options[Optimize`Optimize]);

    Check[
      MapThread[
        If[#1,#1,Message[assignfn::args,Apply[Sequence,#2]]]&,
        {types,errmsgs}
      ]; optlist,      (* Return list of option values. *)
      $Failed,         (* Option of wrong type. *)
      assignfn::args   (* Check only for these messages. *)
    ]
  ];    (* End of OptionTest. *)



(* C assignment format. *)

SetAttributes[CAssign,HoldFirst];

Options[CAssign]:= {
AssignBreak->{Options[$Output,PageWidth][[1,2]]-2,"\\\n"},
AssignCase->Default, AssignEnd->";", AssignFortranNumbers->False, AssignHyperbolic->False,
AssignIndent->"", AssignIndex->0, AssignLabel->"", AssignMaxSize->Infinity,
AssignOptimize->False, AssignPrecision->Ceil[$MachinePrecision]-1,
AssignRange->False, AssignReplace->{" "->""}, AssignTemporary->{"t",Array},
AssignToArray->{}, AssignToFile->"", AssignTrig->True, AssignZero->True};

CAssign[lhs_:"",expr_?(!OptionQ[#]&),opts___?OptionQ]:=
  Module[{optvals},
    optvals /; 
      And[
        (optvals = OptionTest[expr,GetShape[lhs],CAssign,FilterOptions[CAssign,opts]])=!=$Failed,
        optvals = CMain[lhs,expr,optvals,{FilterOptions[Optimize`Optimize,opts]}];
        True
      ]
  ];


(* Perform assignments and code translation. Output resulting list as a 
 column and avoid string delimiters "". *)

SetAttributes[CMain,HoldFirst];

CMain[lhs_,expr_,{linbrk_,acase_,aend_,fnumsQ_,hypQ_,indent_,index_,albl_,
amxsz_,optQ_,prec_,rangeQ_,arep_,tvar_,atoarry_,atofile_,trigQ_,zeroQ_},optopts_]:=
Block[{$RecursionLimit=Infinity},
  Block[atoarry,

    AssignTemporaryIndex = 0;

(* Format C Arrays. *)

    Map[ (Format[#[i__],CForm]:=HoldForm[Part[#,i]])&, atoarry ];

    ColumnForm[
      Flatten[
        CommonAssign[
          Makelhs[lhs,CForm],
          RangeTest[
            CDefs[
              MyN[expr,prec,atoarry,CMain],
            trigQ,hypQ,optQ,prec,atoarry,optopts],
          prec,CForm,rangeQ],
          CForm,
          " = ",acase,aend,tvar,atofile,zeroQ,
          indent,index,albl,linbrk,amxsz,arep
        ]
      ]
    ] //OutputForm
  ]
]; (* End of CMain.*)



(* Define rules for C translation. *)

(* C expression head replacement. *)

SetAttributes[ApplyCDefs,Listable];

ApplyCDefs[expr_]:= CRH[Map[CRH,expr,-2]];

Literal[CRH[ArcTan[x_,y_]]]:= atan2[y,x];

(* Nest logical operators. *)

Literal[CRH[Equal[x_,y_,z__]]]:= Apply[CD[And], Map[CD[Equal][x,#]&,{y,z}] ];
Literal[CRH[e:Unequal[x_,y_,z__]]]:=
  Apply[CD[And],
    Flatten[ Table[Map[CD[Unequal][e[[i]],#]&,Drop[{x,y,z},i]],{i,Length[e]-1}] ]
  ];
Literal[CRH[(h:Greater|GreaterEqual|Less|LessEqual)[x_,y__,z_]]]:=
  Apply[CD[And],MapThread[CD[h],{{x,y},{y,z}}]];
Literal[CRH[Inequality[x_,op_,y_,z__]]]:= CD[And][CD[op][x,y],CRH[Inequality[y,z]]];
Literal[CRH[Inequality[x_,op_,y_]]]:= CD[op][x,y];

(* Recover minus sign. *)

Literal[CRH[Times[-1.,x__]]]:= CD[Times][-1,x];

(* Replace heads in remaining expressions. *)

Literal[CRH[expr_]]:= Operate[CD,expr];

(* Legal C functions. *)

cfuns = {abs,acos,AddTo,asin,atan,atan2,ceil,cos,cosh,Decrement,
  div,DivideBy,exp,fabs,floor,fmod,frexp,Increment,labs,ldexp,ldiv,
  log,log10,mod,modf,pow,Power,PreIncrement,PreDecrement,rand,srand,
  sin,sinh,sqrt,SubtractFrom,tan,tanh,TimesBy};

ANSIC[funct_]:=
 If[MemberQ[cfuns,funct],
   funct,
   Message[AssignFunction::undef,funct,"C"]; funct
 ];

(* Add C definitions. *)

SetAttributes[CDefs,{HoldAll}];

CDefs[expr_,trigQ_,hypQ_,optQ_,prec_,atoarry_,{optopts___}]:=
  Block[{Csc,Cot,Sec,ArcCsc,ArcCot,ArcSec,Csch,Coth,Sech,
    ArcCsch,ArcCoth,ArcSech,acosh,asinh,atanh,CD,pow},
    With[{one=N[1,prec],two=N[2,prec]},
      Module[{optexpr},

(* Handled correctly by CForm. *)

        CD[Times]=Times; CD[Plus]=Plus; CD[Equal]=Equal; CD[Unequal]=Unequal;
        CD[Greater]=Greater; CD[Less]=Less; CD[GreaterEqual]=GreaterEqual;
        CD[LessEqual]=LessEqual; CD[Or]=Or; CD[And]=And; CD[Not]=Not;

(* Needs additional rules. *)

        CD[Power]=Power;

(* Numeric. *)

        CD[Abs]=abs; CD[Conjugate]=conjg; CD[Floor]=floor; CD[Max]=max;
        CD[Min]=min; CD[Mod]=mod;  CD[Random]=rand; CD[Round]=ceil;
        CD[Sign]=sign; CD[Sqrt]=sqrt;

(* Trigonometric related. *)

        CD[ArcCos]=acos; CD[ArcCosh]=acosh; CD[ArcSin]=asin;
        CD[ArcSinh]=asinh; CD[ArcTan]=atan; CD[ArcTanh]=atanh;
        CD[Cos]=cos; CD[Cosh]=cosh; CD[Sin]=sin; CD[Sinh]=sinh;
        CD[Tan]=tan; CD[Tanh]=tanh; CD[Log]=log; CD[exp]=exp;

(* Numbers. *)

        CD[Complex]=Complex; CD[Rational]=Rational;

(* Arrays. *)

        Map[ (CD[#]=#)&, atoarry];

(* Legal C function? Only check head once. *)

        CD[x_]:= CD[x]=ANSIC[x];

(* Add format rules. *)

        If[trigQ,
          Csc[x_]:= Evaluate[one/CD[Sin][x]];
          Cot[x_]:= Evaluate[one/CD[Tan][x]];
          Sec[x_]:= Evaluate[one/CD[Cos][x]];
          ArcCsc[x_]:= Evaluate[CD[ArcSin][one/x]];
          ArcCot[x_]:= Evaluate[CD[ArcTan][one/x]];
          ArcSec[x_]:= Evaluate[CD[ArcCos][one/x]];
        ];

        If[hypQ,
          Csch[x_]:= Evaluate[one/CD[Sinh][x]];
          Coth[x_]:= Evaluate[one/CD[Tanh][x]];
          Sech[x_]:= Evaluate[one/CD[Cosh][x]];
          ArcCsch[x_]:= Evaluate[CD[ArcSinh][one/x]];
          ArcCoth[x_]:= Evaluate[CD[ArcTanh][one/x]];
          ArcSech[x_]:= Evaluate[CD[ArcCosh][one/x]];
          CD[ArcCosh][x_]:= Evaluate[CD[Log][x+CD[Sqrt][x^2-one]]];
          CD[ArcSinh][x_]:= Evaluate[CD[Log][x+CD[Sqrt][x^2+one]]];
          CD[ArcTanh][x_]:= Evaluate[CD[Log][(one+x)/(one-x)]/two];
          CD[ArcTanh][(x_:one)/y_]:= Evaluate[CD[Log][(y+x)/(y-x)]/two];
        ];

(* Apply formatting rules and optimize. *)

        optexpr = If[optQ, AssignOpt[#,optopts], # ]& @ ApplyCDefs[expr];

(* Add remaining formatting rules. These are applied here to avoid
 any conflict with code optimization. *)

        Block[{Power},

(* Rational powers. *)

          Power[x_,Rational[1,2]]:= Evaluate[CD[Sqrt][x]];
          Power[x_,Rational[-1,2]]:= Evaluate[one/CD[Sqrt][x]];

          Power[a_,Rational[b_,c_]]:=
            With[{nb=N[b,prec],nc=N[c,prec]}, pow[a,HoldForm[nb/nc]] ];

(* Remaining powers. *)

          Power[a_,b_?(NumberQ[#]&&#!=-1&)]:= pow[a,N[b,prec]];
          Power[a_,b_?(#=!=-1&)]:= pow[a,b];

          optexpr

        ]
      ]
    ]
  ];  (* End of CDefs. *)



(* Define FORTRAN assignment format. *)

SetAttributes[FortranAssign,HoldFirst];

Options[FortranAssign]:= {
AssignBreak->{If[#>72,72,#]&[-1+Options[$Output,PageWidth][[1,2]]],
  "\n     &  "}, AssignCase->Default, AssignEnd->"",
AssignFortranNumbers->True, AssignHyperbolic->False,
AssignIndent->"        ", AssignIndex->1, AssignLabel->"",
AssignMaxSize->5000, AssignOptimize->False,
AssignPrecision->Ceil[$MachinePrecision]-1, AssignRange->False,
AssignReplace->{" "->""}, AssignTemporary->{"t",Sequence}, AssignToArray->{},
AssignToFile->"", AssignTrig->True, AssignZero->True};


FortranAssign[lhs_:"",expr_?(!OptionQ[#]&),opts___?OptionQ]:=
  Module[{optvals},
    optvals /;
      And[
        (optvals = OptionTest[expr,GetShape[lhs],FortranAssign,
          FilterOptions[FortranAssign,opts]])=!=$Failed,
	    optvals = FMain[lhs,expr,optvals,{FilterOptions[Optimize`Optimize,opts]}];
	    True
	  ]
  ];


SetAttributes[FMain,HoldFirst];

FMain[lhs_,expr_,{linbrk_,acase_,aend_,fnumsQ_,hypQ_,indent_,index_,albl_,
amxsz_,optQ_,prec_,rangeQ_,arep_,tvar_,atoarry_,atofile_,trigQ_,zeroQ_},optopts_]:=
  Block[{d,e,$RecursionLimit=Infinity},
    Module[{newexpr,expsymb,AvoidRule=False},

      AssignTemporaryIndex = 0;

(* Attach rule for formatting real numbers with FortranForm. *)

      If[fnumsQ,
        Unprotect[Real];
        Format[expsymb,FortranForm] = If[prec>8, d, e]; (* Choose exponent. *)

(* Toggle to avoid infinite recursion formatting FORTRAN numbers. *)

        Real/: Format[r_Real,FortranForm]:=
                 (SequenceForm[First[#] 10, expsymb, -1+Last[#]]& @
                   MantissaExponent[r]) /; (AvoidRule=!AvoidRule);
      ];

(* Perform assignments and code translation. *)

      newexpr =
        CommonAssign[
          Makelhs[lhs,FortranForm],
          RangeTest[
            FortranDefs[
              MyN[expr,prec,atoarry,FMain],
            trigQ,hypQ,optQ,prec,atoarry,optopts],
          prec,FortranForm,rangeQ],
          FortranForm,
          " = ",acase,aend,tvar,atofile,zeroQ,
          indent,index,albl,linbrk,amxsz,arep
        ];

(* Remove real number rule. *)

      If[fnumsQ,
        Format[Format`Private`r$_Real,FortranForm]=.;
        Protect[Real];
      ];

(* Output a list as a column and avoid string delimiters "". *)

      ColumnForm[ Flatten[newexpr] ]  //OutputForm

    ]
  ]; (* End of FMain. *)



(* Define rules for FORTRAN translation. *)

(* FORTRAN expression head replacement. *)

SetAttributes[ApplyFortDefs,Listable];

ApplyFortDefs[expr_]:= FRH[Map[FRH,expr,-2]];

Literal[FRH[ArcTan[x_,y_]]]:= atan2[y,x];

(* Nest logical operators. *)

Literal[FRH[Equal[x_,y_,z__]]]:= Apply[FD[And], Map[FD[Equal][x,#]&,{y,z}] ];
Literal[FRH[e:Unequal[x_,y_,z__]]]:=
  Apply[FD[And],
    Flatten[ Table[Map[FD[Unequal][e[[i]],#]&,Drop[{x,y,z},i]],{i,Length[e]-1}] ]
  ];
Literal[FRH[(h:Greater|GreaterEqual|Less|LessEqual)[x_,y__,z_]]]:=
  Apply[FD[And],MapThread[FD[h],{{x,y},{y,z}}]];
Literal[FRH[Inequality[x_,op_,y_,z__]]]:= FD[And][FD[op][x,y],FRH[Inequality[y,z]]];
Literal[FRH[Inequality[x_,op_,y_]]]:= FD[op][x,y];

(* Recover minus sign. *)

Literal[FRH[Times[-1.,x__]]]:= FD[Times][-1,x];

(* Replace heads in remaining expressions. *)

Literal[FRH[expr_]]:= Operate[FD,expr];

(* Legal FORTRAN functions. *)

fortfuns = {abs,acos,aimag,aint,alog,alog10,amax0,amax1,amin0,amin1,amod,
  anint,asin,atan,atan2,cabs,ccos,cexp,char,clog,cmplx,conjg,cos,cosh,
  csin,csqrt,dabs,dacos,dasin,datan,datan2,dble,dcos,dcosh,ddim,dexp,
  dim,dint,dlog,dlog10,dmax1,dmin1,dmod,dnint,dprod,dsign,dsin,dsinh,
  dsqrt,dtan,dtanh,exp,float,iabs,ichar,idim,idint,idnint,ifix,index,
  int,isign,len,lge,lgt,lle,llt,log,log10,max,max0,max1,min,min0,min1,
  mod,nint,real,sign,sin,sinh,sngl,sqrt,tan,tanh};

ANSIF[funct_]:=
 If[MemberQ[fortfuns,funct],
   funct,
   Message[AssignFunction::undef,funct,"FORTRAN"]; funct
 ];

(* Add FORTRAN definitions. *)

SetAttributes[FortranDefs,{HoldAll}];

FortranDefs[expr_,trigQ_,hypQ_,optQ_,prec_,atoarry_,{optopts___}]:=
  Block[{Csc,Cot,Sec,ArcCsc,ArcCot,ArcSec,Csch,Coth,Sech,
    ArcCsch,ArcCoth,ArcSech,acosh,asinh,atanh,FD},
    With[{one=N[1,prec],two=N[2,prec]},
      Module[{optexpr},

(* Handled correctly by FortranForm. *)

        FD[Times]=Times; FD[Plus]=Plus; FD[Equal]=Equal; FD[Unequal]=Unequal;
        FD[Greater]=Greater; FD[Less]=Less; FD[GreaterEqual]=GreaterEqual;
        FD[LessEqual]=LessEqual; FD[Or]=Or; FD[And]=And; FD[Not]=Not;

(* Needs additional rules. *)

        FD[Power]=Power;

(* Numeric. *)

        FD[Abs]=abs; FD[Conjugate]=conjg; FD[Floor]=int; FD[Max]=max;
        FD[Min]=min; FD[Mod]=mod; FD[Re]=real; FD[Round]=nint; FD[Sqrt]=sqrt;

(* Trigonometric related. *)

        FD[ArcCos]=acos; FD[ArcCosh]=acosh; FD[ArcSin]=asin;
        FD[ArcSinh]=asinh; FD[ArcTan]=atan; FD[ArcTanh]=atanh;
        FD[Cos]=cos; FD[Cosh]=cosh; FD[Sin]=sin; FD[Sinh]=sinh;
        FD[Tan]=tan; FD[Tanh]=tanh; FD[Log]=log; FD[exp]=exp;

(* Numbers. *)

        FD[Complex]=Complex; FD[Rational]=Rational;

(* Arrays. *)

        Map[ (FD[#]=#)&, atoarry ];

(* Legal FORTRAN function? Only check head once. *)

        FD[x_]:= FD[x]=ANSIF[x];

(* Add format rules. *)

        If[trigQ,
          Csc[x_]:= Evaluate[one/FD[Sin][x]];
          Cot[x_]:= Evaluate[one/FD[Tan][x]];
          Sec[x_]:= Evaluate[one/FD[Cos][x]];
          ArcCsc[x_]:= Evaluate[FD[ArcSin][one/x]];
          ArcCot[x_]:= Evaluate[FD[ArcTan][one/x]];
          ArcSec[x_]:= Evaluate[FD[ArcCos][one/x]];
        ];

        If[hypQ,
          Csch[x_]:= Evaluate[one/FD[Sinh][x]];
          Coth[x_]:= Evaluate[one/FD[Tanh][x]];
          Sech[x_]:= Evaluate[one/FD[Cosh][x]];
          ArcCsch[x_]:= Evaluate[FD[ArcSinh][one/x]];
          ArcCoth[x_]:= Evaluate[FD[ArcTanh][one/x]];
          ArcSech[x_]:= Evaluate[FD[ArcCosh][one/x]];
          FD[ArcCosh][x_]:= Evaluate[FD[Log][x+FD[Sqrt][x^2-one]]];
          FD[ArcSinh][x_]:= Evaluate[FD[Log][x+FD[Sqrt][x^2+one]]];
          FD[ArcTanh][x_]:= Evaluate[FD[Log][(one+x)/(one-x)]/two];
          FD[ArcTanh][(x_:one)/y_]:= Evaluate[FD[Log][(y+x)/(y-x)]/two];
        ];

(* Apply formatting rules and optimize. *)

        optexpr = If[optQ,AssignOpt[#,optopts],#]& @ ApplyFortDefs[expr];

(* Add remaining formatting rules. These are applied here to avoid
 any conflict with code optimization. *)

        Block[{Power},

(* Rational powers. *)

          Power[x_,Rational[1,2]]:= Evaluate[FD[Sqrt][x]];
          Power[x_,Rational[-1,2]]:= Evaluate[one/FD[Sqrt][x]];

          Power[a_,Rational[b_,c_]]:=
            With[{nb=N[b,prec],nc=N[c,prec]}, HoldForm[a^(nb/nc)] ];

          optexpr
        ]

      ]
    ]
  ];  (* End of FortranDefs. *)




(* Define Maple assignment format. *)

SetAttributes[MapleAssign,HoldAll];

Options[MapleAssign]:= {
AssignBreak->{Options[$Output,PageWidth][[1,2]]-2,"\\\n"},
AssignCase->Default, AssignEnd->";", AssignFortranNumbers->False, AssignHyperbolic->False,
AssignIndent->"", AssignIndex->1, AssignLabel->"", AssignMaxSize->Infinity,
AssignOptimize->False, AssignPrecision->Infinity, AssignRange->False,
AssignReplace->{}, AssignTemporary->{}, AssignToArray->{},
AssignToFile->"", AssignTrig->False, AssignZero->True};


MapleAssign[lhs_:"",expr_?(!OptionQ[#]&),opts___?OptionQ]:=
  Module[{optvals},
    optvals /;
      And[
        (optvals = OptionTest[expr,GetShape[lhs],MapleAssign,opts])=!=$Failed,
	    optvals = MMain[lhs,Unevaluated[expr],optvals]; True
	  ]
  ];


SetAttributes[MMain,HoldFirst];

MMain[lhs_,expr_,{linbrk_,acase_,aend_,fnumsQ_,hypQ_,indent_,index_,albl_,
amxsz_,optQ_,prec_,rangeQ_,arep_,tvar_,atoarry_,atofile_,trigQ_,zeroQ_}]:=
  Block[{$RecursionLimit=Infinity},

    AssignTemporaryIndex = 0;

(* Perform assignments and code translation. *)

    OutputForm[
      Flatten[
        CommonAssign[
          Makelhs[lhs,InputForm],
          MyN[Evaluate[MapleDefs[expr]],prec,atoarry],
          InputForm,
          " := ",acase,aend,tvar,atofile,zeroQ,
          indent,index,albl,linbrk,amxsz,arep
        ]
      ] //ColumnForm
    ]

  ] (* End of MMain. *)



(**** Maple function redefinitions. ****)

SetAttributes[MRH,HoldAll];

ApplyMapleDefs[expr_]:= MapAll[ MRH, Unevaluated[expr] ];

(* Special Mathematica expressions. *)

Literal[MRH[Head[x_]]]:= MapleFun[Head,0,x];
Literal[MRH[Part[x_,y_]]]:= MapleFun[Part,0,x];
Literal[MRH[Part[x_,y_List]]]:= MapleFun[Part,Apply[Sequence,y],x];

Literal[MRH[Replace[x_,rep_]]]:= MapleFun[Replace,rep,x];
Literal[MRH[ReplaceAll[x_,rep_]]]:= MapleFun[ReplaceAll,rep,x];

Literal[MRH[Rule[Trig,True]]]:= trig;
Literal[MRH[Rule[x_,y_]]]:= SequenceForm[MapleArgs["=",x,y]];
Literal[MRH[Equal[x_,y_]]]:= SequenceForm[MapleArgs["=",x,y]];

(* Nest logical operators. *)

Literal[MRH[Equal[x_,y_,z__]]]:=
  Apply[MD[And], Map[SequenceForm[MapleArgs[MD[Equal],x,#]]&,{y,z}] ];
Literal[MRH[e:Unequal[x_,y_,z__]]]:=
  Apply[MD[And],
    Flatten[
      Table[
        Map[SequenceForm[MapleArgs[MD[Unequal],e[[i]],#]]&,Drop[{x,y,z},i]],
      {i,Length[e]-1}]
    ]
  ];
Literal[MRH[Unequal[x_,y_]]]:= SequenceForm[MapleArgs[MD[Unequal],x,y]];
Literal[MRH[(h:Greater|GreaterEqual|Less|LessEqual)[x_,y__,z_]]]:=
  Apply[MD[And],MapThread[SequenceForm[MapleArgs[MD[h],#1,#2]]&,{{x,y},{y,z}}]];
Literal[MRH[(h:Greater|GreaterEqual|Less|LessEqual)[x_,y_]]]:=
  SequenceForm[MapleArgs[MD[h],x,y]];
Literal[MRH[Inequality[x_,op_,y_,z__]]]:=
MD[And][SequenceForm[MapleArgs[MD[op],x,y]],MRH[Inequality[y,z]]];
Literal[MRH[Inequality[x_,op_,y_]]]:= SequenceForm[MapleArgs[MD[op],x,y]];

(* Logicals. *)

MRH[True] = true;  MRH[False] = false;

(* Trigonometric related. *)

Literal[MRH[ArcTan[x_,y_]]]:= Evaluate[MapleFun[MD[ArcTan],y,x]];
Literal[MRH[Power[MRH[E],x_]]]:= Evaluate[MapleFun[MD[Exp],x]];

(* Required for formatting. *)

Literal[MRH[x_SequenceForm|x_OutputForm]]:= x;

(* Arithmetic. *)

Literal[MRH[x_Plus|x_Power|x_Times]]:= x;

Literal[MRH[x_Complex|x_Integer|x_Rational|x_Real]]:= x;

Literal[MRH[Infinity|ComplexInfinity|DirectedInfinity[___]]]:=
  Evaluate[MD[DirectedInfinity]];

(* Lists. *)

Literal[MRH[List[x__]]]:= SequenceForm[OutputForm["["],MapleArgs[",",x],OutputForm["]"]];

(* General function head replacement. *)

special = {And,Or,Complex,D,Factorial,Integrate,Product,Sum,
  Series,Sign,NSolve,Solve};

Literal[MRH[f_[x__]]]:= If[MemberQ[special,f], MD[f][x], MapleFun[f,x] ];

(* Remaining expressions. *)

Literal[MRH[expr_]]:= MD[expr];

(**** End of Maple function redefinitions. ****)


(**** Maple function head redefinitions. ****)

    MD[Abs]=abs; MD[AiryAi]=Ai; MD[AiryBi]=Bi; MD[And]=and; MD[ArcCos]=arccos; MD[ArcCosh]=arccosh;
    MD[ArcCot]=arccot; MD[ArcCoth]=arccoth; MD[ArcCsc]=arccsc; MD[ArcCsch]=arccsch;
    MD[ArcSec]=arcsec; MD[ArcSech]=arcsech; MD[ArcSin]=arcsin; MD[ArcSinh]=arcsinh;
    MD[ArcTan]=arctan; MD[ArcTanh]=arctanh; MD[BernoulliB]=bernoulli;
    MD[Binomial]=binomial; MD[Collect]=collect; MD[Cos]=cos; MD[Cosh]=cosh;
    MD[CosIntegral]=Ci; MD[Cot]=cot; MD[Coth]=coth; MD[Csc]=csc; MD[Csch]=csch;
    MD[D]=diff; MD[Denominator]=denom; MD[DirectedInfinity]=infinity; MD[Erf]=erf;
    MD[Erfc]=erfc; MD[EulerE]=euler; MD[Exp]=exp; MD[Expand]=expand; MD[ExpIntegral]=Ei;
    MD[Factor]=factor; MD[Factorial]=factorial; MD[Gamma]=GAMMA;
    MD[Head]=op; MD[Integrate]=int;
    MD[Length]=nops; MD[Log]=log; MD[LogGamma]=lnGAMMA; MD[Map]=map; MD[Max]=max;
    MD[Min]=min; MD[Mod]=mod; MD[N]=evalf; MD[NIntegrate]=int; MD[NSolve]=fsolve; MD[Not]=not;
    MD[Null]=NULL; MD[Numerator]=num; MD[Or]=or; MD[Part]=op; MD[PolyGamma]=Psi;
    MD[PolyLog]=dilog; MD[Product]=product; MD[Replace]=subs; MD[ReplaceAll]=subs;
    MD[Roots]=RootOf; MD[Round]=round; MD[Sec]=sec; MD[Sech]=sech;
    MD[Series]=series; MD[Sign]=sign; MD[Simplify]=simplify; MD[Sin]=sin;
    MD[Sinh]=sinh; MD[SinIntegral]=Si; MD[Solve]=solve; MD[Sqrt]=sqrt; MD[Sum]=sum;
    MD[Tan]=tan; MD[Tanh]=tanh; MD[Together]=normal;

(* Logical symbols. *)

    MD[Greater]=">"; MD[GreaterEqual]=">="; MD[Less]="<"; MD[LessEqual]="<=";
    MD[Equal]="="; MD[Unequal]="<>";

(* Not yet implemented. *)

    MD[x_]:= x;

(**** End of Maple function head redefinitions. ****)




(**** Rules to convert arguments and functions to Maple form. ****)

(* Recursive argument conversion. *)

MapleArgs[str_,x_]:= x;
MapleArgs[str_,x_,y__]:= Sequence[x,OutputForm[str],MapleArgs[str,y]];

(* Function conversion. *)

MapleFun[f_,x__]:= SequenceForm[MD[f],OutputForm["("],MapleArgs[",",x],OutputForm[")"]];

MapleDerivs[SequenceForm[_,a_,_,b_,_]]:= SequenceForm[a,OutputForm["$"],b];

MapleDerivs[x_]:= x;

MapleRange[SequenceForm[_,x_,_,r2_,_],str_]:=
  SequenceForm[ x,OutputForm["="],MapleArgs[str,1,r2] ];

MapleRange[SequenceForm[_,x_,_,r1_,_,r2_,_],str_]:=
  SequenceForm[ x,OutputForm["="],MapleArgs[str,r1,r2] ];

MapleRange[x__,y_]:= x;

MapleSeries[x_]:=
  Module[{mterms = Map[MapleRange[#,","]&,x],trunc},
    trunc = Max[Map[Part[#,-1]&,mterms]]; (* Max truncation order. *)
    SequenceForm[
      OutputForm["["],
      Apply[Sequence, Insert[ Map[Drop[#,-1]&,mterms], OutputForm["]"], {-1,-2}] ],
      trunc
    ]
  ];

(**** End of rules to convert arguments and functions to Maple form. ****)


(* Define rules for Maple translation. *)

SetAttributes[MapleDefs,{HoldAll}];

MapleDefs[expr_]:=
  Block[{and,or,diff,factorial,int,product,sum,series,sign,fsolve,solve},

(* Add formatting rules. *)

    and[x__]:= SequenceForm[MapleArgs[" and ",x]];
    or[x__]:= SequenceForm[MapleArgs[" or ",x]];
    diff[x_,a_]:= MapleFun[diff,x,MapleDerivs[a]];
    diff[x_,a__,b_]:= diff[ diff[x,a], MapleDerivs[b]];
    factorial[x_]:= SequenceForm[x,OutputForm["!"]];
    int[f_,x_]:= MapleFun[int,f,MapleRange[x,".."]];
    int[f_,x__,y_]:= MapleFun[int,int[f,x], MapleRange[y,".."]];
    product[x_,y_]:= MapleFun[product,x,MapleRange[y,".."]];
    product[x_,y__,z_]:= MapleFun[product,product[x,y],MapleRange[z,".."]];
    sum[x_,y_]:= MapleFun[sum,x,MapleRange[y]];
    sum[x_,y__,z_]:= MapleFun[sum,sum[x,y],MapleRange[z]];
    series[x_,y_]:= MapleFun[series,x,MapleRange[y,","]];
    series[x_,y__,z_]:= MapleFun[mtaylor,x,MapleSeries[{y,z}]];
    sign[x_Complex]:= MapleFun[signum,x];
    sign[x_]:= MapleFun[sign,x];
    fsolve[x__]:= MapleFun[fsolve,x] /. {"["->"{","]"->"}"};
    solve[x__]:= MapleFun[solve,x] /. {"["->"{","]"->"}"};

    ApplyMapleDefs[Unevaluated[expr]]

  ]; (* End of Maple Defs. *)




(* Define assignment for specified format. *)

SetAttributes[Assign,HoldFirst];

Options[Assign]:= {
AssignBreak->{Options[$Output,PageWidth][[1,2]]-1,"\n"},
AssignCase->Default, AssignEnd->"", AssignFortranNumbers->False, AssignHyperbolic->False,
AssignIndent->"", AssignIndex->1, AssignLabel->"", AssignMaxSize->Infinity,
AssignOptimize->False, AssignPrecision->Infinity, AssignRange->False,
AssignReplace->{}, AssignTemporary->{}, AssignToArray->{},
AssignToFile->"", AssignTrig->False, AssignZero->True};


Assign[lhs_:"",expr_?(!OptionQ[#]&),form_?(!OptionQ[#]&),opts___?OptionQ]:=
  Module[{optvals},
    optvals /;
      And[
        (optvals = OptionTest[expr,GetShape[lhs],Assign,opts])=!=$Failed,
	    optvals = AMain[lhs,expr,form,optvals]; True
	  ]
  ];


(* Perform assignments and code translation. *)

SetAttributes[AMain,HoldFirst];

AMain[lhs_,expr_,form_,{linbrk_,acase_,aend_,fnumsQ_,hypQ_,indent_,index_,
albl_,amxsz_,optQ_,prec_,rangeQ_,arep_,tvar_,atoarry_,atofile_,trigQ_,zeroQ_}]:=
Block[{$RecursionLimit=Infinity},

  AssignTemporaryIndex = 0;

  OutputForm[
    Flatten[
      CommonAssign[
        Makelhs[lhs,form],
        MyN[expr,prec,atoarry],
        form,
        " = ",acase,aend,tvar,atofile,zeroQ,
        indent,index,albl,linbrk,amxsz,arep
      ]
    ] //ColumnForm
  ]
]; (* End of AMain. *)



(* Define main assignment formatting. *)

(* This definition ensures compatibility with the code optimization
 package Optimize.m - which returns a list of replacement rules
 and an optimized expression or list of expressions. *)

CommonAssign[lhs_,{optrules:{__Rule},expr_},form_,args__]:=
  Apply[
    CommonAssign[ Join[Makelhs[#1,form],{lhs}],
      Join[#2,{expr}],form,args]&,
    Thread[optrules,Rule]
  ];

(* Assignments to single expressions. *)

CommonAssign[lhs_,expr_,form_,eqstr_,acase_,aend_,tvar_,atofile_,
  zeroQ_,indent_,index_,albl_,linbrk_,amxsz_,arep_]:=
  Module[{outchan,lbllen,redexpr,strings},

(* Remove zero-valued expressions and add array indices to lhs. *)

    redexpr = RemoveZeros[lhs,expr,form,index,!zeroQ];

(* Break up long expressions and convert to strings. *)

    strings =
      Map[BreakExpression[#,eqstr,amxsz,tvar,form]&,redexpr];

(* Apply string replacement rules and indentation/termination strings. *)

    strings =
      Map[
        StringJoin[ indent, #, aend ]&,
        StringReplace[ Flatten[{strings}],arep ]
      ];

(* Attach a label to the first expression. *)

    label = ToString[albl];  lbllen = StringLength[label];

    If[ lbllen =!=0,
      strings[[1]] = StringJoin[label,StringDrop[strings[[1]],lbllen]]
    ];

(* Convert to Lower, or Upper case. *)

    Switch[acase,
      Upper,strings = Map[ToUpperCase,strings],
      Lower,strings = Map[ToLowerCase,strings]
    ];

(* Add continuation characters to break up long lines of code. *)

    If[ListQ[linbrk], strings = Map[BreakLines[#,linbrk]&,strings] ];

(* Output results to a file. *)

    If[atofile=!="",
      outchan = OpenWrite[atofile,FormatType->OutputForm];
      Write[outchan, strings //ColumnForm];
      Close[outchan] ];

    strings    (* Output. *)

  ];  (* End of CommonAssign. *)



(* Add index and delete zero-valued elements. *)

(* Remove zero-valued elements. The variable index is used as
 an offset for the starting index. *)

(* Lists of assignments (remove outer List). *)

RemoveZeros[lhs_List,rhs_List,format_,index_,remzero_]:=
  Apply[Join,MapThread[ RemoveZeros[#1,#2,format,index,remzero]&, {lhs,rhs} ]];

(* No zeros removed. *)

RemoveZeros["",rhs_List,format_,index_,False]:= Thread[{"",Flatten[rhs]}];

RemoveZeros[lhs_,rhs_List,format_,index_,False]:=
  JoinIndices[lhs,FindPositions[rhs,Null],Flatten[rhs],format,index];

RemoveZeros[lhs_,rhs_,format_,index_,False]:= {{lhs,rhs}};

(* Zeros removed. *)

(* Check for expressions with no non-zeros present. *)

AssignZero::continue = "Expression encountered with no non-zero elements.
Continuing with zero assignments.";

RemoveZeros[lhs_,0,format_,index_,True]:=
  (Message[AssignZero::continue]; {{lhs,0}});

RemoveZeros["",rhs_List,format_,index_,True]:=
  Module[{redrhs = Flatten[Delete[rhs,Position[rhs,0,Heads->False]]]},
    If[redrhs==={}, Message[AssignZero::continue]; redrhs = Flatten[rhs]];
    Thread[{"",redrhs}]
  ];

RemoveZeros[lhs_,rhs_List,format_,index_,True]:=
  Module[{redrhs = Flatten[Delete[rhs,Position[rhs,0,Heads->False]]],posntest=0},
    If[redrhs==={},
      Message[AssignZero::continue]; redrhs = Flatten[rhs]; posntest=Null;];
    JoinIndices[lhs,FindPositions[rhs,posntest],redrhs,format,index]
  ];

RemoveZeros[lhs_,rhs_,format_,index_,True]:= {{lhs,rhs}};

(* Find positions of non-null and non-zero elements. *)

FindPositions[expr_,Null]:=
  Position[ Function[x,UnsameQ[x,Null],Listable][expr],True,Heads->False];

FindPositions[expr_,0]:=
  Position[ Function[x,UnsameQ[x,0]&&UnsameQ[x,Null],Listable][expr],True,Heads->False];

(* Join index string to lhs string. *)

JoinIndices[lhs_,posns_List,rhs_,format_,index_]:=
  (MapThread[
    {StringJoin[ lhs, StringIndex[#1,format] ],#2}&,
    {posns+index-1,Flatten[rhs]}
  ]);

(* StringIndex is used to convert the position into the
 appropriate form for an array element. *)

StringIndex[psn_,CForm]:= 
  StringReplace[ ToString[psn] ,{"{"->"[","}"->"]",", "->"]["}];

(* Default is FORTRAN case. *)

StringIndex[psn_,_]:= 
  StringReplace[ ToString[psn] ,{"{"->"(","}"->")"}];



(* Break up large expressions into sub-expressions. *)

BreakUp[subexpr_,maxlen_,temp_]:=
  If[LengthTest[subexpr,maxlen],
    Fragment[subexpr,maxlen,temp],
    subexpr (* else *)
  ];

(* Add temporary variable and sub-expression to list. *)

AddTemp[temp_,subexpr_]:= (parts = Join[parts,{{temp,subexpr}}]; temp);

(* Test used for permissible sub-expression size. *)

LengthTest[expr_,maxlen_]:= ByteCount[expr]>maxlen;

(* Ignore numeric exponents, temporary variables etc. *)

Fragment[subexpr:(_temp|_?AtomQ),maxlen_,tmpvar_]:= subexpr;

(* Binary decomposition of Plus and Times. *)

Fragment[subexpr:(_Plus|_Times),maxlen_,temp_]:=
  If[LengthTest[subexpr,maxlen],
    With[{quo = Quotient[Length[subexpr],2]},
      BreakUp[
        Head[subexpr][
          Fragment[Take[subexpr,quo],maxlen,temp],
          Fragment[Drop[subexpr,quo],maxlen,temp]
        ],
      maxlen,temp]
    ],
    AddTemp[temp[++index],subexpr] (* else *)
  ];

(* n-ary decomposition of remaining functions. *)

Fragment[subexpr_,maxlen_,temp_]:=
  If[LengthTest[subexpr,maxlen],
    BreakUp[Map[Fragment[#,maxlen,temp]&,subexpr],maxlen,temp],
    AddTemp[temp[++index],subexpr] (* else *)
  ];


(* Recursively decompose expression. *)

BreakExpression[args_,eqstr_,Infinity,_,form_]:=
  MyFormat[args,eqstr,form];

Assign::notemp = "No temporary variable was specified. Continuing
with original expression.";

BreakExpression[args_,eqstr_,maxlen_,{}|{"",_},form_]:=
  (Message[Assign::notemp]; MyFormat[args,eqstr,form]);

BreakExpression[{lhs_,expr_},eqstr_,maxlen_,{tvar_,tform_},form_]:=
  Block[{index=0,parts={},$RecursionLimit=Infinity},
    Module[{outexpr,tmp},

(* Array or sequence of temporaries. *)

      If[tform===Array,
        If[form===CForm,
          Format[tmp[i_],form]:=HoldForm[Part[#,i]],
          Format[tmp[i_],form]:=#[i]
        ],
        Format[tmp[i_],form]:= SequenceForm[#,i]
      ]& @ ToExpression[ToString[tvar]];

(* Recursively break up expression and re-use temporary variables. *)

      outexpr = BreakUp[expr,maxlen,tmp];

(* Store maximum number of temporaries introduced. *)

      If[ index>AssignTemporaryIndex, AssignTemporaryIndex=index ];

(* Output list of temp strings and final expression. *)

      {Map[ MyFormat[#,eqstr,form]&, parts ],
       MyFormat[{lhs,outexpr},eqstr,form]}
    ]
  ]; (* End of BreakExpression *)


(* Convert the expression to a string of appropriate form. *)

MyFormat[{"",expr_},_,form_]:= ToString[ expr, FormatType->form ];

MyFormat[{lhs_String,expr_},eqstr_,form_]:=
  StringJoin[ lhs, eqstr, ToString[ expr, FormatType->form ] ];

MyFormat[{lhs_,expr_},eqstr_,form_]:=
  StringJoin[
    ToString[ lhs, FormatType->form ],
    eqstr,
    ToString[ expr, FormatType->form ]
  ];


(* Break up long lines of code and add continuation characters. *)

BreakLines[string_,{lineln_,indstr_}]:= 
  Module[{indlen=StringLength[indstr],linelen,numlines,
    strlen=StringLength[string]},

(* (\n is counted as one character). *)

    linelen = lineln - indlen + 1;

    If[strlen >=linelen,
      numlines = Floor[(strlen - (linelen+indlen))/linelen];
      StringJoin[
        StringTake[string,indlen-1],  (* first part *)
          Table[
            StringTake[
              string,                         (* middle part *)
              {linelen i + indlen,
               linelen (i+1) + indlen-1}
            ]<>indstr
          ,{i,0,numlines}],
          StringTake[
            string                            (* end part *)
          ,((numlines+1) linelen + indlen-1) - strlen]
        ],
        string (* else *)
      ]
    ];  (* End of BreakLines. *)


(* Determine `shape' of lhs lists. Extract elements to avoid evaluation. *)

SetAttributes[GetShape,{Listable,HoldAll}];
GetShape[_]:="";


(* Make lhs into a list of held strings. *)

SetAttributes[Makelhs,{Listable,HoldAll}];
Makelhs[lhs_String,_]:= lhs;
Makelhs[lhs_,InputForm]:= ToString[ HoldForm[lhs] ];
Makelhs[lhs_,form_]:= ToString[ form[HoldForm[lhs]] ];


(* Slightly modified version of N. Arguments to specified
 symbols are temporarily protected from N. *)

SetAttributes[MyN,HoldAll];

(* Protect exponential from N. *)

MyN[expr_,_DirectedInfinity,_,CMain|FMain]:=
  Block[{E}, E /: Power[E,x_]:= exp[x]; expr ];

(* Approximate numeric Exp. *)

MyN[args__,CMain|FMain]:=
  Block[{E}, E /: Power[E,x_?(!NumberQ[#]&)]:= exp[x]; MyN[args] ];

(* Infinite Precision. *)

MyN[expr_,_DirectedInfinity,__]:= expr;

(* Finite precision. *)

MyN[expr_,prec_,{}]:=
  If[prec===Ceil[$MachinePrecision]-1,
    #,
    # //. {r_Real:>SetPrecision[r,prec]}
  ]& @ N[ expr, prec ];

(* Finite Precision, protect array arguments from N. *)

MyN[expr_,prec_,atoarry_]:=
  Block[atoarry,
    SetAttributes[atoarry,NHoldAll];
    MyN[ expr, prec, {} ]
  ];


(* Optimize expressions. *)

AssignOpt[expr_,optopts___?OptionQ]:=
  Check[
    RuleGen[
      Optimize`Optimize[expr,optopts],
      expr
    ],
    Message[AssignOptimize::fail]; expr, (* Else proceed with unoptimized expression. *)
    Optimize`Optimize::args              (* Check only for this message. *)
  ];

(* Check for optimization. *)

RuleGen[{{},expr_},expr_]:= expr;     (* No optimization. *)
RuleGen[optexpr_,_]:= optexpr;        (* Optimization. *)



(* Test range of real and integer numbers. *)

spfortran = {2.^-126,2.^127,HoldForm[-2^31],HoldForm[2^31-1],
             -2^31,2^31-1,"single"};
spc = {2.^-125,2.^128,HoldForm[-2^31],HoldForm[2^31-1],-2^31,2^31-1,"single"};
dpfortran = {2.^-1022,2.^1023,HoldForm[-2^63],HoldForm[2^63-1],
             -2^63,2^63-1,"double"};
dpc = {2.^-1021,2.^1024,HoldForm[-2^63],HoldForm[2^63-1],-2^63,2^63-1,"double"};

SetAttributes[RangeTest,HoldFirst];

RangeTest[expr_,nprec_,form_,False]:= expr;

RangeTest[expr_,nprec_,FortranForm,True]:=
  CheckRange[expr,FortranForm,If[nprec<=8, spfortran, dpfortran]];

RangeTest[expr_,nprec_,CForm,True]:=
  CheckRange[expr,CForm,If[nprec<=8, spc, dpc]];

AssignRange::float = "Expression contains machine numbers outside
the permissible range `1` to `2` for IEEE `3` precision.";

AssignRange::integer = "Expression contains integers outside the
permissible range `1` to `2` which cannot be represented in IEEE
`3` precision and have been converted to floating point numbers.";

CheckRange[expr_,form_,{xrmin_,xrmax_,xihmin_,xihmax_,ximin_,ximax_,prec_}]:=
  Module[{complxQ,intmsg=True,intQ,realQ,rlmsg=True},

    realQ = r_Real?((Abs[#]>xrmax||Abs[#]<xrmin)&&rlmsg&):>
              (If[rlmsg, rlmsg=False; Message[AssignRange::float,xrmin,xrmax,prec]];
                 r);

    intQ = i_Integer?(#>ximax||#<ximin&):>
             (If[intmsg,  intmsg=False; Message[AssignRange::integer,xihmin,xihmax,prec]];
                N[i] /. realQ);

    cmplxQ = Complex[r_,i_]:>Complex[r /. {realQ,intQ},i /. {realQ,intQ}];

    expr /. {realQ,intQ,cmplxQ}
  ]; (* End of CheckRange. *)


End[];  (* End `Private` Context. *)

(* Protect exported symbols. *)

SetAttributes[{Assign,CAssign,FortranAssign,MapleAssign},ReadProtected];

Protect[d,e,abs,acos,acosh,Ai,aimag,aint,alog,alog10,amax0,amax1,amin0,amin1,
amod,and,anint,arccos,arccosh,arccot,arccoth,arccsc,arccsch,arcsec,arcsech,arcsin,
arcsinh,arctan,arctanh,asin,asinh,atan,atan2,atanh,bernoulli,Bi,binomial,cabs,
ccos,ceil,cexp,char,Ci,clog,cmplx,collect,conjg,cos,cosh,cot,coth,csc,csch,csin,
csqrt,dabs,dacos,dasin,datan,datan2,dble,dcos,dcosh,ddim,denom,dexp,dilog,dim,
dint,dlog,dlog10,dmax1,dmin1,dmod,dnint,dprod,dsign,dsin,dsinh,dsqrt,dtan,dtanh,
Ei,erf,erfc,euler,evalf,exp,expand,factor,factorial,false,float,fsolve,GAMMA,iabs,
ichar,idim,idint,idnint,ifix,index,infinity,int,isign,len,lge,lgt,lle,llt,log,log10,lnGAMMA,
map,max,max0,max1,min,min0,min1,mod,mtaylor,nint,not,NULL,num,op,or,pow,Psi,real,RootOf,
round,sec,sech,series,Si,sign,sin,sinh,sngl,solve,sqrt,subs,tan,tanh,true,
Lower,Upper,Assign,AssignBreak,AssignCase,AssignEnd,AssignFortranNumbers,
AssignIndent,AssignHyperbolic,AssignLabel,AssignMaxSize,AssignPrecision,
AssignRange,AssignReplace,AssignTemporary,AssignToArray,AssignToFile,AssignTrig,
CAssign,FortranAssign,MapleAssign];

EndPackage[];    (* End package Context. *)
