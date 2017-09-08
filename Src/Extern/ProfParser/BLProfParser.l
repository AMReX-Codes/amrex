%{
#include <stdio.h>
#include <stdlib.h>

#include <BLProfParser.tab.H>

#define YY_NEVER_INTERACTIVE 1
#define ECHO return (0);       /* for unmatched tokens */

%}

%option noyywrap yylineno

WS  [ \t\r]
ENDL    \n
DIGIT   [0-9]
ALPHA   [a-zA-z]
PLUS    [+]
MINUS    -
DQUOTE  ["]
SQUOTE  [']
SIGN    ({MINUS}|{PLUS})?{WS}*
EXPONENT    ("E"|"e"){WS}*{SIGN}{WS}*{DIGIT}+

DECIMALNUMBER   (({DIGIT}+".")|({DIGIT}*"."{DIGIT}+)){DIGIT}*{WS}*

DOT   "."
COMMA   ,
COLON   :
EQUAL   =
LPAREN  \(
RPAREN  \)
LBRACKET  \[
RBRACKET  \]
SLASH  \/
EOS     '\0'
STAR    \*
POUND \#
POUNDCOMMENT {POUND}+.*

MPITOKEN "MPI"
OMP "OMP"
TIME "time"
DT "DT"
COMMENT {SLASH}{SLASH}+.*

WORD {ALPHA}+
RUNTIME "Run time"

BLPROFVERSION "BLProfVersion"
PHFNAME "phFName"
BLPROFPROC "BLProfProc"
BLPROFDATAFILENAME "bl_prof_D_"
CALCENDTIME "calcEndTime"

COMMPROFVERSION "CommProfVersion"
CALLSTATSPROFVERSION "CallStatsVersion"
NPROCS "NProcs"
COMMSTATSSIZE "CommStatsSize"
NOUTFILES "NOutFiles"
CPDATAPROC "CommProfProc"
NCOMMSTATS "nCommStats"
REGIONNAME "RegionName"
CALLSTATSPROC "CallStatsProc"
FNAME "fName"
TIMEMINMAX "timeMinMax"
NRSS "nRSS"
NTRACESTATS "nTraceStats"
DATAFILE "datafile"
COMMDATAFILENAME "bl_comm_prof_D_"
SEEKPOS "seekpos"
BARRIERNUMBER "bNum"
NAMETAG "nTag"
NAMETAGNAMES "nameTagNames"
REDUCTION "red"
TIMERTIME "timerTime"
TAGRANGE "tagRange"
FINESTLEVEL "FinestLevel"
MAXLEVEL "MaxLevel"
REFRATIO "RefRatio"
PROBDOMAIN "ProbDomain"
COMPUTE "compute"
SERVICE "service"
HEADERFILE "HeaderFile"
COMMHEADERFILENAME "bl_comm_prof_H_"

CSTATSHEADERFILENAME "bl_call_stats_H_"
CSTATSDATAFILENAME "bl_call_stats_D_"
INCLUDEALL "__IncludeAll__"
INCLUDENONE "__IncludeNone__"

CPU "cpu"
SLOT "slot"
CAGE "cage"
CABINET "cabinet"
CAB_POSITION "cab_position"
CAB_ROW "cab_row"
X_COORD "x_coord"
Y_COORD "y_coord"
Z_COORD "z_coord"
PROCESS_SLOTS "process_slots"
PROCESS_SLOTS_FREE "process_slots_free"
PROCESSOR_STATUS_UP "processor_status='up'"
PROCESSOR_STATUS_DOWN "processor_status='down'"
PROCESSOR_TYPE_SERVICE "processor_type='service'"
PROCESSOR_TYPE_COMPUTE "processor_type='compute'"
ALLOC_MODE_BATCH "alloc_mode='batch'"
ALLOC_MODE_OTHER "alloc_mode='other'"

PROCESSOR_ID "processor_id"
OD_ALLOCATOR_ID "od_allocator_id"
NEXT_RED_BLACK_SWITCH "next_red_black_switch"
PROCESSOR_SPEC "processor_spec"
SNULL "null"


QNAME {DQUOTE}.+{DQUOTE}

NAME "name"
INDEX "index"

STEP "STEP"
TIMEGL "TIME"
REGRID "REGRID"
WLB  "with lbase"
LEVEL "Level"
GRIDS "grids"
CELLS "cells"
PCTOD "% of domain"


PARENS ({LPAREN}{RPAREN})
DOTS {DOT}{DOT}+


%%


<<EOF>>     { /* stop here <<EOF>> */ BEGIN(INITIAL); return(ENDOFFILE); }

%{
 /* INITIAL */
%}
<INITIAL>{

{WS}       { /* ignore whitespace */ }
{MPITOKEN} { return(MPITOKEN);   }
{OMP}      { return(OMP);   }
{COMMENT}  { return(COMMENT);   }
{DT}       { return(DT);   }

{BLPROFVERSION}    { return(BLPROFVERSION);   }
{PHFNAME}          { return(PHFNAME);   }
{BLPROFPROC}       { return(BLPROFPROC);   }
{BLPROFDATAFILENAME}{DIGIT}+     { yylval.cValue = strdup(yytext); return(BLPROFDATAFILENAME); }
{CALCENDTIME}      { return(CALCENDTIME);   }

{COMMPROFVERSION}  { return(COMMPROFVERSION);   }
{NPROCS}           { return(NPROCS);   }
{COMMSTATSSIZE}    { return(COMMSTATSSIZE);   }
{NOUTFILES}        { return(NOUTFILES);   }
{FINESTLEVEL}      { return(FINESTLEVEL);   }
{MAXLEVEL}         { return(MAXLEVEL);   }
{REFRATIO}         { return(REFRATIO);   }
{PROBDOMAIN}       { return(PROBDOMAIN);   }
{COMPUTE}          { return(COMPUTE);   }
{SERVICE}          { return(SERVICE);   }
{HEADERFILE}       { return(HEADERFILE);   }
{CPDATAPROC}       { return(CPDATAPROC);   }
{NCOMMSTATS}       { return(NCOMMSTATS);   }
{DATAFILE}         { return(DATAFILE);   }
{SEEKPOS}          { return(SEEKPOS);   }
{BARRIERNUMBER}    { return(BARRIERNUMBER);   }
{NAMETAG}          { return(NAMETAG);   }
{NAMETAGNAMES}     { return(NAMETAGNAMES);   }
{REDUCTION}        { return(REDUCTION);   }
{TIMEMINMAX}       { return(TIMEMINMAX);   }
{TIMERTIME}        { return(TIMERTIME);   }
{TAGRANGE}         { return(TAGRANGE);   }
{NAME}             { return(NAME);   }
{INDEX}            { return(INDEX);   }
{LPAREN}           { return(LPAREN);   }
{RPAREN}           { return(RPAREN);   }
{QNAME}            { yylval.cValue = strdup(yytext); return(QNAME);   }
{COMMDATAFILENAME}{DIGIT}+     { yylval.cValue = strdup(yytext); return(COMMDATAFILENAME);   }
{COMMHEADERFILENAME}{DIGIT}+   { yylval.cValue = strdup(yytext); return(COMMHEADERFILENAME);   }

{CALLSTATSPROFVERSION}  { return(CALLSTATSPROFVERSION);   }
{FNAME}            { return(FNAME);   }
{REGIONNAME}       { return(REGIONNAME);   }
{CALLSTATSPROC}    { return(CALLSTATSPROC);   }
{NRSS}                 { return(NRSS);   }
{NTRACESTATS}      { return(NTRACESTATS);   }
{INCLUDEALL}       { return(INCLUDEALL);   }
{INCLUDENONE}      { return(INCLUDENONE);   }
{CSTATSDATAFILENAME}{DIGIT}+     { yylval.cValue = strdup(yytext); return(CSTATSDATAFILENAME);   }
{CSTATSHEADERFILENAME}{DIGIT}+   { yylval.cValue = strdup(yytext); return(CSTATSHEADERFILENAME);   }

{SIGN}{DIGIT}+                    { yylval.iValue = atol(yytext); return(NUMINTEGER); }
{SIGN}{DECIMALNUMBER}{EXPONENT}?  { yylval.fValue = atof(yytext); return(NUMFLOAT);   }
{SIGN}{DIGIT}+{EXPONENT}          { yylval.fValue = atof(yytext); return(NUMFLOAT);   }


{STEP}       { return(STEP);   }
{TIMEGL}     { return(TIMEGL); }
{REGRID}     { return(REGRID); }
{WLB}        { return(WLB);    }
{LEVEL}      { return(LEVEL);  }
{GRIDS}      { return(GRIDS);  }
{CELLS}      { return(CELLS);  }
{PCTOD}      { return(PCTOD);  }


{RUNTIME}            { return RUNTIME; }

{TIME}{COLON}?   { return(TIME);   }

{POUNDCOMMENT}           { return(POUNDCOMMENT);   }
{CPU}                    { return(CPU);  }
{SLOT}                   { return(SLOT); }
{CAGE}                   { return(CAGE); }
{CABINET}                { return(CABINET); }
{CAB_POSITION}           { return(CAB_POSITION); }
{CAB_ROW}                { return(CAB_ROW); }
{X_COORD}                { return(X_COORD); }
{Y_COORD}                { return(Y_COORD); }
{Z_COORD}                { return(Z_COORD); }
{PROCESS_SLOTS}          { return(PROCESS_SLOTS); }
{PROCESS_SLOTS_FREE}     { return(PROCESS_SLOTS_FREE); }
{PROCESSOR_STATUS_UP}    { return(PROCESSOR_STATUS_UP); }
{PROCESSOR_STATUS_DOWN}  { return(PROCESSOR_STATUS_DOWN); }
{PROCESSOR_TYPE_SERVICE} { return(PROCESSOR_TYPE_SERVICE); }
{PROCESSOR_TYPE_COMPUTE} { return(PROCESSOR_TYPE_COMPUTE); }
{ALLOC_MODE_BATCH}       { return(ALLOC_MODE_BATCH); }
{ALLOC_MODE_OTHER}       { return(ALLOC_MODE_OTHER); }
{PROCESSOR_ID}           { return(PROCESSOR_ID); }
{OD_ALLOCATOR_ID}        { return(OD_ALLOCATOR_ID); }
{NEXT_RED_BLACK_SWITCH}  { return(NEXT_RED_BLACK_SWITCH); }
{PROCESSOR_SPEC}         { return(PROCESSOR_SPEC); }
{SNULL}                  { return(SNULL); }

{ENDL}     { return DELENDLINE; }
{PLUS}     { return PLUS;       }
{MINUS}    { return MINUS;      }
{EQUAL}    { return EQUAL; }
{COLON}    { return COLON; }
{LBRACKET} { return LBRACKET; }
{RBRACKET} { return RBRACKET; }
{COMMA}    { return COMMA; }

{WORD}     { yylval.cValue = strdup(yytext);  return(WORD);   }
{EOS}      { /* stop here EOS */       BEGIN(INITIAL); return(ENDOFFILE);  }
{DOTS}     { return DOTS; }
.          { /* yylval.cValue = strdup(yytext); */  return UNKNOWN; }


}

%{
 /* end INITIAL  */
%}

%%


/* need this because we cannot directly access elements of YY_BUFFER_STATE externally */
char *GetYYBufPos(YY_BUFFER_STATE ybs) {
  return ybs->yy_buf_pos;
}

