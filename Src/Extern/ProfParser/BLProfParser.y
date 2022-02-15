/* ---------------------------------------------------------- */
/*   BLProfParser.y                                           */
/* ---------------------------------------------------------- */
%{
#include <stdio.h>
#include <AMReX_BLProfStats.H>
#include <AMReX_CommProfStats.H>
#include <AMReX_RegionsProfStats.H>


#define theBLPptr ((BLProfStats *) outpptr)

#define CCCOMMENT SLASH(/)
#define SLASH(s) /##s

#define vout0 if(theBLPptr->Verbose() >= 0) cout
#define vout1 if(theBLPptr->Verbose() >= 1) cout
#define vout2 if(theBLPptr->Verbose() >= 2) cout

int yyerror(void *outpptr, const char *s);
extern int yylex();
extern int yylineno;
extern char *yytext;

static int lineCount = 0;

%}

%parse-param {void *outpptr}

%code requires {

#include <iostream>
#include <vector>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
using std::ostream;

  struct IVec {
    int iv[BL_SPACEDIM];
    int &operator[](int ii) { return iv[ii]; }
    int get(int ii) { return iv[ii]; }
    std::ostream &operator<<(std::ostream &os) {
      os << '(' AMREX_D_TERM(<< iv[0], << ',' << iv[1], << ',' << iv[2]) << ')';
      return os;
    }
  };

void CPIV(amrex::IntVect &bliv, IVec &ivec);

}

%union {
  long   iValue;
  double fValue;
  char  *cValue;
  IVec   ivValue;
  std::vector<int> *vintptr;
}

/* declare tokens */
%token <iValue> NUMINTEGER
%token <fValue> NUMFLOAT

%token <cValue> WS ENDL ENDOFFILE
%token <cValue> DIGIT ALPHA PLUS MINUS SIGN EXPONENT
%token <cValue> COMMA EQUAL COLON EOS LPAREN RPAREN SLASH LINE
%token <cValue> LBRACKET RBRACKET DOTS
%token <cValue> PLUSSIGN MINUSSIGN DQUOTE SQUOTE DELCOMMA DELENDLINE
%token <cValue> MPITOKEN OMP WORD UNKNOWN TIME
%token <cValue> DT COMMENT RUNTIME

%token <cValue> BLPROFVERSION PHFNAME BLPROFPROC BLPROFDATAFILENAME
%token <cValue> CALCENDTIME

%token <cValue> COMMPROFVERSION NPROCS COMMSTATSSIZE CPDATAPROC NCOMMSTATS
%token <cValue> DATAFILE COMMDATAFILENAME SEEKPOS PROCNAME
%token <cValue> BARRIERNUMBER NAME QNAME INDEX NAMETAG NAMETAGNAMES
%token <cValue> TIMEMINMAX TIMERTIME REDUCTION TAGRANGE
%token <cValue> STEP TIMEGL REGRID WLB LEVEL GRIDS CELLS PCTOD
%token <cValue> FINESTLEVEL MAXLEVEL REFRATIO PROBDOMAIN COMPUTE SERVICE
%token <cValue> NOUTFILES HEADERFILE COMMHEADERFILENAME

%token <cValue> CALLSTATSPROFVERSION CSTATSHEADERFILENAME CSTATSDATAFILENAME
%token <cValue> REGIONNAME CALLSTATSPROC FNAME INCLUDEALL INCLUDENONE
%token <cValue> NRSS NTRACESTATS
%token <cValue> POUND POUNDCOMMENT CPU SLOT CAGE CABINET CAB_POSITION CAB_ROW
%token <cValue> X_COORD Y_COORD Z_COORD PROCESS_SLOTS PROCESS_SLOTS_FREE
%token <cValue> PROCESSOR_STATUS_UP PROCESSOR_STATUS_DOWN
%token <cValue> PROCESSOR_TYPE_SERVICE PROCESSOR_TYPE_COMPUTE
%token <cValue> ALLOC_MODE_BATCH ALLOC_MODE_OTHER
%token <cValue> PROCESSOR_ID OD_ALLOCATOR_ID NEXT_RED_BLACK_SWITCH PROCESSOR_SPEC
%token <cValue> SNULL


%type <cValue> infile line lines word words error serviceorcompute allocmode

%type <ivValue> intvect

%token <IVec> iv3d

%%

infile: { /* empty */
          vout0 << "infile0::lineCount = " << lineCount << endl;
          return 0;
        }
  |     lines {
          vout0 << "infile1::lineCount = " << lineCount << endl;
          return 0;
        }
;


lines:  line {
          ++lineCount;
        }
  |     lines line {
          ++lineCount;
        }
;


line: DELENDLINE {
        vout2 << "aDELENDLINE" << endl;
      }
  |   knownline DELENDLINE {
        vout2 << "kDELENDLINE" << endl;
      }
  |   words DELENDLINE {
        vout2 << "wDELENDLINE" << endl;
      }
  |   ENDOFFILE {
        vout2 << "infile2::lineCount = " << lineCount << endl;
        vout0 << "ENDOFFILE:  lineCount = " << lineCount << endl;
        return 0;
      }
  |   error DELENDLINE {
        vout2 << "eDELENDLINE:  " << endl;
      }
;


knownline: MPITOKEN WORD WORD NUMINTEGER MPITOKEN WORD {
             vout1 << "kMPI:  " << $4 << endl;
           }

  |        MPITOKEN WORD WORD NUMINTEGER WORD {
             vout1 << "kMPI:  " << $4 << endl;
           }

  |        words MPITOKEN WORD WORD NUMINTEGER words {
             vout1 << "kMPI:  " << $5 << endl;
           }

  |        OMP WORD WORD NUMINTEGER OMP WORD {
             vout1 << "kOMP:  " << $4 << endl;
           }

  |        RUNTIME EQUAL NUMFLOAT {
             vout0 << "RUNTIME:  " << $3 << endl;
           }

  |        BLPROFVERSION NUMINTEGER {
             vout0 << "BLPROFVERSION = " << $2 << endl;
	     theBLPptr->SetBLPVersion($2);
           }

  |        PHFNAME QNAME
           {
             vout0 << "PHFNAME = " << $2 << endl;
             theBLPptr->AddFunctionName($2);
             free($2);
           }

  |        BLPROFPROC NUMINTEGER
	   DATAFILE BLPROFDATAFILENAME
	   SEEKPOS NUMINTEGER
	   {
             vout0 << "BLPROFPROC = " << $2 << endl;
             vout0 << "DATAFILE   = " << $4 << endl;
             vout0 << "SEEKPOS    = " << $6 << endl;
	     theBLPptr->InitBLProfDataBlock($2, $4, $6);
	     free($4);
           }

  |        CALCENDTIME NUMFLOAT {
             vout0 << "CALCENDTIME = " << $2 << endl;
	     theBLPptr->AddCalcEndTime($2);
           }

  |        COMMPROFVERSION NUMINTEGER {
             vout0 << "COMMPROFVERSION = " << $2 << endl;
	     theBLPptr->SetCPVersion($2);
           }

  |        NPROCS NUMINTEGER {
             vout0 << "NPROCS = " << $2 << endl;
	     theBLPptr->SetNProcs($2);
           }

  |        NOUTFILES NUMINTEGER {
             vout0 << "NOUTFILES = " << $2 << endl;
	     theBLPptr->SetNOutFiles($2);
           }

  |        COMMSTATSSIZE NUMINTEGER {
             vout0 << "COMMSTATSSIZE = " << $2 << endl;
	     theBLPptr->SetCSSize($2);
           }

  |        CPDATAPROC NUMINTEGER
           NCOMMSTATS NUMINTEGER
	   DATAFILE COMMDATAFILENAME
	   SEEKPOS NUMINTEGER
	   {
             vout0 << "CPDATAPROC = " << $2 << endl;
             vout0 << "NCOMMSTATS = " << $4 << endl;
             vout0 << "DATAFILE   = " << $6 << endl;
             vout0 << "SEEKPOS    = " << $8 << endl;
	     theBLPptr->InitCommDataBlock($2, $4, $6, $8);
	     free($6);
           }

  |        CPDATAPROC NUMINTEGER
           NCOMMSTATS NUMINTEGER
	   DATAFILE COMMDATAFILENAME
	   SEEKPOS NUMINTEGER
	   WORD
	   {
             vout0 << "CPDATAPROC = " << $2 << endl;
             vout0 << "NCOMMSTATS = " << $4 << endl;
             vout0 << "DATAFILE   = " << $6 << endl;
             vout0 << "SEEKPOS    = " << $8 << endl;
             vout0 << "PROCNAME   = " << $9 << endl;
	     theBLPptr->InitCommDataBlock($2, $4, $6, $8, $9);
	     free($6);
	     free($9);
           }

  |        CPDATAPROC NUMINTEGER
           NCOMMSTATS NUMINTEGER
	   DATAFILE COMMDATAFILENAME
	   SEEKPOS NUMINTEGER
	   WORD NUMINTEGER
	   {
             vout0 << "CPDATAPROC   = " << $2 << endl;
             vout0 << "NCOMMSTATS   = " << $4 << endl;
             vout0 << "COMMDATAFILE = " << $6 << endl;
             vout0 << "SEEKPOS      = " << $8 << endl;
             vout0 << "NODENAME     = " << $9 << endl;
             vout0 << "NODENUMBER   = " << $10 << endl;
	     theBLPptr->InitCommDataBlock($2, $4, $6, $8, $9, $10);
	     free($6);
           }

  |        HEADERFILE COMMHEADERFILENAME
	   {
             vout0 << "HEADERFILE   = " << $2 << endl;
	     theBLPptr->AddCommHeaderFileName($2);
	     free($2);
           }

  |        BARRIERNUMBER NUMINTEGER QNAME NUMINTEGER {
             vout0 << "BARRIERNUMBER = " << $2 << "  NAME = " << $3 << "  INDEX = " << $4 << endl;
	     theBLPptr->AddBarrier($2, $3, $4);
	     free($3);
           }

  |        NAMETAG NUMINTEGER NUMINTEGER {
             vout0 << "NAMETAGINDEX = " << $2 << "  INDEX = " << $3 << endl;
	     theBLPptr->AddNameTag($2, $3);
           }

  |        NAMETAGNAMES QNAME {
             vout0 << "NAMETAGNAMES = " << $2 << endl;
	     theBLPptr->AddNameTagName($2);
	     free($2);
           }

  |        REDUCTION NUMINTEGER NUMINTEGER {
             vout0 << "REDUCTION = " << $2 << "  INDEX = " << $3 << endl;
	     theBLPptr->AddReduction($2, $3);
           }

  |        TIMEMINMAX NUMFLOAT NUMFLOAT {
             vout0 << "TIMEMINMAX = " << $2 << "  " << $3 << endl;
	     theBLPptr->AddTimeMinMax($2, $3);
           }

  |        TIMERTIME NUMFLOAT {
             vout0 << "TIMERTIME = " << $2 << endl;
	     theBLPptr->AddTimerTime($2);
           }
  |        STEP EQUAL NUMINTEGER TIMEGL EQUAL NUMFLOAT COLON REGRID WLB EQUAL NUMINTEGER {
             vout0 << "STEP = " << $3 << "  TIME = " << $6 << "  WLB = " << $11 << endl;
	     //theBLPptr->AddTimerTime($2);
           }

  |        STEP EQUAL NUMINTEGER TIMEGL EQUAL NUMINTEGER COLON REGRID WLB EQUAL NUMINTEGER {
             vout0 << "STEP = " << $3 << "  TIME = " << $6 << "  WLB = " << $11 << endl;
	     //theBLPptr->AddTimerTime($2);
           }

  |        LEVEL NUMINTEGER NUMINTEGER GRIDS NUMINTEGER CELLS NUMFLOAT PCTOD {
             vout0 << "L = " << $2 << "  G = " << $3 << "  C = " << $5 << "  %D = " << $7 << endl;
	     theBLPptr->AddGridLevel($2, $3);
           }

  |        LEVEL NUMINTEGER NUMINTEGER GRIDS NUMINTEGER CELLS NUMINTEGER PCTOD {
             vout0 << "L = " << $2 << "  G = " << $3 << "  C = " << $5 << "  %D = " << $7 << endl;
	     theBLPptr->AddGridLevel($2, $3);
           }

  |        TIMEGL EQUAL NUMFLOAT COLON REGRID WLB EQUAL NUMINTEGER {
             vout0 << "TIME = " << $3 << "  WLB = " << $8 << endl;
	     //theBLPptr->AddTimerTime($2);
           }

  |        TIMEGL EQUAL NUMINTEGER COLON REGRID WLB EQUAL NUMINTEGER {
             vout0 << "TIME = " << $3 << "  WLB = " << $8 << endl;
	     //theBLPptr->AddTimerTime($2);
           }

  |        NUMINTEGER COLON
           LPAREN
             LPAREN NUMINTEGER COMMA NUMINTEGER COMMA NUMINTEGER RPAREN
             LPAREN NUMINTEGER COMMA NUMINTEGER COMMA NUMINTEGER RPAREN
             LPAREN NUMINTEGER COMMA NUMINTEGER COMMA NUMINTEGER RPAREN
	   RPAREN
	   NUMINTEGER NUMINTEGER NUMINTEGER
	   COLON COLON
	   NUMINTEGER
           {
             vout0 << "GLLEV = " << $1;
	     vout0 << "  Box = ((" << $5 << "," << $7 << "," << $9 << ") (";
	     vout0 << $12 << "," << $14 << "," << $16 << ") (";
	     vout0 << $19 << "," << $21 << "," << $23 << "))";
	     vout0 << "  len = " << $26 << " " << $27 << " " << $28;
	     vout0 << "  proc = " << $31;
	     vout0 << endl;
	     theBLPptr->AddGrid3D($1,                  // level
	                          $5, $7, $9,          // loend
	                          $12 , $14 , $16,     // hiend
				  $19 , $21 , $23,     // centering
				  $26 , $27 , $28,     // npoints
				  $31);                // proc
           }

  |        TAGRANGE NUMINTEGER NUMINTEGER {
             vout0 << "TAGRANGE = " << $2 << "  " << $3 << endl;
	     theBLPptr->AddTagRange($2, $3);
           }

  |        FINESTLEVEL NUMINTEGER {
             vout0 << "FINESTLEVEL = " << $2 << endl;
	     theBLPptr->AddFinestLevel($2);
           }

  |        MAXLEVEL NUMINTEGER {
             vout0 << "MAXLEVEL = " << $2 << endl;
	     theBLPptr->AddMaxLevel($2);
           }

  |        REFRATIO NUMINTEGER intvect {
#if (BL_SPACEDIM == 2)
#else
             vout0 << "REFRATIO[lev] = " << $2 << " >> " << $3[0] << " << " << endl;
	     amrex::IntVect iv;
             CPIV(iv, $3);
	     theBLPptr->AddRefRatio($2, iv);
#endif
           }

  |        PROBDOMAIN NUMINTEGER LPAREN intvect intvect intvect RPAREN {
#if (BL_SPACEDIM == 2)
	     //cout << "PROBDOMAIN (spacedim = 2)" << endl;
#else
	     amrex::IntVect ivlo, ivhi, ivc;
             CPIV(ivlo, $4);
             CPIV(ivhi, $5);
             CPIV(ivc,  $6);
	     amrex::Box pd(ivlo, ivhi, ivc);
             //vout0 << "PROBDOMAIN[" << $2 << "] = " << $4 << " " << $5
	           //<< "  " << $6 << endl;
	     theBLPptr->AddProbDomain($2, pd);
#endif
           }

  |        NUMINTEGER COMPUTE
           NUMINTEGER NUMINTEGER NUMINTEGER
	   NUMINTEGER NUMINTEGER
	   NUMINTEGER NUMINTEGER NUMINTEGER
           {
             cout << $1 << " COMPUTE " << $7 << " " << $8 << " " << $9 << " " << $10 << endl;
	     theBLPptr->AddTopoCoord($1, $7, $8, $9, $10);
           }

  |        NUMINTEGER SERVICE
           NUMINTEGER NUMINTEGER NUMINTEGER
	   NUMINTEGER NUMINTEGER
	   NUMINTEGER NUMINTEGER NUMINTEGER
           {
             cout << $1 << " SERVICE " << $7 << " " << $8 << " " << $9 << " " << $10 << endl;
	     theBLPptr->AddTopoCoord($1, $7, $8, $9, $10, true);
           }

  |        CALLSTATSPROFVERSION NUMINTEGER
           {
             vout1 << "CALLSTATSPROFVERSION = " << $2 << endl;
             theBLPptr->SetCSVersion($2);
           }

  |        REGIONNAME QNAME NUMINTEGER
           {
             vout1 << "REGIONNAME = " << $2 << endl;
             vout1 << "REGIONNUMBER = " << $3 << endl;
             theBLPptr->AddRegionName($2, $3);
             free($2);
           }


  |        HEADERFILE CSTATSHEADERFILENAME
           {
             vout1 << "CSTATSHEADERFILENAME   = " << $2 << endl;
             theBLPptr->AddCStatsHeaderFileName($2);
             free($2);
           }

  |        CALLSTATSPROC NUMINTEGER
           NRSS NUMINTEGER
           NTRACESTATS NUMINTEGER
           DATAFILE CSTATSDATAFILENAME
           SEEKPOS NUMINTEGER
           {
             vout1 << "CALLSTATSPROC = " << $2 << endl;
             vout1 << "NRSS                 = " << $4 << endl;
             vout1 << "NTRACESTATS          = " << $6 << endl;
             vout1 << "DATAFILE             = " << $8 << endl;
             vout1 << "SEEKPOS              = " << $10 << endl;
             theBLPptr->InitCStatsDataBlock($2, $4, $6, $8, $10);
             free($8);
           }

  |        FNAME QNAME NUMINTEGER
           {
             vout1 << "FNAME = " << $2 << "  ";
             vout1 << "FNUMBER = " << $3 << endl;
             theBLPptr->AddFunctionName($2, $3);
             free($2);
           }

  |        PLUS QNAME NUMINTEGER
           {
             vout1 << "PLUS = " << $2 << "  " << $3 << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::ON, $2, $3);
             free($2);
           }

  |        MINUS QNAME NUMINTEGER
           {
             vout1 << "MINUS = " << $2 << "  " << $3 << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::OFF, $2, $3);
             free($2);
           }

  |        COLON QNAME NUMINTEGER
           {
             vout1 << "COLON = " << $2 << "  " << $3 << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::UNDEFINED, $2, $3);
             free($2);
           }

  |        INCLUDEALL
           {
             vout1 << "INCLUDEALL" << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::INCLUDEALL);
           }

  |        INCLUDENONE
           {
             vout1 << "INCLUDENONE" << endl;
             theBLPptr->SetFilter(RegionsProfStats::FilterStatus::INCLUDENONE);
           }

  |        POUNDCOMMENT
           {
             vout0 << "POUNDCOMMENT :: " << endl;
           }

  |        CPU                   EQUAL NUMINTEGER COMMA
           SLOT                  EQUAL NUMINTEGER COMMA
           CAGE                  EQUAL NUMINTEGER COMMA
           CABINET               EQUAL SNULL      COMMA
           CAB_POSITION          EQUAL NUMINTEGER COMMA
           CAB_ROW               EQUAL NUMINTEGER COMMA
           X_COORD               EQUAL NUMINTEGER COMMA
           Y_COORD               EQUAL NUMINTEGER COMMA
           Z_COORD               EQUAL NUMINTEGER COMMA
           PROCESS_SLOTS         EQUAL NUMINTEGER COMMA
           PROCESS_SLOTS_FREE    EQUAL NUMINTEGER COMMA
           processorstatus                        COMMA
           serviceorcompute                       COMMA
           allocmode                              COMMA
           PROCESSOR_ID          EQUAL NUMINTEGER COMMA
           OD_ALLOCATOR_ID       EQUAL NUMINTEGER COMMA
           NEXT_RED_BLACK_SWITCH EQUAL SNULL      COMMA
           PROCESSOR_SPEC        EQUAL SNULL
           {
	     vout0 << "_here 0:" << endl;
             vout0 << "CPU = " << $3 << endl;
             vout0 << "SLOT = " << $7 << endl;
             vout0 << "CAGE = " << $11 << endl;
             vout0 << "CAB_POSITION = " << $19 << endl;
             vout0 << "CAB_ROW = " << $23 << endl;
             vout0 << "X_COORD = " << $27 << endl;
             vout0 << "Y_COORD = " << $31 << endl;
             vout0 << "Z_COORD = " << $35 << endl;
             vout0 << "PROCESS_SLOT = " << $39 << endl;
             vout0 << "PROCESS_SLOT_FREE = " << $43 << endl;
             vout0 << "PROCESSOR_ID = " << $53 << endl;
             vout0 << "OD_ALLOCATOR_ID = " << $57 << endl;
	     vout0 << endl;
             theBLPptr->AddEdisonPID($27, $31, $35, $19, $23, $11,
	                             $7, $3, $53);
           }



;


serviceorcompute: PROCESSOR_TYPE_SERVICE
                  {
                    vout1 << "PROCESSOR_TYPE_SERVICE" << endl;
		  }
  |               PROCESSOR_TYPE_COMPUTE
                  {
                    vout1 << "PROCESSOR_TYPE_COMPUTE" << endl;
		  }
;

allocmode: ALLOC_MODE_BATCH
                  {
                    vout1 << "ALLOC_MODE_BATCH" << endl;
		  }
  |               ALLOC_MODE_OTHER
                  {
                    vout1 << "ALLOC_MODE_OTHER" << endl;
		  }
;

processorstatus: PROCESSOR_STATUS_UP
                  {
                    vout1 << "PROCESSOR_STATUS_UP" << endl;
		  }
  |               PROCESSOR_STATUS_DOWN
                  {
                    vout1 << "PROCESSOR_STATUS_DOWN" << endl;
		  }
;


intvect: LPAREN NUMINTEGER COMMA NUMINTEGER RPAREN
         {
	 }
  |      LPAREN NUMINTEGER COMMA NUMINTEGER COMMA NUMINTEGER RPAREN
         {
#if (BL_SPACEDIM == 2)
#else
	   $$[0] = $2;
	   $$[1] = $4;
	   $$[2] = $6;
	   amrex::IntVect iv;
	   iv[0] = $2;
	   iv[1] = $4;
	   iv[2] = $6;
	   //vout0 << "++++ $$ = (" << $2 << "," << $4 << "," << $6 << ")" << endl;
#endif
	 }
;


/*
box: LPAREN intvect intvect RPAREN
	 {
	 }
  |      LPAREN intvect intvect intvect RPAREN
         {
	   vout0 << "(" << $2 << " " << $3 << " " << $4 << ")";
	 }
;
*/


words: word { }
  | words word { }
;


word: WORD {
        vout2<< "wWORD = " << $1 << endl;
	free($1);
      }
  |   NUMINTEGER {
        vout2 << "wNUMINTEGER = " << $1 << endl;
      }
  |   NUMFLOAT {
        vout2 << "wNUMFLOAT = " <<  $1 << endl;
      }
  |   UNKNOWN {
        vout2 << "wUNKNOWN" << endl;
      }
  |   EQUAL {
        vout2 << "wEQUAL" << endl;
      }
  |   COLON {
        vout2 << "wCOLON" << endl;
      }
  |   COMMA {
        vout2 << "wCOMMA" << endl;
      }
  |   LBRACKET {
        vout2 << "wLBRACKET" << endl;
      }
  |   RBRACKET {
        vout2 << "wRBRACKET" << endl;
      }
  |   TIME {
        vout2 << "wTIME" << endl;
      }
  |   DOTS {
        vout2 << "wDOTS" << endl;
      }
  |   DT {
        vout2 << "wDT" << endl;
      }
  |   COMMENT {
        vout1 << "COMMENT:   " << endl;
      }
;


%%

int yyerror(void * /*outpptr*/, const char * /*s*/) {
  cerr << "*** Unrecognized output at line " << yylineno << "  ::  " << yytext << endl;
  return -1;
}



void CPIV(amrex::IntVect &bliv, IVec &ivec) {
#if (BL_SPACEDIM == 2)
  amrex::ignore_unused(bliv, ivec);
#else
  bliv[0] = ivec[0];
  bliv[1] = ivec[1];
  bliv[2] = ivec[2];
#endif
}

