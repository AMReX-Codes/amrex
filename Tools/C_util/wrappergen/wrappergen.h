#ifndef _WRAPPERGEN_H_
#define _WRAPPERGEN_H_

#include "expandingList.h"

#define MAX_IDENT_LEN 256

typedef struct fn_def_struct {
  char *name;
  char **argTypePrefix, **argTypeSuffix;
  char **argNames;
  int  nargs;
  char *returnType;
  xpandList wrapperdefs;
    /* list of integer indices pointing to the wrappers used on this fn */
} fn_def;

#ifndef STRDUP
#define STRDUP(str) strcpy( (char *) malloc( strlen( str ) + 1 ), (str) )
#endif

#ifndef STR_RANGE_DUP
#define STR_RANGE_DUP( cpy, start, end ) { \
  (cpy) = (char *)malloc( (end) - (start) + 1 ); \
  strncpy( (cpy), (start), (end) - (start) ); \
  (cpy)[(end) - (start)] = '\0'; }
#endif

typedef struct variable_ {
  char *typePrefix, *typeSuffix;
  char *rqName;  /* requested name */
} variable;

typedef struct fileinfo_ {
  char *name, *str;    /* str is the string to parse */
  int filenum, lineno;
    /* lineno should be set to the line number in the file that str
       starts on */
} fileinfo;

typedef struct wrapperdef_ {
  char *nameEscape, *prefixCode, *suffixCode;
  /* prefix/suffix code - code to go before/after the call */
  /* if no {{callfn}}, suffixCode will be null */
  variable *vars;
  int nvars, prefixLineNo, suffixLineNo, firstLine;
  fileinfo finfo;
    /* finfo - set *name and filenum so we know who to blame
       if this wrapper is goofy */
    /* when this wrapper is being written out, fill in lineno and string
       and just pass &winfo[i].finfo */
} wrapperdef;


typedef struct wrapperinfo_ {
  xpandList wrapperDefs;
} wrapperinfo;


typedef struct replacement_ {
  char *from, *to;
} replacement;


typedef struct rpcinfo_ {
  fn_def *fn_list;
  xpandList rpc;
  int n_fn;
} rpcinfo;

extern int using_cpp;

/* I got this trick from the Tcl implementation */
#ifdef _ANSI_ARGS_
#undef _ANSI_ARGS_
#endif

#ifdef __STDC__
#define _ANSI_ARGS_(x) x
#else
#define _ANSI_ARGS_(x) ()
#endif

void WriteWrappers _ANSI_ARGS_(( FILE *outf, char **wrapperFiles,
				 int nwrapperFiles, fn_def *fn_list,
				 int n_fn ));

void ReadWrapperFile _ANSI_ARGS_(( FILE *outf, char *fileName, int filenum,
				   fn_def *fn_list, int n_fn,
				   wrapperinfo *winfo ));

char* ReadFileIntoString _ANSI_ARGS_(( const char* inf ));


void ProcessString _ANSI_ARGS_(( FILE *outf,
				 fileinfo *finfo,
				 rpcinfo *rinfo,
				 wrapperinfo *winfo ));

/* either substitute with for{each,all}fn (outf), or define a new
   wrapper with fn[all] (winfo) */
/* escStartLine is set to the first line of the start of the escape */
void ProcessEscape _ANSI_ARGS_(( FILE *outf,
				 fileinfo *finfo,
				 rpcinfo *rinfo,
				 wrapperinfo *winfo,
				 char **escBodyList, int escBodyLen,
				 char *escBody, int escStartLine ));

/* finfo->lineno set to the first line of what is to be read */
/* is returned set to whatever is after what was read */
/* escStartLine tells what line the escape starts on */
int ReadUntilMatch _ANSI_ARGS_(( fileinfo *finfo,
                                 char *start,
				 char *end, char **escbody,
				 int escStartLine ));

int ReadUntilEscape _ANSI_ARGS_(( fileinfo *finfo,
				  char **preceding, char ***escBodyList,
				  int *escBodyLen, char **escBodyLiteral,
				  int *escStartLine ));

int CountNewlines _ANSI_ARGS_(( char *start, char *end ));

int RegisterVarType _ANSI_ARGS_(( char *type, xpandList varTypes ));

/* makes a copy of the string, and freeing will be difficult */
void ListizeString _ANSI_ARGS_(( char *str, char ***list, int *len ));

int IsReservedName _ANSI_ARGS_(( char *name ));

void OutChar _ANSI_ARGS_(( int c, int where, void *outputForm ));

void DoForEach _ANSI_ARGS_(( FILE *outf,
			     fileinfo *finfo,
			     rpcinfo *rinfo,
			     char **argv,
			     int argc, char *escBody, int startLine,
			     char *body ));

void DoForAll _ANSI_ARGS_(( FILE *outf,
			    fileinfo *finfo,
			    rpcinfo *rinfo,
			    char **argv,
			    int argc, char *escBody, int startLine,
			    char *body ));

void DoFn _ANSI_ARGS_(( fileinfo *finfo,
		        rpcinfo *rinfo,
		        wrapperinfo *winfo,
		        char **argv,
			int argc,
                        char *body,
		        int startingLine ));

void DoFnAll _ANSI_ARGS_(( fileinfo *finfo,
			   rpcinfo *rinfo,
			   wrapperinfo *winfo,
			   char **argv,
			   int argc,
			   char *body,
			   int startingLine ));

void ReadFnDef _ANSI_ARGS_(( fileinfo *finfo,
			     rpcinfo *rinfo,
			     wrapperinfo *winfo,
			     char **argv,
			     int argc,
			     char *body,
			     int startingLine,
			     int allFn ));

void WriteFunctionCalls _ANSI_ARGS_(( FILE *outf,
				      fn_def *fn_list,
				      int n_fn,
				      wrapperinfo *winfo ));

int IsFnInList _ANSI_ARGS_(( char *fn, fn_def *fn_list, int n_fn ));

char ***CreateUniqueVarNames _ANSI_ARGS_(( wrapperdef *wrapperList,
					   int nwrappers ));

void ReadVardecl _ANSI_ARGS_(( fileinfo *finfo, int startLine,
			       char *body, wrapperinfo *winfo,
			       xpandList vars ));

int ReadVardeclBasetype _ANSI_ARGS_(( char *filename, int lineno,
				      char *body, char **basetype,
				      char **end ));

int ReadVardeclVarname _ANSI_ARGS_(( char **readPt, char **varPrefix,
				     char **varName, char **varSuffix ));

void CheckForHiddenArgs _ANSI_ARGS_(( fn_def *fn_list, int fn_num,
				      wrapperinfo *winfo, int wrapperNum ));

int IsUnique _ANSI_ARGS_(( char *str, wrapperdef* wrapperList, int nwrappers,
			   char ***others,
			   int wrapperNum ));

void PrintWrapperCode _ANSI_ARGS_(( FILE *outf,
				    fn_def *fn_list,
				    int n_fn,
				    wrapperinfo *winfo,
				    char ***varNames,
				    int fn_num,
				    int wrapperNumIdx ));

#endif
