#ifndef BL_PARMPARSE_H
#define BL_PARMPARSE_H

//
// $Id: ParmParse.H,v 1.20 2001-07-24 18:12:35 car Exp $
//

#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <Array.H>

//
// ParmParse class implements a simple database for the storage and
// retrieval of command-line and input-file arguments.  The entries are
// stored in a static table in (name,value_list) pairs.
//
// The format of the input file is a series of DEFINITIONS.
//
// A DEFINITION is of the form  <name> = <value> <value> ...
// The equal sign is important since the list of values can span multiple
// lines.
//
// Comments in an input file include all text from a '#' character to the
// end of the line.  Here is an example input file:
//
//   niter = 100                # niter is an integer
//   title = "Double Wammy"     # example of a string with spaces
//   cell_size = 0.5 0.75       # cell spacing in each dimension
//   plot.var = Density 1 10    # a list of values
//   plot.var = Energy  5 12    # another list of values
//   bigarray = 1 2 3 4 5 6 7 8 # first part of array
//              9 10 11 12      # continuation of bigarray
//   test = apple "boy blue" 10 20 30 40
//   FILE = prob_file           # insert contents of this "prob_file" here
//
// The "FILE = <filename>" definition is special.  Rather than just
// adding this entry to the database, it reads the contents of <filename>
// into the database.
//
// ParmParse stores all entries in a static table which is built the
// first time a ParmParse object is constructed (usually in main()).
// Subsequent invocations have access to this table.
// A ParmParse constructor has an optional "prefix" argument that will
// limit the searches to only those entries of the table with this prefix
// in name.  For example:
//     ParmParse pp("plot");
// will find only those entries with name given by "plot.<string>".
//
// All values in the table are stored as strings.  For example, the
// values of "cell_size" in the above input file are stored as the
// strings "0.5" and "0.75".  These strings can be returned as either
// string of numeric values by the query functions.
// Character strings with spaces must be delimited by double quotes
// in the input file but the quotes are stripped before they are entered
// into the table.  For example, 'title' in the above input file has a
// single value, the string 'Double Wammy' (without the quotes).
// Each value in the list associated with a definition can be referred to
// by its index number.  The index numbers start at 0 just like an array
// in the C programming language.  Consider the definition of "test" in
// the above input file.  The first value 'apple'is a string with index
// 0.  The second value 'boy blue' is a string with index 1.  The
// remaining four values are integers indexed 2, 3, 4, and 5.
//
// For a string value to represent an integer or float it must fit the
// following regular experssion:
//   Sign    ::= '+' | '-'
//   Digit   ::= '0' | '1' | ... | '9'
//   Integer ::= [Sign]Digit+
//   Exp     ::= ('e'|'E')Integer
//   Float   ::= ( Integer[.Digit*][Exp] | [Integer].Digit+[Exp] )
//
// Where '+' indicates one or more occurences, '*' represents zero or
// more occurences, '|' means one or the other and '[]' represents zero
// or one occurence.
//
// Note that floats and doubles have the same string representation and
// that the FORTRAN "double" exponent format is not supported.
// That is, 1.0d+3 is not a valid representation of a floating point
// number but that 1.0e+3 is acceptable.
//
// There are a host of functions allowing the user to query the database
// and retrieve values.  Here are some general rules about the names of
// the member functions:
//
// * Functions with the string "get" in their names attempt to get a
//   value or an array of values from the table.  They generate a
//   run-time error if they are not successful.
//
// * Functions with the string "query" in their names attempt to get a
//   value or an array of values from the table.  They return the value 1
//   (true) if they are successful and 0 (false) if not.
//
// * Functions with the string "kth" in their names get values from
//   the Kth entry with the given name in the database.  This is
//   necessary since there may be multiple definitions with the same
//   name in the database.
//
// * Functions without the string "kth" in their names get values from
//   the last entry with the given name in the database.  Note that the
//   definitions from the command line are appended to the database table
//   and hence will be the last entries.
//
// * Functions with the string "arr" in their names get an Array of
//   values from the given entry in the table.  The array argument is
//   resized (if necessary) to hold all the values requested.
//
// * Functions without the string "arr" in their names get single
//   values from the given entry in the table.
//
// The following is a code sample showing how to use ParmParse:
//
// main(int argc, char **argv)
// {
//     char* in_file_name = argv[1];
//     ParmParse::Initialize(argc-2, argv+2, in_file_name);
//
//     // Query table for value of "niter".  If not in table
//     // then set to default value
//     if (!pp.query("niter",niter)) niter = 20;
//
//     // read array of cell sizes if in table
//     Array<float> dx;
//     if (nx=pp.countval("cell_size")) {
//        // get nx values starting at index 0 and store in dx.
//        // dx is automatically resized here.
//        pp.getarr("cell_size",dx,0,nx);
//     }
//     ParmParse::Finalize();
// }
//
// void do_graphics()
// {
//    //
//    // Will only query entries with the "plot" prefix:
//    //
//    ParmParse pp("plot");
//    //
//    // Read all variables with "plot.var" keyword.
//    //
//    std::string var_name;
//    Array<int> range;
//    int num = pp.countname("var");
//    for (int k = 0; k < num; k++)
//    {
//       //
//       // Element 0 in list is a string.
//       //
//       pp.getkth("var",k,var_name,0);
//       //
//       // Elements 1 and 2 are integers.
//       // Note that "range" will be resized to hold 2 elements.
//       //
//       pp.getktharr("var",k,range,1,2);
//       cout << "variable = " << var_name << "lo, hi = ",
//            << range[0] << " " << range[1] << endl;
//    }
// }
// -----------------------------------------------------------------
// -----------------------  END COMMENTS ---------------------------
// -----------------------------------------------------------------


//@Man:
//@Memo: Parse Parameters From Command Line and Input Files
/*@Doc:

  The ParmParse class is used to interpret parameters passed in to a program
  from the command line and an arbitrary collection of input files.  The
  parameters are stored in static table that can be queried by any object
  of type ParmParse.  A parameter is a  "definition".  A definition is
  of the form "<name> = <value><value>...<value>".  It is stored in the table
  as a name, value-list pair.

  In the following example, niter is a definition with the single
  integer value 10; name is a definition with the string value "big
  code" and dx is a definition with the two floating point values 0.5
  and 0.75 and iv is an IntVect(5,4)

  prog niter = 10 name = "big code" dx = 0.5 0.75 iv=(5,4)

  The ParmParse class has two constructors.  The first is responsible for
  building the table and is usually called by the main routine of an
  application.  It has arguments for the command line argc and argv parameters,
  as well as an optional filename argument for reading definitions from an
  input file.  The table is built by reading the input file first (if it
  exists) with the command line arguments added to the end of the table.
  The order of a definition in the table is significant, so command line
  parameters can be used to override definitions in the input file.  A
  definition of the explicit form: FILE=<filename> is not added to the table
  but is a directive to include the named file at that point in the table.

  The second constructor is generally used by other classes in the code.  It
  permits access to the table via a large collection of query functions.
  Both constructors have an optional prefix argument that narrows the search
  to entries in the table with the same prefix.  For example, let PlanR be a
  ParmParse object with code prefix "ope".  PlanR.get("val",v) will look for
  an entry in the parameter list of the form: ope.val==<value>, and will
  reject all entries not starting with the correct code prefix.

  The query functions search the table for definition names that match a given
  string (and prefix) and return values from the corresponding value list.
  The values can be returned as ints, Array<int>s, floats, Array<float>s,
  doubles, Array<double>s, std::strings, or Array<aSring>s.  All values in the
  table are stored as std::string objects, but if an int, float, or double is
  requested, the translation is done automatically.  In the previous example,
  the value of niter could be returned as either an std::string, an int, a double,
  or a float.  The values of dx can be returned as std::strings, floats, or
  doubles, but the value of name can be returned only as an std::string.

  Comments in an input file include all text from a `\#' character to the
  end of the line.  Here is a sample input file:

     niter = 100

     title = "Double Wammy"

     cell\_size = 0.5 0.75

     plot.var = Density 1 10

     plot.var = Energy  5 12

     bigarray = 1 2 3 4 5 6 7 8

                9 10 11 12

     aBox    = ((0,0) (5,5))		

     test = apple "boy blue" 10 20 30 40

     FILE = prob\_file
*/

class Box;
class IntVect;

class ParmParse
{
public:

    enum { LAST = -1, FIRST = 0 };

    /*@ManDoc: Construct an additional ParmParse object sharing the same
               internal table as any other such objects in existence.  If
               prefix is specified, load this string as the code prefix
               for this particular ParmParse object.
    */
    
    explicit ParmParse (const char* prefix = 0);
    
    //
    //@ManDoc: Returns true if name is in table.
    //
    bool contains (const char* name) const;

    /*@ManDoc: Returns the number of values associated with nth occurence of
               name (prepended with the prefix) in the table.  n == -1 implies
               the last occurence.
    */
    int countval (const char* name,
                  int         n = LAST) const;

    /*@ManDoc: Returns the number of times the given name (prepended with
               prefix) appears in the table.
    */
    int countname (const char* name) const;

    /*@ManDoc: Returns the number of times the given name (prepended with
               prefix) appears in the table.
    */
    int countname (const std::string& name) const;
    //
    //@ManDoc: Write the contents of the table in ASCII to the ostream.
    //
    static void dumpTable (std::ostream& os);

    /*@ManDoc: Get the ival'th value of kth occurrence of the requested name.
               If successful, the value is converted to a bool and stored
               in reference ref.  If the kth occurrence does not exist or
               ival'th value does not exist, or if the printed representation
               of the value cannot be converted to an bool, an error message is
               output and the program halts.   Note that ival == 0 is the first
               value in the list.  ParmParse converte the value 'true', and non-zero
	       integers or floats to bool(true), and bool(false) for 'false'
	       or zero integer or float values.
    */
    void getkth (const char* name,
                 int         k,
                 bool&        ref,
                 int         ival = FIRST) const;
    //
    //@ManDoc: Same as getkth() but searches for the last occurrence of name.
    //
    void get (const char* name,
              bool&        ref,
              int         ival = FIRST) const;

    /*@ManDoc: Similar to getkth() but returns 0 if there is no kth occurrence
               of name.  If successful, it returns 1 and stores the value in
               ref.  If the kth occurrence exists, but ival'th value of that
               occurrence does not, or if there is a type mismatch, then the
               program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  bool&       ref,
                  int         ival = FIRST) const;
    //
    //@ManDoc: Same as querykth() but searches for the last occurrence of name.
    //
    int query (const char* name,
               bool&        ref,
               int         ival = FIRST) const;

    /*@ManDoc: Get the ival'th value of kth occurrence of the requested name.
               If successful, the value is converted to an int and stored
               in reference ref.  If the kth occurrence does not exist or
               ival'th value does not exist, or if the printed representation
               of the value cannot be converted to an int, an error message is
               output and the program halts.   Note that ival == 0 is the first
               value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 int&        ref,
                 int         ival = FIRST) const;

    //@ManDoc: Same as getkth() but searches for the last occurrence of name.
    //
    void get (const char* name,
              int&        ref,
              int         ival = FIRST) const;

    /*@ManDoc: Similar to getkth() but returns 0 if there is no kth occurrence
               of name.  If successful, it returns 1 and stores the value in
               ref.  If the kth occurrence exists, but ival'th value of that
               occurrence does not, or if there is a type mismatch, then the
               program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  int&        ref,
                  int         ival = FIRST) const;
    //
    //@ManDoc: Same as querykth() but searches for the last occurrence of name.
    //
    int query (const char* name,
               int&        ref,
               int         ival = FIRST) const;

    /*@ManDoc: Get the ival'th value of kth occurrence of the requested name.
               If successful, the value is converted to a float and stored
               in reference ref.  If the kth occurrence does not exist or
               ival'th value does not exist, or if the printed representation
               of the value cannot be converted to a float, an error message
               is output and the program halts.   Note that ival == 0 is the
               first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 float&      ref,
                 int         ival = FIRST) const;
    //
    //@ManDoc: Same as getkth() but searches for the last occurrence of name.
    //
    void get (const char* name,
              float&      ref,
              int         ival = FIRST) const;

    /*@ManDoc: Similar to getkth() but returns 0 if there is no kth occurrence
               of name.  If successful, it returns 1 and stores the value in
               ref.  If the kth occurrence exists, but ival'th value of that
               occurrence does not, or if there is a type mismatch, then the
               program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  float&      ref,
                  int         ival = FIRST) const;
    //
    //@ManDoc: Same as querykth() but searches for the last occurrence of name.
    //
    int query (const char* name,
               float&      ref,
               int         ival = FIRST) const;

    /*@ManDoc: Get the ival'th value of kth occurrence of the requested name.
               If successful, the value is converted to a double and stored
               in reference ref.  If the kth occurrence does not exist or
               ival'th value does not exist, or if the printed representation
               of the value cannot be converted to a double, an error message
               is output and the program halts.   Note that ival = 0 is the
               first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 double&     ref,
                 int         ival = FIRST) const;
    //
    //@ManDoc: Same as getkth() but searches for the last occurrence of name.
    //
    void get (const char* name,
              double&     ref,
              int         ival = FIRST) const;

    /*@ManDoc: Similar to getkth() but returns 0 if there is no kth occurrence
               of name.  If successful, it returns 1 and stores the value in
               ref.  If the kth occurrence exists, but ival'th value of that
               occurrence does not, or if there is a type mismatch, then the
               program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  double&     ref,
                  int         ival = FIRST) const;
    //
    //@ManDoc: Same as querykth() but searches for the last occurrence of name.
    //
    int query (const char* name,
               double&     ref,
               int         ival = FIRST) const;

    /*@ManDoc: Get the ival'th value of kth occurrence of the requested name.
               If successful, the value is converted to a std::string and stored
               in reference ref.  If the kth occurrence does not exist or
               ival'th value does not exist, or if the printed representation
               of the value cannot be converted to a std::string, an error message
               is output and the program halts.   Note that ival = 0 is the
               first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 std::string&    ref,
                 int         ival = FIRST) const;
    //
    //@ManDoc: Same as getkth() but searches for the last occurrence of name.
    //
    void get (const char* name,
              std::string&    ref,
              int         ival = FIRST) const;

    /*@ManDoc: Similar to getkth() but returns 0 if there is no kth occurrence
               of name.  If successful, it returns 1 and stores the value in
               ref.  If the kth occurrence exists, but ival'th value of that
               occurrence does not, or if there is a type mismatch, then the
               program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  std::string&   ref,
                  int         ival = FIRST) const;
    //
    //@ManDoc: Same as querykth() but searches for the last occurrence of name.
    //
    int query (const char* name,
               std::string&    ref,
               int         ival = FIRST) const;

    /*@ManDoc: Get the ival'th value of kth occurrence of the requested name.
               If successful, the value is converted to an IntVect and stored
               in reference ref.  If the kth occurrence does not exist or
               ival'th value does not exist, or if the printed representation
               of the value cannot be converted to a IntVect, an error message
               is output and the program halts.   Note that ival = 0 is the
               first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 IntVect&    ref,
                 int         ival = FIRST) const;
    //
    //@ManDoc: Same as getkth() but searches for the last occurrence of name.
    //
    void get (const char* name,
              IntVect&    ref,
              int         ival = FIRST) const;

    /*@ManDoc: Similar to getkth() but returns 0 if there is no kth occurrence
               of name.  If successful, it returns 1 and stores the value in
               ref.  If the kth occurrence exists, but ival'th value of that
               occurrence does not, or if there is a type mismatch, then the
               program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  IntVect&   ref,
                  int         ival = FIRST) const;
    //
    //@ManDoc: Same as querykth() but searches for the last occurrence of name.
    //
    int query (const char* name,
               IntVect&    ref,
               int         ival = FIRST) const;

    /*@ManDoc: Get the ival'th value of kth occurrence of the requested name.
               If successful, the value is converted to a Box and stored
               in reference ref.  If the kth occurrence does not exist or
               ival'th value does not exist, or if the printed representation
               of the value cannot be converted to a Box, an error message
               is output and the program halts.   Note that ival = 0 is the
               first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 Box&        ref,
                 int         ival = FIRST) const;
    //
    //@ManDoc: Same as getkth() but searches for the last occurrence of name.
    //
    void get (const char* name,
              Box&        ref,
              int         ival = FIRST) const;

    /*@ManDoc: Similar to getkth() but returns 0 if there is no kth occurrence
               of name.  If successful, it returns 1 and stores the value in
               ref.  If the kth occurrence exists, but ival'th value of that
               occurrence does not, or if there is a type mismatch, then the
               program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  Box&        ref,
                  int         ival = FIRST) const;
    //
    //@ManDoc: Same as querykth() but searches for the last occurrence of name.
    //
    int query (const char* name,
               Box&        ref,
               int         ival = FIRST) const;

    /*@ManDoc: Gets an Array<int> of num\_val values from kth occurrence of
               given name.  If successful, the values are converted to an int
               and stored in the Array<int> object ref.  ref is resized (if
               necessary) to hold num\_val values.  The value in the list
               indexed by start\_ix is copied into Array<int>[0], Array<int>[1]
               holds start\_ix+1, etc.  If the kth occurrence does not exist
               or there are fewer than start\_ix + num\_val values associated
               with the kth occurrence, or if some of the values cannot be
               converted to an int, an error message is reported and the
               program halts.
    */
    void getktharr (const char* name,
                    int         k,
                    Array<int>& ref,
                    int         start_ix,
                    int         num_val) const;
    //
    //@ManDoc: Same as getktharr() but searches for last occurrence of name.
    //
    void getarr (const char* name,
                 Array<int>& ref,
                 int         start_ix,
                 int         num_val) const;
    //
    //@ManDoc: queryktharr() is to querykth() as getktharr() is to getkth().
    //
    int queryktharr (const char* name,
                     int         k,
                     Array<int>& ref,
                     int         start_ix,
                     int         num_val) const;
    //
    //@ManDoc: Same as queryktharr() but searches for last occurrence of name.
    //
    int queryarr (const char* name,
                  Array<int>& ref,
                  int         start_ix,
                  int         num_val) const;

    /*@ManDoc: Gets an Array<float> of num\_val values from kth occurrence of
               given name.  If successful, the values are converted to a float
               and stored in the Array<float> object ref.  ref is resized (if
               necessary) to hold num\_val values.  The value in the list
               indexed by start\_ix is copied into Array<float>[0],
               Array<float>[1] holds start\_ix+1, etc.  If the kth occurrence
               does not exist or there are fewer than start\_ix + num\_val
               values associated with the kth occurrence, or if some of the
               values cannot be converted to a float, an error message is
               reported and the program halts.
    */
    void getktharr (const char*   name,
                    int           k,
                    Array<float>& ref,
                    int           start_ix,
                    int           num_val) const;
    //
    //@ManDoc: Same as getktharr() but searches for last occurrence of name.
    //
    void getarr (const char*   name,
                 Array<float>& ref,
                 int           start_ix,
                 int           num_val) const;
    //
    //@ManDoc: queryktharr() is to querykth() as getktharr() is to getkth().
    //
    int queryktharr (const char*   name,
                     int           k,
                     Array<float>& ref,
                     int           start_ix,
                     int           num_val) const;
    //
    //@ManDoc: Same as queryktharr() but searches for last occurrence of name.
    //
    int queryarr (const char*   name,
                  Array<float>& ref,
                  int           start_ix,
                  int           num_val) const;

    /*@ManDoc: Gets an Array<double> of num\_val values from kth occurrence of
               given name.  If successful, the values are converted to a double
               and stored in the Array<double> object ref.  ref is resized (if
               necessary) to hold num\_val values.  The value in the list
               indexed by start\_ix is copied into Array<double>[0],
               Array<double>[1] holds start\_ix+1, etc.  If the kth occurrence
               does not exist or there are fewer than start\_ix + num\_val
               values associated with the kth occurrence, or if some of the
               values cannot be converted to a double, an error message is
               reported and the program halts.
    */
    void getktharr (const char*    name,
                    int            k,
                    Array<double>& ref,
                    int            start_ix,
                    int            num_val) const;
    //
    //@ManDoc: Same as getktharr() but searches for last occurrence of name.
    //
    void getarr (const char*    name,
                 Array<double>& ref,
                 int            start_ix,
                 int            num_val) const;
    //
    //@ManDoc: queryktharr() is to querykth() as getktharr() is to getkth().
    //
    int queryktharr (const char*    name,
                     int            k,
                     Array<double>& ref,
                     int            start_ix,
                     int            num_val) const;
    //
    //@ManDoc: Same as queryktharr() but searches for last occurrence of name.
    //
    int queryarr (const char*    name,
                  Array<double>& ref,
                  int            start_ix,
                  int            num_val) const;

    /*@ManDoc: Gets an Array<std::string> of num\_val values from kth occurrence of
               given name.  If successful, the values are converted to an
               std::string and stored in the Array<std::string> object ref.  ref is
               resized (if necessary) to hold num\_val values.  The value in
               the list indexed by start\_ix is copied into Array<std::string>[0],
               Array<std::string>[1] holds start\_ix+1, etc.  If the kth occurrence
               does not exist or there are fewer than start\_ix + num\_val
               values associated with the kth occurrence, or if some of the
               values cannot be converted to an std::string, an error message is
               reported and the program halts.
    */
    void getktharr (const char*     name,
                    int             k,
                    Array<std::string>& ref,
                    int             start_ix,
                    int             num_val) const;
    //
    //@ManDoc: Same as getktharr() but searches for last occurrence of name.
    //
    void getarr (const char*     name,
                 Array<std::string>& ref,
                 int             start_ix,
                 int             num_val) const;
    //
    //@ManDoc: queryktharr() is to querykth() as getktharr() is to getkth().
    //
    int queryktharr (const char*     name,
                     int             k,
                     Array<std::string>& ref,
                     int             start_ix,
                     int             num_val) const;
    //
    //@ManDoc: Same as queryktharr() but searches for last occurrence of name.2
    //
    int queryarr (const char*     name,
                  Array<std::string>& ref,
                  int             start_ix,
                  int             num_val) const;

    /*@ManDoc: Gets an Array<IntVect> of num\_val values from kth occurrence of
               given name.  If successful, the values are converted to an
               IntVect and stored in the Array<IntVect> object ref.  ref is
               resized (if necessary) to hold num\_val values.  The value in
               the list indexed by start\_ix is copied into Array<IntVect>[0],
               Array<IntVect>[1] holds start\_ix+1, etc.  If the kth occurrence
               does not exist or there are fewer than start\_ix + num\_val
               values associated with the kth occurrence, or if some of the
               values cannot be converted to an IntVect, an error message is
               reported and the program halts.
    */
    void getktharr (const char*     name,
                    int             k,
                    Array<IntVect>& ref,
                    int             start_ix,
                    int             num_val) const;
    //
    //@ManDoc: Same as getktharr() but searches for last occurrence of name.
    //
    void getarr (const char*     name,
                 Array<IntVect>& ref,
                 int             start_ix,
                 int             num_val) const;
    //
    //@ManDoc: queryktharr() is to querykth() as getktharr() is to getkth().
    //
    int queryktharr (const char*     name,
                     int             k,
                     Array<IntVect>& ref,
                     int             start_ix,
                     int             num_val) const;
    //
    //@ManDoc: Same as queryktharr() but searches for last occurrence of name.2
    //
    int queryarr (const char*     name,
                  Array<IntVect>& ref,
                  int             start_ix,
                  int             num_val) const;

    /*@ManDoc: Gets an Array<Box> of num\_val values from kth occurrence of
               given name.  If successful, the values are converted to an
               Box and stored in the Array<Box> object ref.  ref is
               resized (if necessary) to hold num\_val values.  The value in
               the list indexed by start\_ix is copied into Array<Box>[0],
               Array<Box>[1] holds start\_ix+1, etc.  If the kth occurrence
               does not exist or there are fewer than start\_ix + num\_val
               values associated with the kth occurrence, or if some of the
               values cannot be converted to an Box, an error message is
               reported and the program halts.
    */
    void getktharr (const char*     name,
                    int             k,
                    Array<Box>& ref,
                    int             start_ix,
                    int             num_val) const;
    //
    //@ManDoc: Same as getktharr() but searches for last occurrence of name.
    //
    void getarr (const char*     name,
                 Array<Box>& ref,
                 int             start_ix,
                 int             num_val) const;
    //
    //@ManDoc: queryktharr() is to querykth() as getktharr() is to getkth().
    //
    int queryktharr (const char*     name,
                     int             k,
                     Array<Box>& ref,
                     int             start_ix,
                     int             num_val) const;
    //
    //@ManDoc: Same as queryktharr() but searches for last occurrence of name.2
    //
    int queryarr (const char*     name,
                  Array<Box>& ref,
                  int             start_ix,
                  int             num_val) const;

    /*@ManDoc: Construct an initial ParmParse object from the argc and argv
               passed in to main().  An error will be signalled if another
               ParmParse object currently exists.  If parfile is specified,
               read the parameters in from that file first and then append
               those derived from argv to the table.  If prefix is specified,
               load this string as the code prefix for this particular
               ParmParse object.
    */
    static void Initialize(int         argc,
			   char**      argv,
			   const char* parfile);


    /*@ManDoc: The destructor.  The internal static table will only be deleted
               if there are no other ParmParse objects in existence.
    */
    static void Finalize();
private:
    //
    // Prefix used in keyword search.
    //
    std::string thePrefix;
};

#endif /*BL_PARMPARSE_H*/
