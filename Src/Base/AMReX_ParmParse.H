
#ifndef BL_PARMPARSE_H
#define BL_PARMPARSE_H

#include <stack>
#include <string>
#include <iosfwd>
#include <vector>
#include <list>
#include <array>

namespace amrex {

class Box;
class IntVect;

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
// strings or numeric values by the query functions.
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
//     Vector<float> dx;
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
//    Vector<int> range;
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


/**
* \brief  Parse Parameters From Command Line and Input Files
*
*  The ParmParse class is used to interpret parameters passed in to a program
*  from the command line and an arbitrary collection of input files.  The
*  parameters are stored in static table that can be queried by any object
*  of type ParmParse.  A parameter is a  "definition".  A definition is
*  of the form "\<name\> = \<value\>\<value\>...\<value\>".  It is stored in the table
*  as a name, value-list pair.
*
*  In the following example, niter is a definition with the single
*  integer value 10; name is a definition with the string value "big
*  code" and dx is a definition with the two floating point values 0.5
*  and 0.75 and iv is an IntVect(5,4)
*
*  prog niter = 10 name = "big code" dx = 0.5 0.75 iv=(5,4)
*
*  The ParmParse class has two constructors.  The first is responsible for
*  building the table and is usually called by the main routine of an
*  application.  It has arguments for the command line argc and argv parameters,
*  as well as an optional filename argument for reading definitions from an
*  input file.  The table is built by reading the input file first (if it
*  exists) with the command line arguments added to the end of the table.
*  The order of a definition in the table is significant, so command line
*  parameters can be used to override definitions in the input file.  A
*  definition of the explicit form: FILE=\<filename\> is not added to the table
*  but is a directive to include the named file at that point in the table.
*
*  The second constructor is generally used by other classes in the code.  It
*  permits access to the table via a large collection of query functions.
*  Both constructors have an optional prefix argument that narrows the search
*  to entries in the table with the same prefix.  For example, let PlanR be a
*  ParmParse object with code prefix "ope".  PlanR.get("val",v) will look for
*  an entry in the parameter list of the form: ope.val==\<value\>, and will
*  reject all entries not starting with the correct code prefix.
*
*  The query functions search the table for definition names that match a given
*  string (and prefix) and return values from the corresponding value list.
*  The values can be returned as ints, Vector\<int\>s, floats, Vector\<float\>s,
*  doubles, Vector\<double\>s, std::strings, or Vector\<aSring\>s.  All values in the
*  table are stored as std::string objects, but if an int, float, or double is
*  requested, the translation is done automatically.  In the previous example,
*  the value of niter could be returned as either an std::string, an int, a double,
*  or a float.  The values of dx can be returned as std::strings, floats, or
*  doubles, but the value of name can be returned only as an std::string.
*
*  Comments in an input file include all text from a # character to the
*  end of the line.  Here is a sample input file:
*
*     niter = 100
*
*     title = "Double Wammy"
*
*     cell_size = 0.5 0.75
*
*     plot.var = Density 1 10
*
*     plot.var = Energy  5 12
*
*     bigarray = 1 2 3 4 5 6 7 8
*
*                9 10 11 12
*
*     aBox    = ((0,0) (5,5))
*
*     test = apple "boy blue" 10 20 30 40
*
*     FILE = prob_file
*/

class ParmParse
{
public:
    class Frame;
    class Record;

    enum { LAST = -1, FIRST = 0, ALL = -1 };
    /**
    * \brief Construct an additional ParmParse object sharing the same
    * internal table as any other such objects in existence.  If
    * prefix is specified, load this string as the code prefix
    * for this particular ParmParse object.
    */
    explicit ParmParse (const std::string& prefix = std::string());
    //! Returns true if name is in table.
    bool contains (const char* name) const;
    /**
    * \brief Returns the number of values associated with nth occurence of
    * name (prepended with the prefix) in the table.  n == -1 implies
    * the last occurence.
    */
    int countval (const char* name,
                  int         n = LAST) const;
    /**
    * \brief Returns the number of times the given name (prepended with
    * prefix) appears in the table.
    */
    int countname (const std::string& name) const;
    /**
    * \brief Returns the number of records of the given name (prepended with prefix) appears
    * in the table
    */
    int countRecords (const std::string& name) const;
    //! Returns the nth record of the given name, or zero if none found.
    Record getRecord (const std::string& name, int n = LAST) const;
    //! Write the contents of the table in ASCII to the ostream.
    static void dumpTable (std::ostream& os, bool prettyPrint = false);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a bool and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to an bool, an error message is
    * output and the program halts.   Note that ival == 0 is the first
    * value in the list.  ParmParse converte the value 'true', and non-zero
    * integers or floats to bool(true), and bool(false) for 'false'
    * or zero integer or float values.
    */
    void getkth (const char* name,
                 int         k,
                 bool&        ref,
                 int         ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (const char* name,
              bool&        ref,
              int         ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  bool&       ref,
                  int         ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (const char* name,
               bool&        ref,
               int         ival = FIRST) const;
    //! Add a key 'name'with value 'ref' to the end of the PP table.
    void add (const char* name,
              const bool  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to an int and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to an int, an error message is
    * output and the program halts.   Note that ival == 0 is the first
    * value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 int&        ref,
                 int         ival = FIRST) const;

    //! Same as getkth() but searches for the last occurrence of name.
    void get (const char* name,
              int&        ref,
              int         ival = FIRST) const;
    /**
    * Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  int&        ref,
                  int         ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (const char* name,
               int&        ref,
               int         ival = FIRST) const;
    //! Add a key 'name'with value 'ref' to the end of the PP table.
    void add (const char* name,
              const int  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to an int and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to an int, an error message is
    * output and the program halts.   Note that ival == 0 is the first
    * value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 long&       ref,
                 int         ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (const char* name,
              long&       ref,
              int         ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  long&       ref,
                  int         ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (const char* name,
               long&        ref,
               int         ival = FIRST) const;
    //! Add a key 'name'with value 'ref' to the end of the PP table.
    void add (const char* name,
              const long  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a float and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a float, an error message
    * is output and the program halts.   Note that ival == 0 is the
    * first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 float&      ref,
                 int         ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (const char* name,
              float&      ref,
              int         ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  float&      ref,
                  int         ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (const char* name,
               float&      ref,
               int         ival = FIRST) const;
    //! Add a key 'name'with value 'ref' to the end of the PP table.
    void add (const char* name,
              const float  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a double and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a double, an error message
    * is output and the program halts.   Note that ival = 0 is the
    * first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 double&     ref,
                 int         ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (const char* name,
              double&     ref,
              int         ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  double&     ref,
                  int         ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (const char* name,
               double&     ref,
               int         ival = FIRST) const;
    //! Add a key 'name'with value 'ref' to the end of the PP table.
    void add (const char* name,
              const double  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a std::string and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a std::string, an error message
    * is output and the program halts.   Note that ival = 0 is the
    * first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 std::string&    ref,
                 int         ival = FIRST) const;

    //! Same as getkth() but searches for the last occurrence of name.
    void get (const char* name,
              std::string&    ref,
              int         ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  std::string&   ref,
                  int         ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (const char* name,
               std::string&    ref,
               int         ival = FIRST) const;
    //! Add a key 'name'with value 'ref' to the end of the PP table.
    void add (const char* name,
              const std::string&  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to an IntVect and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a IntVect, an error message
    * is output and the program halts.   Note that ival = 0 is the
    * first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 IntVect&    ref,
                 int         ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (const char* name,
              IntVect&    ref,
              int         ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  IntVect&   ref,
                  int         ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (const char* name,
               IntVect&    ref,
               int         ival = FIRST) const;
    //! Add a key 'name'with value 'ref' to the end of the PP table.
    void add (const char* name,
              const IntVect&  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a Box and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a Box, an error message
    * is output and the program halts.   Note that ival = 0 is the
    * first value in the list.
    */
    void getkth (const char* name,
                 int         k,
                 Box&        ref,
                 int         ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (const char* name,
              Box&        ref,
              int         ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (const char* name,
                  int         k,
                  Box&        ref,
                  int         ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (const char* name,
               Box&        ref,
               int         ival = FIRST) const;
    //! Add a key 'name'with value 'ref' to the end of the PP table.
    void add (const char* name,
              const Box&  val);
    /**
    * \brief Gets an std::vector\<int\> of num_val values from kth occurrence of
    * given name.  If successful, the values are converted to an int
    * and stored in the std::vector\<int\> object ref.  ref is resized (if
    * necessary) to hold num_val values.  The value in the list
    * indexed by start_ix is copied into std::vector\<int\>[0], std::vector\<int\>[1]
    * holds start_ix+1, etc.  If the kth occurrence does not exist
    * or there are fewer than start_ix + num_val values associated
    * with the kth occurrence, or if some of the values cannot be
    * converted to an int, an error message is reported and the
    * program halts.
    */
    void getktharr (const char* name,
                    int         k,
                    std::vector<int>& ref,
                    int         start_ix = FIRST,
                    int         num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (const char* name,
                 std::vector<int>& ref,
                 int         start_ix = FIRST,
                 int         num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (const char* name,
                     int         k,
                     std::vector<int>& ref,
                     int         start_ix = FIRST,
                     int         num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (const char* name,
                  std::vector<int>& ref,
                  int         start_ix = FIRST,
                  int         num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ptr' to the end of the PP table.
    void addarr (const char* name,
                const std::vector<int>&  ptr);

    /**
    * \brief Gets an std::vector\<long\> of num_val values from kth occurrence of
    * given name.  If successful, the values are converted to a long
    * and stored in the std::vector\<long\> object ref.  ref is resized (if
    * necessary) to hold num_val values.  The value in the list
    * indexed by start_ix is copied into std::vector\<long\>[0], std::vector\<long\>[1]
    * holds start_ix+1, etc.  If the kth occurrence does not exist
    * or there are fewer than start_ix + num_val values associated
    * with the kth occurrence, or if some of the values cannot be
    * converted to a long, an error message is reported and the
    * program halts.
    */
    void getktharr (const char* name,
                    int         k,
                    std::vector<long>& ref,
                    int         start_ix = FIRST,
                    int         num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (const char* name,
                 std::vector<long>& ref,
                 int         start_ix = FIRST,
                 int         num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (const char* name,
                     int         k,
                     std::vector<long>& ref,
                     int         start_ix = FIRST,
                     int         num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (const char* name,
                  std::vector<long>& ref,
                  int         start_ix = FIRST,
                  int         num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ptr' to the end of the PP table.
    void addarr (const char* name,
                const std::vector<long>&  ptr);

    /**
    * \brief Gets an std::vector\<float\> of num_val values from kth occurrence of
    * given name.  If successful, the values are converted to a float
    * and stored in the std::vector\<float\> object ref.  ref is resized (if
    * necessary) to hold num_val values.  The value in the list
    * indexed by start_ix is copied into std::vector\<float\>[0],
    * std::vector\<float\>[1] holds start_ix+1, etc.  If the kth occurrence
    * does not exist or there are fewer than start_ix + num_val
    * values associated with the kth occurrence, or if some of the
    * values cannot be converted to a float, an error message is
    * reported and the program halts.
    */
    void getktharr (const char*   name,
                    int           k,
                    std::vector<float>& ref,
                    int           start_ix = FIRST,
                    int           num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (const char*   name,
                 std::vector<float>& ref,
                 int           start_ix = FIRST,
                 int           num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (const char*   name,
                     int           k,
                     std::vector<float>& ref,
                     int           start_ix = FIRST,
                     int           num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (const char*   name,
                  std::vector<float>& ref,
                  int           start_ix = FIRST,
                  int           num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ptr' to the end of the PP table.
    void addarr (const char* name,
                const std::vector<float>&  ptr);
    /**
    * \brief Gets an std::vector\<double\> of num_val values from kth occurrence of
    * given name.  If successful, the values are converted to a double
    * and stored in the std::vector\<double\> object ref.  ref is resized (if
    * necessary) to hold num_val values.  The value in the list
    * indexed by start_ix is copied into std::vector\<double\>[0],
    * std::vector\<double\>[1] holds start_ix+1, etc.  If the kth occurrence
    * does not exist or there are fewer than start_ix + num_val
    * values associated with the kth occurrence, or if some of the
    * values cannot be converted to a double, an error message is
    * reported and the program halts.
    */
    void getktharr (const char*    name,
                    int            k,
                    std::vector<double>& ref,
                    int            start_ix = FIRST,
                    int            num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (const char*    name,
                 std::vector<double>& ref,
                 int            start_ix = FIRST,
                 int            num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (const char*    name,
                     int            k,
                     std::vector<double>& ref,
                     int            start_ix = FIRST,
                     int            num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (const char*    name,
                  std::vector<double>& ref,
                  int            start_ix = FIRST,
                  int            num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ptr' to the end of the PP table.
    void addarr (const char* name,
                const std::vector<double>&  ptr);
    /**
    * \brief Gets an std::vector<std::string> of num_val values from kth occurrence of
    * given name.  If successful, the values are converted to an
    * std::string and stored in the std::vector<std::string> object ref.  ref is
    * resized (if necessary) to hold num_val values.  The value in
    * the list indexed by start_ix is copied into std::vector<std::string>[0],
    * std::vector<std::string>[1] holds start_ix+1, etc.  If the kth occurrence
    * does not exist or there are fewer than start_ix + num_val
    * values associated with the kth occurrence, or if some of the
    * values cannot be converted to an std::string, an error message is
    * reported and the program halts.
    */
    void getktharr (const char*     name,
                    int             k,
                    std::vector<std::string>& ref,
                    int             start_ix = FIRST,
                    int             num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (const char*     name,
                 std::vector<std::string>& ref,
                 int             start_ix = FIRST,
                 int             num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (const char*     name,
                     int             k,
                     std::vector<std::string>& ref,
                     int             start_ix = FIRST,
                     int             num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.2
    int queryarr (const char*     name,
                  std::vector<std::string>& ref,
                  int             start_ix = FIRST,
                  int             num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ptr' to the end of the PP table.
    void addarr (const char* name,
                const std::vector<std::string>&  ptr);
    /**
    * \brief Gets an std::vector\<IntVect\> of num_val values from kth occurrence of
    * given name.  If successful, the values are converted to an
    * IntVect and stored in the std::vector\<IntVect\> object ref.  ref is
    * resized (if necessary) to hold num_val values.  The value in
    * the list indexed by start_ix is copied into std::vector\<IntVect\>[0],
    * std::vector\<IntVect\>[1] holds start_ix+1, etc.  If the kth occurrence
    * does not exist or there are fewer than start_ix + num_val
    * values associated with the kth occurrence, or if some of the
    * values cannot be converted to an IntVect, an error message is
    * reported and the program halts.
    */
    void getktharr (const char*     name,
                    int             k,
                    std::vector<IntVect>& ref,
                    int             start_ix = FIRST,
                    int             num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (const char*     name,
                 std::vector<IntVect>& ref,
                 int             start_ix = FIRST,
                 int             num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (const char*     name,
                     int             k,
                     std::vector<IntVect>& ref,
                     int             start_ix = FIRST,
                     int             num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.2
    int queryarr (const char*     name,
                  std::vector<IntVect>& ref,
                  int             start_ix = FIRST,
                  int             num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ptr' to the end of the PP table.
    void addarr (const char* name,
                const std::vector<IntVect>&  ptr);
    /**
    * \brief Gets an std::vector\<Box\> of num_val values from kth occurrence of
    * given name.  If successful, the values are converted to an
    * Box and stored in the std::vector\<Box\> object ref.  ref is
    * resized (if necessary) to hold num_val values.  The value in
    * the list indexed by start_ix is copied into std::vector\<Box\>[0],
    * std::vector\<Box\>[1] holds start_ix+1, etc.  If the kth occurrence
    * does not exist or there are fewer than start_ix + num_val
    * values associated with the kth occurrence, or if some of the
    * values cannot be converted to an Box, an error message is
    * reported and the program halts.
    */
    void getktharr (const char*     name,
                    int             k,
                    std::vector<Box>& ref,
                    int             start_ix = FIRST,
                    int             num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (const char*     name,
                 std::vector<Box>& ref,
                 int             start_ix = FIRST,
                 int             num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (const char*     name,
                     int             k,
                     std::vector<Box>& ref,
                     int             start_ix = FIRST,
                     int             num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.2
    int queryarr (const char*     name,
                  std::vector<Box>& ref,
                  int             start_ix = FIRST,
                  int             num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ptr' to the end of the PP table.
    void addarr (const char* name,
                const std::vector<Box>&  ptr);

    template <typename T, std::size_t N>
    void get (const char* name, std::array<T,N>& ref) const {
        std::vector<T> v;
        this->getarr(name, v);
        AMREX_ALWAYS_ASSERT(v.size() >= N);
        for (std::size_t i = 0; i < N; ++i) {
            ref[i] = v[i];
        }
    }

    template <typename T, std::size_t N>
    int query (const char* name, std::array<T,N>& ref) const {
        std::vector<T> v;
        int exist = this->queryarr(name, v);
        if (exist) {
            AMREX_ALWAYS_ASSERT(v.size() >= N);
            for (std::size_t i = 0; i < N; ++i) {
                ref[i] = v[i];
            }
        }
        return exist;
    }

    /**
    * \brief Construct an initial ParmParse object from the argc and argv
    * passed in to main().  An error will be signalled if another
    * ParmParse object currently exists.  If parfile is specified,
    * read the parameters in from that file first and then append
    * those derived from argv to the table.  If prefix is specified,
    * load this string as the code prefix for this particular
    * ParmParse object.
    */
    static void Initialize(int         argc,
			   char**      argv,
			   const char* parfile);
    /**
    * \brief The destructor.  The internal static table will only be deleted
    * if there are no other ParmParse objects in existence.
    */
    static void Finalize();

    struct PP_entry;
    typedef std::list<PP_entry> Table;
    static void appendTable(ParmParse::Table& tab);
    const Table& table() const {return m_table;}

private:
    //
    // Don't allow one ParmParse object to be assigned to another.
    // Given that we use a const reference for "m_table" there's
    // really no way to make assignment work.  It really doesn't make
    // much sense anyway.  There's really only one "table" used by ParmParse,
    // a copy of which is maintained by each MPI process.  I'm adding this
    // solely so that the error message the compiler generates is a little
    // clearer than the one generated by the default compiler-generated
    // assignment operator.
    //
    ParmParse& operator= (const ParmParse& rhs);

    friend class Frame;
    friend class Record;

    explicit ParmParse (const Table& tbl);
    //
    //! Set/Get the prefix.
    std::string getPrefix() const;
    std::string setPrefix(const std::string& str);
    void pushPrefix(const std::string& str);
    void popPrefix();
    std::string prefixedName (const std::string& str) const;
    //
    //! Prefix used in keyword search.
    std::stack<std::string> m_pstack;
    const Table& m_table;
};

struct ParmParse::PP_entry
{
    PP_entry (const std::string&            name,
              const std::list<std::string>& vals);
    PP_entry (const std::string& name,
              const std::string& vals);
    PP_entry (const std::string& name,
	      const Table& table);
    PP_entry (const PP_entry& pe);
    PP_entry& operator= (const PP_entry& pe);
    ~PP_entry ();
    std::string print() const;

    std::string              m_name;
    std::vector<std::string> m_vals;
    Table*                   m_table;
    mutable bool             m_queried;
};


class ParmParse::Frame
{
public:
    Frame (ParmParse& pp, const std::string& str);
    ~Frame ();
    void push(const std::string& str);
    void pop();
    std::string getPrefix() const;
private:
    ParmParse& m_pp;
    int        m_np;
};

class ParmParse::Record
{
public:
    const ParmParse* operator->() const;
    const ParmParse& operator* () const;
private:
    friend class ParmParse;
    explicit Record (const ParmParse& pp);
    ParmParse m_pp;
};

std::ostream& operator<< (std::ostream& os, const ParmParse::PP_entry& pp);

}

#endif /*BL_PARMPARSE_H*/
