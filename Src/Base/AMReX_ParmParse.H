#ifndef AMREX_PARMPARSE_H_
#define AMREX_PARMPARSE_H_
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Enum.H>
#include <AMReX_INT.H>
#include <AMReX_IntVect.H>
#include <AMReX_IParser.H>
#include <AMReX_Parser.H>
#include <AMReX_TypeTraits.H>

#include <array>
#include <iosfwd>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

namespace amrex {

template<int dim>
class BoxND;
using Box = BoxND<AMREX_SPACEDIM>;
template<int dim>
class IntVectND;
using IntVect = IntVectND<AMREX_SPACEDIM>;
template<int dim>
class RealVectND;
using RealVect = RealVectND<AMREX_SPACEDIM>;

namespace ppdetail {
    template <class T, class Enable = void>
    struct ArithmeticOptional_TT : std::false_type {};

    template <class T>
    struct ArithmeticOptional_TT<T, std::enable_if_t<std::is_arithmetic_v<T>>>
        : std::true_type
    {
        using value_type = T;
    };

    template <class T>
    struct ArithmeticOptional_TT<std::optional<T>,
                                 std::enable_if_t<std::is_arithmetic_v<T>>>
        : std::true_type
    {
        using value_type = T;
    };

    template <class T>
    inline constexpr bool IsArithmeticOptional_v = ArithmeticOptional_TT<T>::value;

    template <class T>
    using underlying_type_t = typename ArithmeticOptional_TT<T>::value_type;
}

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
/*
   niter = 100                # niter is an integer
   title = "Double Wammy"     # example of a string with spaces
   cell_size = 0.5 0.75       # cell spacing in each dimension
   plot.var = Density 1 10    # a list of values
   plot.var = Energy  5 12    # another list of values
   bigarray = 1 2 3 4 5 6 7 8 \
              9 10 11 12      # continuation of bigarray
   multi_line_string = "This is a
                        multi-line string."
   my_2d_table = \
           # col 1      2      3
             {{ 11.0,  12.0,  13.0 }   # row 1
              { 21.0,  22.0,  23.0 }   # row 2
              { 31.0,  32.0,  33.0 }   # row 3
              { 41.0,  42.0,  43.0 } } # row 4
   test = apple "boy blue" 10 20 30 40
   FILE = prob_file           # insert contents of this "prob_file" here
*/
// For values spanning multiple lines except for table, one must use '\' at
// the end of a line
// for continuation, otherwise it's a runtime error. Note that there must be
// at least one space before the continuation character `\`.  Multiple lines
// inside a pair of double quotes are considered a single string containing
// '\n's. The "FILE = <filename>" definition is special.  Rather than just
// adding this entry to the database, it reads the contents of <filename>
// into the database.
// For CI/CD workflows and out-of-source tests, the environment variable
// AMREX_INPUTS_FILE_PREFIX can be set to prefix every FILE = <filename>
// with a custom path.
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
// following regular expression:
//   Sign    ::= '+' | '-'
//   Digit   ::= '0' | '1' | ... | '9'
//   Integer ::= [Sign]Digit+
//   Exp     ::= ('e'|'E')Integer
//   Float   ::= ( Integer[.Digit*][Exp] | [Integer].Digit+[Exp] )
//
// Where '+' indicates one or more occurrences, '*' represents zero or
// more occurrences, '|' means one or the other and '[]' represents zero
// or one occurrence.
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
*  ```
*  executable niter = 10 name = "big code" dx = 0.5 0.75 iv=(5,4)
*  ```
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
* ```
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
*     bigarray = 1 2 3 4 5 6 7 8 \
*                9 10 11 12
*
*     multi_line_string = "This is a
*                          multi-line string."
*
*     aBox    = ((0,0) (5,5))
*
*     test = apple "boy blue" 10 20 30 40
*
*     FILE = prob_file
* ```
*
*  Preprocessing of AMREX_SPACEDIM is supported. It supports `#if`, `#elif`,
*  `#else`, and `#endif`.  The condition must be `AMREX_SPACEDIM op D`, where op
*  is >, <, >=, <=, or ==, and D is 1, 2, or 3. The parentheses around the
*  condition are optional. Some examples are shown below.
*
* \code
*    #if (AMREX_SPACEDIM == 1)
*    n_cell = 256
*    #elif (AMREX_SPACEDIM == 2)
*    n_cell = 128 128
*    #else
*    n_cell = 64 64 64
*    #endif
*
*    #if AMREX_SPACEDIM >= 2
*    t = 0.5
*    #else
*    t = 1.5
*    #endif
* \endcode
*
*  Math expression is supported for integers and reals. For example
*
* ```
*    n_cell = 128
*    amrex.n_cell = n_cell*2 8 16**2
* ```
*
*  becomes
*
* ```
*    n_cell = 128
*    amrex.n_cell = 256 8 256
* ```
*
*  More details can be found at https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse
*/
class ParmParse
{
public:
    enum { LAST = -1, FIRST = 0, ALL = -1 };
    /**
    * \brief Construct an additional ParmParse object sharing the same
    * internal table as any other such objects in existence.  If prefix is
    * specified, load this string as the code prefix for this particular
    * ParmParse object. If parser_prefix is specified, it will be used as
    * prefixed in math expression evaluations.
    */
    explicit ParmParse (std::string prefix = std::string(),
                        std::string parser_prefix = std::string());

    //! Returns true if name is in table.
    [[nodiscard]] bool contains (std::string_view name) const;
    /**
    * \brief Returns the number of values associated with nth occurrence of
    * name (prepended with the prefix) in the table.  n == -1 implies
    * the last occurrence.
    */
    [[nodiscard]] int countval (std::string_view name, int n = LAST) const;
    /**
    * \brief Returns the number of times the given name (prepended with
    * prefix) appears in the table.
    */
    [[nodiscard]] int countname (std::string_view name) const;
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a bool and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to an bool, an error message is
    * output and the program halts.   Note that ival == 0 is the first
    * value in the list.  ParmParse converts the value 'true', and non-zero
    * integers or floats to bool(true), and bool(false) for 'false'
    * or zero integer or float values.
    */
    void getkth (std::string_view name,
                 int              k,
                 bool&            ref,
                 int              ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              bool&            ref,
              int              ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  bool&            ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               bool&            ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, bool val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to an int and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to an int, an error message is
    * output and the program halts.  Note that ival == 0 is the first
    * value in the list.
    */
    void getkth (std::string_view name,
                 int              k,
                 int&             ref,
                 int              ival = FIRST) const;

    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              int&             ref,
              int              ival = FIRST) const;
    /**
    * Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  int&             ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               int&             ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, int val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to an int and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to an int, an error message is
    * output and the program halts.  Note that ival == 0 is the first
    * value in the list.
    */
    void getkth (std::string_view name,
                 int              k,
                 long&            ref,
                 int              ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              long&            ref,
              int              ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  long&            ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               long&            ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, long val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to an int and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to an int, an error message is
    * output and the program halts.   Note that ival == 0 is the first
    * value in the list.
    */
    void getkth (std::string_view name,
                 int              k,
                 long long&       ref,
                 int              ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              long long&       ref,
              int              ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  long long&       ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               long long&       ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, long long  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a float and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a float, an error message
    * is output and the program halts.   Note that ival == 0 is the
    * first value in the list.
    */
    void getkth (std::string_view name,
                 int              k,
                 float&           ref,
                 int              ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              float&           ref,
              int              ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  float&           ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               float&           ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, float  val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a double and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a double, an error message
    * is output and the program halts.   Note that ival = 0 is the
    * first value in the list.
    */
    void getkth (std::string_view name,
                 int         k,
                 double&     ref,
                 int         ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              double&          ref,
              int              ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  double&          ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               double&          ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, double val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a std::string and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a std::string, an error message
    * is output and the program halts.   Note that ival = 0 is the
    * first value in the list.
    */
    void getkth (std::string_view name,
                 int              k,
                 std::string&     ref,
                 int              ival = FIRST) const;

    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              std::string&     ref,
              int              ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  std::string&     ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               std::string&     ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, const std::string& val);


    //! Similar to get() but store the whole line including empty spaces in
    //! the middle if present as a single string. For example, `foo = a b c`
    //! is treated by this function as a single string "a b c", whereas
    //! get() will give only "a" and getarr() will store the results in
    //! `std::vectcor<std::string>`. Note this does not preserve the number
    //! of empty spaces, because of how ParmParse parses the line. For
    //! example, `a b  c` with two spaces between `b` and `c` will become "a
    //! b c".
    void getline (std::string_view name, std::string& ref) const;

    //! Similar to query() but store the whole line including empty spaces
    //! in the middle if present as a single string. For example, `foo = a b
    //! c` is treated by this function as a single string "a b c", whereas
    //! query() will give only "a" and queryarr() will store the results in
    //! `std::vectcor<std::string>`. Note this does not preserve the number
    //! of empty spaces, because of how ParmParse parses the line. For
    //! example, `a b  c` with two spaces between `b` and `c` will become "a
    //! b c".
    int queryline (std::string_view name, std::string& ref) const;

    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to an IntVect and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a IntVect, an error message
    * is output and the program halts.   Note that ival = 0 is the
    * first value in the list.
    */
    void getkth (std::string_view name,
                 int              k,
                 IntVect&         ref,
                 int              ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              IntVect&         ref,
              int              ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  IntVect&         ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               IntVect&         ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, const IntVect& val);
    /**
    * \brief Get the ival'th value of kth occurrence of the requested name.
    * If successful, the value is converted to a Box and stored
    * in reference ref.  If the kth occurrence does not exist or
    * ival'th value does not exist, or if the printed representation
    * of the value cannot be converted to a Box, an error message
    * is output and the program halts.   Note that ival = 0 is the
    * first value in the list.
    */
    void getkth (std::string_view name,
                 int              k,
                 Box&             ref,
                 int              ival = FIRST) const;
    //! Same as getkth() but searches for the last occurrence of name.
    void get (std::string_view name,
              Box&             ref,
              int              ival = FIRST) const;
    /**
    * \brief Similar to getkth() but returns 0 if there is no kth occurrence
    * of name.  If successful, it returns 1 and stores the value in
    * ref.  If the kth occurrence exists, but ival'th value of that
    * occurrence does not, or if there is a type mismatch, then the
    * program signals an error and halts.
    */
    int querykth (std::string_view name,
                  int              k,
                  Box&             ref,
                  int              ival = FIRST) const;
    //! Same as querykth() but searches for the last occurrence of name.
    int query (std::string_view name,
               Box&             ref,
               int              ival = FIRST) const;
    //! Add a key 'name' with value 'val' to the end of the PP table.
    void add (std::string_view name, const Box& val);
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
    void getktharr (std::string_view  name,
                    int               k,
                    std::vector<int>& ref,
                    int               start_ix = FIRST,
                    int               num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (std::string_view  name,
                 std::vector<int>& ref,
                 int               start_ix = FIRST,
                 int               num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (std::string_view  name,
                     int               k,
                     std::vector<int>& ref,
                     int               start_ix = FIRST,
                     int               num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (std::string_view  name,
                  std::vector<int>& ref,
                  int               start_ix = FIRST,
                  int               num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ref' to the end of the PP table.
    void addarr (std::string_view name, const std::vector<int>& ref);

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
    void getktharr (std::string_view   name,
                    int                k,
                    std::vector<long>& ref,
                    int                start_ix = FIRST,
                    int                num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (std::string_view   name,
                 std::vector<long>& ref,
                 int                start_ix = FIRST,
                 int                num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (std::string_view   name,
                     int                k,
                     std::vector<long>& ref,
                     int                start_ix = FIRST,
                     int                num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (std::string_view   name,
                  std::vector<long>& ref,
                  int                start_ix = FIRST,
                  int                num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ref' to the end of the PP table.
    void addarr (std::string_view name, const std::vector<long>& ref);

    /**
    * \brief Gets an std::vector\<long long\> of num_val values from kth occurrence of
    * given name.  If successful, the values are converted to a long long
    * and stored in the std::vector\<long long\> object ref.  ref is resized (if
    * necessary) to hold num_val values.  The value in the list
    * indexed by start_ix is copied into std::vector\<long long\>[0], std::vector\<long long\>[1]
    * holds start_ix+1, etc.  If the kth occurrence does not exist
    * or there are fewer than start_ix + num_val values associated
    * with the kth occurrence, or if some of the values cannot be
    * converted to a long long, an error message is reported and the
    * program halts.
    */
    void getktharr (std::string_view        name,
                    int                     k,
                    std::vector<long long>& ref,
                    int                     start_ix = FIRST,
                    int                     num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (std::string_view        name,
                 std::vector<long long>& ref,
                 int                     start_ix = FIRST,
                 int                     num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (std::string_view        name,
                     int                     k,
                     std::vector<long long>& ref,
                     int                     start_ix = FIRST,
                     int                     num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (std::string_view        name,
                  std::vector<long long>& ref,
                  int                     start_ix = FIRST,
                  int                     num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ref' to the end of the PP table.
    void addarr (std::string_view name, const std::vector<long long>& ref);

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
    void getktharr (std::string_view    name,
                    int                 k,
                    std::vector<float>& ref,
                    int                 start_ix = FIRST,
                    int                 num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (std::string_view    name,
                 std::vector<float>& ref,
                 int                 start_ix = FIRST,
                 int                 num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (std::string_view    name,
                     int                 k,
                     std::vector<float>& ref,
                     int                 start_ix = FIRST,
                     int                 num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (std::string_view    name,
                  std::vector<float>& ref,
                  int                 start_ix = FIRST,
                  int                 num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ref' to the end of the PP table.
    void addarr (std::string_view name, const std::vector<float>& ref);
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
    void getktharr (std::string_view     name,
                    int                  k,
                    std::vector<double>& ref,
                    int                  start_ix = FIRST,
                    int                  num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (std::string_view     name,
                 std::vector<double>& ref,
                 int                  start_ix = FIRST,
                 int                  num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (std::string_view     name,
                     int                  k,
                     std::vector<double>& ref,
                     int                  start_ix = FIRST,
                     int                  num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.
    int queryarr (std::string_view     name,
                  std::vector<double>& ref,
                  int                  start_ix = FIRST,
                  int                  num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ref' to the end of the PP table.
    void addarr (std::string_view name, const std::vector<double>& ref);
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
    void getktharr (std::string_view          name,
                    int                       k,
                    std::vector<std::string>& ref,
                    int                       start_ix = FIRST,
                    int                       num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (std::string_view          name,
                 std::vector<std::string>& ref,
                 int                       start_ix = FIRST,
                 int                       num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (std::string_view          name,
                     int                       k,
                     std::vector<std::string>& ref,
                     int                       start_ix = FIRST,
                     int                       num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.2
    int queryarr (std::string_view          name,
                  std::vector<std::string>& ref,
                  int                       start_ix = FIRST,
                  int                       num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ref' to the end of the PP table.
    void addarr (std::string_view name, const std::vector<std::string>& ref);
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
    void getktharr (std::string_view      name,
                    int                   k,
                    std::vector<IntVect>& ref,
                    int                   start_ix = FIRST,
                    int                   num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (std::string_view      name,
                 std::vector<IntVect>& ref,
                 int                   start_ix = FIRST,
                 int                   num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (std::string_view      name,
                     int                   k,
                     std::vector<IntVect>& ref,
                     int                   start_ix = FIRST,
                     int                   num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.2
    int queryarr (std::string_view      name,
                  std::vector<IntVect>& ref,
                  int                   start_ix = FIRST,
                  int                   num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ref' to the end of the PP table.
    void addarr (std::string_view name, const std::vector<IntVect>& ref);
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
    void getktharr (std::string_view  name,
                    int               k,
                    std::vector<Box>& ref,
                    int               start_ix = FIRST,
                    int               num_val = ALL) const;
    //! Same as getktharr() but searches for last occurrence of name.
    void getarr (std::string_view  name,
                 std::vector<Box>& ref,
                 int               start_ix = FIRST,
                 int               num_val = ALL) const;
    //! queryktharr() is to querykth() as getktharr() is to getkth().
    int queryktharr (std::string_view  name,
                     int               k,
                     std::vector<Box>& ref,
                     int               start_ix = FIRST,
                     int               num_val = ALL) const;
    //! Same as queryktharr() but searches for last occurrence of name.2
    int queryarr (std::string_view  name,
                  std::vector<Box>& ref,
                  int               start_ix = FIRST,
                  int               num_val = ALL) const;
    //! Add a key 'name' with vector of values 'ref' to the end of the PP table.
    void addarr (std::string_view name, const std::vector<Box>& ref);

    /*
     * \brief Query IntVect from array
     *
     * This reads IntVect from an array (e.g., `8 16 8`), not the format
     * using parentheses (e.g., `(8,16,8)`).
     */
    int queryarr (std::string_view name, IntVect& ref) const;

    /*
     * \brief Get IntVect from array
     *
     * This reads IntVect from an array (e.g., `8 16 8`), not the format
     * using parentheses (e.g., `(8,16,8)`).
     */
    void getarr (std::string_view name, IntVect& ref) const;

    //! Query RealVect from array
    int queryarr (std::string_view name, RealVect& ref) const;

    //! Get RealVect from array
    void getarr (std::string_view name, RealVect& ref) const;

    template <typename T, std::size_t N>
    void get (std::string_view name, std::array<T,N>& ref) const {
        std::vector<T> v;
        this->getarr(name, v);
        AMREX_ALWAYS_ASSERT(v.size() >= N);
        for (std::size_t i = 0; i < N; ++i) {
            ref[i] = v[i];
        }
    }

    template <typename T, std::size_t N>
    int query (std::string_view name, std::array<T,N>& ref) const {
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
     * \brief If `name` is found, the value in the ParmParse database will
     * be stored in the `ref` argument.  If not, the value in `ref` will be
     * added to the ParmParse database.  The return value indicates if it
     * existed previously.
     */
    template <typename T, std::enable_if_t<!IsStdVector<T>::value, int> = 0>
    int queryAdd (std::string_view name, T& ref) {
        int exist = this->query(name, ref);
        if (!exist) {
            this->add(name, ref);
        }
        return exist;
    }

    int queryAdd (std::string_view name, std::string& ref) {
        int exist = this->query(name, ref);
        if (!exist && !ref.empty()) {
            this->add(name, ref);
        }
        return exist;
    }

    /**
     * \brief If `name` is found, the value in the ParmParse database will
     * be stored in the `ref` argument.  If not, the value in `ref` will be
     * added to the ParmParse database.  The return value indicates if it
     * existed previously.
     *
     * If `name` is found, then the ref argument will be reallocated
     * (and resized) according to the number of values in the inputs.
     */
    template <typename T>
    int queryAdd (std::string_view name, std::vector<T>& ref) {
        std::vector<T> empty;
        int exist = this->queryarr(name, empty);
        if (exist) {
            ref = std::move(empty);
        }
        if (!exist && !ref.empty()) {
            this->addarr(name, ref);
        }
        return exist;
    }

    /**
     * \brief If `name` is found, the value in the ParmParse database will
     * be stored in the `ref` argument.  If not, the value in `ref` will be
     * added to the ParmParse database.  The return value indicates if it
     * existed previously.
     */
    template <typename T>
    int queryAdd (std::string_view name, std::vector<T>& ref, int num_val) {
        int exist = this->queryarr(name, ref, 0, num_val);
        if (!exist) {
            this->addarr(name, ref);
        }
        return exist;
    }

    /**
     * \brief If `name` is found, the value in the ParmParse database will
     * be stored in the `ref` argument.  If not, the value in `ref` will be
     * added to the ParmParse database.  The return value indicates if it
     * existed previously.
     */
    template <typename T, std::size_t N>
    int queryAdd (std::string_view name, std::array<T,N>& ref) {
        std::vector<T> v;
        int exist = this->queryarr(name, v);
        if (exist) {
            AMREX_ALWAYS_ASSERT(v.size() >= N);
            for (std::size_t i = 0; i < N; ++i) {
                ref[i] = v[i];
            }
        } else {
            v.resize(N);
            for (std::size_t i = 0; i < N; ++i) {
                v[i] = ref[i];
            }
            this->addarr(name, v);
        }
        return exist;
    }

    /**
     * \brief Query with Parser. If `name` is found, this uses amrex::Parser
     * to parse the entire list of empty space separated values as a single
     * scalar. The return value indicates whether it's found. Note that
     * queryWithParser will be used recursively for unresolved symbols.
     */
    int queryWithParser (std::string_view name, bool& ref) const;
    int queryWithParser (std::string_view name, int& ref) const;
    int queryWithParser (std::string_view name, long& ref) const;
    int queryWithParser (std::string_view name, long long& ref) const;
    int queryWithParser (std::string_view name, float& ref) const;
    int queryWithParser (std::string_view name, double& ref) const;

    /**
     * \brief Query with Parser. The return value indicates whether it's
     * found. Note that queryWithParser will be used for unresolved
     * symbols. If the number of elements in the input does not equal to
     * `nvals`, it's a runtime error.
     */
    int queryarrWithParser (std::string_view name, int nvals, bool* ptr) const;
    int queryarrWithParser (std::string_view name, int nvals, int* ptr) const;
    int queryarrWithParser (std::string_view name, int nvals, long* ptr) const;
    int queryarrWithParser (std::string_view name, int nvals, long long* ptr) const;
    int queryarrWithParser (std::string_view name, int nvals, float* ptr) const;
    int queryarrWithParser (std::string_view name, int nvals, double* ptr) const;
    template <typename T, std::enable_if_t<std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>,int> = 0>
    int queryarrWithParser (std::string_view name, int nvals, std::vector<T>& ref) const
    {
        if (this->contains(name)) {
            if (int(ref.size()) < nvals) { ref.resize(nvals); }
            return this->queryarrWithParser(name, nvals, ref.data());
        } else {
            return 0;
        }
    }

    /**
     * \brief Query with Parser. If `name` is found, this uses amrex::Parser
     * to parse the entire list of empty space separated values as a single
     * scalar. If not, the value in `ref` will be added to the ParmParse
     * database. The return value indicates whether it's found.
     */
    template <typename T, std::enable_if_t<std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>,int> = 0>
    int queryAddWithParser (std::string_view name, T& ref)
    {
        int exist = this->queryWithParser(name, ref);
        if (!exist) {
            this->add(name, ref);
        }
        return exist;
    }

    /**
     * \brief Get with Parser. If `name` is found, this uses amrex::Parser
     * to parse the entire list of empty space separated values as a single
     * scalar. If not, it's a runtime error.
     */
    template <typename T, std::enable_if_t<std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>,int> = 0>
    void getWithParser (std::string_view name, T& ref) const
    {
        int exist = this->queryWithParser(name, ref);
        if (!exist) {
            amrex::Error(std::string("ParmParse::getWithParser: failed to get ")+std::string(name));
        }
    }

    /**
     * \brief Get with Parser. If `name` is not found, it's a runtime
     * error. If the number of elements does not equal to `nvals`, it's also
     * a runtime error.
     */
    template <typename T, std::enable_if_t<std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>,int> = 0>
    void getarrWithParser (std::string_view name, int nvals, T* ptr) const
    {
        int exist = this->queryarrWithParser(name, nvals, ptr);
        if (!exist) {
            amrex::Error(std::string("ParmParse::getarrWithParser: failed to get ")+std::string(name));
        }
    }

    /**
     * \brief Get with Parser. If `name` is not found, it's a runtime
     * error. If the number of elements does not equal to `nvals`, it's also
     * a runtime error.
     */
    template <typename T, std::enable_if_t<std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>,int> = 0>
    void getarrWithParser (std::string_view name, int nvals, std::vector<T>& ref) const
    {
        int exist = this->queryarrWithParser(name, nvals, ref);
        if (!exist) {
            amrex::Error(std::string("ParmParse::getarrWithParser: failed to get ")+std::string(name));
        }
    }

    /*
     * \brief Evaluate given string as math expression
     *
     * For unknown symbols, ParmParse database will be queried.
     */
    template <typename T, std::enable_if_t<std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>,int> = 0>
    T eval (std::string const& expr) const
    {
        if constexpr (std::is_integral_v<T>) {
            auto const parser = this->makeIParser(expr, {});
            auto const exe = parser.compileHost<0>();
            return static_cast<T>(exe()); // In the future, we might add safety check.
        } else {
            auto const parser = this->makeParser(expr, {});
            auto const exe = parser.compileHost<0>();
            return static_cast<T>(exe());
        }
    }

    /*
     * \brief Query two names.
     *
     * This function queries with `new_name` first. If it's not found, it
     * will try again with `old_name`.
     */
    template <typename T>
    int query (const char* new_name, const char* old_name, T& ref)
    {
        return (this->query(new_name, ref) ||
                this->query(old_name, ref));
    }

    /**
     * \brief Get using two names.
     *
     * This function queries with `new_name` first, If it's not found, it
     * will try again with `old_name`. It's an error if neither name is
     * found.
     */
    template <typename T>
    void get (const char* new_name, const char* old_name, T& ref)
    {
        auto exist = this->query(new_name, old_name, ref);
        if (!exist) {
            amrex::ErrorStream() << "ParmParse::get failed to find "
                                 << new_name << " and " << old_name << '\n';
            ParmParse::dumpTable(amrex::ErrorStream());
            amrex::Abort();
        }
    }

    /**
     * \brief. Query enum value using given name.
     *
     * Here T is an enum class defined by AMREX_ENUM. The return value
     * indicates if `name` is found. An exception is thrown, if the found
     * string associated with the name cannot be converted to an enumerator
     * (i.e., the string does not match any names in the definition of T).
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    int query (std::string_view name, T& ref, int ival = FIRST) const
    {
        std::string s;
        int exist = this->query(name, s, ival);
        if (exist) {
            try {
                ref = amrex::getEnum<T>(s);
            } catch (...) {
                if (amrex::Verbose() > 0 ) {
                    amrex::Print() << "amrex::ParmParse::query (input name: "
                                   << this->prefixedName(name) << "):\n";
                }
                throw;
            }
        }
        return exist;
    }


    //! Add a key 'name' with value 'val' to the end of the ParmParse
    //! database. Here T is an AMREX_ENUM type.
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    void add (std::string_view name, T const& val)
    {
        this->add(name, amrex::getEnumNameString(val));
    }

    /**
     * \brief. Get enum value using given name.
     *
     * Here T is an enum class defined by AMREX_ENUM. It's a runtime error,
     * if `name` is not found. An exception is thrown, if the found string
     * associated with the name cannot be converted to an enumerator (i.e.,
     * the string does not match any names in the definition of T).
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    void get (std::string_view name, T& ref, int ival = FIRST) const
    {
        std::string s;
        this->get(name, s, ival);
        try {
            ref = amrex::getEnum<T>(s);
        } catch (...) {
            if (amrex::Verbose() > 0 ) {
                amrex::Print() << "amrex::ParmParse::get (input name:  "
                               << this->prefixedName(name) << "):\n";
            }
            throw;
        }
    }

    //! Query an array of enum values using given name.
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    int queryarr (std::string_view name,
                  std::vector<T>& ref,
                  int start_ix = FIRST,
                  int num_val = ALL) const
    {
        std::vector<std::string> s;
        int exist = this->queryarr(name, s, start_ix, num_val);
        if (exist) {
            ref.resize(s.size());
            for (std::size_t i = 0; i < s.size(); ++i) {
                try {
                    ref[i] = amrex::getEnum<T>(s[i]);
                }  catch (...) {
                    if (amrex::Verbose() > 0 ) {
                        amrex::Print() << "amrex::ParmParse::queryarr (input name:  "
                                       << this->prefixedName(name) << "):\n";
                    }
                    throw;
                }
            }
        }
        return exist;
    }

    //! Get an array of enum values using given name.
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    void getarr (std::string_view name,
                 std::vector<T>& ref,
                 int start_ix = FIRST,
                 int num_val = ALL) const
    {
        std::vector<std::string> s;
        this->getarr(name, s, start_ix, num_val);
        ref.resize(s.size());
        for (std::size_t i = 0; i < s.size(); ++i) {
            try {
                ref[i] = amrex::getEnum<T>(s[i]);
            }  catch (...) {
                if (amrex::Verbose() > 0 ) {
                    amrex::Print() << "amrex::ParmParse::getarr (input name:  "
                                   << this->prefixedName(name) << "):\n";
                }
                throw;
            }
        }
    }

    /**
     * \brief. Query enum value using given name.
     *
     * Here T is an enum class defined by AMREX_ENUM. The return value
     * indicates if `name` is found. An exception is thrown, if the found
     * string associated with the name cannot be case-insensitively
     * converted to an enumerator (i.e., the found string, not `name`, does
     * not case-insensitively match any names in the definition of T). If
     * there are multiple matches, the first one is used.
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    int query_enum_case_insensitive (std::string_view name, T& ref, int ival = FIRST) const
    {
        std::string s;
        int exist = this->query(name, s, ival);
        if (exist) {
            try {
                ref = amrex::getEnumCaseInsensitive<T>(s);
            } catch (...) {
                if (amrex::Verbose() > 0) {
                    amrex::Print() << "amrex::ParmParse::query_enum_case_insensitive (input name: "
                                   << this->prefixedName(name) << "):\n";
                }
                throw;
            }
        }
        return exist;
    }

    /**
     * \brief. Get enum value using given name.
     *
     * Here T is an enum class defined by AMREX_ENUM. It's a runtime error,
     * if `name` is not found. An exception is thrown, if the found string
     * associated with the name cannot be case-insensitively converted to an
     * enumerator (i.e., the found string, not `name`, does not
     * case-insensitively match any names in the definition of T). If there
     * are multiple matches, the first one is used.
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    void get_enum_case_insensitive (std::string_view name, T& ref, int ival = FIRST) const
    {
        int exist = this->query_enum_case_insensitive(name, ref, ival);
        if (!exist) {
            std::string msg("get_enum_case_insensitive(\"");
            msg.append(name).append("\",").append(amrex::getEnumClassName<T>())
                .append("&) failed.");
            amrex::Abort(msg);
        }
    }

    /**
     * \brief. Query enum value using given name.
     *
     * Here T is an enum class defined by AMREX_ENUM. The return value
     * indicates if `name` is found. An exception is thrown, if the found
     * string associated with the name cannot be case-insensitively
     * converted to an enumerator (i.e., the found string, not `name`, does
     * not case-insensitively match any names in the definition of T). If
     * there are multiple matches, the first one is used. Characters in
     * `ignores` will be ignored as if they don't exist in the value part of
     * ParamParse entries (e.g., `name = value`).
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    int query_enum_sloppy (std::string_view name, T& ref, std::string_view const& ignores,
                           int ival = FIRST) const
    {
        std::string s;
        int exist = this->query(name, s, ival);
        if (exist) {
            try {
                s.erase(std::remove_if(s.begin(), s.end(),
                                       [&] (auto const& c) {
                                           return ignores.find(c) != std::string_view::npos; }),
                        s.end());
                ref = amrex::getEnumCaseInsensitive<T>(s);
            } catch (...) {
                if (amrex::Verbose() > 0) {
                    amrex::Print() << "amrex::ParmParse::query_enum_sloppy (input name: "
                                   << this->prefixedName(name) << "):\n";
                }
                throw;
            }
        }
        return exist;
    }

    /**
     * \brief. Get enum value using given name.
     *
     * Here T is an enum class defined by AMREX_ENUM. It's a runtime error,
     * if `name` is not found. An exception is thrown, if the found string
     * associated with the name cannot be case-insensitively converted to an
     * enumerator (i.e., the found string, not `name`, does not
     * case-insensitively match any names in the definition of T). If there
     * are multiple matches, the first one is used. Characters in `ignores`
     * will be ignored as if they don't exist in the value part of
     * ParamParse entries (e.g., `name = value`).
     */
    template <typename T, typename ET = amrex_enum_traits<T>,
              std::enable_if_t<ET::value,int> = 0>
    void get_enum_sloppy (std::string_view name, T& ref, std::string_view const& ignores,
                          int ival = FIRST) const
    {
        int exist = this->query_enum_sloppy(name, ref, ignores, ival);
        if (!exist) {
            std::string msg("get_enum_sloppy(\"");
            msg.append(name).append("\",").append(amrex::getEnumClassName<T>())
                .append("&) failed.");
            amrex::Abort(msg);
        }
    }

    /**
     * \brief Query T with Parser, but treat the number as double precision
     * during parsing.
     *
     * The final result is cast to T. It may result in a runtime error if
     * the conversion is not safe. T is either arithmetic type or
     * std::optional of arithmetic type.
     */
    template <typename T, std::enable_if_t<ppdetail::IsArithmeticOptional_v<T>, int> = 0>
    int queryAsDouble (std::string_view name, T& ref) const
    {
        using value_type = ppdetail::underlying_type_t<T>;
        double dref;
        int exist = queryWithParser(name, dref);
        if (exist) {
            if (std::is_integral_v<value_type>) {
                dref = std::round(dref);
            }
            auto vref = static_cast<value_type>(dref);
            if constexpr (std::is_integral_v<value_type> && !std::is_same_v<value_type,bool>) {
                if (static_cast<double>(vref) != dref) {
                    amrex::Abort("ParmParse:: queryAsDouble is not safe");
                }
            }
            ref = vref;
        }
        return exist;
    }

    /**
     * \brief Query T array with Parser, but treat the number as double
     * precision during parsing.
     *
     * The final result is cast to T's. It may result in a runtime error if
     * the conversion is not safe. T is either arithmetic type or
     * std::optional of arithmetic type.
     */
    template <typename T, std::enable_if_t<ppdetail::IsArithmeticOptional_v<T>, int> = 0>
    int queryarrAsDouble (std::string_view name, int nvals, T* ptr) const
    {
        using value_type = ppdetail::underlying_type_t<T>;
        std::vector<double> dref(nvals);
        int exist = queryarrWithParser(name, nvals, dref.data());
        if (exist) {
            for (int i = 0; i < nvals; ++i) {
                if (std::is_integral_v<value_type>) {
                    dref[i] = std::round(dref[i]);
                }
                auto vref = static_cast<value_type>(dref[i]);
                if constexpr (std::is_integral_v<value_type> && !std::is_same_v<value_type,bool>) {
                    if (static_cast<double>(vref) != dref[i]) {
                        amrex::Abort("ParmParse:: queryarrAsDouble is not safe");
                    }
                }
                ptr[i] = vref;
            }
        }
        return exist;
    }

    /**
     * \brief Get T with Parser, but treat the number as double precision
     * during parsing.
     *
     * The final result is cast to T. It may result in a runtime error if
     * the conversion is not safe. T is either arithmetic type or
     * std::optional of arithmetic type.
     */
    template <typename T, std::enable_if_t<ppdetail::IsArithmeticOptional_v<T>, int> = 0>
    void getAsDouble (std::string_view name, T& ref) const
    {
        int exist = this->queryAsDouble(name, ref);
        if (!exist) {
            amrex::Error(std::string("ParmParse::getAsDouble: failed to get ")+std::string(name));
        }
    }

    /**
     * \brief Get T array with Parser, but treat the number as double
     * precision during parsing.
     *
     * The final result is cast to T's. It may result in a runtime error if
     * the conversion is not safe. T is either arithmetic type or
     * std::optional of arithmetic type.
     */
    template <typename T, std::enable_if_t<ppdetail::IsArithmeticOptional_v<T>, int> = 0>
    void getarrAsDouble (std::string_view name, int nvals, T* ptr) const
    {
        int exist = this->queryarrAsDouble(name, nvals, ptr);
        if (!exist) {
            amrex::Error(std::string("ParmParse::getarrAsDouble: failed to get ")+std::string(name));
        }
    }

    //! Remove given name from the table.
    int remove (std::string_view name);

    //! Make Parser using given string `func` as function body and `vars` as
    //! variable names. Constants known to ParmParse will be set. It's a
    //! runtime error, if there are unknown symbols in `func`.
    [[nodiscard]] Parser makeParser (std::string const& func,
                                     Vector<std::string> const& vars) const;

    //! Make IParser using given string `func` as function body and `vars`
    //! as variable names. Constants known to ParmParse will be set. It's a
    //! runtime error, if there are unknown symbols in `func`.
    [[nodiscard]] IParser makeIParser (std::string const& func,
                                       Vector<std::string> const& vars) const;

    /**
     * \brief Query vector of vector. The return value indicates whether
     * it's found. The table (i.e., vector of vector) is in the format of
     * `{{a00, a01, a02...} {a10, a11, a12...} ...}`.
     */
    int querytable (std::string_view name, std::vector<std::vector<double>>& ref) const;
    int querytable (std::string_view name, std::vector<std::vector<float>>& ref) const;
    int querytable (std::string_view name, std::vector<std::vector<int>>& ref) const;

    /**
     * \brief Get vector of vector. It's an error if it is not found.  The
     * table (i.e., vector of vector) is in the format of
     * `{{a00, a01, * a02...} {a10, a11, a12...} ...}`.
     */
    template <typename T>
    void gettable (std::string_view name, std::vector<std::vector<T>>& ref) const
    {
        if (this->querytable(name, ref) == 0) {
            amrex::ErrorStream() << "ParmParse::gettable: " << name
                                 << " not found in database\n";
            ParmParse::dumpTable(amrex::ErrorStream());
            amrex::Abort();
        }
    }

    /**
    * \brief Construct an initial ParmParse object from the argc and argv
    * passed in to main().  An error will be signalled if another
    * ParmParse object currently exists.  If parfile is specified,
    * read the parameters in from that file first and then append
    * those derived from argv to the table.
    */
    static void Initialize (int argc, char** argv, const char* parfile);
    static void Initialize (int argc, char** argv, const std::string& parfile) {
        return Initialize(argc, argv, parfile.c_str());
    }
    /**
    * \brief The destructor.  The internal static table will only be deleted
    * if there are no other ParmParse objects in existence.
    */
    static void Finalize ();

    //! Set prefix used by math expression Parser
    static void SetParserPrefix (std::string a_prefix);

    static int Verbose ();
    static void SetVerbose (int v);

    //! Write the contents of the table in ASCII to the ostream.
    static void dumpTable (std::ostream& os, bool prettyPrint = false);

    //! Write all entries in the table in a pretty way to the ostream. If there are
    //! duplicates, only the last one is printed.
    static void prettyPrintTable (std::ostream& os);

    //! Write unused inputs in a pretty way to the ostream. If there are
    //! duplicates, only the last one is printed.
    static void prettyPrintUnusedInputs (std::ostream& os);

    //! Write used inputs in a pretty way to the ostream. If there are
    //! duplicates, only the last one is printed.
    static void prettyPrintUsedInputs (std::ostream& os);

    //! Add keys and values from a file to the end of the PP table.
    static void addfile (std::string const& filename);

    static bool QueryUnusedInputs ();

    //! Any unused [prefix.]* parameters?
    [[nodiscard]] static bool hasUnusedInputs (const std::string& prefix = std::string());

    //! Returns unused [prefix.]* parameters.
    [[nodiscard]] static std::vector<std::string> getUnusedInputs (const std::string& prefix = std::string());

    //! Returns [prefix.]* parameters.
    [[nodiscard]] static std::set<std::string> getEntries (const std::string& prefix = std::string());

    struct PP_entry {
        // There can be multiple occurrences for a given name (e.g.,
        // multiple lines starting with `foo =` in inputs. For each
        // occurrence, there can be multiple values. Thus, the use of
        // vector<vector<std::string>>.
        std::vector<std::vector<std::string>> m_vals;
        mutable Long                          m_count = 0;
        std::variant<
            std::string*,
            bool*,
            int*,
            long*,
            long long*,
            amrex::IntVect*,
            amrex::Box*,
            float*,
            double*
        > m_typehint = static_cast<std::string*>(nullptr);
    };
    using Table = std::unordered_map<std::string, PP_entry>;

    [[nodiscard]] const Table& table() const {return *m_table;}

    //! keyword for files to load
    static std::string const FileKeyword;

    static std::string ParserPrefix;

    [[nodiscard]] std::string const& getPrefix () const;

    [[nodiscard]] std::string prefixedName (std::string_view str) const;

protected:

    std::string m_prefix; // Prefix used in keyword search
    std::string m_parser_prefix; // Prefix used by Parser
    Table* m_table;
};

}

#endif /* AMREX_PARMPARSE_H_ */
