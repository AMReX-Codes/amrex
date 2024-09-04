#include <AMReX_ParmParse.H>
#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_RealVect.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <typeinfo>
#include <type_traits>

extern "C" void amrex_init_namelist (const char*);
extern "C" void amrex_finalize_namelist ();

namespace amrex {

namespace {
    bool initialized = false;
    ParmParse::Table g_table;
    std::vector<std::set<std::string>> g_parser_recursive_symbols;
    namespace pp_detail {
        int verbose = -1;
    }
}

std::string const ParmParse::FileKeyword = "FILE";
std::string       ParmParse::ParserPrefix;

ParmParse::ParmParse (std::string prefix, std::string parser_prefix)
    : m_prefix(std::move(prefix)),
      m_parser_prefix(std::move(parser_prefix)),
      m_table(&g_table)
{}

namespace
{

std::string pp_to_pretty_string (std::string const& name,
                                 std::vector<std::string> const& vals)
{
    std::stringstream ss;
    ss << name << " =";
    for (auto const& v : vals) {
        ss << " " << v;
    }
    return ss.str();
}

std::string pp_to_string (std::string const& name,
                          std::vector<std::string> const& vals)
{
    std::stringstream ss;
    ss << name << "(nvals = " << vals.size() << ") " << " :: [";
    for (std::size_t i = 0; i < vals.size(); ++i) {
        ss << vals[i];
        if ( i < vals.size()-1 ) { ss << ", "; }
    }
    ss << "]";
    return ss.str();
}

enum PType
{
    pDefn,
    pEQ_sign,
    pValue,
    pEOF
};

template <class T>
bool
isT (const std::string& str, T& val)
{
    std::istringstream s(str);
    s >> val;
    if ( s.fail() ) { return false; }
    std::string left;
    std::getline(s, left);
    if ( !left.empty() ) { return false; }
    return true;
}

template <typename T, std::enable_if_t<std::is_floating_point_v<T>,int> = 0>
bool
is_floating_point (const std::string& str, T& val)
{
    if (str == "nan") {
        val = std::numeric_limits<T>::quiet_NaN();
        return true;
    } else if (str == "inf") {
        val = std::numeric_limits<T>::infinity();
        return true;
    } else if (str == "-inf") {
        val = -std::numeric_limits<T>::infinity();
        return true;
    } else {
        return isT(str, val);
    }
}

template <class T>
bool
is (const std::string& str, T& val)
{
    return isT(str, val);
}

template <>
bool
is (const std::string& str, float& val)
{
    return is_floating_point(str, val);
}

template <>
bool
is (const std::string& str, double& val)
{
    return is_floating_point(str, val);
}

template <>
bool
is (const std::string& str, std::string& val)
{
    val = str;
    return true;
}

template <>
bool
is (const std::string& str, bool& val)
{
    auto const lo_str = amrex::toLower(str);
    if ( lo_str == "true" || lo_str == "t" )
    {
        val = true;
        return true;
    }
    if ( lo_str == "false" || lo_str == "f" )
    {
        val = false;
        return true;
    }
    int int_val;
    if ( isT(str, int_val) )
    {
        val = int_val != 0;
        return true;
    }
    double dbl_val;
    if ( isT(str, dbl_val) )
    {
        val = dbl_val != 0;
        return true;
    }
    return false;
}

template <class T> const char* tok_name(const T&) { return typeid(T).name(); }
template <class T> const char* tok_name(std::vector<T>&) { return tok_name(T());}

//
// Simple lexical analyser.
//

enum lexState
{
    START,
    STRING,
    QUOTED_STRING,
    IDENTIFIER,
    LIST
};

const char* const
state_name[] =
{
   "START",
   "STRING",
   "QUOTED_STRING",
   "IDENTIFIER",
   "LIST"
};

int
eat_garbage (const char*& str)
{
    int num_linefeeds = 0;
    for (;;)
    {
        if ( *str == 0 ) { break; } // NOLINT
        else if ( *str == '#' )
        {
            while ( *str && *str != '\n' )
            {
                str++;
            }
            continue;
        }
        else if ( std::isspace(*str) )
        {
            if (*str == '\n') { ++num_linefeeds; }
            str++;
        }
        else if ( *str == '\\' ) // '\' followed by a line break is continuation to next line
        {
            // Unfortunately, a line break has three variants, \r, \n, and \r\n.
            if (*(str+1) == '\n') {
                str += 2;
            } else if (*(str+1) == '\r') {
                if (*(str+2) == '\n') {
                    str += 3;
                } else {
                    str += 2;
                }
            } else {
                break;
            }
        }
        else
        {
            break;
        }
    }
    return num_linefeeds;
}

PType
getToken (const char*& str, std::string& ostr, int& num_linefeeds)
{
   //
   // Eat white space and comments.
   //
   num_linefeeds = eat_garbage(str);
   //
   // Check for end of file.
   //
   if ( *str == 0 )
   {
       return pEOF;
   }
   //
   // Start token scan.
   //
   lexState state = START;
   int      pcnt  = 0; // Tracks nested parens
   while (true)
   {
       char ch = *str;
       if ( ch == 0 )
       {
           amrex::Error("ParmParse::getToken: EOF while parsing");
       }
       switch (state)
       {
       case START:
           if ( ch == '=' )
           {
               ostr += ch; str++;
               return pEQ_sign;
           }
           else if ( ch == '"' )
           {
               str++;
               state = QUOTED_STRING;
           }
           else if ( ch == '(' )
           {
               ostr += ch; str++; pcnt = 1;
               state = LIST;
           }
           else if ( std::isalpha(ch) )
           {
               ostr += ch; str++;
               state = IDENTIFIER;
           }
           else
           {
               ostr += ch; str++;
               state = STRING;
           }
           break;
       case IDENTIFIER:
           if ( std::isalnum(ch) || ch == '_' || ch == '.' || ch == '[' || ch == ']' || ch == '+' || ch == '-' )
           {
               ostr += ch; str++;
           }
           else if ( std::isspace(ch) || ch == '=' )
           {
               return pDefn;
           }
           else
           {
               ostr += ch; str++;
               state = STRING;
           }
           break;
       case LIST:
           if ( ch == '(' )
           {
               ostr += ch; str++; pcnt++;
           }
           else if ( ch == ')' )
           {
               ostr += ch; str++; pcnt--;
               if ( pcnt == 0 )
               {
                   return pValue;
               }
           }
           else
           {
               ostr += ch; str++;
           }
           break;
       case STRING:
           if ( std::isspace(ch) || ch == '=' )
           {
               return pValue;
           }
           else
           {
               ostr += ch; str++;
           }
           break;
       case QUOTED_STRING:
           if ( ch == '"' )
           {
               str++;
               return pValue;
           }
           else
           {
               ostr += ch; str++;
           }
           break;
       default:
           amrex::ErrorStream() << "ParmParse::getToken(): invalid string = " << ostr << '\n'
                                << "STATE = " << state_name[state]
                                << ", next char = " << ch << '\n'
                                << ", rest of input = \n" << str << '\n';
           amrex::Abort();
       }
   }
}

//
// Return the index of the n'th occurrence of a parameter name,
// except if n==-1, return the index of the last occurrence.
// Return 0 if the specified occurrence does not exist.
//
std::vector<std::string> const*
ppindex (const ParmParse::Table& table, int n, const std::string& name)
{
    auto found = table.find(name);
    if (found == table.cend()) { return nullptr; }

#ifdef AMREX_USE_OMP
#pragma omp atomic update
#endif
    ++(found->second.m_count);

    if (n == ParmParse::LAST) {
        return &(found->second.m_vals.back());
    } else {
        return &(found->second.m_vals[n]);
    }
}

void bldTable (const char*& str, ParmParse::Table& tab);

bool isTrue(std::smatch const& sm)
{
    const std::string op = sm[1].str();
    const int dim = std::stoi(sm[2].str());
    if (op == "<") {
        return AMREX_SPACEDIM < dim;
    } else if (op == ">") {
        return AMREX_SPACEDIM > dim;
    } else if (op == "==") {
        return AMREX_SPACEDIM == dim;
    } else if (op == "<=") {
        return AMREX_SPACEDIM <= dim;
    } else if (op == ">=") {
        return AMREX_SPACEDIM >= dim;
    } else {
        return false;
    }
}

void
read_file (const char* fname, ParmParse::Table& tab)
{
    //
    // Space for input file if it exists.
    //
    if ( fname != nullptr && fname[0] != 0 )
    {
#ifdef AMREX_USE_MPI
        if (ParallelDescriptor::Communicator() == MPI_COMM_NULL)
        {
            throw std::runtime_error("read_file: AMReX must be initialized");
        }
#endif

        Vector<char> fileCharPtr;
        std::string filename = fname;
        ParallelDescriptor::ReadAndBcastFile(filename, fileCharPtr);

        std::istringstream is(fileCharPtr.data());
        std::ostringstream os_cxx(std::ios_base::out);
        std::ostringstream os_fortran(std::ios_base::out);
        bool fortran_namelist = false;
        std::regex if_regex("^\\s*#\\s*if\\s+\\(?\\s*AMREX_SPACEDIM\\s*(>|<|==|>=|<=)\\s*([1-3])\\s*\\)?\\s*$"); // NOLINT
        std::regex elif_regex("^\\s*#\\s*elif\\s+\\(?\\s*AMREX_SPACEDIM\\s*(>|<|==|>=|<=)\\s*([1-3])\\s*\\)?\\s*$"); // NOLINT
        std::regex else_regex("^\\s*#\\s*else\\s*$"); // NOLINT
        std::regex endif_regex("^\\s*#\\s*endif\\s*$"); // NOLINT
        std::vector<bool> valid_region;  // Keep this block or not?
        std::vector<bool> has_true;      // Has previous if/elif ever been true?
        for (std::string line; std::getline(is, line); ) {
            std::smatch sm;
            if (std::regex_match(line, sm, if_regex)) {
                bool r = isTrue(sm);
                valid_region.push_back(r);
                has_true.push_back(r);
                continue;
            } else if (std::regex_match(line, sm, elif_regex)) {
                if (has_true.back() == false) {
                    // If none of the previous if/elif is true
                    bool r = isTrue(sm);
                    valid_region.back() = r;
                    has_true.back() = r;
                } else {
                    // If any of the previous if/elif is true
                    valid_region.back() = false;
                }
                continue;
            } else if (std::regex_match(line, sm, else_regex)) {
                if (has_true.back() == false) {
                    // If none of the previous if/elif is true,
                    valid_region.back() = true;
                } else {
                    valid_region.back() = false;
                }
                continue;
            } else if (std::regex_match(line, sm, endif_regex)) {
                valid_region.pop_back();
                has_true.pop_back();
                continue;
            }

            if (std::find(std::begin(valid_region), std::end(valid_region), false)
                != std::end(valid_region)) {
                continue;
            }

            auto r = std::find_if(std::begin(line), std::end(line),
                                  [](int c) -> bool { return !std::isspace(c); });
            if (fortran_namelist) { // already inside fortran namelist
                // os_fortran << line << "\n";
                // pgi and ibm do not like `\n`.  We strip comments for them too.
                os_fortran << line.substr(0, line.find('!')) << " ";
                if (r != std::end(line) && *r == '/') {
                    fortran_namelist = false; // end of Fortran namelist
                }
            } else if (r != std::end(line) && *r == '&') {
                // os_fortran << line << "\n";
                // pgi and ibm do not like `\n`.  We strip comments for them too.
                os_fortran << line.substr(0, line.find('!')) << " ";
                fortran_namelist = true;  // begin of Fortran namelist
            } else {
                os_cxx << line << "\n";
            }
        }

        std::string filestring_cxx = os_cxx.str();
        const char* b = filestring_cxx.c_str();
        bldTable(b, tab);

#if !defined(BL_NO_FORT)
        std::string filestring_fortran = os_fortran.str();
        amrex_init_namelist(filestring_fortran.c_str());
#endif
    }
}

void
addDefn (std::string& def, std::vector<std::string>& val, ParmParse::Table& tab)
{
    //
    // Check that defn exists.
    //
    if ( def.empty() )
    {
        val.clear();
        return;
    }
    //
    // Check that it has values.
    //
    if ( val.empty() )
    {
        amrex::ErrorStream() << "ParmParse::addDefn(): no values for definition " << def << "\n";
        amrex::Abort();
    }
    //
    // Check if this defn is a file include directive.
    //
    if ( def == ParmParse::FileKeyword && val.size() == 1 )
    {
        //
        // Read file and add to this table.
        //
        const char* fname = val.front().c_str();
        read_file(fname, tab);
    }
    else
    {
        tab[def].m_vals.push_back(val);
    }
    val.clear();
    if ( def != ParmParse::FileKeyword ) {
        def = std::string();
    }
}

void
bldTable (const char*& str, ParmParse::Table& tab)
{
    std::string              cur_name;
    std::vector<std::string> cur_list;
    std::vector<int>         cur_linefeeds;

    for (;;)
    {
        std::string tokname;
        int num_linefeeds;

        PType token = getToken(str, tokname, num_linefeeds);

        switch (token)
        {
        case pEOF:
        {
            if (std::accumulate(cur_linefeeds.begin(), cur_linefeeds.end(), int(0)) > 0)
            {
                std::string error_message("ParmParse: Multiple lines in ");
                error_message.append(cur_name).append(" =");
                for (auto const& x : cur_list) {
                    error_message.append(" ").append(x);
                }
                error_message.append(". Must use \\ for line continuation.");
                amrex::Abort(error_message);
            }
            addDefn(cur_name,cur_list,tab);
            return;
        }
        case pEQ_sign:
        {
            if ( cur_name.empty() )
            {
                amrex::Abort("ParmParse::bldTable() EQ with no current defn");
            }
            if ( !cur_list.empty() )
            {
                //
                // Read one too far, remove last name on list.
                //
                auto tmp_str = cur_list.back();
                cur_list.pop_back();
                cur_linefeeds.pop_back();
                if (std::accumulate(cur_linefeeds.begin(), cur_linefeeds.end(), int(0)) > 0)
                {
                    std::string error_message("ParmParse: Multiple lines in ");
                    error_message.append(cur_name).append(" =");
                    for (auto const& x : cur_list) {
                        error_message.append(" ").append(x);
                    }
                    error_message.append(". Must use \\ for line continuation.");
                    amrex::Abort(error_message);
                }
                addDefn(cur_name,cur_list,tab);
                cur_name = std::move(tmp_str);
            }
            cur_linefeeds.clear();
            break;
        }
        case pDefn:
        {
            if ( cur_name.empty() )
            {
                cur_name = tokname;
                break;
            }
            //
            // Otherwise, fall through, this may be a string.
            //
            AMREX_FALLTHROUGH;
        }
        case pValue:
        {
            if ( cur_name.empty() )
            {
                std::string msg("ParmParse::bldTable(): value with no defn: ");
                msg += tokname;
                amrex::Abort(msg.c_str());
            }
            cur_list.push_back(tokname);
            cur_linefeeds.push_back(num_linefeeds);
            break;
        }
        } // switch (token)
    }
}

template <typename T>
bool pp_parser (const ParmParse::Table& table, const std::string& parser_prefix,
                const std::string& name, const std::string& val, T& ref,
                bool use_querywithparser);

template <class T>
bool
squeryval (const ParmParse::Table& table,
           const std::string&      parser_prefix,
           const std::string&      name,
           T&                      ref,
           int                     ival,
           int                     occurrence)
{
    //
    // Get last occurrence of name in table.
    //
    auto const* def = ppindex(table, occurrence, name);
    if ( def == nullptr )
    {
        return false;
    }
    //
    // Does it have ival values?
    //
    if ( ival >= static_cast<int>(def->size()) )
    {
        amrex::ErrorStream() << "ParmParse::queryval no value number"
                             << ival << " for ";
        if ( occurrence ==  ParmParse::LAST )
        {
            amrex::ErrorStream() << "last occurrence of ";
        }
        else
        {
            amrex::ErrorStream() << " occurrence " << occurrence << " of ";
        }
        amrex::ErrorStream() << name << '\n' << pp_to_string(name,*def) << '\n';
        amrex::Abort();
    }

    const std::string& valname = (*def)[ival];

    bool ok = is(valname, ref);
    if ( !ok )
    {
        if constexpr (std::is_same_v<T,int> ||
                      std::is_same_v<T,long> ||
                      std::is_same_v<T,long long> ||
                      std::is_same_v<T,float> ||
                      std::is_same_v<T,double>)
        {
            if (pp_parser(table, parser_prefix, name, valname, ref, false)) {
                return true;
            }
        } else {
            amrex::ignore_unused(parser_prefix);
        }

        amrex::ErrorStream() << "ParmParse::queryval type mismatch on value number "
                             << ival << " of " << '\n';
        if ( occurrence == ParmParse::LAST )
        {
            amrex::ErrorStream() << " last occurrence of ";
        }
        else
        {
            amrex::ErrorStream() << " occurrence number " << occurrence << " of ";
        }
        amrex::ErrorStream() << name << '\n';
        amrex::ErrorStream() << " Expected an \""
                             << tok_name(ref)
                             << "\" type  which can't be parsed from the string \""
                             << valname << "\"\n"
                             << pp_to_string(name,*def) << '\n';
        amrex::Abort();
    }
    return true;
}

template <class T>
void
sgetval (const ParmParse::Table& table,
         const std::string&      parser_prefix,
         const std::string&      name,
         T&                      ref,
         int                     ival,
         int                     occurrence)
{
    if ( squeryval(table, parser_prefix, name,ref,ival,occurrence) == 0 )
    {
        amrex::ErrorStream() << "ParmParse::getval ";
        if ( occurrence >= 0 )
        {
            amrex::ErrorStream() << "occurrence number "
                                 << occurrence
                                 << " of ";
        }

        amrex::ErrorStream() << "ParmParse::getval(): "
                             << name
                             << " not found in table"
                             << '\n';
        ParmParse::dumpTable(amrex::ErrorStream());
        amrex::Abort();
    }
}

template <class T>
bool
squeryarr (const ParmParse::Table& table,
           const std::string&      parser_prefix,
           const std::string&      name,
           std::vector<T>&         ref,
           int                     start_ix,
           int                     num_val,
           int                     occurrence)
{
    //
    // Get last occurrence of name in table.
    //
    auto const* def = ppindex(table,occurrence, name);
    if ( def == nullptr )
    {
        return false;
    }
    //
    // Does it have sufficient number of values and are they all
    // the same type?
    //
    if ( num_val == ParmParse::ALL )
    {
        num_val = static_cast<int>(def->size());
    }

    if ( num_val == 0 ) { return true; }

    int stop_ix = start_ix + num_val - 1;
    if ( static_cast<int>(ref.size()) <= stop_ix )
    {
        ref.resize(stop_ix + 1);
    }
    if ( stop_ix >= static_cast<int>(def->size()) )
    {
        amrex::ErrorStream() << "ParmParse::queryarr too many values requested for";
        if ( occurrence == ParmParse::LAST )
        {
            amrex::ErrorStream() << " last occurrence of ";
        }
        else
        {
            amrex::ErrorStream() << " occurrence " << occurrence << " of ";
        }
        amrex::ErrorStream() << name << '\n' << pp_to_string(name,*def) << '\n';
        amrex::Abort();
    }
    for ( int n = start_ix; n <= stop_ix; n++ )
    {
        const std::string& valname = (*def)[n];
        bool ok = is(valname, ref[n]);
        if ( !ok )
        {
            if constexpr (std::is_same_v<T,int> ||
                          std::is_same_v<T,long> ||
                          std::is_same_v<T,long long> ||
                          std::is_same_v<T,float> ||
                          std::is_same_v<T,double>)
            {
                if (pp_parser(table, parser_prefix, name, valname, ref[n], false)) {
                    continue;
                }
            } else {
                amrex::ignore_unused(parser_prefix);
            }

            amrex::ErrorStream() << "ParmParse::queryarr type mismatch on value number "
                                 <<  n << " of ";
            if ( occurrence == ParmParse::LAST )
            {
                amrex::ErrorStream() << " last occurrence of ";
            }
            else
            {
                amrex::ErrorStream() << " occurrence number " << occurrence << " of ";
            }
            amrex::ErrorStream() << name << '\n';
            amrex::ErrorStream() << " Expected an \""
                                 << tok_name(ref)
                                 << "\" type which can't be parsed from the string \""
                                 << valname << "\"\n"
                                 << pp_to_string(name,*def) << '\n';
            amrex::Abort();
        }
    }
    return true;
}

template <class T>
void
sgetarr (const ParmParse::Table& table,
         const std::string&      parser_prefix,
         const std::string&      name,
         std::vector<T>&         ref,
         int                     start_ix,
         int                     num_val,
         int                     occurrence)
{
    if ( squeryarr(table,parser_prefix,name,ref,start_ix,num_val,occurrence) == 0 )
    {
        amrex::ErrorStream() << "ParmParse::sgetarr ";
        if ( occurrence >= 0 )
        {
            amrex::ErrorStream() << "occurrence number " << occurrence << " of ";
        }
        amrex::ErrorStream() << "ParmParse::sgetarr(): "
                             << name
                             << " not found in table"
                             << '\n';
        ParmParse::dumpTable(amrex::ErrorStream());
        amrex::Abort();
    }
}

template <class T>
void
saddval (const std::string& name, const T& ref)
{
    std::stringstream val;
    val << std::setprecision(17) << ref;

    auto& entry = g_table[name];
    entry.m_vals.emplace_back(std::vector<std::string>{val.str()});
    ++entry.m_count;
}

template <class T>
void
saddarr (const std::string& name, const std::vector<T>& ref)
{
    std::vector<std::string> arr;
    arr.reserve(ref.size());
    for (auto const& item : ref) {
        std::stringstream val;
        val << std::setprecision(17) << item;
        arr.push_back(val.str());
    }

    auto& entry = g_table[name];
    entry.m_vals.emplace_back(std::move(arr));
    ++entry.m_count;
}

// Initialize ParmParse.
void
ppinit (int argc, char** argv, const char* parfile, ParmParse::Table& table)
{
    if ( parfile != nullptr )
    {
        read_file(parfile, table);
    }

    if ( argc > 0 )
    {
        std::string argstr;
        const char SPACE = ' ';
        for ( int i = 0; i < argc; i++ )
        {
            argstr += argv[i];
            argstr += SPACE;
        }
        ParmParse::Table arg_table;
        const char* b = argstr.c_str();
        bldTable(b, arg_table);
        //
        // Append arg_table to end of existing table.
        //
        for (auto& [name, arg_entry] : arg_table) {
            auto& src = arg_entry.m_vals;
            auto& dst = table[name].m_vals;
            std::move(std::begin(src), std::end(src), std::back_inserter(dst));
        }
    }
    initialized = true;
}

bool unused_table_entries_q (const ParmParse::Table& table,
                             const std::string& prefix = std::string())
{
    if (prefix.empty()) {
        return std::any_of(table.begin(), table.end(),
                           [] (auto const& x) -> bool {
                               return x.second.m_count == 0;
                           });
    } else {
        auto s = prefix + '.';
        return std::any_of(table.begin(), table.end(),
                           [&] (auto const& x) -> bool {
                               return x.second.m_count == 0
                                   && x.first.substr(0,s.size()) == s;
                           });
    }
}

void pp_print_unused (const std::string& pfx, const ParmParse::Table& table)
{
    std::vector<std::string> sorted_names;
    sorted_names.reserve(table.size());
    for (auto const& [name, entry] : table) {
        if (entry.m_count == 0) {
            sorted_names.push_back(name);
        }
    }
    std::sort(sorted_names.begin(), sorted_names.end());

    for (auto const& name : sorted_names) {
        auto const& entry = table.at(name);
        for (auto const& vals : entry.m_vals) {
            amrex::AllPrint() << pfx << "::" << pp_to_string(name, vals) << '\n';
        }
    }
}

template <class T>
bool squeryWithParser (const ParmParse::Table& table,
                       const std::string&      parser_prefix,
                       const std::string&      name,
                       T&                      ref);

template <typename T, typename PARSER_t = std::conditional_t<std::is_integral_v<T>,
                                                             IParser, Parser>>
PARSER_t
pp_make_parser (std::string const& func, Vector<std::string> const& vars,
                ParmParse::Table const& table, std::string const& parser_prefix,
                bool use_querywithparser)
{
    using value_t =  std::conditional_t<std::is_integral_v<T>, long long, double>;

    std::vector<std::string> prefixes;
    prefixes.reserve(3);
    prefixes.emplace_back();
    if (! parser_prefix.empty()) {
        prefixes.emplace_back(parser_prefix+".");
    }
    if (! ParmParse::ParserPrefix.empty()) {
        prefixes.emplace_back(ParmParse::ParserPrefix+".");
    }

    PARSER_t parser(func);

    auto symbols = parser.symbols();
    for (auto const& var : vars) {
        symbols.erase(var);
    }

    for (auto const& s : symbols) {
        value_t v = 0;
        bool r = false;
        for (auto const& pf : prefixes) {
            if (use_querywithparser) {
                r = squeryWithParser(table, parser_prefix, pf+s, v);
            } else {
                r = squeryval(table, parser_prefix, pf+s, v,
                              ParmParse::FIRST, ParmParse::LAST);
            }
            if (r) { break; }
        }
        if (r == false) {
            amrex::Error("ParmParse: failed to parse " + func);
        }
        parser.setConstant(s, v);
    }
    if (!vars.empty()) {
        parser.registerVariables(vars);
    }

    return parser;
}

template <typename T>
bool pp_parser (const ParmParse::Table& table, const std::string& parser_prefix,
                const std::string& name, const std::string& val, T& ref,
                bool use_querywithparser)
{
    auto& recursive_symbols = g_parser_recursive_symbols[OpenMP::get_thread_num()];
    if (auto found = recursive_symbols.find(name); found != recursive_symbols.end()) {
        amrex::Error("ParmParse: recursive reference to "+name+" is not allowed");
        return false;
    } else {
        recursive_symbols.insert(name);
    }

    auto parser = pp_make_parser<T>(val, {}, table, parser_prefix, use_querywithparser);
    auto exe = parser.template compileHost<0>();
    ref = static_cast<T>(exe());

    recursive_symbols.erase(name);
    return true;
}

}  // End of unnamed namespace.

std::string
ParmParse::prefixedName (const std::string_view& str) const
{
    AMREX_ASSERT( ! str.empty() );

    if (m_prefix.empty()) {
        return std::string(str);
    } else {
        std::string r = m_prefix + '.';
        r.append(str);
        return r;
    }
}

void
ParmParse::addfile (std::string const& filename) {
#ifdef AMREX_USE_MPI
    if (ParallelDescriptor::Communicator() == MPI_COMM_NULL)
    {
        throw std::runtime_error("ParmParse::addfile: AMReX must be initialized");
    }
#endif

    auto file = FileKeyword;
    std::vector<std::string> val{{filename}};
    addDefn(file, val, g_table);
}

void
ParmParse::Initialize (int         argc,
                       char**      argv,
                       const char* parfile)
{
    if ( initialized )
    {
        amrex::Error("ParmParse::Initialize(): already initialized!");
    }

    g_parser_recursive_symbols.resize(OpenMP::get_max_threads());

    ppinit(argc, argv, parfile, g_table);

    amrex::ExecOnFinalize(ParmParse::Finalize);
}

bool
ParmParse::QueryUnusedInputs ()
{
    if ( ParallelDescriptor::IOProcessor() && unused_table_entries_q(g_table))
    {
        if (ParmParse::Verbose()) {
            amrex::OutStream() << "Unused ParmParse Variables:\n";
            pp_print_unused("  [TOP]", g_table);
            amrex::OutStream() << '\n';
        }
        return true;
    }
    return false;
}

bool
ParmParse::hasUnusedInputs (const std::string& prefix)
{
    return unused_table_entries_q(g_table, prefix);
}

std::vector<std::string>
ParmParse::getUnusedInputs (const std::string& prefix)
{
    std::vector<std::string> sorted_names;
    const std::string prefixdot = prefix.empty() ? std::string() : prefix+".";
    for (auto const& [name, entry] : g_table) {
        if (entry.m_count == 0 &&
            name.substr(0,prefixdot.size()) == prefixdot)
        {
            sorted_names.push_back(name);
        }
    }
    std::sort(sorted_names.begin(), sorted_names.end());

    std::vector<std::string> r;
    for (auto const& name : sorted_names) {
        auto const& entry = g_table[name];
        for (auto const& vals : entry.m_vals) {
            std::string tmp(name);
            tmp.append(" =");
            for (auto const& v : vals) {
                tmp += " " + v;
            }
            r.emplace_back(std::move(tmp));
        }
    }

    return r;
}

std::set<std::string>
ParmParse::getEntries (const std::string& prefix)
{
    std::set<std::string> r;
    const std::string prefixdot = prefix.empty() ? std::string() : prefix+".";
    for (auto const& [name, entry] : g_table) {
        if (name.substr(0,prefixdot.size()) == prefixdot) {
            r.insert(name);
        }
    }
    return r;
}

int
ParmParse::Verbose ()
{
    if (pp_detail::verbose < 0) {
        pp_detail::verbose = std::max(amrex::Verbose(),0);
        ParmParse pp("amrex.parmparse");
        if (! pp.query("verbose", "v", pp_detail::verbose)) {
            pp.add("verbose", pp_detail::verbose);
        }
    }
    return pp_detail::verbose;
}

void
ParmParse::SetVerbose (int v)
{
    pp_detail::verbose = v;
}

void
ParmParse::Finalize ()
{
    if ( ParallelDescriptor::IOProcessor() && unused_table_entries_q(g_table))
    {
        if (ParmParse::Verbose()) {
            amrex::OutStream() << "Unused ParmParse Variables:\n";
            pp_print_unused("  [TOP]", g_table);
            amrex::OutStream() << '\n';
        }
        if (amrex::system::abort_on_unused_inputs) {
            amrex::Abort("ERROR: unused ParmParse variables.");
        }
    }
    g_table.clear();

#if !defined(BL_NO_FORT)
    amrex_finalize_namelist();
#endif

    g_parser_recursive_symbols.clear();

    pp_detail::verbose = -1;
    initialized = false;
}

void
ParmParse::SetParserPrefix (std::string a_prefix)
{
    ParmParse::ParserPrefix = std::move(a_prefix);
}

void
ParmParse::dumpTable (std::ostream& os, bool prettyPrint)
{
    std::vector<std::string> sorted_names;
    sorted_names.reserve(g_table.size());
    for (auto const& [name, entry] : g_table) {
        sorted_names.push_back(name);
    }
    std::sort(sorted_names.begin(), sorted_names.end());

    for (auto const& name : sorted_names) {
        auto const& entry = g_table[name];
        if (prettyPrint && entry.m_count > 0) {
            for (auto const& vals : entry.m_vals) {
                os << pp_to_pretty_string(name, vals) << '\n';
            }
        }
        else {
            for (auto const& vals : entry.m_vals) {
                os << pp_to_string(name, vals) << '\n';
            }
        }
    }
}

void
ParmParse::prettyPrintTable (std::ostream& os)
{
    std::vector<std::string> sorted_names;
    sorted_names.reserve(g_table.size());
    for (auto const& [name, entry] : g_table) {
        sorted_names.push_back(name);
    }
    std::sort(sorted_names.begin(), sorted_names.end());

    for (auto const& name : sorted_names) {
        auto const& entry = g_table[name];
        std::vector<std::string> value_string;
        std::unordered_map<std::string,int> count;
        for (auto const& vals : entry.m_vals) {
            value_string.emplace_back(pp_to_pretty_string(name, vals));
            ++count[value_string.back()];
        }
        for (auto const& s : value_string) {
            if (--count[s] == 0) {
                os << s << '\n';
            }
        }
    }
}

int
ParmParse::countval (const char* name,
                     int         n) const
{
    //
    // First find n'th occurrence of name in table.
    //
    auto const* def = ppindex(*m_table, n, prefixedName(name));
    return def == nullptr ? 0 : static_cast<int>(def->size());
}

// BOOL
void
ParmParse::getkth (const char* name,
                   int         k,
                   bool&       ref,
                   int         ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name,
                bool&       ref,
                int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     bool&       ref,
                     int         ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name,
                  bool&       ref,
                  int         ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, // NOLINT(readability-make-member-function-const)
                const bool  val)
{
    saddval(prefixedName(name),val);
}

// INT
void
ParmParse::getkth (const char* name, int k, int& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, int& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, int& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, int& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const int val) // NOLINT(readability-make-member-function-const)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<int>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<int>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<int>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<int>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<int>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}


// LONG
void
ParmParse::getkth (const char* name, int k, long& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, long& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, long& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, long& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, // NOLINT(readability-make-member-function-const)
                const long  val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<long>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<long>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<long>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<long>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<long>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}

// long long
void
ParmParse::getkth (const char* name, int k, long long& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, long long& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, long long& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, long long& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const long long val) // NOLINT(readability-make-member-function-const)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<long long>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<long long>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<long long>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<long long>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<long long>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}

// FLOAT
void
ParmParse::getkth (const char* name, int k, float& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, float& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, float& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, float& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const float val) // NOLINT(readability-make-member-function-const)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<float>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<float>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<float>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<float>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<float>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}



// DOUBLE
void
ParmParse::getkth (const char* name, int k, double& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, double& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, double& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, double& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const double val) // NOLINT(readability-make-member-function-const)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<double>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<double>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<double>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<double>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<double>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}



// STRING
void
ParmParse::getkth (const char* name, int k, std::string& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, std::string& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, std::string& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, std::string& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const std::string& val) // NOLINT(readability-make-member-function-const)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<std::string>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<std::string>& ref,
                   int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<std::string>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<std::string>& ref,
                     int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<std::string>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}



// INTVECT
void
ParmParse::getkth (const char* name, int k, IntVect& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, IntVect& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, IntVect& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, IntVect& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const IntVect& val) // NOLINT(readability-make-member-function-const)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<IntVect>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<IntVect>& ref,
                   int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<IntVect>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<IntVect>& ref,
                     int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<IntVect>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}

// BOX
void
ParmParse::getkth (const char* name, int k, Box& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, Box& ref, int ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, Box& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, Box& ref, int ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const Box& val) // NOLINT(readability-make-member-function-const)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<Box>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<Box>& ref,
                   int start_ix, int num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<Box>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<Box>& ref,
                     int start_ix, int num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<Box>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}


int
ParmParse::queryarr (const char* name, IntVect& ref) const
{
    std::vector<int> v;
    int exist = this->queryarr(name, v);
    if (exist) {
        AMREX_ALWAYS_ASSERT(v.size() == AMREX_SPACEDIM);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) { ref[i] = v[i]; }
    }
    return exist;
}

void
ParmParse::getarr (const char* name, IntVect& ref) const
{
    std::vector<int> v;
    this->getarr(name, v);
    AMREX_ALWAYS_ASSERT(v.size() == AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) { ref[i] = v[i]; }
}

int
ParmParse::queryarr (const char* name, RealVect& ref) const
{
    std::vector<Real> v;
    int exist = this->queryarr(name, v);
    if (exist) {
        AMREX_ALWAYS_ASSERT(v.size() == AMREX_SPACEDIM);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) { ref[i] = v[i]; }
    }
    return exist;
}

void
ParmParse::getarr (const char* name, RealVect& ref) const
{
    std::vector<Real> v;
    this->getarr(name, v);
    AMREX_ALWAYS_ASSERT(v.size() == AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) { ref[i] = v[i]; }
}

//
// Return number of occurrences of parameter name.
//

int
ParmParse::countname (const std::string& name) const
{
    auto pname = prefixedName(name);
    auto found = m_table->find(pname);
    if (found != m_table->cend()) {
        return static_cast<int>(found->second.m_vals.size());
    } else {
        return 0;
    }
}

//
// Return true if name in table.
//

bool
ParmParse::contains (const char* name) const
{
    auto pname = prefixedName(name);
    auto found = m_table->find(pname);
    if (found != m_table->cend()) {
#ifdef AMREX_USE_OMP
#pragma omp atomic update
#endif
        ++(found->second.m_count);
        return true;
    } else {
        return false;
    }
}

int
ParmParse::remove (const char* name)
{
    auto const pname = prefixedName(name);
    auto n = m_table->erase(pname);
    return static_cast<int>(n);
}

namespace {
template <class T>
bool squeryWithParser (const ParmParse::Table& table,
                       const std::string&      parser_prefix,
                       const std::string&      name,
                       T&                      ref)
{
    std::vector<std::string> vals;
    bool exist = squeryarr(table, parser_prefix, name, vals,
                           ParmParse::FIRST, ParmParse::ALL, ParmParse::LAST);
    if (!exist) { return false; }

    std::string combined_string;
    for (auto const& v : vals) {
        combined_string.append(v);
    }
    return pp_parser(table, parser_prefix, name, combined_string, ref, true);
}

template <class T>
bool squeryarrWithParser (const ParmParse::Table& table,
                          const std::string&      parser_prefix,
                          const std::string&      name,
                          int                     nvals,
                          T*                      ref)
{
    std::vector<std::string> vals;
    bool exist = squeryarr(table, parser_prefix, name, vals,
                           ParmParse::FIRST, ParmParse::ALL, ParmParse::LAST);
    if (!exist) { return false; }

    AMREX_ALWAYS_ASSERT(int(vals.size()) == nvals);
    for (int ival = 0; ival < nvals; ++ival) {
        bool r = pp_parser(table, parser_prefix, name, vals[ival], ref[ival], true);
        if (!r) { return false; }
    }
    return true;
}
}

int
ParmParse::queryWithParser (const char* name, int& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (const char* name, long& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (const char* name, long long& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (const char* name, float& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (const char* name, double& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryarrWithParser (const char* name, int nvals, int* ref) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ref);
}

int
ParmParse::queryarrWithParser (const char* name, int nvals, long* ref) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ref);
}

int
ParmParse::queryarrWithParser (const char* name, int nvals, long long* ref) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ref);
}

int
ParmParse::queryarrWithParser (const char* name, int nvals, float* ref) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ref);
}

int
ParmParse::queryarrWithParser (const char* name, int nvals, double* ref) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ref);
}

Parser
ParmParse::makeParser (std::string const& func,
                       Vector<std::string> const& vars) const
{
    return pp_make_parser<double>(func, vars, *m_table, m_parser_prefix, true);
}

IParser
ParmParse::makeIParser (std::string const& func,
                        Vector<std::string> const& vars) const
{
    return pp_make_parser<long long>(func, vars, *m_table, m_parser_prefix, true);
}

}
