#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>
#include <AMReX_BLFort.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iterator>
#include <list>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

extern "C" void amrex_init_namelist (const char*);
extern "C" void amrex_finalize_namelist ();

namespace amrex {

#ifdef AMREX_XSDK
static bool finalize_verbose = false;
#else
static bool finalize_verbose = true;
#endif

std::string const ParmParse::FileKeyword = "FILE";

//
// Used by constructor to build table.
//
ParmParse::PP_entry::PP_entry (std::string name, const std::list<std::string>& vals)
    :
    m_name(std::move(name)),
    m_vals(vals.begin(), vals.end()),
    m_table(nullptr),
    m_queried(false)
{
}

ParmParse::PP_entry::PP_entry (std::string name, const std::string& val)
    :
    m_name(std::move(name)),
    m_vals({val}),
    m_table(nullptr),
    m_queried(false)
{
}

ParmParse::PP_entry::PP_entry (std::string name, const std::list<PP_entry>& table)
    :
    m_name(std::move(name)),
    m_table(new Table(table)),
    m_queried(false)
{
}

ParmParse::PP_entry::PP_entry (const PP_entry& pe)
    : m_name(pe.m_name),
      m_vals(pe.m_vals),
      m_table(nullptr),
      m_queried(pe.m_queried)
{
    if ( pe.m_table )
    {
        m_table = new Table(*pe.m_table);
    }
}

ParmParse::PP_entry::~PP_entry ()
{
    delete m_table;
}

ParmParse::PP_entry&
ParmParse::PP_entry::operator= (const PP_entry& pe)
{
    if ( &pe == this ) return *this;
    m_name = pe.m_name;
    m_vals = pe.m_vals;
    m_table = nullptr;
    m_queried = pe.m_queried;
    if ( pe.m_table )
    {
        m_table = new Table(*pe.m_table);
    }
    return *this;
}

std::string
ParmParse::PP_entry::print () const {
    std::stringstream t;
    t << m_name << " = ";
    int n = static_cast<int>(m_vals.size());
    for ( int i = 0; i < n; i++)
    {
        t << m_vals[i];
        if ( i < n-1 ) t << " ";
    }
    return t.str();
}

std::ostream&
operator<< (std::ostream& os, const ParmParse::PP_entry& pp)
{
    os << pp.m_name << "(nvals = " << pp.m_vals.size() << ") " << " :: [";
    int n = static_cast<int>(pp.m_vals.size());
    for ( int i = 0; i < n; i++ )
    {
        os << pp.m_vals[i];
        if ( i < n-1 ) os << ", ";
    }
    os << "]";

    if ( !os )
    {
        amrex::Error("write on ostream failed");
    }
    return os;
}

namespace
{
enum PType
{
    pDefn,
    pValue,
    pEQ_sign,
    pOpenBracket,
    pCloseBracket,
    pEOF
};

template <class T>
bool
isT (const std::string& str, T& val)
{
    std::istringstream s(str);
    s >> val;
    if ( s.fail() ) return false;
    std::string left;
    std::getline(s, left);
    if ( !left.empty() ) return false;
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
    if ( str == "true" || str == "t" )
    {
        val = true;
        return true;
    }
    if ( str == "false" || str == "f" )
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

ParmParse::Table g_table;

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

void
eat_garbage (const char*& str)
{
    for (;;)
    {
        if ( *str == 0 ) break; // NOLINT
        else if ( *str == '#' )
        {
            while ( *str && *str != '\n' )
            {
                str++;
            }
            continue;
        }
        else if ( isspace(*str) )
        {
            str++;
        }
        else
        {
            break;
        }
    }
}

PType
getToken (const char*& str,
          std::string& ostr)
{
#define ERROR_MESS                                                             \
   amrex::ErrorStream() << "ParmParse::getToken(): invalid string = " << ostr << '\n'; \
   amrex::ErrorStream() << "STATE = " << state_name[state]              \
                        << ", next char = " << ch << '\n';              \
   amrex::ErrorStream() << ", rest of input = \n" << str << '\n';       \
   amrex::Abort()
   //
   // Eat white space and comments.
   //
   eat_garbage(str);
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
           else if ( ch == '{' )
           {
               str++;
               return pOpenBracket;
           }
           else if ( ch == '}' )
           {
               str++;
               return pCloseBracket;
           }
           else if ( isalpha(ch) )
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
           if ( isalnum(ch) || ch == '_' || ch == '.' || ch == '[' || ch == ']' || ch == '+' || ch == '-' )
           {
               ostr += ch; str++;
           }
           else if ( isspace(ch) || ch == '=' )
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
           if ( isspace(ch) || ch == '=' )
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
           ERROR_MESS;
       }
   }
#undef ERROR_MESS
}


//
// Keyword aware string comparison.
//


bool
ppfound (const std::string& keyword,
         const ParmParse::PP_entry& pe,
         bool recordQ)
{
    return (recordQ == (pe.m_table!=nullptr)) && (keyword == pe.m_name);
}

//
// Return the index of the n'th occurrence of a parameter name,
// except if n==-1, return the index of the last occurrence.
// Return 0 if the specified occurrence does not exist.
//

const ParmParse::PP_entry*
ppindex (const ParmParse::Table& table,
         int         n,
         const std::string& name,
         bool recordQ)
{
    const ParmParse::PP_entry* fnd = nullptr;

    if ( n == ParmParse::LAST )
    {
        //
        // Search from back of list.
        //
        for (auto li = table.crbegin(), REnd = table.crend(); li != REnd; ++li)
        {
            if ( ppfound(name, *li, recordQ) )
            {
                fnd = &*li;
                break;
            }
        }
    }
    else
    {
        for (auto const& li : table)
        {
            if ( ppfound(name, li, recordQ) )
            {
                fnd = &li;
                if ( --n < 0 )
                {
                    break;
                }
            }
        }
        if ( n >= 0)
        {
            fnd = nullptr;
        }
    }

    if ( fnd )
    {
        //
        // Found an entry; mark all occurrences of name as used.
        //
        for (const auto& li : table)
        {
            if ( ppfound(name, li, recordQ) )
            {
                li.m_queried = true;
            }
        }
    }
    return fnd;
}

void
bldTable (const char*& str, std::list<ParmParse::PP_entry>& tab);

namespace {
    bool isTrue(std::smatch const& sm) {
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
}

void
read_file (const char* fname, std::list<ParmParse::PP_entry>& tab)
{
    //
    // Space for input file if it exists.
    //
    if ( fname != nullptr && fname[0] != 0 )
    {
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
addDefn (std::string&         def,
         std::list<std::string>&   val,
         std::list<ParmParse::PP_entry>& tab)
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
        tab.emplace_back(def,val);
    }
    val.clear();
    if ( def != ParmParse::FileKeyword )
        def = std::string();
}

void
addTable (std::string& def,
          ParmParse::Table& val,
          std::list<ParmParse::PP_entry>& tab)
{
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
        amrex::ErrorStream() << "ParmParse::addTable(): no values for Table " << def << "\n";
        amrex::Abort();
    }
    tab.emplace_back(def, val);
    val.clear();
    def = std::string();
}

void
bldTable (const char*&                    str,
          std::list<ParmParse::PP_entry>& tab)
{
    std::string            cur_name;
    std::list<std::string> cur_list;
    ParmParse::Table       cur_table;
    std::string            tmp_str;

    for (;;)
    {
        std::string tokname;

        PType token = getToken(str,tokname);

        switch (token)
        {
        case pCloseBracket:
            if ( !cur_name.empty() && cur_list.empty() )
            {
                amrex::Abort("ParmParse::bldTable() defn with no list");
            }
            AMREX_FALLTHROUGH;
        case pEOF:
            addDefn(cur_name,cur_list,tab);
            return;
        case pOpenBracket:
            if ( cur_name.empty() )
            {
                amrex::Abort("ParmParse::bldTabe() '{' with no blocknamne");
            }
            if ( !cur_list.empty() )
            {
                tmp_str = cur_list.back();
                cur_list.pop_back();
                addDefn(cur_name, cur_list, tab);
                cur_name = tmp_str;
            }
            bldTable(str, cur_table);
            addTable(cur_name, cur_table, tab);
            break;
        case pEQ_sign:
            if ( cur_name.empty() )
            {
                amrex::Abort("ParmParse::bldTable() EQ with no current defn");
            }
            if ( !cur_list.empty() )
            {
                tmp_str = cur_list.back();
                cur_list.pop_back();
                addDefn(cur_name,cur_list,tab);
                cur_name = tmp_str;
            }
            //
            // Read one too far, remove last name on list.
            //
            break;
        case pDefn:
            if ( cur_name.empty() )
            {
                cur_name = tokname;
                break;
            }
            //
            // Otherwise, fall through, this may be a string.
            //
            AMREX_FALLTHROUGH;
        case pValue:
            if ( cur_name.empty() )
            {
                std::string msg("ParmParse::bldTable(): value with no defn: ");
                msg += tokname;
                amrex::Abort(msg.c_str());
            }
            cur_list.push_back(tokname);
            break;
        }
    }
}

namespace
{
template <class T>
bool
squeryval (const ParmParse::Table& table,
           const std::string& name,
           T&           ref,
           int          ival,
           int          occurrence)
{
    //
    // Get last occurrence of name in table.
    //
    const ParmParse::PP_entry* def = ppindex(table, occurrence, name, false);
    if ( def == nullptr )
    {
        return false;
    }
    //
    // Does it have ival values?
    //
    if ( ival >= static_cast<int>(def->m_vals.size()) )
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
        amrex::ErrorStream() << def->m_name << '\n' << *def << '\n';
        amrex::Abort();
    }

    const std::string& valname = def->m_vals[ival];

    bool ok = is(valname, ref);
    if ( !ok )
    {
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
        amrex::ErrorStream() << def->m_name << '\n';
        amrex::ErrorStream() << " Expected an \""
                  << tok_name(ref)
                  << "\" type  which can't be parsed from the string \""
                  << valname << "\"\n"
                  << *def << '\n';
        amrex::Abort();
    }
    return true;
}

template <class T>
void
sgetval (const ParmParse::Table& table,
         const std::string& name,
         T&           ref,
         int          ival,
         int          occurrence)
{
    if ( squeryval(table, name,ref,ival,occurrence) == 0 )
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
           const std::string& name,
           std::vector<T>&    ref,
           int          start_ix,
           int          num_val,
           int          occurrence)
{
    //
    // Get last occurrence of name in table.
    //
    const ParmParse::PP_entry *def = ppindex(table,occurrence, name, false);
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
        num_val = static_cast<int>(def->m_vals.size());
    }

    if ( num_val == 0 ) return true;

    int stop_ix = start_ix + num_val - 1;
    if ( static_cast<int>(ref.size()) <= stop_ix )
    {
        ref.resize(stop_ix + 1);
    }
    if ( stop_ix >= static_cast<int>(def->m_vals.size()) )
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
        amrex::ErrorStream() << def->m_name << '\n' << *def << '\n';
        amrex::Abort();
    }
    for ( int n = start_ix; n <= stop_ix; n++ )
    {
        const std::string& valname = def->m_vals[n];
        bool ok = is(valname, ref[n]);
        if ( !ok )
        {
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
            amrex::ErrorStream() << def->m_name << '\n';
            amrex::ErrorStream() << " Expected an \""
                      << tok_name(ref)
                      << "\" type which can't be parsed from the string \""
                      << valname << "\"\n"
                      << *def << '\n';
            amrex::Abort();
        }
    }
    return true;
}

template <class T>
void
sgetarr (const ParmParse::Table& table, const std::string& name,
         std::vector<T>& ref, int start_ix, int num_val, int occurrence)
{
    if ( squeryarr(table,name,ref,start_ix,num_val,occurrence) == 0 )
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
    ParmParse::PP_entry entry(name,val.str());
    entry.m_queried=true;
    g_table.push_back(entry);
}


template <class T>
void
saddarr (const std::string& name, const std::vector<T>& ref)
{
    std::list<std::string> arr;
    for (auto const& item : ref) {
        std::stringstream val;
        val << std::setprecision(17) << item;
        arr.push_back(val.str());
    }
    ParmParse::PP_entry entry(name,arr);
    entry.m_queried=true;
    g_table.push_back(entry);
}

}

//
// Initialize ParmParse.
//

namespace {
    bool initialized = false;
}

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
        std::list<ParmParse::PP_entry> arg_table;
        const char* b = argstr.c_str();
        bldTable(b, arg_table);
        //
        // Append arg_table to end of existing table.
        //
        g_table.splice(table.end(), arg_table);
    }
    initialized = true;
}

}  // End of unnamed namespace.

std::string
ParmParse::prefixedName (const std::string& str) const
{
    if ( str.empty() )
    {
        amrex::Error("ParmParse::prefixedName: has empty name");
    }
    if ( !m_pstack.top().empty())
    {
        return m_pstack.top() + '.' + str;
    }
    return str;
}

void
ParmParse::addfile (std::string const& filename) {
    auto l = std::list<std::string>{filename};
    auto file = FileKeyword;
    addDefn(file,
            l,
            g_table);
}

void
ParmParse::pushPrefix (const std::string& str)
{
    std::string s(str);
    if ( !s.empty() )
    {
        if ( !m_pstack.top().empty() )
        {
            s = m_pstack.top() + "." + s;
        }
        m_pstack.push(s);
    }
}

void
ParmParse::popPrefix ()
{
    if ( m_pstack.size() <= 1 )
    {
        amrex::Error("ParmParse::popPrefix: stack underflow");
    }
    m_pstack.pop();
}

std::string
ParmParse::getPrefix() const
{
    return m_pstack.top();
}

ParmParse::ParmParse (const std::string& prefix)
    :
    m_table(g_table)
{
    m_pstack.push(prefix);
}

ParmParse::ParmParse (Table& a_table)
    : m_table(a_table)
{
    m_pstack.push("");
}

ParmParse::Frame::Frame (ParmParse& pp, const std::string& pfix)
    :
    m_pp(pp)
{
    push(pfix);
    BL_ASSERT( m_np == 1 );
}

ParmParse::Frame::~Frame ()
{
    BL_ASSERT( m_np > 0 );
    while ( m_np )
    {
        pop();
    }
    BL_ASSERT( m_np == 0 );
}

void
ParmParse::Frame::push (const std::string& str)
{
    m_pp.pushPrefix(str);
    m_np++;
}

void
ParmParse::Frame::pop ()
{
    BL_ASSERT( m_np > 0);
    m_pp.popPrefix();
    m_np--;
}

std::string
ParmParse::Frame::getPrefix () const
{
    return m_pp.getPrefix();
}

void
ParmParse::appendTable(ParmParse::Table& tab)
{
  g_table.splice(g_table.end(), tab);
}

static
bool
unused_table_entries_q (const ParmParse::Table& table, const std::string& prefix = std::string())
{
    for (auto const& li : table) // NOLINT(readability-use-anyofallof)
    {
        if ( li.m_table )
        {
            if ( !li.m_queried )
            {
                if (prefix.empty()) {
                    return true;
                } else {
                    if (li.m_name.substr(0,prefix.size()+1) == prefix+".") {
                        return true;
                    }
                }
            }
            else
            {
                if (unused_table_entries_q(*li.m_table, prefix)) return true;
            }
        }
        else if ( !li.m_queried )
        {
            if (prefix.empty()) {
                return true;
            } else {
                if (li.m_name.substr(0,prefix.size()+1) == prefix+".") {
                    return true;
                }
            }
        }
    }
    return false;
}

static
void
finalize_table (const std::string& pfx, const ParmParse::Table& table)
{
    for (auto const& li : table)
    {
        if ( li.m_table )
        {
            if ( !li.m_queried )
            {
                if (finalize_verbose) {
                    amrex::AllPrint() << "Record " << li.m_name << std::endl;
                }
            }
            else
            {
                finalize_table(pfx + "::" + li.m_name, *li.m_table);
            }
        }
        else if ( !li.m_queried )
        {
            if (finalize_verbose) {
                amrex::AllPrint() << pfx << "::" << li << std::endl;
            }
        }
    }
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

    ppinit(argc, argv, parfile, g_table);

    amrex::ExecOnFinalize(ParmParse::Finalize);
}

bool
ParmParse::QueryUnusedInputs ()
{
    if ( ParallelDescriptor::IOProcessor() && unused_table_entries_q(g_table))
    {
      finalize_verbose = amrex::system::verbose;
      if (finalize_verbose) amrex::OutStream() << "Unused ParmParse Variables:\n";
      finalize_table("  [TOP]", g_table);
      if (finalize_verbose) amrex::OutStream() << std::endl;
      return true;
    }
    return false;
}

bool
ParmParse::hasUnusedInputs (const std::string& prefix)
{
    return unused_table_entries_q(g_table, prefix);
}

static
void
get_entries_under_prefix (std::vector<std::string>& found_entries,
                         const ParmParse::Table& table,
                         const std::string& prefix,
                         const bool only_unused = false,
                         const bool add_values = false)
{
    const std::string prefixdot = prefix.empty() ? std::string() : prefix+".";
    for (auto const& entry : table) {
        if ((! only_unused) || (only_unused && ! entry.m_queried)) {
            if (entry.m_name.substr(0,prefixdot.size()) == prefixdot) {
                std::string tmp(entry.m_name);
                if (add_values) {
                    tmp.append(" =");
                    for (auto const& v : entry.m_vals) {
                        tmp += " " + v;
                    }
                }
                found_entries.emplace_back(std::move(tmp));
            }
        }

        if (entry.m_table) {
            get_entries_under_prefix(found_entries, table, prefix,
                                     only_unused, add_values);
        }
    }
}

std::vector<std::string>
ParmParse::getUnusedInputs (const std::string& prefix)
{
    std::vector<std::string> r;
    get_entries_under_prefix(r, g_table, prefix, true, true);
    return r;
}

std::set<std::string>
ParmParse::getEntries (const std::string& prefix)
{
    std::vector<std::string> r;
    get_entries_under_prefix(r, g_table, prefix, false, false);
    return std::set<std::string>(r.begin(), r.end());
}

void
ParmParse::Finalize ()
{
    if ( ParallelDescriptor::IOProcessor() && unused_table_entries_q(g_table))
    {
      finalize_verbose = amrex::system::verbose;
      if (finalize_verbose) amrex::OutStream() << "Unused ParmParse Variables:\n";
      finalize_table("  [TOP]", g_table);
      if (finalize_verbose) amrex::OutStream() << std::endl;
      //
      // First loop through and delete all queried entries.
      //
      if (amrex::system::abort_on_unused_inputs) amrex::Abort("ERROR: unused ParmParse variables.");
    }
    g_table.clear();

#if !defined(BL_NO_FORT)
    amrex_finalize_namelist();
#endif

    initialized = false;
}

void
ParmParse::dumpTable (std::ostream& os, bool prettyPrint)
{
    for (auto const& li : g_table)
    {
        if(prettyPrint && li.m_queried) {
            os << li.print() << std::endl;
        }
        else {
            os << li << std::endl;
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
    const PP_entry* def = ppindex(m_table, n, prefixedName(name), false);
    return def == nullptr ? 0 : static_cast<int>(def->m_vals.size());
}

// BOOL
void
ParmParse::getkth (const char* name,
                   int         k,
                   bool&       ref,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name,
                bool&       ref,
                int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     bool&       ref,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name,
                  bool&       ref,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name,
                const bool  val)
{
    saddval(prefixedName(name),val);
}

// INT
void
ParmParse::getkth (const char* name, int k, int& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, int& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, int& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, int& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const int val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<int>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<int>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<int>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<int>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<int>& ref)
{
    saddarr(prefixedName(name),ref);
}


// LONG
void
ParmParse::getkth (const char* name, int k, long& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, long& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, long& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, long& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name,
                const long  val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<long>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<long>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<long>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<long>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<long>& ref)
{
    saddarr(prefixedName(name),ref);
}

// long long
void
ParmParse::getkth (const char* name, int k, long long& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, long long& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, long long& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, long long& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const long long val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<long long>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<long long>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<long long>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<long long>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<long long>& ref)
{
    saddarr(prefixedName(name),ref);
}

// FLOAT
void
ParmParse::getkth (const char* name, int k, float& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, float& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, float& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, float& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const float val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<float>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<float>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<float>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<float>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<float>& ref)
{
    saddarr(prefixedName(name),ref);
}



// DOUBLE
void
ParmParse::getkth (const char* name, int k, double& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, double& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, double& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, double& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const double val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<double>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<double>& ref, int start_ix,
                   int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<double>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<double>& ref, int start_ix,
                     int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<double>& ref)
{
    saddarr(prefixedName(name),ref);
}



// STRING
void
ParmParse::getkth (const char* name, int k, std::string& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, std::string& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, std::string& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, std::string& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const std::string& val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<std::string>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<std::string>& ref,
                   int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<std::string>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<std::string>& ref,
                     int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<std::string>& ref)
{
    saddarr(prefixedName(name),ref);
}



// INTVECT
void
ParmParse::getkth (const char* name, int k, IntVect& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, IntVect& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, IntVect& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, IntVect& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const IntVect& val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<IntVect>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<IntVect>& ref,
                   int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<IntVect>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<IntVect>& ref,
                     int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<IntVect>& ref)
{
    saddarr(prefixedName(name),ref);
}

// BOX
void
ParmParse::getkth (const char* name, int k, Box& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (const char* name, Box& ref, int ival) const
{
    sgetval(m_table, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (const char* name, int k, Box& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (const char* name, Box& ref, int ival) const
{
    return squeryval(m_table, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (const char* name, const Box& val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (const char* name, int k, std::vector<Box>& ref,
                      int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name, std::vector<Box>& ref,
                   int start_ix, int num_val) const
{
    sgetarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name, int k, std::vector<Box>& ref,
                        int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char* name, std::vector<Box>& ref,
                     int start_ix, int num_val) const
{
    return squeryarr(m_table, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (const char* name, const std::vector<Box>& ref)
{
    saddarr(prefixedName(name),ref);
}


//
// Return number of occurrences of parameter name.
//

int
ParmParse::countname (const std::string& name) const
{
    int cnt = 0;
    for (auto const& li : m_table)
    {
        if ( ppfound(prefixedName(name), li, false) )
        {
            cnt++;
        }
    }
    return cnt;
}

int
ParmParse::countRecords (const std::string& name) const
{
    int cnt = 0;
    for (auto const& li : m_table)
    {
        if ( ppfound(prefixedName(name), li, true) )
        {
            cnt++;
        }
    }
    return cnt;
}

//
// Return true if name in table.
//

bool
ParmParse::contains (const char* name) const
{
    for (auto const& li : m_table)
    {
       if ( ppfound(prefixedName(name), li, false))
       {
           //
           // Found an entry; mark all occurrences of name as used.
           //
           for (auto& lli : m_table)
           {
               if ( ppfound(prefixedName(name), lli, false) )
               {
                   lli.m_queried = true;
               }
           }
           return true;
       }
    }
    return false;
}

int
ParmParse::remove (const char* name)
{
    int r = 0;
    for (auto it = m_table.begin(); it != m_table.end(); ) {
        if (ppfound(prefixedName(name), *it, false)) {
            it = m_table.erase(it);
            ++r;
        } else {
            ++it;
        }
    }
    return r;
}

ParmParse::Record
ParmParse::getRecord (const std::string& name, int n) const
{
    const PP_entry* pe = ppindex(m_table, n, prefixedName(name), true);
    if ( pe == nullptr )
    {
        amrex::ErrorStream() << "ParmParse::getRecord: record " << name << " not found" << std::endl;
        amrex::Abort();
        return Record(ParmParse());
    } else {
        return Record(ParmParse(*pe->m_table));
    }
}

//
//
//

ParmParse::Record::Record ( const ParmParse& pp )
    : m_pp(pp)
{
}

const ParmParse*
ParmParse::Record::operator-> () const
{
    return &m_pp;
}

const ParmParse&
ParmParse::Record::operator* () const
{
    return m_pp;
}

}
