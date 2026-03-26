#include <AMReX_ParmParse.H>
#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_RealVect.H>
#include <AMReX_String.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
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
    std::vector<std::set<std::string>> g_parser_recursive_symbols(1);
    std::string g_toml_table_key;
    /// \cond DOXYGEN_IGNORE
    namespace pp_detail {
        int verbose = -1;
    }
    /// \endcond
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

bool is_toml_array (std::string const& token);

std::string pp_to_pretty_string (std::string const& name,
                                 std::vector<std::string> const& vals,
                                 std::vector<ParmParse::QuoteType> const* quotes,
                                 ParmParse::PP_entry const* entry)
{
    std::stringstream ss;
    ss << name << " =";
    for (std::size_t i = 0; i < vals.size(); ++i) {
        auto const& v = vals[i];
        if (is_toml_array(v)) {
            ss << " " << v.substr(5); // length of $$ARR is 5.
        } else {
            ParmParse::QuoteType quote = ParmParse::QuoteType::None;
            if (quotes != nullptr && i < quotes->size()) {
                quote = (*quotes)[i];
            }
            switch (quote) {
            case ParmParse::QuoteType::Double:
                ss << " \"" << v << "\"";
                break;
            case ParmParse::QuoteType::Triple:
                ss << R"( """)" << v << R"(""")";
                break;
            case ParmParse::QuoteType::None:
            default:
                ss << " " << v;
                break;
            }
        }
    }
    if (entry && entry->m_parsed && ! entry->m_last_vals.empty()) {
        int min_col = 36;
        int pad = min_col - static_cast<int>(ss.str().size());
        if (pad > 0) {
            ss << std::string(pad, ' ');
        }
        ss.precision(17);
        ss << "    #";
        for (auto const& x : entry->m_last_vals) {
            std::visit([&] (auto&& arg)
            {
                ss << " " << arg;
            }, x);
        }
    }
    return ss.str();
}

std::string pp_to_string (std::string const& name,
                          std::vector<std::string> const& vals)
{
    std::stringstream ss;
    ss << name << "(nvals = " << vals.size() << ") " << " :: [";
    if (vals.size() == 1 && is_toml_array(vals[0])) {
        ss << vals[0].substr(5); // length of $$ARR is 5.
    } else {
        for (std::size_t i = 0; i < vals.size(); ++i) {
            ss << vals[i];
            if ( i < vals.size()-1 ) { ss << ", "; }
        }
    }
    ss << "]";
    return ss.str();
}

enum class PType
{
    Defn,
    EQ_sign,
    Value,
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
is_literal_bool (const std::string& str, T& val)
{
    auto const lo_str = amrex::toLower(str);
    if ( lo_str == "true" || lo_str == "t" ) {
        val = static_cast<T>(1);
        return true;
    } else if ( lo_str == "false" || lo_str == "f" ) {
        val = static_cast<T>(0);
        return true;
    } else {
        return false;
    }
}

template <class T>
bool
is (const std::string& str, T& val)
{
    if constexpr (std::is_integral_v<T>) {
        if (is_literal_bool(str, val)) {
            return true;
        }
        if (isT(str, val)) {
            return true;
        }
        // Treat 123., 123.0, 123.00 etc. as integer.
        auto dec = str.find('.');
        if (dec == std::string::npos) {
            return false;
        }
        if (dec+1 == str.size()) {
            std::string stripped = str;
            stripped.pop_back();
            return isT(stripped, val);
        }
        auto begin_it = str.begin() + static_cast<std::ptrdiff_t>(dec+1);
        auto end_it = str.end();
        if (!std::all_of(begin_it, end_it, [] (char c) { return c == '0'; })) {
            return false;
        }
        std::string stripped = str;
        stripped.erase(dec);
        return isT(stripped, val);
    } else if constexpr (std::is_floating_point_v<T>) {
        if (is_literal_bool(str, val)) {
            return true;
        }
        return is_floating_point(str, val);
    } else {
        return isT(str, val);
    }
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
    if (is_literal_bool(str, val)) {
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

enum class lexState
{
    START,
    STRING,
    QUOTED_STRING,
    TRIPLY_QUOTED_STRING,
    IDENTIFIER,
    LIST,
    INITIALIZER,
    ARRAY
};

int
eat_garbage (const char*& str, bool* newline_from_comment = nullptr)
{
    int num_linefeeds = 0;
    if (newline_from_comment) { *newline_from_comment = false; }
    for (;;)
    {
        if ( *str == 0 ) { break; } // NOLINT
        else if ( *str == '#' )
        {
            while ( *str && *str != '\n' )
            {
                str++;
            }
            if (*str == '\n') {
                if (newline_from_comment) { *newline_from_comment = true; }
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

void eat_comment (const char*& str)
{
    if ( *str == '#' )
    {
        while ( *str && *str != '\n' )
        {
            str++;
        }
    }
}

PType
getToken (const char*& str, std::string& ostr, int& num_linefeeds,
          bool& newline_from_comment, ParmParse::QuoteType& quote_type)
{
   //
   // Eat white space and comments.
   //
   num_linefeeds = eat_garbage(str, &newline_from_comment);
   quote_type = ParmParse::QuoteType::None;
    //
    // Check for end of file.
    //
    if ( *str == 0 )
    {
        return PType::pEOF;
    }
    //
    // Start token scan.
    //
    lexState state = lexState::START;
    int      pcnt  = 0; // Tracks nested parens
    int      cbcnt = 0; // Tracks nested curly braces
    int      sbcnt = 0; // Tracks nested square brackets
    bool     array_in_string = false;
    bool     array_escape = false;
    while (true)
    {
        char ch = *str;
        if ( ch == 0 )
        {
            amrex::Error("ParmParse::getToken: EOF while parsing");
        }
        switch (state)
        {
        case lexState::START:
            if ( ch == '=' )
            {
                ostr += ch; str++;
                return PType::EQ_sign;
            }
            else if ( ch == '"' )
            {
                if (*(str+1) == '"' && *(str+2) == '"') {
                    str += 3;
                    if ((*str) == '\n') {
                        // A newline immediately following the opening
                        // delimiter will be trimmed.
                        ++str;
                    }
                    state = lexState::TRIPLY_QUOTED_STRING;
                    quote_type = ParmParse::QuoteType::Triple;
                } else {
                    str++;
                    state = lexState::QUOTED_STRING;
                    quote_type = ParmParse::QuoteType::Double;
                }
            }
            else if ( ch == '(' )
            {
                ostr += ch; str++; pcnt = 1;
                state = lexState::LIST;
            }
            else if ( ch == '{' )
            {
                ostr += ch; str++; cbcnt = 1;
                state = lexState::INITIALIZER;
            }
            else if ( ch == '[' )
            {
                ostr += "$$ARR[";
                str++; sbcnt = 1;
                array_in_string = false;
                array_escape = false;
                state = lexState::ARRAY;
            }
            else if ( std::isalpha(ch) )
            {
                ostr += ch; str++;
                state = lexState::IDENTIFIER;
            }
            else
            {
                ostr += ch; str++;
                state = lexState::STRING;
            }
            break;
        case lexState::IDENTIFIER:
            if ( std::isalnum(ch) || ch == '_' || ch == '.' || ch == '[' || ch == ']' || ch == '+' || ch == '-' )
            {
                ostr += ch; str++;
            }
            else if ( std::isspace(ch) || ch == '=' )
            {
                return PType::Defn;
            }
            else
            {
                ostr += ch; str++;
                state = lexState::STRING;
            }
            break;
        case lexState::LIST:
            eat_comment(str);
            ch = *str;
            if ( ch == '(' )
            {
                ostr += ch; str++; pcnt++;
            }
            else if ( ch == ')' )
            {
                ostr += ch; str++; pcnt--;
                if ( pcnt == 0 && cbcnt == 0 && sbcnt == 0 )
                {
                    return PType::Value;
                }
            }
            else
            {
                ostr += ch; str++;
            }
            break;
        case lexState::INITIALIZER:
            eat_garbage(str);
            ch = *str;
            if ( ch == '{' )
            {
                ostr += ch; str++; cbcnt++;
            }
            else if ( ch == '}' )
            {
                ostr += ch; str++; cbcnt--;
                if ( cbcnt == 0 && pcnt == 0 && sbcnt == 0 )
                {
                    return PType::Value;
                }
            }
            else
            {
                ostr += ch; str++;
            }
            break;
        case lexState::ARRAY:
        {
            if (!array_in_string) {
                eat_garbage(str);
                ch = *str;
            } else {
                ch = *str;
            }

            if (array_in_string) {
                ostr += ch; str++;
                if (array_escape) {
                    array_escape = false;
                } else if (ch == '\\') {
                    array_escape = true;
                } else if (ch == '"') {
                    array_in_string = false;
                }
            }
            else if ( ch == '"' )
            {
                array_in_string = true;
                array_escape = false;
                ostr += ch; str++;
            }
            else if ( ch == '[' )
            {
                ostr += ch; str++; sbcnt++;
            }
            else if ( ch == ']' )
            {
                ostr += ch; str++; sbcnt--;
                if ( sbcnt == 0 && pcnt == 0 && cbcnt == 0)
                {
                    return PType::Value;
                }
            }
            else
            {
                ostr += ch; str++;
            }
            break;
        }
        case lexState::STRING:
            if ( std::isspace(ch) || ch == '=' )
            {
                return PType::Value;
            }
            else
            {
                ostr += ch; str++;
            }
            break;
        case lexState::TRIPLY_QUOTED_STRING:
            if ( (ch == '"') && (*(str+1) == '"') && (*(str+2) == '"') )
            {
                str += 3;
                return PType::Value;
            }
            else
            {
                ostr += ch; str++;
            }
            break;
        case lexState::QUOTED_STRING:
            if ( ch == '"' )
            {
                str++;
                return PType::Value;
            }
            else
            {
                ostr += ch; str++;
            }
            break;
        default:
            amrex::ErrorStream() << "ParmParse::getToken(): invalid string = " << ostr << '\n'
                                 << "STATE = " << static_cast<int>(state)
                                 << ", next char = " << ch << '\n'
                                 << ", rest of input = \n" << str << '\n';
            amrex::Abort();
        }
    }
}

std::string is_valid_table_key (std::string const& str)
{
    if (str.size() >= 8 && str.substr(0,6) == "$$ARR[" && str.back() == ']') {
        auto key = str.substr(6, str.size()-7);
        bool r = std::isalpha(key[0]);
        for (std::size_t i = 1; i < key.size() && r; ++i) {
            char ch = key[i];
            r = std::isalnum(ch) || ch == '_' || ch == '.' || ch == '-' || ch == '"';
        }
        if (r) { return key; }
    }
    return {};
}

//
// Return the index of the n'th occurrence of a parameter name,
// except if n==-1, return the index of the last occurrence.
// Return nullptr if the specified occurrence does not exist.
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
        if(found->second.m_vals.size() < (std::size_t)n + 1) {
            return nullptr;
        }
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
    } else if (op == "!=") {
        return AMREX_SPACEDIM != dim;
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
        std::string filename = fname;

        // optional prefix to search files in
        char const *amrex_inputs_file_prefix_c = std::getenv("AMREX_INPUTS_FILE_PREFIX");
        if (amrex_inputs_file_prefix_c != nullptr) {
            // we expect a directory path as the prefix: append a trailing "/" if missing
            auto amrex_inputs_file_prefix = std::string(amrex_inputs_file_prefix_c);
            if (amrex_inputs_file_prefix.back() != '/') {
                amrex_inputs_file_prefix += "/";
            }
            filename = amrex_inputs_file_prefix + filename;
        }

#ifdef AMREX_USE_MPI
        if (ParallelDescriptor::Communicator() == MPI_COMM_NULL)
        {
            throw std::runtime_error("read_file: AMReX must be initialized");
        }
#endif

        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(filename, fileCharPtr);

        std::istringstream is(fileCharPtr.data());
        std::ostringstream os_cxx(std::ios_base::out);
        std::ostringstream os_fortran(std::ios_base::out);
        bool fortran_namelist = false;
        std::regex if_regex("^\\s*#\\s*if\\s+\\(?\\s*AMREX_SPACEDIM\\s*(>|<|==|>=|<=|!=)\\s*([1-3])\\s*\\)?\\s*$"); // NOLINT
        std::regex elif_regex("^\\s*#\\s*elif\\s+\\(?\\s*AMREX_SPACEDIM\\s*(>|<|==|>=|<=|!=)\\s*([1-3])\\s*\\)?\\s*$"); // NOLINT
        std::regex else_regex("^\\s*#\\s*else\\s*$"); // NOLINT
        std::regex endif_regex("^\\s*#\\s*endif\\s*$"); // NOLINT
        std::regex if_gpu_regex("^\\s*#\\s*ifdef\\s+AMREX_USE_GPU\\s*$"); // NOLINT
        std::regex if_not_gpu_regex("^\\s*#\\s*ifndef\\s+AMREX_USE_GPU\\s*$"); // NOLINT
        std::vector<bool> valid_region;  // Keep this block or not?
        std::vector<bool> has_true;      // Has previous if/elif ever been true?
        for (std::string line; std::getline(is, line); ) {
            std::smatch sm;
            if (std::regex_match(line, if_gpu_regex)) {
#ifdef AMREX_USE_GPU
                bool r = true;
#else
                bool r = false;
#endif
                valid_region.push_back(r);
                has_true.push_back(r);
                continue;
            } else if (std::regex_match(line, if_not_gpu_regex)) {
#ifdef AMREX_USE_GPU
                bool r = false;
#else
                bool r = true;
#endif
                valid_region.push_back(r);
                has_true.push_back(r);
                continue;
            } else if (std::regex_match(line, sm, if_regex)) {
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
addDefn (std::string& def, std::vector<std::string>& val,
         std::vector<ParmParse::QuoteType>& val_quotes, ParmParse::Table& tab)
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
    AMREX_ALWAYS_ASSERT(val.size() == val_quotes.size());
    //
    // Check if this defn is a file include directive.
    //
    if ( def == ParmParse::FileKeyword && val.size() == 1 )
    {
        // We need to provide a clean environment for included file.
        auto prev_toml_table_key = std::exchange(g_toml_table_key, std::string{});
        //
        // Read file and add to this table.
        //
        const char* fname = val.front().c_str();
        read_file(fname, tab);
        g_toml_table_key = std::move(prev_toml_table_key);
    }
    else
    {
        std::string key;
        if (g_toml_table_key.empty()) {
            key = def;
        } else {
            key.append(g_toml_table_key).append(".").append(def);
        }
        tab[key].m_vals.push_back(val);
        tab[key].m_quotes.push_back(val_quotes);
    }
    val.clear();
    val_quotes.clear();
    def = std::string();
}

void
bldTable (const char*& str, ParmParse::Table& tab)
{
    std::string              cur_value;
    std::vector<std::string> cur_list;
    std::vector<ParmParse::QuoteType> cur_quotes;
    std::vector<int>         cur_linefeeds;

    for (;;)
    {
        std::string tokvalue;
        int num_linefeeds;
        bool newline_from_comment = false;
        ParmParse::QuoteType tok_quote = ParmParse::QuoteType::None;

        PType toktype = getToken(str, tokvalue, num_linefeeds, newline_from_comment, tok_quote);

        switch (toktype)
        {
        case PType::pEOF:
        {
            if (std::accumulate(cur_linefeeds.begin(), cur_linefeeds.end(), int(0)) > 0)
            {
                std::string error_message("ParmParse: Multiple lines in ");
                error_message.append(cur_value).append(" =");
                for (auto const& x : cur_list) {
                    error_message.append(" ").append(x);
                }
                error_message.append(". Must use \\ for line continuation.");
                amrex::Abort(error_message);
            }
            addDefn(cur_value,cur_list,cur_quotes,tab);
            return;
        }
        case PType::EQ_sign:
        {
            if ( !cur_list.empty() )
            {
                //
                // Read one too far, remove last name on list.
                //
                auto tmp_str = cur_list.back();
                cur_list.pop_back();
                cur_quotes.pop_back();
                cur_linefeeds.pop_back();
                if (std::accumulate(cur_linefeeds.begin(), cur_linefeeds.end(), int(0)) > 0)
                {
                    std::string error_message("ParmParse: Multiple lines in ");
                    error_message.append(cur_value).append(" =");
                    for (auto const& x : cur_list) {
                        error_message.append(" ").append(x);
                    }
                    error_message.append(". Must use \\ for line continuation.");
                    amrex::Abort(error_message);
                }
                addDefn(cur_value,cur_list,cur_quotes,tab);
                cur_value = std::move(tmp_str);
            }
            cur_linefeeds.clear();
            break;
        }
        case PType::Defn:
        {
            if ( cur_value.empty() )
            {
                cur_value = std::move(tokvalue);
                break;
            }
            //
            // Otherwise, fall through, this may be a string.
            //
            AMREX_FALLTHROUGH;
        }
        case PType::Value:
        {
            auto table_key = is_valid_table_key(tokvalue);
            bool table_header_on_newline = (num_linefeeds > 0) || newline_from_comment;
            if (cur_value.empty() && cur_list.empty() && !table_key.empty()) {
                g_toml_table_key = std::move(table_key);
            } else if ((table_header_on_newline) && !cur_list.empty() && !table_key.empty()) {
                addDefn(cur_value,cur_list,cur_quotes,tab);
                g_toml_table_key = std::move(table_key);
            } else {
                cur_list.push_back(std::move(tokvalue));
                cur_quotes.push_back(tok_quote);
                cur_linefeeds.push_back(num_linefeeds);
            }
            break;
        }
        } // switch (toktype)
    }
}

template <typename T>
bool pp_parser (const ParmParse::Table& table, const std::string& parser_prefix,
                const std::string& name, const std::string& val, T& ref,
                bool use_querywithparser);

template <typename T>
void pp_entry_set_last_val (ParmParse::PP_entry const& entry, int ival, T ref, bool parsed)
{
#ifdef AMREX_USE_OMP
#pragma omp single nowait
#endif
    {
        if (ival >= int(entry.m_last_vals.size())) {
            entry.m_last_vals.resize(ival+1);
        }
        entry.m_last_vals[ival] = ref;
        if (parsed) { entry.m_parsed = true; }
    }
}

template <typename T>
void read_array_1d (std::vector<T>& ref, std::string const& str)
{
    ref.clear();
    std::istringstream is(str);
    auto throw_parse_error = [&str]() {
        throw std::runtime_error("ParmParse: failed to parse array element in " + str);
    };
    T v{};
    is.ignore(100000, '[');
    if (!(is >> v)) {
        throw_parse_error();
    }
    ref.push_back(v);
    while (true) {
        is >> std::ws;
        auto nc = is.peek();
        if (nc == ',') {
            is.ignore(1, ',');
            is >> std::ws;
            nc = is.peek();
            if (nc == ']') { return; }
            if (!(is >> v)) {
                throw_parse_error();
            }
            ref.push_back(v);
            continue;
        } else {
            break;
        }
    }
}

bool is_escaped_quote (std::string const& str, std::string::size_type pos)
{
    // An odd number of backslashes means `"` is escaped.
    std::string::size_type bs = 0;
    for (auto i = pos; i > 0 && str[i-1] == '\\'; --i) { ++bs; }
    return (bs % 2) == 1;
}

std::size_t find_next_unquoted (std::string const& str, std::size_t start, char target)
{
    bool in_string = false;
    for (std::size_t i = start; i < str.size(); ++i) {
        char c = str[i];
        if (c == '"' && !is_escaped_quote(str, i)) {
            in_string = !in_string;
            continue;
        }
        if (!in_string && c == target) {
            return i;
        }
    }
    return std::string::npos;
}

void read_array_1d (std::vector<std::string>& ref, std::string const& str)
{
    ref.clear();
    std::string::size_type pos = str.find('[');
    if (pos == std::string::npos) { return; }
    while (true) {
        pos = str.find('"', pos+1);
        if (pos != std::string::npos) {
            auto open_pos = pos;
            while (true) {
                pos = str.find('"', pos+1);
                if (pos != std::string::npos) {
                    if (!is_escaped_quote(str, pos)) {
                        ref.push_back(str.substr(open_pos+1, pos-(open_pos+1)));
                        break;
                    }
                } else {
                    amrex::ErrorStream() << "ParmParse: unmatched quotes in string array\n";
                    amrex::Abort();
                    return;
                }
            }
        } else {
            break;
        }
    }
}

bool balanced_brackets (std::string const& str, std::string::size_type ib,
                        std::string::size_type ie)
{
    if (ie-ib >= 3) {
        if (str[ib] == '[' && str[ib+1] == '"' &&
            str[ie] == ']' && str[ie-1] == '"') {
            return !is_escaped_quote(str, ie-1);
        } else {
            return false;
        }
    } else {
        return false;
    }
}

template <typename T>
void read_array_2d (std::vector<std::vector<T>>& ref, std::string const& str)
{
    ref.clear();
    std::string::size_type pos = find_next_unquoted(str, 0, '[');
    if (pos == std::string::npos) { return; }
    for (int row_index = 0; row_index < 1000000; ++row_index) { // NOLINT
        pos = find_next_unquoted(str, pos+1, '[');
        if (pos != std::string::npos) {
            auto open_pos = pos;
            while (true) {
                pos = find_next_unquoted(str, pos+1, ']');
                if (pos != std::string::npos) {
                    if constexpr (std::is_same_v<T,std::string>) {
                        if (!balanced_brackets(str, open_pos, pos)) {
                            continue; // continue the searching for ']'
                        }
                    }
                    ref.resize(row_index+1);
                    try {
                        read_array_1d(ref[row_index], str.substr(open_pos, pos-open_pos+1));
                    } catch (...) {
                        throw;
                    }
                    break;
                } else {
                    throw std::runtime_error("ParmParse: unmatched [] in nested arrays\n");
                }
            }
        } else {
            break;
        }
    }
}

bool is_toml_1d_array (std::string const& token);

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
    // Get specified occurrence of name in table.
    //
    auto const* def = ppindex(table, occurrence, name);
    if ( def == nullptr )
    {
        return false;
    }

    auto const& entry = table.at(name);

#ifdef AMREX_USE_OMP
#pragma omp single nowait
#endif
    {
        using T_ptr = std::decay_t<T>*;
        entry.m_typehint = static_cast<T_ptr>(nullptr);
    }

    //
    // Handle TOML array: stored as single token "$$ARR[v0,v1,...]"
    //
    if (!(def->empty()) && is_toml_1d_array((*def)[0])) {
        std::vector<T> toml_vals;
        read_array_1d(toml_vals, (*def)[0]);
        if (ival >= static_cast<int>(toml_vals.size())) {
            amrex::ErrorStream() << "ParmParse::queryval no value number "
                                 << ival << " for ";
            if ( occurrence ==  ParmParse::LAST ) {
                amrex::ErrorStream() << "last occurrence of ";
            } else {
                amrex::ErrorStream() << " occurrence " << occurrence << " of ";
            }
            amrex::ErrorStream() << name << '\n' << pp_to_string(name,*def) << '\n';
            amrex::Abort();
        }
        ref = toml_vals[ival];
        return true;
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

    constexpr bool is_integral_floating = (std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>);

    if (is(valname, ref)) {
        if constexpr (is_integral_floating) {
            pp_entry_set_last_val(entry, ival, ref, false);
        }
    } else {
        if constexpr (is_integral_floating) {
            if (pp_parser(table, parser_prefix, name, valname, ref, false)) {
                pp_entry_set_last_val(entry, ival, ref, true);
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
                             << " not found in database"
                             << '\n';
        ParmParse::dumpTable(amrex::ErrorStream());
        amrex::Abort();
    }
}

// Checks if token matches $$ARR[...]
bool is_toml_array (std::string const& token)
{
    return token.size() >= 7 && token.compare(0,6,"$$ARR[") == 0 &&
        token.back() == ']';
}

// Checks if token matches $$ARR[[...]]
bool is_toml_2d_array (std::string const& token)
{
    auto sz = token.size();
    return sz >= 9 && token.compare(0,7,"$$ARR[[") == 0 &&
        token.compare(sz-2,2,"]]") == 0;
}

// Checks if token matches $$ARR[...] but not $$ARR[[...]]
bool is_toml_1d_array (std::string const& token)
{
    return is_toml_array(token) && !is_toml_2d_array(token);
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

    bool const toml_array = !(def->empty()) && is_toml_array((*def)[0]);
    std::vector<T> toml_vals;
    if (toml_array) {
        read_array_1d(toml_vals, (*def)[0]);
    }

    auto const& entry = table.at(name);

#ifdef AMREX_USE_OMP
#pragma omp single nowait
#endif
    {
        using T_ptr = std::decay_t<T>*;
        entry.m_typehint = static_cast<T_ptr>(nullptr);
    }

    //
    // Does it have sufficient number of values and are they all
    // the same type?
    //
    int available = toml_array ? static_cast<int>(toml_vals.size())
                               : static_cast<int>(def->size());
    if ( num_val == ParmParse::ALL )
    {
        num_val = available;
    }

    if ( num_val == 0 ) { return true; }

    constexpr bool is_integral_floating = (std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>);

    int stop_ix = start_ix + num_val - 1;
    if ( static_cast<int>(ref.size()) <= stop_ix )
    {
        ref.resize(stop_ix + 1);
    }
    if ( stop_ix >= available )
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
    if (toml_array) {
        for (int n = start_ix; n <= stop_ix; ++n) {
            ref[n] = toml_vals[n];
            if constexpr (is_integral_floating) {
                pp_entry_set_last_val(entry, n, ref[n], true);
            }
        }
        return true;
    }

    for ( int n = start_ix; n <= stop_ix; n++ )
    {
        const std::string& valname = (*def)[n];
        if (is(valname, ref[n])) {
            if constexpr (is_integral_floating) {
                pp_entry_set_last_val(entry, n, ref[n], false);
            }
        } else {
            if constexpr (is_integral_floating) {
                if (pp_parser(table, parser_prefix, name, valname, ref[n], false)) {
                    pp_entry_set_last_val(entry, n, ref[n], true);
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
                             << " not found in database"
                             << '\n';
        ParmParse::dumpTable(amrex::ErrorStream());
        amrex::Abort();
    }
}

template <class T>
std::string to_toml_value (const T& ref)
{
    using TT = std::remove_reference_t<T>;
    std::stringstream ss;
    if constexpr (std::is_floating_point_v<TT>) {
        ss << std::setprecision(std::numeric_limits<TT>::max_digits10);
    } else if constexpr (std::is_same_v<TT,bool>) {
        ss << std::boolalpha;
    }
    ss << ref;
    std::string s = ss.str();
    if constexpr (std::is_floating_point_v<TT>) {
        const std::regex digits_only(R"([+-]?\d+)");
        if (std::regex_match(s, digits_only)) {
            s += ".0";
        }
    }
    return s;
}

template <class T>
void
saddval (const std::string& name, const T& ref)
{
    std::string s = to_toml_value(ref);
    auto& entry = g_table[name];
    entry.m_vals.emplace_back(1, std::move(s));
    auto qt = std::is_same_v<std::remove_reference_t<T>,std::string>
        ? ParmParse::QuoteType::Double : ParmParse::QuoteType::None;
    entry.m_quotes.emplace_back(1, qt);
    ++entry.m_count;
    using T_ptr = std::decay_t<T>*;
    entry.m_typehint = static_cast<T_ptr>(nullptr);
}

template <class T>
void
saddarr (const std::string& name, const std::vector<T>& ref)
{
    std::vector<std::string> arr;
    arr.reserve(ref.size());
    for (auto const& item : ref) {
        arr.push_back(to_toml_value(item));
    }

    auto& entry = g_table[name];
    auto arr_size = arr.size();
    entry.m_vals.emplace_back(std::move(arr));
    auto qt = std::is_same_v<std::remove_reference_t<T>,std::string>
        ? ParmParse::QuoteType::Double : ParmParse::QuoteType::None;
    entry.m_quotes.emplace_back(arr_size, qt);
    ++entry.m_count;
    using T_ptr = std::decay_t<T>*;
    entry.m_typehint = static_cast<T_ptr>(nullptr);
}

// Initialize ParmParse.
void
ppinit (int argc, char** argv, const char* parfile, ParmParse::Table& table)
{
    g_toml_table_key.clear();

    // Check environment first
    if (char const* env = std::getenv("AMREX_DEFAULT_INIT")) {
        std::string env_s = std::string(env) + '\n';
        char const* s = env_s.c_str();
        bldTable(s, table);
    }

    g_toml_table_key.clear();

    if ( parfile != nullptr )
    {
        read_file(parfile, table);
    }

    g_toml_table_key.clear();

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
            auto& src_quotes = arg_entry.m_quotes;
            auto& dst_quotes = table[name].m_quotes;
            std::move(std::begin(src_quotes), std::end(src_quotes), std::back_inserter(dst_quotes));
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

template <typename T, typename PARSER_t = std::conditional_t<std::is_integral_v<T>
                                                             && !std::is_same_v<bool,T>,
                                                             IParser, Parser>>
PARSER_t
pp_make_parser (std::string const& func, Vector<std::string> const& vars,
                ParmParse::Table const& table, std::string const& parser_prefix,
                bool use_querywithparser)
{
    using value_t =  std::conditional_t<std::is_integral_v<T> && !std::is_same_v<bool,T>,
                                        long long, double>;

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

    bool recursive = false;
    auto& recursive_symbols = g_parser_recursive_symbols[OpenMP::get_thread_num()];

    for (auto const& s : symbols) {
        value_t v = 0;
        bool r = false;
        for (auto const& pf : prefixes) {
            std::string pfs = pf + s;
            if (auto found = recursive_symbols.find(pfs); found != recursive_symbols.end()) {
                recursive = true;
                continue;
            }
            if (use_querywithparser) {
                r = squeryWithParser(table, parser_prefix, pfs, v);
            } else {
                r = squeryval(table, parser_prefix, pfs, v,
                              ParmParse::FIRST, ParmParse::LAST);
            }
            if (r) { break; }
        }
        if (r == false) {
            std::string msg("ParmParse: failed to parse "+func);
            if (recursive) {
                msg.append(" due to recursive symbol ").append(s);
            } else {
                msg.append(" due to unknown symbol ").append(s);
            }
            amrex::Error(msg);
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

std::string const&
ParmParse::getPrefix () const
{
    return m_prefix;
}

std::string
ParmParse::prefixedName (std::string_view str) const
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
    // this is required because we will BCast the file content in sub-function calls
    if (ParallelDescriptor::Communicator() == MPI_COMM_NULL)
    {
        throw std::runtime_error("ParmParse::addfile: AMReX must be initialized");
    }
#endif

    // check the file exists and give a user-friendly error
    if (ParallelDescriptor::IOProcessor())
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            FileExists(filename),
            "ParmParse::addfile: file does not exist: " + filename
        );
    }

    g_toml_table_key.clear();

    // add the file
    auto file = FileKeyword;
    std::vector<std::string> val(1, filename);
    std::vector<ParmParse::QuoteType> val_quotes(1, ParmParse::QuoteType::None);
    addDefn(file, val, val_quotes, g_table);

    g_toml_table_key.clear();
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
    g_parser_recursive_symbols.resize(1);

    pp_detail::verbose = -1;
    initialized = false;
}

void
ParmParse::SetParserPrefix (std::string a_prefix)
{
    ParmParse::ParserPrefix = std::move(a_prefix);
}

// dumpTable is a diagnostic view and its output is not intended to be fed back to ParmParse.
// Use prettyPrintTable when you need canonical, re-readable ParmParse syntax.
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
            for (std::size_t i = 0; i < entry.m_vals.size(); ++i) {
                auto const& vals = entry.m_vals[i];
                auto const* quotes = (i < entry.m_quotes.size())
                    ? &(entry.m_quotes[i]) : nullptr;
                os << pp_to_pretty_string(name, vals, quotes, nullptr) << '\n';
            }
        }
        else {
            for (auto const& vals : entry.m_vals) {
                os << pp_to_string(name, vals) << '\n';
            }
        }
    }
}

namespace {

enum class PPFlag { all, unused, used };

void pretty_print_table (std::ostream& os, PPFlag pp_flag)
{
    std::vector<std::string> sorted_names;
    sorted_names.reserve(g_table.size());
    for (auto const& [name, entry] : g_table) {
        bool to_print;
        if (pp_flag == PPFlag::used) {
            to_print = (entry.m_count > 0);
        } else if (pp_flag == PPFlag::unused) {
            to_print = (entry.m_count == 0);
        } else {
            to_print = true;
        }
        if (to_print) { sorted_names.push_back(name); }
    }
    std::sort(sorted_names.begin(), sorted_names.end());

    for (auto const& name : sorted_names) {
        auto const& entry = g_table[name];
        if (! entry.m_vals.empty()) {
            auto const idx = entry.m_vals.size() - 1;
            auto const& val = entry.m_vals[idx];
            auto const* quotes = (idx < entry.m_quotes.size())
                ? &(entry.m_quotes[idx]) : nullptr;
            os << pp_to_pretty_string(name, val, quotes, &entry) << '\n';
        }
    }
}

}

// prettyPrintTable emits valid ParmParse syntax that can be parsed again.
void
ParmParse::prettyPrintTable (std::ostream& os)
{
    pretty_print_table(os, PPFlag::all);
}

void
ParmParse::prettyPrintUnusedInputs (std::ostream& os)
{
    pretty_print_table(os, PPFlag::unused);
}

void
ParmParse::prettyPrintUsedInputs (std::ostream& os)
{
    pretty_print_table(os, PPFlag::used);
}

int
ParmParse::countval (std::string_view name,
                     int              n) const
{
    //
    // First find n'th occurrence of name in table.
    //
    auto const* def = ppindex(*m_table, n, prefixedName(name));
    if (def == nullptr) { return 0; }

    if (!(def->empty()) && is_toml_array((*def)[0])) {
        auto const& token = (*def)[0];
        std::size_t count = 0;
        if (is_toml_2d_array(token)) {
            std::string::size_type pos = find_next_unquoted(token, 0, '[');
            while (true) {
                pos = find_next_unquoted(token, pos+1, '[');
                if (pos != std::string::npos) {
                    pos = find_next_unquoted(token, pos+1, ']');
                    if (pos != std::string::npos) {
                        ++count;
                    } else {
                        throw std::runtime_error("ParmParse: unmatched [] in nested arrays\n");
                    }
                } else {
                    break;
                }
            }
        } else {
            std::string::size_type pos = find_next_unquoted(token, 0, '[');
            while (true) {
                pos = find_next_unquoted(token, pos+1, ',');
                if (pos != std::string::npos) {
                    if (pos+1 < token.size() && token[pos+1] != ']') {
                        // Note that we don't need to worry about ", ]",
                        // because getToken eats garbage.
                        ++count;
                    }
                } else {
                    break;
                }
            }
            if (token != "$$ARR[]") {
                ++count; // unless it's an empty array, increase by 1 because we counted commas, not really the number of elements
            }
        }
        return static_cast<int>(count);
    } else {
        return static_cast<int>(def->size());
    }
}

// BOOL
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   bool&            ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                bool&            ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     bool&            ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  bool&            ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view name, // NOLINT(readability-make-member-function-const)
                bool       val)
{
    saddval(prefixedName(name),val);
}

// INT
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   int&             ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                int&             ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     int&             ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  int&             ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view name, // NOLINT(readability-make-member-function-const)
                int        val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (std::string_view  name,
                      int               k,
                      std::vector<int>& ref,
                      int               start_ix,
                      int               num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (std::string_view  name,
                   std::vector<int>& ref,
                   int               start_ix,
                   int               num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (std::string_view  name,
                        int               k,
                        std::vector<int>& ref,
                        int               start_ix,
                        int               num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (std::string_view  name,
                     std::vector<int>& ref,
                     int               start_ix,
                     int               num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (std::string_view        name, // NOLINT(readability-make-member-function-const)
                   const std::vector<int>& ref)
{
    saddarr(prefixedName(name),ref);
}


// LONG
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   long&            ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                long&            ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     long&            ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  long&            ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view name, // NOLINT(readability-make-member-function-const)
                long       val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (std::string_view   name,
                      int                k,
                      std::vector<long>& ref,
                      int                start_ix,
                      int                num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (std::string_view   name,
                   std::vector<long>& ref,
                   int                start_ix,
                   int                num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (std::string_view   name,
                        int                k,
                        std::vector<long>& ref,
                        int                start_ix,
                        int                num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (std::string_view   name,
                     std::vector<long>& ref,
                     int                start_ix,
                     int                num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (std::string_view name, const std::vector<long>& ref) // NOLINT(readability-make-member-function-const)
{
    saddarr(prefixedName(name),ref);
}

// long long
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   long long&       ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                long long&       ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     long long&       ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  long long&       ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view name, // NOLINT(readability-make-member-function-const)
                long long  val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (std::string_view        name,
                      int                     k,
                      std::vector<long long>& ref,
                      int                     start_ix,
                      int                     num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (std::string_view        name,
                   std::vector<long long>& ref,
                   int                     start_ix,
                   int                     num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (std::string_view        name,
                        int                     k,
                        std::vector<long long>& ref,
                        int                     start_ix,
                        int                     num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

int
ParmParse::queryarr (std::string_view        name,
                     std::vector<long long>& ref,
                     int                     start_ix,
                     int                     num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (std::string_view              name, // NOLINT(readability-make-member-function-const)
                   const std::vector<long long>& ref)
{
    saddarr(prefixedName(name),ref);
}

// FLOAT
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   float&           ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                float&           ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     float&           ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  float&           ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view name, // NOLINT(readability-make-member-function-const)
                float      val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (std::string_view    name,
                      int                 k,
                      std::vector<float>& ref,
                      int                 start_ix,
                      int                 num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (std::string_view    name,
                   std::vector<float>& ref,
                   int                 start_ix,
                   int                 num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (std::string_view    name,
                        int                 k,
                        std::vector<float>& ref,
                        int                 start_ix,
                        int                 num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (std::string_view    name,
                     std::vector<float>& ref,
                     int                 start_ix,
                     int                 num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (std::string_view          name, // NOLINT(readability-make-member-function-const)
                   const std::vector<float>& ref)
{
    saddarr(prefixedName(name),ref);
}



// DOUBLE
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   double&          ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                double&          ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     double&          ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  double&          ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view name, // NOLINT(readability-make-member-function-const)
                double     val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (std::string_view     name,
                      int                  k,
                      std::vector<double>& ref,
                      int                  start_ix,
                      int                  num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (std::string_view     name,
                   std::vector<double>& ref,
                   int                  start_ix,
                   int                  num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (std::string_view     name,
                        int                  k,
                        std::vector<double>& ref,
                        int                  start_ix,
                        int                  num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (std::string_view     name,
                     std::vector<double>& ref,
                     int                  start_ix,
                     int                  num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (std::string_view           name, // NOLINT(readability-make-member-function-const)
                   const std::vector<double>& ref)
{
    saddarr(prefixedName(name),ref);
}



// STRING
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   std::string&     ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                std::string&     ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     std::string&     ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  std::string&     ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view   name, // NOLINT(readability-make-member-function-const)
                const std::string& val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (std::string_view          name,
                      int                       k,
                      std::vector<std::string>& ref,
                      int                       start_ix,
                      int                       num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (std::string_view          name,
                   std::vector<std::string>& ref,
                   int                       start_ix,
                   int                       num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (std::string_view          name,
                        int                       k,
                        std::vector<std::string>& ref,
                        int                       start_ix,
                        int                       num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (std::string_view          name,
                     std::vector<std::string>& ref,
                     int                       start_ix,
                     int                       num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (std::string_view                name, // NOLINT(readability-make-member-function-const)
                   const std::vector<std::string>& ref)
{
    saddarr(prefixedName(name),ref);
}



// INTVECT
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   IntVect&         ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                IntVect&         ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     IntVect&         ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  IntVect&         ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view name, // NOLINT(readability-make-member-function-const)
                const IntVect&   val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (std::string_view      name,
                      int                   k,
                      std::vector<IntVect>& ref,
                      int                   start_ix,
                      int                   num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (std::string_view      name,
                   std::vector<IntVect>& ref,
                   int                   start_ix,
                   int                   num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (std::string_view      name,
                        int                   k,
                        std::vector<IntVect>& ref,
                        int                   start_ix,
                        int                   num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (std::string_view      name,
                     std::vector<IntVect>& ref,
                     int                   start_ix,
                     int                   num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (std::string_view            name, // NOLINT(readability-make-member-function-const)
                   const std::vector<IntVect>& ref)
{
    saddarr(prefixedName(name),ref);
}

// BOX
void
ParmParse::getkth (std::string_view name,
                   int              k,
                   Box&             ref,
                   int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

void
ParmParse::get (std::string_view name,
                Box&             ref,
                int              ival) const
{
    sgetval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

int
ParmParse::querykth (std::string_view name,
                     int              k,
                     Box&             ref,
                     int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival,k);
}

int
ParmParse::query (std::string_view name,
                  Box&             ref,
                  int              ival) const
{
    return squeryval(*m_table,m_parser_prefix, prefixedName(name),ref,ival, LAST);
}

void
ParmParse::add (std::string_view name, // NOLINT(readability-make-member-function-const)
                const Box&       val)
{
    saddval(prefixedName(name),val);
}

void
ParmParse::getktharr (std::string_view  name,
                      int               k,
                      std::vector<Box>& ref,
                      int               start_ix,
                      int               num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val,k);
}

void
ParmParse::getarr (std::string_view  name,
                   std::vector<Box>& ref,
                   int               start_ix,
                   int               num_val) const
{
    sgetarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (std::string_view  name,
                        int               k,
                        std::vector<Box>& ref,
                        int               start_ix,
                        int               num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix, num_val,k);
}

int
ParmParse::queryarr (std::string_view  name,
                     std::vector<Box>& ref,
                     int               start_ix,
                     int               num_val) const
{
    return squeryarr(*m_table,m_parser_prefix, prefixedName(name),ref,start_ix,num_val, LAST);
}

void
ParmParse::addarr (std::string_view        name, // NOLINT(readability-make-member-function-const)
                   const std::vector<Box>& ref)
{
    saddarr(prefixedName(name),ref);
}


int
ParmParse::queryarr (std::string_view name,
                     IntVect&         ref) const
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
ParmParse::getarr (std::string_view name, IntVect& ref) const
{
    std::vector<int> v;
    this->getarr(name, v);
    AMREX_ALWAYS_ASSERT(v.size() == AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) { ref[i] = v[i]; }
}

int
ParmParse::queryarr (std::string_view name, RealVect& ref) const
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
ParmParse::getarr (std::string_view name, RealVect& ref) const
{
    std::vector<Real> v;
    this->getarr(name, v);
    AMREX_ALWAYS_ASSERT(v.size() == AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) { ref[i] = v[i]; }
}

void
ParmParse::getline (std::string_view name, std::string& ref) const
{
    std::vector<std::string> tmp;
    getarr(name, tmp);
    ref = amrex::join(tmp, ' ');
}

int
ParmParse::queryline (std::string_view name, std::string& ref) const
{
    std::vector<std::string> tmp;
    auto r = queryarr(name, tmp);
    if (r) {
        ref = amrex::join(tmp, ' ');
    }
    return r;
}

//
// Return number of occurrences of parameter name.
//

int
ParmParse::countname (std::string_view name) const
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
ParmParse::contains (std::string_view name) const
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
ParmParse::remove (std::string_view name)
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

    constexpr bool is_integral_floating = (std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>);

    auto const& entry = table.at(name);

    if (pp_parser(table, parser_prefix, name, combined_string, ref, true))
    {
#ifdef AMREX_USE_OMP
#pragma omp single nowait
#endif
        {
            using T_ptr = std::decay_t<T>*;
            entry.m_typehint = static_cast<T_ptr>(nullptr);
        }

        if constexpr (is_integral_floating) {
            pp_entry_set_last_val(entry, 0, ref, true);
        }

        return true;
    } else {
        return false;
    }
}

template <class T>
bool squeryarrWithParser (const ParmParse::Table& table,
                          const std::string&      parser_prefix,
                          const std::string&      name,
                          int                     nvals,
                          T*                      ptr)
{
    std::vector<std::string> vals;
    bool exist = squeryarr(table, parser_prefix, name, vals,
                           ParmParse::FIRST, ParmParse::ALL, ParmParse::LAST);
    if (!exist) { return false; }

    constexpr bool is_integral_floating = (std::is_same_v<T,bool> ||
                                           std::is_same_v<T,int> ||
                                           std::is_same_v<T,long> ||
                                           std::is_same_v<T,long long> ||
                                           std::is_same_v<T,float> ||
                                           std::is_same_v<T,double>);

    auto const& entry = table.at(name);

    AMREX_ALWAYS_ASSERT(int(vals.size()) == nvals);
    for (int ival = 0; ival < nvals; ++ival) {
        bool r = pp_parser(table, parser_prefix, name, vals[ival], ptr[ival], true);
        if (r) {
            if constexpr (is_integral_floating) {
                pp_entry_set_last_val(entry, ival, ptr[ival], true);
            }
        } else {
            return false;
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp single nowait
#endif
    {
        using T_ptr = std::decay_t<T>*;
        entry.m_typehint = static_cast<T_ptr>(nullptr);
    }

    return true;
}
}

int
ParmParse::queryWithParser (std::string_view name, bool& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (std::string_view name, int& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (std::string_view name, long& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (std::string_view name, long long& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (std::string_view name, float& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryWithParser (std::string_view name, double& ref) const
{
    return squeryWithParser(*m_table,m_parser_prefix,prefixedName(name),ref);
}

int
ParmParse::queryarrWithParser (std::string_view name, int nvals, bool* ptr) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ptr);
}

int
ParmParse::queryarrWithParser (std::string_view name, int nvals, int* ptr) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ptr);
}

int
ParmParse::queryarrWithParser (std::string_view name, int nvals, long* ptr) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ptr);
}

int
ParmParse::queryarrWithParser (std::string_view name, int nvals, long long* ptr) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ptr);
}

int
ParmParse::queryarrWithParser (std::string_view name, int nvals, float* ptr) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ptr);
}

int
ParmParse::queryarrWithParser (std::string_view name, int nvals, double* ptr) const
{
    return squeryarrWithParser(*m_table,m_parser_prefix,prefixedName(name),nvals,ptr);
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

namespace
{
template <typename T>
std::vector<T> read_table_row (std::istream& is)
{
    std::vector<T> r;
    is >> std::ws;
    char c;
    T v{};
    is >> c;
    if (c == '{') {
        is >> v;
        r.push_back(v);
        while (true) {
            is >> std::ws;
            auto nc = is.peek();
            if (nc == ',') {
                is.ignore(10000, ',');
                is >> v;
                r.push_back(v);
                continue;
            } else {
                break;
            }
        }
        is.ignore(100000,  '}');
    } else {
        amrex::Error("ParmParse::querytable: read_table_row expected \'{\'");
    }
    if (is.fail()) {
        amrex::Error("ParmParse::querytable read_table_row failed to read table");
    }
    return r;
}

template <typename T>
void read_table_2d (std::vector<std::vector<T>>& ref, std::string const& str)
{
    std::istringstream is(str);
    is >> std::ws;
    char c;
    is >> c;
    if (c == '{') {
        for (int row_index = 0; row_index < 1000000; ++row_index) {
            if (auto row = read_table_row<T>(is); ! row.empty()) {
                if (row_index == 0) { ref.clear(); }
                ref.emplace_back(std::move(row));
            } else {
                break;
            }
            is >> std::ws;
            auto nc = is.peek();
            if (nc == ',') {
                is >> c; // skip optional ','
                is >> std::ws;
                nc = is.peek();
            }
            if (nc == '}') { break; }
        }
        is.ignore(100000,  '}');
    } else {
        amrex::Error("ParmParse::querytable: read_table_2d expected \'{\'");
    }
    if (is.fail()) {
        amrex::Error("ParmParse::querytable read_table_2d failed to read table");
    }
}
}

int ParmParse::querytable (std::string_view name, std::vector<std::vector<double>>& ref) const
{
    std::string table_s;
    int r = query(name, table_s);
    if (r) {
        read_table_2d(ref, table_s);
    }
    return r;
}

int ParmParse::querytable (std::string_view name, std::vector<std::vector<float>>& ref) const
{
    std::string table_s;
    int r = query(name, table_s);
    if (r) {
        read_table_2d(ref, table_s);
    }
    return r;
}

int ParmParse::querytable (std::string_view name, std::vector<std::vector<int>>& ref) const
{
    std::string table_s;
    int r = query(name, table_s);
    if (r) {
        read_table_2d(ref, table_s);
    }
    return r;
}

int ParmParse::queryarr (std::string_view name, std::vector<std::vector<double>>& ref) const
{
    std::string arr_s;
    int r = query(name, arr_s);
    if (r) {
        read_array_2d(ref, arr_s);
    }
    return r;
}

int ParmParse::queryarr (std::string_view name, std::vector<std::vector<float>>& ref) const
{
    std::string arr_s;
    int r = query(name, arr_s);
    if (r) {
        read_array_2d(ref, arr_s);
    }
    return r;
}

int ParmParse::queryarr (std::string_view name, std::vector<std::vector<int>>& ref) const
{
    std::string arr_s;
    int r = query(name, arr_s);
    if (r) {
        read_array_2d(ref, arr_s);
    }
    return r;
}

int ParmParse::queryarr (std::string_view name, std::vector<std::vector<std::string>>& ref) const
{
    std::string arr_s;
    int r = query(name, arr_s);
    if (r) {
        read_array_2d(ref, arr_s);
    }
    return r;
}

}
