//
// $Id: ParmParse.cpp,v 1.41 2002-03-22 23:29:24 car Exp $
//
#include <winstd.H>

#include <typeinfo>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <vector>
#include <list>

#include <BoxLib.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Box.H>
#include <IntVect.H>

//
// Used by constructor to build table.
//
ParmParse::PP_entry::PP_entry (const std::string& name,
                    const std::list<std::string>& vals)
    :
    m_name(name),
    m_table(0),
    m_queried(false)
{
    for ( std::list<std::string>::const_iterator li = vals.begin(); li != vals.end(); ++li )
    {
	m_vals.push_back(*li);
    }
}

ParmParse::PP_entry::PP_entry (const std::string& name,
                    const std::list<PP_entry>&   table)
    :
    m_name(name),
    m_table(new Table(table)),
    m_queried(false)
{
}

ParmParse::PP_entry::PP_entry (const PP_entry& pe)
    : m_name(pe.m_name),
      m_vals(pe.m_vals),
      m_table(0),
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
    m_table = 0;
    m_queried = pe.m_queried;
    if ( pe.m_table )
    {
	m_table = new Table(*pe.m_table);
    }
    return *this;
}

std::ostream&
operator<< (std::ostream& os, const ParmParse::PP_entry& pp)
{
    os << pp.m_name << "(nvals = " << pp.m_vals.size() << ") " << " :: [";
    int n = pp.m_vals.size();
    for ( int i = 0; i < n; i++ )
    {
	os << pp.m_vals[i];
	if ( i < n-1 ) os << ", ";
    }
    os << "]";

    if ( !os )
    {
        BoxLib::Error("write on ostream failed");
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
is (const std::string& str, T& val)
{
    std::istringstream s(str);
    s >> val;
    if ( !s ) return false;
    return true;
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
    if ( is(str, int_val) )
    {
        val = int_val != 0;
        return true;
    }
    double dbl_val;
    if ( is(str, dbl_val) )
    {
        val = dbl_val != 0;
        return true;
    }
    return false;
}

ParmParse::Table g_table;
typedef std::list<ParmParse::PP_entry>::iterator list_iterator;
typedef std::list<ParmParse::PP_entry>::const_iterator const_list_iterator;

template <class T> const char* tok_name(const T&) { return typeid(T).name(); }
template <class T> const char* tok_name(std::vector<T>&) { return tok_name(T());}

#if 0
template <> const char* tok_name(const bool&)        { return "bool";        }
template <> const char* tok_name(const int&)         { return "int";         }
template <> const char* tok_name(const float&)       { return "float";       }
template <> const char* tok_name(const double&)      { return "double";      }
template <> const char* tok_name(const std::string&) { return "std::string"; }
template <> const char* tok_name(const Box&)         { return "Box";         }
template <> const char* tok_name(const IntVect&)     { return "IntVect";     }
#endif

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
	if ( *str == 0 ) break;
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
	  char*       ostr)
{
#define ERROR_MESS 							\
   ostr[k++] = '\0';							\
   std::cerr << "ParmParse::getToken(): invalid string = " << ostr << '\n'; \
   std::cerr << "STATE = " << state_name[state]				\
             << ", next char = " << ch << '\n';				\
   std::cerr << ", rest of input = \n" << str << '\n';		\
   BoxLib::Abort()
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
   int k = 0;			// index of output string
   int pcnt = 0;		// Tracks nested parens
   while (1)
   {
       char ch = *str;
       if ( ch == 0 )
       {
	   BoxLib::Error("ParmParse::getToken: EOF while parsing");
       }
       switch (state)
       {
       case START:
           if ( ch == '=' )
           {
               ostr[k++] = ch; str++;
               ostr[k++] = 0;
               return pEQ_sign;
           }
           else if ( ch == '"' )
           {
               str++;
               state = QUOTED_STRING;
           }
	   else if ( ch == '(' )
	   {
	       ostr[k++] = ch; str++; pcnt = 1;
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
               ostr[k++] = ch; str++;
               state = IDENTIFIER;
           }
           else
           {
               ostr[k++] = ch; str++;
               state = STRING;
           }
           break;
       case IDENTIFIER:
           if ( isalnum(ch) || ch == '_' || ch == '.' || ch == '[' || ch == ']' )
           {
               ostr[k++] = ch; str++;
           }
           else if ( isspace(ch) || ch == '=' )
           {
               ostr[k++] = 0;
               return pDefn;
           }
           else
           {
               ostr[k++] = ch; str++;
               state = STRING;
           }
           break;
       case LIST:
	   if ( ch == '(' )
	   {
	       ostr[k++] = ch; str++; pcnt++;
	   }
	   else if ( ch == ')' )
	   {
	       ostr[k++] = ch; str++; pcnt--;
	       if ( pcnt == 0 )
	       {
		   ostr[k++] = 0;
		   return pValue;
	       }
	   }
	   else
	   {
	       ostr[k++] = ch; str++;
	   }
	   break;
       case STRING:
           if ( isspace(ch) || ch == '=' )
           {
               ostr[k++] = 0;
               return pValue;
           }
           else
           {
               ostr[k++] = ch; str++;
           }
           break;
       case QUOTED_STRING:
           if ( ch == '"' )
           {
               str++;
               ostr[k++] = 0;
               return pValue;
           }
           else
           {
               ostr[k++] = ch; str++;
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

static
bool
ppfound (const std::string& keyword,
	 const ParmParse::PP_entry& pe,
	 bool recordQ)
{
    return (recordQ == (pe.m_table!=0)) && (keyword == pe.m_name);
}

//
// Return the index of the n'th occurence of a parameter name,
// except if n==-1, return the index of the last occurence.
// Return 0 if the specified occurence does not exist.
//

const ParmParse::PP_entry*
ppindex (const ParmParse::Table& table,
	 int         n,
	 const std::string& name,
	 bool recordQ)
{
    const ParmParse::PP_entry* fnd = 0;

    if ( n == ParmParse::LAST )
    {
        //
        // Search from back of list.
        //
        for (std::list<ParmParse::PP_entry>::const_reverse_iterator li = table.rbegin(); li != table.rend(); ++li)
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
        for ( const_list_iterator li=table.begin(); li != table.end(); ++li )
        {
            if ( ppfound(name, *li, recordQ) )
            {
                fnd = &*li;
                if ( --n < 0 )
		{
                    break;
		}
            }
        }
        if ( n >= 0)
	{
            fnd = 0;
	}
    }

    if ( fnd )
    {
        //
        // Found an entry; mark all occurences of name as used.
        //
        for ( const_list_iterator li=table.begin(); li != table.end(); ++li )
	{
            if ( ppfound(name, *li, recordQ) )
	    {
                li->m_queried = true;
	    }
	}
    }
    return fnd;
}

void
bldTable (const char*& str, std::list<ParmParse::PP_entry>& tab);

static void
read_file (const char*           fname,
	   std::list<ParmParse::PP_entry>& tab)
{
    //
    // Space for input file if it exists.
    // Note: on CRAY, access requires (char*) not (const char*).
    //
    if ( fname != 0 && fname[0] != 0 )
    {
        FILE* pffd = fopen(fname, "rb");
        if ( pffd == 0 )
        {
            std::cerr << "ParmParse::read_file(): couldn't open \""
                      << fname
                      << "\"";
            BoxLib::Abort();
        }
        //
        // Get the length.
        //
        fseek(pffd, 0, 2);
        int pflen = (int)ftell(pffd);
        rewind(pffd);
        char* str = new char[pflen+1];
        memset(str,0,pflen+1);
        int nread = fread(str, 1, pflen, pffd);
        if ( !(nread == pflen) )
        {
            std::cerr << "ParmParse::read_file(): fread() only "
                      << nread
                      << " bytes out of "
                      << pflen
                      << " from "
                      << fname;
            BoxLib::Abort();
        }
        fclose(pffd);
	const char* b = str;
        bldTable(b, tab);
        delete [] str;
    }
}

void
addDefn (std::string&         def,
	 std::list<std::string>&   val,
	 std::list<ParmParse::PP_entry>& tab)
{
    static const std::string FileKeyword("FILE");
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
        std::cerr << "ParmParse::addDefn(): no values for definition " << def << "\n";
        BoxLib::Abort();
    }
    //
    // Check if this defn is a file include directive.
    //
    if ( def == FileKeyword && val.size() == 1 )
    {
        //
        // Read file and add to this table.
        //
        const char* fname = val.front().c_str();
        read_file(fname, tab);
    }
    else
    {
        tab.push_back(ParmParse::PP_entry(def,val));
    }
    val.clear();
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
        std::cerr << "ParmParse::addTable(): no values for Table " << def << "\n";
	BoxLib::Abort();
    }
    tab.push_back(ParmParse::PP_entry(def, val));
    val.clear();
    def = std::string();
}

void
bldTable (const char*&           str,
	  std::list<ParmParse::PP_entry>& tab)
{
    std::string       cur_name;
    std::list<std::string> cur_list;
    ParmParse::Table  cur_table;
    std::string       tmp_str;

    PType    token;
    const int SCRATCH_STR_LEN  = 200;
    char      tokname[SCRATCH_STR_LEN];

    for (;;)
    {
	token = getToken(str, tokname);

	switch (token)
	{
	case pCloseBracket:
	    if ( !cur_name.empty() && cur_list.empty() )
	    {
		BoxLib::Abort("ParmParse::bldTable() defn with no list");
	    }
	case pEOF:
	    addDefn(cur_name,cur_list,tab);
	    return;
	case pOpenBracket:
	    if ( cur_name.empty() )
	    {
		BoxLib::Abort("ParmParse::bldTabe() '{' with no blocknamne");
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
		BoxLib::Abort("ParmParse::bldTable() EQ with no current defn");
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
	case pValue:
	    if ( cur_name.empty() )
	    {
		tokname[SCRATCH_STR_LEN-1] = 0;
		std::string msg("ParmParse::bldTable(): value with no defn: ");
		msg += tokname;
		BoxLib::Abort(msg.c_str());
	    }
	    cur_list.push_back(tokname);
	    break;
	}
    }
}

template <class T>
static
bool
squeryval (const ParmParse::Table& table,
	   const std::string& name,
	   T&           ptr,
	   int          ival,
	   int          occurence)
{
    //
    // Get last occurrance of name in table.
    //
    const ParmParse::PP_entry* def = ppindex(table, occurence, name, false);
    if ( def == 0 )
    {
        return false;
    }
    //
    // Does it have ival values?
    //
    if ( ival >= def->m_vals.size() )
    {
        std::cerr << "ParmParse::queryval no value number"
                  << ival << " for ";
        if ( occurence ==  ParmParse::LAST )
	{
            std::cerr << "last occurence of ";
	}
        else
	{
            std::cerr << " occurence " << occurence << " of ";
	}
        std::cerr << def->m_name << '\n' << *def << '\n';
        BoxLib::Abort();
    }

    const std::string& valname = def->m_vals[ival];

    bool ok = is(valname, ptr);
    if ( !ok )
    {
        std::cerr << "ParmParse::queryval type mismatch on value number "
                  << ival << " of " << '\n';
        if ( occurence == ParmParse::LAST )
	{
            std::cerr << " last occurence of ";
	}
        else
	{
            std::cerr << " occurence number " << occurence << " of ";
	}
        std::cerr << def->m_name << '\n';
        std::cerr << " Expected an \""
                  << tok_name(ptr)
                  << "\" type  which can't be parsed from the string \""
                  << valname << "\"\n"
		  << *def << '\n';
        BoxLib::Abort();
    }
    return true;
}

template <class T>
static
void
sgetval (const ParmParse::Table& table,
	 const std::string& name,
	 T&           ptr,
	 int          ival,
	 int          occurence)
{
    if ( squeryval(table, name,ptr,ival,occurence) == 0 )
    {
        std::cerr << "ParmParse::getval ";
        if ( occurence >= 0 )
	{
            std::cerr << "occurence number "
                      << occurence
                      << " of ";
	}

        std::cerr << "ParmParse::getval(): "
                  << name
                  << " not found in table"
                  << '\n';
        ParmParse::dumpTable(std::cerr);
        BoxLib::Abort();
    }
}

template <class T>
static
bool
squeryarr (const ParmParse::Table& table,
	   const std::string& name,
	   std::vector<T>&    ptr,
	   int          start_ix,
	   int          num_val,
	   int          occurence)
{
    //
    // Get last occurrance of name in table.
    //
    const ParmParse::PP_entry *def = ppindex(table,occurence, name, false);
    if ( def == 0 )
    {
        return false;
    }
    //
    // Does it have sufficient number of values and are they all
    // the same type?
    //
    if ( num_val == ParmParse::ALL )
    {
	num_val = def->m_vals.size();
    }

    if ( num_val == 0 ) return true;

    int stop_ix = start_ix + num_val - 1;
    if ( ptr.size() < stop_ix )
    {
        ptr.resize(stop_ix + 1);
    }
    if ( stop_ix >= def->m_vals.size() )
    {
        std::cerr << "ParmParse::queryarr too many values requested for";
        if ( occurence == ParmParse::LAST )
	{
            std::cerr << " last occurence of ";
	}
        else
	{
            std::cerr << " occurence " << occurence << " of ";
	}
        std::cerr << def->m_name << '\n' << *def << '\n';
        BoxLib::Abort();
    }
    for ( int n = start_ix; n <= stop_ix; n++ )
    {
	const std::string& valname = def->m_vals[n];
	bool ok = is(valname, ptr[n]);
	if ( !ok )
	{
	    std::cerr << "ParmParse::queryarr type mismatch on value number "
		      <<  n << " of ";
	    if ( occurence == ParmParse::LAST )
	    {
		std::cerr << " last occurence of ";
	    }
	    else
	    {
		std::cerr << " occurence number " << occurence << " of ";
	    }
	    std::cerr << def->m_name << '\n';
	    std::cerr << " Expected an \""
		      << tok_name(ptr)
		      << "\" type which can't be parsed from the string \""
		      << valname << "\"\n"
		      << *def << '\n';
	    BoxLib::Abort();
	}
    }
    return true;
}

template <class T>
static
void
sgetarr (const ParmParse::Table& table,
	 const std::string&  name,
	 std::vector<T>&           ptr,
	 int          start_ix,
	 int          num_val,
	 int          occurence)
{
    if ( squeryarr(table,name,ptr,start_ix,num_val,occurence) == 0 )
    {
        std::cerr << "ParmParse::sgetarr ";
        if ( occurence >= 0 )
	{
            std::cerr << "occurence number " << occurence << " of ";
	}
        std::cerr << "ParmParse::sgetarr(): "
                  << name
                  << " not found in table"
                  << '\n';
        ParmParse::dumpTable(std::cerr);
        BoxLib::Abort();
    }
}

//
// Initialize ParmParse.
//

bool initialized = false;

void
ppinit (int argc, char** argv, const char* parfile, ParmParse::Table& table)
{
    if ( parfile != 0 )
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
	BoxLib::Error("ParmParse::prefixedName: has empty name");
    }
    if ( !m_pstack.top().empty())
    {
        return m_pstack.top() + '.' + str;
    }
    return str;
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
	BoxLib::Error("ParmParse::popPrefix: stack underflow");
    }
    m_pstack.pop();
}

std::string
ParmParse::getPrefix() const
{
    return m_pstack.top();
}

void
ParmParse::Initialize (int         argc,
		       char**      argv,
		       const char* parfile)
{
    if ( initialized )
    {
	BoxLib::Error("ParmParse::Initialize(): already initialized!");
    }
    ppinit(argc, argv, parfile, g_table);
}

ParmParse::ParmParse (const std::string& prefix)
    :
    m_table(g_table)
{
    m_pstack.push(prefix);
}

ParmParse::ParmParse (const Table& table)
    : m_table(table)
{
    m_pstack.push("");
}

ParmParse::Frame::Frame (ParmParse& pp, const std::string& pfix)
    :
    m_pp(pp), m_np(0)
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

namespace
{
bool
unused_table_entries_q (const ParmParse::Table& table)
{
    for ( const_list_iterator li = table.begin(); li != table.end(); ++li )
    {
	if ( li->m_table )
	{
	    if ( !li->m_queried )
	    {
		return true;
	    }
	    else
	    {
		return unused_table_entries_q(*li->m_table);
	    }
	}
	else if ( !li->m_queried )
	{
	    return true;
	}
    }
    return false;
}

void
finalize_table (std::string pfx, const ParmParse::Table& table)
{
    for ( const_list_iterator li = table.begin(); li != table.end(); ++li )
    {
	if ( li->m_table )
	{
	    if ( !li->m_queried )
	    {
		std::cout << "Record " << li->m_name << std::endl;
	    }
	    else
	    {
		finalize_table(pfx + "::" + li->m_name, *li->m_table);
	    }
	}
	else if ( !li->m_queried )
	{
	    std::cout << pfx << "::" << *li << std::endl;
	}
    }
}
}


void
ParmParse::Finalize ()
{
    if ( ParallelDescriptor::IOProcessor() && unused_table_entries_q(g_table))
    {
	std::cout << "Unused ParmParse Variables:\n";
	finalize_table("[TOP]", g_table);
	std::cout << "done.\n";
	//
	// First loop through and delete all queried entries.
	//
    }
    g_table.clear();
}

void
ParmParse::dumpTable (std::ostream& os)
{
    for ( const_list_iterator li=g_table.begin(); li != g_table.end(); ++li )
    {
	os << *li << std::endl;
    }
}

int
ParmParse::countval (const char* name,
                     int         n) const
{
    //
    // First find n'th occurance of name in table.
    //
    const PP_entry* def = ppindex(m_table, n, prefixedName(name), false);
    return def == 0 ? 0 : def->m_vals.size();
}

// BOOL
void
ParmParse::getkth (const char* name,
                   int         k,
                   bool&        ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
ParmParse::get (const char* name,
                bool&        ptr,
                int ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     bool&        ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  bool&        ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

// INT
void
ParmParse::getkth (const char* name,
                   int         k,
                   int&        ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
ParmParse::get (const char* name,
                int&        ptr,
                int ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     int&        ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  int&        ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
ParmParse::getktharr (const char* name,
                      int         k,
                      std::vector<int>& ptr,
                      int         start_ix,
                      int         num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name,
                   std::vector<int>& ptr,
                   int         start_ix,
                   int         num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char* name,
                        int         k,
                        std::vector<int>& ptr,
                        int         start_ix,
                        int         num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

int
ParmParse::queryarr (const char* name,
                     std::vector<int>& ptr,
                     int         start_ix,
                     int         num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

// FLOAT
void
ParmParse::getkth (const char* name,
                   int         k,
                   float&      ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
ParmParse::get (const char* name,
                float&      ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     float&      ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  float&      ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
ParmParse::getktharr (const char*   name,
                      int           k,
                      std::vector<float>& ptr,
                      int           start_ix,
                      int           num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
ParmParse::getarr (const char*   name,
                   std::vector<float>& ptr,
                   int           start_ix,
                   int           num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char*   name,
                        int           k,
                        std::vector<float>& ptr,
                        int           start_ix,
                        int           num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char*   name,
                     std::vector<float>& ptr,
                     int           start_ix,
                     int           num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

// DOUBLE
void
ParmParse::getkth (const char* name,
                   int         k,
                   double&     ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
ParmParse::get (const char* name,
                double&     ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     double&     ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  double&     ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
ParmParse::getktharr (const char*    name,
                      int            k,
                      std::vector<double>& ptr,
                      int            start_ix,
                      int            num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
ParmParse::getarr (const char*    name,
                   std::vector<double>& ptr,
                   int            start_ix,
                   int            num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char*    name,
                        int            k,
                        std::vector<double>& ptr,
                        int            start_ix,
                        int            num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char*    name,
                     std::vector<double>& ptr,
                     int            start_ix,
                     int            num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

// STRING
void
ParmParse::getkth (const char* name,
                   int         k,
                   std::string&    ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
ParmParse::get (const char* name,
                std::string&    ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     std::string&    ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  std::string&    ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
ParmParse::getktharr (const char*     name,
                      int             k,
                      std::vector<std::string>& ptr,
                      int             start_ix,
                      int             num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
ParmParse::getarr (const char*     name,
                   std::vector<std::string>& ptr,
                   int             start_ix,
                   int             num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char*     name,
                        int             k,
                        std::vector<std::string>& ptr,
                        int             start_ix,
                        int             num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char*     name,
                     std::vector<std::string>& ptr,
                     int             start_ix,
                     int             num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

// INTVECT
void
ParmParse::getkth (const char* name,
                   int         k,
                   IntVect&    ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
ParmParse::get (const char* name,
                IntVect&    ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     IntVect&    ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  IntVect&    ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
ParmParse::getktharr (const char*     name,
                      int             k,
                      std::vector<IntVect>& ptr,
                      int             start_ix,
                      int             num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
ParmParse::getarr (const char*     name,
                   std::vector<IntVect>& ptr,
                   int             start_ix,
                   int             num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char*     name,
                        int             k,
                        std::vector<IntVect>& ptr,
                        int             start_ix,
                        int             num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char*     name,
                     std::vector<IntVect>& ptr,
                     int             start_ix,
                     int             num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

// BOX
void
ParmParse::getkth (const char* name,
                   int         k,
                   Box&    ptr,
                   int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival,k);
}

void
ParmParse::get (const char* name,
                Box&    ptr,
                int         ival) const
{
    sgetval(m_table, prefixedName(name),ptr,ival, LAST);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     Box&    ptr,
                     int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  Box&    ptr,
                  int         ival) const
{
    return squeryval(m_table, prefixedName(name),ptr,ival, LAST);
}

void
ParmParse::getktharr (const char*     name,
                      int             k,
                      std::vector<Box>& ptr,
                      int             start_ix,
                      int             num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val,k);
}

void
ParmParse::getarr (const char*     name,
                   std::vector<Box>& ptr,
                   int             start_ix,
                   int             num_val) const
{
    sgetarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}

int
ParmParse::queryktharr (const char*     name,
                        int             k,
                        std::vector<Box>& ptr,
                        int             start_ix,
                        int             num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix, num_val,k);
}

int
ParmParse::queryarr (const char*     name,
                     std::vector<Box>& ptr,
                     int             start_ix,
                     int             num_val) const
{
    return squeryarr(m_table, prefixedName(name),ptr,start_ix,num_val, LAST);
}


//
// Return number of occurences of parameter name.
//

int
ParmParse::countname (const std::string& name) const
{
    int cnt = 0;
    for ( const_list_iterator li=m_table.begin(); li != m_table.end(); ++li )
    {
	if ( ppfound(prefixedName(name), *li, false) )
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
    for ( const_list_iterator li=m_table.begin(); li != m_table.end(); ++li )
    {
	if ( ppfound(prefixedName(name), *li, true) )
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
    for ( const_list_iterator li = m_table.begin(); li != m_table.end(); ++li )
    {
       if ( ppfound(prefixedName(name), *li, false))
       {
           //
           // Found an entry; mark all occurences of name as used.
           //
           for ( const_list_iterator lli = m_table.begin(); lli != m_table.end(); ++lli )
	   {
               if ( ppfound(prefixedName(name), *lli, false) )
	       {
                   lli->m_queried = true;
	       }
	   }
           return true;
       }
    }
    return false;
}

ParmParse::Record
ParmParse::getRecord (const std::string& name, int n) const
{
    const PP_entry* pe = ppindex(m_table, n, prefixedName(name), true);
    if ( pe == 0 )
    {
	std::cerr << "ParmParse::getRecord: record " << name << " not found" << std::endl;
	BoxLib::Abort();
    }
    return Record(ParmParse(*pe->m_table));
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
