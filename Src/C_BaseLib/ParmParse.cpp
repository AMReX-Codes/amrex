//
// $Id: ParmParse.cpp,v 1.23 2001-07-23 18:15:21 lijewski Exp $
//

#include <iostream>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>

#include <List.H>
#include <BoxLib.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>

namespace
{
enum PPType
{
    ppDefn,
    ppOption,
    ppBool,
    ppInt,
    ppFloat,
    ppDouble,
    ppString,
    ppEQ_sign,
    ppEOF
};

struct PP_entry
{
    PP_entry (std::string&       name,
              PPType         typ,
              List<std::string>& vals);

    Array<std::string> val;
    std::string        defname;
    PPType         deftype;
    bool           queried;

    void dump (std::ostream& os) const;
};

bool
isInteger (const std::string& str,
	   int&               val)
{
    //
    // Start token scan.
    //
    char* endp = 0;
    val = int(::strtol(str.c_str(), &endp, 10));
    return *endp == 0;
}

bool
isDouble (const std::string& str,
	  double&        val)
{
   char* endp = 0;
   val = ::strtod(str.c_str(), &endp);
   return *endp == 0;
}
 

bool
isBoolean (const std::string& str,
	   bool&          val)
{
    if ( str == "true"  )
    {
        val = true;
        return true;
    }
    if ( str == "false" )
    {
        val = false;
        return true;
    }
    int int_val;
    if ( isInteger(str, int_val) )
    {
        val = int_val != 0;
        return true;
    }
    double dbl_val;
    if ( isDouble(str, dbl_val) )
    {
        val = dbl_val != 0;
        return true;
    }
    return false;
}

int             num_obj = 0;
List<PP_entry*> table;

const char* const
tok_name[] =
{
   "NAME",
   "OPTION",
   "INTEGER",
   "FLOAT",
   "DOUBLE",
   "STRING",
   "EQ_SIGN",
   "EOF"
};

//
// Simple lexical analyser.
//

enum lexState
{
    START,
    MINUS,
    SIGN,
    OPTION,
    STRING,
    QUOTED_STRING,
    INTEGER,
    START_FRACTION,
    FRACTION,
    START_EXP,
    SIGNED_EXP,
    EXP,
    PREFIX,
    SUFFIX,
    STAR
};

const char* const
state_name[] =
{
   "START",
   "MINUS",
   "SIGN",
   "OPTION",
   "STRING",
   "QUOTED_STRING",
   "INTEGER",
   "START_FRACTION",
   "FRACTION",
   "START_EXP",
   "SIGNED_EXP",
   "EXP",
   "PREFIX",
   "SUFFIX",
   "STAR"
};

void
eat_garbage (const char* str,
             int&        i,
             int         len)
{
    for (;;)
    {
        if (i < len && str[i] == '#')
        {
            while (i < len && str[i] != '\n')
                i++;
        }
        else if (i < len && isspace(str[i]))
            i++;
        else
            break;
    }
}

PPType
getToken (const char* str,
                     int&        i,
                     int         slen,
                     char*       ostr)
{
#define ERROR_MESS \
   ostr[k++] = '\0'; \
   std::cerr << "ParmParse::getToken(): invalid string = " << ostr << '\n'; \
   std::cerr << "STATE = " << state_name[state] \
             << ", next char = " << ch << '\n'; \
   std::cerr << ", rest of input = \n" << (str+i) << '\n'; \
   BoxLib::Abort()
   //
   // Eat white space and comments.
   //
   eat_garbage(str,i,slen);
   //
   // Check for end of file.
   //
   if (i >= slen || str[i] == '\0')
       return ppEOF;
   //
   // Start token scan.
   //
   lexState state = START;
   int k = 0;
   while (1)
   {
       char ch = str[i];
       switch (state)
       {
       case START:
           if (ch == '=')
           {
               ostr[k++] = ch; i++;
               ostr[k++] = 0;
               return ppEQ_sign;
           }
           else if (ch == '"')
           {
               i++;
               state = QUOTED_STRING;
           }
           else if (ch == '*')
           {
               ostr[k++] = ch; i++;
               state = STAR;
           }
           else if (isalpha(ch) || ch == '_')
           {
               ostr[k++] = ch; i++;
               state = PREFIX;
           }
           else if (ch == '-')
           {
               ostr[k++] = ch; i++;
               state = MINUS;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case MINUS:
           if (isalpha(ch) || ch == '_')
           {
               k--;
               ostr[k++] = ch; i++;
               state = OPTION;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case OPTION:
           if (isalnum(ch) || ch=='_' || ch=='.')
           {
               ostr[k++] = ch; i++;
           }
           else if (isspace(ch))
           {
               ostr[k++] = 0;
               return ppOption;
           }
           else
           {
               ERROR_MESS;
           }
           break;
       case STAR:
           if (ch == '.')
           {
               ostr[k++] = ch; i++;
               state = SUFFIX;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case PREFIX:
           if (isalnum(ch) || ch == '_')
           {
               ostr[k++] = ch; i++;
           }
           else if (ch == '.')
           {
               ostr[k++] = ch; i++;
               state = SUFFIX;
           }
           else if (isspace(ch) || ch == '=')
           {
               ostr[k++] = 0;
               return ppDefn;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case SUFFIX:
           if (isalnum(ch) || ch == '_')
           {
               ostr[k++] = ch; i++;
           }
           else if (isspace(ch) || ch == '=')
           {
               ostr[k++] = 0;
               return ppDefn;
           }
           else
           {
               ostr[k++] = ch; i++;
               state = STRING;
           }
           break;
       case QUOTED_STRING:
           if (ch == '"')
           {
               i++;
               ostr[k++] = 0;
               return ppString;
           }
           else
           {
               ostr[k++] = ch; i++;
           }
           break;
       case STRING:
           if (isspace(ch) || ch == '=')
           {
               ostr[k++] = 0;
               return ppString;
           }
           else
           {
               ostr[k++] = ch; i++;
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
ppfound (const char*    keyword,
         const std::string& key,
         const std::string& prefix)
{
    //
    // Return true if key==keyword || key == prefix.keyword.
    //
    if (!prefix.empty())
    {
        std::string tmp(prefix);
        tmp += '.';
        tmp += keyword;
        return (key == tmp);
    }
    else
    {
        return (key == keyword);
    }
}

//
// Return the index of the n'th occurence of a parameter name,
// except if n==-1, return the index of the last occurence.
// Return 0 if the specified occurence does not exist.
//

const PP_entry*
ppindex (int         n,
	 const char* name,
	 const std::string& thePrefix)
{
    PP_entry* fnd = 0;

    if (n < 0)
    {
        //
        // Search from back of list.
        //
        for (ListIterator<PP_entry*> li = table.last(); li; --li)
        {
            if (ppfound(name,li()->defname,thePrefix))
            {
                fnd = li();
                break;
            }
        }
    }
    else
    {
        for (ListIterator<PP_entry*> li(table); li; ++li)
        {
            if (ppfound(name,li()->defname,thePrefix))
            {
                fnd = li();
                if (--n < 0)
                    break;
            }
        }
        if (n >= 0)
            fnd = 0;
    }

    if (fnd)
    {
        //
        // Found an entry; mark all occurences of name as used.
        //
        for (ListIterator<PP_entry*> li(table); li; ++li)
            if (ppfound(name,li()->defname,thePrefix))
                li()->queried = true;
    }

    return fnd;
}

void
bldTable (const char*      str,
	  int              lenstr,
	  List<PP_entry*>& tab);

static void
read_file (const char*      fname,
                      List<PP_entry*>& tab)
{
    //
    // Space for input file if it exists.
    // Note: on CRAY, access requires (char*) not (const char*).
    //
    if (fname != 0 && fname[0] != 0)
    {
        FILE* pffd = fopen(fname, "rb");
        if (pffd == 0)
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
        if (!(nread == pflen))
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
        bldTable(str,pflen+1,tab);
        delete [] str;
    }
}

void
addDefn (std::string&         def,
                    List<std::string>&   val,
                    List<PP_entry*>& tab)
{
    static const std::string FileKeyword("FILE");
    //
    // Check that defn exists.
    //
    if (def.length() == 0)
    {
        val.clear();
        return;
    }
    //
    // Check that it has values.
    //
    if (val.isEmpty())
    {
        std::cerr << "ParmParse::addDefn(): no values for definition " << def;
        BoxLib::Abort();
    }
    //
    // Check if this defn is a file include directive.
    //
    if (def == FileKeyword && val.length() == 1)
    {
        //
        // Read file and add to this table.
        //
        const char* fname = val.firstElement().c_str();
        read_file(fname, tab);
    }
    else
    {
        PP_entry* pp = new PP_entry(def,ppDefn,val);
        tab.append(pp);
    }
    val.clear();
    def = std::string();
}

void
bldTable (const char*      str,
                     int              lenstr,
                     List<PP_entry*>& tab)
{
   std::string       cur_name;
   List<std::string> cur_list;
   std::string       tmp_str;
   PP_entry      *pp;

   int       i = 0;
   PPType    token;
   const int SCRATCH_STR_LEN  = 200;
   char      tokname[SCRATCH_STR_LEN];

   while (true)
   {
      token = getToken(str,i,lenstr,tokname);

      switch (token)
      {
      case ppEOF:
          addDefn(cur_name,cur_list,tab);
          return;
      case ppOption:
          addDefn(cur_name,cur_list,tab);
          tmp_str = tokname;
          pp = new PP_entry(tmp_str,ppOption,cur_list);
          tab.append(pp);
          break;
      case ppEQ_sign:
          if (cur_name.length() == 0)
              BoxLib::Abort("ParmParse::bldTable() EQ with no current defn");
          if (cur_list.isEmpty())
              //
              // First time we see equal sign, just ignore.
              //
              break;
          //
          // Read one too far, remove last name on list.
          //
          tmp_str = cur_list.lastElement();
          cur_list.removeLast();
          addDefn(cur_name,cur_list,tab);
          cur_name = tmp_str;
          break;
      case ppDefn:
          if (cur_name.length() == 0)
          {
              cur_name = tokname;
              break;
          }
          //
          // Otherwise, fall through, this may be a string.
          //
      case ppBool:
      case ppInt:
      case ppFloat:
      case ppDouble:
      case ppString:
          if (cur_name.length() == 0)
          {
              tokname[SCRATCH_STR_LEN-1] = 0;
              std::string msg("ParmParse::bldTable(): value with no defn: ");
              msg += tokname;
              BoxLib::Abort(msg.c_str());
          }
          cur_list.append(tokname);
          break;
      }
   }
}


void
PP_entry::dump (std::ostream& os) const
{
    static const char TokenInitial[] = { 'N','O','I','F','D','S','=','E' };

    char tmp[200];
    long nval = val.length();
    sprintf(tmp,
            "(%c,%1d) %15s :: ",
            TokenInitial[deftype],
            int(nval),
            defname.c_str());
    os << tmp;
    for (int i = 0; i < nval; i++)
    {
       os << " ("
          << TokenInitial[ppString]
          << ','
          << val[i]
          << ')';
    }
    os << '\n';

    if (os.fail())
        BoxLib::Error("PP_entry::dump(ostream&) failed");
}

static
int
squeryval (const char*  name,
	  const std::string& thePrefix,
	  const PPType type,
	  void*        ptr,
	  int          ival,
	  int          occurence)
{
    //
    // Get last occurrance of name in table.
    //
    const PP_entry *def = ppindex(occurence,name, thePrefix);
    if (def == 0)
        return 0;
    //
    // Does it have ival values?
    //
    if (ival >= def->val.length())
    {
        std::cerr << "ParmParse::queryval no value number"
                  << ival << " for ";
        if (occurence < 0)
            std::cerr << "last occurence of ";
        else
            std::cerr << " occurence " << occurence << " of ";
        std::cerr << def->defname << '\n';
        def->dump(std::cerr);
        BoxLib::Abort();
    }

    const std::string& valname = def->val[ival];

    bool ok = false;
    double val_dbl;
    //
    // Retrieve value.
    //
    switch (type)
    {
    case ppBool:
        ok = isBoolean(valname, *(bool*)ptr);
        break;
    case ppInt:
        ok = isInteger(valname,*(int*)ptr);
        break;
    case ppFloat:
        ok = isDouble(valname,val_dbl);
        if (ok)
            *(float*)ptr = (float) val_dbl;
        break;
    case ppDouble:
        ok = isDouble(valname,val_dbl);
        *(double*)ptr = val_dbl;
        break;
    case ppString:
        ok = true;
        *(std::string*)ptr = valname;
        break;
    default:
        BoxLib::Abort("ParmParse::queryval invalid type");
    }
    if (!ok)
    {
        std::cerr << "ParmParse::queryval type mismatch on value number "
                  << ival << " of " << '\n';
        if (occurence < 0)
            std::cerr << " last occurence of ";
        else
            std::cerr << " occurence number " << occurence << " of ";
        std::cerr << def->defname << '\n';
        std::cerr << " Expected "
                  << tok_name[type]
                  << " value = "
                  << valname << '\n';
        def->dump(std::cerr);
        BoxLib::Abort();
    }
    return 1;
}

static
void
sgetval (const char*  name,
	const std::string& thePrefix,
	const PPType type,
	void*        ptr,
	int          ival,
	int          occurence)
{
    if (squeryval(name,thePrefix,type,ptr,ival,occurence) == 0)
    {
        std::cerr << "ParmParse::getval ";
        if (occurence >= 0)
            std::cerr << "occurence number "
                      << occurence
                      << " of ";

        if (thePrefix.size())
            std::cerr << thePrefix << '.';

        std::cerr << "ParmParse::getval(): "
                  << name
                  << " not found in table"
                  << '\n';
        ParmParse::dumpTable(std::cerr);
        BoxLib::Abort();
    }
}

static
int
squeryarr (const char*  name,
		     const std::string& thePrefix,
                     const PPType type,
                     void*        ptr,
                     int          start_ix,
                     int          num_val,
                     int          occurence)
{
    //
    // Get last occurrance of name in table.
    //
    const PP_entry *def = ppindex(occurence,name, thePrefix);
    if (def == 0)
        return 0;
    //
    // Does it have sufficient number of values and are they all
    // the same type?
    //
    int stop_ix = start_ix + num_val - 1;
    if (stop_ix >= def->val.length())
    {
        std::cerr << "ParmParse::queryarr too many values requested for";
        if (occurence < 0)
            std::cerr << " last occurence of ";
        else
            std::cerr << " occurence " << occurence << " of ";
        std::cerr << def->defname << '\n';
        def->dump(std::cerr);
        BoxLib::Abort();
    }

    for (int n = start_ix; n <= stop_ix; n++)
    {
       const std::string& valname = def->val[n];
       //
       // Retrieve value.
       //
       bool ok = false;
       double val_dbl;
       switch (type)
       {
       case ppBool:
           ok = isBoolean(valname,*(bool*)ptr);
           ptr = (bool*)ptr+1;
           break;
       case ppInt:
           ok = isInteger(valname,*(int*)ptr);
           ptr = (int*)ptr+1;
           break;
       case ppFloat:
           ok = isDouble(valname,val_dbl);
           if (ok)
               *(float*)ptr = (float) val_dbl;
           ptr = (float*)ptr+1;
           break;
       case ppDouble:
           ok = isDouble(valname,*(double*)ptr);
           ptr = (double*)ptr+1;
           break;
       case ppString:
           ok = true;
           *(std::string*)ptr = valname;
           ptr = (std::string*)ptr+1;
           break;
       default:
           BoxLib::Error("ParmParse::get invalid type");
       }
       if (!ok)
       {
           std::cerr << "ParmParse::queryarr type mismatch on value number "
                <<  n << " of ";
           if (occurence < 0)
               std::cerr << " last occurence of ";
           else
               std::cerr << " occurence number " << occurence << " of ";
           std::cerr << def->defname << '\n';
           std::cerr << " Expected "
                     << tok_name[type]
                     << " value = "
                     << valname << '\n';
           def->dump(std::cerr);
           BoxLib::Abort();
       }
    }

    return 1;
}

static
void
sgetarr (const char*  name,
	const std::string& thePrefix,
	const PPType type,
	void*        ptr,
	int          start_ix,
	int          num_val,
	int          occurence)
{
    if (squeryarr(name,thePrefix,type,ptr,start_ix,num_val,occurence) == 0)
    {
        std::cerr << "ParmParse::sgetarr ";
        if (occurence >= 0)
            std::cerr << "occurence number " << occurence << " of ";
        if (thePrefix.size())
            std::cerr << thePrefix << '.';
        std::cerr << "ParmParse::sgetarr(): "
                  << name
                  << " not found in table"
                  << '\n';
        ParmParse::dumpTable(std::cerr);
        BoxLib::Abort();
    }
}

//
// Used by constructor to build table.
//
static void ppinit (int argc, char** argv, const char* parfile);

PP_entry::PP_entry (std::string&          name,
                    PPType typ,
                    List<std::string>&    vals)
    :
    val(vals.length()),
    defname(name),
    deftype(typ),
    queried(false)
{
   ListIterator<std::string> li(vals);
   for (int i = 0; li; i++, ++li)
      val[i] = vals[li];
}

//
// Initialize ParmParse.
//

void
ppinit (int argc, char** argv, const char* parfile)
{
    if (parfile != 0)
       read_file(parfile,table);

    if (argc > 0)
    {
        std::string argstr;
        const char SPACE = ' ';
        for (int i = 0; i < argc; i++)
        {
            argstr += argv[i];
            argstr += SPACE;
        }
        List<PP_entry*> arg_table;
        bldTable(argstr.c_str(), argstr.length()+1, arg_table);
        //
        // Append arg_table to end of existing table.
        //
        table.catenate(arg_table);
    }
}
}

ParmParse::ParmParse (int         argc,
                      char**      argv,
                      const char* prefix,
                      const char* parfile)
{
    if (table.length() > 0)
       BoxLib::Error("ParmParse::ParmParse(): table already built");
    num_obj++;
    if (prefix != 0)
       thePrefix = prefix;
    ppinit(argc, argv, parfile);
}

ParmParse::ParmParse (const char* prefix)
{
    num_obj++;
    if (prefix != 0)
       thePrefix = prefix;
}

ParmParse::~ParmParse ()
{
    if (--num_obj == 0)
    {
        //
        // First loop through and delete all queried entries.
        //
        for (ListIterator<PP_entry*> li(table); li; )
        {
            if (li()->queried)
            {
                delete table[li];

                ListIterator<PP_entry*> tmp = li++;
                
                table.remove(tmp);
            }
            else
            {
                ++li;
            }
        }

        if (!table.isEmpty())
        {
            //
            // Now spit out unused entries.
            //
            if (ParallelDescriptor::IOProcessor())
                std::cout << "\nUnused ParmParse Variables: ";

            for (ListIterator<PP_entry*> li(table); li; ++li)
            {
                if (ParallelDescriptor::IOProcessor())
                    std::cout << li()->defname << ' ';

                delete table[li];
            }

            if (ParallelDescriptor::IOProcessor())
                std::cout << '\n' << '\n';

            table.clear();
        }
    }
}

void
ParmParse::dumpTable (std::ostream& os)
{
   for (ListIterator<PP_entry*> li(table); li; ++li)
      li()->dump(os);
}

int
ParmParse::countval (const char* name,
                     int         n)
{
    //
    // First find n'th occurance of name in table.
    //
    const PP_entry* def = ppindex(n,name,thePrefix);
    return def == 0 ? 0 : def->val.length();
}


void
ParmParse::getkth (const char* name,
                   int         k,
                   bool&        ptr,
                   int         ival)
{
    sgetval(name,thePrefix,ppBool,&ptr,ival,k);
}

void
ParmParse::get (const char* name,
                bool&        ptr,
                int ival)
{
    sgetval(name,thePrefix,ppBool,&ptr,ival,-1);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     bool&        ptr,
                     int         ival)
{
    return squeryval(name,thePrefix,ppBool,&ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  bool&        ptr,
                  int         ival)
{
    return squeryval(name,thePrefix,ppBool,&ptr,ival,-1);
}

void
ParmParse::getkth (const char* name,
                   int         k,
                   int&        ptr,
                   int         ival)
{
    sgetval(name,thePrefix,ppInt,&ptr,ival,k);
}

void
ParmParse::get (const char* name,
                int&        ptr,
                int ival)
{
    sgetval(name,thePrefix,ppInt,&ptr,ival,-1);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     int&        ptr,
                     int         ival)
{
    return squeryval(name,thePrefix,ppInt,&ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  int&        ptr,
                  int         ival)
{
    return squeryval(name,thePrefix,ppInt,&ptr,ival,-1);
}

void
ParmParse::getkth (const char* name,
                   int         k,
                   float&      ptr,
                   int         ival)
{
    sgetval(name,thePrefix,ppFloat,&ptr,ival,k);
}

void
ParmParse::get (const char* name,
                float&      ptr,
                int         ival)
{
    sgetval(name,thePrefix,ppFloat,&ptr,ival,-1);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     float&      ptr,
                     int         ival)
{
    return squeryval(name,thePrefix,ppFloat,&ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  float&      ptr,
                  int         ival)
{
    return squeryval(name,thePrefix,ppFloat,&ptr,ival,-1);
}

void
ParmParse::getkth (const char* name,
                   int         k,
                   double&     ptr,
                   int         ival)
{
    sgetval(name,thePrefix,ppDouble,&ptr,ival,k);
}

void
ParmParse::get (const char* name,
                double&     ptr,
                int         ival)
{
    sgetval(name,thePrefix,ppDouble,&ptr,ival,-1);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     double&     ptr,
                     int         ival)
{
    return squeryval(name,thePrefix,ppDouble,&ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  double&     ptr,
                  int         ival)
{
    return squeryval(name,thePrefix,ppDouble,&ptr,ival,-1);
}

void
ParmParse::getkth (const char* name,
                   int         k,
                   std::string&    ptr,
                   int         ival)
{
    sgetval(name,thePrefix,ppString,&ptr,ival,k);
}

void
ParmParse::get (const char* name,
                std::string&    ptr,
                int         ival)
{
    sgetval(name,thePrefix,ppString,&ptr,ival,-1);
}

int
ParmParse::querykth (const char* name,
                     int         k,
                     std::string&    ptr,
                     int         ival)
{
    return squeryval(name,thePrefix,ppString,&ptr,ival,k);
}

int
ParmParse::query (const char* name,
                  std::string&    ptr,
                  int         ival)
{
    return squeryval(name,thePrefix,ppString,&ptr,ival,-1);
}

void
ParmParse::getktharr (const char* name,
                      int         k,
                      Array<int>& ptr,
                      int         start_ix,
                      int         num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    sgetarr(name,thePrefix,ppInt,&(ptr[0]),start_ix,num_val,k);
}

void
ParmParse::getarr (const char* name,
                   Array<int>& ptr,
                   int         start_ix,
                   int         num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    sgetarr(name,thePrefix,ppInt,&(ptr[0]),start_ix,num_val,-1);
}

int
ParmParse::queryktharr (const char* name,
                        int         k,
                        Array<int>& ptr,
                        int         start_ix,
                        int         num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    return squeryarr(name,thePrefix,ppInt,&(ptr[0]),start_ix,num_val,k);
}

int
ParmParse::queryarr (const char* name,
                     Array<int>& ptr,
                     int         start_ix,
                     int         num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    return squeryarr(name,thePrefix,ppInt,&(ptr[0]),start_ix,num_val,-1);
}

void
ParmParse::getktharr (const char*   name,
                      int           k,
                      Array<float>& ptr,
                      int           start_ix,
                      int           num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    sgetarr(name,thePrefix,ppFloat,&(ptr[0]),start_ix,num_val,k);
}

void
ParmParse::getarr (const char*   name,
                   Array<float>& ptr,
                   int           start_ix,
                   int           num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    sgetarr(name,thePrefix,ppFloat,&(ptr[0]),start_ix,num_val,-1);
}

int
ParmParse::queryktharr (const char*   name,
                        int           k,
                        Array<float>& ptr,
                        int           start_ix,
                        int           num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    return squeryarr(name,thePrefix,ppFloat,&(ptr[0]),start_ix, num_val,k);
}

int
ParmParse::queryarr (const char*   name,
                     Array<float>& ptr,
                     int           start_ix,
                     int           num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    return squeryarr(name,thePrefix,ppFloat,&(ptr[0]),start_ix,num_val,-1);
}

void
ParmParse::getktharr (const char*    name,
                      int            k,
                      Array<double>& ptr,
                      int            start_ix,
                      int            num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    sgetarr(name,thePrefix,ppDouble,&(ptr[0]),start_ix,num_val,k);
}

void
ParmParse::getarr (const char*    name,
                   Array<double>& ptr,
                   int            start_ix,
                   int            num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    sgetarr(name,thePrefix,ppDouble,&(ptr[0]),start_ix,num_val,-1);
}

int
ParmParse::queryktharr (const char*    name,
                        int            k,
                        Array<double>& ptr,
                        int            start_ix,
                        int            num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    return squeryarr(name,thePrefix,ppDouble,&(ptr[0]),start_ix, num_val,k);
}

int
ParmParse::queryarr (const char*    name,
                     Array<double>& ptr,
                     int            start_ix,
                     int            num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    return squeryarr(name,thePrefix,ppDouble,&(ptr[0]),start_ix,num_val,-1);
}

void
ParmParse::getktharr (const char*     name,
                      int             k,
                      Array<std::string>& ptr,
                      int             start_ix,
                      int             num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    sgetarr(name,thePrefix,ppString,&(ptr[0]),start_ix,num_val,k);
}

void
ParmParse::getarr (const char*     name,
                   Array<std::string>& ptr,
                   int             start_ix,
                   int             num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    sgetarr(name,thePrefix,ppString,&(ptr[0]),start_ix,num_val,-1);
}

int
ParmParse::queryktharr (const char*     name,
                        int             k,
                        Array<std::string>& ptr,
                        int             start_ix,
                        int             num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    return squeryarr(name,thePrefix,ppString,&(ptr[0]),start_ix, num_val,k);
}

int
ParmParse::queryarr (const char*     name,
                     Array<std::string>& ptr,
                     int             start_ix,
                     int             num_val)
{
    if (ptr.length() < num_val)
        ptr.resize(num_val);
    return squeryarr(name,thePrefix,ppString,&(ptr[0]),start_ix,num_val,-1);
}

int
ParmParse::countname (const std::string& name)
{
    return countname(name.c_str());
}

//
// Return number of occurences of parameter name.
//

int
ParmParse::countname (const char* name)
{
    int cnt = 0;
    for (ListIterator<PP_entry*> li(table); li; ++li)
       if (ppfound(name,li()->defname,thePrefix))
           cnt++;
    return cnt;
}

//
// Return true if name in table.
//

bool
ParmParse::contains (const char* name)
{
    for (ListIterator<PP_entry*> li(table); li; ++li)
    {
       if (ppfound(name,li()->defname,thePrefix))
       {
           //
           // Found an entry; mark all occurences of name as used.
           //
           for (ListIterator<PP_entry*> lli(table); lli; ++lli)
               if (ppfound(name,lli()->defname,thePrefix))
                   lli()->queried = true;
           return true;
       }
    }
    return false;
}

