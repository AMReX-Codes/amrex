//BL_COPYRIGHT_NOTICE

//
// $Id: ParmParse.cpp,v 1.4 1997-09-24 22:06:45 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
using std::cin;
using std::cout;
using std::cerr;
#else
#include <iostream.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#endif

#include <BoxLib.H>
#include <ParmParse.H>

int             ParmParse::xargc   = -1;
int             ParmParse::num_obj = 0;
char**          ParmParse::xargv   = 0;
List<PP_entry*> ParmParse::table;

static const char* tok_name[] =
{
   "NAME", "OPTION", "INTEGER", "FLOAT", "DOUBLE", "STRING", "EQ_SIGN", "EOF"
};

PP_entry::PP_entry (aString&          name,
                    ParmParse::PPType typ,
                    List<aString>&    vals)
          : defname(name),
            deftype(typ),
            val(vals.length())
{
   ListIterator<aString> li(vals);
   for (int i = 0; li; i++, ++li)
      val[i] = vals[li];
}

void
PP_entry::dump (ostream& os) const
{
    static const char TokenInitial[] =
    {
        'N', 'O', 'I', 'F', 'D', 'S', '=', 'E'
    };

    char tmp[200];
    long nval = val.length();
    sprintf(tmp,
            "(%c,%1d) %15s :: ",
            TokenInitial[deftype],
            int(nval),
            defname.c_str());
    os << tmp;
    for (int i = 0; i < nval; i++)
       os << " (" << TokenInitial[ParmParse::ppString] << ',' << val[i] << ')';
    os << '\n';

    if (os.fail())
        BoxLib::Error("PP_entry::dump(ostream&) failed");
}

ParmParse::ParmParse (int         argc,
                      char**      argv,
                      const char* prefix,
                      const char* parfile)
{
    if (table.length() > 0)
       BoxLib::Abort("ParmParse::ParmParse(): table already built");
    num_obj++;
    xargc = argc;
    xargv = argv;
    if (prefix != 0)
       thePrefix = prefix;
    ppinit(parfile);
}

ParmParse::ParmParse (const char* prefix)
{
    num_obj++;
    if (xargc < 0)
        BoxLib::Error("ParmParse::ParmParse(): class not properly initialized");
    if (prefix != 0)
       thePrefix = prefix;
}

void ParmParse::dumpTable (ostream& os)
{
   for (ListIterator<PP_entry*> li(table); li; ++li)
      li()->dump(os);
}

//
// Initialize ParmParse.
//

void
ParmParse::ppinit (const char* parfile)
{
    if (xargc < 0)
    {
        cerr << "ParmParse::ppinit(): xargc="
             << xargc
             << " not initiated in main";
        BoxLib::Abort();
    }

    if (parfile != 0)
       read_file(parfile,table);

    if (xargc > 0)
    {
        aString argstr;
        const char SPACE = ' ';
        for (int i = 0; i < xargc; i++)
        {
            argstr += xargv[i];
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

ParmParse::~ParmParse ()
{
   if (--num_obj == 0)
   {
      for (ListIterator<PP_entry*> li(table); li; ++li)
         delete table[li];
      table.clear();
   }
}

//
// Keyword aware string comparison.
//

static
bool
ppfound (const char*    keyword,
         const aString& key,
         const aString& prefix)
{
    //
    // Return true if key==keyword || key == prefix.keyword.
    //
    if (!prefix.isNull())
    {
        aString tmp(prefix);
        tmp += '.';
        tmp += keyword;
        return (key == tmp);
    }
    else
        return (key == keyword);
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
           return true;
    }
    return false;
}

//
// Return the index of the n'th occurence of a parameter name,
// except if n==-1, return the index of the last occurence.
// Return 0 if the specified occurence does not exist.
//

const PP_entry*
ParmParse::ppindex (int         n,
                    const char* name) const
{
    const PP_entry* fnd = 0;
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
                n--;
                if (n < 0)
                    break;
            }
        }
        if (n >= 0)
            fnd = 0;
    }
    return fnd;
}

void
ParmParse::getval (const char*  name,
                   const PPType type,
                   void*        ptr,
                   int          ival,
                   int          occurence)
{

    if (queryval(name,type,ptr,ival,occurence) == 0)
    {
        cerr << "ParmParse::getval ";
        if (occurence >= 0)
            cerr << "occurence number " << occurence << " of ";
        if (!thePrefix.isNull())
            cerr << thePrefix << '.';
        cerr << "ParmParse::getval(): " << name
             << " not found in table"   << '\n';
        dumpTable(cerr);
        BoxLib::Abort();
    }
}

int
ParmParse::queryval (const char*  name,
                     const PPType type,
                     void*        ptr,
                     int          ival,
                     int          occurence)
{
    //
    // Get last occurrance of name in table.
    //
    const PP_entry *def = ppindex(occurence,name);
    if (def == 0)
        return 0;
    //
    // Does it have ival values?
    //
    if (ival >= def->val.length())
    {
        cerr << "ParmParse::queryval no value number"
             << ival << " for ";
        if (occurence < 0)
            cerr << "last occurence of ";
        else
            cerr << " occurence " << occurence << " of ";
        cerr << def->defname << '\n';
        def->dump(cerr);
        abort();
    }

    const aString& valname = def->val[ival];

    int ok;
    double val_dbl;
    //
    // Retrieve value.
    //
    switch (type)
    {
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
        *(aString*)ptr = valname;
        break;
    default:
        BoxLib::Abort("ParmParse::queryval invalid type");
    }
    if (!ok)
    {
        cerr << "ParmParse::queryval type mismatch on value number "
             << ival << " of " << '\n';
        if (occurence < 0)
            cerr << " last occurence of ";
        else
            cerr << " occurence number " << occurence << " of ";
        cerr << def->defname << '\n';
        cerr << " Expected " << tok_name[type]
             << " value = " << valname << '\n';
        def->dump(cerr);
        abort();
    }
    return 1;
}

void
ParmParse::getarr (const char*  name,
                   const PPType type,
                   void*        ptr,
                   int          start_ix,
                   int          num_val,
                   int          occurence)
{
    if (queryarr(name,type,ptr,start_ix,num_val,occurence) == 0)
    {
        cerr << "ParmParse::getarr ";
        if (occurence >= 0)
            cerr << "occurence number " << occurence << " of ";
        if (!thePrefix.isNull())
            cerr << thePrefix << '.';
        cerr << "ParmParse::getarr(): " << name
             << " not found in table"   << '\n';
        dumpTable(cerr);
        BoxLib::Abort();
    }
}

int
ParmParse::queryarr (const char*  name,
                     const PPType type,
                     void*        ptr,
                     int          start_ix,
                     int          num_val,
                     int          occurence)
{
    //
    // Get last occurrance of name in table.
    //
    const PP_entry *def = ppindex(occurence,name);
    if (def == 0)
        return 0;
    //
    // Does it have sufficient number of values and are they all
    // the same type?
    //
    int stop_ix = start_ix + num_val - 1;
    if (stop_ix >= def->val.length())
    {
        cerr << "ParmParse::queryarr too many values requested for";
        if (occurence < 0)
            cerr << " last occurence of ";
        else
            cerr << " occurence " << occurence << " of ";
        cerr << def->defname << '\n';
        def->dump(cerr);
        abort();
    }

    for (int n = start_ix; n <= stop_ix; n++)
    {
       const aString& valname = def->val[n];
       //
       // Retrieve value.
       //
       int ok = false;
       double val_dbl;
       switch (type)
       {
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
           *(aString*)ptr = valname;
           ptr = (aString*)ptr+1;
           break;
       default:
           cerr << "ParmParse::get invalid type" << '\n';
           abort();
       }
       if (!ok)
       {
           cerr << "ParmParse::queryarr type mismatch on value number "
                <<  n << " of ";
           if (occurence < 0)
               cerr << " last occurence of ";
           else
               cerr << " occurence number " << occurence << " of ";
           cerr << def->defname << '\n';
           cerr << " Expected " << tok_name[type]
                << " value = " << valname << '\n';
           def->dump(cerr);
           abort();
       }
    }

    return 1;
}

void
ParmParse::read_file (const char*      fname,
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
            cerr << "ParmParse::read_file(): couldn't open \""
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
        if (str == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
        memset(str,0,pflen+1);
        int nread = fread(str, 1, pflen, pffd);
        if (!(nread == pflen))
        {
            cerr << "ParmParse::read_file(): fread() only "
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

static
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

static const char* state_name[] =
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

ParmParse::PPType
ParmParse::getToken (const char* str,
                     int&        i,
                     int         slen,
                     char*       ostr)
{
#define ERROR_MESS \
   ostr[k++] = '\0'; \
   cerr << "ParmParse::getToken(): invalid string = " << ostr << '\n'; \
   cerr << "STATE = " << state_name[state] \
        << ", next char = " << ch << '\n'; \
   cerr << ", rest of input = \n" << (str+i) << '\n'; \
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

void
ParmParse::bldTable (const char*      str,
                     int              lenstr,
                     List<PP_entry*>& tab)
{
   aString       cur_name;
   List<aString> cur_list;
   aString       tmp_str;
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
          if (pp == 0)
              BoxLib::OutOfMemory(__FILE__, __LINE__);
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
      case ppInt:
      case ppFloat:
      case ppDouble:
      case ppString:
          if (cur_name.length() == 0)
              BoxLib::Abort("ParmParse::bldTable(): value with no defn");
          cur_list.append(tokname);
          break;
      }
   }
}

void
ParmParse::addDefn (aString&         def,
                    List<aString>&   val,
                    List<PP_entry*>& tab)
{
    static const aString FileKeyword("FILE");
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
        cerr << "ParmParse::addDefn(): no values for definition " << def;
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
        if (pp == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);
        tab.append(pp);
    }
    val.clear();
    def = aString();
}
