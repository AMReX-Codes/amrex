//
// $Id: aString.cpp,v 1.11 2001-07-17 23:02:30 lijewski Exp $
//

#include <cctype>

#include <BLassert.H>
#include <aString.H>

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

StringRep::StringRep (int _len)
{
    bufferlength = _len;
    s = new char [bufferlength];
}

StringRep::~StringRep ()
{
    delete [] s;
    s = 0;
}

bool
aString::isNull () const
{
    return len == 0;
}

int
aString::length () const
{
    return len;
}

double
aString::toDouble () const
{
    return len == 0 ? 0 : atof(p->s);
}

int
aString::toInteger () const
{
    return len == 0 ? 0 : atoi(p->s);
}

long
aString::toLong () const
{
    return len == 0 ? 0 : atol(p->s);
}

const char*
aString::c_str () const
{
    return p->s;
}

char
aString::operator[] (int index) const
{
    BL_ASSERT(index >=0 && index < len);
    return p->s[index];
}

aString
operator+ (const aString& left,
           const aString& right)
{
    aString result(left);
    return result += right;
}

bool
operator< (const aString& left,
           const aString& right)
{
    return ::strcmp(left.c_str(), right.c_str()) < 0;
}

bool
operator<= (const aString& left,
            const aString& right)
{
    return ::strcmp(left.c_str(), right.c_str()) <= 0;
}

bool
operator!= (const aString& left,
            const aString& right)
{
    return ::strcmp(left.c_str(), right.c_str()) != 0;
}

bool
operator== (const aString& left,
            const aString& right)
{
    return ::strcmp(left.c_str(), right.c_str()) == 0;
}

bool
operator>= (const aString& left,
            const aString& right)
{
    return ::strcmp(left.c_str(), right.c_str()) >= 0;
}

bool
operator>  (const aString& left,
            const aString& right)
{
    return ::strcmp(left.c_str(), right.c_str()) > 0;
}

bool
operator< (const aString& left,
           const char*    right)
{
    return ::strcmp(left.c_str(), right) < 0;
}

bool
operator<= (const aString& left,
            const char*    right)
{
    return ::strcmp(left.c_str(), right) <= 0;
}

bool
operator!= (const aString& left,
            const char*    right)
{
    return ::strcmp(left.c_str(), right) != 0;
}

bool
operator== (const aString& left,
            const char*    right)
{
    return ::strcmp(left.c_str(), right) == 0;
}

bool
operator>= (const aString& left,
            const char*    right)
{
    return ::strcmp(left.c_str(), right) >= 0;
}

bool
operator>  (const aString& left,
            const char*    right)
{
    return ::strcmp(left.c_str(), right) > 0;
}

bool
operator< (const char*    left,
           const aString& right)
{
    return ::strcmp(left, right.c_str()) < 0;
}

bool
operator<= (const char*    left,
            const aString& right)
{
    return ::strcmp(left, right.c_str()) <= 0;
}

bool
operator!= (const char*    left,
            const aString& right)
{
    return ::strcmp(left, right.c_str()) != 0;
}

bool
operator== (const char*    left,
            const aString& right)
{
    return ::strcmp(left, right.c_str()) == 0;
}

bool
operator>= (const char*    left,
            const aString& right)
{
    return ::strcmp(left, right.c_str()) >= 0;
}

bool
operator>  (const char*    left,
            const aString& right)
{
    return ::strcmp(left, right.c_str()) > 0;
}

void
StringRep::resize (int n)
{
    if (n > bufferlength)
    {
        char* ns = new char[n];
        ::memcpy(ns,s,bufferlength);
        bufferlength = n;
        delete [] s;
        s = ns;
    }
}

aString::aString ()
    : p(new StringRep(1)),
      len(0)
{
    p->s[0] = 0;
}

aString::aString (char c)
    : p(new StringRep(2)),
      len(1)
{
    p->s[0] = c;
    p->s[1] = 0;
    if (c == '\0')
        len = 0;
}

aString::aString (int size)
    : p(new StringRep(size+1)),
      len(0)
{
    ::memset(p->s,'\0',p->bufferlength);
}

aString::aString (const char* initialtext)
{
    BL_ASSERT(initialtext != 0);
    len = ::strlen(initialtext);
    p = new StringRep(len + 1);
    ::memcpy(p->s,initialtext,len+1);
}

aString::aString (const aString& initialstring)
    : p(initialstring.p),
      len(initialstring.len)
{}

aString&
aString::operator= (const aString& rhs)
{
    p   = rhs.p;
    len = rhs.len;
    return *this;
}

aString&
aString::operator+= (const aString& val)
{
    copyModify();
    int clen = length() + val.length();
    p->resize(clen+1);
    ::memcpy(&(p->s[len]),val.p->s, val.length()+1);
    len = clen;
    return *this;
}

aString&
aString::operator+= (const char* s)
{
    BL_ASSERT(s != 0);
    copyModify();
    int slen = ::strlen(s);
    int clen = length() + slen;
    p->resize(clen+1);
    ::memcpy(&(p->s[len]),s, slen+1);
    len = clen;
    return *this;
}

aString&
aString::operator+= (char c)
{
    if (!(c == '\0'))
    {
        copyModify();
        p->resize(len+2);
        p->s[len++] = c;
        p->s[len]   = 0;
    }
    return *this;
}

char&
aString::operator[] (int index)
{
    BL_ASSERT(index >= 0 && index < len);
    copyModify();
    return p->s[index];
}

void
aString::copyModify ()
{
    if (!p.unique())
    {
        StringRep* np = new StringRep(len+1);
        ::memcpy(np->s,p->s,len+1);
        p = np;
    }
}

aString&
aString::toUpper ()
{
    copyModify();
    for (char *pp = p->s; pp != 0; pp++)
        *pp = toupper(*pp);
    return *this;
}

aString&
aString::toLower ()
{
    copyModify();
    for (char *pp = p->s; pp != 0; pp++)
        *pp = tolower(*pp);
    return *this;
}

std::istream&
operator>> (std::istream& is,
            aString&      str)
{
    const int BufferSize = 128;
    char buf[BufferSize + 1];
    int index = 0;
    //
    // Nullify str.
    //
    str = "";
    //
    // Eat leading whitespace.
    //
    char c;
    do { is.get(c); } while (is.good() && isspace(c));
    buf[index++] = c;
    //
    // Read until next whitespace.
    //
    while (is.get(c) && !isspace(c))
    {
        buf[index++] = c;
        if (index == BufferSize)
        {
            buf[BufferSize] = 0;
            str += buf;
            index = 0;
        }
    }
    is.putback(c);
    buf[index] = 0;
    str += buf;
    if (is.fail())
        BoxLib::Abort("operator>>(istream&,aString&) failed");
    return is;
}

std::ostream&
operator<< (std::ostream&  out,
            const aString& str)
{
    out.write(str.c_str(), str.len);
    if (out.fail())
        BoxLib::Error("operator<<(ostream&,aString&) failed");
    return out;
}

std::istream&
aString::getline (std::istream& is)
{
    char      c;
    const int BufferSize = 256;
    char      buf[BufferSize + 1];
    int       index = 0;

    *this = "";
    //
    // Get those characters.
    // We read the newline but don't add it to the string.
    //
    while (is.get(c))
    {
        if (c == '\n')
            break;

        buf[index++] = c;

        if (index == BufferSize)
        {
            buf[BufferSize] = 0;
            *this += buf;
            index = 0;
        }
    }
    buf[index] = 0;
    *this += buf;

    if (!(is || is.eof()))
        BoxLib::Abort("aString::getline(istream&) failed");

    return is;
}

std::vector<aString>
aString::tokenize (const aString& separators) const
{
    std::vector<char*> ptr;
    //
    // Make copy of line that we can modify.
    //
    char* line = new char[length()+1];

    (void) strcpy(line, c_str());

    char* token = 0;

    if (!((token = strtok(line, separators.c_str())) == 0))
    {
        ptr.push_back(token);
        while (!((token = strtok(0, separators.c_str())) == 0))
            ptr.push_back(token);
    }

    std::vector<aString> tokens(ptr.size());

    for (unsigned int i = 1; i < ptr.size(); i++)
    {
        char* p = ptr[i];

        while (strchr(separators.c_str(), *(p-1)) != 0)
            *--p = 0;
    }

    for (unsigned int i = 0; i < ptr.size(); i++)
        tokens[i] = ptr[i];

    delete line;

    return tokens;
}

#ifdef BL_NAMESPACE
}
#endif

