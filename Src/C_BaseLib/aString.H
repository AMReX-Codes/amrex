#ifndef BL_ASTRING_H
#define BL_ASTRING_H

//
// $Id: aString.H,v 1.16 2001-07-20 17:01:45 car Exp $
//

#include <cstddef>
#include <string>
#include <iostream>

class aString;


//
//@Man:
//@Memo: A String Class
/*@Doc:

  The class aString is used to store and manipulate character strings.  It
  has an efficient underlying storage mechanism and some useful string
  manipulation operations.

  The aString class is implemented using a character array and reference
  count. Two aStrings may reference the same underlying character array with
  a reference count of two.  When an aString copy constructor or copy
  operator is applied the reference count on the underlying character array
  is incremented but the actual string is not copied.  That is, copying an
  aString is an inexpensive operation.  When an aString is destructed, the
  reference count is decremented.  The underlying character array is deleted
  only when the reference count goes to zero.  Any operator that modifies an
  aString will make its own copy of the character array before the
  modification, unless it's the sole owner of the character array in the
  aString.
*/

class aString
    : public std::string
{
public:
    //
    //@ManDoc: Constructs an empty string.
    //
    aString ();

    /*@ManDoc: Constructs an aString containing the single character c.
               If c is the null character then this will be the empty string.
    */
    explicit aString (char c);

    /*@ManDoc: Constructs an empty aString but allocates enough space to hold
               a character string of length len.  This may be useful when
               reading in very long lines with the `getline' function; i.e.
               it may be slightly more efficient to allocate a very large
               aString once and then use `getline' than to use `getline'
               on an empty string.  In general though, you won't notice any
               difference and you really shouldn't use this constructor.
    */
    aString(const std::string& str);
    
    explicit aString (size_type len);

    /*@ManDoc: Constructs an aString initialized to the character string s.
               It is an error is `s' is the null string.
    */
    aString (const char* s);
    /*@ManDoc: Converts all characters in this aString to upper case.
               Returns a reference to self to allow for operator chaining.
    */
    aString& toUpper ();

    /*@ManDoc: Converts all characters in this aString to lower case.
               Returns a reference to self to allow for operator chaining.
    */
    aString& toLower ();

    /*@ManDoc: Read the next line from the input stream strm and store it
               as the value of the string.  The delimiter for a line of text
               is the newline character.  The newline is extracted from the
               istream, but it is NOT stored in the aString.  There is no
               limit to the length of string that can be extracted.
    */
    std::istream& getline (std::istream& strm);

    /*@ManDoc: Returns the number of characters stored in this aString.
               This does not include the terminating null character.
    */
    size_type length () const;

    //
    //@ManDoc: Returns true if this is the null string.
    //
    bool isNull () const;

    /*@ManDoc: Returns a reference to the kth character in the string.  An
               error occurs if k < 0 or k >= length().
    */
    char& operator [] (size_type k);

    /*@ManDoc: Returns kth character in the string.  An error occurs
               if k < 0 or k >= length().
    */
    char operator[] (size_type k) const;


    /*@ManDoc: Returns the value of the string as a double.  In the case
               when the string is empty or when the initial characters
               in the string are strictly alphabetic characters, this
               returns zero.  It's equivalent to `atof(c\_str())'.
    */
    double toDouble () const;

    /*@ManDoc: Returns the value of the string as a integer.  In the case
               when the string is empty or when the initial characters
               in the string are strictly alphabetic characters, this
               returns zero.  It's equivalent to `atoi(c\_str())'.
    */
    int toInteger () const;

    /*@ManDoc: Returns the value of the string as a long.  In the case
               when the string is empty or when the initial characters
               in the string are strictly alphabetic characters, this
               returns zero.  It's equivalent to `atol(c\_str())'.
    */
    long toLong () const;

};

#endif /*BL_ASTRING_H*/
