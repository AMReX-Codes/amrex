
#ifndef BL_MISC_H
#define BL_MISC_H

//
// $Id: Misc.H,v 1.4 2000-10-02 20:52:36 lijewski Exp $
//

#include <SPACE.H>

/*@ManDoc:
  Returns a copy of the minimum of the two values.  Will work with
  any type T that supplies a T::operator<(const T, const T), or an
  equivalent construct, as well as a copy constructor.   In
  particular, it will work with any integral or floating-point
  intrinsic type.
*/
#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif

template <class T>
inline
T
Min (const T& a,
     const T& b)
{
    return a < b ? a : b;
}

/*@ManDoc:
  Returns a copy of the maximum of the two values.  Will work with
  any type T that supplies a T::operator>(const T, const T), or an
  equivalent construct, as well as a copy constructor.   In
  particular, it will work with any integral or floating-point
  intrinsic type.
*/

template <class T>
inline
T
Max (const T& a,
     const T& b)
{
    return a > b ? a : b;
}

/*@ManDoc:
  Returns a copy of the absolute value of the value.  Will work with
  any type T that can be compared against zero with
  T::operator>(const T, const T), or an equivalent construct, as
  well as a copy constructor, and an T::operator-() or equivalent.
  In particular, it will work with any integral or floating-point
  intrinsic type.
*/

template <class T>
inline
T
Abs (const T& a)
{
    return a > 0 ? a : -a;
}

/*@ManDoc:
  Swaps the two values.  Type T must have a copy constructor and
  an assignment operator.
*/

template <class T>
inline
void
Swap (T& a, T& b)
{
    T tmp = a;
    a = b;
    b = tmp;
}

#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif

#endif /*BL_MISC_H*/
