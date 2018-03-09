#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif
#include <iostream>
#include "LSProblem.H"

#include "NamespaceHeader.H"

// Null Constructor
LSProblem<1>::LSProblem()
{
}

LSProblem<1>::LSProblem(const LSProblem<1> & a_lsProblem)
{
}

// Destructor
LSProblem<1>::~LSProblem()
{
}

int LSProblem<1>::recursiveCount(const int & a_degreeP)
{
  return 1;
}

void LSProblem<1>::setNumMonomials()
{
  m_numP = 1;
  m_numPLess1 = 1;
}

// Equals operator
void LSProblem<1>::operator=(const LSProblem & a_lSProblem)
{
}

void LSProblem<1>::print(ostream & a_out)const
{

}

ostream& operator<<(ostream      & a_out,
                    LSProblem<1> & a_lSProblem)
{
  a_lSProblem.print(a_out);
  return a_out;
}

#include "NamespaceFooter.H"
