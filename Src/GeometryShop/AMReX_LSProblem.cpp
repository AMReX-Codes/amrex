#include "AMReX_LSProblem.H"

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


