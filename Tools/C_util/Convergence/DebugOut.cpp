/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#include "DebugOut.H"
#include "AMReX_Print.H"
#include "AMReX_BoxArray.H"
#include "AMReX_parstream.H"
#include "AMReX_PlotFileUtil.H"
#include <iomanip>
#include <string>
#include <sstream>

using namespace amrex;
using namespace std;

string convertIntDD(int number)
{
  stringstream ss;//create a stringstream
  ss << number;//add number to the stream
  return ss.str();//return a string with the contents of the stream
}

void printDBL(const BoxArray* a_fabPtr)
{
  amrex::Print() << "BoxArray contains:"  << *a_fabPtr << endl;
    
}
void printBA(const BoxArray* a_fabPtr)
{
  amrex::Print() << "BoxArray contains:"  << *a_fabPtr << endl;
    
}
void printFAB(const FArrayBox* a_fabPtr)
{
  const FArrayBox& fab = *a_fabPtr;
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
  {
    amrex::Print() << "\t" << bit() ;
    for (int ivar = 0; ivar < fab.nComp(); ivar++)
    {
      amrex::Print() << "\t" << fab(bit(),ivar);
    }
    amrex::Print() << "\n";
  }
}

void printBFR(const BaseFab<Real>* a_fabPtr)
{
  const BaseFab<Real>& fab = *a_fabPtr;
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
  {
    amrex::Print() << "\t" << bit() ;
    for (int ivar = 0; ivar < fab.nComp(); ivar++)
    {
      amrex::Print() << "\t" << fab(bit(),ivar);
    }
    amrex::Print() << "\n";
  }
}
void printBFI(const BaseFab<int>* a_fabPtr)
{
  const BaseFab<int>& fab = *a_fabPtr;
  BoxIterator bit(fab.box());
  for (bit.reset(); bit.ok(); ++bit)
  {
    amrex::Print() << "\t" << bit() ;
    for (int ivar = 0; ivar < fab.nComp(); ivar++)
    {
      amrex::Print() << "\t" << fab(bit(),ivar);
    }
    amrex::Print() << "\n";
  }
}

void printBL(const BoxArray* a_dblInPtr)
{
  amrex::Print() << "BoxArray: ";
  for(int ibox = 0; ibox < a_dblInPtr->size(); ibox++)
  {
    amrex::Print() << (*a_dblInPtr)[ibox] << "   ";
  }
  amrex::Print() << "\n";
}


void printBox(const Box* a_boxPtr)
{
  amrex::Print() << "Box:" << *a_boxPtr << "\n";
}

void 
viewMF(const MultiFab* a_data)
{
#if 0
  Vector<string> names(a_data->nComp());
  for(int icomp = 0; icomp < a_data->nComp(); icomp++)
  {
    names[icomp] = string("var_") + convertIntDD(icomp);
  }
  string filename("debug_file.plt");
  Geometry geom(a_data->getDomain());
  Real time = 0; int level_step = 0;  

  WriteSingleLevelPlotfile(filename, *a_data, names, geom, time, level_step);

  string command = "visit -o " + filename + string("/Header");
  int ret = std::system(command.c_str());
  amrex::Print() << "data output to " << filename << ".  Visit was called and got return value " << ret << endl;
#endif
}


void 
maxMinMF(const MultiFab* a_data)
{
  amrex::Print() << "max min for multifab "  <<endl;
  for(int icomp = 0; icomp < a_data->nComp(); icomp++)
  {
    Real maxVal = a_data->max(icomp);
    Real minVal = a_data->min(icomp);
    amrex::Print() << "comp = " << icomp << ", max = " << maxVal << ", min = " << minVal << endl;
  }
}


void 
printMF(const MultiFab* a_data)
{
  amrex::Print() << "data for multifab "  <<endl;
  int ibox = 0;
  for(MFIter mfi(*a_data); mfi.isValid(); ++mfi)
  {
    amrex::Print() << "box " << ibox << ":";
    const FArrayBox& fab = (*a_data)[mfi];
    printFAB(&fab);
  }
}


void 
printMFEdge(const MultiFab* a_data)
{
  amrex::Print() << "data for multifab "  <<endl;
  int ibox = 0;
  for(MFIter mfi(*a_data); mfi.isValid(); ++mfi)
  {
    amrex::Print() << "box " << ibox << ":";
    Box bx = mfi.validbox();
    Box sidebox = adjCellLo(bx, 0, 3);
    sidebox.shift(0,3);
    const FArrayBox& fab = (*a_data)[mfi];
    for(BoxIterator bit(sidebox); bit.ok(); ++bit)
    {
      amrex::Print() << bit()  << "\t" ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
      {
        amrex::Print() << "\t" << fab(bit(),ivar);
      }
      amrex::Print() << "\n";
    }
  }
}

void 
printFABEdge(const FArrayBox* a_data, const Box* bxptr)
{
  amrex::Print() << "data for fab "  <<endl;

  const FArrayBox& fab = (*a_data);
  Box bx;
  if(bxptr != NULL)
  {
    bx = *bxptr;
  }
  else
  {
    bx = fab.box();
  }

  Box sidebox = adjCellLo(bx, 0, 3);
  sidebox.shift(0,3);
  for(BoxIterator bit(sidebox); bit.ok(); ++bit)
  {
    amrex::Print() << bit()  << "\t" ;
    for (int ivar = 0; ivar < fab.nComp(); ivar++)
    {
      amrex::Print() << "\t" << fab(bit(),ivar);
    }
    amrex::Print() << "\n";
  }
}
