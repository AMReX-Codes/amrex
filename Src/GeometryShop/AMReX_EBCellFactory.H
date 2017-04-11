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


#ifndef _EBCELLFACTORY_H_
#define _EBCELLFACTORY_H_

#include <cmath>
#include <cstdlib>
#include <vector>
#include "AMReX_EBISLayout.H"
#include "AMReX_EBCellFAB.H"


/// Factory class to produce EBCellFABs.
/**
   Factory class to produce EBCellFABs.
   This is needed for LevelData construction.
 */
namespace amrex
{
class EBCellFactory
  : public FabFactory<EBCellFAB>
{
public:
  /// factory function.
  /**
      Creates a new baseebcellfab object
      and returns a pointer to it.  Responsiblitly
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual EBCellFAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const
    {
      shared_ptr<FabArray<EBGraph> > graphs = m_ebisl.getAllGraphs();
      shared_ptr<FabArray<EBData> >  data   = m_ebisl.getAllData();
      EBISBox ebisBox((*graphs)[box_index], (*data)[box_index]);
      return new EBCellFAB(ebisBox, box, ncomps);
    }
  ///
  /**
     create the factory with an ispace.  calls full define function
  */
  EBCellFactory(const EBISLayout& a_ebisl)
    {
      m_ebisl = a_ebisl;
    }

  ///
  virtual ~EBCellFactory();

private:
  EBISLayout m_ebisl;

  //disallowed constructors and all that
  EBCellFactory(const EBCellFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }
  EBCellFactory()
  {
    MayDay::Error("disallowed constructor");
  }
  void operator=(const EBCellFactory& a_inputs)
  {
    MayDay::Error("disallowed constructor");
  }

};

#endif
