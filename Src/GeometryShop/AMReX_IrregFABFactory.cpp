
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
#include "AMReX_IrregFABFactory.H"
#include "AMReX_FabArray.H"

namespace amrex
{

  IrregFAB*
  IrregFABFactory::
  create (const Box& box, int ncomps, const FabInfo& info, int box_index) const
  {
    EBGraph& graph       = (*m_graphs)[box_index];
    IntVectSet ivs;
    if(m_useSets)
    {
      ivs   = (*m_sets)[box_index];
    }
    else
    {
      Box restbox = box;
      restbox &= graph.getDomain();
      ivs = graph.getIrregCells(restbox);
    }

    return new  IrregFAB(ivs, graph, ncomps);
  }
}
