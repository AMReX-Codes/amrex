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

#ifndef _EBFLUXFACTORY_H___
#define _EBFLUXFACTORY_H___

#include <cmath>
#include <cstdlib>
#include "AMReX_EBISLayout.H"
#include "AMReX_EBFluxFAB.H"
#include <AMReX_FabFactory.H>

namespace amrex
{

/// Factory class to produce BaseEBFluxFABs.
/**
   Factory class to produce BaseEBFluxFABs.
   This is needed for LevelData construction.
*/
  class EBFluxFactory
    : public FabFactory<EBFluxFAB>
  {
  public:
    /// factory function.
    /**
       Creates a new baseebcellfab object
       and returns a pointer to it.  Responsiblitly
       for calling operator 'delete' on this pointer is passed to the user.
    */
    virtual EBFluxFAB* create(const Box& box, int ncomps, const FabInfo& info, int box_index) const override;

    virtual EBFluxFactory* clone () const override {
        return new EBFluxFactory(*this);
    }

    ///
    /**
       create the factory with an ispace.  calls full define function
    */
    EBFluxFactory(const EBISLayout& a_ebisl);

    ///
    virtual ~EBFluxFactory();

  private:
    EBISLayout m_ebisl;

  };
}

#endif
