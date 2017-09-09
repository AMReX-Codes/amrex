
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

#ifndef _EBDebugDump___H_
#define _EBDebugDump___H_

#include "AMReX_EBDebugOut.H"

namespace amrex
{

  void dontCallThisFunctionEver()
  {

    dumpFAB(NULL);

    dumpDBL(NULL);

    dumpBFR(NULL);

    dumpBFI(NULL);
 
    dumpBA(NULL);

    dumpIVS(NULL);

    dumpBox(NULL);

    dumpEBFAB(NULL);

    dumpIVFAB(NULL);
    
    dumpEBFace(NULL);


    dumpVVoFs( NULL);

    dumpVFaces(NULL);

    dumpEBLevel(NULL);

    dumpVoFSten(NULL);

    dumpFaceSten(NULL);
  }
  
}

#endif

