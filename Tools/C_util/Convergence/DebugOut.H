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

#ifndef _DebugOut_H_
#define _DebugOut_H_

#include "AMReX_Box.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_FabArray.H"
#include "AMReX_MultiFab.H"


///bunch of functions to call from the debugger
using namespace amrex;

///
extern void printFAB(const FArrayBox* a_fabPtr);

///
extern void printDBL(const BoxArray* a_fabPtr);

///
extern void printBFR(const BaseFab<Real>* a_fabPtr);

///
extern void printBFI(const BaseFab<int>* a_fabPtr);

///
/**
   Print inputs to standard out.  For use inside debuggers.
*/
extern void printBA(const BoxArray* a_dblInPtr);

///
/**
   Print inputs to standard out.  For use inside debuggers.
*/
extern void printBox(const Box* a_boxPtr);

///
extern void viewMF(const MultiFab* a_data);

///
extern void maxMinMF(const MultiFab* a_data);

extern void printMF(const MultiFab* a_data);
extern void printMFEdge(const MultiFab* a_data);
extern void printFABEdge(const FArrayBox* a_data, const Box* bx = NULL);

#endif
