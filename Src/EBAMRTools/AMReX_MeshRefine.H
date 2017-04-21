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


#ifndef _MeshRefine_Wrapper_
#define _MeshRefine_Wrapper_

#include <cmath>
#include <cstdlib>
#include <vector>
#include <AMReX_Box.H>
#include <AMReX_IntVectSet.H>
#include <AMReX_BoxArray.H>


/// function to make grids
/**
   This makes a set of grids from a set of tags.   
   This is not as general as Chombo's meshrefine (it has
   none of the stuff necessary for lbase != 0).
   This is useful for  writing tests, however.  
*/
namespace amrex
{
  extern void
  MeshRefine(std::vector<BoxArray>           &   a_grids,
             const std::vector<IntVectSet>   &   a_tags,
             const std::vector<int>          &   a_refRat,
             const int                       &   a_maxLev,
             const int                       &   a_blockingFactor,
             const int                       &   a_properNestingRadius,
             const int                       &   a_maxGridSize,
             const Real                      &   a_gridEfficiency,
             const Box                       &   a_coarsestDomain);
}
#endif
