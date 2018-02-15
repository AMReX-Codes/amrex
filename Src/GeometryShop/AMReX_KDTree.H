 /*

LICENSE NOTICE

This source code is a part of the Project X software library.  Project X
solves partial differential equations
in multiple dimensions using an adaptive, discontinuous Galerkin finite
element method, and has been
specifically designed to solve the compressible Euler and Navier-Stokes
equations.

Copyright © 2003-2007 Massachusetts Institute of Technology



This library is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser
General Public License as published by the Free Software Foundation;
either version 2.1 of the License,
or (at your option) any later version.



This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser
General Public License in lgpl.txt for more details.



You should have received a copy of the GNU Lesser General Public License
along with this library; if not, write
to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
Boston, MA 02111-1307 USA.

This software was developed with the assistance of Government
sponsorship related to Government Contract
Number F33615-03-D-3306.  The Government retains certain rights to this
software under that Contract.


For more information about Project X, please contact:


David L. Darmofal
Department of Aeronautics & Astronautics
Massachusetts Institute of Technology
77 Massachusetts Ave, Room 37-401
Cambridge, MA 02139

Phone: (617) 258-0743
FAX:   (617) 258-5143

E-MAIL: darmofal@mit.edu
URL:    http://raphael.mit.edu

*/

/*
This file is part of ``kdtree'', a library for working with kd-trees.
Copyright (C) 2007-2009 John Tsiombikas <nuclear@member.fsf.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/


#ifndef _KDTREE_H_
#define _KDTREE_H_

#ifndef USE_LIST_NODE_ALLOCATOR
#define USE_LIST_NODE_ALLOCATOR 1
#endif

#define RESTRICT __restrict__

#include <string>
using std::string;

#include "AMReX_KDStruct.H"



  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  /*******************KDTree GENERAL INFORMATION*******************/
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  /* PLEASE be familiar with binary trees, KDTrees,
     and Linked Lists before using these functions.
     (See Wikipedia for more info) */

  /* This file defines functions for interacting with KDTree and ListHead */
  /* (which is a LinkedList) data structures defined in KDStruct.h */

  /* This library follows an abstract, object-oriented design methodology. */

  /* KDTree and ListHead allow for each node to reference any type of data, */
  /* as long as every node stores the same type.  That is, each node */
  /* has a field: void *data; */
  /* If you want KDTree/ListHead to free their data when destroyed, */
  /* you MUST call DataDestructor with a function pointer toan appropriate */
  /* free() method for the data. */

  /* KDTrees provide a fast way of searching N k-dimensional data (points).*/
  /* This library allows the user to find the closest data point to a given */
  /* query point in log^k(N) time.  This library also supports searching for */
  /* all data points within a given distance of query point. */

  /* To use KDTree: */
  /* 1) KDCreate to make a KDTree for k-dimensional data */
  /* 1a) KDDataDestructor to set the data destructor, if needed. */
  /* 2) KDInsert the data points (must be done one at a time) */
  /*    NOTE: KDInsert creates a COPY of the input data points */
  /* 2a) KDInsert2d and KDInsert3d are special cases for 2d and 3d data. */
  /* 3a) KDNearestNeighbor to search for the nearest data point to a query point */
  /* 3b) KDNearestRange to search for all data points within some distance of a query */
  /*     Note that NearestRange creates and returns a KDResult list, which must
         be RELEASED when useage is complete! */
  /* 4) Call KDFree */

  /* NOTE: a KDTree with dim<<1 is the same as a Binary Tree. */
  /* Feel free to use KDTree for this purpose, but remember that a delete */
  /* function DOES NOT currently exist! */

  /* FUTURE: KDResult will be removed, and KDTree will use ListHead objects
     to return the results of KDNearestRange. */

  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  /*******************ListHead GENERAL INFORMATION*****************/
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/

  /* ListHead objects "own" LinkedLists (LList).  For maximum abstraction, the */
  /* user CANNOT interact with ListNodes directly.  Instead they must */
  /* use the LList.* functions (described below). ListHead objects provide */
  /* way to iterate over ListNodes and insert/delete nodes. */

  /* As with KDTree, LinkedLists allow for arbitrary data to be stored at */
  /* each ListNode.  Again, LListDataDestructor MUST be called if you want */
  /* LList to free these data when the LList is destroyed */

  /* To use Linked Lists: */
  /* 1) LListCreate to create a ListHead.  The linked list owned by  */
  /*    the ListHead can be ordered by a Real key or unordered.   */
  /* 1a) LListDataDestructor if needed (see KDTree for more details) */
  /* INSERT: */
  /* 2) LListInsert to insert nodes with optional data.  LListInsert will
        insert in order if the associated ListHead is ordered.  Otherwise
        it will insert at the current iterator position. */
  /* 2a) LListInsertFirst inserts a node at the front of the list.  If the
         list is ordered, the node will be inserted in order. */
  /* ITERATORS: */
  /* Since users CANNOT interact with the LinkedList directly, an iterator is*/
  /* provided (as in Java) to allow users to navigate the list. */
  /* 3a) LListRewind sets the iterator to the first node (or to NULL if list is empty) */
  /* 3b) LListNext increments the iterator. */
  /* 3c) LListIterIsNull determines if the iterator is pointing to a LList node or not. */
  /* 3d) LListAtLast determines if the iterator is pointing to the LAST node. */
  /* 3e) LListItem returns the data stored at the node currently pointed to by the list iterator */
  /* 3f) LListKey returns the key at the current node. */
  /* NOTE: the list iterator in ListHead starts off NULL */

  /* 4) Call LListFree to clear the entire LinkedList and destroy the
     ListHead object. */

  /****************************************************************/
  /****************************************************************/
  /****************************************************************/
  /****************************************************************/





/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/************************KDTree FUNCTIONS ****************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/






  /*
typedef struct _kdtree KDTree;
struct _kdtree;
typedef struct _kdres KDResult;
struct _kdres;
typedef struct _linkedlisthead ListHead;
struct _linkedlisthead;
  */

/******************************************************************/
  KDTree *KDCreate(int const k, int *ierr);
/* create a kd-tree for "k"-dimensional data */

/******************************************************************/
int KDFree(KDTree *tree);
/* free the KDTree (calls KDClear, then destroys the head) */

/******************************************************************/
int KDClear(KDTree *tree);
/* remove all the elements from the tree */
/* this does not destroy the head, so the KDTree can be re-filled with */
/* new information */


/******************************************************************/
  void KDSetGlobalData(KDTree *tree, void *data);
/* if called with non-null 2nd argument, tree's globalData pointer
   will be set.  This functionality is usually used in conjunction
   with KDSetGlobalDataDestructor, but that is not required.

   Does not check that tree->globalData is NULL; i.e., this will
   overwrite any existing pointer
 */

/******************************************************************/
  int KDSetGlobalDataDestructor(KDTree *tree, void (*destr)(void*));
/* if called with non-null 2nd argument, then KDTree
   guarantees that destr(tree->globalData) will be called by KDFree().

   tree->globalDestrFlag must be -1 before this call.
   tree->globalDestrFlag will be set to 1.
 */

  void KDGetGlobalData(KDTree *tree, void** data);

/******************************************************************/
int KDSetDataDestructor(KDTree *tree, void (*destr)(void*));
/* if called with non-null 2nd argument, the function provided
 * will be called on data pointers (see KDinsert) when nodes
 * are to be removed from the tree.

 tree->globalDestrFlag must be -1 before this call.
 tree->globalDestrFlag will be set to 0.
 */

/******************************************************************/
  int KDInsert(KDTree *tree, Real const *pos, void *data);
/* insert a node, specifying its position, and optional data */

/******************************************************************/
  int KDTreeStatistics(KDTree const *tree);
  /*
    Prints statistics about the structure of the tree:
    <>max/min/avg depth printed to stdout
    <>depth of each leaf printed to "leafDepths.txt"
   */

  int KDTreePrint(KDTree const *tree);

  int KDNodePrint(KDNode const *node,
                  const string &prefix);

/******************************************************************/
  int KDSearch(KDTree const *tree, Real const *pos, int *foundFlag);
  /* Searches *tree to determine whether a KDNode with coordinate *pos
     has already been inserted.
     *foundflag is set to 1 if *pos is pre-existing and 0 otherwise.
  */

/******************************************************************/
  int KDExhaustiveSearch(KDTree const *tree, Real const *pos, int *numFound);
  /* Searches the entire KDTree exhaustively for KDNodes holding position pos.
     FOR UNIT TESTING ONLY.  DO NOT USE THIS FUNCTION.
  */


/******************************************************************/
int KDNearestNeighbor(KDTree const *tree, Real const* const queryPoint,
                        Real *resultPoint, void **resultData,
                        Real *bestDistSquared, int const approx);


/******************************************************************/
extern int
KDBBoxSearch(KDTree const *tree, Real const * RESTRICT xbb,
               int *pnData, void ***resultData);
/*
  PURPOSE: returns all tree data items within a bounding box

  INPUT:
   tree :      tree on which to perform search
   xbb :       bounding box ( xmin, xmax, ymin, ymax .... )

  OUTPUT:
   nData :      number of tree items inside bounding box
   resultData : array of data items corresponding to bounding box (reallocated)

  RETURN: error code
*/


/******************************************************************/
  KDResult *KDNearestRange(KDTree *kd,
                             Real const *query,
                             Real const range, int *ierr);
/* Find any nearest nodes from the specified point within a range.
 *
 * This function returns a pointer to a result set, which can be manipulated
 * by the KDResult_* functions.
 * The returned pointer can be null as an indication of an error. Otherwise
 * a valid result set is always returned which may contain 0 or more elements.
 * The result set must be deallocated with KDResultFree, after use.
 */


/******************************************************************/
void KDResultFree(KDResult *set);
/* frees a result set returned by KDnearest_range() */

/******************************************************************/
int KDResultSize(KDResult const *set);
/* returns the size of the result set (in elements) */

/******************************************************************/
void KDResultRewind(KDResult *set);
/* rewinds the result set iterator */

/******************************************************************/
int KDResultEnd(KDResult const *set);
/* returns non-zero if the set iterator reached the end after the last element */

/******************************************************************/
int KDResultNext(KDResult *set);
/* advances the result set iterator, returns non-zero on success, zero if
 * there are no more elements in the result set.
 */

/******************************************************************/
void *KDResultItem(KDResult const *set, Real *pos);
/* returns the data pointer (can be null) of the current result set item
 * and optionally sets its position to the pointers(s) if not null.
 */


/******************************************************************/
void *KDResultItemData(KDResult const *set);
/* equivalent to KDResultItem(set, 0) */




/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/************************LinkedList FUNCTIONS ************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/




/******************************************************************/
  int LListInsert(ListHead *head, void *data,
                    Real const key, int const useListIter);
/*
Inserts a node with data and key.  If head->ordered is 0, then
the node will be inserted in order by key (increasing order).

If useListIter is 0, behavior is the same as LListInsertFirst.

If useListIter is nonzero, unordered insertion places the node AFTER
the node that list->listIter points to.  Returns -1 if listIter == NULL.
Ordered insertion attempts to insert AFTER listIter; if it can't (b/c such
an insertion is out of order OR listIter == NULL), then it will search
and insert at the correct location.
*/

/******************************************************************/
  int LListInsertFirst(ListHead *head, void *data,
                         Real const key);
/*
unordered: Inserts a node with data and key at the head of head->llist.
ordered: Begins at the head and searches forward until the correct position
  for this key is found.
 */

/******************************************************************/
  ListHead *LListCreate(int const ordered, int *ierr);
/* create a linked list object.  if ordered == 1, the list will ensure
   that LListNodes are inserted in order.  Else no ordering is guaranteed.*/

/******************************************************************/
  void LListFree(ListHead *head);
  /*Frees the linked list associated with this ListHead. */

  /*****************************************************************/
  /*  int LListIndexOf(ListHead *head, void (*isEqual)(void*,void*), void* data1, int goThereFlag, int indexOf);*/
  /* Using the comparator isEqual(void* data1, void* data2), where *data1
     is input and *data2 corresponds to the data in some LListNode,
     IndexOf returns the first occurence of data1 in the list.

     If there are no occurences, IndexOf = -1.

     If goThereFlag = 1, the list iterator will point to the node at position
     indexOf when LListIndexOf completes.  If indexOf = -1, the list iterator will
     be returned to the first node.
     If goThereFlag = 0, the list iterator will not be changed.
  */

/******************************************************************/
  int LListSize(ListHead const *head);
  /*Returns the size of the LinkedList */

/******************************************************************/
  void LListRewind(ListHead *head);
  /*Rewinds the iterator one step.*/

/******************************************************************/
  int LListIterIsNull(ListHead const *head);
  /*Returns non-zero if the LinkedList iterator has reached the end of the list.*/

/******************************************************************/
  int LListAtLast(ListHead const *head);
  /*Returns non-zero if the LinkedList iterator is at the last node of the list.*/


/******************************************************************/
  int LListNext(ListHead *head);
  /*
    Advances the LinkedList iterator.  Returns non-zero on success; returns zero if there
    are no more elements (i.e. if LListEnd would return nonzero).
  */

/******************************************************************/
  void *LListItem(ListHead const *head);
  /*
    Returns the data pointer of the LListNode pointed to by head->listIter.
    Returns NULL if listIter points to nothing.
  */

/******************************************************************/
  Real LListKey(ListHead const *head);
  /* Returns the value of the key for the node pointed to by listIter. */

/******************************************************************/
void LListDataDestructor(ListHead *head, void (*destr)(void*));
/* if called with non-null 2nd argument, the function provided
 * will be called on data pointers (see KDInsert) when nodes
 * are to be removed from the tree.
 */

#ifdef USE_LIST_NODE_ALLOCATOR
  void LListFinalize( void );
  void KDTreeFinalize( void );
#endif

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/************************Stack FUNCTIONS *****************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

  void *StackPop(ListHead *head, int *ierr);
  /* Pop the first node off of the stack.
     Throws an error if the stack is empty. */

  /*******************************************************************/
  int StackPush(ListHead *head, void *data);
  /* Push a node onto the stack.
     The node key is NOT set, so it could take any value. */

  /*******************************************************************/
  void *StackPopWithKey(ListHead *head, Real *key, int *ierr);
  /* Pop the first node off of the stack, also returning the key.
     Throws an error if the stack is empty.

     WARNING: it is NOT SAFE to call this function with a stack
     built with StackPush b/c the key is NOT set! */

  /*******************************************************************/
  int StackPushWithKey(ListHead *head, void *data, Real const key);
  /* Push a node onto the stack with a specified key. */

  /*******************************************************************/
  ListHead *StackCreate(int *ierr, int const prealloc);
  /* Create a Stack object. */

  /*******************************************************************/
  void StackFree(ListHead *head);
  /* Free the Stack object, even if the stack is non-empty.
     This will free all stack members with a warning. */



#endif  /* _KDTREE_H_ */
