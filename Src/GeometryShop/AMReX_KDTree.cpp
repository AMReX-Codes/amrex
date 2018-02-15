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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits>

#include "PXStuff.H"

/* Structures */
#include "AMReX_KDStruct.H"
#include "AMReX_KDTree.H"
#include "AMReX_parstream.H"


#include "AMReX_REAL.H"


#define SQ(x)                   ((x) * (x))


static int KDClearRec(KDNode *node, void (*destr)(void*), char globalDestrFlag);
static int KDClearIter(KDNode *node, void (*destr)(void*), char globalDestrFlag);
static int ResultListInsert(ResultNode *list, KDNode *item, Real dist_sq);
static void ClearResults(KDResult *set);

static void ClearLList(ListHead *head);

#ifdef USE_LIST_NODE_ALLOCATOR
static int AllocResultNode(ResultNode **);
static int AllocLListNode(LListNode **);
static void FreeResultNode(ResultNode*);
static void FreeLListNode(LListNode*);
#else
static int AllocResultNode((ResultNode **)node)
{
  //return Allocate(1,sizeof(ResultNode),(void **)&(node) );
  *node = (ResultNode*)malloc(sizeof(ResultNode));
  return PX_NO_ERROR;
}
static int AllocLListNode((LListNode **)node)
{
  //return Allocate(1,sizeof(LListNode),(void **)&(*node) );
  *node = (LListNode *)malloc(sizeof(LListNode));
  return PX_NO_ERROR;
}
#define FreeResultNode(n)               free(n)

static void FreeLListNode(LListNode *node)
{
  free( node );
}
#endif


/******************************************************************/
KDTree *KDCreate(int const k, int *ierr)
{
  KDTree *tree; /* the tree to be created */

  if (k > 255)
  {
    *ierr = PX_CODE_FLOW_ERROR;
    return NULL;
  }
  //dimension > 255 not supported b/c nodes store dimension as a uchar

  //*ierr = Allocate(1,sizeof(KDTree),(void **)&(tree) ); /*allocate space for the struc*/
  //if (*ierr !=PX_NO_ERROR) return tree;
  tree = (KDTree *)malloc(sizeof(KDTree));

  (*ierr) = PX_NO_ERROR;
  /* initialize */
  tree->dim = k;
  tree->root = NULL;
  tree->destr = NULL;
  tree->globalData = NULL;
  tree->globalDestrFlag = -1;
  return tree;
}

/******************************************************************/
int KDFree(KDTree *tree)
{
  // int ierr;
  if (tree) 
{
    /* Make sure all nodes are freed */
    PXErrorReturn( KDClear(tree) );
    /* free KDTree structure */
    if (tree->globalDestrFlag == 1 && tree->globalData != NULL)
      tree->destr(tree->globalData);

    free( tree );
  }
  return PX_NO_ERROR;
}

/******************************************************************/
static int KDClearRec(KDNode *node, void (*destr)(void*), char globalDestrFlag)
{
  if (node == NULL) return PX_NO_ERROR; /*base case: stop if node is NULL */

  /*recursive calls on subtrees */
  KDClearRec(node->left, destr, globalDestrFlag);
  KDClearRec(node->right, destr, globalDestrFlag);

  /* call data destroy function pointer if needed */
  if (globalDestrFlag == 0 && destr != NULL && node->data != NULL) 
{
    (*destr)(node->data);
  }
  /* free the node struct */
  free( node );

  return PX_NO_ERROR;
}

/******************************************************************/
static int KDClearIter(KDNode *node, void (*destr)(void*), char globalDestrFlag)
{
  /* Clears all nodes of a KDTree below *node, including *node. */

  int ierr; /*error code*/
  KDNode *curNode; /* current node for deletion */
  ListHead *stack; /* manual program stack: used in place of the system
                      stack to avoid recursion */

  if (NULL == node)
    return PX_NO_ERROR; /* Tree is empty, nothing to do */

  stack = StackCreate(&ierr, 0);
  PXErrorReturn(ierr);

  PXErrorReturn( StackPush(stack, (void *) node) ); /* Push root node onto the stack */

  while (LListSize(stack) != 0)
  {
    curNode = (KDNode *) StackPop(stack, &ierr);
    PXErrorReturn(ierr);

    /* Push curNode's children */
    if (curNode->left != NULL)
      PXErrorReturn( StackPush(stack, (void *) curNode->left) );

    if (curNode->right != NULL)
      PXErrorReturn( StackPush(stack, (void *) curNode->right) );

    /* free curNode */
    if (globalDestrFlag == 0 && destr != NULL && curNode->data != NULL)
    {
      (*destr)(curNode->data);
    }
    /* free the node struct */
    free( curNode );
  }

  StackFree(stack);

  return PX_NO_ERROR;
}

/******************************************************************/
int KDClear(KDTree *tree)
{
  // int ierr;
  /* Free the subtrees*/
  PXErrorReturn( KDClearIter(tree->root, tree->destr, tree->globalDestrFlag) );
  /* Set root node to null */
  tree->root = NULL;
  return PX_NO_ERROR;
}

/******************************************************************/
void KDSetGlobalData(KDTree *tree, void *data)
{
  tree->globalData = data;
}

void KDGetGlobalData(KDTree *tree, void** data)
{
  *data = tree->globalData;
}

/******************************************************************/
int KDSetGlobalDataDestructor(KDTree *tree, void (*destr)(void*))
{
  /* Sets the destructor function for the data held by each KDNode.
     This is NULL by default and can be left that way if the KDNodes
     do not hold data. */
  if (tree->globalDestrFlag != -1) return PX_CODE_FLOW_ERROR;
  if (destr == NULL) return PX_NO_ERROR;
  tree->destr = destr;
  tree->globalDestrFlag = 1;

  return PX_NO_ERROR;
}

/******************************************************************/
int KDSetDataDestructor(KDTree *tree, void (*destr)(void*))
{
  /* Sets the destructor function for the data held by each KDNode.
     This is NULL by default and can be left that way if the KDNodes
     do not hold data. */
  if (tree->globalDestrFlag != -1) return PX_CODE_FLOW_ERROR;
  if (destr == NULL) return PX_NO_ERROR;
  tree->destr = destr;
  tree->globalDestrFlag = 0;

  return PX_NO_ERROR;
}


/******************************************************************/
static inline int InsertIter(KDTree *tree, Real const *pos, void *data)
{
  int curDir=-1;              // current normal direction for inserted node's dividing
  // int newDir;                 // new normal direction for the inserted node's dividing
  int d, dim;                 // #of space dimensions for this tree
  KDNode *curNode = NULL;     // current node for traversal
  KDNode *newNode;            // new node to be allocated
  KDNode **nextNode;          // the nextNode to be searched
  Real *RESTRICT curNodePos = NULL;
  enum PXE_Boolean equalFlag = PXE_False; // flag to check nodes are equal

  dim = tree->dim;

  nextNode = &(tree->root);
  while (NULL != (*nextNode))
  {

    curNode    = *nextNode;
    curDir     = curNode->dir;
    curNodePos = curNode->pos;

    if (pos[curDir] < curNodePos[curDir]) 
{
      nextNode = &(curNode->left);
    }
    else if (pos[curDir] > curNodePos[curDir])
    {
      nextNode = &(curNode->right);
    }
    else
    {
      /* pos[node->dir] == node->pos[node->dir], so we must check for the
         possibility that these two points are identical. */
      equalFlag = PXE_True;
      for (d = 0; d<dim; d++)
      {
        if (pos[d] != curNodePos[d])
        {
          /* break if an inequality is detected */
          equalFlag = PXE_False;
          break;
        }
      }
      if ( equalFlag ) {
        // printf("Warning node with this coordinate already exists in tree, not adding it\n");
        /* Node is already in the tree. Do not reinsert it.*/
        if (tree->globalDestrFlag == 0 && tree->destr != NULL && data != NULL)
          tree->destr(data);

        return PX_NO_ERROR;
      }
      else 
      {
        /* if the points are not the same, insert left (convention) */
        nextNode = &(curNode->left);
      }
    } /* close else statement for possible node equality */
  } /* close loop to search for inserted node's position */

  /* If (*nextNode) is NULL, that means we have reached a null pointer
     from curNode->left OR curNode->right (i.e. a leaf node).  Thus
     the search is complete.

     KEY IDEA: nextNode points to the ADDRESS of the location where the
     new node should be inserted!

     Now, insert the input *pos and *data in the place pointed to by
     nextNode. */

  /* curNode is now the node to be inserted */
  //PXErrorReturn( Allocate( sizeof(KDNode), 1, (void **)&(newNode) ) );
  newNode = (KDNode*)malloc(sizeof(KDNode));

  /*copy the point from the input--do not allow pointer chains */
  newNode->data = data;
  newNode->dir  = (curDir + 1) % dim;
  newNode->left = newNode->right = NULL;
  memcpy(newNode->pos, pos, dim*sizeof(Real));
  if ( curNode != NULL ) 
{
    memcpy(newNode->xmin, curNode->xmin, dim*sizeof(Real));
    memcpy(newNode->xmax, curNode->xmax, dim*sizeof(Real));

    if ( pos[curDir] > curNodePos[curDir] ) newNode->xmin[curDir] = curNodePos[curDir];
    else                                    newNode->xmax[curDir] = curNodePos[curDir];
  }
  else 
  {
    for ( d = 0; d < dim; d++) newNode->xmin[d] = -DBL_MAX;
    for ( d = 0; d < dim; d++) newNode->xmax[d] =  DBL_MAX;
  }


  *nextNode = newNode; /* put newNode into its position as determined by
                          the search. */
  return PX_NO_ERROR;

}

/******************************************************************/
int KDInsert(KDTree *tree, Real const *pos, void *data)
{
  /*Wrapper for InsertIter */
  return InsertIter(tree, pos, data);
}


/******************************************************************/
int KDTreeStatistics(KDTree const *tree)
{
  int ierr; /*error code*/
  int numNodes = 0, numLeaves = 0, leafCtr = 0, numSingletons = 0;
  int maxStackSize = 0;

  KDNode *curNode; /* current node for deletion */
  ListHead *stack; /* manual program stack: used in place of the system
                      stack to avoid recursion */

  if (NULL == tree || NULL == tree->root)
  {
    printf("KDTree is empty, nothing to do.\n");
    return PX_NO_ERROR; /* Tree is empty, nothing to do */
  }

  stack = StackCreate(&ierr, 0);
  PXErrorReturn(ierr);

  PXErrorReturn( StackPush(stack, (void *) tree->root) ); /* Push root node onto the stack */

  /* exhaustively visit every node; in-order traversal */
  /* determine number of leaves & number of nodes */
  while (LListSize(stack) != 0)
  {
    if (LListSize(stack) > maxStackSize)
      maxStackSize = LListSize(stack);

    curNode = (KDNode *) StackPop(stack, &ierr);
    PXErrorReturn(ierr);
    numNodes += 1;

    if (curNode->left == NULL && curNode->right == NULL)
    {
      numLeaves += 1;
    }
    else
    {
      if ( (curNode->right != NULL && curNode->left == NULL) ||
           (curNode->left != NULL && curNode->right == NULL))
        numSingletons += 1;

      /* Push curNode's children */
      if (curNode->right != NULL)
        PXErrorReturn( StackPush(stack, (void *) curNode->right) );

      if (curNode->left != NULL)
        PXErrorReturn( StackPush(stack, (void *) curNode->left) );
    }

  }

  Real depth = 0.0;
  int *leafDepths;
  //PXErrorReturn( Allocate( numLeaves, sizeof(int), (void **)&(leafDepths) ) );
  leafDepths = (int*)malloc(numLeaves*sizeof(int));

  PXErrorReturn( StackPushWithKey(stack, (void *) tree->root, depth) ); /* Push root node onto the stack */

  /* exhaustively visit every node; in-order traversal */
  /* gather leaf depth info */
  while (LListSize(stack) != 0)
  {
    curNode = (KDNode *) StackPopWithKey(stack, &depth, &ierr);
    PXErrorReturn(ierr);

    if (curNode->left == NULL && curNode->right == NULL)
    {
      leafDepths[leafCtr] = (int) depth;
      leafCtr += 1;
    }
    else
    {
      /* Push curNode's children */
      if (curNode->right != NULL)
        PXErrorReturn( StackPushWithKey(stack, (void *) curNode->right, depth + 1.0) );

      if (curNode->left != NULL)
        PXErrorReturn( StackPushWithKey(stack, (void *) curNode->left, depth + 1.0) );
    }

  }
  StackFree(stack);

  /* compute statistics and such */
  int minDepth = leafDepths[0];
  int maxDepth = leafDepths[0];
  int sumDepth = 0;
  int i, curDepth;

//  FILE *leafDepthFile;
//  if ( (leafDepthFile = fopen("leafDepths.txt","w"))==NULL){
//    printf("failure opening file leafDepths.txt\n");
//    PXErrorReturn(PX_READWRITE_ERROR);
//  }

  for (i=0; i<numLeaves; i++)
  {
    curDepth = leafDepths[i];
    sumDepth += curDepth;
    if (minDepth > curDepth)
      minDepth = curDepth;
    else if (maxDepth < curDepth)
      maxDepth = curDepth;

//    fprintf(leafDepthFile, "%d\n",curDepth);
  }

  char diag[10000];
  sprintf(diag,"\n numNodes = %d, numSingletons = %d, numLeaves = %d, maxStackSize = %d\n", numNodes, numSingletons, numLeaves,maxStackSize);
  pout() << diag;
  sprintf(diag,"maxDepth = %d, minDepth = %d, avg = %.6E\n",maxDepth, minDepth, ((Real) sumDepth)/((Real) numLeaves));
  pout() << diag;

  free(leafDepths);
//  fclose(leafDepthFile);

  return PX_NO_ERROR;
}



/******************************************************************/
static inline int isEqual(int const dim, Real const *pos1, Real const *pos2)
{
  /* Returns 1 if arrays pos1 and pos2 are element-wise equal.
     Returns 0 otherwise. */
  int dimCtr; /*iterate over dimensions */
  for (dimCtr = 0; dimCtr<dim; dimCtr++)
  {
    if (pos1[dimCtr] != pos2[dimCtr])
    {
      /* stop if an inequality is detected */
      return 0;
    }
  }
  return 1;
}


/******************************************************************/
int KDSearch(KDTree const *tree, Real const *pos, int *foundFlag)
{
  /* This search is conducted pseudo-recursively (i.e. using a while loop,
     which is guaranteed to terminate b/c the KDTree has finite size and
     each loop step goes down one level.

     The code is slightly modified from KDNNSBase, which is functionality used by
     KDNearestNeighbor, and is not available to users.
  */
  int dim = tree->dim; /*number of dimensions of points in the tree*/
  int curDir = 0; /*normal direction of the current hyperplane we're splitting on, e.g. x, y, z in 3d*/
  KDNode *tempNode = tree->root; /*temp variable to hold tree->root */
  Real *RESTRICT tempPos = tempNode->pos;

  if (!tempNode)
  { /*make sure tree is not empty*/
    (*foundFlag) = 0;   /* Tree is empty; clearly we didn't find it.*/
    return PX_NO_ERROR;
  }

  /* base case, only 1 node */
  if (tempNode->left == NULL && tempNode->right == NULL)
  { /*if tree has only 1 node*/
    /* printf("Tree has only 1 node.\n");*/
    (*foundFlag) = isEqual(dim, pos, tempPos);
    return PX_NO_ERROR;
  }

  (*foundFlag) = 0;
  for (;;)
  {
    tempPos = tempNode->pos;
    (*foundFlag) = isEqual(dim, pos, tempPos);
    if ((*foundFlag) ==  1)
      return PX_NO_ERROR;

    curDir = tempNode->dir;

    if (pos[curDir] < tempPos[curDir])
    {
      /*go left*/
      tempNode = tempNode->left;
    }
    else
    {
      /*go right*/
      tempNode = tempNode->right;
    }

    if (tempNode == NULL)
      break;
  }
  /* No node matches *pos, return successfully */
  return PX_NO_ERROR;

}

/******************************************************************/
static int KDExhaustiveSearchRecursive(KDNode *node, Real *pos, int dim, int *numFound)
{
  int ierr = PX_NO_ERROR;

  /* error code*/
  Real *nodePos = node->pos;

  if (isEqual(dim, pos, nodePos) == 1)
    (*numFound)++;

  if (node->left != NULL)
    ierr = KDExhaustiveSearchRecursive(node->left, pos, dim, numFound);

  if (ierr != PX_NO_ERROR)
  {
    (*numFound) = 0;
    PXErrorReturn(ierr);
  }

  if (node->right != NULL)
    ierr = KDExhaustiveSearchRecursive(node->right, pos, dim, numFound);

  if (ierr != PX_NO_ERROR)
    (*numFound) = 0;

  PXErrorReturn(ierr);

  return ierr;
}



/******************************************************************/
static int KDExhaustiveSearchIter(KDNode const *node, Real const *pos, int const dim, int *numFound)
{
  int ierr; /*error code*/
  Real *nodePos;
  KDNode *curNode; /* current node for deletion */
  ListHead *stack; /* manual program stack: used in place of the system
                      stack to avoid recursion */

  if (NULL == node)
  {
    (*numFound) = 0;
    return PX_NO_ERROR; /* Tree is empty, nothing to do */
  }

  stack = StackCreate(&ierr, 0);
  PXErrorReturn(ierr);

  PXErrorReturn( StackPush(stack, (void *) node) ); /* Push root node onto the stack */

  while (LListSize(stack) != 0)
  {
    curNode = (KDNode *) StackPop(stack, &ierr);
    PXErrorReturn(ierr);

    /* Push curNode's children */
    if (curNode->left != NULL)
      PXErrorReturn( StackPush(stack, (void *) curNode->left) );

    if (curNode->right != NULL)
      PXErrorReturn( StackPush(stack, (void *) curNode->right) );

    /* check curNode for equality & increment if needed */
    nodePos = curNode->pos;
    if (isEqual(dim, pos, nodePos) == 1)
      (*numFound)++;

  }

  StackFree(stack);

  return PX_NO_ERROR;
}


/******************************************************************/
int KDExhaustiveSearch(KDTree const *tree, Real const *pos, int *numFound)
{
  /* Searches the entire KDTree exhaustively for KDNodes holding position pos.
     FOR UNIT TESTING ONLY.  DO NOT USE THIS FUNCTION.

     The search is preorder, i.e. root, left, right).
  */
  KDNode *tempNode = tree->root; /*temp variable to hold tree->root */
  (*numFound) = 0;

  if (tempNode == NULL)
    return PX_NO_ERROR; /* tree is empty, we found nothing */

  if (tempNode->left == NULL && tempNode->right == NULL)
  { /*if tree has only 1 node*/
    /* printf("Tree has only 1 node.\n");*/
    if (isEqual(tree->dim, pos, tempNode->pos) == 1)
      (*numFound)++;
    return PX_NO_ERROR;
  }
  return KDExhaustiveSearchIter(tempNode, pos, tree->dim, numFound);

}
/* [A-Z][A-Za-z]*(.*) -> \& */


/******************************************************************/
static inline Real distanceSquared(int const dim, Real const *pos1, Real const *pos2)
{
  int i; /*iterate over dimensions*/
  Real distSq = 0; /*variable to accumulate the distance */

  /*calculate square distance between pos1 and pos2 */
  for (i=0; i<dim; i++) 
{
    distSq += SQ(pos1[i] - pos2[i]);
  }
  return distSq;
}

/******************************************************************/
static int KDNNSRecursiveSearch
(KDNode *node, Real const * RESTRICT queryPoint, Real **resultPoint, void **resultData,
 Real *bestDistSquared, int dim)
{
  KDNode *nearNode, *farNode; /* the left and right children of this node.
                                 If the query point falls strictly to the left of this node's
                                 dividing line, then near is the left child. */
  Real *nodePos;
  Real *xmin, *xmax;       /* bounds of box region corresponding to node */
  int d, curDir; /*normal direction of the current hyperplane we're splitting on, e.g. x, y, z in 3d*/
  // int ierr; /*error flag*/
  Real curDistSquared;   /*square dist from the point owned by this node to query point   */

  if (node == NULL)
    return PX_NO_ERROR;

  /* Check whether the box corresponding to this node is possible */
  xmin = node->xmin;
  xmax = node->xmax;
  curDistSquared = 0.0;
  for ( d = 0; d < dim; d++) 
{
    if      ( queryPoint[d] < xmin[d] ) curDistSquared += ( queryPoint[d] - xmin[d] )*( queryPoint[d] - xmin[d] );
    else if ( queryPoint[d] > xmax[d] ) curDistSquared += ( queryPoint[d] - xmax[d] )*( queryPoint[d] - xmax[d] );
  }
  if ( curDistSquared > *bestDistSquared ) return PX_NO_ERROR;

  /* Check current code */
  curDir  = node->dir;
  nodePos = node->pos;
  curDistSquared = distanceSquared(dim, nodePos, queryPoint);
  if (curDistSquared < *bestDistSquared)
  {
    *bestDistSquared = curDistSquared;
    *resultPoint = nodePos;
    *resultData = node->data;
  }

  if (queryPoint[curDir] < nodePos[curDir])
  {
    nearNode = node->left;
    farNode = node->right;
  }
  else
  {
    nearNode = node->right;
    farNode = node->left;
  }

  /* Search both left and right, but search closest */
  PXErrorReturn( KDNNSRecursiveSearch(nearNode, queryPoint, resultPoint, resultData, bestDistSquared, dim) );
  PXErrorReturn( KDNNSRecursiveSearch( farNode, queryPoint, resultPoint, resultData, bestDistSquared, dim) );

  return PX_NO_ERROR;
}




/******************************************************************/
/* Nearest Neighbor Search
   worst: O(k*N^(1-1/k)) where N = # pts and k = dimension
   expected: O(log2(N))
*/

static int KDNNSBase(KDTree const *tree, Real const *queryPoint,
                     Real **resultPoint, void **resultData, Real *bestDistSquared,
                     int const approx)
{
  CH_TIME("KDNNSBase");

  int curDir = 0; /*normal direction of the current hyperplane we're splitting on, e.g. x, y, z in 3d*/
  Real curDistSquared = 0; /*square dist from the point owned by this node to query point */
  KDNode *tempNode = tree->root; /*temp variable to hold tree->root */
  Real *RESTRICT nodePos;

  if (!tempNode) /*make sure tree is not empty*/
    return PX_BAD_INPUT;
  /* error, tree empty*/

  if (tempNode->left == NULL && tempNode->right == NULL)
  { /*if tree has only 1 node*/
    //printf("Tree has only 1 node.\n");
    *resultPoint = tempNode->pos;
    *bestDistSquared = distanceSquared(tree->dim, *resultPoint, queryPoint);
    *resultData = tempNode->data;
    return PX_NO_ERROR;
  }

  /* Pretending that our task is to insert queryPoint into the tree,
   * find where it would go.
   * This information is used to get a reasonable initial guess of
   * bestDistSquared, which cheapens the recursive search.
   * The value resulting here is used if an approximate solution is requested.*/
  for (;;)
  {
    curDir = tempNode->dir;
    nodePos = tempNode->pos;
    curDistSquared = distanceSquared(tree->dim, nodePos, queryPoint);
    if (curDistSquared < *bestDistSquared)
    {
      *bestDistSquared = curDistSquared;
      *resultPoint = nodePos;
      *resultData = tempNode->data;
      /*printf("Result point updated: (%.3f, %.3f, %.3f)\n",
        (*resultPoint)[0],(*resultPoint)[1],(*resultPoint)[2]);
      */
    }

    if (queryPoint[curDir] < nodePos[curDir])
    {
      /*go left*/
      tempNode = tempNode->left;
    }
    else
    {
      /*go right*/
      tempNode = tempNode->right;
    }

    if (tempNode == NULL)
      break;
  }

  if (approx == 1)
    return PX_NO_ERROR;

  {
    /* Now we must recursively search the tree for points. */
    CH_TIME("KDNNSBase_KDNNSRecursiveSearch");
    return KDNNSRecursiveSearch(tree->root, queryPoint, resultPoint, resultData,
                                bestDistSquared, tree->dim);
  }
}

/******************************************************************/
int KDNearestNeighbor(KDTree const *tree, Real const *queryPoint,
                      Real *resultPoint, void **resultData, Real *pbestDistSquared,
                      int const approx)
{
  CH_TIME("KDNearestNeighbor");
  
  // int ierr; /*error code*/
  // int d;
  Real bestDistSquared;             // best distance squared
  Real *resultPointInternal = NULL; /* internal representation of resultPoint
                                     * b/c we do not want the user to have any
                                     * kind of access to the data in the tree nodes*/


  /* Sanity-check the inputs */
  bestDistSquared = DBL_MAX;
  if ( pbestDistSquared != NULL ) 
{
    if (( *pbestDistSquared >= 0.0 ) && ( resultPoint != NULL ) && (resultData != NULL)) 
{
      /* At lest all inputs are consistent to have previous guess */
      bestDistSquared = distanceSquared(tree->dim, queryPoint, resultPoint);

      /* Check bestDistSquared matches input */
      if (Abs(bestDistSquared - (*pbestDistSquared)) > 1e-15)
        bestDistSquared = DBL_MAX;
    }
    (*pbestDistSquared) = bestDistSquared;
  }


  /*Pass execution on to the base search. */
  PXErrorReturn( KDNNSBase(tree, queryPoint, &resultPointInternal, resultData,
                           &bestDistSquared, approx) );


  if ((resultPoint != NULL) && (resultPointInternal != NULL) ) 
{
    /* update the output resultPoint. */
    memcpy(resultPoint, resultPointInternal, tree->dim*sizeof(Real));
  }

  if ( pbestDistSquared != NULL ) 
{
    (*pbestDistSquared) = bestDistSquared;
  }

  return PX_NO_ERROR;
}


/******************************************************************/
static int
KDListNodesRecursive( ListHead *list, KDNode *node )
{
  if ( node == NULL ) return PX_NO_ERROR;

  /* insert node */
  PXErrorReturn( LListInsert( list, node, 0, 0) );


  /* Insert children */
  PXErrorReturn( KDListNodesRecursive( list, node->left  ) );
  PXErrorReturn( KDListNodesRecursive( list, node->right ) );

  return PX_NO_ERROR;
}


/******************************************************************/
static int
KDBBoxSearchRecursive(ListHead *list, KDNode *node, int const dim,
                      Real const * RESTRICT xbb)
{
  Real *nodePos;           // coordinate of node
  Real *xmin, *xmax;       // bounds of box region corresponding to node
  int d;                      // dimension
  enum PXE_Boolean flag;

  if (node == NULL)
    return PX_NO_ERROR;

  xmin = node->xmin;
  xmax = node->xmax;

  /* Check whether the box corresponding to this node is possible */
  flag  = PXE_False;
  for ( d = 0; d < dim; d++) 
{
    if ( xbb[2*d+1] < xmin[d] ) 
    {
      flag = PXE_True; 
      break; 
    }
    if ( xbb[2*d  ] > xmax[d] ) 
    { 
      flag = PXE_True; 
      break; 
    }
  }
  if ( flag ) 
  {
    /* bounding boxes don't overlap */
    return PX_NO_ERROR;
  }

  /* Check whether the box corresponding to this node contained within bounding box */
  flag = PXE_True;
  for ( d = 0; d < dim; d++) 
  {
    if      ( xbb[2*d  ] >= xmin[d] ) 
    { 
      flag = PXE_False; 
      break; 
    }
    else if ( xbb[2*d+1] <= xmax[d] ) 
    { 
      flag = PXE_False; 
      break; 
    }
  }
  if ( flag ) 
  {
    /* Return this node and all children */
    PXErrorReturn( KDListNodesRecursive( list, node ) );
    return PX_NO_ERROR;
  }

  /* Check whether current node is in Bounding Box */
  nodePos = node->pos;
  flag = PXE_True;
  for ( d = 0; d < dim; d++) 
  {
    if      ( xbb[2*d  ] > nodePos[d] ) 
    { 
      flag = PXE_False; 
      break; 
    }
    else if ( xbb[2*d+1] < nodePos[d] ) 
    { 
      flag = PXE_False; 
      break; 
    }
  }
  if ( flag ) 
  {
    /* Add this node */
    PXErrorReturn( LListInsert( list, node, 0, 0) );
  }

  /* Check children */
  PXErrorReturn( KDBBoxSearchRecursive( list, node->left , dim, xbb ) );
  PXErrorReturn( KDBBoxSearchRecursive( list, node->right, dim, xbb ) );

  return PX_NO_ERROR;
}





/******************************************************************/
/* Bounding box search
   worst: O(k*N^(1-1/k)) where N = # pts and k = dimension
   expected: O(log2(N))
*/

int KDBBoxSearch(KDTree const *tree, Real const * RESTRICT xbb,
                 int *pnData, void ***resultData)
{
  int ierr;
  int i, nData;
  KDNode *node;
  // Real *temp;
  ListHead *list;


  /* Create linked list to store tree nodes */
  list = LListCreate(0, &ierr); PXErrorReturn(ierr);


  PXErrorReturn( KDBBoxSearchRecursive( list, tree->root, tree->dim, xbb ) );

  nData = LListSize( list );
  //PXErrorReturn( ReAllocate( nData, sizeof(void *), (void **) resultData) );
  *resultData = (void **)realloc(*resultData,nData*sizeof(void*));

  /* Set output */
  (*pnData) = nData;
  for ( i = 0; i < nData; i++) 
  {
    node = (KDNode *)StackPop( list, &ierr ); PXErrorReturn(ierr);
    (*resultData)[i] = node->data;
  }

  LListFree(list);

  return PX_NO_ERROR;
}









/******************************************************************/
static int FindNearest(KDNode *node, Real const *query,
                       Real const range, ResultNode *list,
                       int const ordered, int const dim, int *ierr)
{
  Real distSq;/*square dist from the point owned by this node to query point */
  Real dx; /* distance to the dividing hyperplane owned by this node */
  Real *nodePos = node->pos;
  int ret; /*tracks how many points are added to the Result list by each recursive search */
  int addedRes = 0; /*The total number of points found to be in range */

  if (node == NULL)
  {
    (*ierr) = PX_NO_ERROR;
    return 0; /*base case*/
  }

  distSq = distanceSquared(dim, nodePos, query);
  /*If the query is already known to be in range, then add it to the KDResult list*/
  if (distSq <= SQ(range)) 
  {
    (*ierr) = ResultListInsert(list, node, ordered ? distSq : -1.0);
    if (*ierr != PX_NO_ERROR)
      return -1;

    addedRes = 1;
  }

  dx = query[node->dir] - nodePos[node->dir];

  /*search the "near" child of this TreeNode */
  ret = FindNearest(dx <= 0.0 ? node->left : node->right, query, range, list, ordered, dim, ierr);
  if (ret >= 0 && Abs(dx) < range) 
  {
    /*If we cannot rule out the "far" child, search it recursively as well.
      We can only rule it out if its dividing hyperplane does not intersect
      the circle of radius "range" centered at "query" */
    addedRes += ret;
    ret = FindNearest(dx <= 0.0 ? node->right : node->left, query, range, list, ordered, dim,ierr);
  }

  if (*ierr != PX_NO_ERROR)
  {
    return -1;
  }
  addedRes += ret;

  return addedRes;
}

/******************************************************************/
KDResult *KDNearestRange(KDTree *kd, Real const *query,
                         Real const range, int *ierr)
{
  int ret; /*variable to store the size of the KDResult object */
  KDResult *rset; /* the KDResult object used to track the output.
                   * It holds a linked list which contain the points
                   * that are within range of the query point. */

  /*allocate KDResult struct */
  //*ierr = Allocate(1,sizeof(KDResult),(void **)&(rset) );
  //if (*ierr != PX_NO_ERROR) return rset;
  rset = (KDResult*)malloc(sizeof(KDResult));


  *ierr = AllocResultNode(&(rset->rlist));
  if (*ierr != PX_NO_ERROR)
  {
    KDResultFree(rset);
    return NULL;
  }
  //    free( rset );
  //  return 0;

  rset->rlist->next = 0;
  rset->tree = kd;

  /*Perform the search for neighbors within distance range */
  ret = FindNearest(kd->root, query, range, rset->rlist, 0, kd->dim,ierr);
  if (*ierr != PX_NO_ERROR)
  {
    KDResultFree(rset);
    return NULL;
  }
  rset->size = ret; /*set size*/
  KDResultRewind(rset); /*rewind iterator to start */
  return rset;
}

/******************************************************************/
void KDResultFree(KDResult *rset)
{
  ClearResults(rset);
  FreeResultNode(rset->rlist);
  free( rset );
}

/******************************************************************/
int KDResultSize(KDResult const *set)
{
  return (set->size);
}

/******************************************************************/
void KDResultRewind(KDResult *rset)
{
  rset->riter = rset->rlist->next;
}

/******************************************************************/
int KDResultEnd(KDResult const *rset)
{
  return rset->riter == 0;
}

/******************************************************************/
int KDResultNext(KDResult *rset)
{
  rset->riter = rset->riter->next;
  return rset->riter != 0;
}

/******************************************************************/
void *KDResultItem(KDResult const *rset, Real *pos)
{
  if (rset->riter) 
  {
    if (pos) 
    {
      memcpy(pos, rset->riter->pos, rset->tree->dim * sizeof *pos);
    }
    return rset->riter->data;
  }
  return 0;
}


/******************************************************************/
void *KDResultItemData(KDResult const *set)
{
  return KDResultItem(set, 0);
}


/* ------------------ LinkedList functions ------------------- */
/******************************************************************/
int LListInsert(ListHead *head, void *data,
                Real const key, int const useListIter)
{
  // int ierr;
  LListNode *insertNode=NULL; /* node to be inserted*/
  LListNode *tempNode; /* pointer to step through the LinkedList*/

  if (useListIter == 0)
  {
    return LListInsertFirst(head, data, key);
  }
  if (head->listIter == NULL)
  {
    if (head->ordered == 0)
    {
      return PX_BAD_INPUT; /* error: cannot insert node after a null node!*/
    }else
      return LListInsertFirst(head, data, key);
  }

  tempNode = head->listIter;
  if (head->ordered != 0) 
  {
    /* Insert after listIter but check that ordering is not violated*/
    if (tempNode->key < key && (tempNode->next == NULL || tempNode->next->key >= key))
    {
      /* insert */


      PXErrorReturn( AllocLListNode(&insertNode) );
      insertNode->data = data;
      insertNode->key = key;

      insertNode->next = tempNode->next;
      tempNode->next = insertNode;
      head->size += 1;
      return PX_NO_ERROR;
    }
    else
    {
      /* Inserting here would violate the ordering.  Find the correct spot.*/
      return LListInsertFirst(head, data, key);
    }
  }
  else
  {
    /* insert after listIter */


    PXErrorReturn( AllocLListNode(&insertNode) );
    insertNode->data = data;
    insertNode->key = key;

    insertNode->next = tempNode->next;
    tempNode->next = insertNode;
    head->size += 1;
    return PX_NO_ERROR;
  }
}

/******************************************************************/
int LListInsertFirst(ListHead *head, void *data,
                     Real const key)
{
  // int ierr;
  LListNode *insertNode=NULL; /* node to be inserted */
  LListNode *tempNode; /* pointer to step through the LinkedList */


  PXErrorReturn( AllocLListNode(&insertNode) );
  insertNode->data = data;
  insertNode->key = key;

  tempNode = head->llist;
  if (head->ordered != 0) 
  {
    while (tempNode->next != NULL && tempNode->next->key < key) 
    {
      tempNode = tempNode->next;
    }
  }
  insertNode->next = tempNode->next;
  tempNode->next = insertNode;
  head->size += 1;
  return PX_NO_ERROR;
}


/******************************************************************/
ListHead *LListCreate(int const ordered, int *ierr)
{
  ListHead *head=NULL; /* holds the ListHead object for the ListHead to be created */

  if (ordered != 0 && ordered != 1)
  {
    (*ierr) = PX_BAD_INPUT;
    return NULL;
  }

  /* allocate space for ListHead struct */
  //*ierr = Allocate(1,sizeof(ListHead),(void **)&(head) );
  //if (*ierr != PX_NO_ERROR) return head;
  head = (ListHead*)malloc(sizeof(ListHead));

  /* allocate the head-node of this ListHead's LinkedList */
  *ierr = AllocLListNode(&(head->llist));
  if (*ierr != PX_NO_ERROR)
  {
    free(head );
    return NULL;
  }

  (*ierr) = PX_NO_ERROR;
  /* set pointers/fields appropriately */
  head->llist->next = NULL;
  head->llist->data = NULL;

  head->size = 0;
  head->listIter = NULL;
  head->ordered = ordered;
  head->destr = NULL;
  return head;
}


/******************************************************************/
void LListFree(ListHead *head)
{
  ClearLList(head);
  FreeLListNode(head->llist);
  free( head );
}

/******************************************************************/
int LListSize(ListHead const *head)
{
  return (head->size);
}

/******************************************************************/
void LListRewind(ListHead *head)
{
  head->listIter = head->llist->next;
}

/******************************************************************/
int LListIterIsNull(ListHead const *head)
{
  return head->listIter == NULL;
}

/******************************************************************/
int LListAtLast(ListHead const *head)
{
  return head->listIter->next == NULL;
}

/******************************************************************/
int LListNext(ListHead *head)
{
  head->listIter = head->listIter->next;
  return head->listIter != NULL;
}

/******************************************************************/
void *LListItem(ListHead const *head)
{
  if (head->listIter != NULL) 
  {
    return head->listIter->data;
  }
  return NULL;
}

/******************************************************************/
Real LListKey(ListHead const *head)
{
  if (head->listIter != NULL) 
  {
    return head->listIter->key;
  }
  printf("ERROR: Asked for the key of a NULL node!");
  return std::numeric_limits<double>::quiet_NaN();
}

/******************************************************************/
void LListDataDestructor(ListHead *head, void (*destr)(void*))
{
  head->destr = destr;
}


/* ---- static helpers ---- */
/* special list node allocators, if needed */
#ifdef USE_LIST_NODE_ALLOCATOR
static ResultNode *free_nodes = NULL;

//#ifdef USE_LIST_NODE_ALLOCATOR
static LListNode *free_listnodes = NULL;


/******************************************************************/
static int AllocResultNode(ResultNode **node)
{
  // int ierr; /* error code */
  int ierr;;

  /*if there are no free nodes available, allocate a new node*/
  if (!free_nodes) 
  {
    //ierr = Allocate(1,sizeof(ResultNode),(void **)&(*node) );
    *node = (ResultNode*)malloc(sizeof(ResultNode));
    ierr = PX_NO_ERROR;
  } 
  else 
  {
    /*otherwise, grab an un-used node from the free nodes list.*/
    *node = free_nodes;
    free_nodes = free_nodes->next;
    (*node)->next = 0;
    ierr = PX_NO_ERROR;
  }


  return ierr;
}

/******************************************************************/
static void FreeResultNode(ResultNode *node)
{
  /*nodes are not immediately deleted: instead they are added to
    a list of free nodes (free_listnodes) which are used to
    quickly allocate additional nodes later. */
  node->next = free_nodes;
  free_nodes = node;

}
//#endif        /* list node allocator or not */


/******************************************************************/
static int AllocLListNode(LListNode **node)
{
  //int ierr; /* error code*/
  int ierr; //LListNode *node; /*pointer to a newly allocated node, if needed */

  /* obtain mutex lock if needed */

  /*if there are no free nodes available, allocate a new node*/
  if (free_listnodes == NULL) 
  {
    //ierr = Allocate(1,sizeof(LListNode),(void **)&(*node) );
    *node = (LListNode*)malloc(sizeof(LListNode));
    ierr = PX_NO_ERROR;
  } 
  else 
  {
    /*otherwise, grab an un-used node from the free nodes list.*/
    *node = free_listnodes;
    free_listnodes = free_listnodes->next;
    (*node)->next = NULL;
    ierr = PX_NO_ERROR;
  }

  /* unlock */

  return ierr;
}

/******************************************************************/
static void FreeLListNode(LListNode *node)
{
  /*if needed, obtain mutex lock*/

  /*nodes are not immediately deleted: instead they are added to
    a list of free nodes (free_listnodes) which are used to
    quickly allocate additional nodes later. */
  node->next = free_listnodes;
  free_listnodes = node;

  /*release lock*/
}

/******************************************************************/
void LListFinalize()
{
  LListNode *node;

  node = free_listnodes;
  while (NULL != free_listnodes)
  {
    node = free_listnodes;
    free_listnodes = free_listnodes->next;
    free( node );
  }
  /* at completion, free_listnodes will be NULL */
}

/******************************************************************/
void KDTreeFinalize()
{
  ResultNode *node;

  node = free_nodes;
  while (NULL != free_nodes)
  {
    node = free_nodes;
    free_nodes = free_nodes->next;
    free( node );
  }

}

#endif  /* list node allocator or not */


/******************************************************************/
/* inserts the item. if dist_sq is >= 0, then do an ordered insert */
static int ResultListInsert(ResultNode *list, KDNode *item, Real dist_sq)
{
  // int ierr;
  ResultNode *rnode=NULL; /* pointer to the node to be allocated */

  PXErrorReturn( AllocResultNode(&rnode) );

  //rnode->item = item;
  rnode->data = item->data;
  rnode->pos = item->pos;
  rnode->dist_sq = dist_sq;

  /*if dist_sq is non-negative, insert this node in sorted order.
    else place it at the head the list. */
  if (dist_sq >= 0.0) 
  {
    while (list->next && list->next->dist_sq > dist_sq) 
    {
      list = list->next;
    }
  }
  rnode->next = list->next;
  list->next = rnode;
  return PX_NO_ERROR;
}

/******************************************************************/
static void ClearResults(KDResult *rset)
{
  ResultNode *tmp; /* temp node for pointer dancing during deletion */
  ResultNode *node = rset->rlist->next;/*pointer to the first node*/

  while (node) 
  {
    tmp = node;
    node = node->next;
    FreeResultNode(tmp);
  }

  rset->rlist->next = 0;
}

/******************************************************************/
static void ClearLList(ListHead *head)
{
  LListNode *tmp; /* temp node for pointer dancing during deletion */
  LListNode *node = head->llist->next; /*pointer to the first node*/

  while (node != NULL) 
  {
    tmp = node;
    node = node->next;
    if (head->destr != NULL && tmp->data != NULL)
    {
      head->destr(tmp->data);
    }
    FreeLListNode(tmp);
  }

  head->llist->next = NULL;
}


/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/************************Stack FUNCTIONS *****************************/
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

/******************************************************************/
void *StackPop(ListHead *head, int *ierr)
{
  LListNode *tempNode; /* Node to be removed */
  void *data;
  tempNode = head->llist->next; /* set tempNode to the first node */
  if (NULL == tempNode)
  {
    printf("Tried to pop an empty stack!\n");
    // (*ierr) = PXError(PX_CODE_FLOW_ERROR);
    PXError((*ierr),PX_CODE_FLOW_ERROR);
    return NULL;
  }

  /* unlink tempNode from the Stack */
  head->llist->next = tempNode->next;
  data = tempNode->data;

  FreeLListNode(tempNode);
  head->size -= 1;
  return data;
}

/******************************************************************/
int StackPush(ListHead *head, void *data)
{
  // int ierr;
  LListNode *insertNode=NULL;

  PXErrorReturn( AllocLListNode(&insertNode) );
  insertNode->data = data;
  /* insertNode->key = 0.0; */
  insertNode->next = head->llist->next;
  head->llist->next = insertNode;
  head->size += 1;
  return PX_NO_ERROR;

}

/******************************************************************/
void *StackPopWithKey(ListHead *head, Real *key, int *ierr)
{
  LListNode *tempNode; /* Node to be removed */
  void *data;
  tempNode = head->llist->next; /* set tempNode to the first node */
  if (NULL == tempNode)
  {
    printf("Tried to pop an empty stack!\n");
    // (*ierr) = PXError(PX_CODE_FLOW_ERROR);
    PXError((*ierr),PX_CODE_FLOW_ERROR);
    return NULL;
  }

  /* unlink tempNode from the Stack */
  head->llist->next = tempNode->next;
  data = tempNode->data;
  *key = tempNode->key;

  FreeLListNode(tempNode);
  head->size -= 1;
  return data;
}

/******************************************************************/
int StackPushWithKey(ListHead *head, void *data, Real const key)
{
  // int ierr;
  LListNode *insertNode=NULL;

  PXErrorReturn( AllocLListNode(&insertNode) );
  insertNode->data = data;
  insertNode->key = key;
  insertNode->next = head->llist->next;
  head->llist->next = insertNode;
  head->size += 1;
  return PX_NO_ERROR;

}


/******************************************************************/
ListHead *StackCreate(int *ierr, int const prealloc)
{
  return LListCreate(0, ierr);

  ListHead *head;
  LListNode *preallocNodes;
  head = LListCreate(0,ierr);
  if ((*ierr) != PX_NO_ERROR)
  {
    return NULL;
  }
  //PXErrorReturn(*ierr);
  if (prealloc != 0)
  {
    printf("WARNING: Preallocation not yet supported!\n");

    //*ierr = Allocate(prealloc,sizeof(LListNode), (void **)&(preallocNodes) );
    preallocNodes = (LListNode*)malloc(prealloc*sizeof(LListNode));

  }
  //return head;
}

/******************************************************************/
void StackFree(ListHead *head)
{
  if (head->size != 0)
  {
    printf("Warning: free-ing non-empty stack!\n");
  }
  ClearLList(head);
  FreeLListNode(head->llist);
  free( head );
}



