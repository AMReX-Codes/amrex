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


#ifndef _KDSTRUCT_H_
#define _KDSTRUCT_H_

#include "AMReX_REAL.H"

/*!
  \file   PXKDStruct.h

   This header file defines the structures used to create KDTrees
   and LinkedLists in ProjectX.
*/


#define KD_MAX_DIM 4

typedef struct _kdnode KDNode;
struct _kdnode {
  void *data;
  struct _kdnode *left, *right; /* negative/positive side */
  Real pos[KD_MAX_DIM];
  Real xmin[KD_MAX_DIM];
  Real xmax[KD_MAX_DIM];
  unsigned char dir;
};

typedef struct _resultnode ResultNode;
struct _resultnode {
  //KDNode *item;
  void *data; //ptr to data from the kdnode
  Real *pos; //ptr to coordinates of kdnode
  Real dist_sq;
  struct _resultnode *next;
};


typedef struct _kdtree KDTree;
struct _kdtree {
  KDNode *root;
  int dim;
  char globalDestrFlag;
  /* value: -1 => destr == NULL == globalData, no node-based data to free
             0 => destr != NULL, globalData == NULL, destr called on each node's data ptr
             1 => destr != NULL, globalData != NULL, destr called on globalData ptr (only once)
  */
  void *globalData;
  void (*destr)(void*); //function called to release memory
};

typedef struct _kdres KDResult;
struct _kdres {
        KDTree *tree;
        ResultNode *rlist, *riter;
        int size;
};

typedef struct _linkedlistnode LListNode;
struct _linkedlistnode {
  void *data;
  Real key;
  struct _linkedlistnode *next;
};

typedef struct _linkedlisthead ListHead;
struct _linkedlisthead {
  LListNode *llist, *listIter;
  int size;
  int ordered;
  void (*destr)(void*);
};


#endif
