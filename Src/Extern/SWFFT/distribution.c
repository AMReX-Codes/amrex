/*
 *                 Copyright (C) 2017, UChicago Argonne, LLC
 *                            All Rights Reserved
 *
 *           Hardware/Hybrid Cosmology Code (HACC), Version 1.0
 *
 * Salman Habib, Adrian Pope, Hal Finkel, Nicholas Frontiere, Katrin Heitmann,
 *      Vitali Morozov, Jeffrey Emberson, Thomas Uram, Esteban Rangel
 *                        (Argonne National Laboratory)
 *
 *  David Daniel, Patricia Fasel, Chung-Hsing Hsu, Zarija Lukic, James Ahrens
 *                      (Los Alamos National Laboratory)
 *
 *                               George Zagaris
 *                                 (Kitware)
 *
 *                            OPEN SOURCE LICENSE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *   1. Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer. Software changes,
 *      modifications, or derivative works, should be noted with comments and
 *      the author and organization's name.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the names of UChicago Argonne, LLC or the Department of Energy
 *      nor the names of its contributors may be used to endorse or promote
 *      products derived from this software without specific prior written
 *      permission.
 *
 *   4. The software and the end-user documentation included with the
 *      redistribution, if any, must include the following acknowledgment:
 *
 *     "This product includes software produced by UChicago Argonne, LLC under
 *      Contract No. DE-AC02-06CH11357 with the Department of Energy."
 *
 * *****************************************************************************
 *                                DISCLAIMER
 * THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND. NEITHER THE
 * UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR 
 * UCHICAGO ARGONNE, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, 
 * EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE
 * ACCURARY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, DATA, APPARATUS,
 * PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
 * PRIVATELY OWNED RIGHTS.
 *
 * *****************************************************************************
 */

#include <assert.h>
#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "distribution_c.h"

#ifndef USE_SLAB_WORKAROUND
#define USE_SLAB_WORKAROUND 0
#endif

enum {
  REDISTRIBUTE_1_TO_3,
  REDISTRIBUTE_3_TO_1,
  REDISTRIBUTE_2_TO_3,
  REDISTRIBUTE_3_TO_2
};

//#define DEBUG_CONDITION (self == 0 || self == 1)
#define DEBUG_CONDITION false

// return comma or period depending on position in a list
static inline char *separator(int i, int n)
{
  return i == (n - 1) ? "." : ", ";
}


//Go from rank of processor to its cartesian coords, and vica versa. 
//Assumes the ranks increment in the z dimension then y then x.
void Coord_cube(int myrank, 
		int coord[], 
		distribution_t *d)
{
  coord[0]=myrank/(d->process_topology_3.nproc[1]*d->process_topology_3.nproc[2]);
  coord[1]=(myrank%(d->process_topology_3.nproc[1]*d->process_topology_3.nproc[2]))/(d->process_topology_3.nproc[2]);
  coord[2]=(myrank%(d->process_topology_3.nproc[1]*d->process_topology_3.nproc[2]))%(d->process_topology_3.nproc[2]);
  return;
}

void Rank_cube(int * myrank, 
	       int coord[], 
	       distribution_t *d)
{
  *myrank = coord[2] + (d->process_topology_3.nproc[2])*(coord[1] + d->process_topology_3.nproc[1]*coord[0]);
  return;
}


/*
  The subsequent member functions are used to look up and ranks of x,y, and z 
  pencils from their coordinates, and vica versa.
  The ordering of the ranks is such that pencils will be going through cubes 
  with the same rank sequencially. (since the cubes ranks are already 
  determined and can not be changed, these routines figure out which ranks 
  the pencils should be assigned so that there is no communication hangs.)
*/
void Coord_x_pencils(int myrank, 
		     int coord[], 
		     distribution_t *d)
{
  // asserts only one processor in x_direction
  assert(d->process_topology_2_x.nproc[0] == 1);
  //since x_pencils only have one processor in the x_direction.
  coord[0]=0;
  int num_pen_in_cube_col=d->process_topology_2_x.nproc[1]/d->process_topology_3.nproc[1];
  int num_pen_in_cube_row=d->process_topology_2_x.nproc[2]/d->process_topology_3.nproc[2];
  int num_cubes=(d->process_topology_3.nproc[2]*d->process_topology_3.nproc[1]);
  
/*
  the x_pencil ranks increment in each cube sequencially, after reaching the 
  last cube the second slot in the first cube is the next rank, and then the 
  process repeats. num_repeats, is the number of times this repetition had to 
  have occured to increment to the current rank.
*/
  int num_repeats=myrank/(num_cubes);
  
  //now subtract the difference of how many repetitions, to find the lowest 
  //rank in the cube it resides. 
  int low_rank=myrank-num_repeats*num_cubes;
  
  //find the y and z coords of the low_rank, then adjust coords for ranks 
  //that repeated around the cube.
  coord[1] = (low_rank/d->process_topology_3.nproc[2])*num_pen_in_cube_col 
    + num_repeats%num_pen_in_cube_col;
  coord[2] = (low_rank%d->process_topology_3.nproc[2])*num_pen_in_cube_row + num_repeats/num_pen_in_cube_col;
    
  return;
}

void Rank_x_pencils(int * myrank, 
		    int coord[], 
		    distribution_t *d)
{
  int num_pen_in_cube_col=d->process_topology_2_x.nproc[1]/d->process_topology_3.nproc[1];
  int num_pen_in_cube_row=d->process_topology_2_x.nproc[2]/d->process_topology_3.nproc[2];
  if(num_pen_in_cube_col == 0)
    fprintf(stderr,"num_cube_col%d ", 
	    d->process_topology_2_x.nproc[1]/d->process_topology_3.nproc[1]);
  if(num_pen_in_cube_row ==0)
    fprintf(stderr,"num_cube_row%d ", d->process_topology_3.nproc[2]);
  assert(num_pen_in_cube_col !=0 && num_pen_in_cube_row !=0);
  int alpha = coord[1]%num_pen_in_cube_col;
  int num_cubes = (d->process_topology_3.nproc[2]*d->process_topology_3.nproc[1]);
  int beta = coord[2]%num_pen_in_cube_row;
  *myrank = (alpha*num_cubes) 
    + ((coord[1]/num_pen_in_cube_col)*d->process_topology_3.nproc[2]) 
    + (beta*(num_cubes)*num_pen_in_cube_col) + coord[2]/num_pen_in_cube_row;
  return;
}

void Coord_y_pencils(int myrank, 
		     int coord[], 
		     distribution_t *d)
{
  // asserts only one processor in y_direction
  assert(d->process_topology_2_y.nproc[1] == 1);
  //since y_pencils only have one processor in the y_direction.
  coord[1] = 0;
  int num_pen_in_cube_row = d->process_topology_2_y.nproc[2]/d->process_topology_3.nproc[2];
  int alpha = myrank%(d->process_topology_2_y.nproc[2]);
  coord[0] = myrank/d->process_topology_2_y.nproc[2];
  
  coord[2] = (alpha/d->process_topology_3.nproc[2]) 
    + (alpha%d->process_topology_3.nproc[2])*num_pen_in_cube_row;
  
  return;
}

void Rank_y_pencils(int * myrank, 
		    int coord[], 
		    distribution_t *d)
{
  int num_pen_in_cube_col = d->process_topology_2_y.nproc[0]/d->process_topology_3.nproc[0];
  int num_pen_in_cube_row = d->process_topology_2_y.nproc[2]/d->process_topology_3.nproc[2];
  //WHY ARE THESE COMMENTED OUT?
  //if(num_pen_in_cube_col ==0)fprintf(stderr,"num_cube_col%d ", d->process_topology_2_y.nproc[1]/d->process_topology_3.nproc[1]);
  //if(num_pen_in_cube_row ==0)fprintf(stderr,"num_cube_row%d ", d->process_topology_3.nproc[2]);
  assert(num_pen_in_cube_col !=0 && num_pen_in_cube_row !=0);
  int beta = coord[2]%num_pen_in_cube_row;
  *myrank = coord[0]*d->process_topology_2_y.nproc[2] 
    + beta*d->process_topology_3.nproc[2] 
    + coord[2]/num_pen_in_cube_row;
  return;
}

void Coord_z_pencils(int myrank, 
		     int coord[], 
		     distribution_t *d)
{
  // asserts only one processor in z_direction
  assert(d->process_topology_2_z.nproc[2] == 1);
  //since z_pencils only have one processor in the z_direction.
  coord[2] = 0;
  int num_pen_in_cube_col = d->process_topology_2_z.nproc[1]/d->process_topology_3.nproc[1];
  int num_pen_in_cube_row = d->process_topology_2_z.nproc[0]/d->process_topology_3.nproc[0];
  int num_pen_in_cube = d->process_topology_3.nproc[2];
  int alpha = myrank/(d->process_topology_2_z.nproc[1]*num_pen_in_cube_row);
  coord[0] = alpha*num_pen_in_cube_row + (myrank%num_pen_in_cube)/num_pen_in_cube_col;
  coord[1] = ((myrank%(d->process_topology_2_z.nproc[1]*num_pen_in_cube_row))/num_pen_in_cube)*num_pen_in_cube_col + myrank%num_pen_in_cube_col;
  
  return;
}

void Rank_z_pencils(int * myrank, 
		    int coord[], 
		    distribution_t *d)
{
  int num_pen_in_cube_col = d->process_topology_2_z.nproc[1]/d->process_topology_3.nproc[1];
  int num_pen_in_cube_row = d->process_topology_2_z.nproc[0]/d->process_topology_3.nproc[0];
  int num_pen_in_cube = d->process_topology_3.nproc[2];
  if(num_pen_in_cube_col == 0)
    fprintf(stderr,"num_cube_col%d ", 
	    d->process_topology_2_z.nproc[1]/d->process_topology_3.nproc[1]);
  if(num_pen_in_cube_row == 0)
    fprintf(stderr,"num_cube_row%d ", d->process_topology_3.nproc[2]);
  assert(num_pen_in_cube_col !=0 && num_pen_in_cube_row !=0);
  int alpha = coord[1]%num_pen_in_cube_col;
  int beta = coord[0]%num_pen_in_cube_row;
  *myrank = alpha 
    + ((coord[1]/num_pen_in_cube_col)*num_pen_in_cube) 
    + (beta*num_pen_in_cube_col) 
    + (coord[0]/num_pen_in_cube_row)*d->process_topology_2_z.nproc[1]*num_pen_in_cube_row;
  return;
}


// create 1-, 2- and 3-d cartesian data distributions comm MPI Communicator
void distribution_init(MPI_Comm comm, 
		       const int n[], 
		       const int Ndims[], 
		       distribution_t *d, 
               const int* rmap,
		       bool debug)
{
/* 
   As of 09/06/2011 The MPI function MPI_Dims_create is used to come up with 
   the most evenly distributed number of processors for the 3d distribution. 
   Since this can actually vary between machines, we should later write our own
   prime factorization function that does that for us. For the 2d distribution 
   pencils, Dims_create is also used, but the code then checks if the pencils 
   it outputs fits inside the 3d cuboids that were created. If it does not, it 
   tries swapping the dimensions of the pencils, and if they still do not fit, 
   it takes the 3d cubes dimensions of processors (np1,np2,np3) and (for 
   Z-pencils for example) makes pencils of the form (np1*np3,np2,1), or 
   (np1,np2*np3,1) which fit inside the cubes. However, even then, it is not 
   ensured that this will work since for example, if np1*np3 is bigger then 
   the number of points in one dimension (Ng) then there are not enough points 
   for each processor to have at least one point in that dimension. So the code
   checks this and asserts three variables check_x_dims check_y_dims, and 
   check_z_dims, which will assert if these kinda errors happen (as well as 
   checking errors coming from picking the total number of processors and Ng 
   in a way where the cubes will not fit for any orientation (like 100 procs 
   and Ng=101!)). Curretly the fix to these errors is to pick better values 
   for Ng and the total number of processors that work, however when we do 
   have our own prime factorization method, then that method could also make 
   pencils that fit inside the proper distribution (and we would not need so 
   many checks). In the mean time, to pick these "better" values for Ng, the 
   user should pick values such that: Ng % np1, Ng % np2, and Ng % np3 all 
   equal zero, and that np1*np2, np2*np3, and np3*np1 are all less then Ng.
   (in other words, the cubes created fit inside the number of grid points, 
   and the number of pencils created is not more then the number of points 
   in a dimension (Ng)).
*/
  d->parent = comm;

  int nproc;//num_processors
  int self; //rank
  int ndim = 3;
  int period[3];

  MPI_Comm_rank(comm, &self);
  MPI_Comm_size(comm, &nproc);

  // Construct the rankmap[grid] --> rank map (converts grids to ranks) from the input 
  // (if none is provided then we assume the trivial map {0, 1, ..., nproc}).
  d->rankmap  = (int *) malloc(sizeof(int)*nproc);
  if(rmap) for(int i=0; i<nproc; i++) d->rankmap[i] = rmap[i];
  else for(int i=0; i<nproc; i++) d->rankmap[i] = i;

  // Construct the inverse gridmap[rank] --> grid map (converts ranks to grids)
  d->gridmap = (int *) malloc(sizeof(int)*nproc);
  for(int i=0; i<nproc; i++) {
    for(int j=0; j<nproc; j++) {
      if(i == d->rankmap[j]) { d->gridmap[i] = j; break; }
    }
  }

  // Map this rank to the correct grid
  self = d->gridmap[self]; 

  if (!self) 
    printf("Initializing redistribution using a %s layout on %d ranks.\n",
#ifdef PENCIL
	   "pencil"
#else
	   "slab"
#endif
	   ,nproc);

  d->debug = debug;
  for (int i = 0; i < 3; ++i)
    d->n[i] = n[i];

  // set up process grid with 1d decomposition (SLABs)
  d->process_topology_1.nproc[0] = 0;
  d->process_topology_1.nproc[1] = 1; // don't distribute outer dimensions
  d->process_topology_1.nproc[2] = 1; // don't distribute outer dimensions
  period[0] = period[1] = period[2] = 1;
  //process_topology_1.nproc is filled with number of processors in each dim
  MPI_Dims_create(nproc, ndim, d->process_topology_1.nproc); 

  if(self == 0) {
    printf("distribution 1D: [%d:%d:%d]\n",
	   d->process_topology_1.nproc[0],
	   d->process_topology_1.nproc[1],
	   d->process_topology_1.nproc[2]);
    fflush(stdout);
  }

  if (d->debug && 0 == self) {
    fprintf(stderr, "Process grids:\n");
    fprintf(stderr, "  1d: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", 
	      d->process_topology_1.nproc[i], 
	      separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }
  //creates the new communicator
  MPI_Cart_create(comm, ndim, d->process_topology_1.nproc, period, 0, 
		  &d->process_topology_1.cart);
  //gets .self (is coordinate)
  MPI_Cart_get(d->process_topology_1.cart, ndim, d->process_topology_1.nproc, 
	       d->process_topology_1.period, d->process_topology_1.self);
  //calculates the local dimensions (number of points in each dimension)
  d->process_topology_1.n[0] = n[0] / d->process_topology_1.nproc[0];
  d->process_topology_1.n[1] = n[1] / d->process_topology_1.nproc[1];
  d->process_topology_1.n[2] = n[2] / d->process_topology_1.nproc[2];



  // set up process grid with 3d decomposition (CUBE)
  d->process_topology_3.nproc[0] = 0;
  d->process_topology_3.nproc[1] = 0;
  d->process_topology_3.nproc[2] = 0;
  period[0] = period[1] = period[2] = 1;
  Custom3D_Dims_create(Ndims, nproc, ndim, d->process_topology_3.nproc);

  if(self == 0) {
    printf("distribution 3D: [%d:%d:%d]\n",
	   d->process_topology_3.nproc[0],
	   d->process_topology_3.nproc[1],
	   d->process_topology_3.nproc[2]);
    fflush(stdout);
  }
  
  if (d->debug && 0 == self) {
    fprintf(stderr, "  3d: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", 
	      d->process_topology_3.nproc[i], 
	      separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }

  MPI_Cart_create(comm, ndim, d->process_topology_3.nproc, period, 0, 
		  &d->process_topology_3.cart);
  //finds cartesian coordinate of this current rank
  Coord_cube(self,d->process_topology_3.self,d);

  if(debug){
/*
  this debug statment checks to see if the way coordinates found by 
  calculation matchs MPI's coord system (MPI might differ between machines 
  so this is why the code calculates the coord system itself, however with 
  debug on, can check if it matches MPI(even tho it is not enforced to match 
  it.)).
*/
    int prev_coord[3];
    prev_coord[0]=d->process_topology_3.self[0];
    prev_coord[1]=d->process_topology_3.self[1];
    prev_coord[2]=d->process_topology_3.self[2];
    MPI_Cart_get(d->process_topology_3.cart, ndim, 
		 d->process_topology_3.nproc, 
		 d->process_topology_3.period, 
		 d->process_topology_3.self);
    for(int i=0; i < 3; i++)
      if(prev_coord[i] != d->process_topology_3.self[i])
	abort();
  }
  assert(n[0]%d->process_topology_3.nproc[0] == 0);
  assert(n[0]%d->process_topology_3.nproc[1] == 0);
  assert(n[0]%d->process_topology_3.nproc[2] == 0);
  
  //set local dimensions
  d->process_topology_3.n[0] = n[0] / d->process_topology_3.nproc[0];
  d->process_topology_3.n[1] = n[1] / d->process_topology_3.nproc[1];
  d->process_topology_3.n[2] = n[2] / d->process_topology_3.nproc[2];



  // set up process grid with 2d decomposition (z_PENCILs )
  d->process_topology_2_z.nproc[0] = 0;
  d->process_topology_2_z.nproc[1] = 0;
  d->process_topology_2_z.nproc[2] = 1; // don't distribute outer dimension 
  period[0] = period[1] = period[2] = 1;
  MPI_Dims_create(nproc, ndim, d->process_topology_2_z.nproc);
  d->process_topology_2_z.n[0] = n[0] / d->process_topology_2_z.nproc[0];
  d->process_topology_2_z.n[1] = n[1] / d->process_topology_2_z.nproc[1];
  d->process_topology_2_z.n[2] = n[2] / d->process_topology_2_z.nproc[2];
  //variable used to ensure that pencils created fit inside the cuboids, 
  //if not the code will assert out.
  bool check_z_dims=false; 
  if(d->process_topology_2_z.n[0] != 0 
     && d->process_topology_2_z.n[1] != 0 
     && d->process_topology_2_z.n[2] != 0)
  {// protects from dividing by zero.
    check_z_dims = ((d->process_topology_3.n[0]) % (d->process_topology_2_z.n[0]) == 0) 
      && ((d->process_topology_3.n[1]) % (d->process_topology_2_z.n[1]) == 0) 
      && (n[0] % (d->process_topology_2_z.nproc[0]) == 0) 
      && (n[0] % (d->process_topology_2_z.nproc[1]) == 0);
    
    if(self==0 && debug && !check_z_dims)
      fprintf(stderr,"Need to fix Z PENCILS z_procs(%d,%d,%d) 3d.ns(%d,%d,%d) 2d_z.ns(%d,%d,%d).... \n", 
	      d->process_topology_2_z.nproc[0],
	      d->process_topology_2_z.nproc[1],
	      d->process_topology_2_z.nproc[2],
	      d->process_topology_3.n[0],
	      d->process_topology_3.n[1],
	      d->process_topology_3.n[2],
	      d->process_topology_2_z.n[0],
	      d->process_topology_2_z.n[1],
	      d->process_topology_2_z.n[2]);
   
    //try swaping pencil dimensions if current setup pencil dimensions dont 
    //fit inside the cubes.
    if(!(check_z_dims) 
       && ((d->process_topology_3.n[0]) % (d->process_topology_2_z.n[1]) == 0) 
       && ((d->process_topology_3.n[1]) % (d->process_topology_2_z.n[0]) == 0))
    {

      if(self==0 && debug)
	fprintf(stderr,"Swaping Z pencils in initialization  (%d,%d,%d).... \n", 
		d->process_topology_2_z.nproc[0],
		d->process_topology_2_z.nproc[1],
		d->process_topology_2_z.nproc[2]);
      int temp=d->process_topology_2_z.nproc[0];
      d->process_topology_2_z.nproc[0] = d->process_topology_2_z.nproc[1];
      d->process_topology_2_z.nproc[1] = temp;
      d->process_topology_2_z.nproc[2] = d->process_topology_2_z.nproc[2];
      
      d->process_topology_2_z.n[0] = n[0] / d->process_topology_2_z.nproc[0];
      d->process_topology_2_z.n[1] = n[1] / d->process_topology_2_z.nproc[1];
      d->process_topology_2_z.n[2] = n[2] / d->process_topology_2_z.nproc[2];
      check_z_dims = ((d->process_topology_3.n[0]) % (d->process_topology_2_z.n[0]) == 0) 
	&& ((d->process_topology_3.n[1]) % (d->process_topology_2_z.n[1]) == 0)
	&& (n[0] % (d->process_topology_2_z.nproc[0]) == 0) 
	&& (n[0] % (d->process_topology_2_z.nproc[1]) == 0);
    }
  } else {
    check_z_dims=false;
  }
  /*
    if that did not work, make a pencil that does if inside the 3d cuboids by 
    taking the cuboids dimensions (np1,np2,np3) and making pencils 
    (np1,np2*np3,1), or (np1*np3,np2,1) on the most evenly distributed 
    dimensions
  */
  if(!check_z_dims){
    if(self==0 && debug)
      fprintf(stderr,"MAKING Z PENCILS FIT zprocs(%d,%d,%d) z.ns(%d,%d,%d).... \n", 
	      d->process_topology_2_z.nproc[0],
	      d->process_topology_2_z.nproc[1],
	      d->process_topology_2_z.nproc[2],
	      d->process_topology_2_z.n[0],
	      d->process_topology_2_z.n[1],
	      d->process_topology_2_z.n[2]);
    
    d->process_topology_2_z.nproc[2]=1;
    if(d->process_topology_3.n[0]>d->process_topology_3.n[1])
    {
      d->process_topology_2_z.nproc[1]=d->process_topology_3.nproc[1]*d->process_topology_3.nproc[2];
      d->process_topology_2_z.nproc[0]=d->process_topology_3.nproc[0];
      if((n[0] % (d->process_topology_2_z.nproc[0]) != 0) 
	 || (n[0] % (d->process_topology_2_z.nproc[1]) != 0))
      {
	d->process_topology_2_z.nproc[0]=d->process_topology_3.nproc[0]*d->process_topology_3.nproc[2];
	d->process_topology_2_z.nproc[1]=d->process_topology_3.nproc[1];
      }
    } else {
      d->process_topology_2_z.nproc[0]=d->process_topology_3.nproc[0]*d->process_topology_3.nproc[2];
      d->process_topology_2_z.nproc[1]=d->process_topology_3.nproc[1];
      if((n[0] % (d->process_topology_2_z.nproc[0]) != 0) 
	 || (n[0] % (d->process_topology_2_z.nproc[1]) != 0))
      {
	d->process_topology_2_z.nproc[1]=d->process_topology_3.nproc[1]*d->process_topology_3.nproc[2];
	d->process_topology_2_z.nproc[0]=d->process_topology_3.nproc[0];
      }
    }
    d->process_topology_2_z.n[0] = n[0] / d->process_topology_2_z.nproc[0];
    d->process_topology_2_z.n[1] = n[1] / d->process_topology_2_z.nproc[1];
    d->process_topology_2_z.n[2] = n[2] / d->process_topology_2_z.nproc[2];
    if(self==0 && debug)
      fprintf(stderr,"MAKING Z PENCILS FIT AFTER zprocs(%d,%d,%d) z.ns(%d,%d,%d)...\n", 
	      d->process_topology_2_z.nproc[0],
	      d->process_topology_2_z.nproc[1],
	      d->process_topology_2_z.nproc[2],
	      d->process_topology_2_z.n[0],
	      d->process_topology_2_z.n[1],
	      d->process_topology_2_z.n[2]);
    if(d->process_topology_2_z.n[0] != 0 
       && d->process_topology_2_z.n[1] != 0 
       && d->process_topology_2_z.n[2] != 0)
    {// protects from dividing by zero.
      check_z_dims=((d->process_topology_3.n[0]) % (d->process_topology_2_z.n[0]) == 0) 
	&& ((d->process_topology_3.n[1]) % (d->process_topology_2_z.n[1]) == 0)
	&& (n[0] % (d->process_topology_2_z.nproc[0]) == 0) 
	&& (n[0] % (d->process_topology_2_z.nproc[1]) == 0);
    } else {
      check_z_dims=false;
    }
  }
    
  if (d->debug && 0 == self) {
    fprintf(stderr, "  2d_z: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", 
	      d->process_topology_2_z.nproc[i], 
	      separator(i, ndim));
    }
    fprintf(stderr, "\n");
  } 
  if(!check_z_dims && debug && (self==0)){
    FILE * outfile;
    outfile= fopen("error.data","a");
    fprintf(outfile,"Z DIMS FAILS:(%d,%d,%d) (%d,%d,%d) \n",
	    d->process_topology_2_z.nproc[0],
	    d->process_topology_2_z.nproc[1],
	    d->process_topology_2_z.nproc[2], 
	    d->process_topology_3.nproc[0],
	    d->process_topology_3.nproc[1],
	    d->process_topology_3.nproc[2]);
  }
  assert(check_z_dims);
/*
  if this happends, it is because the dimensions were chosen incorrectly. 
  Either to many processors for the number of points in one dimenison (could 
  not do at least 1 point per processor), or the methods above could 
  not make a distribution of pencils that fit in the cubiods, which would 
  happen if the user gave numbers that wouldent work (we require the number 
  of processors in each dimension of the cuboid must be modulo the number of 
  points in that dimension, otherwise, this error will happen).
*/
  MPI_Cart_create(comm, 
		  ndim, 
		  d->process_topology_2_z.nproc, 
		  period, 
		  0, 
		  &d->process_topology_2_z.cart);
  //find the cartesian coord of the current rank (for the z_pencil)
  Coord_z_pencils(self,d->process_topology_2_z.self,d);

  if(self == 0) {
    printf("distribution 2z: [%d:%d:%d]\n",
	   d->process_topology_2_z.nproc[0],
	   d->process_topology_2_z.nproc[1],
	   d->process_topology_2_z.nproc[2]);
    fflush(stdout);
  }



  // set up process grid with 2d decomposition (x_PENCILs)
  d->process_topology_2_x.nproc[0] = 1; // don't distribute outer dimension
  d->process_topology_2_x.nproc[1] = 0;
  d->process_topology_2_x.nproc[2] = 0;
  period[0] = period[1] = period[2] = 1;
  MPI_Dims_create(nproc, ndim, d->process_topology_2_x.nproc);
  d->process_topology_2_x.n[0] = n[0] / d->process_topology_2_x.nproc[0];
  d->process_topology_2_x.n[1] = n[1] / d->process_topology_2_x.nproc[1];
  d->process_topology_2_x.n[2] = n[2] / d->process_topology_2_x.nproc[2];
  //variable used to ensure that pencils created fit inside the cuboids, 
  //if not the code will assert out.
  bool check_x_dims = false;
  if(d->process_topology_2_x.n[0] != 0 
     && d->process_topology_2_x.n[1] != 0 
     && d->process_topology_2_x.n[2] != 0)
  {// protects from dividing by zero.
    check_x_dims = ((d->process_topology_3.n[2]) % (d->process_topology_2_x.n[2]) == 0) 
      && ((d->process_topology_3.n[1]) % (d->process_topology_2_x.n[1]) == 0) 
      && (n[0] % (d->process_topology_2_x.nproc[2]) == 0) 
      && (n[0] % (d->process_topology_2_x.nproc[1]) == 0);
    if(self==0 && debug && !check_x_dims)
      fprintf(stderr,"Need to fix X PENCILS x_procs(%d,%d,%d) 3d.ns(%d,%d,%d) 2d_x.ns(%d,%d,%d)...\n", 
	      d->process_topology_2_x.nproc[0],
	      d->process_topology_2_x.nproc[1],
	      d->process_topology_2_x.nproc[2],
	      d->process_topology_3.n[0],
	      d->process_topology_3.n[1],
	      d->process_topology_3.n[2],
	      d->process_topology_2_x.n[0],
	      d->process_topology_2_x.n[1],
	      d->process_topology_2_x.n[2]);

    //try swaping pencil dimensions if current setup does not have pencils 
    //that fit inside cubes.
    if(!(check_x_dims) 
       && ((d->process_topology_3.n[2]) % (d->process_topology_2_x.n[1]) == 0) 
       && ((d->process_topology_3.n[1]) % (d->process_topology_2_x.n[2]) == 0))
    {
      if(self==0 && debug)
	fprintf(stderr,"Swaping X pencils in initialization .... \n");
      d->process_topology_2_x.nproc[0] = d->process_topology_2_x.nproc[0];
      int temp = d->process_topology_2_x.nproc[1];
      d->process_topology_2_x.nproc[1] = d->process_topology_2_x.nproc[2];
      d->process_topology_2_x.nproc[2] = temp;
   
      d->process_topology_2_x.n[0] = n[0] / d->process_topology_2_x.nproc[0];
      d->process_topology_2_x.n[1] = n[1] / d->process_topology_2_x.nproc[1];
      d->process_topology_2_x.n[2] = n[2] / d->process_topology_2_x.nproc[2];
      check_x_dims = ((d->process_topology_3.n[2]) % (d->process_topology_2_x.n[2]) == 0) 
	&& ((d->process_topology_3.n[1]) % (d->process_topology_2_x.n[1]) == 0)
	&& (n[0] % (d->process_topology_2_x.nproc[2]) == 0) 
	&& (n[0] % (d->process_topology_2_x.nproc[1]) == 0);
    } 
  } else{
    check_x_dims=false;
  }
  /*
    if that did not work, make a pencil that does by taking the cuboid 
    (np1,np2,np3) and making pencils of the form (1,np2*np1,np3) or 
    (1,np2*np1,np3) depending on the most even distribution it can.
  */
  if(!check_x_dims){
    if(self==0 && debug)
      fprintf(stderr,"MAKING X PENCILS FIT xprocs(%d,%d,%d) x.ns(%d,%d,%d)...\n",
	      d->process_topology_2_x.nproc[0],
	      d->process_topology_2_x.nproc[1],
	      d->process_topology_2_x.nproc[2],
	      d->process_topology_2_x.n[0],
	      d->process_topology_2_x.n[1],
	      d->process_topology_2_x.n[2]);

    d->process_topology_2_x.nproc[0] = 1;
    if(d->process_topology_3.nproc[2] > d->process_topology_3.nproc[1])
    {
      d->process_topology_2_x.nproc[1] = d->process_topology_3.nproc[1]*d->process_topology_3.nproc[0];
      d->process_topology_2_x.nproc[2] = d->process_topology_3.nproc[2];
      if((n[0] % (d->process_topology_2_x.nproc[2]) != 0) 
	 || (n[0] % (d->process_topology_2_x.nproc[0]) != 0))
      {
	d->process_topology_2_x.nproc[2]=d->process_topology_3.nproc[2]*d->process_topology_3.nproc[0];
	d->process_topology_2_x.nproc[1]=d->process_topology_3.nproc[1];
      }

    } else {
      d->process_topology_2_x.nproc[2] = d->process_topology_3.nproc[2]*d->process_topology_3.nproc[0];
      d->process_topology_2_x.nproc[1] = d->process_topology_3.nproc[1];
      if((n[0] % (d->process_topology_2_x.nproc[2]) != 0) 
	 || (n[0] % (d->process_topology_2_x.nproc[0]) != 0))
      {
	d->process_topology_2_x.nproc[1]=d->process_topology_3.nproc[1]*d->process_topology_3.nproc[0];
	d->process_topology_2_x.nproc[2]=d->process_topology_3.nproc[2];
      }
    }
    d->process_topology_2_x.n[0] = n[0] / d->process_topology_2_x.nproc[0];
    d->process_topology_2_x.n[1] = n[1] / d->process_topology_2_x.nproc[1];
    d->process_topology_2_x.n[2] = n[2] / d->process_topology_2_x.nproc[2];
    if(self==0 && debug)
      fprintf(stderr,"MAKING X PENCILS FIT AFTER xprocs(%d,%d,%d) x.ns(%d,%d,%d)...\n",
	      d->process_topology_2_x.nproc[0],
	      d->process_topology_2_x.nproc[1],
	      d->process_topology_2_x.nproc[2],
	      d->process_topology_2_x.n[0],
	      d->process_topology_2_x.n[1],
	      d->process_topology_2_x.n[2]);
    if(d->process_topology_2_x.n[0] != 0 
       && d->process_topology_2_x.n[1] != 0 
       && d->process_topology_2_x.n[2] != 0)
    {// protects from dividing by zero.
      check_x_dims = ((d->process_topology_3.n[2]) % (d->process_topology_2_x.n[2]) == 0) 
	&& ((d->process_topology_3.n[1]) % (d->process_topology_2_x.n[1]) == 0)
	&& (n[0] % (d->process_topology_2_x.nproc[2]) == 0) 
	&& (n[0] % (d->process_topology_2_x.nproc[1]) == 0);
    } else {
      check_x_dims=false;
    }  
  }
   
  if (d->debug && 0 == self) {
    fprintf(stderr, "  2d_x: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", 
	      d->process_topology_2_x.nproc[i], 
	      separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }
  if(!check_x_dims && debug && (self==0)){
    FILE * outfile;
    outfile= fopen("error.data","a");
    fprintf(outfile,"X DIMS FAILS:(%d,%d,%d) (%d,%d,%d) \n",
	    d->process_topology_2_x.nproc[0],
	    d->process_topology_2_x.nproc[1],
	    d->process_topology_2_x.nproc[2], 
	    d->process_topology_3.nproc[0],
	    d->process_topology_3.nproc[1],
	    d->process_topology_3.nproc[2]);
  }
  assert(check_x_dims);
/*
  if this happends, it is because the dimensions were chosen incorrectly. 
  Either to many processors for the number of points in one dimenison (could 
  not do at least 1 point per processor), or the methods above could not make 
  a distribution of pencils that fit in the cubiods, which would happen if the 
  user gave numbers that wouldent work (we require the number of processors in 
  each dimension of the cuboid must be modulo the number of points in that 
  dimension, otherwise, this error will happen).
*/
  MPI_Cart_create(comm, 
		  ndim, 
		  d->process_topology_2_x.nproc, 
		  period, 
		  0, 
		  &d->process_topology_2_x.cart);
  Coord_x_pencils(self, d->process_topology_2_x.self, d);

  if(self == 0) {
    printf("distribution 2x: [%d:%d:%d]\n",
	   d->process_topology_2_x.nproc[0],
	   d->process_topology_2_x.nproc[1],
	   d->process_topology_2_x.nproc[2]);
    fflush(stdout);
  }
  


  // set up process grid with 2d decomposition (y_PENCILs)
  d->process_topology_2_y.nproc[0] = 0;
  d->process_topology_2_y.nproc[1] = 1; // don't distribute outer dimension
  d->process_topology_2_y.nproc[2] = 0;
  period[0] = period[1] = period[2] = 1;
  MPI_Dims_create(nproc, ndim, d->process_topology_2_y.nproc);
  d->process_topology_2_y.n[0] = n[0] / d->process_topology_2_y.nproc[0];
  d->process_topology_2_y.n[1] = n[1] / d->process_topology_2_y.nproc[1];
  d->process_topology_2_y.n[2] = n[2] / d->process_topology_2_y.nproc[2];
  //variable used to ensure that pencils created fit inside the cuboids, 
  //if not the code will assert out.
  bool check_y_dims=false;
  if(d->process_topology_2_y.n[0] != 0 
     && d->process_topology_2_y.n[1] != 0 
     && d->process_topology_2_y.n[2] != 0)
  {// protects from dividing by zero.
    check_y_dims = (((d->process_topology_3.n[2]) % (d->process_topology_2_y.n[2]) == 0) 
		    && ((d->process_topology_3.n[0]) % (d->process_topology_2_y.n[0]) == 0) 
		    && (n[0] % (d->process_topology_2_y.nproc[2]) == 0) 
		    && (n[0] % (d->process_topology_2_y.nproc[0]) == 0));
    if(self==0 && debug && !check_y_dims)
      fprintf(stderr,"Need to fix Y PENCILS y_procs(%d,%d,%d) 3d.ns(%d,%d,%d) 2d_y.ns(%d,%d,%d)...\n",
	      d->process_topology_2_y.nproc[0],
	      d->process_topology_2_y.nproc[1],
	      d->process_topology_2_y.nproc[2],
	      d->process_topology_3.n[0],
	      d->process_topology_3.n[1],
	      d->process_topology_3.n[2],
	      d->process_topology_2_y.n[0],
	      d->process_topology_2_y.n[1],
	      d->process_topology_2_y.n[2]);
    //try swaping pencil dimensions if the current dimension of the pencils 
    //does not fit inside the cubes.
    if(!(check_y_dims) 
       && ((d->process_topology_3.n[2]) % (d->process_topology_2_y.n[0]) == 0) 
       && ((d->process_topology_3.n[0]) % (d->process_topology_2_y.n[2]) == 0))
    {
      if(self==0 && debug)
	fprintf(stderr,"Swaping Y pencils in initialization .... \n");
      
      int temp = d->process_topology_2_y.nproc[0];
      d->process_topology_2_y.nproc[0] = d->process_topology_2_y.nproc[2];
      d->process_topology_2_y.nproc[2] = temp;
      d->process_topology_2_y.nproc[1] = d->process_topology_2_y.nproc[1];
      
      d->process_topology_2_y.n[0] = n[0] / d->process_topology_2_y.nproc[0];
      d->process_topology_2_y.n[1] = n[1] / d->process_topology_2_y.nproc[1];
      d->process_topology_2_y.n[2] = n[2] / d->process_topology_2_y.nproc[2];
      check_y_dims = (((d->process_topology_3.n[2]) % (d->process_topology_2_y.n[2]) == 0) 
		      && ((d->process_topology_3.n[0]) % (d->process_topology_2_y.n[0]) == 0) 
		      && (n[0] % (d->process_topology_2_y.nproc[2]) == 0) 
		      && (n[0] % (d->process_topology_2_y.nproc[0]) == 0));
    }
  } else {
    check_y_dims = false;
  }
/*
  if that did not work, make a pencil that does by taking the cuboid 
  (np1,np2,np3) and making pencils of the form (np1,1,np3*np2) or 
  (np1*np2,1,np3) depending on the most even distribution it can.
*/
  if(!check_y_dims){
    if(self==0 && debug)
      fprintf(stderr,"MAKING Y PENCILS FIT yprocs(%d,%d,%d) y.ns(%d,%d,%d)...\n", 
	      d->process_topology_2_y.nproc[0],
	      d->process_topology_2_y.nproc[1],
	      d->process_topology_2_y.nproc[2],
	      d->process_topology_2_y.n[0],
	      d->process_topology_2_y.n[1],
	      d->process_topology_2_y.n[2]);
    
    d->process_topology_2_y.nproc[1]=1;
    if(d->process_topology_3.nproc[2] > d->process_topology_3.nproc[0])
    {
      d->process_topology_2_y.nproc[0] = d->process_topology_3.nproc[0]*d->process_topology_3.nproc[1];
      d->process_topology_2_y.nproc[2] = d->process_topology_3.nproc[2];
      if((n[0] % (d->process_topology_2_y.nproc[2]) != 0) 
	 || (n[0] % (d->process_topology_2_y.nproc[0]) != 0))
      {
	d->process_topology_2_y.nproc[2] = d->process_topology_3.nproc[2]*d->process_topology_3.nproc[1];
	d->process_topology_2_y.nproc[0] = d->process_topology_3.nproc[0];
      }
    } else {
      d->process_topology_2_y.nproc[2] = d->process_topology_3.nproc[2]*d->process_topology_3.nproc[1];
      d->process_topology_2_y.nproc[0] = d->process_topology_3.nproc[0];
      if((n[0] % (d->process_topology_2_y.nproc[2]) != 0) 
	 || (n[0] % (d->process_topology_2_y.nproc[0]) != 0))
      {
	d->process_topology_2_y.nproc[0] = d->process_topology_3.nproc[0]*d->process_topology_3.nproc[1];
	d->process_topology_2_y.nproc[2] = d->process_topology_3.nproc[2];
      }
    }
    
    d->process_topology_2_y.n[0] = n[0] / d->process_topology_2_y.nproc[0];
    d->process_topology_2_y.n[1] = n[1] / d->process_topology_2_y.nproc[1];
    d->process_topology_2_y.n[2] = n[2] / d->process_topology_2_y.nproc[2];
    if(self==0 && debug)
      fprintf(stderr,"MAKING Y PENCILS FIT AFTER yprocs(%d,%d,%d) y.ns(%d,%d,%d)...\n",
	      d->process_topology_2_y.nproc[0],
	      d->process_topology_2_y.nproc[1],
	      d->process_topology_2_y.nproc[2],
	      d->process_topology_2_y.n[0],
	      d->process_topology_2_y.n[1],
	      d->process_topology_2_y.n[2]);
    if(d->process_topology_2_y.n[0] != 0 && d->process_topology_2_y.n[1] != 0 
       && d->process_topology_2_y.n[2] != 0)
    {// protects from dividing by zero.
      check_y_dims = (((d->process_topology_3.n[2]) % (d->process_topology_2_y.n[2]) == 0) 
		      && ((d->process_topology_3.n[0]) % (d->process_topology_2_y.n[0]) == 0) 
		      && (n[0] % (d->process_topology_2_y.nproc[2]) == 0) 
		      && (n[0] % (d->process_topology_2_y.nproc[0]) == 0));
    } else {
      check_y_dims=false;
    }
  }
   
  if (d->debug && 0 == self) {
    fprintf(stderr, "  2d_y: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", 
	      d->process_topology_2_y.nproc[i], 
	      separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }
  if(!check_y_dims && debug && (self==0)){
    FILE * outfile;
    outfile = fopen("error.data","a");
    fprintf(outfile,"Y DIMS FAILS:(%d,%d,%d) (%d,%d,%d) \n",
	    d->process_topology_2_y.nproc[0],
	    d->process_topology_2_y.nproc[1],
	    d->process_topology_2_y.nproc[2], 
	    d->process_topology_3.nproc[0],
	    d->process_topology_3.nproc[1],
	    d->process_topology_3.nproc[2]);
  }
  assert(check_y_dims);
/*
  if this happends, it is because the dimensions were chosen incorrectly. 
  Either to many processors for the number of points in one dimenison (could 
  not do at least 1 point per processor), or the methods above could 
  not make a distribution of pencils that fit in the cubiods, which would 
  happen if the user gave numbers that wouldent work (we require the number of 
  processors in each dimension of the cuboid must be modulo the number of 
  points in that dimension, otherwise, this error will happen).
*/
  MPI_Cart_create(comm, 
		  ndim, 
		  d->process_topology_2_y.nproc, 
		  period, 
		  0, 
		  &d->process_topology_2_y.cart);
  //find the cartesian coord of the current rank (for the y_pencil)
  Coord_y_pencils(self,d->process_topology_2_y.self,d);

  if(self == 0) {
    printf("distribution 2y: [%d:%d:%d]\n",
	   d->process_topology_2_y.nproc[0],
	   d->process_topology_2_y.nproc[1],
	   d->process_topology_2_y.nproc[2]);
    fflush(stdout);
  }


  
  if (d->debug) {
    int myrank_cube;
    Rank_cube(&myrank_cube,d->process_topology_3.self,d);
    int myrank_x;
    Rank_x_pencils(&myrank_x,d->process_topology_2_x.self,d);
    int myrank_y;
    Rank_y_pencils(&myrank_y,d->process_topology_2_y.self,d);
    int myrank_z;
    Rank_z_pencils(&myrank_z,d->process_topology_2_z.self,d);
    if(myrank_z != self 
       || myrank_y != self 
       || myrank_x != self 
       || myrank_cube != self)
      abort(); //means ranks were calculated wrong.
    if (0 == self) {
      fprintf(stderr, "Process map:\n");
    }
    for (int p = 0; p < nproc; ++p) {
      MPI_Barrier(comm);
      if (p == self) {
	fprintf(stderr, "  %d: 1d = (%d, %d, %d), 2d_x = (%d, %d, %d) rank is= %d,2d_y = (%d, %d, %d) rank is= %d,2d_z = (%d, %d, %d) rank is= %d, 3d = (%d, %d, %d). rank is= %d\n",
		self,
		d->process_topology_1.self[0], 
		d->process_topology_1.self[1], 
		d->process_topology_1.self[2],
		d->process_topology_2_x.self[0], 
		d->process_topology_2_x.self[1], 
		d->process_topology_2_x.self[2],
		myrank_x,
		d->process_topology_2_y.self[0], 
		d->process_topology_2_y.self[1], 
		d->process_topology_2_y.self[2],
		myrank_y,
		d->process_topology_2_z.self[0], 
		d->process_topology_2_z.self[1], 
		d->process_topology_2_z.self[2],
		myrank_z,
		d->process_topology_3.self[0], 
		d->process_topology_3.self[1], 
		d->process_topology_3.self[2],
		myrank_cube);
      }
    }
  }

  //allocate size of buffers used to hold pencil chunks of data in the 
  //distribution routines for 3d to 1d and vica versa.
  int buff_z_chunk = d->process_topology_2_z.n[0]*d->process_topology_2_z.n[1]*d->process_topology_3.n[2];
  int buff_y_chunk = d->process_topology_2_y.n[0]*d->process_topology_2_y.n[2]*d->process_topology_3.n[1];
  int buff_x_chunk = d->process_topology_2_x.n[1]*d->process_topology_2_x.n[2]*d->process_topology_3.n[0];
  int buff_size = 0;
  if(buff_z_chunk > buff_y_chunk){
    buff_size=buff_z_chunk;
  } else {
    buff_size=buff_y_chunk;
  }
  if(buff_x_chunk > buff_size)
    buff_size = buff_x_chunk;
  
  d->d2_chunk=(complex_t *) malloc(sizeof(complex_t)*buff_size);
  d->d3_chunk=(complex_t *) malloc(sizeof(complex_t)*buff_size);
}

// Use MPI_Dims_create to create a 3D decomposition, or use a user-provided decomposition
// if it is appropriate (it has the required number of processors)
void Custom3D_Dims_create(const int Ndims[], int nproc, int ndims, int dims[])
{

  int check = 1;
  for(int i=0; i<ndims; i++) check *= Ndims[i];

  if(check == nproc) {
    for(int i=0; i<ndims; i++) dims[i] = Ndims[i];
  }
  else {
    MPI_Dims_create(nproc, ndims, dims);
  } 

}

// create 1-, 2- and 3-d cartesian data distributions with explicitly
// provided dimension lists
void distribution_init_explicit(MPI_Comm comm, 
				const int n[], 
                                int nproc_1d[],
                                int nproc_2d_x[],
                                int nproc_2d_y[],
                                int nproc_2d_z[],
                                int nproc_3d[],
                                distribution_t *d, 
				bool debug)
{
  d->parent = comm;

  int nproc;
  int self;
  int ndim = 3;
  int period[3];
  
  MPI_Comm_rank(comm, &self);
  MPI_Comm_size(comm, &nproc);
  
  if (!self) printf("Initializing redistribution using a %s layout on %d ranks.\n",
#ifdef PENCIL
		    "pencil"
#else
		    "slab"
#endif
		    ,nproc);
  
  d->debug = debug;
  for (int i = 0; i < 3; ++i)
    d->n[i] = n[i];
  
  // check supplied dimension lists are valid
  assert(nproc_1d[0] == nproc);
  assert(nproc_1d[1] == 1);
  assert(nproc_1d[2] == 1);
  
  assert(nproc_2d_x[1] * nproc_2d_x[2] == nproc);
  assert(nproc_2d_x[0] == 1);
  
  assert(nproc_2d_y[0] * nproc_2d_y[2] == nproc);
  assert(nproc_2d_y[1] == 1);
  
  assert(nproc_2d_z[0] * nproc_2d_z[1] == nproc);
  assert(nproc_2d_z[2] == 1);
  
  assert(nproc_3d[0] * nproc_3d[1] * nproc_3d[2]== nproc);
  
  // set up process grid with 1d decomposition (SLABs)
  period[0] = period[1] = period[2] = 1;
  MPI_Cart_create(comm, ndim, nproc_1d, period, 0, &d->process_topology_1.cart);
  MPI_Cart_get(d->process_topology_1.cart, ndim, d->process_topology_1.nproc, d->process_topology_1.period, d->process_topology_1.self);
  d->process_topology_1.n[0] = n[0] / d->process_topology_1.nproc[0];
  d->process_topology_1.n[1] = n[1] / d->process_topology_1.nproc[1];
  d->process_topology_1.n[2] = n[2] / d->process_topology_1.nproc[2];
  if (d->debug && 0 == self) {
    fprintf(stderr, "Process grids:\n");
    fprintf(stderr, "  1d: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", d->process_topology_1.nproc[i], separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }
  
  // set up process grid with 3d decomposition (CUBE)
  period[0] = period[1] = period[2] = 1;
  MPI_Cart_create(comm, ndim, nproc_3d, period, 0, &d->process_topology_3.cart);
  Coord_cube(self,d->process_topology_3.self,d);
  if (d->debug && 0 == self) {
    fprintf(stderr, "  3d: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", d->process_topology_3.nproc[i], separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }
  if(debug){
    int prev_coord[3];
    prev_coord[0]=d->process_topology_3.self[0];
    prev_coord[1]=d->process_topology_3.self[1];
    prev_coord[2]=d->process_topology_3.self[2];
    MPI_Cart_get(d->process_topology_3.cart, ndim, d->process_topology_3.nproc, d->process_topology_3.period, d->process_topology_3.self);
    for(int i=0; i < 3; i++){
      if(prev_coord[i] != d->process_topology_3.self[i])abort();//Cube coordinates calculated wrong!
    }
  }
  d->process_topology_3.n[0] = n[0] / d->process_topology_3.nproc[0];
  d->process_topology_3.n[1] = n[1] / d->process_topology_3.nproc[1];
  d->process_topology_3.n[2] = n[2] / d->process_topology_3.nproc[2];
  
  // set up process grid with 2d_x decomposition (X dim Pencils)
  period[0] = period[1] = period[2] = 1;
  MPI_Cart_create(comm, ndim, nproc_2d_x, period, 0, &d->process_topology_2_x.cart);
  d->process_topology_2_x.nproc[0]=nproc_2d_x[0];
  d->process_topology_2_x.nproc[1]=nproc_2d_x[1];
  d->process_topology_2_x.nproc[2]=nproc_2d_x[2];
  d->process_topology_2_x.n[0] = n[0] / d->process_topology_2_x.nproc[0];
  d->process_topology_2_x.n[1] = n[1] / d->process_topology_2_x.nproc[1];
  d->process_topology_2_x.n[2] = n[2] / d->process_topology_2_x.nproc[2];
  
  bool check_x_dims=((d->process_topology_3.n[2]) % (d->process_topology_2_x.n[2]) == 0) && ((d->process_topology_3.n[1]) % (d->process_topology_2_x.n[1]) == 0) && (n[0] % (d->process_topology_2_x.nproc[2]) == 0) && (n[0] % (d->process_topology_2_x.nproc[0]) == 0);
  if(!check_x_dims && debug && (self==0)){
    FILE * outfile;
    outfile= fopen("error.data","a");
    fprintf(outfile,"X DIMS FAILS:(%d,%d,%d) (%d,%d,%d) \n",d->process_topology_2_x.nproc[0],d->process_topology_2_x.nproc[1],d->process_topology_2_x.nproc[2], d->process_topology_3.nproc[0],d->process_topology_3.nproc[1],d->process_topology_3.nproc[2]);
  }
  assert(check_x_dims);//if this happends, it is because the dimensions were chosen incorrectly. Either to many processors for the number of points in one dimenison (could not do at least 1 point per processor), or the methods above could 
  //not make a distribution of pencils that fit in the cubiods, which would happen if the user gave numbers that wouldent work (we require the number of processors in each dimension of the cuboid must be modulo the number of points 
  //in that dimension, otherwise, this error will happen).
  Coord_x_pencils(self,d->process_topology_2_x.self,d);
  
  if (d->debug && 0 == self) {
    fprintf(stderr, "  2d_x: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", d->process_topology_2_x.nproc[i], separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }
  
  // set up process grid with 2d_y decomposition (Y dim Pencils)
  period[0] = period[1] = period[2] = 1;
  MPI_Cart_create(comm, ndim, nproc_2d_y, period, 0, &d->process_topology_2_y.cart);
  d->process_topology_2_y.nproc[0]=nproc_2d_y[0];
  d->process_topology_2_y.nproc[1]=nproc_2d_y[1];
  d->process_topology_2_y.nproc[2]=nproc_2d_y[2];
  d->process_topology_2_y.n[0] = n[0] / d->process_topology_2_y.nproc[0];
  d->process_topology_2_y.n[1] = n[1] / d->process_topology_2_y.nproc[1];
  d->process_topology_2_y.n[2] = n[2] / d->process_topology_2_y.nproc[2];
  
  
  bool check_y_dims=(((d->process_topology_3.n[2]) % (d->process_topology_2_y.n[2]) == 0) && ((d->process_topology_3.n[0]) % (d->process_topology_2_y.n[0]) == 0) && (n[0] % (d->process_topology_2_y.nproc[2]) == 0) && (n[0] % (d->process_topology_2_y.nproc[0]) == 0));
  if(!check_y_dims && debug && (self==0)){
    FILE * outfile;
    outfile= fopen("error.data","a");
    fprintf(outfile,"Y DIMS FAILS:(%d,%d,%d) (%d,%d,%d) \n",d->process_topology_2_y.nproc[0],d->process_topology_2_y.nproc[1],d->process_topology_2_y.nproc[2], d->process_topology_3.nproc[0],d->process_topology_3.nproc[1],d->process_topology_3.nproc[2]);
  }
  assert(check_y_dims);//if this happends, it is because the dimensions were chosen incorrectly. Either to many processors for the number of points in one dimenison (could not do at least 1 point per processor), or the methods above could 
  //not make a distribution of pencils that fit in the cubiods, which would happen if the user gave numbers that wouldent work (we require the number of processors in each dimension of the cuboid must be modulo the number of points 
  //in that dimension, otherwise, this error will happen).
  Coord_y_pencils(self,d->process_topology_2_y.self,d);
  
  if (d->debug && 0 == self) {
    fprintf(stderr, "  2d_y: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", d->process_topology_2_y.nproc[i], separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }
  
  // set up process grid with 2d_z decomposition (Z dim pencils)
  period[0] = period[1] = period[2] = 1;
  MPI_Cart_create(comm, ndim, nproc_2d_z, period, 0, &d->process_topology_2_z.cart);
  d->process_topology_2_z.nproc[0]=nproc_2d_z[0];
  d->process_topology_2_z.nproc[1]=nproc_2d_z[1];
  d->process_topology_2_z.nproc[2]=nproc_2d_z[2];
  d->process_topology_2_z.n[0] = n[0] / d->process_topology_2_z.nproc[0];
  d->process_topology_2_z.n[1] = n[1] / d->process_topology_2_z.nproc[1];
  d->process_topology_2_z.n[2] = n[2] / d->process_topology_2_z.nproc[2];
  
  
  bool check_z_dims=((d->process_topology_3.n[0]) % (d->process_topology_2_z.n[0]) == 0) && ((d->process_topology_3.n[1]) % (d->process_topology_2_z.n[1]) == 0) && (n[0] % (d->process_topology_2_z.nproc[0]) == 0) && (n[0] % (d->process_topology_2_z.nproc[1]) == 0);
  if(!check_z_dims && debug && (self==0)){
    FILE * outfile;
    outfile= fopen("error.data","a");
    fprintf(outfile,"Z DIMS FAILS:(%d,%d,%d) (%d,%d,%d) \n",d->process_topology_2_z.nproc[0],d->process_topology_2_z.nproc[1],d->process_topology_2_z.nproc[2], d->process_topology_3.nproc[0],d->process_topology_3.nproc[1],d->process_topology_3.nproc[2]);
  }
  assert(check_z_dims);//if this happends, it is because the dimensions were chosen incorrectly. Either to many processors for the number of points in one dimenison (could not do at least 1 point per processor), or the methods above could 
  //not make a distribution of pencils that fit in the cubiods, which would happen if the user gave numbers that wouldent work (we require the number of processors in each dimension of the cuboid must be modulo the number of points 
  //in that dimension, otherwise, this error will happen).
  Coord_z_pencils(self,d->process_topology_2_z.self,d);
  
  if (d->debug && 0 == self) {
    fprintf(stderr, "  2d_z: ");
    for (int i = 0; i < ndim; ++i) {
      fprintf(stderr, "%d%s", d->process_topology_2_z.nproc[i], separator(i, ndim));
    }
    fprintf(stderr, "\n");
  }
  //assert that all pencils fit in the cuboid.
  
  if (d->debug) {
    int myrank_cube;
    Rank_cube(&myrank_cube,d->process_topology_3.self,d);
    int myrank_z;
    Rank_z_pencils(&myrank_z,d->process_topology_2_z.self,d);
    int myrank_y;
    Rank_y_pencils(&myrank_y,d->process_topology_2_y.self,d);
    int myrank_x;
    Rank_x_pencils(&myrank_x,d->process_topology_2_x.self,d);
    if(myrank_z != self || myrank_y != self || myrank_x != self || myrank_cube != self)abort(); //means ranks were calculated wrong.
    if (0 == self) {
      fprintf(stderr, "Process map:\n");
    }
    for (int p = 0; p < nproc; ++p) {
      MPI_Barrier(comm);
      if (p == self) {
	fprintf(stderr,
		"  %d: 1d = (%d, %d, %d), 2d_x = (%d, %d, %d) rank (%d), 2d_y = (%d, %d, %d) rank (%d), 2d_z = (%d, %d, %d) rank (%d), 3d = (%d, %d, %d) rank (%d).\n",
		self,
		d->process_topology_1.self[0], d->process_topology_1.self[1], d->process_topology_1.self[2],
		d->process_topology_2_x.self[0], d->process_topology_2_x.self[1], d->process_topology_2_x.self[2],myrank_x,
		d->process_topology_2_y.self[0], d->process_topology_2_y.self[1], d->process_topology_2_y.self[2],myrank_y,
		d->process_topology_2_z.self[0], d->process_topology_2_z.self[1], d->process_topology_2_z.self[2],myrank_z,
		d->process_topology_3.self[0], d->process_topology_3.self[1], d->process_topology_3.self[2],myrank_cube);
      }
    }
  }
}




///
// clean up the data distribution
//   d    distribution descriptor
///
void distribution_fini(distribution_t *d)
{
  MPI_Comm_free(&d->process_topology_1.cart);
  MPI_Comm_free(&d->process_topology_2_x.cart);
  MPI_Comm_free(&d->process_topology_2_y.cart);
  MPI_Comm_free(&d->process_topology_2_z.cart);
  MPI_Comm_free(&d->process_topology_3.cart);
  free(d->d2_chunk);
  free(d->d3_chunk);
  free(d->gridmap);
  free(d->rankmap);
}


///
// check that the dimensions, n, of an array are commensurate with the
// process grids of this distribution
//   n    (global) grid dimensions
//   d    distribution descriptor
///
void distribution_assert_commensurate(distribution_t *d)
{
  for (int i = 0; i < 3; ++i) {
#if defined(PENCIL)
    assert(0 == (d->n[i] % d->process_topology_2_x.nproc[i]));
    assert(0 == (d->n[i] % d->process_topology_2_y.nproc[i]));
    assert(0 == (d->n[i] % d->process_topology_2_z.nproc[i]));
#else
    assert(0 == (d->n[i] % d->process_topology_1.nproc[i]));
#endif
    assert(0 == (d->n[i] % d->process_topology_3.nproc[i]));
  }
}


// forward declarations
static void redistribute(const complex_t *, complex_t *, distribution_t *, int);
static void redistribute_2_and_3(const complex_t *, complex_t *, distribution_t *, int, int);
static void redistribute_slab(const complex_t *, complex_t *, distribution_t *, int);


///
// redistribute a 1-d to a 3-d data distribution
//   a    input
//   b    ouput
//   d    distribution descriptor
///
void distribution_1_to_3(const complex_t *a,
                         complex_t *b,
                         distribution_t *d)
{
  if (USE_SLAB_WORKAROUND) {
    redistribute_slab(a, b, d, REDISTRIBUTE_1_TO_3);
  } else {
    redistribute(a, b, d, REDISTRIBUTE_1_TO_3);
  }
}


///
// redistribute a 3-d to a 1-d data distribution
//   a    input
//   b    ouput
//   d    distribution descriptor
///
void distribution_3_to_1(const complex_t *a,
                         complex_t *b,
                         distribution_t *d)
{
  if (USE_SLAB_WORKAROUND) {
    redistribute_slab(a, b, d, REDISTRIBUTE_3_TO_1);
  } else {
    redistribute(a, b, d, REDISTRIBUTE_3_TO_1);
  }
}


///
// redistribute between 1- and 3-d distributions.
//   a    input
//   b    ouput
//   d    distribution descriptor
//   dir  direction of redistribution
//
// This actually does the work.
///
static void redistribute(const complex_t *a,
                         complex_t *b,
                         distribution_t *d,
                         int direction)
{
  int remaining_dim[3];
  MPI_Comm subgrid_cart;
  int subgrid_self;
  int subgrid_nproc;
  
  // exchange data with processes in a 2-d slab of 3-d subdomains
  
  remaining_dim[0] = 0;
  remaining_dim[1] = 1;
  remaining_dim[2] = 1;
  MPI_Cart_sub(d->process_topology_3.cart, remaining_dim, &subgrid_cart);
  MPI_Comm_rank(subgrid_cart, &subgrid_self);
  MPI_Comm_size(subgrid_cart, &subgrid_nproc);
  
  for (int p = 0; p < subgrid_nproc; ++p) {
    int d1_peer = (subgrid_self + p) % subgrid_nproc;
    int d3_peer = (subgrid_self - p + subgrid_nproc) % subgrid_nproc;
    int coord[2];
    int sizes[3];
    int subsizes[3];
    int starts[3];
    MPI_Datatype d1_type;
    MPI_Datatype d3_type;
    
    MPI_Cart_coords(subgrid_cart, d1_peer, 2, coord);
    if (0) {
      int self;
      MPI_Comm_rank(MPI_COMM_WORLD, &self);
      fprintf(stderr, "%d: d1_peer, d1_coord, d3_peer = %d, (%d, %d), %d\n",
	      self, d1_peer, coord[0], coord[1], d3_peer);
    }
    
    // create dataypes representing a subarray in the 1- and 3-d distributions
    
    sizes[0] = d->process_topology_1.n[0];
    sizes[1] = d->process_topology_1.n[1];
    sizes[2] = d->process_topology_1.n[2];
    subsizes[0] = d->process_topology_1.n[0];
    subsizes[1] = d->process_topology_3.n[1];
    subsizes[2] = d->process_topology_3.n[2];
    starts[0] = 0;
    starts[1] = coord[0] * d->process_topology_3.n[1];
    starts[2] = coord[1] * d->process_topology_3.n[2];
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE_COMPLEX, &d1_type);
    MPI_Type_commit(&d1_type);
    
    sizes[0] = d->process_topology_3.n[0];
    sizes[1] = d->process_topology_3.n[1];
    sizes[2] = d->process_topology_3.n[2];
    subsizes[0] = d->process_topology_1.n[0];
    subsizes[1] = d->process_topology_3.n[1];
    subsizes[2] = d->process_topology_3.n[2];
    starts[0] = d3_peer * d->process_topology_1.n[0];
    starts[1] = 0;
    starts[2] = 0;
    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE_COMPLEX, &d3_type);
    MPI_Type_commit(&d3_type);
    
    // exchange data
    
    if (direction == REDISTRIBUTE_3_TO_1) {
      MPI_Sendrecv((char *) a, 1, d3_type, d3_peer, 0,
		   (char *) b, 1, d1_type, d1_peer, 0,
		   subgrid_cart, MPI_STATUS_IGNORE);
    } else if (direction == REDISTRIBUTE_1_TO_3) {
      MPI_Sendrecv((char *) a, 1, d1_type, d1_peer, 0,
		   (char *) b, 1, d3_type, d3_peer, 0,
		   subgrid_cart, MPI_STATUS_IGNORE);
    } else {
      abort();
    }
    
    // free datatypes
    
    MPI_Type_free(&d1_type);
    MPI_Type_free(&d3_type);
  }
  
  MPI_Comm_free(&subgrid_cart);
}


///
// redistribute between 1- and 3-d distributions.
//   a    input
//   b    ouput
//   d    distribution descriptor
//   dir  direction of redistribution
//
// This actually does the work, using slabs of subarrays to work
// around an issue in Open MPI with large non-contiguous datatypes.
///
static void redistribute_slab(const complex_t *a,
                              complex_t *b,
                              distribution_t *d,
                              int direction)
{
  int remaining_dim[3];
  MPI_Comm subgrid_cart;
  int subgrid_self;
  int subgrid_nproc;
  ptrdiff_t d1_slice = d->process_topology_1.n[1] * d->process_topology_1.n[2] * sizeof(complex_t);
  ptrdiff_t d3_slice = d->process_topology_3.n[1] * d->process_topology_3.n[2] * sizeof(complex_t);
  
  // exchange data with processes in a 2-d slab of 3-d subdomains
  
  remaining_dim[0] = 0;
  remaining_dim[1] = 1;
  remaining_dim[2] = 1;
  MPI_Cart_sub(d->process_topology_3.cart, remaining_dim, &subgrid_cart);
  MPI_Comm_rank(subgrid_cart, &subgrid_self);
  MPI_Comm_size(subgrid_cart, &subgrid_nproc);
  
  for (int p = 0; p < subgrid_nproc; ++p) {
    int coord[2];
    int d1_peer = (subgrid_self + p) % subgrid_nproc;
    int d3_peer = (subgrid_self - p + subgrid_nproc) % subgrid_nproc;
    
    MPI_Cart_coords(subgrid_cart, d1_peer, 2, coord);
    if (0) {
      int self;
      MPI_Comm_rank(MPI_COMM_WORLD, &self);
      fprintf(stderr, "%d: d1_peer, d1_coord, d3_peer = %d, (%d, %d), %d\n",
	      self, d1_peer, coord[0], coord[1], d3_peer);
    }
    
    for (int slice = 0; slice < d->process_topology_1.n[0]; ++slice) {
      int sizes[2];
      int subsizes[2];
      int starts[2];
      MPI_Datatype d1_type;
      MPI_Datatype d3_type;
      ptrdiff_t d1_offset = slice * d1_slice;
      ptrdiff_t d3_offset = (slice + d3_peer * d->process_topology_1.n[0]) * d3_slice;
      
      // create subarray dataypes representing the slice subarray in the 1- and 3-d distributions
      
      sizes[0] = d->process_topology_1.n[1];
      sizes[1] = d->process_topology_1.n[2];
      subsizes[0] = d->process_topology_3.n[1];
      subsizes[1] = d->process_topology_3.n[2];
      starts[0] = coord[0] * d->process_topology_3.n[1];
      starts[1] = coord[1] * d->process_topology_3.n[2];
      MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE_COMPLEX, &d1_type);
      MPI_Type_commit(&d1_type);
      
      MPI_Type_contiguous(d->process_topology_3.n[1] * d->process_topology_3.n[2],
			  MPI_DOUBLE_COMPLEX,
			  &d3_type);
      MPI_Type_commit(&d3_type);
      
      // exchange data
      
      if (direction == REDISTRIBUTE_3_TO_1) {
	MPI_Sendrecv((char *) a + d3_offset, 1, d3_type, d3_peer, 0,
		     (char *) b + d1_offset, 1, d1_type, d1_peer, 0,
		     subgrid_cart, MPI_STATUS_IGNORE);
      } else if (direction == REDISTRIBUTE_1_TO_3) {
	MPI_Sendrecv((char *) a + d1_offset, 1, d1_type, d1_peer, 0,
		     (char *) b + d3_offset, 1, d3_type, d3_peer, 0,
		     subgrid_cart, MPI_STATUS_IGNORE);
      } else {
	abort();
      }
      
      // free datatypes
      
      MPI_Type_free(&d1_type);
      MPI_Type_free(&d3_type);
    }
  }
  
  MPI_Comm_free(&subgrid_cart);
}


///
// redistribute a 2-d to a 3-d data distribution
//   a    input
//   b    ouput
//   d    distribution descriptor
///
void distribution_2_to_3(const complex_t *a,
                         complex_t *b,
                         distribution_t *d, int z_dim)
{
  redistribute_2_and_3(a, b, d, REDISTRIBUTE_2_TO_3, z_dim);
}


///
// redistribute a 3-d to a 2-d data distribution
//   a    input
//   b    ouput
//   d    distribution descriptor
///
void distribution_3_to_2(const complex_t *a,
                         complex_t *b,
                         distribution_t *d, int z_dim)
{
  redistribute_2_and_3(a, b, d, REDISTRIBUTE_3_TO_2, z_dim);
}


///
// redistribute between 2- and 3-d distributions.
//   a    input
//   b    ouput
//   d    distribution descriptor
//   dir  direction of redistribution
//
// This actually does the work.
///
static void redistribute_2_and_3(const complex_t *a,
                                 complex_t *b,
                                 distribution_t *d,
                                 int direction,
				 int z_dim)
{
  int self = d->process_topology_1.self[0];
  int npeers;
  int me=0;//determines which processor to print
  bool print_me=false; //prints info on proccessor whose rank = me.
  bool print_mess=false;//prints communication sends and recieves without actually doing the comms(intended to debug comm hangs).
  bool print_result=false /*true*/;//prints a line in a file called "passed.data" which happends if the code runs completly.
  assert(z_dim==0||z_dim==1||z_dim==2);
  int x_dim=0,y_dim=0;
  //x_dim, y_dim and z_dim are the dimensions of the x,y,z axis of the pencil with respect to the original axis(where index 2 is into the grid, 1 is vertical translation and 0 is horizontal).
  switch(z_dim){
    case 0: x_dim=1; y_dim=2; 
      if((self == me) && print_me)fprintf(stderr, "DOING X PENCILS!...\n"); break;
    case 1: x_dim=2; y_dim=0;
      if((self == me && print_me))fprintf(stderr, "DOING Y PENCILS!...\n"); break;
    case 2: x_dim=0; y_dim=1;
      if((self == me && print_me))fprintf(stderr, "DOING Z PENCILS!...\n"); break;
    default: assert("incorrect inputed dimension");
  }
  
  // assuming dimensions are all commensurate, then the number of
  // peers to exchange with is the number of processes in the z_dimension
  // direction in the 3d distribution
  npeers = d->process_topology_3.nproc[z_dim]; //picked last direction (lets say into the grid)
  
  // book-keeping for the processor translation in the x-y plane
  int p0 = 0;
  int p1 = 0;
  int p1max = 0;
  
  MPI_Request req1=MPI_REQUEST_NULL;
  MPI_Request req2=MPI_REQUEST_NULL;
  
  int pencil_sizes[3];
  int cube_sizes[3];
  int subsizes[3];
  
  
  cube_sizes[x_dim] = d->process_topology_3.n[x_dim];
  cube_sizes[y_dim] = d->process_topology_3.n[y_dim]; 
  cube_sizes[z_dim] = d->process_topology_3.n[z_dim];
  
  //set varibles used to calculate the subarrays of each pencil and cube.
  switch(z_dim){
    case 0: 
      p1max = d->process_topology_2_x.nproc[x_dim] / d->process_topology_3.nproc[x_dim] - 1; 
      //find out the size of the chunk you need to use (stored in subsizes), and set sizes to the local size of the pencil.
      //The x and y dimensions of the subchunck will be the dimensions of the pencil (since the code asserts at the beginning that all pencils fit inside the 3d cuboid.)
      //The z dimension will be the dimension of the cuboid, since this will always be <= to the z_dim of the pencil.
      pencil_sizes[x_dim] = d->process_topology_2_x.n[x_dim];
      pencil_sizes[y_dim] = d->process_topology_2_x.n[y_dim];  
      pencil_sizes[z_dim] = d->process_topology_2_x.n[z_dim]; 
      subsizes[x_dim] = d->process_topology_2_x.n[x_dim];
      subsizes[y_dim] = d->process_topology_2_x.n[y_dim];   
      break;
    case 1: 
      p1max = d->process_topology_2_y.nproc[x_dim] / d->process_topology_3.nproc[x_dim] - 1; 
      pencil_sizes[x_dim] = d->process_topology_2_y.n[x_dim];
      pencil_sizes[y_dim] = d->process_topology_2_y.n[y_dim];  
      pencil_sizes[z_dim] = d->process_topology_2_y.n[z_dim]; 
      subsizes[x_dim] = d->process_topology_2_y.n[x_dim];
      subsizes[y_dim] = d->process_topology_2_y.n[y_dim];   
      break;
    case 2: 
      p1max = d->process_topology_2_z.nproc[y_dim] / d->process_topology_3.nproc[y_dim] - 1; 
      pencil_sizes[x_dim] = d->process_topology_2_z.n[x_dim];
      pencil_sizes[y_dim] = d->process_topology_2_z.n[y_dim];  
      pencil_sizes[z_dim] = d->process_topology_2_z.n[z_dim]; 
      subsizes[x_dim] = d->process_topology_2_z.n[x_dim];
      subsizes[y_dim] = d->process_topology_2_z.n[y_dim];   
      break;
  }
  subsizes[z_dim] = d->process_topology_3.n[z_dim];
  int chunk_size=subsizes[0]*subsizes[1]*subsizes[2];//size of data chunks that will be communicated between pencil and cube distributions.
  
  //set variables that will be used to find pencils chunks
  int pencil_dims[3]={0,0,0};// size of entire pencil in its local coord system 
  int local_sizes[3]={0,0,0}; //size of chunck in its local coord system.
  if(z_dim==2){
    local_sizes[0]=subsizes[0];
    local_sizes[1]=subsizes[1];
    local_sizes[2]=subsizes[2];
    pencil_dims[0]=d->process_topology_2_z.n[0];//pencil dims in grid coord system (where index 2 is the z direction).
    pencil_dims[1]=d->process_topology_2_z.n[1];
    pencil_dims[2]=d->process_topology_2_z.n[2];
  }
  else if(z_dim==1){
    
    local_sizes[0]=subsizes[0];
    local_sizes[1]=subsizes[2];
    local_sizes[2]=subsizes[1];
    pencil_dims[0]=d->process_topology_2_y.n[0];
    pencil_dims[1]=d->process_topology_2_y.n[2];
    pencil_dims[2]=d->process_topology_2_y.n[1];
  }
  else if(z_dim==0){
    local_sizes[0]=subsizes[2];
    local_sizes[1]=subsizes[1];
    local_sizes[2]=subsizes[0];
    pencil_dims[0]=d->process_topology_2_x.n[2];
    pencil_dims[1]=d->process_topology_2_x.n[1];
    pencil_dims[2]=d->process_topology_2_x.n[0];
  }
  
  if((self == me) && print_me)fprintf(stderr, "%d, %d, %d, %d Dimensions!...\n", x_dim,y_dim,z_dim, p1max);
  
  // communicate with our peers
  for (int p = 0; p < npeers; ++p) {
    if((self == me) && print_me)fprintf(stderr, "%d, %d, %d Made it beg-for!...\n", self,p, npeers);
    
    int d2_coord[3];
    int d2_peer;
    int d2_peer_coord[3];
    int d3_coord[3];
    int d3_peer;
    int d3_peer_coord[3];
    int recv_peer;
    int send_peer;
    int d2_array_start[3];
    int d3_array_start[3];
    //turn the processor coordinate into one specified by the number of data points in each dimension.
    for (int i = 0; i < 3; ++i) {
      switch(z_dim){
	case 0: d2_coord[i]  = d->process_topology_2_x.self[i] * d->process_topology_2_x.n[i]; break;
	case 1: d2_coord[i]  = d->process_topology_2_y.self[i] * d->process_topology_2_y.n[i]; break;
	case 2: d2_coord[i]  = d->process_topology_2_z.self[i] * d->process_topology_2_z.n[i]; break;
      }
    }
    //over every iteration of the loop, transverse down the pencil (since it will be divided in chunks whose coordinates will only differ in the z_dimension.
    d2_coord[z_dim] += p * d->process_topology_3.n[z_dim]; 
    
    
    if((self == me) && print_me)fprintf(stderr, "%d, %d, %d Coord!...\n", d2_coord[0],d2_coord[1],d2_coord[2]);
    
    
    //d2_array_start is the starting index of the chunk in the pencils local coordinates.
    d2_array_start[0] = d2_coord[x_dim] % pencil_sizes[x_dim]; 
    d2_array_start[1] = d2_coord[y_dim] % pencil_sizes[y_dim]; 
    d2_array_start[2] = d2_coord[z_dim] % pencil_sizes[z_dim]; 
    
    if (DEBUG_CONDITION || ((self== me) && print_me)) {
      fprintf(stderr,
	      "%d: pencil_sizes=(%d,%d,%d), cube_sizes=(%d,%d,%d), subsizes=(%d,%d,%d),d2_coord=(%d,%d,%d), d2_array_start=(%d,%d,%d) \n",
	      self,
	      pencil_sizes[0], pencil_sizes[1], pencil_sizes[2],
	      cube_sizes[0], cube_sizes[1], cube_sizes[2],
	      subsizes[0], subsizes[1], subsizes[2],
	      d2_coord[0], d2_coord[1], d2_coord[2],
	      d2_array_start[0],d2_array_start[1],d2_array_start[2]);
    }
    
    
    //if making cuboids from pencils, right here we need to fill the d2_chunk array with the data that later needs to be sent to a cuboid.
        //The array is a chunk of the pencil and is why we needed to calculate the starting index for the array in the local coordinates.
    if(direction == REDISTRIBUTE_2_TO_3){	
      int64_t ch_indx=0;
      int dims_size=pencil_dims[0]*pencil_dims[1]*pencil_dims[2];
      for(int i0=d2_array_start[0];i0<d2_array_start[0]+local_sizes[0];i0++){
	for(int i1=d2_array_start[1];i1<d2_array_start[1]+local_sizes[1];i1++){
	  for(int i2=d2_array_start[2];i2<d2_array_start[2]+local_sizes[2];i2++){
	    int64_t local_indx=pencil_dims[2]*(pencil_dims[1]*i0+i1) + i2;
	    assert(local_indx < dims_size);
	    assert(ch_indx <chunk_size && ch_indx >= 0 && local_indx>=0 && local_indx < dims_size);
	    d->d2_chunk[ch_indx]=a[local_indx];
	    ch_indx++;
	  }
	}
      }
      
      if((self == me) && print_me)fprintf(stderr, "%d, %d, %d, pencil_dims!...\n", pencil_dims[0],pencil_dims[1],pencil_dims[2]);
    }
    
    // what peer in the 3d distribution owns this subarray? 
    for (int i = 0; i < 3; ++i) {
      d3_peer_coord[i] = d2_coord[i] / d->process_topology_3.n[i];
    }
    if((self == me) && print_me)fprintf(stderr, "%d, %d, %d Cube that hits pencil coord!...\n",d3_peer_coord[0],d3_peer_coord[1],d3_peer_coord[2]);
    //find the rank of this peer.
    switch(z_dim){
      case 0: MPI_Cart_rank(d->process_topology_3.cart, d3_peer_coord, &d3_peer); break;
      case 1: MPI_Cart_rank(d->process_topology_3.cart, d3_peer_coord, &d3_peer); break;
      case 2: MPI_Cart_rank(d->process_topology_3.cart, d3_peer_coord, &d3_peer); break;
    }
    if((self == me) && print_me)fprintf(stderr, "%d, %d, Made it half way!...\n", self,p);
    if((self == me) && print_me)fprintf(stderr, "%d, %d, PEER!...\n", self,d3_peer);
    
    //By here in the for loop, we have broken the pencil into a chunk and found which cuboid it resides; over every iteration, the for-loop will break up the pencil in the z_dimension.
    //From here on we do the opposite. We divide the cuboid into chunks (that are the same size as the ones in the pencil), and determine which pencils own these chunks.
    
    
    // what is the coordinate of my pth subarray in the 3d distribution?
    for (int i = 0; i < 3; ++i) {
      switch(z_dim){
	case 0: d3_coord[i]  = d->process_topology_3.self[i] * d->process_topology_3.n[i]; break;
	case 1: d3_coord[i]  = d->process_topology_3.self[i] * d->process_topology_3.n[i]; break;
	case 2: d3_coord[i]  = d->process_topology_3.self[i] * d->process_topology_3.n[i]; break;
      }
    }
    
    //now unlike above, we dont need to iterate in the z_dim, because for each processor its subarrays inward dimension is already set by the cubes z_dim.
    //Instead, each iteration of the for-loop will look at different subarrays whose locations in the cuboid differ by local x and y coords.
    
    switch(z_dim){
      //p1 is a place holder for the first translation . The outside for-loop will increment the coord in that direction, say x_dim, 
      //and keep doing so until all of the chunks in that dimension are calculated. Then it will increment p0 in the other dimension (in this example the y) 
      //and repeat until all of the subchunks in the x and y dimensions are calculated.
      //are found. 
      //Note: p0 and p1 will increment different dimensions depending of whether it is using the x y or z pencils, this is because the set up of the coordinate system for each 
      //pencil is different and to ensure that no communications hang up later, the directions coded below are unique for each type of pencil.
      case 0:
	d3_coord[y_dim] += p0 * d->process_topology_2_x.n[y_dim]; 
	d3_coord[x_dim] += p1 * d->process_topology_2_x.n[x_dim]; 
	break;
      case 1:
	d3_coord[y_dim] += p0 * d->process_topology_2_y.n[y_dim]; 
	d3_coord[x_dim] += p1 * d->process_topology_2_y.n[x_dim]; 
	break;
      case 2:
	d3_coord[x_dim] += p0 * d->process_topology_2_z.n[x_dim]; 
	d3_coord[y_dim] += p1 * d->process_topology_2_z.n[y_dim]; 
	break;
    }
    if (p1 == p1max) {
      p0++;
      p1 = 0;
    } else {
      p1++;
    }
    // create a dataype for my pth subrarray in the 3d distribution
    
    
    //d3_array_start holds the starting index of the chunk in the cubes local coordinates(note the cubes local coord system is actually the same as the grids global coord system, by set up)
    
    d3_array_start[x_dim] = d3_coord[x_dim] % cube_sizes[x_dim]; 
    d3_array_start[y_dim] = d3_coord[y_dim] % cube_sizes[y_dim]; 
    d3_array_start[z_dim] = d3_coord[z_dim] % cube_sizes[z_dim]; 
    
    //make starting point so that it coincides with the starting point of the pencil from the pencils coordinate system. (for z_pencils nothing needs to be changed, since it already
    //has the coordinate system of the grid, however, the x and y pencils have different starting points of the subchunk in their coord systems.)
    if(z_dim==0 || z_dim ==1){
      d3_array_start[2]=d3_array_start[2]+subsizes[2]-1;
    }
    if(print_me && (self==me))fprintf(stderr,"D3_array_start is (%d,%d,%d) and subsizes is (%d,%d,%d) \n",d3_array_start[0],d3_array_start[1],d3_array_start[2],subsizes[0],subsizes[1],subsizes[2]);
    
    
    //If sending cube chunks to pencils, need to fill those chunks with data here. The chunks are filled in the order 
    //such that when the pencil recieves the chunk, in its local array indexing, it assumes that the array is already 
    //filled such that it is contiguous. Therefore, complicated for-loops below fill the array in the cubes local indexing to match what the pencil will
    //expect. 
    if(direction == REDISTRIBUTE_3_TO_2){
      int64_t ch_indx=0;
      int dims_size=cube_sizes[0]*cube_sizes[1]*cube_sizes[2];
      if((self == me) && print_me)fprintf(stderr, "%d, %d, MAKE 3D Chunk...\n", self,d3_peer);
      switch(z_dim){
	case 0:
	  for(int i2=d3_array_start[y_dim];i2>d3_array_start[y_dim]-subsizes[y_dim];i2--){//perhaps y_dim
	    for(int i1=d3_array_start[x_dim];i1<d3_array_start[x_dim]+subsizes[x_dim];i1++){//perhaps x_dim
	      for(int i0=d3_array_start[z_dim];i0<d3_array_start[z_dim]+subsizes[z_dim];i0++){//perhaps z_dim
		int64_t local_indx=d->process_topology_3.n[2]*(d->process_topology_3.n[1]*i0+i1) + i2;
		assert(local_indx < dims_size);
		assert(ch_indx <chunk_size && ch_indx >= 0 && local_indx>=0 && local_indx < dims_size);
		d->d3_chunk[ch_indx]=a[local_indx];
		ch_indx++;
	      }
	    }
	  }
	  break;
	case 1:
	  for(int i0=d3_array_start[y_dim];i0<d3_array_start[y_dim]+subsizes[y_dim];i0++){
	    for(int i2=d3_array_start[x_dim];i2>d3_array_start[x_dim]-subsizes[x_dim];i2--){
	      for(int i1=d3_array_start[z_dim];i1<d3_array_start[z_dim]+subsizes[z_dim];i1++){
		int64_t local_indx=d->process_topology_3.n[2]*(d->process_topology_3.n[1]*i0+i1) + i2;
		assert(local_indx < dims_size);
		assert(ch_indx <chunk_size && ch_indx >= 0 && local_indx>=0 && local_indx < dims_size);
		d->d3_chunk[ch_indx]=a[local_indx];
		ch_indx++;
	      }
	    }
	  }
	  
	  break;
	case 2:
	  for(int i0=d3_array_start[x_dim];i0<d3_array_start[x_dim]+subsizes[x_dim];i0++){
	    for(int i1=d3_array_start[y_dim];i1<d3_array_start[y_dim]+subsizes[y_dim];i1++){
	      for(int i2=d3_array_start[z_dim];i2<d3_array_start[z_dim]+subsizes[z_dim];i2++){
		int64_t local_indx=d->process_topology_3.n[2]*(d->process_topology_3.n[1]*i0+i1) + i2;
		assert(local_indx < dims_size);
		assert(ch_indx <chunk_size && ch_indx >= 0 && local_indx>=0 && local_indx < dims_size);
		d->d3_chunk[ch_indx]=a[local_indx];
		ch_indx++;
	      }
	    }
	  }
	  
	  break;
      }
    }
    
    if (DEBUG_CONDITION || ((self == me) && print_me)) {
      fprintf(stderr,
	      "%d: pencil_sizes=(%d,%d,%d), cube_sizes=(%d,%d,%d), subsizes=(%d,%d,%d), d3_coord=(%d,%d,%d), d3_array_start=(%d,%d,%d) \n",
	      self,
	      pencil_sizes[0], pencil_sizes[1], pencil_sizes[2],
	      cube_sizes[0], cube_sizes[1], cube_sizes[2],
	      subsizes[0], subsizes[1], subsizes[2],
	      d3_coord[0], d3_coord[1], d3_coord[2],
	      d3_array_start[0],d3_array_start[1],d3_array_start[2]);
    }
    
    // what peer in the 2d distribution owns this subarray?
    for (int i = 0; i < 3; ++i) {
      switch(z_dim){
	case 0:
	  d2_peer_coord[i] = d3_coord[i] / d->process_topology_2_x.n[i];
	  break;
	case 1:
	  d2_peer_coord[i] = d3_coord[i] / d->process_topology_2_y.n[i];
	  break;
	case 2:
	  d2_peer_coord[i] = d3_coord[i] / d->process_topology_2_z.n[i];
	  break;
      }
    }
    d2_peer_coord[z_dim] = 0;//since these are pencils, there is no two pencils in this direction.
    if((self == me) && print_me)fprintf(stderr, "%d, %d, %d PENCIL that hits chunk!...\n",d2_peer_coord[0],d2_peer_coord[1],d2_peer_coord[2]);
    switch(z_dim){
      //find its rank
      case 0:
        Rank_x_pencils(&d2_peer,d2_peer_coord,d);
	break;
      case 1:
        Rank_y_pencils(&d2_peer,d2_peer_coord,d);
	break;
      case 2:
	Rank_z_pencils(&d2_peer,d2_peer_coord,d);
	break;
    }
    if((self == me) && print_me)fprintf(stderr, "%d, %d, %d Made it before comm!...\n", self,p, npeers);
    
    // record the communication to be done in a schedule. Make sure to map each grid to the correct rank
    if (direction == REDISTRIBUTE_3_TO_2) {
      recv_peer = d->rankmap[d3_peer];
      send_peer = d->rankmap[d2_peer];
    } else if (direction == REDISTRIBUTE_2_TO_3) {
      recv_peer = d->rankmap[d2_peer];
      send_peer = d->rankmap[d3_peer];
    } else {
      abort();
    }
    //comunication of the chunks:
    //if print_mess boolean is set to true, then the code runs without sending any messages, and is used to test which messages would be sent in the entire run.
    //(designed to debug comm hangups, if they occur).
    
    if(direction == REDISTRIBUTE_3_TO_2){
      
      if((self == me) && print_mess)fprintf(stderr, " I am %d, making request to recieve from %d...\n", self,recv_peer);
      if(!print_mess)MPI_Irecv((void *) d->d2_chunk, chunk_size, MPI_DOUBLE_COMPLEX, recv_peer, 0, d->process_topology_1.cart, &req1);
      
      if((self == me) && print_mess)fprintf(stderr, " I am %d, making request to send to %d...\n", self,send_peer);
      if(!print_mess)MPI_Isend((void *) d->d3_chunk, chunk_size, MPI_DOUBLE_COMPLEX, send_peer, 0, d->process_topology_1.cart, &req2);
      
      if((self == me) && print_mess)fprintf(stderr, " I am %d, waiting to recieve from %d...\n", self,recv_peer);
      //fprintf(stderr, " I am %d, waiting to recieve from %d...\n", self,recv_peer);
      if(!print_mess)MPI_Wait(&req1,MPI_STATUS_IGNORE);
      
      //if((self == me || self == 1 || self == 2 || self == 3) && print_me)fprintf(stderr, " I am %d, waiting to send to %d...\n", self,send_peer);
      //fprintf(stderr, " I am %d, waiting to send to %d...\n", self,send_peer);
      if(self==me && print_mess)fprintf(stderr, " I am %d, waiting to send to %d...\n", self,send_peer);
      if(!print_mess)MPI_Wait(&req2,MPI_STATUS_IGNORE);
      
      //fill the local array with the received chunk.
      int64_t ch_indx=0;
      int dims_size=pencil_dims[0]*pencil_dims[1]*pencil_dims[2];
      if(self==me && print_me)fprintf(stderr,"REAL SUBSIZES (%d,%d,%d)\n",subsizes[x_dim],subsizes[y_dim],subsizes[z_dim]);
      if(self==me && print_me)fprintf(stderr,"PENCIL DIMENSION VS. local sizes (%d,%d,%d) vs (%d,%d,%d)\n",pencil_dims[0],pencil_dims[1],pencil_dims[2],local_sizes[0],local_sizes[1],local_sizes[2]);
      if(self==me && print_me)fprintf(stderr,"DIM_2_ARRAY_START (%d,%d,%d) \n",d2_array_start[0],d2_array_start[1],d2_array_start[2]);
      for(int i0=d2_array_start[0];i0<d2_array_start[0]+local_sizes[0];i0++){
	for(int i1=d2_array_start[1];i1<d2_array_start[1]+local_sizes[1];i1++){
	  for(int i2=d2_array_start[2];i2<d2_array_start[2]+local_sizes[2];i2++){
	    int64_t local_indx=pencil_dims[2]*(pencil_dims[1]*i0+i1) + i2;
	    //if(self==me)fprintf(stderr,"local_indx = %d ",local_indx);
	    //if(local_indx >= dims_size)fprintf(stderr,"WOW, in third for, dims is (%d), we are %d and my rank is %d",dims_size,local_indx,self);
	    assert(local_indx < dims_size);
	    assert(ch_indx <chunk_size && ch_indx >= 0 && local_indx>=0 && local_indx < dims_size);
	    b[local_indx]=d->d2_chunk[ch_indx];
	    //if((p==0 || p==1 || p==2 || p==3 || p==4 || p==5) && self==me)fprintf(stderr,"(%f,%f) ",real(d->d2_chunk[ch_indx]),imag(d->d2_chunk[ch_indx]));
	    ch_indx++;
	  }
                        	}
      }
      //     if((p==0 ||p==1 || p==2 || p==3 || p==4 || p==5) && self==me)fprintf(stderr,"P is %d \n",p);
      
    } 
    else if (direction == REDISTRIBUTE_2_TO_3) {
      
      if((self == me) && print_mess)fprintf(stderr, " I am %d, making request to recieve from %d...\n", self,recv_peer);
      if(!print_mess)MPI_Irecv((void *) d->d3_chunk, chunk_size, MPI_DOUBLE_COMPLEX, recv_peer, 0, d->process_topology_1.cart, &req1);
      
      if((self == me) && print_mess)fprintf(stderr, " I am %d, making request to send to %d...\n", self,send_peer);
      if(!print_mess)MPI_Isend((void *) d->d2_chunk, chunk_size, MPI_DOUBLE_COMPLEX, send_peer, 0, d->process_topology_1.cart, &req2);
      
      if((self == me) && print_mess)fprintf(stderr, " I am %d, waiting to recieve from %d...\n", self,recv_peer);
      if(!print_mess)MPI_Wait(&req1,MPI_STATUS_IGNORE);
      
      if((self == me) && print_mess)fprintf(stderr, " I am %d, waiting to send to %d...\n", self,send_peer);
      if(!print_mess)MPI_Wait(&req2,MPI_STATUS_IGNORE);
      int64_t ch_indx=0;
      int dims_size=(d->process_topology_3.n[2])*(d->process_topology_3.n[1])*(d->process_topology_3.n[0]);
      if(z_dim==0){
	//fill the local array with the received chunk.
	
	for(int i2=d3_array_start[y_dim];i2>d3_array_start[y_dim]-subsizes[y_dim];i2--){
	  for(int i1=d3_array_start[x_dim];i1<d3_array_start[x_dim]+subsizes[x_dim];i1++){
	    for(int i0=d3_array_start[z_dim];i0<d3_array_start[z_dim]+subsizes[z_dim];i0++){
	      int64_t local_indx=d->process_topology_3.n[2]*(d->process_topology_3.n[1]*i0+i1) + i2;
	      //if(local_indx >= dims_size)fprintf(stderr,"WOW, in fourth for, dims is (%d), we are %d and my rank is %d",dims_size,local_indx,self);
	      assert(local_indx < dims_size);
	      assert(ch_indx <chunk_size && ch_indx >= 0 && local_indx>=0 && local_indx < dims_size);
	      b[local_indx]=d->d3_chunk[ch_indx];
	      //                         if(p==3 && self==me)fprintf(stderr,"(%f,%f) ",real(d->d3_chunk[ch_indx]),imag(d->d3_chunk[ch_indx]));
	      ch_indx++;
	    }
	  }
	}
      }
      else if(z_dim==1){
	for(int i0=d3_array_start[y_dim];i0<d3_array_start[y_dim]+subsizes[y_dim];i0++){
	  for(int i2=d3_array_start[x_dim];i2>d3_array_start[x_dim]-subsizes[x_dim];i2--){
	    for(int i1=d3_array_start[z_dim];i1<d3_array_start[z_dim]+subsizes[z_dim];i1++){
	      int64_t local_indx=d->process_topology_3.n[2]*(d->process_topology_3.n[1]*i0+i1) + i2;
	      //if(local_indx >= dims_size)fprintf(stderr,"WOW, in fourth for, dims is (%d), we are %d and my rank is %d",dims_size,local_indx,self);
	      assert(local_indx < dims_size);
	      assert(ch_indx <chunk_size && ch_indx >= 0 && local_indx>=0 && local_indx < dims_size);
	      b[local_indx]=d->d3_chunk[ch_indx];
	      //                             if(p==0 && self==me)fprintf(stderr,"(%f,%f) ",real(d->d3_chunk[ch_indx]),imag(d->d3_chunk[ch_indx]));
	      ch_indx++;
	    }
	  }
	}
	
      }
      else if(z_dim==2){
	for(int i0=d3_array_start[x_dim];i0<d3_array_start[x_dim]+subsizes[x_dim];i0++){
	  for(int i1=d3_array_start[y_dim];i1<d3_array_start[y_dim]+subsizes[y_dim];i1++){
	    for(int i2=d3_array_start[z_dim];i2<d3_array_start[z_dim]+subsizes[z_dim];i2++){
	      int64_t local_indx=d->process_topology_3.n[2]*(d->process_topology_3.n[1]*i0+i1) + i2;
	      assert(local_indx < dims_size);
	      assert(ch_indx <chunk_size && ch_indx >= 0 && local_indx>=0 && local_indx < dims_size);
	      b[local_indx]=d->d3_chunk[ch_indx];
	      //                   if(p==1 && self==me)fprintf(stderr,"(%f,%f) ",real(d->d3_chunk[ch_indx]),imag(d->d3_chunk[ch_indx]));
	      ch_indx++;
	    }
	  }
	}
	
      }
      else{
	abort();
      }
    }
    
    if (DEBUG_CONDITION) {
      fprintf(stderr,
	      "%d: npeers,p,p0,p1,p1max=(%d,%d,%d,%d,%d), "
	      "d3_coord=(%d,%d,%d), d2_peer_coord=(%d,%d,%d), "
	      "d2_coord=(%d,%d,%d), d3_peer_coord=(%d,%d,%d), "
	      "recv_peer=%d, send_peer=%d\n",
	      self,
	      npeers, p, p0, p1, p1max,
	      d3_coord[0], d3_coord[1], d3_coord[2],
	      d2_peer_coord[0], d2_peer_coord[1], d2_peer_coord[2],
	      d2_coord[0], d2_coord[1], d2_coord[2],
	      d3_peer_coord[0], d3_peer_coord[1], d3_peer_coord[2],
	      recv_peer, send_peer);
    }
    
    if((self == me) && print_me)fprintf(stderr, "%d, %d, %d Made it end-for!...\n", self,p, npeers);
  }
  
  //if((self == me) && print_me)fprintf(outfile, "   Made it all the way! for z_dim =(%d) and num_proc = (%d)...\n", z_dim, d->process_topology_1.nproc[0]);
  if((self == me) && print_result){
    FILE * outfile;
    outfile= fopen("passed.data","a");
    if (outfile) fprintf(outfile, "   Made it all the way! for z_dim =(%d) and num_proc = (%d)...\n", z_dim, d->process_topology_1.nproc[0]);
    if (outfile) fclose(outfile);
  }
//    fprintf(stderr, "%d, Made it all the way!...\n", self);
}
