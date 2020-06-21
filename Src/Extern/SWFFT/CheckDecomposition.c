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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#include <mpi.h>

#include "distribution.h"

static inline const char *separator(int i, int n)
{
  return i == (n - 1) ? "." : ", ";
}

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  int ndim = 3, period[3], self = 0, nproc, n[3], debug = 1;
  int explicit3d = 0, np3d = 0;

  distribution_t dist;
  distribution_t *d = &dist;
  d->debug = 1;

  if(argc < 5) {
    fprintf(stderr,"\n");
    fprintf(stderr,"USAGE: %s <ngx> <ngy> <ngz> <Nproc> [nx ny nz]\n",argv[0]);
    fprintf(stderr,"\n");
    fprintf(stderr,"Required: ng? = number of global grid vertexes in each dimension\n");
    fprintf(stderr,"Required: Nproc = total number of MPI ranks\n");
    fprintf(stderr,"Optional: n? = number of MPI ranks in each dimension for 3D Cartesian communicator if setting explicitly; Nproc = nx*ny*nz enforced\n");
    fprintf(stderr,"\n");
    exit(-1);
  }

  for(int i=0; i<ndim; i++) {
    n[i] = atoi(argv[i+1]);;
    d->n[i] = n[i];
    //d->padding[i] = 0;
  }

  nproc = atoi(argv[4]);

  if(argc >= 8) {
    explicit3d = 1;
    np3d = 1;
    for(int i=0; i<ndim; i++) {
      d->process_topology_3.nproc[i] = atoi(argv[i+5]);
      np3d *= d->process_topology_3.nproc[i];
    }
    if(np3d != nproc) {
      fprintf(stderr,"ERROR: %d * %d * %d = %d != %d\n",
	      d->process_topology_3.nproc[0],
	      d->process_topology_3.nproc[1],
	      d->process_topology_3.nproc[2],
	      np3d, nproc);
      exit(-1);
    }
  }



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

  //calculates the local dimensions (number of points in each dimension)
  d->process_topology_1.n[0] = n[0] / d->process_topology_1.nproc[0];
  d->process_topology_1.n[1] = n[1] / d->process_topology_1.nproc[1];
  d->process_topology_1.n[2] = n[2] / d->process_topology_1.nproc[2];
  


  // set up process grid with 3d decomposition (CUBE)
  if(!explicit3d) {
    d->process_topology_3.nproc[0] = 0;
    d->process_topology_3.nproc[1] = 0;
    d->process_topology_3.nproc[2] = 0;
    period[0] = period[1] = period[2] = 1;
    MPI_Dims_create(nproc, ndim, d->process_topology_3.nproc);
  }

  if(self == 0) {
    printf("distribution 3D: [%d:%d:%d]\n",
	   d->process_topology_3.nproc[0],
	   d->process_topology_3.nproc[1],
	   d->process_topology_3.nproc[2]);
    fflush(stdout);
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
      fprintf(stderr,"Need to fix Z PENCILS z_procs(%d,%d,%d) 3d.ns(%d,%d,%d) 2d_z.ns(%d,%d,%d)\n", 
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
	fprintf(stderr,"Swaping Z pencils in initialization  (%d,%d,%d)\n", 
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
//  if that did not work, make a pencil that does if inside the 3d cuboids by 
//  taking the cuboids dimensions (np1,np2,np3) and making pencils 
//  (np1,np2*np3,1), or (np1*np3,np2,1) on the most evenly distributed 
//  dimensions
  if(!check_z_dims){
    if(self==0 && debug)
      fprintf(stderr,"MAKING Z PENCILS FIT zprocs(%d,%d,%d) z.ns(%d,%d,%d)\n", 
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
      fprintf(stderr,"MAKING Z PENCILS FIT AFTER zprocs(%d,%d,%d) z.ns(%d,%d,%d)\n", 
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
  //assert(check_z_dims);
  if(!check_z_dims)
    fprintf(stderr,"assert(check_z_dims) would have failed.\n");

//  if this happens, it is because the dimensions were chosen incorrectly. 
//  Either to many processors for the number of points in one dimenison (could 
//  not do at least 1 point per processor), or the methods above could 
//  not make a distribution of pencils that fit in the cubiods, which would 
//  happen if the user gave numbers that wouldent work (we require the number 
//  of processors in each dimension of the cuboid must be modulo the number of 
//  points in that dimension, otherwise, this error will happen).

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
      && (n[0] % (d->process_topology_2_x.nproc[0]) == 0);
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
	&& (n[0] % (d->process_topology_2_x.nproc[0]) == 0);
    } 
  } else{
    check_x_dims=false;
  }
//    if that did not work, make a pencil that does by taking the cuboid 
//    (np1,np2,np3) and making pencils of the form (1,np2*np1,np3) or 
//    (1,np2*np1,np3) depending on the most even distribution it can.
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
	&& (n[0] % (d->process_topology_2_x.nproc[0]) == 0);
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
  //assert(check_x_dims);
  if(!check_x_dims)
    fprintf(stderr,"assert(check_x_dims) would have failed.\n");
//  if this happens, it is because the dimensions were chosen incorrectly. 
//  Either to many processors for the number of points in one dimenison (could 
//  not do at least 1 point per processor), or the methods above could not make 
//  a distribution of pencils that fit in the cubiods, which would happen if the 
//  user gave numbers that wouldent work (we require the number of processors in 
//  each dimension of the cuboid must be modulo the number of points in that 
//  dimension, otherwise, this error will happen).

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
//  if that did not work, make a pencil that does by taking the cuboid 
//  (np1,np2,np3) and making pencils of the form (np1,1,np3*np2) or 
//  (np1*np2,1,np3) depending on the most even distribution it can.
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
  //assert(check_y_dims);
  if(!check_y_dims)
    fprintf(stderr,"assert(check_y_dims) would have failed.\n");

//  if this happens, it is because the dimensions were chosen incorrectly. 
//  Either to many processors for the number of points in one dimenison (could 
//  not do at least 1 point per processor), or the methods above could 
//  not make a distribution of pencils that fit in the cubiods, which would 
//  happen if the user gave numbers that wouldent work (we require the number of 
//  processors in each dimension of the cuboid must be modulo the number of 
//  points in that dimension, otherwise, this error will happen).

  if(self == 0) {
    printf("distribution 2y: [%d:%d:%d]\n",
	   d->process_topology_2_y.nproc[0],
	   d->process_topology_2_y.nproc[1],
	   d->process_topology_2_y.nproc[2]);
    fflush(stdout);
  }



  MPI_Finalize();
  return 0;
}
