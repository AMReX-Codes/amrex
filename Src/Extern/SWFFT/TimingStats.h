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

#ifndef HACC_TIMINGSTATS_H
#define HACC_TIMINGSTATS_H

#include <math.h>

#include <mpi.h>

// lightweight timing statistics from MPI_Wtime() calls
// C header only, no static variables
// prints maximum, average/mean, minimum, and stddev

#ifdef __cplusplus
extern "C" {
#endif

inline
void printTimingStats(MPI_Comm comm,        // comm for MPI_Allreduce()
		      const char *preamble, // text at beginning of line
		      double dt)            // delta t in seconds
{
  int myrank, nranks;
  double max, min, sum, avg, var, stdev;

  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nranks);

  MPI_Allreduce(&dt, &max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&dt, &min, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&dt, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
  avg = sum/nranks;

  dt -= avg;
  dt *= dt;
  MPI_Allreduce(&dt, &var, 1, MPI_DOUBLE, MPI_SUM, comm);
  var *= 1.0/nranks;
  stdev = sqrt(var);

  if(myrank==0) {
    printf("%s  max %.3es  avg %.3es  min %.3es  dev %.3es\n",
	   preamble, max, avg, min, stdev);
  }

  MPI_Barrier(comm);

  return;
}

#ifdef __cplusplus
}
#endif

#endif // HACC_TIMINGSTATS_H
