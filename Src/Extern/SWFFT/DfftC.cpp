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

///
// Give C linkage to C++ Dfft class so that Fortran can access its functions. 
///

#include "complex-type.h"
#include "Distribution.hpp"
#include "Dfft.hpp"

extern "C" {

  hacc::Dfft* Dfft__new(hacc::Distribution &dist) {
    return new hacc::Dfft(dist);
  }

  void Dfft__makePlans(hacc::Dfft* This, complex_t *forward_output, complex_t *forward_scratch, 
                      complex_t *backward_input, complex_t *backward_scratch, unsigned int flags) {
    This->makePlans(forward_output, forward_scratch, backward_input, backward_scratch, flags);
  }

  void Dfft__forward(hacc::Dfft* This, complex_t const *in) {
    This->forward(in);
  }

  void Dfft__backward(hacc::Dfft* This, complex_t *out) {
    This->backward(out);
  }

  size_t Dfft__global_size(hacc::Dfft* This) {
    return This->global_size();
  }

  size_t Dfft__local_size(hacc::Dfft* This) {
    return This->local_size();
  }

  void Dfft__global_ng(hacc::Dfft* This, int n[3]) {
    for(size_t i = 0; i < 3; ++i) n[i] = This->global_ng(i);
  }

  void Dfft__self_rspace(hacc::Dfft* This, int n[3]) {
    for(size_t i = 0; i < 3; ++i) n[i] = This->self_rspace(i);
  } 

  void Dfft__nproc_rspace(hacc::Dfft* This, int n[3]) {
    for(size_t i = 0; i < 3; ++i) n[i] = This->nproc_rspace(i);
  }

  void Dfft__local_ng_rspace(hacc::Dfft* This, int n[3]) {
    for(size_t i = 0; i < 3; ++i) n[i] = This->local_ng_rspace(i);
  }

  void Dfft__self_kspace(hacc::Dfft* This, int n[3]) {
    for(size_t i = 0; i < 3; ++i) n[i] = This->self_kspace(i);
  }

  void Dfft__nproc_kspace(hacc::Dfft* This, int n[3]) {
    for(size_t i = 0; i < 3; ++i) n[i] = This->nproc_kspace(i);
  }

  void Dfft__local_ng_kspace(hacc::Dfft* This, int n[3]) {
    for(size_t i = 0; i < 3; ++i) n[i] = This->local_ng_kspace(i);
  }

  MPI_Fint Dfft__parent_comm(hacc::Dfft* This) {
    MPI_Comm comm = This->parent_comm();
    return MPI_Comm_c2f(comm);
  }

  MPI_Fint Dfft__cartcomm_rspace(hacc::Dfft* This) {
    MPI_Comm comm = This->cartcomm_rspace();
    return MPI_Comm_c2f(comm);
  }

  MPI_Fint Dfft__cartcomm_kspace(hacc::Dfft* This) {
    MPI_Comm comm = This->cartcomm_kspace();
    return MPI_Comm_c2f(comm);
  }

  void Dfft__delete(hacc::Dfft* This) {
    delete This;
  }

}

