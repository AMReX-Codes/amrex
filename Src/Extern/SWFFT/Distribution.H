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

// I think this should work for global {ngx, ngy, ngz}

#ifndef HACC_DISTRIBUTION_HPP
#define HACC_DISTRIBUTION_HPP

#include <vector>

///
// Distribution / partition / decomposition of data
//
// A C++ wrapper around distribution.h
///

#include "complex-type.h"
#include "distribution_c.h"

namespace hacc {

class Distribution {

public:


  //
  // constructors
  //
  
  // standard setup

  Distribution(MPI_Comm comm, 
	       int const n[], 
	       bool debug = false)
    : m_comm(comm), m_rmap(NULL), m_debug(debug)
  {
    int Ndims[3] = { 0, 0, 0 };
    initialize(comm, n, Ndims);
  }

  Distribution(MPI_Comm comm, 
	       int ng, 
	       bool debug = false)
    : m_comm(comm), m_rmap(NULL), m_debug(debug)
  {
    int n[3] = { ng, ng, ng };
    int Ndims[3] = { 0, 0, 0 };
    initialize(comm, n, Ndims);
  }

  // custom setup with 3d decomposition and grid to rank map

  Distribution(MPI_Comm comm,
           int const n[],
           int const Ndims[],
           int *rmap,
           bool debug = false)
    : m_comm(comm), m_rmap(rmap), m_debug(debug)
  {
    initialize(comm, n, Ndims);
  }

  Distribution(MPI_Comm comm,
           int ng,
           int const Ndims[],
           int *rmap,
           bool debug = false)
    : m_comm(comm), m_rmap(rmap), m_debug(debug)
  {
    int n[3] = { ng, ng, ng };
    initialize(comm, n, Ndims);
  }

  //
  // destructor
  //

  virtual ~Distribution() 
  {
    distribution_fini(&m_d);
  }


  // initialization
  
  void initialize(MPI_Comm comm, int const n[], int const Ndims[]) {
    int flag;
    MPI_Initialized(&flag);
    if (flag == 0) {
      MPI_Init(0, 0);
    }

    distribution_init(comm, n, Ndims, &m_d, m_rmap, m_debug);
  }


  // redistribution

  void redistribute_1_to_3(const complex_t *a, complex_t *b) {
    distribution_1_to_3(a, b, &m_d);
  }

  void redistribute_3_to_1(const complex_t *a, complex_t *b) {
    distribution_3_to_1(a, b, &m_d);
  }

  void redistribute_2_to_3(const complex_t *a, complex_t *b, int axis) {
    distribution_2_to_3(a, b, &m_d, axis);
  }

  void redistribute_3_to_2(const complex_t *a, complex_t *b, int axis) {
    distribution_3_to_2(a, b, &m_d, axis);
  }


  // grid sizes

  size_t local_size() const {
    size_t size = 1;
    for (int i = 0; i < 3; ++i) {
      size *= (m_d.n[i] / m_d.process_topology_3.nproc[i]);
    }
    return size;
  }

  size_t global_size() const {
    size_t size = 1;
    for (int i = 0; i < 3; ++i) {
      size *= m_d.n[i];
    }
    return size;
  }

  int global_ng(int i) const { return m_d.n[i];}
  int local_ng_1d(int i) const { return m_d.process_topology_1.n[i];}
  int local_ng_2d_x(int i) const { return m_d.process_topology_2_x.n[i];}
  int local_ng_2d_y(int i) const { return m_d.process_topology_2_y.n[i];}
  int local_ng_2d_z(int i) const { return m_d.process_topology_2_z.n[i];}
  int local_ng_3d(int i) const { return m_d.process_topology_3.n[i];}

  int const (& global_ng() const)[3] { return m_d.n;}
  int const (& local_ng_1d() const)[3] { return m_d.process_topology_1.n;}
  int const (& local_ng_2d_x() const)[3] { return m_d.process_topology_2_x.n;}
  int const (& local_ng_2d_y() const)[3] { return m_d.process_topology_2_y.n;}
  int const (& local_ng_2d_z() const)[3] { return m_d.process_topology_2_z.n;}
  int const (& local_ng_3d() const)[3] { return m_d.process_topology_3.n;}


  // numbers of ranks

  int nproc() const { return m_d.process_topology_1.nproc[0];}

  int nproc_1d(int i) const { return m_d.process_topology_1.nproc[i];}
  int nproc_2d_x(int i) const { return m_d.process_topology_2_x.nproc[i];}
  int nproc_2d_y(int i) const { return m_d.process_topology_2_y.nproc[i];}
  int nproc_2d_z(int i) const { return m_d.process_topology_2_z.nproc[i];}
  int nproc_3d(int i) const { return m_d.process_topology_3.nproc[i];}

  int const (& nproc_1d() const)[3] { return m_d.process_topology_1.nproc;}
  int const (& nproc_2d_x() const)[3] { return m_d.process_topology_2_x.nproc;}
  int const (& nproc_2d_y() const)[3] { return m_d.process_topology_2_y.nproc;}
  int const (& nproc_2d_z() const)[3] { return m_d.process_topology_2_z.nproc;}
  int const (& nproc_3d() const)[3] { return m_d.process_topology_3.nproc;}


  // rank location

  int self() const { return m_d.process_topology_1.self[0];}

  int self_1d(int i) const { return m_d.process_topology_1.self[i];}
  int self_2d_x(int i) const { return m_d.process_topology_2_x.self[i];}
  int self_2d_y(int i) const { return m_d.process_topology_2_y.self[i];}
  int self_2d_z(int i) const { return m_d.process_topology_2_z.self[i];}
  int self_3d(int i) const { return m_d.process_topology_3.self[i];}

  int const (& self_1d() const)[3] { return m_d.process_topology_1.self;}
  int const (& self_2d_x() const)[3] { return m_d.process_topology_2_x.self;}
  int const (& self_2d_y() const)[3] { return m_d.process_topology_2_y.self;}
  int const (& self_2d_z() const)[3] { return m_d.process_topology_2_z.self;}
  int const (& self_3d() const)[3] { return m_d.process_topology_3.self;}


  // communicators

  MPI_Comm cart_1d() const { return m_d.process_topology_1.cart;}
  MPI_Comm cart_2d_x() const { return m_d.process_topology_2_x.cart;}
  MPI_Comm cart_2d_y() const { return m_d.process_topology_2_y.cart;}
  MPI_Comm cart_2d_z() const { return m_d.process_topology_2_z.cart;}
  MPI_Comm cart_3d() const { return m_d.process_topology_3.cart;}

  MPI_Comm parent_comm() const { return m_comm;}


  // 

  int rank_2d_x(int c[]) {
    int r;
    Rank_x_pencils(&r, c, &m_d);
    return r;
  }

  int rank_2d_y(int c[]) {
    int r;
    Rank_y_pencils(&r, c, &m_d);
    return r;
  }

  int rank_2d_z(int c[]) {
    int r;
    Rank_z_pencils(&r, c, &m_d);
    return r;
  }


  // 

  void coords_2d_x(int r, int c[]) { Coord_x_pencils(r, c, &m_d);}
  void coords_2d_y(int r, int c[]) { Coord_y_pencils(r, c, &m_d);}
  void coords_2d_z(int r, int c[]) { Coord_z_pencils(r, c, &m_d);}


public:
  // This is public for now until we refactor the Solver* classes
  // to use the C++ interface.
  distribution_t m_d;

protected:
  MPI_Comm m_comm;
  int* m_rmap;
  bool m_debug;
};

} // namespace hacc
#endif // HACC_DISTRIBUTION_HPP
