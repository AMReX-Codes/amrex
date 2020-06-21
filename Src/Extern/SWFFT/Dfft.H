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

#ifndef HACC_DFFT_HPP
#define HACC_DFFT_HPP

#ifdef ESSL_FFTW
#include <fftw3_essl.h>
#else
#include <fftw3.h>
#endif

#include "complex-type.h"
#include "Distribution.H"
#include "Error.h"

// DFFT_TIMING
// 0 = no info
// 1 = summary info
// 2 = fine grain info
#ifndef DFFT_TIMING
#define DFFT_TIMING 0
#endif

#if DFFT_TIMING
#include "TimingStats.h"
#endif

#define FFTW_ADDR(X) reinterpret_cast<fftw_complex*>(&(X)[0])

namespace hacc {

class Dfft {

public:

  // 
  int global_ng(int i) const { return d.global_ng(i); }
  int const (& global_ng() const)[3] { return d.global_ng(); }
  size_t global_size() const { return d.global_size(); }
  MPI_Comm parent_comm() const { return d.parent_comm(); }



  // 
  size_t local_size() const { return d.local_size(); }



  // rank location in r-space
  int self_rspace(int i) const { return d.self_3d(i);}
  int const (& self_rspace() const)[3] { return d.self_3d();}

  // number of ranks in r-space
  int nproc_rspace(int i) const { return d.nproc_3d(i);}
  int const (& nproc_rspace() const)[3] { return d.nproc_3d();}

  // local grid dimensions in r-space
  int local_ng_rspace(int i) const { return d.local_ng_3d(i);}
  int const (& local_ng_rspace() const)[3] { return d.local_ng_3d();}

  // 3D cartesian communicator in r-space
  MPI_Comm cartcomm_rspace() const { return d.cart_3d(); }



  // rank location in k-space
  int self_kspace(int i) const { return d.self_2d_z(i);}
  int const (& self_kspace() const)[3] { return d.self_2d_z();}

  // number of ranks in r-space
  int nproc_kspace(int i) const { return d.nproc_2d_z(i);}
  int const (& nproc_kspace() const)[3] { return d.nproc_2d_z();}

  // local grid dimensions in k-space
  int local_ng_kspace(int i) const { return d.local_ng_2d_z(i);}
  int const (& local_ng_kspace() const)[3] { return d.local_ng_2d_z();}

  // 3D cartesian communicator in k-space
  MPI_Comm cartcomm_kspace() const { return d.cart_2d_z(); }



  void forward(complex_t const *in) {

    if (PlansMade != true) Error() << "Dfft buffers not set";

#if DFFT_TIMING
    double start, stop;
    double tdist=0.0, tdfft=0.0;
    start = MPI_Wtime();
#endif

    distribution_3_to_2(&in[0], &m_fs[0], &d.m_d, 0);   // in --> fs

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  fd32", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    fftw_execute(m_plan_f_x);                           // fs --> fo

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdfft += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  f1dx", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_2_to_3(&m_fo[0], &m_fs[0], &d.m_d, 0); // fo --> fs

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  fd23", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_3_to_2(&m_fs[0], &m_fo[0], &d.m_d, 1); // fs --> fo

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  fd32", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    fftw_execute(m_plan_f_y);                           // fo --> fs

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdfft += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  f1dy", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_2_to_3(&m_fs[0], &m_fo[0], &d.m_d, 1); // fs --> fo

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  fd23", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_3_to_2(&m_fo[0], &m_fs[0], &d.m_d, 2); // fo --> fs

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  fd32", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    fftw_execute(m_plan_f_z);                           // fs --> fo

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdfft += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  f1dz", stop-start);
#endif
    printTimingStats(d.parent_comm(), "DFFT fdist", tdist);
    printTimingStats(d.parent_comm(), "DFFT fdfft", tdfft);
#endif

  }



  void forward(float const *in, size_t ghost0, size_t ghost1) {

    if (PlansMade != true) Error() << "Dfft buffers not set";

#if DFFT_TIMING
    double start = MPI_Wtime();
#endif

    // copy from overloaded in to m_fo
    const int *local_dim_r = local_ng_rspace();
    int indexc = 0;
    int indexf = ghost0*(ghost0 + local_dim_r[2] + ghost1)*(ghost0 + local_dim_r[1] + ghost1);    
    for(int local_r0=0; local_r0 < local_dim_r[0]; ++local_r0) {
      indexf += ghost0*(ghost0 + local_dim_r[2] + ghost1);
      for(int local_r1=0; local_r1 < local_dim_r[1]; ++local_r1) {
        indexf += ghost0;
        for(int local_r2=0; local_r2 < local_dim_r[2]; ++local_r2) {
          m_fo[indexc] = in[indexf];
          indexf++;
          indexc++;
        } // 2
        indexf += ghost1;
      } // 1
      indexf += ghost1*(ghost0 + local_dim_r[2] + ghost1);
    } // 0

#if DFFT_TIMING
    double stop = MPI_Wtime();
    printTimingStats(d.parent_comm(), "DFFT fcopy", stop-start);
#endif

    // foward FFT from m_fo into m_fo
    forward(&m_fo[0]);
  }



  void forward(float const *in, size_t ghost) {
    return forward(in, ghost, ghost+1);
  }



  void backward(complex_t *out) {

    if (PlansMade != true) Error() << "Dfft buffers not set";

#if DFFT_TIMING
    double start, stop;
    double tdist=0.0, tdfft=0.0;
    start = MPI_Wtime();
#endif

    fftw_execute(m_plan_b_z);                           // bi --> bs

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdfft += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  b1dz", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_2_to_3(&m_bs[0], &m_bi[0], &d.m_d, 2); // bs --> bi

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  bd23", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_3_to_2(&m_bi[0], &m_bs[0], &d.m_d, 1); // bi --> bs

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  bd32", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    fftw_execute(m_plan_b_y);                           // bs --> bi

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdfft += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  b1dy", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_2_to_3(&m_bi[0], &m_bs[0], &d.m_d, 1); // bi --> bs

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  bd23", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_3_to_2(&m_bs[0], &m_bi[0], &d.m_d, 0); // bs --> bi

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  bd32", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    fftw_execute(m_plan_b_x);                           // bi --> bs

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdfft += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  b1dx", stop-start);
#endif
    MPI_Barrier(d.parent_comm());
    start = MPI_Wtime();
#endif

    distribution_2_to_3(&m_bs[0], &out[0], &d.m_d, 0);  // bs --> out

#if DFFT_TIMING
    stop = MPI_Wtime();
    tdist += stop-start;
#if DFFT_TIMING > 1
    printTimingStats(d.parent_comm(), "DFFT  bd23", stop-start);
#endif
    printTimingStats(d.parent_comm(), "DFFT bdist", tdist);
    printTimingStats(d.parent_comm(), "DFFT bdfft", tdfft);
#endif

  }



  void backward(float *out, size_t ghost0, size_t ghost1) {

    // backward FFT from m_bi into m_bi
    backward(&m_bi[0]);

#if DFFT_TIMING
    double start = MPI_Wtime();
#endif

    // copy from m_bi to overloaded out
    const int *local_dim_r = local_ng_rspace();
    int indexc = 0;
    int indexf = ghost0*(ghost0 + local_dim_r[2] + ghost1)*(ghost0 + local_dim_r[1] + ghost1);
    for(int local_r0=0; local_r0 < local_dim_r[0]; ++local_r0) {
      indexf += ghost0*(ghost0 + local_dim_r[2] + ghost1);
      for(int local_r1=0; local_r1 < local_dim_r[1]; ++local_r1) {
        indexf += ghost0;
        for(int local_r2=0; local_r2 < local_dim_r[2]; ++local_r2) {
          out[indexf] = std::real(m_bi[indexc]);
          indexf++;
          indexc++;
        } // 2
        indexf += ghost1;
      } // 1
      indexf += ghost1*(ghost0 + local_dim_r[2] + ghost1);
    } // 0

#if DFFT_TIMING
    double stop = MPI_Wtime();
    printTimingStats(d.parent_comm(), "DFFT bcopy", stop-start);
#endif

  }



  void backward(float *out, size_t ghost) {
    return backward(out, ghost, ghost+1);
  }



  Dfft(Distribution &dist)
    : d(dist), PlansMade(false)
  {}



  // arrays given at class initialization are hard-wired to the fft plans
  Dfft(Distribution &dist,
       complex_t *forward_output,
       complex_t *forward_scratch,
       complex_t *backward_input,
       complex_t *backward_scratch,
       unsigned int flags = FFTW_MEASURE)
    : d(dist), PlansMade(false)
  {
    makePlans(forward_output, 
	      forward_scratch, 
	      backward_input,
	      backward_scratch,
	      flags);
  }



  // array rules:
  // forward/backward method calls over-write these arrays
  // forward method has separate input array argument
  // backward method has separate output array argument
  // forward_scratch != forward_output && backward_scratch != backward_input
  // no required relationship between forward and backward arrays

  // actually, can these transforms be done in-place?
  // for now leave as-is, but test later

  // FFTW_MEASURE = 0
  void makePlans(complex_t *forward_output,
		 complex_t *forward_scratch,
		 complex_t *backward_input,
		 complex_t *backward_scratch,
		 unsigned int flags = FFTW_MEASURE)
  {
    if(forward_output == forward_scratch) 
      Error() << "Dfft::setBuffers() forward_output == forward_scratch";
    if(backward_input == backward_scratch) 
      Error() << "Dfft::setBuffers() backward_input == backward_scratch";

    m_fo = forward_output;
    m_fs = forward_scratch;
    m_bi = backward_input;
    m_bs = backward_scratch;

#if DFFT_TIMING
    double start = MPI_Wtime();
#endif

    // fs --> fo
    m_plan_f_x = fftw_plan_many_dft(1, // rank
				    &(d.m_d.process_topology_2_x.n[0]), // const int *n,
				    d.m_d.process_topology_2_x.n[1] * d.m_d.process_topology_2_x.n[2], // howmany
                                    FFTW_ADDR(m_fs),
                                    NULL, // const int *inembed,
                                    1, // int istride,
                                    d.m_d.process_topology_2_x.n[0], // int idist,
                                    FFTW_ADDR(m_fo),
                                    NULL, // const int *onembed,
                                    1, // int ostride,
                                    d.m_d.process_topology_2_x.n[0], // int odist,
                                    FFTW_FORWARD, // int sign,
                                    flags); // unsigned flags

    // fo --> fs
    m_plan_f_y = fftw_plan_many_dft(1, // rank
                                    &(d.m_d.process_topology_2_y.n[1]), // const int *n,
                                    d.m_d.process_topology_2_y.n[0] * d.m_d.process_topology_2_y.n[2], // howmany
                                    FFTW_ADDR(m_fo),
                                    NULL, // const int *inembed,
                                    1, // int istride,
                                    d.m_d.process_topology_2_y.n[1], // int idist,
                                    FFTW_ADDR(m_fs),
                                    NULL, // const int *onembed,
                                    1, // int ostride,
                                    d.m_d.process_topology_2_y.n[1], // int odist,
                                    FFTW_FORWARD, // int sign,
                                    flags); // unsigned flags

    // fs --> fo
    m_plan_f_z = fftw_plan_many_dft(1, // rank
                                    &(d.m_d.process_topology_2_z.n[2]), // const int *n,
                                    d.m_d.process_topology_2_z.n[1] * d.m_d.process_topology_2_z.n[0], // howmany
                                    FFTW_ADDR(m_fs),
                                    NULL, // const int *inembed,
                                    1, // int istride,
                                    d.m_d.process_topology_2_z.n[2], // int idist,
                                    FFTW_ADDR(m_fo),
                                    NULL, // const int *onembed,
                                    1, // int ostride,
                                    d.m_d.process_topology_2_z.n[2], // int odist,
                                    FFTW_FORWARD, // int sign,
                                    flags); // unsigned flags

    // bi --> bs
    m_plan_b_z = fftw_plan_many_dft(1, // rank
                                    &(d.m_d.process_topology_2_z.n[2]), // const int *n,
                                    d.m_d.process_topology_2_z.n[1] * d.m_d.process_topology_2_z.n[0], // howmany
                                    FFTW_ADDR(m_bi),
                                    NULL, // const int *inembed,
                                    1, // int istride,
                                    d.m_d.process_topology_2_z.n[2], // int idist,
                                    FFTW_ADDR(m_bs),
                                    NULL, // const int *onembed,
                                    1, // int ostride,
                                    d.m_d.process_topology_2_z.n[2], // int odist,
                                    FFTW_BACKWARD, // int sign,
                                    flags); // unsigned flags

    // bs --> bi
    m_plan_b_y = fftw_plan_many_dft(1, // rank
                                    &(d.m_d.process_topology_2_y.n[1]), // const int *n,
                                    d.m_d.process_topology_2_y.n[0] * d.m_d.process_topology_2_y.n[2], // howmany
                                    FFTW_ADDR(m_bs),
                                    NULL, // const int *inembed,
                                    1, // int istride,
                                    d.m_d.process_topology_2_y.n[1], // int idist,
                                    FFTW_ADDR(m_bi),
                                    NULL, // const int *onembed,
                                    1, // int ostride,
                                    d.m_d.process_topology_2_y.n[1], // int odist,
                                    FFTW_BACKWARD, // int sign,
                                    flags); // unsigned flags

    // bi --> bs
    m_plan_b_x = fftw_plan_many_dft(1, // rank
                                    &(d.m_d.process_topology_2_x.n[0]), // const int *n,
                                    d.m_d.process_topology_2_x.n[1] * d.m_d.process_topology_2_x.n[2], // howmany
                                    FFTW_ADDR(m_bi),
                                    NULL, // const int *inembed,
                                    1, // int istride,
                                    d.m_d.process_topology_2_x.n[0], // int idist,
                                    FFTW_ADDR(m_bs),
                                    NULL, // const int *onembed,
                                    1, // int ostride,
                                    d.m_d.process_topology_2_x.n[0], // int odist,
                                    FFTW_BACKWARD, // int sign,
                                    flags); // unsigned flags

#if DFFT_TIMING
    double stop = MPI_Wtime();
    printTimingStats(d.parent_comm(), "DFFT  init", stop-start);
#endif

    PlansMade = true;
  }



  ~Dfft() {
    if(PlansMade == true) {
      fftw_destroy_plan(m_plan_f_x);
      fftw_destroy_plan(m_plan_f_y);
      fftw_destroy_plan(m_plan_f_z);
      fftw_destroy_plan(m_plan_b_z);
      fftw_destroy_plan(m_plan_b_y);
      fftw_destroy_plan(m_plan_b_x);
    }
  }



  Distribution & get_d() {return d;}

protected:

  Distribution &d;
  bool PlansMade;
  complex_t *m_fo;
  complex_t *m_fs;
  complex_t *m_bi;
  complex_t *m_bs;
  fftw_plan m_plan_f_x;
  fftw_plan m_plan_f_y;
  fftw_plan m_plan_f_z;
  fftw_plan m_plan_b_z;
  fftw_plan m_plan_b_y;
  fftw_plan m_plan_b_x;
};

} // namespace hacc

#endif // HACC_DFFT_HPP
