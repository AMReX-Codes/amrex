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

#ifndef HACC_ALIGNEDALLOCATOR_H
#define HACC_ALIGNEDALLOCATOR_H

#include <stddef.h>
#include <stdlib.h>

#include <new>

namespace hacc {
template <typename T, size_t N>
class AlignedAllocator
{
public:
  typedef T value_type;
  typedef T *pointer;
  typedef T &reference;
  typedef const T *const_pointer;
  typedef const T &const_reference;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  template <typename U>
  struct rebind {
  	typedef AlignedAllocator<U, N> other;
  };

public:
  AlignedAllocator() throw() {};
  AlignedAllocator(const AlignedAllocator&) throw() {};

  template <typename U, size_t M>
  AlignedAllocator(const AlignedAllocator<U, M>&) throw() {};

public:
  ~AlignedAllocator() throw () {};

public:
  pointer address(reference x) const { return &x; }
  const_pointer address (const_reference x) const { return &x; }

  size_type max_size() const throw() { return size_t(-1) / sizeof(T); }

  void construct(pointer p, const_reference val) { ::new ((void*)p) T(val); }
  void destroy(pointer p) { ((T*)p)->~T(); }

public:
  pointer allocate(size_type n,
                   const void * /*hint*/ = 0)
  {
    pointer p;
    if (posix_memalign((void **) &p, N, n*sizeof(T)) != 0) {
      throw std::bad_alloc();
    }

    return p;
  }

  void deallocate(pointer p, size_type n)
  {
    free((void *) p);
  }
};
} // namespace hacc

#endif // HACC_ALIGNEDALLOCATOR_H

