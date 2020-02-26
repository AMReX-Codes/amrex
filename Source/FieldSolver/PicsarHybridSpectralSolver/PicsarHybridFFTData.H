/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PICSAR_HYBRID_FFTDATA_H_
#define WARPX_PICSAR_HYBRID_FFTDATA_H_

#ifdef WARPX_USE_PSATD
#   include <fftw3.h>


// FFTData is a stuct containing a 22 pointers to auxiliary arrays
// 1-11: padded arrays in real space (required by FFTW); 12-22: arrays in spectral space
struct FFTData
{
    static constexpr int N = 22;
    void* data[N] = { nullptr };

    ~FFTData () {
        for (int i = 0; i < N; ++i) { // The memory is allocated with fftw_alloc.
            fftw_free(data[i]);
            data[i] = nullptr;
        }
    }

    FFTData () = default;

    FFTData (FFTData && rhs) noexcept {
        for (int i = 0; i < N; ++i) {
            data[i] = rhs.data[i];
            rhs.data[i] = nullptr;
        }
    }

    FFTData (FFTData const&) = delete;
    void operator= (FFTData const&) = delete;
    void operator= (FFTData&&) = delete;
};

#endif // WARPX_USE_PSATD
#endif // WARPX_PICSAR_HYBRID_FFTDATA_H_
