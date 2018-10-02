#ifndef AMREX_RESTRICT_H_
#define AMREX_RESTRICT_H_

#ifdef __cplusplus

#ifdef _WIN32
#define AMREX_RESTRICT __restrict
#else
#define AMREX_RESTRICT __restrict__
#endif

#else

#define AMREX_RESTRICT restrict

#endif

#endif
