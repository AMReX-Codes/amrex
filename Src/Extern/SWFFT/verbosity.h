#ifndef VERBOSITY_H
#define VERBOSITY_H

// This provides a wrapper for amrex::Verbosity() to control the level of SWFFT output.

#ifdef __cplusplus
extern "C"
#endif
int verbosity ();

#endif
