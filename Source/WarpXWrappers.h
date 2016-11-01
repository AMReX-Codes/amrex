#ifndef WARPX_WRAPPERS_H_
#define WARPX_WRAPPERS_H_

#ifdef __cplusplus
extern "C" {
#endif

  void boxlib_init (int argc, char*** argv);

  void boxlib_finalize ();

  void warpx_init ();

  void warpx_finalize ();

  void warpx_evolve (int numsteps);  // -1 means the inputs parameter will be used.

#ifdef __cplusplus
}
#endif

#endif
