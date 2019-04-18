# WarpX

## Overview

**Warning: This is an alpha release of WarpX. The code is still in active development. Robustness and performance may fluctuate at this stage. The input and output formats may evolve.**

WarpX is an advanced electromagnetic Particle-In-Cell code.
It supports many features including Perfectly-Matched Layers (PML), mesh refinement, and the boosted-frame technique.

## Documentation

In order to learn how to install and run the code, please see the online documentation: https://ecp-warpx.github.io/index.html

## Testing the code

The code can be tested by running
```
./run_tests.sh
```
from the root folder of WarpX (after downloading the sources of `amrex` and `picsar`, as explained in the documentation).

The tests can be influenced by environment variables:
- `export WARPX_TEST_DIM=3` or `export WARPX_TEST_DIM=2` in order to select only the tests that correspond to this dimension
- `export WARPX_TEST_ARCH=CPU` or `export WARPX_TEST_ARCH=GPU` in order to run the tests on CPU or GPU respectively.
- `export WARPX_TEST_COMMIT=...` in order to test a specific commit.