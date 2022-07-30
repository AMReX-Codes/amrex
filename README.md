<div align="center">
<img src="https://github.com/amrex-codes/amrex-codes.github.io/blob/main/images/AMReX_logo_small_banner_500.png" alt="AMReX Logo">

<p align="center">
  <a href="https://doi.org/10.21105/joss.01370">
  <img src="http://joss.theoj.org/papers/10.21105/joss.01370/status.svg" alt="Citing">
  </a>
  <a href="https://doi.org/10.5281/zenodo.2555438">
  <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.2555438.svg" alt="DOI">
  </a>
  <img src="https://github.com/AMReX-codes/amrex/workflows/cmake/badge.svg?branch=development" alt="CI: CMake on Development">
</p>


<p>
A software framework for massively parallel block-structured adaptive mesh
refinement applications.
</p>

[Overview](#Overview) -
[Features](#Features) -
[Documentation](#Documentation) -
[Gallery](#Gallery) -
[Get Help](#get-help) -
[Contribute](#Contribute) -
[License](#License) -
[Citation](#Citation)

</div>

## Overview

AMReX is a software framework designed to accelerate scientific discovery for
applications solving partial differential equations on block-structured meshes. Its
massively parallel adaptive mesh refinement (AMR) algorithms focus computational
resources and allow scalable performance on heterogeneous architectures so that
scientists can efficiently resolve details in large simulations.
AMReX is developed at [LBNL](https://www.lbl.gov/), [NREL](https://www.nrel.gov/),
and [ANL](https://www.anl.gov/) as part of the Block-Structured
AMR Co-Design Center in the Department of Energy's Exascale Computing Project.

More information is available at the [AMReX website](https://amrex-codes.github.io/).

## Features

- C++ and Fortran interfaces
- Support for cell-centered, face-centered, edge-centered, and nodal data
- Support for hyperbolic, parabolic, and elliptic solves on a hierarchical adaptive grid structure
- Optional subcycling in time for time-dependent PDEs
- Support for particles
- Embedded boundary description of irregular geometry
- Parallelization via flat MPI, OpenMP, hybrid MPI/OpenMP, or MPI/MPI
- GPU Acceleration with CUDA (NVidia), HIP (AMD) or SYCL (Intel) backends
- Parallel I/O
- Plotfile format supported by Amrvis, VisIt, ParaView and yt
- Built-in profiling tools

## Documentation

Four types of documentation are available:
- [User's Guide](https://amrex-codes.github.io/amrex/docs_html/) -- For more information about AMReX features and functions
- [Example Codes](https://amrex-codes.github.io/amrex/tutorials_html/#example-codes) -- The fastest way to start your own project
- [Guided Tutorials](https://amrex-codes.github.io/amrex/tutorials_html/GuidedTutorials.html) -- Learn basic AMReX topics in a progressive way
- [Technical Reference](https://amrex-codes.github.io/amrex/doxygen/) -- Conveniently searchable code documentation via Doxygen

## Gallery

AMReX supports several Exascale Computing Project software applications, such as
ExaSky, WarpX, Pele(Combustion), Astro, and MFiX-Exa. AMReX has also been used
in a wide variety of other scientific simulations, some of which, can be seen
in our application [gallery](https://amrex-codes.github.io/amrex/gallery.html).

<div align="center">
<img src="https://github.com/amrex-codes/amrex-codes.github.io/blob/main/images/gallery_small.gif" alt="Gallery Slideshow">
</div>

## Get Help

You can also view questions
and ask your own on our [GitHub Discussions](https://github.com/AMReX-Codes/amrex/discussions) page.  
To obtain additional help, simply post an issue.

## Contribute

We are always happy to have users contribute to the AMReX source code. To
contribute, issue a pull request against the development branch.
Any level of changes are welcomed: documentation, bug fixes, new test problems,
new solvers, etc. For more details on how to contribute to AMReX, please see
[CONTRIBUTING.md](CONTRIBUTING.md).

## License

AMReX Copyright (c) 2017, The Regents of the University of California,
through Lawrence Berkeley National Laboratory and the Alliance for
Sustainable Energy, LLC., through National Renewable Energy Laboratory
(subject to receipt of any required approvals from the U.S. Dept. of
Energy).  All rights reserved.

If you have questions about your rights to use or distribute this
software, please contact Berkeley Lab's Innovation & Partnerships
Office at IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the
U.S. Department of Energy and the U.S. Government consequently retains
certain rights. As such, the U.S. Government has been granted for
itself and others acting on its behalf a paid-up, nonexclusive,
irrevocable, worldwide license in the Software to reproduce,
distribute copies to the public, prepare derivative works, and perform
publicly and display publicly, and to permit other to do so.

License for AMReX can be found at [LICENSE](LICENSE).

## Citation

To cite AMReX, please use [![Citing](http://joss.theoj.org/papers/10.21105/joss.01370/status.svg)](https://doi.org/10.21105/joss.01370)

```
@article{AMReX_JOSS,
  doi = {10.21105/joss.01370},
  url = {https://doi.org/10.21105/joss.01370},
  year = {2019},
  month = may,
  publisher = {The Open Journal},
  volume = {4},
  number = {37},
  pages = {1370},
  author = {Weiqun Zhang and Ann Almgren and Vince Beckner and John Bell and Johannes Blaschke and Cy Chan and Marcus Day and Brian Friesen and Kevin Gott and Daniel Graves and Max Katz and Andrew Myers and Tan Nguyen and Andrew Nonaka and Michele Rosso and Samuel Williams and Michael Zingale},
  title = {{AMReX}: a framework for block-structured adaptive mesh refinement},
  journal = {Journal of Open Source Software}
}
```
