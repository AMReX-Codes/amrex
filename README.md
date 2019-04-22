[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2555438.svg)](https://doi.org/10.5281/zenodo.2555438)

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

License for AMReX can be found at [license.txt](license.txt).

## Development Model

Development generally follows the following ideas:

  * New features are committed to the `development` branch.

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

    If a change is critical, we can cherry-pick the commit from
    `development` to `master`.

  * Bug fixes, questions and contributions of new features are welcome!

       * Bugs should be reported through GitHub issues
       * We suggest asking questions through GitHub issues as well
       * *Any contributions of new features that have the potential
         to change answers should be done via pull requests.*
         A pull request should be generated from your fork of
         amrex and target the `development` branch.

         If there are a number of small commits making up the PR, we may
         wish to squash commits upon merge to have a clean history.
         *Please ensure that your PR title and first post are descriptive,
         since these will be used for a squashed commit message.*

  * On the first workday of each month, we perform a merge of
    `development` into `master`.  For this merge to take place, we
    need to be passing the regression tests.

    To accommodate this need, we close the merge window into
    `development` a few days before the merge day.  While the merge
    window is closed, only bug fixes should be pushed into
    `development`.  Once the merge from `development` -> `master` is
    done, the merge window reopens.

## Core Developers

People who make a number of substantive contributions will be named
"core developers" of AMReX.  The criteria for becoming a core
developer are flexible, but generally involve one of the following:

  * 100 non-trivial commits to `amrex/Src/`  *and/or*

  * addition of a new algorithm / module  *and/or*

  * substantial input into the code design process or testing

If a core developer is inactive for multiple years, we may reassess their
status as a core developer.

The current list of core developers is: Ann Almgren (LBNL), Vince Beckner, John Bell (LBNL), Johannes Blaschke (LBNL), Cy Chan (LBNL), Marcus Day (LBNL), Brian Friesen (NERSC), Kevin Gott (NERSC), Daniel Graves (LBNL), Max Katz (NVIDIA), Andrew Myers (LBNL), Tan Nguyen (LBNL), Andrew Nonaka (LBNL), Michele Rosso (LBNL), Sam Williams (LBNL), Weiqun Zhang (LBNL), Michael Zingale (Stonybrook University).

