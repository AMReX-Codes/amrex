# CTest dashboard settings, read automatically by ctest from the top of the
# source tree. This does not affect an ordinary `ctest` run, only submissions.
#
# Tests are registered by include(CTest) in the top-level CMakeLists.txt when
# AMReX_ENABLE_TESTS=ON. To build, test and publish results in one step:
#
#   ctest -D Experimental       # ad hoc run, from a branch or a laptop
#   ctest -D Nightly            # scheduled run, grouped by the start time below
#
# Results appear at https://my.cdash.org/index.php?project=AMReX

set(CTEST_PROJECT_NAME AMReX)

# Nightly runs are grouped into the day that begins at this time, so a run that
# starts late still lands in the right column. Must match the nightly start
# time configured for the project on CDash.
set(CTEST_NIGHTLY_START_TIME "00:00:00 UTC")

set(CTEST_SUBMIT_URL https://my.cdash.org/submit.php?project=AMReX)
set(CTEST_DROP_SITE_CDASH TRUE)
