
# Migration From BoxLib

To help C++ BoxLib users migrate to AMReX, we have provided a set of
tools in `Tools/Migration`.  We recommend you go through the following
steps in the order listed.  You probably should perform the migration
in a branch.  For each step, it is recommended that you start with a
fresh clone of your latest migration branch.  It is assume that the
environment variable `AMREX_HOME` is set to the AMReX directory. 

**The scripts used in this process typically perform search and
  replace in the directory tree rooted at the current directory.  So
  do not run them outside your application code.**

## Step 0

Make sure your code work with the [latest version of BoxLib on
github](https://github.com/BoxLib-Codes/BoxLib). 

## Step 1

AMReX `migration/1-amrex_home` branch should be used in this step.  In
this step, we replace `BOXLIB_HOME` with `AMREX_HOME` in the GNU Make
system by running the following command in your application directory.

    $AMREX_HOME/Tools/Migration/step-1-amrex_home/amrex_home.sh



