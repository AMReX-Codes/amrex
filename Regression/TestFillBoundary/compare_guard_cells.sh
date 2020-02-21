#!/bin/bash

# Copyright 2020 Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


# This script compares runs WarpX with (i) default guard cell exchanges and
# (ii) exchanging all guard cells, and check that they match to a certain
# tolerance level (typically 2.e-15). For a large number of configurations
# (PML, moving window, nox=2 etc.) it runs both versions and applies the
# AMReX fcompare utility to check if they match.
#
# This script can run in three modes:
# - Specify commit
#   ./compare_guard_cells.sh commit <commit hash>
# - specify branch
#   ./compare_guard_cells.sh branch https://github.com/ECP-WarpX/WarpX.git dev
# - Just use the current git head. Typically if you want
#   to explore a few commits manually.
#   ./compare_guard_cells.sh manual
#
# In both cases, it does the following:
# - git clone WarpX, AMReX and PICSAR (optional)
# - Compile (optional)
# - Compile fcompare (optional)
# - Loop over many configurations and compare the results
#
# The list of all simulations is stored in ran.txt
# The list of all comparisons that failed (and the fcompare output) is stored in failed.txt
#
# To use this script, create an empty directory, copy both this script and the
# inputs.2d file in an EMPTY directory (or it may mess up with your WarpX repo and other config)
# and execute it as shown above

# exit if error
set -e

# Choose your options
DO_GIT_CLONE=true
COMPILE=true
COMPILE_fcompare=true
CLEAN=true
DIM=2d
REL_TOL=2e-15
plotfile=plt00010
OMP_NUM_THREADS=1

# Some general variables
CURR_DIR=$(pwd)
WRPX_DIR=$CURR_DIR/WarpX
FCOMPARE=$CURR_DIR/fcompare.ex
EXECUTABLE=$CURR_DIR/main2d.gnu.TPROF.MPI.OMP.ex

# Read input options
MODE=$1 # "branch" or "commit" or "manual"
if [ $MODE ==  "commit" ]; then
    COMMIT=$2
elif [ $MODE == "branch" ]; then
    FORK=$2
    BRANCH=$3
elif [ $MODE == "manual" ]; then
    echo "manual mode: this script does not do any git checkout"
else
    echo "ERROR: MODE must be branch or commit or manual"
    exit
fi

if [ $DO_GIT_CLONE == true ]; then
    # Clone all three repos
    rm -rf WarpX picsar amrex
    git clone https://github.com/ECP-WarpX/WarpX.git
    git clone https://bitbucket.org/berkeleylab/picsar.git
    git clone --branch development https://github.com/AMReX-Codes/amrex.git
fi

# compile fcompare
if [ $COMPILE_fcompare == true ]; then
    # compile fcompare, and copy executable
    cd $CURR_DIR/amrex/Tools/Plotfile
    make -j 4 DIM=2
    mv $CURR_DIR/amrex/Tools/Plotfile/fcompare.gnu.ex $FCOMPARE
fi

# Compile executable to test
if [ $COMPILE == true ]; then
    cd $WRPX_DIR
    # Compile test, and copy executable
    if [ $MODE == "commit" ]; then
        git checkout $COMMIT
    fi
    if [ $MODE == "branch" ]; then
        git remote add tst $FORK
        git fetch tst
        git checkout $BRANCH
        git pull tst $BRANCH
    fi
    echo "git status for the test before compiling:"
    echo "-----------------------------------------"
    git status
    if [ $CLEAN == true ]; then
        make realclean
    fi
    make -j 4 DIM=2
    mv $WRPX_DIR/Bin/main2d.gnu.TPROF.MPI.OMP.ex $EXECUTABLE
fi

cd $CURR_DIR
if [ $CLEAN == true ]; then
    # Clean previous runs and create result directories
    rm -f Backtrace*
    rm -rf ref tst
    rm -rf OUTPUT ERROR
    rm -f ran.txt failed.txt
fi
mkdir -p OUTPUT/ref OUTPUT/tst ERROR/ref ERROR/tst

# From now on, do not exit if error (fcompare errors)
set +e

# Loop over many configurations possible
for NCI_CORR in 0 1; do
    for FILTER in 0 1; do
        for MOVING_WINDOW in 0 1; do
            for NOX in 1 3; do
                for SOLVER in yee ckc; do
                    for PML in 0 1; do
                        for MAXLEV in 0 1; do
                            for DIVE in 0 1; do
                                for SUBCYCLING in 0 1; do
                                    for NODAL in 0 1; do
                                        for GATHER in energy-conserving momentum-conserving; do
                                            # Exception: no-MR + subcycling (meaningless)
                                            if [ $MAXLEV == 0 ] && [ $SUBCYCLING == 1 ]; then continue; fi
                                            # Exception: PML + nodal (not implemented)
                                            if [ $PML == 1 ] && [ $NODAL == 1 ]; then continue; fi
                                            # Exception: divE cleaning + nodal (not implemented)
                                            if [ $DIVE == 1 ] && [ $NODAL == 1 ]; then continue; fi
                                            # Exception: NCI corrector + momentum-conserving gather
                                            if [ $NCI_CORR == 1 ] && [ $GATHER == momentum-conserving ]; then continue; fi
                                            # Exceptions: PSATD does not work with NCI corrector, subcycling,
                                            if [[ $EXECUTABLE == *"PSATD"* ]]; then
                                                if  [ $NCI_CORR == 1 ] || [ $SOLVER == ckc ] || [ $SUBCYCLING == 1 ] || \
                                                        [ $NODAL == 0 ] || [ $GATHER == momentum-conserving ] || \
                                                        [ $CURDEPO == esirkepov ]; then continue
                                                fi
                                            fi
                                            RUNNAME=NCI_CORR.$NCI_CORR.FILTER.$FILTER.MOVING_WINDOW.$MOVING_WINDOW.NOX.$NOX.SOLVER.$SOLVER.PML.$PML.MAXLEV.$MAXLEV.DIVE.$DIVE.SUBCYCLING.$SUBCYCLING.NODAL.$NODAL.GATHER.$GATHER
                                            echo $RUNNAME >> ran.txt
                                            for WHICH in ref tst; do
                                                rm -rf $WHICH
                                                SAFE_GUARD_CELLS=0
                                                if [ $WHICH == ref ]; then SAFE_GUARD_CELLS=1; fi
                                                # Run simulation
                                                sleep .01
                                                echo $SAFE_GUARD_CELLS
                                                mpirun -np 2 $EXECUTABLE inputs.$DIM \
                                                       particles.use_fdtd_nci_corr=$NCI_CORR \
                                                       warpx.use_filter=$FILTER \
                                                       warpx.do_moving_window=$MOVING_WINDOW \
                                                       interpolation.nox=$NOX interpolation.noy=$NOX interpolation.noz=$NOX \
                                                       algo.maxwell_fdtd_solver=$SOLVER \
                                                       amr.plot_file=$WHICH/plt \
                                                       warpx.do_pml=$PML \
                                                       amr.max_level=$MAXLEV \
                                                       warpx.do_dive_cleaning=$DIVE \
                                                       warpx.do_subcycling=$SUBCYCLING \
                                                       warpx.do_nodal=$NODAL \
                                                       algo.field_gathering=$GATHER \
                                                       warpx.safe_guard_cells=$SAFE_GUARD_CELLS \
                                                       >  ./OUTPUT/$WHICH/$RUNNAME.txt \
                                                       2> ./ERROR/$WHICH/$RUNNAME.txt
                                                # MPI is sometimes unhappy without this pause (??)
                                                sleep .01
                                            done
                                            # variable plot_agree is empty if and only if plotfiles disagree
                                            plot_agree=`$FCOMPARE -r $REL_TOL ./ref/$plotfile ./tst/$plotfile | grep AGREE`
                                            # if plotfiles disagree, print different and exit
                                            if [[ -z $plot_agree ]] ; then
                                                echo "!! ^   COMPARISON failed     !!" >> ran.txt
                                                echo $RUNNAME >> failed.txt
                                                $FCOMPARE -r $REL_TOL ./ref/$plotfile ./tst/$plotfile >> failed.txt
                                            fi
                                            rm -rf ref/* tst/*
                                        done
                                    done
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done
