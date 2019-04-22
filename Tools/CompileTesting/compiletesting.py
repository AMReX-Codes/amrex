#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import shlex
import subprocess
import calendar
import time
import argparse

def compiletesting(arg_string):

    parser = argparse.ArgumentParser(description="compile tests")
    parser.add_argument("--redo_failed", action="store_true")
    parser.add_argument("--typecheck", action="store_true")
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--make_flags", type=str, default="")
    if not arg_string is None:
        args = parser.parse_args(arg_string)
    else:
        args = parse.parse_args()

    if args.redo_failed:
        test_list = []
        f = open("failed_tests", 'r')
        for line in f.readlines():
            test_list.append(line[:-1])
        f.close()
    elif args.full:
        test_list = ['Tutorials/Basic/HelloWorld_C',
                     'Tutorials/Basic/HelloWorld_F',
                     'Tutorials/Basic/HeatEquation_EX1_C/Exec',
                     'Tutorials/Basic/HeatEquation_EX1_F',
                     'Tutorials/Amr/Advection_AmrCore/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_F/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_octree_F/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_AmrLevel/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_AmrLevel/Exec/UniformVelocity',
                     'Tutorials/EB/CNS/Exec/Combustor',
                     'Tutorials/EB/CNS/Exec/Pulse',
                     'Tutorials/EB/LevelSet/Exec',
                     'Tutorials/LinearSolvers/ABecLaplacian_C',
                     'Tutorials/LinearSolvers/ABecLaplacian_F',
                     'Tutorials/Particles/CellSortedParticles',
#                     'Tutorials/Particles/ElectromagneticPIC',
#                     'Tutorials/Particles/ElectrostaticPIC',
                     'Tutorials/Particles/NeighborList',
                     'OldTutorials/DataServicesTest0',
#                     'OldTutorials/MultiColor_C',
                     'OldTutorials/MultiFabTests_C',
                     'OldTutorials/MultiGrid_C',
                     'OldTutorials/Tiling_C',
                     'OldTutorials/Tiling_Heat_C',
                     'OldTutorials/WaveEquation_C',
                     'Tests/BBIOBenchmark',
                     'Tests/C_BaseLib',
                     'Tests/IOBenchmark',
                     'Tests/NoFort',
                     'Tests/LinearSolvers/CellEB',
                     'Tests/LinearSolvers/CellEB2',
                     'Tests/LinearSolvers/C_CellMG',
                     'Tests/LinearSolvers/ComparisonTest',
                     'Tests/LinearSolvers/C_TensorMG',
                     'Tests/MKDir']

    else:
        test_list = ['Tutorials/Amr/Advection_AmrCore/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_F/Exec/SingleVortex',
                     'Tutorials/EB/CNS/Exec/Pulse',
                     'Tutorials/LinearSolvers/ABecLaplacian_C',
                     'Tutorials/LinearSolvers/ABecLaplacian_F',
                     'Tutorials/Particles/NeighborList',
                     'Tests/NoFort',
                     'Tests/LinearSolvers/CellEB2',
                     'Tests/LinearSolvers/ComparisonTest']

    print("Test List: ", test_list)

    TOP = os.getcwd()

    start_time = calendar.timegm(time.gmtime())

    failed_tests = []
    for test in test_list:
        print("Compile", test)
        os.chdir(os.path.join(TOP,test))
        
        command = "make realclean"
        outfile = "makerealclean.ou"
        run(command, outfile)
        
        command = "make -j4 " + args.make_flags
        if args.typecheck:
            command += " typecheck"
        outfile = "make.ou"
        run(command, outfile)        

        test_success = False
        if args.typecheck:
            f = open(outfile, 'r')
            for line in f.readlines():
                if "functions checked, 0 error" in line:
                    test_success = True
                    break
        else:
            os.chdir(os.path.join(TOP,test))
            for file in os.listdir('./'):
                if file.endswith(".ex") or file.endswith(".exe"):
                    t = os.path.getmtime(file)
                    test_success = t > start_time

        if test_success:
            print("    success")
        else:
            print("    failed")
            failed_tests.append(test)

    print("")

    os.chdir(TOP)
    if failed_tests:
        print("Failed tests: ", failed_tests)
        f = open("failed_tests", 'w')
        for t in failed_tests:
            f.write(t+"\n")
        f.close()
    else:
        print("Compile tests passed.")

def run(command, outfile=None):

    # shlex.split will preserve inner quotes
    p0 = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    stdout0, stderr0 = p0.communicate()
    rc = p0.returncode
    p0.stdout.close()
    p0.stderr.close()

    if outfile:
        try: cf = open(outfile, 'wb')
        except IOError:
            print("ERROR: unable to open file for writing")
        else:
            cf.write(stdout0)
            cf.write(stderr0)
            cf.close()
    else:
        print("    ", stdout0)
        if stderr0:
            print("    ", stderr0)

    return rc

if __name__ == "__main__":
    compiletesting(sys.argv[1:])

