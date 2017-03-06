#!/usr/bin/env python

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
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--test_pgas", action="store_true")
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
        test_list = ['Tutorials/AMR_Adv_C/Exec/SingleVortex',
                     'Tutorials/AMR_Adv_C/Exec/UniformVelocity',
                     'Tutorials/AMR_Adv_C_v2/Exec/SingleVortex',
                     'Tutorials/AMR_Adv_F/Exec/SingleVortex',
                     'Tutorials/AMR_Adv_F/Exec/UniformVelocity',
                     'Tutorials/DataServicesTest0',
                     'Tutorials/Exp_CNS_NoSpec_F',
                     'Tutorials/GettingStarted_C',
                     'Tutorials/GettingStarted_F',
                     'Tutorials/HeatEquation_EX1_C',
                     'Tutorials/HeatEquation_EX1_F',
                     'Tutorials/HeatEquation_EX2_C',
                     'Tutorials/HeatEquation_EX2_F',
                     'Tutorials/HeatEquation_EX3_F',
                     'Tutorials/HeatEquation_EX4_F',
                     'Tutorials/HeatEquation_EX5_F',
                     'Tutorials/HelloWorld_C',
                     'Tutorials/HelloWorld_CF',
                     'Tutorials/MultiColor_C',
                     'Tutorials/MultiFabTests_C',
                     'Tutorials/MultiGrid_C',
                     'Tutorials/MultiGrid_F',
                     'Tutorials/PIC_C',
                     'Tutorials/Random_F',
                     'Tutorials/Sidecar_EX1',
                     'Tutorials/Tiling_C',
                     'Tutorials/Tiling_Heat_C',
                     'Tutorials/Tiling_Heat_F',
                     'Tutorials/TwoGrid_PIC_C',
                     'Tutorials/WaveEquation_C',
                     'Tutorials/WaveEquation_F',
                     'Tests/BBIOBenchmark',
                     'Tests/C_BaseLib',
                     'Tests/FillBoundaryComparison',
                     'Tests/IOBenchmark',
                     'Tests/LinearSolvers/C_CellMG',
                     'Tests/LinearSolvers/ComparisonTest',
                     'Tests/LinearSolvers/C_TensorMG',
                     'Tests/LinearSolvers/F_MG',
                     'Tests/MKDir',
                     'MiniApps/AMR_Adv_Diff_F90',
                     'MiniApps/FillBoundary',
                     'MiniApps/MultiGrid_C',
                     'MiniApps/SMC']
        if (args.test_pgas):
            test_list += ['Tutorials/PGAS_HEAT', 'MiniApps/PGAS_SMC']
    else:
        test_list = ['Tutorials/AMR_Adv_C/Exec/SingleVortex',
                     'Tutorials/AMR_Adv_C_v2/Exec/SingleVortex',
                     'Tutorials/AMR_Adv_F/Exec/SingleVortex',
                     'Tutorials/PIC_C',
                     'Tests/LinearSolvers/ComparisonTest']

    print "Test List: ", test_list

    TOP = os.getcwd()

    start_time = calendar.timegm(time.gmtime())

    for test in test_list:
        print "Compile", test
        os.chdir(os.path.join(TOP,test))
        
        command = "make realclean"
        outfile = "makerealclean.ou"
        run(command, outfile)
        
        command = "make -j4 " + args.make_flags
        outfile = "make.ou"
        run(command, outfile)        

        os.chdir(TOP)

    print ""

    failed_tests = []
    for test in test_list:
        test_success = False
        os.chdir(os.path.join(TOP,test))
        for file in os.listdir('./'):
            if file.endswith(".ex") or file.endswith(".exe"):
                t = os.path.getmtime(file)
                test_success = t > start_time
        if not test_success:
            failed_tests.append(test)

    os.chdir(TOP)
    if failed_tests:
        print "Failed tests: ", failed_tests
        f = open("failed_tests", 'w')
        for t in failed_tests:
            f.write(t+"\n")
        f.close()
    else:
        print "Compile tests passed."

def run(command, outfile=None):

    # shlex.split will preserve inner quotes
    p0 = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    stdout0, stderr0 = p0.communicate()
    rc = p0.returncode
    p0.stdout.close()
    p0.stderr.close()

    if outfile:
        try: cf = open(outfile, 'w')
        except IOError:
            print "ERROR: unable to open file for writing"
        else:
            cf.write(stdout0)
            cf.write(stderr0)
            cf.close()
    else:
        print "    ", stdout0
        if stderr0:
            print "    ", stderr0

    return rc

if __name__ == "__main__":
    compiletesting(sys.argv[1:])

