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
    parser.add_argument("--test_pgas", action="store_true",default=False)
    parser.add_argument("--test_boxlib", action="store_true",default=False)
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
                     'Tutorials/Basic/HeatEquation_EX1_C',
                     'Tutorials/Basic/HeatEquation_EX1_F',
                     'Tutorials/Amr/Advection_AmrCore/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_F/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_octree_F/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_AmrLevel/Exec/SingleVortex',
                     'Tutorials/Amr/Advection_AmrLevel/Exec/UniformVelocity',
                     'OldTutorials/DataServicesTest0',
                     'OldTutorials/GettingStarted_C',
                     'OldTutorials/HeatEquation_EX2_C',
                     'OldTutorials/MultiColor_C',
                     'OldTutorials/MultiFabTests_C',
                     'OldTutorials/MultiGrid_C',
#                     'OldTutorials/PIC_C',
                     'OldTutorials/Sidecar_EX1',
                     'OldTutorials/Tiling_C',
                     'OldTutorials/Tiling_Heat_C',
#                     'OldTutorials/TwoGrid_PIC_C',
                     'OldTutorials/WaveEquation_C',
                     'Tests/BBIOBenchmark',
                     'Tests/C_BaseLib',
                     'Tests/IOBenchmark',
                     'Tests/LinearSolvers/C_CellMG',
                     'Tests/LinearSolvers/ComparisonTest',
                     'Tests/LinearSolvers/C_TensorMG',
                     'Tests/MKDir',
                     'MiniApps/FillBoundary',
                     'MiniApps/MultiGrid_C']
        if (args.test_pgas):
            test_list += ['Tests/FillBoundaryComparison',
                          'OldTutorials/PGAS_HEAT',
                          'MiniApps/PGAS_SMC']
        if (args.test_boxlib):
            test_list += ['OldTutorials/AMR_Adv_F/Exec/SingleVortex',
                          'OldTutorials/AMR_Adv_F/Exec/UniformVelocity',
                          'OldTutorials/Exp_CNS_NoSpec_F',
                          'OldTutorials/GettingStarted_F',
                          'OldTutorials/HeatEquation_EX1_F',
                          'OldTutorials/HeatEquation_EX2_F',
                          'OldTutorials/HeatEquation_EX3_F',
                          'OldTutorials/HeatEquation_EX4_F',
                          'OldTutorials/HeatEquation_EX5_F',
                          'OldTutorials/MultiGrid_F',
                          'OldTutorials/Random_F',
                          'OldTutorials/Tiling_Heat_F',
                          'OldTutorials/WaveEquation_F',
                          'Tests/LinearSolvers/F_MG',
                          'MiniApps/SMC']

    else:
        test_list = ['Tutorials/Amr/Advection_AmrCore/Exec/SingleVortex',
                     'OldTutorials/AMR_Adv_CF/Exec/SingleVortex',
#                     'OldTutorials/PIC_C',
                     'Tests/LinearSolvers/ComparisonTest']
        if (args.test_boxlib):
            test_list += ['Tutorials/Basic/HeatEquation_EX1_F',
                          'OldTutorials/AMR_Adv_F/Exec/SingleVortex']

    print "Test List: ", test_list

    TOP = os.getcwd()

    start_time = calendar.timegm(time.gmtime())

    failed_tests = []
    for test in test_list:
        print "Compile", test
        os.chdir(os.path.join(TOP,test))
        
        command = "make realclean"
        outfile = "makerealclean.ou"
        run(command, outfile)
        
        command = "make -j4 " + args.make_flags
        outfile = "make.ou"
        run(command, outfile)        

        test_success = False
        os.chdir(os.path.join(TOP,test))
        for file in os.listdir('./'):
            if file.endswith(".ex") or file.endswith(".exe"):
                t = os.path.getmtime(file)
                test_success = t > start_time
        if test_success:
            print ("    success")
        else:
            print ("    failed")
            failed_tests.append(test)

    print ""

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

