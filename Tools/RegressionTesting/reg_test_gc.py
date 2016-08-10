#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import string
import argparse

import test_util
import params
import test_report as report

def reg_test_gc(argv):
    usage = """
    ./reg_test_gc [--before|-b 2000-00-00]
       testfile.ini
    """

    if len(sys.argv) == 1:
        print usage
        sys.exit(2)
        
    try:
        opts, next = getopt.getopt(argv[1:], "b:",
                                   ["before="])

    except getopt.GetoptError:
        print "invalid calling sequence"
        print usage
        sys.exit(2)

    # defaults
    gcdate = ""
    
    for o, a in opts:
        if o == "--before" or o == "-b" :
            gcdate = a

    try:
        testFile = next[0]

    except IndexError:
        print "ERROR: a test file was not specified"
        print usage
        sys.exit(2)

    if not gcdate:
        print "ERROR: date was not specified"
        print usage
        sys.exit(2)
            
    gcd = valid_date(gcdate)
    if gcd == '':
        print "ERROR: invalid date", gcdate
        print usage
        sys.exit(2)


    workdir = os.getcwd()

    print "loading ", testFile

    args=test_util.get_args([testFile])

    suite, testList = params.load_params(args)
    activeTestList = [t.name for t in testList]

    benchmarkTestList = [t for t in testList if not (t.compileTest or t.restartTest)]
    benchmarkNotFound = {}
    for t in benchmarkTestList:
        benchmarkNotFound[t.name] = ''


    ### clean up the web dir
    print "\ncleaning ", suite.webTopDir

    os.chdir(suite.webTopDir)
    validDirs = []
    for d in os.listdir(suite.webTopDir):
        if (d.startswith("20") and os.path.isdir(d)):
            statusFile = d + '/' + d + '.status'
            if (os.path.isfile(statusFile)):
                validDirs.append(d)
    validDirs.sort()
    validDirs.reverse()
    
    latestBMDate = {}

    for d in validDirs:
        bmtests = benchmarkNotFound.keys()
        if d >= gcd and bmtests:
            if isBenchmarkDir(d):
                for t in bmtests:
                    if findBenchmark(d,t):
                        del benchmarkNotFound[t]
                        latestBMDate[t] = d
        else:
            if isBenchmarkDir(d) and bmtests:
                found = False
                for t in bmtests:
                    if findBenchmark(d,t):
                        found = True
                        del benchmarkNotFound[t]
                        latestBMDate[t] = d
                if not found:
                    rmDir(d)
            else:
                rmDir(d)


    ### clean up the test dir
    testDirs = os.path.join(suite.testTopDir,suite.suiteName+"-tests")
    print "\ncleaning ", testDirs

    os.chdir(testDirs)
    validDirs = []
    for d in os.listdir(testDirs):
        if (d.startswith("20") and os.path.isdir(d)):
            validDirs.append(d)
    validDirs.sort()
    validDirs.reverse()

    for d in validDirs:
        if d < gcd:
            tests = [t for t in os.listdir(d) if os.path.isdir(os.path.join(d,t))]
            found = False
            for t in tests:
                if t in latestBMDate.keys() and latestBMDate[t] == d:
                    found = True
                    break
            if not found:
                rmDir(d)
    
    print "\ncreating suite report..."
    report.report_all_runs(suite, activeTestList)

    print "\nGarbage cleaning finished."


def valid_date(gcdate):
    try:
        y,m,d = gcdate.split("-")
    except ValueError:
        return ''

    try:
        yi = int(y)
        if yi > 2099 or yi < 2000:
            return ''
    except ValueError:
        return ''

    try:
        mi = int(m)
        if mi > 12 or mi < 1:
            print 'm='+m+'!'
            return ''
    except ValueError:
        return ''

    try:
        di = int(d)
        if di > 31 or di < 1:
            return ''
    except ValueError:
        return ''

    return '-'.join([y,m.zfill(2),d.zfill(2)])


def isBenchmarkDir(d):
    f = open(os.path.join(d,d+'.status'), 'r')
    line = f.readline()
    f.close()
    if string.find(line, 'BENCHMARKS UPDATED') == -1:
        return False
    else:
        return True


def findBenchmark(d, t):
    try:
        f = open(os.path.join(d,t+'.status'), 'r')
    except:
        return False
    line = f.readline()
    f.close()
    if string.find(line, 'benchmarks updated') == -1:
        return False
    else:
        return True


def rmDir(d):
    print '  deleting', d
    shutil.rmtree(d)
    

if __name__== "__main__":
    reg_test_gc(sys.argv)
