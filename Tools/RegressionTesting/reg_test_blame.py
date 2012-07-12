#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import datetime
import time
import string
import tarfile
import subprocess
import smtplib
import email
import getpass
import socket
import testnew as reg_test

class testObj2:
    def __init__(self, name):
        self.name = name
        self.passed = True
        self.lastFailure = ''
        self.lastSuccess = ''
        self.BoxLibHeads = ['','']
        self.thisCodeHeads = ['','']

def reg_test_blame(argv):
    usage = """
    ./reg_test_blame [--single_test test]
       testfile.ini
    """

    if len(sys.argv) == 1:
        print usage
        sys.exit(2)
        
    try:
        opts, next = getopt.getopt(argv[1:], "",
                                   ["single_test="])

    except getopt.GetoptError:
        print "invalid calling sequence"
        print usage
        sys.exit(2)

    # defaults
    single_test = ""
    
    for o, a in opts:
        if o == "--single_test":
            single_test = a
            
    try:
        testFile = next[0]

    except IndexError:
        print "ERROR: a test file was not specified"
        print usage
        sys.exit(2)

    origdir = os.getcwd()

    #--------------------------------------------------------------------------
    # read in the test information
    #--------------------------------------------------------------------------
    reg_test.bold("loading " + testFile)

    suite, testList = reg_test.LoadParams(testFile)

    if (len(testList) == 0):
        reg_test.fail("No valid tests defined")

    if not single_test == "":
        found = 0
        for obj in testList:
            if (obj.name == single_test):
                found = 1
                newTestList = [obj]
                break

        if (not found):
            fail("ERROR: %s is not a valid test" % (single_test))
        else:
            testList = newTestList

    webDirs = buildWebDirs(suite)

    failedTestList = []
    for obj in testList:
        to2 = check_test_status(obj, suite, webDirs)
        if not to2.passed:
            failedTestList.append(to2)

    if not failedTestList:
        reg_test.bold("All tests are successful. No one to blame.")
        sys.exit(0)
    else:
        print "\n"
        reg_test.bold("List of failed test:")
        for obj in failedTestList:
            print "  %s " % obj.name

    badCommitList = []
    for obj in failedTestList:
        print "\n"
        reg_test.bold("Working on "+obj.name)
        bad_commits, oldestFailed = find_someone_to_blame(obj, suite, testFile, origdir)
        badCommitList.append((bad_commits, oldestFailed))

    print "\n"
    reg_test.bold("reporting...")
    for bad, obj in zip(badCommitList,failedTestList) :
        coms = bad[0]
        oldestFailed = bad[1]
        print "\n"
        print obj.name, "first failed at BoxLib", ', '.join(coms[0][:-1])
        if suite.sourceTree != "BoxLib":
            print "     plus", suite.suiteName, ', '.join(coms[1][:-1])
        if oldestFailed:
            print "   The strang thing is the above test was successful last time with the same source code."
        else:
            if suite.sourceTree == "BoxLib" or coms[0][-1] > coms[1][-1] :
                print "Regtester blames BoxLib!"
            else:
                print "Regtester blames", suite.suiteName+"!"


def buildWebDirs(suite):
    os.chdir(suite.webTopDir)
    validDirs = []
    for file in os.listdir(suite.webTopDir):
        if (file.startswith("20") and os.path.isdir(file)):
            # look for the status file
            statusFile = file + '/' + file + '.status'
            if (os.path.isfile(statusFile)):
                validDirs.append(file)

    validDirs.sort()
    validDirs.reverse()
    return validDirs

def check_test_status(tobj, suite, webDirs):
    tobj2 = testObj2(tobj.name)
    os.chdir(suite.webTopDir)
    for d in webDirs:
        sFile = d+'/'+tobj.name+'.status'
        if os.path.isfile(sFile):
            sf = open(sFile, 'r')
            testPassed = 0
            for line in sf.readlines():
                if (string.find(line, "PASSED") >= 0):
                    testPassed = 1
                    break
            if testPassed and tobj2.lastFailure == '':
                tobj2.passed = True
                return tobj2
            elif testPassed:
                tobj2.lastSuccess = d
                f = open(d+'/git.BoxLib.HEAD', 'r')
                tobj2.BoxLibHeads[0] = f.readline().rstrip('\n').strip()
                f.close()
                f = open(d+'/git.'+suite.suiteName+'.HEAD', 'r')
                tobj2.thisCodeHeads[0] = f.readline().rstrip('\n').strip()
                f.close()
                break
            else:
                tobj2.passed = False
                if tobj2.lastFailure == '':
                    tobj2.lastFailure = d
                f = open(d+'/git.BoxLib.HEAD', 'r')
                tobj2.BoxLibHeads[1] = f.readline().rstrip('\n')
                f.close()
                f = open(d+'/git.'+suite.suiteName+'.HEAD', 'r')
                tobj2.thisCodeHeads[1] = f.readline().rstrip('\n')
                f.close()
    return tobj2

def buildGITCommits(topDir, heads):
    commits = []
    os.chdir(topDir)
    prog = ["git", "log", heads[0]+'~..'+heads[1], '--format="%h|%cd|%cn|%ce|%ct"']
    p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()
    for line in stdout.split('\n'):
        if line:
            words = line.strip('"').split('|')
            commits.append((words[0], words[1], words[2], words[3], int(words[4])))
    return commits    

def find_someone_to_blame(to2, suite, testFile, origdir):
    boxLibCommits = buildGITCommits(suite.boxLibDir, to2.BoxLibHeads)
    LboxLib = len(boxLibCommits)

    if suite.sourceTree != "BoxLib":
        thisCodeCommits = buildGITCommits(suite.sourceDir, to2.thisCodeHeads)
        LthisCode = len(thisCodeCommits)

    if suite.sourceTree == "BoxLib":
        allCommits = []
        for bc in boxlibcommits:
            allCommits.append((bc,))
    else:
        ib = 0
        ic = 0
        allCommits = [(boxLibCommits[ib], thisCodeCommits[ic])]
        while ib < LboxLib-1 or ic < LthisCode-1 :
            tbnext = 0
            tcnext = 0
            if ib < LboxLib-1:
                tbnext = boxLibCommits[ib+1][-1]
            if ic < LthisCode-1:
                tcnext = thisCodeCommits[ic+1][-1]
            if tbnext > tcnext:
                ib = ib+1
            else:
                ic = ic+1
            allCommits.append((boxLibCommits[ib], thisCodeCommits[ic]))

    iSuccess = len(allCommits)-1
    iFail = 0 
    
    if iFail >= iSuccess:
        reg_test.fail("This does not make sense: iFail >= iSuccess")

    while iFail != iSuccess-1:
        iTry = (iSuccess + iFail) / 2
        status = run_old_version(allCommits[iTry], to2, suite, testFile, origdir)
        if status == 0:
            iSuccess = iTry
        else:
            iFail = iTry
    
    oldestFailed = False
    if iSuccess == len(allCommits)-1: 
        # make sure the last succsssful test still works. 
        # If it doesn't work anymore, maybe we should blame things like upgraded compiler.
        status = run_old_version(allCommits[iSuccess], to2, suite, testFile, origdir)
        if status > 0:
            oldestFailed = True

    return allCommits[iFail], oldestFailed

def run_old_version(coms, test2, suite, testFile, origdir):

    os.chdir(origdir)

    bhash = coms[0][0]
    if suite.sourceTree == "BoxLib":
        chash = coms[0][0]
    else:
        chash = coms[1][0]        
    print " \n"
    print "  Try BoxLib", coms[0][0:-1]
    if suite.sourceTree != "BoxLib":
        print "    plus", suite.suiteName, coms[1][0:-1]
    print "\n"

    testargv = ['', '--use_old_version', '--boxLibGitHash', bhash, '--sourceGitHash', chash, '--single_test', test2.name, testFile]
    return reg_test.testSuite(testargv)


if __name__== "__main__":
    reg_test_blame(sys.argv)
    
