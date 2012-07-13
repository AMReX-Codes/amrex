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
        self.sourceHeads = ['','']
        self.extSrcHeads = ['','']

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


    #--------------------------------------------------------------------------
    # if we are doing a single test, remove all other tests
    #--------------------------------------------------------------------------
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

    # build list of failed tests
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
        print obj.name, "first failed at"
        for src in coms.keys():
            print "      ", src, ', '.join(coms[src][:-1])
        if oldestFailed:
            print "   The strang thing is the above test was successful last time with the same source code."
        else:
            last_commit_time = 0
            failed_src = ''
            for src in coms.keys():
                if coms[src][-1] > last_commit_time:
                    last_commit_time = coms[src][-1]
                    failed_src = src
            print "\n"
            reg_test.bold("Regtester blames "+failed_src+"!")


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
                tobj2.BoxLibHeads[0] = f.readline().rstrip('\n')
                f.close()
                if suite.sourceTree != "BoxLib":
                    f = open(d+'/git.'+suite.srcName+'.HEAD', 'r')
                    tobj2.sourceHeads[0] = f.readline().rstrip('\n')
                    f.close()
                if suite.useExtSrc:
                    f = open(d+'/git.'+suite.extSrcName+'.HEAD', 'r')
                    tobj2.extSrcHeads[0] = f.readline().rstrip('\n')
                    f.close()
                break
            else:
                tobj2.passed = False
                tobj2.lastFailure = d
                f = open(d+'/git.BoxLib.HEAD', 'r')
                tobj2.BoxLibHeads[1] = f.readline().rstrip('\n')
                f.close()
                if suite.sourceTree != "BoxLib":
                    f = open(d+'/git.'+suite.srcName+'.HEAD', 'r')
                    tobj2.sourceHeads[1] = f.readline().rstrip('\n')
                    f.close()
                if suite.useExtSrc:
                    f = open(d+'/git.'+suite.extSrcName+'.HEAD', 'r')
                    tobj2.extSrcHeads[1] = f.readline().rstrip('\n')
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

    codes = []
    commits = {}

    codes = ['BoxLib']
    commits['BoxLib'] = buildGITCommits(suite.boxLibDir, to2.BoxLibHeads)

    if suite.sourceTree != "BoxLib":
        codes.append(suite.srcName)
        commits[suite.srcName] = buildGITCommits(suite.sourceDir, to2.sourceHeads)

    if suite.useExtSrc:
        codes.append(suite.extSrcName)
        commits[suite.extSrcName] = buildGITCommits(suite.extSrcDir, to2.extSrcHeads)

    # merge commits so that we have a linear history
    icommits = {}
    for code in codes:
        icommits[code] = 0

    linCommits = [{}]
    for code in codes:
        icom = icommits[code]
        linCommits[0][code] = commits[code][icom]

    def still_more_commits(icommits):
        for code in codes:
            icom = icommits[code]
            if icom < len(commits[code]) - 1:
                return True
        return False

    while still_more_commits(icommits):
        tnext = {}
        for code in codes:
            icom = icommits[code]
            if icom < len(commits[code]) - 1:
                tnext[code] = commits[code][icom][-1]
            else:
                tnext[code] = 0
                
        whichcode = ''
        tlatest = 0
        for code in codes:
            if tnext[code] > tlatest:
                tlatest = tnext[code]
                whichcode = code

        icommits[whichcode] += 1

        tempcommits = {}
        for code in codes:
            icom = icommits[code]
            tempcommits[code] = commits[code][icom]

        linCommits.append(tempcommits)

    iSuccess = len(linCommits)-1
    iFail = 0 
    
    if iFail >= iSuccess:
        reg_test.fail("This does not make sense: iFail >= iSuccess")

    while iFail != iSuccess-1:
        iTry = (iSuccess + iFail) / 2
        status = run_old_version(linCommits[iTry], to2, suite, testFile, origdir)
        if status == 0:
            iSuccess = iTry
        else:
            iFail = iTry
    
    oldestFailed = False
    if iSuccess == len(linCommits)-1: 
        # make sure the last succsssful test still works. 
        # If it doesn't work anymore, maybe we should blame things like upgraded compiler.
        status = run_old_version(linCommits[iSuccess], to2, suite, testFile, origdir)
        if status > 0:
            oldestFailed = True

    return linCommits[iFail], oldestFailed

def run_old_version(coms, test2, suite, testFile, origdir):

    os.chdir(origdir)

    bhash = coms['BoxLib'][0]

    if suite.sourceTree != "BoxLib":
        mhash = coms[suite.srcName][0]
    else:
        mhash = ""

    if suite.useExtSrc:
        ehash = coms[suite.extSrcName][0]
    else:
        ehash = ""
        
    print " \n"
    print "  Try ...." 
    for code in coms.keys():
        print "   ", code, coms[code][0:-1] 
    print "\n"

    testargv = ['', '--do_temp_run', '--boxLibGitHash', bhash, '--sourceGitHash', mhash, '--extSrcGitHash', ehash, '--single_test', test2.name, testFile]
    return reg_test.testSuite(testargv)


if __name__== "__main__":
    reg_test_blame(sys.argv)
    
