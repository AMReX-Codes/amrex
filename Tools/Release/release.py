#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import readline
import subprocess
import ppCleanup as ppc
    
def release(argv):

    usage = """
    /path/to/release.py CodeName CLEANWORD_FILE
    """
    if len(argv) != 3:
        print usage
        sys.exit(2)

    code = argv[1]
    cwfile = argv[2]

    if code.lower() == 'castro':
        code = 'Castro'
        private_src = 'gamera:/usr/local/gitroot/Castro'
        public_src = 'gamera:/usr/local/releases/Castro'
    else:
        sys.exit(1)

    if not os.path.exists(cwfile):
        print "clean word file:",cwfile,"does't exist."
        sys.exit(1)


    workdir = os.getcwd()

    my_private_dir = workdir+'/private/'
    my_public_dir = workdir+'/public/'

    my_private_git = os.path.join(my_private_dir, code)
    my_public_git  = os.path.join( my_public_dir, code)

    if os.path.exists(my_private_git):
        print "\nTry to clone the private git repo to", my_private_git
        print "But", my_private_git,"already exists."
        sys.exit(1)
    elif not os.path.exists(my_private_dir):
        os.makedirs(my_private_dir)

    doGitClone(my_private_dir, private_src)

    if os.path.exists(my_public_git):
        print "\nTry to clone the public git repo to", my_public_git
        print "But", my_public_git,"already exists."
        sys.exit(1)
    elif not os.path.exists(my_public_dir):
        os.makedirs(my_public_dir)

    doGitClone(my_public_dir, public_src)

    print "\nrsync from private to public ..."
    systemCall("rsync -ac --delete --exclude '.git' "
               +os.path.normpath(my_private_git)+"/ "
               +os.path.normpath(my_public_git))

    print "\nCleaning up IFDEFs ..."
    for root, dirs, files in os.walk(my_public_git):
        for name in files:
            f = os.path.join(root,name)
            fext = os.path.splitext(f)[1]
            if fext in ['.f','.f90','.F','.F90','.c','.cpp','.h','.H']:
                ftmp = f+'.gitorig'
                os.rename(f,ftmp)
#                systemCall("~/mygitrepo/BoxLib/Tools/ppCleanup/ppCleanup.py -c "+
#                           "~/mygitrepo/BoxLib/Tools/ppCleanup/cleanWords.txt "+
#                           " -o "+f+" "+ftmp)    
                ppc.ppCleanup(cwfile, ftmp, f)
                os.remove(ftmp)
        if '.git' in dirs:
            dirs.remove('.git')

    doGitCommit(my_public_git)

    last_tag = getGitTag(my_public_git)
    print "\nThe last release is called", last_tag
    while True:
        next_tag = raw_input("\nWhat's the next release called? ")
        next_tag = next_tag.strip()
        uin = raw_input("\nIt will be called '"+next_tag+"'. Is that right? (y or n) ")
        if uin == 'y':
            break
        else:
            print "\nLet's try again."

    doGitTag(my_public_git, next_tag)
    doGitTag(my_private_git, next_tag)

    print "\nYou now need to make an important choice!  There are three things need to be done."
    print "   (1) Push the new public release using 'git push' in", my_public_git+"."
    print "   (2) Push the new tag using 'git push --tags' in", my_public_git+"."
    print "   (3) Push the new tag using 'git push --tags' in", my_private_git+"."
    while True:
        uin = raw_input("\nDo you want these to be done by the script? (y or n) ")
        if uin == 'y':
            break
        elif uin == 'n':
            break
        else:
            print "What did you type?", uin+"?"
            
    if uin == 'y':
        print "\nOK. The script will do it for you."
        doGitPush(my_public_git)
        doGitPushTag(my_public_git)
        doGitPushTag(my_private_git)
        print "\nCongratulations! ", next_tag, "has been released!"
    else:
        print "\nOK. You are going to do it yourself."


def doGitClone(d, s):
    d0 = os.getcwd()
    os.chdir(d)
    prog = ["git", "clone", s]
    print "\n Cloning", s
    p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    p.stdin.close()
    os.chdir(d0)
    return


def doGitCommit(d):
    d0 = os.getcwd()
    os.chdir(d)

    prog = ["git", "status", "-s"]
    p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    p.stdin.close()    

    ulist = []

    lines = stdout.splitlines()
    for line in lines:
        words = line.split()
        if words[0] == "??": # untracked
            ulist.append(words[1])

    prog = ["git", "add"]
    prog.extend(ulist)
    p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    p.stdin.close()    

    print "\nCommiting changes"
    prog = ["git", "commit", "-a"]
    p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    p.stdin.close()    

    os.chdir(d0)
    return


def getGitTag(d):
    d0 = os.getcwd()
    os.chdir(d)
    prog = ["git", "for-each-ref", "--format", "'%(refname) %(taggerdate)'", "refs/tags"]
    p = subprocess.Popen(prog, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    os.chdir(d0)
    lastline = '\n'.join(stdout.splitlines()[-1:])
    words = lastline.split()
    return words[0].split('/')[-1]


def doGitTag(d,t):
    d0 = os.getcwd()
    os.chdir(d)
    print "\nTagging in", d
    prog = ["git", "tag", "-f", "-m", t, t]
    p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    p.stdin.close()
    os.chdir(d0)
    return


def doGitPush(d):
    d0 = os.getcwd()
    os.chdir(d)
    prog = ["git", "push"]
    print "\n Pushing in", d
    p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    p.stdin.close()
    os.chdir(d0)
    return


def doGitPushTag(d):
    d0 = os.getcwd()
    os.chdir(d)
    prog = ["git", "push", "--tags"]
    print "\n Pushing tags in", d
    p = subprocess.Popen(prog, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    p.stdout.close()
    p.stderr.close()
    p.stdin.close()
    os.chdir(d0)
    return
    

def systemCall(string):    
    status = os.system('bash -c "' + string + '"')
    return status


if __name__== "__main__":
    release(sys.argv)
