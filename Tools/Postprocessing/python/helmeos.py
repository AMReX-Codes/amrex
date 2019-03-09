# this module is intented to be imported by another python program
#
# it builds a python wrapper around the fortran source code helmeos.f90 using
# the f2py framework.  The easiest way to get this module to work is to make
# sure that 1) your PYTHONPATH can find this module for importing, and 2) 
# setting a FPARALLEL environment variable
#  e.g. (for bash)
#   export FPARALLEL=~/MAESTRO/fParallel
#

import os
import sys
import subprocess
import numpy

# figure out if we need to build the helmeos python wrapper
_built = True
try:
    import fhelmEOS
except ImportError:
    _built = False

# is the helmeos initialized?
_initialized = False

# relative to fParallel
helmEOS_path = '/extern/EOS/helmeos/'

# useful strings
helmEOS = 'helmeos.f90'
helmTable = 'helm_table.dat'
vectorEOS = 'vector_eos_f90.dek'

# useful factor
ZERO = 0.0e0

# EOS input modes
input_rt = 1
input_rh = 2
input_tp = 3
input_rp = 4
input_re = 5
input_ps = 6


# this method calls f2py to build the shared object file
# it needs to know, or it tries to make a naive guess, of the location of
# the fParallel directory
def _build_fhelmEOS(fParallel_path=None):
    global helmEOS_path

    # let's look for the fParallel directory if it is not specified
    if not fParallel_path:
        # check environment
        if "FPARALLEL" in os.environ.keys():
            fParallel_path = os.environ["FPARALLEL"]
        # check current directory
        if os.path.isdir('fParallel'): fParallel_path = './fParallel'
        # check one dir up
        elif os.path.isdir('../fParallel'): fParallel_path = '../fParallel'

    if not fParallel_path or not os.path.isdir(fParallel_path):
        print "Could not find fParallel location as specified: %s" % (
            fParallel_path)
        sys.exit(1)

    fParallel_path = os.path.abspath(fParallel_path)
    print "fParallel found at: %s" % fParallel_path

    helmEOS_path = fParallel_path + helmEOS_path

    # open and read the helmeos.f90 file if possible
    try: 
        ifh = open(helmEOS_path+helmEOS)
    except IOError:
        print "Could not locate %s at expected location %s" % (
            helmeos,helmEOS_path)
        sys.exit(1)

    fileContents = ifh.readlines()
    ifh.close()

    # fix some things so we don't have to worry about modules
    for lineno, line in enumerate(fileContents):
        # remove the bl_error and bl_warn stuff
        if line.find("use bl_error") >=0: 
            fileContents[lineno] = '!' + line
        elif line.find("bl_error") >= 0:
            txt = line.split('(')[1].split(')')[0]
            fileContents[lineno] = "print *, " + txt + "\n stop \n"
        elif line.find("bl_warn") >= 0:
            txt = line.split('(')[1].split(')')[0]
            fileContents[lineno] = "print *, " + txt + "\n"
    
    # dump the updated helmeos.f90 to a temporary, local file
    try:
        ofh = open(helmEOS,'w')
    except IOError:
        print "Couldn't open file %s for writing." % helmEOS
        sys.exit(1)

    ofh.writelines(fileContents)
    ofh.close()

    # grab the vector_eos_f90.dek file
    try:
        subprocess.check_call(["cp", helmEOS_path + vectorEOS, "."])
    except subprocess.CalledProcessError:
        print "Couldn't copy the needed vector_eos_f90.dek file"
        sys.exit(1)

    # let's try building it
    print "Building the shared object file..."
    cmd = ["f2py", "-m", "fhelmEOS", "-c", "--quiet", helmEOS]
    try:
        subprocess.check_call(cmd,
                              stdout=open('/dev/null','w'),
                              stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print "Couldn't build the python wrapper for %s" % helmEOS
        sys.exit(1)

    # remove the un-needed local helmeos.f90 and include file
    try: 
        subprocess.check_call(["rm", helmEOS, vectorEOS])
    except subprocess.CalledProcessError:
        print "Couldn't cleanup the local %s and %s files!" % (
            helmEOS, vectorEOS)

    # import the module to make sure things are happy
    import fhelmEOS


# this just creates a symbolic link to the helm_table needed for an EOS call
def _link_to_helmTable():
    if not os.path.exists(helmTable):
        print "Grabbing helm_table.dat..."
        os.symlink(helmEOS_path+helmTable,helmTable)


# this wraps the _build_fhelmEOS method and checks to see if we already
# have a shared object file already built.  If we do, then no need to waste
# cycles building another one - unless of course you want to with rebuild=True
def fwrap(fParallel=None,rebuild=False):
    global _built

    if rebuild: _built = False

    if not _built: _build_fhelmEOS(fParallel)
    else: 
        print ""
        print "Great - wrapper already built."
        print "If you want to force a rebuild, call fwrap(rebuild=True)"
        print ""

    _link_to_helmTable()

    _built = True

# this class is essentially a data storage container for all the interesting 
# thermodynamic properties
# 
# a call to the eos will return an object of this type
class EOS_dataType(object):

    # initialize some important variables
    # depending on the input_type for the EOS, some of these won't be needed
    # as initial data, but will be calculated
    def __init__(self, 
                 input_type, 
                 abar, zbar,
                 den, temp,
                 h, p, e, s):

        self.input_type = input_type      # see EOS input modes at top
        self.abar = abar                  # average molar mass
        self.zbar = zbar                  # average charge

        self.den = den                    # density
        self.temp = temp                  # temperature
        self.h = h                        # specific enthalpy
        self.p = p                        # pressure
        self.e = e                        # specific internal energy
        self.s = s                        # specific entropy

    # we overwrite this method just for prettier printing
    def __str__(self):
        string = ""
        keys = self.__dict__.keys()
        keys.sort()
        for k in keys:
            string += "%s = %s\n" % (k,self.__dict__[k])

        return string


# this is the standard eos routine that calls the fortran helmeos code
#
#
def eos(input_type=None, 
        abar=None, zbar=None,
        den=ZERO, temp=ZERO,
        h=ZERO, p=ZERO, e=ZERO, s=ZERO):

    import fhelmEOS

    global _initialized

    # number of variables returned from helmeos
    nRetQuant = 23

    # initialize if this is the first call
    if not _initialized: 
        fhelmEOS.helmeos_init()
        _initialized = True

    # sanity checks
    if not input_type: 
        print "eos: Must specify input type!"
        sys.exit(1)
    if not abar or not zbar:
        print "eos: Must specify abar and zbar!"
        sys.exit(1)

    # build an EOS_dataType storage container
    myEOSData = EOS_dataType(input_type,abar,zbar,den,temp,h,p,e,s)

    # for now, this only works with the input_rt input_type
    # this will need fixin'
    if input_type != input_rt:
        print "ERROR: only support input_type = input_rt!"
        sys.exit(1)

    # call the EOS and store the output to a numpy array
    data = numpy.zeros(nRetQuant,numpy.float64)
    data[:] = fhelmEOS.helmeos(do_coulomb = True,
                               npts = 1,
                               temp_row = temp,
                               den_row = den,
                               abar_row = abar,
                               zbar_row = zbar)

    # data[0] = eosfail in fortran; check for errors
    if int(data[0]) != 0:
        print "Error in eos: eosfail = %s" % int(data[0])
        sys.exit(1)

    # parse the returned quantities from helmeos and store them in our
    # EOS_dataType object
    myEOSData.e = float(data[1])
    myEOSData.p = data[2]
    myEOSData.cv = data[3]
    myEOSData.cp = data[4]
    xne = data[5]
    xnp = data[6]
    myEOSData.eta = data[7]
    pele = data[8]
    ppos = data[9]
    myEOSData.dpdr = data[10]
    myEOSData.dpdT = data[11]
    myEOSData.dpda = data[12]
    myEOSData.dpdz = data[13]
    myEOSData.dedr = data[14]
    myEOSData.dedT = data[15]
    myEOSData.deda = data[16]
    myEOSData.dedz = data[17]
    myEOSData.gam1 = data[18]
    cs = data[19] # includes relativistic corrections, which we don't need
    myEOSData.s = data[20]
    myEOSData.dsdr = data[21]
    myEOSData.dsdT = data[22]

    # derived quantities
    myEOSData.ne = xne + xnp
    myEOSData.pele = pele + ppos
    myEOSData.h = myEOSData.e + myEOSData.p/myEOSData.den
    myEOSData.cs = numpy.sqrt(myEOSData.gam1*myEOSData.p/myEOSData.den)

    # return the EOS_dataType object
    return myEOSData


# a simple example using some sample data; see the test_helmeos.py script
# for a more thorough example
if __name__ == "__main__":
    fwrap(sys.argv[1])
    EOSData = eos(input_rt,
                  abar=1.91937636223,
                  zbar=1.41657199235,
                  den=691846.534,
                  temp=774872038.0)

    print EOSData
