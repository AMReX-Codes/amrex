"""
 a module for parsing MAESTRO particle files

  Particle histories in MAESTRO are stored in a set of files, e.g.
  timestamp_#, where the # in the file name is the number of the CPU
  (actually MPI process) with particles.  As particles advect, they
  may jump to different CPUs, so a given particle's history may
  initially be written out in one file and then switch to a different
  file.  This module sorts the particle data to be in proper time
  order.

  This module should be imported in your own python script, for example:

    >>> import parseparticles
    >>> myparticles = parseparticles.parseParticleFile(myParticleFiles)

    where myParticleFiles is a list of the timestamp files containing
    the particle data.

    This returns a dictionary whose keys are internal ids (comprised
    of the particle's id and its origin CPU) and values are particle
    objects, which contain all of the associated information about a
    particle's history


  Usually we don't care about the original particle ID information, so
  we just strip out the particle objects into a list:

    >>> particles = myparticles.values()

    Here, each item in the list is a particle object.  


  Each particle object contains data describing the particle
  identity, as well as a history.  The history itself is a list of
  particleInstance objects.  Each particleInstance object stores
  the state of the particle at an instance in time.  For example:

    >>> time = particles[n].history[i].t 

    is the time associated with particleInstance i of particle n.
    There is a unique particleInstance for particle n for every
    timestep of this particle's existence.

    Coordinate information is accessed as:

    >>> x = particles[n].history[i].xyz[0]

    for the x-coordinate (xyz[1] would be y, etc.).

    State data is stored in a numpy array.  The state names
    corresponding to each data element are stored in the main particle
    object and can be looked-up via the getVarIndex() method.  For
    example, to store the density of particle n at timestep i into the
    dens variable, we would do:

    >>> idens = particles[n].getVarIndex("density")
    >>> dens = particles[n].history[i].data[idens]

    Any field that was stored in the original timestamp file would be
    present in the data array (and the name should be available
    via getVarIndex())


  The total number of particles can be found simply as len(particles),
  and the number of particleInstances for a given particle can be
  found as len(particles[n].history)
""" 

import sys
import string
import numpy


# a global dictionary to store the particle objects.  Each unique
# particle stored with its own key:value pair -- the key is the unique
# particle ID and the value is the particle object itself.  The
# particle objects in turn contain many particleInstance objects in
# the history.
_particleDict = {}


# this works for < 1000 MPI processes, but should ultimately be something 
# dynamic to work for all cases
_idFactor = 1000

#------------------------------------------------------------------------
class particleInstance(object):
    """
       a particleInstance object represents a single instance of a
       particle, and will carry the position, time, and associated
       data.
    """
    __slots__ = ['xyz', 't', 'data']    # memory management

    def __init__(self, xyz, t, data):
        self.xyz = xyz
        self.t = t
        self.data = data



    def __str__(self):
        """
           string representation of the particleInstance for printing
        """
        string = "particle pos:  (%g, %g, %g) \n" %  \
            (self.xyz[0], self.xyz[1], self.xyz[2]) + \
            "         time: %g \n" % (self.t) 
        
        return string

    
    
    def value(self):
        """
           return the value of a particleInstance for comparison 
           purposes.  The value is simply the time.
        """
        return self.t

    def __cmp__(self, other):
        return cmp(self.value(), other.value())



#------------------------------------------------------------------------
class particle(object):
    """
       a particle object stores a history of a single particle, each
       element of which is a particleInstance object.
    """

    __slots__ = ["pid", "originCPU", "dim", "finalized", 
                 "history", "dataNames", "numInstances"]

    def __new__(cls, pid=None, originCPU=None, *args, **kwargs):
        """
           Create a new particle instance in the global _particleDict
           dictionary.
        """

        # before building a new particle, make sure it isn't already
        # built by checking the global _particleDict
        id = pid*_idFactor + originCPU
        if id not in _particleDict:
            obj = object.__new__(cls)
            _particleDict[id] = obj
        return _particleDict[id]

    

    def __init__(self, pid, originCPU, dim, dataNames):
        """
           initialize a particle object
        """
        
        # a MAESTRO particle is identified by 2 numbers, the pid and
        # the CPU that it was created on.  Together, these uniquely
        # identify the particle.
        self.pid = pid
        self.originCPU = originCPU

        # dimensionality of particle data
        self.dim = dim

        # finalized is 1 when we have finished adding data and sorted
        # the history in time-order
        self.finalized = 0   

        # the history list will store instances of the particle at
        # different times.
        self.history = []

        # the dataNames will store the names of the data associated
        # with each particle instance
        self.dataNames = list(dataNames)

        # keep track of the number of particle history instances we've
        # stored
        self.numInstances = 0
    


    def addInstance(self, xyz=[-1.0,-1.0,-1.0],t=0.0, 
                    dataValues=[]):
        """
           add a particleInstance object to the particle history to
           store a single instance (in time) of this particle
        """

        if (self.finalized == 1):
            print "ERROR: particle data finalized before addInstance call"
            sys.exit(2)


        xyz = numpy.array(xyz)
        dataValues = numpy.array(dataValues)


        # add this particle instance to the particle history
        self.history.append(particleInstance(xyz,t,dataValues))
            
        self.numInstances += 1



    def finalize(self):
        """
           sort the particle histories in time order (since they may
           have been spread across multiple processors / files.
  
           finalize() should only be called after all particle data
           has been read in
        """

        # sort the history by time -- the particleInstance objects use
        # time as their value for comparison, so sort() should handle
        # this nicely for us.
        self.history.sort()

        self.finalized = 1



    def flatten(self):
        """
           return numpy arrays with the coordinates and time for all
           instances of this particle
        """

        if (not self.finalized == 1):
            print "ERROR: particle data not finalized before flatten call"
            sys.exit(2)

        coords = numpy.zeros((self.dim, self.numInstances), dtype=numpy.float64)
        time   = numpy.zeros((self.numInstances), dtype=numpy.float64)

        n = 0
        while (n < len(self.history)):
            coords[:,n] = self.history[n].xyz[:self.dim]
            
            time[n] = self.history[n].t

            n += 1

        return coords, time


    def getVarIndex(self, varname):
        """
           returns the index into the history[].data array corresponding
           to the variable varname
        """

        index = self.dataNames.index(varname)
        
        return index



    def __str__(self):
        string = "particle ID:           %d \n" % (self.pid) + \
                 "particle origin CPU:   %d \n" % (self.originCPU) + \
                 "dimensionality:        %d \n" % (self.dim) + \
                 "num. stored instances: %d \n" % (self.numInstances) + \
                 "associated data: \n"

        for item in self.dataNames:
            string += "  %s" % (item)
        
        string += "\n"

        return string


#------------------------------------------------------------------------
def ensure_list(obj):
    """
       make sure an object is a list; used, for example, when a single
       file is passed into the parseParticleFile method
    """

    import types

    if obj == None:
        return [obj]
    if not isinstance(obj, types.ListType):
        return [obj]
    return obj



#-----------------------------------------------------------------------------
def parseParticleFile(maestroParticleFiles):
    """
       read in all the particle data from a file and add each particle
       instance to the _particleDict dictionary, grouping all of the
       history for a unique particle (identified by the pid and
       originCPU) together in a particle object in the dictionary
    """

    haveHeader = 0

    for maestroParticleFile in ensure_list(maestroParticleFiles):

        # read the file line by line
        mf = open(maestroParticleFile, "r")

        for line in mf:

            # skip blank lines
            if line.isspace(): continue

            # look for a header information
            if (line.startswith("#")):

                # store old header, if it exists
                if (haveHeader == 1):
                    oldHeader = list(dataNames)  # list() makes a copy

                
                fields = string.split(line[1:])

                # make sure we know what we are doing -- the first 2
                # fields should be the particle ID and origin CPU
                if (fields[0] == "part-ID" and fields[1] == "origin-CPU"):
                    ipid = 0
                    ioCPU = 1

                else:
                    print "ERROR: particle file columns not in expected order"
                    sys.exit(2)


                # the next fields should be x, y, and z, depending on the
                # dimensionality
                if (fields[2] == "x" and fields[3] == "y" and 
                    fields[4] == "z"):
                    dim = 3
                    ix = 2; iy = 3; iz = 4
                
                elif (fields[2] == "x" and fields[3] == "y"):
                    dim = 2
                    ix = 2; iy = 3

                elif (fields[2] == "x"):
                    dim = 1
                    ix = 2

                else:
                    print "ERROR: particle file columns not in expected order"
                    sys.exit(2)


                # then comes time
                if (fields[2 + dim] == "time"):
                    it = 2 + dim
                
                else:
                    print "ERROR: particle file columns not in expected order"
                    sys.exit(2)
            

                # everything else is associated data
                if (len(fields) > 3 + dim):
                    idata = 3 + dim
                    ndata = len(fields) - idata
                    dataNames = fields[idata:]
                else:
                    ndata = 0


                # did the header change?
                if (haveHeader == 1):
                    if (not oldHeader == dataNames):
                        print "ERROR: header changed while parsing files"
                        sys.exit(2)


                # done with the header
                haveHeader = 1
                numData = len(dataNames)
                continue


            # normal particle data -- first find the required particle data
            fields = map(float,string.split(line))

            pid = int(fields[ipid])
            originCPU = int(fields[ioCPU])

            id = pid*_idFactor + originCPU

            xyz = fields[ix:ix+dim]
            time = fields[it]
            dataValues = fields[idata:]

            if (not len(dataValues) == numData):
                print "ERROR: number of data values not equal to number of names"
                sys.exit(2)


            if not id in _particleDict: 
                particle(pid, originCPU, dim, dataNames)

            _particleDict[id].addInstance(xyz=xyz,
                                          t=time,
                                          dataValues=dataValues)


    # finish up by making sure that each particle's history is
    # in proper time-order
    for id in _particleDict:
        _particleDict[id].finalize()


    return _particleDict
