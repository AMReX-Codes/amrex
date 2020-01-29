# Copyright 2017 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""call back operations
=====================

These are the functions which allow installing user created functions so that
they are called at various places along the time step.

For each call back type, the following three functions are defined.
 - install___: Installs a function to be called at that specified time
 - uninstall___: Uninstalls the function (so it won't be called anymore)
 - isinstalled___: Checks if the function is installed

These functions all take a function or instance method as an argument. Note that
if an instance method is used, an extra reference to the method's object is saved.

The install can be done using a decorator, which has the prefix "callfrom". See example below.

Functions can be called at the following times:
 - afterinit <installafterinit>: immediately after the init is complete
 - beforeEsolve <installbeforeEsolve>: before the solve for E fields
 - afterEsolve <installafterEsolve>: after the solve for E fields
 - beforedeposition <installbeforedeposition>: before the particle deposition (for charge and/or current)
 - afterdeposition <installafterdeposition>: after particle deposition (for charge and/or current)
 - beforestep <installbeforestep>: before the time step
 - afterstep <installafterstep>: after the time step
 - particlescraper <installparticlescraper>: just after the particle boundary conditions are applied
                                             but before lost particles are processed
 - particleloader <installparticleloader>: at the time that the standard particle loader is called
 - particleinjection <installparticleinjection>: called when particle injection happens, after the position
                                                  advance and before deposition is called, allowing a user defined
                                                  particle distribution to be injected each time step
 - appliedfields <installappliedfields>: allows directly specifying any fields to be applied to the particles
                                         during the advance

To use a decorator, the syntax is as follows. This will install the function myplots to be called after each step.

@callfromafterstep
def myplots():
  ppzx()

This is equivalent to the following:

def myplots():
  ppzx()
installafterstep(myplots)

"""
from __future__ import generators

import types
import copy
import time
import ctypes
import sys
import numpy

from ._libwarpx import libwarpx

class CallbackFunctions(object):
    """
    Class to handle call back function lists.

    Note that for functions passed in that are methods of a class instance,
    a full reference of the instance is saved. This extra reference means
    that the object will not actually deleted if the user deletes the
    original reference. This is good since the user does not need to keep
    the reference to the object (for example it can be created using a local
    variable in a function). It may be bad if the user thinks an object was
    deleted, but it actually isn't since it had (unkown to the user)
    installed a method in one of the call back lists
    """

    def __init__(self,name=None,lcallonce=0):
        self.funcs = []
        self.time = 0.
        self.timers = {}
        self.name = name
        self.lcallonce = lcallonce

    def __call__(self,*args,**kw):
        "Call all of the functions in the list"
        tt = self.callfuncsinlist(*args,**kw)
        self.time = self.time + tt
        if self.lcallonce: self.funcs = []

    def clearlist(self):
        self.funcs = []

    def __nonzero__(self):
        "Returns True if functions are installed, otherwise False"
        return self.hasfuncsinstalled()

    def __len__(self):
        "Returns number of functions installed"
        return len(self.funcs)

    def hasfuncsinstalled(self):
        "Checks if there are any functions installed"
        return len(self.funcs) > 0

    def _getmethodobject(self,func):
        "For call backs that are methods, returns the method's instance"
        return func[0]

    def callbackfunclist(self):
        "Generator returning callable functions from the list"
        funclistcopy = copy.copy(self.funcs)
        for f in funclistcopy:
            if isinstance(f,list):
                object = self._getmethodobject(f)
                if object is None:
                    self.funcs.remove(f)
                    continue
                result = getattr(object,f[1])
            elif isinstance(f,basestring):
                import __main__
                if f in __main__.__dict__:
                    result = __main__.__dict__[f]
                    # --- If the function with the name is found, then replace the
                    # --- name in the list with the function.
                    self.funcs[self.funcs.index(f)] = result
                else:
                    continue
            else:
                result = f
            if not callable(result):
                print("\n\nWarning: a call back was found that is not callable.")
                if self.name is not None:
                    print("For %s"%self.name)
                print("Only callable objects can be installed.")
                print("It is possible that the callable's name has been overwritten")
                print("by something not callable. This can happen during restart")
                print("if a function name had later been used as a variable name.")
                print(self.name)
                if isinstance(f,basestring):
                    print("The name of the call back is %s"%f)
                print("\n\n")
                continue
            yield result

    def installfuncinlist(self,f):
        "Check if the specified function is installed"
        if isinstance(f,types.MethodType):
            # --- If the function is a method of a class instance, then save a full
            # --- reference to that instance and the method name.
            finstance = f.im_self
            fname = f.__name__
            self.funcs.append([finstance,fname])
        elif callable(f):
            # --- If a function had already been installed by name, then skip the install.
            # --- This is problematic, since no warning message is given, but it is unlikely
            # --- to arise under normal circumstances.
            # --- The purpose of this check is to avoid redundant installation of functions
            # --- during a restore from a dump file. Without the check, functions that had been
            # --- installed via a decorator would be installed an extra time since the source
            # --- of the function contains the decoration (which is activated when the source
            # --- is exec'd).
            if f.__name__ not in self.funcs:
                self.funcs.append(f)
        else:
            self.funcs.append(f)

    def uninstallfuncinlist(self,f):
        "Uninstall the specified function"
        # --- An element by element search is needed
        # --- f can be a function or method object, or a name (string).
        # --- Note that method objects can not be removed by name.
        funclistcopy = copy.copy(self.funcs)
        for func in funclistcopy:
            if f == func:
                self.funcs.remove(f)
                return
            elif isinstance(func,list) and isinstance(f,types.MethodType):
                object = self._getmethodobject(func)
                if f.im_self is object and f.__name__ == func[1]:
                    self.funcs.remove(func)
                    return
            elif isinstance(func,basestring):
                if f.__name__ == func:
                    self.funcs.remove(func)
                    return
            elif isinstance(f,basestring):
                if isinstance(func,basestring): funcname = func
                elif isinstance(func,list): funcname = None
                else:                        funcname = func.__name__
                if f == funcname:
                    self.funcs.remove(func)
                    return
        raise Exception('Warning: no such function had been installed')

    def isinstalledfuncinlist(self,f):
        "Checks if the specified function is installed"
        # --- An element by element search is needed
        funclistcopy = copy.copy(self.funcs)
        for func in funclistcopy:
            if f == func:
                return 1
            elif isinstance(func,list) and isinstance(f,types.MethodType):
                object = self._getmethodobject(func)
                if f.im_self is object and f.__name__ == func[1]:
                    return 1
            elif isinstance(func,basestring):
                if f.__name__ == func:
                    return 1
        return 0

    def callfuncsinlist(self,*args,**kw):
        "Call the functions in the list"
        bb = time.time()
        for f in self.callbackfunclist():
            #barrier()
            t1 = time.time()
            f(*args,**kw)
            #barrier()
            t2 = time.time()
            # --- For the timers, use the function (or method) name as the key.
            self.timers[f.__name__] = self.timers.get(f.__name__,0.) + (t2 - t1)
        aa = time.time()
        return aa - bb

#=============================================================================

# --- Now create the actual instances.
_afterinit = CallbackFunctions('afterinit')
_beforeEsolve = CallbackFunctions('beforeEsolve')
_afterEsolve = CallbackFunctions('afterEsolve')
_beforedeposition = CallbackFunctions('beforedeposition')
_afterdeposition = CallbackFunctions('afterdeposition')
_particlescraper = CallbackFunctions('particlescraper')
_particleloader = CallbackFunctions('particleloader')
_beforestep = CallbackFunctions('beforestep')
_afterstep = CallbackFunctions('afterstep')
_afterrestart = CallbackFunctions('afterrestart',lcallonce=1)
_particleinjection = CallbackFunctions('particleinjection')
_appliedfields = CallbackFunctions('appliedfields')

# --- Create the objects that can be called from C.
# --- Note that each of the CFUNCTYPE instances need to be saved
_CALLBACK_FUNC_0 = ctypes.CFUNCTYPE(None)
_c_afterinit = _CALLBACK_FUNC_0(_afterinit)
libwarpx.warpx_set_callback_py_afterinit(_c_afterinit)
_c_beforeEsolve = _CALLBACK_FUNC_0(_beforeEsolve)
libwarpx.warpx_set_callback_py_beforeEsolve(_c_beforeEsolve)
_c_afterEsolve = _CALLBACK_FUNC_0(_afterEsolve)
libwarpx.warpx_set_callback_py_afterEsolve(_c_afterEsolve)
_c_beforedeposition = _CALLBACK_FUNC_0(_beforedeposition)
libwarpx.warpx_set_callback_py_beforedeposition(_c_beforedeposition)
_c_afterdeposition = _CALLBACK_FUNC_0(_afterdeposition)
libwarpx.warpx_set_callback_py_afterdeposition(_c_afterdeposition)
_c_particlescraper = _CALLBACK_FUNC_0(_particlescraper)
libwarpx.warpx_set_callback_py_particlescraper(_c_particlescraper)
_c_particleloader = _CALLBACK_FUNC_0(_particleloader)
libwarpx.warpx_set_callback_py_particleloader(_c_particleloader)
_c_beforestep = _CALLBACK_FUNC_0(_beforestep)
libwarpx.warpx_set_callback_py_beforestep(_c_beforestep)
_c_afterstep = _CALLBACK_FUNC_0(_afterstep)
libwarpx.warpx_set_callback_py_afterstep(_c_afterstep)
_c_afterrestart = _CALLBACK_FUNC_0(_afterrestart)
libwarpx.warpx_set_callback_py_afterrestart(_c_afterrestart)
_c_particleinjection = _CALLBACK_FUNC_0(_particleinjection)
libwarpx.warpx_set_callback_py_particleinjection(_c_particleinjection)
_c_appliedfields = _CALLBACK_FUNC_0(_appliedfields)
libwarpx.warpx_set_callback_py_appliedfields(_c_appliedfields)

#=============================================================================
def printcallbacktimers(tmin=1.,lminmax=False,ff=None):
    """Prints timings of installed functions.
    - tmin=1.: only functions with time greater than tmin will be printed
    - lminmax=False: If True, prints the min and max times over all processors
    - ff=None: If given, timings will be written to the file object instead of stdout
    """
    if ff is None: ff = sys.stdout
    for c in [_afterinit,_beforeEsolve,_afterEsolve,
              _beforedeposition,_afterdeposition,
              _particlescraper,
              _particleloader,
              _beforestep,_afterstep,
              _afterrestart,
              _particleinjection,
              _appliedfields]:
        for fname, time in c.timers.items():
            #vlist = numpy.array(gather(time))
            vlist = numpy.array([time])
            #if me > 0: continue
            vsum = numpy.sum(vlist)
            if vsum <= tmin: continue
            vrms = numpy.sqrt(max(0.,numpy.sum(vlist**2)/len(vlist) - (numpy.sum(vlist)/len(vlist))**2))
            npes = 1. # Only works for one processor
            ff.write('%20s %s %10.4f  %10.4f %10.4f'%(c.name,fname,vsum,vsum/npes,vrms))
            if lminmax:
                vmin = numpy.min(vlist)
                vmax = numpy.max(vlist)
                ff.write('  %10.4f  %10.4f'%(vmin,vmax))
            it = libwarpx.warpx_getistep(0)
            if it > 0:
                ff.write('   %10.4f'%(vsum/npes/(it)))
            ff.write('\n')

#=============================================================================
# ----------------------------------------------------------------------------
def callfromafterinit(f):
    installafterinit(f)
    return f
def installafterinit(f):
    "Adds a function to the list of functions called after the init"
    _afterinit.installfuncinlist(f)
def uninstallafterinit(f):
    "Removes the function from the list of functions called after the init"
    _afterinit.uninstallfuncinlist(f)
def isinstalledafterinit(f):
    "Checks if the function is called after a init"
    return _afterinit.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfrombeforeEsolve(f):
    installbeforeEsolve(f)
    return f
def installbeforeEsolve(f):
    "Adds a function to the list of functions called before an E solve"
    _beforeEsolve.installfuncinlist(f)
def uninstallbeforeEsolve(f):
    "Removes the function from the list of functions called before an E solve"
    _beforeEsolve.uninstallfuncinlist(f)
def isinstalledbeforeEsolve(f):
    "Checks if the function is called before an E solve"
    return _beforeEsolve.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterEsolve(f):
    installafterEsolve(f)
    return f
def installafterEsolve(f):
    "Adds a function to the list of functions called after an E solve"
    _afterEsolve.installfuncinlist(f)
def uninstallafterEsolve(f):
    "Removes the function from the list of functions called after an E solve"
    _afterEsolve.uninstallfuncinlist(f)
def isinstalledafterEsolve(f):
    "Checks if the function is called after an E solve"
    return _afterEsolve.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfrombeforedeposition(f):
    installbeforedeposition(f)
    return f
def installbeforedeposition(f):
    "Adds a function to the list of functions called before a particle deposition"
    _beforedeposition.installfuncinlist(f)
def uninstallbeforedeposition(f):
    "Removes the function from the list of functions called before a particle deposition"
    _beforedeposition.uninstallfuncinlist(f)
def isinstalledbeforedeposition(f):
    "Checks if the function is called before a particle deposition"
    return _beforedeposition.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterdeposition(f):
    installafterdeposition(f)
    return f
def installafterdeposition(f):
    "Adds a function to the list of functions called after a particle deposition"
    _afterdeposition.installfuncinlist(f)
def uninstallafterdeposition(f):
    "Removes the function from the list of functions called after a particle deposition"
    _afterdeposition.uninstallfuncinlist(f)
def isinstalledafterdeposition(f):
    "Checks if the function is called after a particle deposition"
    return _afterdeposition.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromparticlescraper(f):
    installparticlescraper(f)
    return f
def installparticlescraper(f):
    "Adds a function to the list of functions called to scrape particles"
    _particlescraper.installfuncinlist(f)
def uninstallparticlescraper(f):
    "Removes the function from the list of functions called to scrape particles"
    _particlescraper.uninstallfuncinlist(f)
def isinstalledparticlescraper(f):
    "Checks if the function is called to scrape particles"
    return _particlescraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromparticleloader(f):
    installparticleloader(f)
    return f
def installparticleloader(f):
    "Adds a function to the list of functions called to load particles"
    _particleloader.installfuncinlist(f)
def uninstallparticleloader(f):
    "Removes the function from the list of functions called to load particles"
    _particleloader.uninstallfuncinlist(f)
def isinstalledparticleloader(f):
    "Checks if the function is called to load particles"
    return _particleloader.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfrombeforestep(f):
    installbeforestep(f)
    return f
def installbeforestep(f):
    "Adds a function to the list of functions called before a step"
    _beforestep.installfuncinlist(f)
def uninstallbeforestep(f):
    "Removes the function from the list of functions called before a step"
    _beforestep.uninstallfuncinlist(f)
def isinstalledbeforestep(f):
    "Checks if the function is called before a step"
    return _beforestep.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterstep(f):
    installafterstep(f)
    return f
def installafterstep(f):
    "Adds a function to the list of functions called after a step"
    _afterstep.installfuncinlist(f)
def uninstallafterstep(f):
    "Removes the function from the list of functions called after a step"
    _afterstep.uninstallfuncinlist(f)
def isinstalledafterstep(f):
    "Checks if the function is called after a step"
    return _afterstep.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterrestart(f):
    raise Exception('restart call back not implemented yet')
    installafterrestart(f)
    return f
def installafterrestart(f):
    "Adds a function to the list of functions called immediately after a restart"
    raise Exception('restart call back not implemented yet')
    _afterrestart.installfuncinlist(f)
def uninstallafterrestart(f):
    "Removes the function from the list of functions called immediately after a restart"
    raise Exception('restart call back not implemented yet')
    _afterrestart.uninstallfuncinlist(f)
def isinstalledafterrestart(f):
    "Checks if the function is called immediately after a restart"
    raise Exception('restart call back not implemented yet')
    return _afterrestart.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromparticleinjection(f):
    installparticleinjection(f)
    return f
def installparticleinjection(f):
    """
    Adds a user defined function that is to be called when particle
    injection happens, after the position advance and before deposition is
    called, allowing a user defined particle distribution to be injected
    each time step"""
    _particleinjection.installfuncinlist(f)
def uninstallparticleinjection(f):
    "Removes the function installed by installparticleinjection"
    _particleinjection.uninstallfuncinlist(f)
def isinstalledparticleinjection(f):
    "Checks if the function is called when particles injection happens"
    return _particleinjection.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromappliedfields(f):
    raise Exception('applied fields call back not implemented yet')
    installappliedfields(f)
    return f
def installappliedfields(f):
    """
    Adds a user defined function which can specify E and B fields which are applied
    to the particles during the particle advance.
    """
    raise Exception('applied fields call back not implemented yet')
    _appliedfields.installfuncinlist(f)
def uninstallappliedfields(f):
    "Removes the function installed by installappliedfields"
    raise Exception('applied fields call back not implemented yet')
    _appliedfields.uninstallfuncinlist(f)
def isinstalledappliedfields(f):
    "Checks if the function is called when which applies fields"
    raise Exception('applied fields call back not implemented yet')
    return _appliedfields.isinstalledfuncinlist(f)

