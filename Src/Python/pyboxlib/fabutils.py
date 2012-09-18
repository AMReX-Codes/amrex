"""Various Fabric (fabfile.org) utilities for launching and managing
BoxLib runs.

These routines assume

"""

import glob
import os

from fabric.api import *
from fabric.colors import *
from fabric.utils import *
from fabric.contrib.files import *

from utils import *

env.mpirun   = 'mpirun'
env.nprocs   = 1
env.nthreads = 1


def runexe(name, rundir, probin):

    if rundir:
        rundir = os.path.join(env.bin, rundir)
    else:
        rundir = env.bin

    with cd(rundir):

        if exists('out'):
            puts(red('skipping: ' + name))
            puts(red(indent("{rundir}/out already exists".format(rundir=rundir))))
            return

        puts(green('running: ' + name))
        puts(green(indent('exe: ' + env.bin + '/' + env.exe)))
        puts(green(indent('cwd: ' + rundir)))
        puts(green(indent('probin: ' + probin)))

        run('{mpirun} -n {nprocs} -env OMP_NUM_THREADS {nthreads} {bin}/{exe} {probin} 2> err | tee out'.format(
                mpirun=env.mpirun, nprocs=env.nprocs, nthreads=env.nthreads, 
                bin=env.bin, exe=env.exe, probin=probin))


class stage(object):

    def __init__(self, stage='staging'):
        self.stage = stage

    def __enter__(self):
        self.mkdir()
        return self

    def __exit__(self, type, value, tb):
        import shutil
        puts(green('syncing stage'))
        local('rsync -auz {src}/* {host}:{dst}'.format(src=self.stage, host=env.host, dst=env.bin))
        shutil.rmtree(self.stage)

    def mkdir(self, *dirnames):
        path = os.path.join(self.stage, *dirnames)
        os.makedirs(path)
        return path + os.sep


@task            
def rsync():

  if env.host == 'localhost':
    return

  exclude = [ '*~', 't', '.git', '*.exe', 'plt*', '*.pyc', '*.so', 'staging' ]
  exclude = map(lambda x: "'" + x + "'", exclude)
  exclude = ' --exclude '.join(exclude)

  for d in env.rsync_dirs:
      command = "rsync -aPz --exclude {exclude} {src}/ {host}:{dst}".format(
          exclude=exclude, host=env.host, src=d, dst=d)
      local(command)


@task
def make(target=''):

  with cd(env.bin):
    run('make %s' % target)
