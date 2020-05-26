import sys
import numpy as np
import vtk.util.numpy_support as vtknp
from vtk import vtkDataObject, vtkCompositeDataSet, vtkMultiBlockDataSet

# default values of control parameters
numBins = 10
meshName = ''
arrayName = ''
arrayCen = vtkDataObject.POINT
outFile = 'hist'
verbose = True

def Initialize():
    r = comm.Get_rank()
    if r == 0:
        if verbose:
            sys.stderr.write( \
                'Initialize numBins=%d meshName=%s arrayName=%s arrayCen=%d outFile=%s\n'%( \
                numBins, meshName, arrayName, arrayCen, outFile))

        # check for valid control parameters
        if not meshName:
            raise RuntimeError('meshName was not set')
        if not arrayName:
            raise RuntimeError('arrayName was not set')

def Execute(adaptor):
    r = comm.Get_rank()

    # get the mesh and array we need
    mesh = adaptor.GetMesh(meshName, True)
    adaptor.AddArray(mesh, meshName, arrayCen, arrayName)

    # force composite data to simplify computations
    if not isinstance(mesh, vtkCompositeDataSet):
        s = comm.Get_size()
        mb = vtkMultiBlockDataSet()
        mb.SetNumberOfBlocks(s)
        mb.SetBlock(r, mesh)
        mesh = mb

    # compute the min and max over local blocks
    mn = sys.float_info.max
    mx = -mn
    it = mesh.NewIterator()
    while not it.IsDoneWithTraversal():
        do = it.GetCurrentDataObject()

        atts = do.GetPointData() if arrayCen == vtkDataObject.POINT \
             else do.GetCellData()

        da = vtknp.vtk_to_numpy(atts.GetArray(arrayName))

        mn = min(mn, np.min(da))
        mx = max(mx, np.max(da))

        it.GoToNextItem()

    # compute global min and max
    mn = comm.allreduce(mn, op=MPI.MIN)
    mx = comm.allreduce(mx, op=MPI.MAX)

    # compute the histogram over local blocks
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        do = it.GetCurrentDataObject()

        atts = do.GetPointData() if arrayCen == vtkDataObject.POINT \
             else do.GetCellData()

        da = vtknp.vtk_to_numpy(atts.GetArray(arrayName))

        h,be = np.histogram(da, bins=numBins, range=(mn,mx))

        hist = hist + h if 'hist' in globals() else h

        it.GoToNextItem()

    # compute the global histogram on rank 0
    h = comm.reduce(hist, root=0, op=MPI.SUM)

    # rank 0 write to disk
    if r == 0:
        t = adaptor.GetDataTime()
        ts = adaptor.GetDataTimeStep()
        fn = '%s_%s_%s_%05d.txt'%(outFile, meshName, arrayName, ts)
        f = open(fn, 'w')
        f.write('step : %d\n'%(ts))
        f.write('time : %0.6g\n'%(t))
        f.write('num bins : %d\n'%(numBins))
        f.write('range : %0.6g %0.6g\n'%(mn, mx))
        f.write('bin edges : ')
        for v in be:
            f.write('%0.6g '%(v))
        f.write('\n')
        f.write('counts : ')
        for v in h:
            f.write('%d '%(v))
        f.write('\n')
        f.close()
        if verbose:
            sys.stderr.write('Execute "%s" written\n'%(fn))

def Finalize():
    r = comm.Get_rank()
    if r == 0 and verbose:
        sys.stderr.write('Finalize\n')
    return 0

