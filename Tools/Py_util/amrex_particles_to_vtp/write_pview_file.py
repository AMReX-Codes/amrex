import numpy as np

def write_paraview_file_structmesh(fname,xp,yp,ccdata,ncdata):

    outfile=open(fname,'w')

    Npx=xp.shape[1]
    Npy=xp.shape[0]
    zero=0
    one=1

    outfile.write("<?xml version=\"1.0\"?>\n")
    outfile.write("<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
    outfile.write("<StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n"%(one,Npx,one,Npy,one,one))

    #piece extent is useful when writing in parallel
    #it corresponds to extents for each processor in global coordinates, eg: 1-2 is proc 0, 3-4 is proc 0
    outfile.write("<Piece Extent=\"%d %d %d %d %d %d\">\n"%(one,Npx,one,Npy,one,one))

    outfile.write("<PointData>\n")
    n_ncdata=ncdata.shape[0]
    if(n_ncdata > 0):
        for ndataset in range(n_ncdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"node_data%d\" format=\"ascii\">\n"%(ndataset))
            for i in range(ncdata.shape[1]):
                outfile.write("%e\t"%(ncdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</PointData>\n")

    outfile.write("<CellData>\n")
    n_ccdata=ccdata.shape[0]
    if(n_ccdata > 0):
        for ndataset in range(n_ccdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"cell_data%d\" format=\"ascii\">\n"%(ndataset))
            for i in range(ccdata.shape[1]):
                outfile.write("%e\t"%(ccdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</CellData>\n")

    outfile.write("<Points>\n")

    outfile.write("<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")
    for j in range(Npy):
        for i in range(Npx):
            outfile.write("%e\t%e\t%e\t"%(xp[j][i],yp[j][i],zero))
    outfile.write("\n</DataArray>\n")

    outfile.write("</Points>\n")
    outfile.write("</Piece>\n")
    outfile.write("</StructuredGrid>\n")
    outfile.write("</VTKFile>")

def write_paraview_file_unst_trimesh(fname,pts,conn,ccdata,ncdata):

    outfile=open(fname,'w')

    Npts=pts.shape[0]
    Ntri=conn.shape[0]
    zero=0
    one=1

    outfile.write("<?xml version=\"1.0\"?>\n")
    outfile.write("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
    outfile.write("<UnstructuredGrid>\n")
    outfile.write("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n"%(Npts,Ntri))

    outfile.write("<PointData>\n")
    n_ncdata=ncdata.shape[0]
    if(n_ncdata > 0):
        for ndataset in range(n_ncdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"Point_data%d\" format=\"ascii\">\n"%(ndataset))
            for i in range(ncdata.shape[1]):
                outfile.write("%e "%(ncdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</PointData>\n")

    outfile.write("<CellData>\n")
    n_ccdata=ccdata.shape[0]
    if(n_ccdata > 0):
        for ndataset in range(n_ccdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"Cell_data%d\" format=\"ascii\">\n"%(ndataset))
            for i in range(ccdata.shape[1]):
                outfile.write("%e "%(ccdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</CellData>\n")

    outfile.write("<Points>\n")
    outfile.write("<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n")
    for i in range(Npts):
        outfile.write("%e\t%e\t%e\t"%(pts[i][0],pts[i][1],zero))
    outfile.write("\n</DataArray>\n")
    outfile.write("</Points>\n")

    outfile.write("<Cells>\n")
    outfile.write("<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n")
    for i in range(Ntri):
        for j in range(3):
            outfile.write("%d "%(conn[i][j]))
    outfile.write("\n</DataArray>\n")

    outfile.write("<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n")
    for i in range(Ntri):
        offs=3*(i+1)
        outfile.write("%d "%(offs))
    outfile.write("\n</DataArray>\n")
    outfile.write("<DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n")
    for i in range(Ntri):
        tri_type=5
        outfile.write("%d "%(tri_type))
    outfile.write("\n</DataArray>\n")

#note, different types are
#triangle - 5
#quad     - 9
#tetrahedron - 10
#hexahedron  - 12
#prism       - 13
#pyramid     - 14

    outfile.write("</Cells>\n")
    outfile.write("</Piece>\n")
    outfile.write("</UnstructuredGrid>\n")
    outfile.write("</VTKFile>\n")

    outfile.close()

def write_paraview_file_cartmesh(fname,dx,prob_lo,N,ncdata,ccdata):
    
    zero=0
    one=1
    outfile=open(fname,'w')
    outfile.write("<?xml version=\"1.0\"?>\n")
    outfile.write("<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
    outfile.write("<RectilinearGrid WholeExtent=\"%d\t%d\t%d\t%d\t%d\t%d\">\n"%(one,N[0],one,N[1],one,one))
    outfile.write("<Piece Extent=\"%d\t%d\t%d\t%d\t%d\t%d\">\n"%(one,N[0],one,N[1],one,one))

    outfile.write("<PointData>\n")
    n_ncdata=ncdata.shape[0]
    if(n_ncdata > 0):
        for ndataset in range(n_ncdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"Point_data%d\" format=\"ascii\">\n"%(ndataset))

            for i in range(ncdata.shape[1]):
                outfile.write("%e "%(ncdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</PointData>\n")

    outfile.write("<CellData>\n")
    n_ccdata=ccdata.shape[0]
    if(n_ccdata > 0):
        for ndataset in range(n_ccdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"Cell_data%d\" format=\"ascii\">\n"%(ndataset))
            for i in range(ccdata.shape[1]):
                outfile.write("%e "%(ccdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</CellData>\n")

    outfile.write("<Coordinates>\n")

    outfile.write("<DataArray type=\"Float32\" Name=\"X\"  format=\"ascii\">\n")
    for i in range(N[0]):
        outfile.write("%e\t"%(prob_lo[0]+i*dx[0]))
    outfile.write("\n</DataArray>\n")

    outfile.write("<DataArray type=\"Float32\" Name=\"Y\"  format=\"ascii\">\n")
    for i in range(N[1]):
        outfile.write("%e\t"%(prob_lo[1]+i*dx[1]))
    outfile.write("\n</DataArray>\n")
    
    outfile.write("<DataArray type=\"Float32\" Name=\"Z\"  format=\"ascii\">\n")
    outfile.write("%e\t"%(0.0))
    outfile.write("\n</DataArray>\n")
    
    outfile.write("</Coordinates>\n")
    outfile.write("</Piece>\n")
    outfile.write("</RectilinearGrid>\n")
    outfile.write("</VTKFile>")

    outfile.close()


def write_paraview_file_particles(fname,pts,ncdata):

    outfile=open(fname,'w')

    Npts=pts.shape[0]
    zero=0
    one=1

    outfile.write("<?xml version=\"1.0\"?>\n")
    outfile.write("<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
    outfile.write("<PolyData>\n")
    outfile.write("<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"%(Npts))

    outfile.write("<PointData>\n")
    n_ncdata=ncdata.shape[0]
    if(n_ncdata > 0):
        for ndataset in range(n_ncdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"Point_data%d\" format=\"ascii\">\n"%(ndataset))
            for i in range(ncdata.shape[1]):
                outfile.write("%e "%(ncdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</PointData>\n")

    outfile.write("<Points>\n")
    outfile.write("<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n")
    for i in range(Npts):
        outfile.write("%e\t%e\t%e\t"%(pts[i][0],pts[i][1],pts[i][2]))
    outfile.write("\n</DataArray>\n")
    outfile.write("</Points>\n")

    outfile.write("</Piece>\n")
    outfile.write("</PolyData>\n")
    outfile.write("</VTKFile>\n")

    outfile.close()
