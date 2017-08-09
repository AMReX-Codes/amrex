#include <iostream>
#include <string.h>

#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H"

const int BLOCKZ=4;
const int BLOCKY=8;
typedef double *** Array3D;

Array3D Array3DMalloc(int nx, int ny, int nz)
{
   Array3D U = (Array3D)malloc(sizeof(double**)*nz + sizeof(double *)*ny*nz +sizeof(double)*nx*ny*nz);

   for(int k=0;k<nz;k++){
       U[k] = ((double**) U) + nz + k*ny;
    }
   double *Ustart = (double *) (U[nz-1] + ny);
   for(int k=0;k<nz;k++)
       for(int j=0;j<ny;j++)
           U[k][j] = Ustart + k*nx*ny + j*nx;
  return U;
}


using namespace amrex;

class Jacobi :public Task{
    private:
	int rankx, ranky, rankz;
	Array3D U, Un, b;
	int iteration;
	double c1;
	double c2;
    public:
	static int nIters;
	static int tx, ty, tz;
	static int nx, ny, nz;
	Jacobi(){
	    iteration=-1;
	    double c1=1.0/6.0;
	    double c2=1.0;
	}

	void initializeData(){
            U = Array3DMalloc(nx, ny, nz);
            Un = Array3DMalloc(nx, ny, nz);
            b = Array3DMalloc(nx, ny, nz);
	    for (int k=0; k<nz; k++)
		for (int j=0; j<ny; j++)
		    for(int i=0; i<nx; i++){
			if(((rankz==0)&&(k==0)) ||((rankz==tz-1)&&(k==nz-1)) ||((ranky==0)&&(j==0)) ||\
				((ranky==ty-1)&&(j==ny-1))||((rankx==0)&&(i==0)) ||((rankx==tx-1)&&(i==nx-1))){
			    U[k][j][i] = 0;
			    Un[k][j][i] = 0;
			}else{
			    U[k][j][i] = 1;
			    Un[k][j][i] = 1;
			}
			b[k][j][i] = 0;
		    }
	}


	void Job(){
	    TaskName taskID= MyName();
	    rankx= taskID[0];
	    ranky= taskID[1];
	    rankz= taskID[2];
	    if(iteration==-1){
		initializeData();
	    }else if(iteration<nIters){
		for (int k0 = 1; k0 < nz - 1; k0+=BLOCKZ) {
		    int k1= k0+BLOCKZ<nz-1?k0+BLOCKZ:nz-1;
		    for (int j0 = 1; j0 < ny - 1; j0+=BLOCKY) {
			int j1= j0+BLOCKY<ny-1?j0+BLOCKY:ny-1;
			for (int k = k0; k < k1; k++) {
			    for (int j = j0; j < j1; j++){
				double *Un0 = &Un[k][j][0];
				double *up = &U[k-1][j][0];
				double *down = &U[k+1][j][0];
				double *east = &U[k][j][0]-1;
				double *west = &U[k][j][1];
				double *north = &U[k][j+1][0];
				double *south = &U[k][j-1][0];
				double *bcentral = &b[k][j][0];
				for (int i = 1; i < nx-1; i++) Un0[i]= c1 * (up[i] + down[i] + east[i] + west[i] + north[i] + south[i] - c2*bcentral[i]);
			    }
			}
		    }
		}
		double ***temp = NULL;
		temp = U;
		U = Un;
		Un = temp;
	    }
            iteration++;
	}
	bool Dependency(){
	    return true; //just for now, I'll update this function later
	}
	void PostCompletion(){
            if(iteration==nIters){
               free(U);
               free(Un);
               free(b);
            }
	}
};

int Jacobi::tx, Jacobi::ty, Jacobi::tz;
int Jacobi::nx, Jacobi::ny, Jacobi::nz;
int Jacobi::nIters;
int main(int argc,char *argv[])
{
    int argCount = 0;
    int verbose=0;
    int rank, nProcs;
    /* Argument list
       -tx: number of tasks in X dimension
       -ty: number of tasks in Y dimension
       -tz: number of tasks in Z dimension
       -nx: number of cells in X dimension
       -ny: number of cells in Y dimension
       -nz: number of cells in Z dimension
       -i : number of iterations
       -v: print out task graph information
     */ 
    while(++argCount <argc) {
	if(!strcmp(argv[argCount], "-tx")) Jacobi::tx= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-ty")) Jacobi::ty= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-tz")) Jacobi::tz= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-nx")) Jacobi::nx = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-ny")) Jacobi::ny = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-nz")) Jacobi::nz = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-i"))  Jacobi::nIters = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-v"))  verbose = true;
    }
    RTS rts;
    rts.RTS_Init(&rank, &nProcs);
    string graphName= "3DJacobi";
    if(verbose && rank==0){
	cout<< "Creating a 3DJacobi Graph with ( "<< Jacobi::tx << ", " << Jacobi::ty <<", " << Jacobi::tz << ") tasks" << endl;
	cout<< "Running the graph with "<< rts.RTS_ProcCount() << " processes" <<endl;
    }
    ArrayGraph<Jacobi, 3> *JacobiGraph= new ArrayGraph<Jacobi, 3>(graphName, PointVect<3>(Jacobi::tx, Jacobi::ty, Jacobi::tz), rank, nProcs);
    rts.RTS_Run(JacobiGraph, false);
    if(verbose && rank==0) printf("Graph execution completed\n");
    rts.RTS_Finalize();
};
