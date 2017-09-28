//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 08-14-2017
#include <iostream>
#include <string.h>
#include <math.h>

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

double  residual(double*** U, const int nx, const int ny, const int nz){

    double c = (1.0/6.0);
    double err = 0;

    for (int k=1; k<nz-1; k++)
	for (int j=1; j<ny-1; j++)
	    for (int i=1; i<nx-1; i++){
		double du =  (U[k-1][j][i] + U[k+1][j][i] + U[k][j-1][i] + U[k][j+1][i] + U[k][j][i-1] + U[k][j][i+1] - 6 * U[k][j][i]);
		err = err +  du * du;
	    }
    return err;
}

using namespace amrex;

double  global_err;

class Jacobi :public Task{
    private:
	int rankx, ranky, rankz;
	Array3D U, Un, b;
	double *bufEast=NULL, *bufWest=NULL, *bufNorth=NULL, *bufSouth=NULL, *bufUp=NULL, *bufDown=NULL;
	int iteration;
	double c1;
	double c2;
	int nx, ny, nz;
    public:
	static int nIters;
	static int tx, ty, tz;
	static int Nx, Ny, Nz;
	Jacobi(){
	    iteration=-1;
	    c1=1.0/6.0;
	    c2=1.0;
	}

	void initializeData(){
	    TaskName taskID= MyName();
	    rankx= taskID[0];
	    ranky= taskID[1];
	    rankz= taskID[2];
	    nx= Nx/tx+2;
	    ny= Ny/ty+2;
	    nz= Nz/tz+2;
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
	    bufEast= (double*)malloc(ny*nz*sizeof(double));
	    bufWest= (double*)malloc(ny*nz*sizeof(double));
	    bufNorth= (double*)malloc(nx*nz*sizeof(double));
	    bufSouth= (double*)malloc(nx*nz*sizeof(double));
	}

	void Job(){
	    if(iteration==-1){
		initializeData();
	    }else if(iteration<=nIters){ //compute stencil from iter 0 to nIters-1, when iteration=nIters we just update ghost cell the from the previous iteration

		if(iteration>0){// iterations 1 to nIters
		    if(rankz >0) {
			Pull(TaskName(rankx, ranky, rankz-1), (char*)&U[0][0][0], nx*ny*sizeof(double));
		    }
		    if(rankz<tz-1) {
			Pull(TaskName(rankx, ranky, rankz+1), (char*)&U[nz-1][0][0], nx*ny*sizeof(double));
		    }
		    if(ranky >0) {
			Pull(TaskName(rankx, ranky-1, rankz), (char*)bufNorth, nx*nz*sizeof(double));
			double* idx;
			double* ptr= (double*)bufNorth;
			for(int z=0; z < nz; z++) {
			    idx = &U[z][0][0];
			    memcpy(idx, ptr, sizeof(double)*nx);
			    ptr += nx;
			}
		    }
		    if(ranky<ty-1) {
			Pull(TaskName(rankx, ranky+1, rankz), (char*)bufSouth, nx*nz*sizeof(double));
			double* idx;
			double* ptr= (double*)bufSouth;
			for(int z=0; z < nz; z++) {
			    idx = &U[z][ny-1][0];
			    memcpy(idx, ptr, sizeof(double)*nx);
			    ptr += nx;
			}
		    }
		    if(rankx >0) {
			Pull(TaskName(rankx-1, ranky, rankz), (char*)bufWest, ny*nz*sizeof(double));
			double* idx;
			double* ptr= (double*)bufWest;
			for(int z=0; z < nz; z++) {
			    for(int y=0; y < ny; y++) {
				idx = &U[z][y][0];
				*idx= *ptr;
				ptr++;
			    }
			}
		    }
		    if(rankx<tx-1) {
			Pull(TaskName(rankx+1, ranky, rankz), (char*)bufEast, ny*nz*sizeof(double));
			double* idx;
			double* ptr= (double*)bufEast;
			for(int z=0; z < nz; z++) {
			    for(int y=0; y < ny; y++) {
				idx = &U[z][y][nx-1];
				*idx= *ptr;
				ptr++;
			    }
			}
		    }
		}


		if(iteration<nIters){ //iterations 0 to nIters-1
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
				    for (int i = 1; i < nx-1; i++){
					Un0[i]= c1 * (up[i] + down[i] + east[i] + west[i] + north[i] + south[i] - c2*bcentral[i]);
				    } 
				}
			    }
			}
		    }
		    double ***temp = NULL;
		    temp = U;
		    U = Un;
		    Un = temp;
		}

		if(iteration<nIters){ //iterations 0 to nIters-1
		    if(rankz >0) {
			Push(TaskName(rankx, ranky, rankz-1), (char*)&U[1][0][0], nx*ny*sizeof(double));
		    }
		    if(rankz<tz-1) {
			Push(TaskName(rankx, ranky, rankz+1), (char*)&U[nz-2][0][0], nx*ny*sizeof(double));
		    }
		    if(ranky >0) {
			double* idx;
			double* ptr= (double*)bufNorth;
			for(int z=0; z < nz; z++) {
			    idx = &U[z][1][0];
			    memcpy(ptr, idx, sizeof(double)*nx);
			    ptr += nx;
			}
			Push(TaskName(rankx, ranky-1, rankz), (char*)bufNorth, nx*nz*sizeof(double));
		    }
		    if(ranky<ty-1) {
			double* idx;
			double* ptr= (double*)bufSouth;
			for(int z=0; z < nz; z++) {
			    idx = &U[z][ny-2][0];
			    memcpy(ptr, idx, sizeof(double)*nx);
			    ptr += nx;
			}
			Push(TaskName(rankx, ranky+1, rankz), (char*)bufSouth, nx*nz*sizeof(double));
		    }
		    if(rankx >0) {
			double* idx;
			double* ptr= (double*)bufWest;
			for(int z=0; z < nz; z++) {
			    for(int y=0; y < ny; y++) {
				idx = &U[z][y][1];
				*ptr= *idx;
				ptr ++;
			    }
			}
			Push(TaskName(rankx-1, ranky, rankz), (char*)bufWest, ny*nz*sizeof(double));
		    }
		    if(rankx<tx-1) {
			double* idx;
			double* ptr= (double*)bufEast;
			for(int z=0; z < nz; z++) {
			    for(int y=0; y < ny; y++) {
				idx = &U[z][y][nx-2];
				*ptr= *idx;
				ptr ++;
			    }
			}
			Push(TaskName(rankx+1, ranky, rankz), (char*)bufEast, ny*nz*sizeof(double));
		    }
		}
	    }
	    iteration++;
	}
	bool Dependency(){
	    if(iteration<=0) return true;
	    else{
		bool satisfied=true;
		TaskName name= MyName();
		if(name[0]>0){ 
		    satisfied= Depend_on(TaskName(name[0]-1, name[1], name[2])); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[0]<tx-1){
		    satisfied= Depend_on(TaskName(name[0]+1, name[1], name[2])); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[1]>0){
		    satisfied= Depend_on(TaskName(name[0], name[1]-1, name[2])); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[1]<ty-1){
		    satisfied= Depend_on(TaskName(name[0], name[1]+1, name[2])); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[2]>0){
		    satisfied= Depend_on(TaskName(name[0], name[1], name[2]-1)); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[2]<tz-1){
		    satisfied= Depend_on(TaskName(name[0], name[1], name[2]+1)); 
		    if(!satisfied) return false; //early cascade
		}
		return true;
	    }
	}
	void PostCompletion(){
	    if(iteration<=nIters){
		KeepTaskAlive();
	    }else{
		double res= residual(Un, nx, ny, nz);
		LocalAtomicAdd(&global_err, res);
		free(U);
		free(Un);
		free(b);
		free(bufEast);
		free(bufWest);
		free(bufNorth);
		free(bufSouth);
		SelfDestroy();
                bool satisfied=true;
                TaskName name= MyName();
	    }
	}
};

int Jacobi::tx=1, Jacobi::ty=1, Jacobi::tz=1;
int Jacobi::Nx=256, Jacobi::Ny=256, Jacobi::Nz=256;
int Jacobi::nIters=10;
int main(int argc,char *argv[])
{
    int argCount = 0;
    int verbose=0;
    int rank, nProcs;
    /* Argument list
       -tx: number of tasks in X dimension
       -ty: number of tasks in Y dimension
       -tz: number of tasks in Z dimension
       -Nx: number of cells in X dimension
       -Ny: number of cells in Y dimension
       -Nz: number of cells in Z dimension
       -i : number of iterations
       -v: print out task graph information
     */ 
    while(++argCount <argc) {
	if(!strcmp(argv[argCount], "-tx")) Jacobi::tx= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-ty")) Jacobi::ty= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-tz")) Jacobi::tz= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-Nx")) Jacobi::Nx = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-Ny")) Jacobi::Ny = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-Nz")) Jacobi::Nz = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-i"))  Jacobi::nIters = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-v"))  verbose = true;
    }
    global_err=0.;
    RTS rts;
    rts.Init();
    rank= rts.MyProc();
    nProcs= rts.ProcCount();
    string graphName= "3DJacobi";
    if(verbose && rank==0){
	cout<< "Creating a 3DJacobi Graph with ( "<< Jacobi::tx << ", " << Jacobi::ty <<", " << Jacobi::tz << ") tasks" << endl;
	cout<< "Running the graph with "<< nProcs << " processes" <<endl;
    }
    double time= -rts.Time();
    rts.Barrier();
    ArrayGraph<Jacobi, 3> *JacobiGraph= new ArrayGraph<Jacobi, 3>(graphName, PointVect<3>(Jacobi::tx, Jacobi::ty, Jacobi::tz), rank, nProcs);
    rts.Barrier();
    rts.Iterate(JacobiGraph);
    double res= global_err;
    double finalErr;
    rts.ReductionSum(&res, &finalErr, 1, 0); //reduce to process 0
    rts.Barrier();
    time +=rts.Time();
    if(rank==0) {
	cout<<"Residual: " << sqrt(finalErr/((double)(Jacobi::Nx+1)*(double)(Jacobi::Ny+1)*(double)(Jacobi::Nz+1))) <<endl;
	cout<<"Graph execution takes "<< time << " seconds"<<endl;
	double gflops = Jacobi::nIters*(double)Jacobi::Nx*(double)Jacobi::Ny*(double)Jacobi::Nz*8/(1.0e9);
	cout<<"GFLOP/S " << gflops/time <<endl;
    }
    delete JacobiGraph;
    rts.Finalize();
};
