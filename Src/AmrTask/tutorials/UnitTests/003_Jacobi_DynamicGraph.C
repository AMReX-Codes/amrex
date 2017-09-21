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
int tx=1, ty=1, tz=1;
int Nx=256, Ny=256, Nz=256;
int nIters=10;

struct domain{
    Array3D U, Un, b;
    double *bufEast=NULL, *bufWest=NULL, *bufNorth=NULL, *bufSouth=NULL;
    double c1;
    double c2;
    int nx, ny, nz;
};

    class DynamicTaskAssociate{
        public:
        TaskName TaskAssociate(TaskName name){
	     assert(name.Dim()==4);
	     TaskName origin=name;
	     origin[3]=0;
	     return origin;
        }
    };


class Jacobi :public Task{
    private:
	domain *dom;
	Array3D U, Un, b;
	double *bufEast, *bufWest, *bufNorth, *bufSouth;
	double c1;
	double c2;
	int nx, ny, nz;

    public:
	Jacobi(domain *domIn){
	    dom = domIn;
	}
	void Job(){
	    TaskName taskID= MyName();
	    int rankx= taskID[0];
	    int ranky= taskID[1];
	    int rankz= taskID[2];
	    int iteration= taskID[3];
	    bufEast=dom->bufEast; bufWest=dom->bufWest; bufNorth=dom->bufNorth; bufSouth=dom->bufSouth;
	    c1=dom->c1;
	    c2=dom->c2;
	    nx=dom->nx; ny=dom->ny; nz=dom->nz;
	    if(iteration%2==1){ U=dom->U; Un=dom->Un;}
	    else {U=dom->Un; Un=dom->U;}
	    b= dom->b;
	    if(iteration>1){ //iterations 2 to nIters+1
		if(rankz >0) {
		    Pull(TaskName(rankx, ranky, rankz-1, iteration-1), (char*)&U[0][0][0], nx*ny*sizeof(double));
		}
		if(rankz<tz-1) {
		    Pull(TaskName(rankx, ranky, rankz+1, iteration-1), (char*)&U[nz-1][0][0], nx*ny*sizeof(double));
		}
		if(ranky >0) {
		    Pull(TaskName(rankx, ranky-1, rankz, iteration-1), (char*)bufNorth, nx*nz*sizeof(double));
		    double* idx;
		    double* ptr= (double*)bufNorth;
		    for(int z=0; z < nz; z++) {
			idx = &U[z][0][0];
			memcpy(idx, ptr, sizeof(double)*nx);
			ptr += nx;
		    }
		}
		if(ranky<ty-1) {
		    Pull(TaskName(rankx, ranky+1, rankz, iteration-1), (char*)bufSouth, nx*nz*sizeof(double));
		    double* idx;
		    double* ptr= (double*)bufSouth;
		    for(int z=0; z < nz; z++) {
			idx = &U[z][ny-1][0];
			memcpy(idx, ptr, sizeof(double)*nx);
			ptr += nx;
		    }
		}
		if(rankx >0) {
		    Pull(TaskName(rankx-1, ranky, rankz, iteration-1), (char*)bufWest, ny*nz*sizeof(double));
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
		    Pull(TaskName(rankx+1, ranky, rankz, iteration-1), (char*)bufEast, ny*nz*sizeof(double));
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


	    if(iteration<= nIters) //iterations 1 to nIters
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

	    if(iteration<=nIters){ //iterations 1 to nIters
		if(rankz >0) {
		    Push(TaskName(rankx, ranky, rankz-1, iteration+1), (char*)&Un[1][0][0], nx*ny*sizeof(double));
		}
		if(rankz<tz-1) {
		    Push(TaskName(rankx, ranky, rankz+1, iteration+1), (char*)&Un[nz-2][0][0], nx*ny*sizeof(double));
		}
		if(ranky >0) {
		    double* idx;
		    double* ptr= (double*)bufNorth;
		    for(int z=0; z < nz; z++) {
			idx = &Un[z][1][0];
			memcpy(ptr, idx, sizeof(double)*nx);
			ptr += nx;
		    }
		    Push(TaskName(rankx, ranky-1, rankz, iteration+1), (char*)bufNorth, nx*nz*sizeof(double));
		}
		if(ranky<ty-1) {
		    double* idx;
		    double* ptr= (double*)bufSouth;
		    for(int z=0; z < nz; z++) {
			idx = &Un[z][ny-2][0];
			memcpy(ptr, idx, sizeof(double)*nx);
			ptr += nx;
		    }
		    Push(TaskName(rankx, ranky+1, rankz, iteration+1), (char*)bufSouth, nx*nz*sizeof(double));
		}
		if(rankx >0) {
		    double* idx;
		    double* ptr= (double*)bufWest;
		    for(int z=0; z < nz; z++) {
			for(int y=0; y < ny; y++) {
			    idx = &Un[z][y][1];
			    *ptr= *idx;
			    ptr ++;
			}
		    }
		    Push(TaskName(rankx-1, ranky, rankz, iteration+1), (char*)bufWest, ny*nz*sizeof(double));
		}
		if(rankx<tx-1) {
		    double* idx;
		    double* ptr= (double*)bufEast;
		    for(int z=0; z < nz; z++) {
			for(int y=0; y < ny; y++) {
			    idx = &Un[z][y][nx-2];
			    *ptr= *idx;
			    ptr ++;
			}
		    }
		    Push(TaskName(rankx+1, ranky, rankz, iteration+1), (char*)bufEast, ny*nz*sizeof(double));
		}
	    }
	}
	bool Dependency(){
	    TaskName name= MyName();
	    int iteration= name[3];
	    if(iteration==1) return true;
	    else{
		bool satisfied=true;
		if(name[0]>0){ 
		    satisfied= Depend_on(TaskName(name[0]-1, name[1], name[2], iteration-1)); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[0]<tx-1){
		    satisfied= Depend_on(TaskName(name[0]+1, name[1], name[2], iteration-1)); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[1]>0){
		    satisfied= Depend_on(TaskName(name[0], name[1]-1, name[2], iteration-1)); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[1]<ty-1){
		    satisfied= Depend_on(TaskName(name[0], name[1]+1, name[2], iteration-1)); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[2]>0){
		    satisfied= Depend_on(TaskName(name[0], name[1], name[2]-1, iteration-1)); 
		    if(!satisfied) return false; //early cascade
		}
		if(name[2]<tz-1){
		    satisfied= Depend_on(TaskName(name[0], name[1], name[2]+1, iteration-1)); 
		    if(!satisfied) return false; //early cascade
		}
		return true;
	    }
	}
	void PostCompletion(){
	    if(MyName()[3]<=nIters){
		Jacobi *jac= new Jacobi(dom);
		TaskName newName= MyName();
		newName[3]++;
		jac->SetName(newName);
		RegisterTask(jac);
	    }else{
		double res= residual(Un, nx, ny, nz);
		LocalAtomicAdd(&global_err, res);
		free(dom->U);
		free(dom->Un);
		free(dom->b);
		free(dom->bufEast);
		free(dom->bufWest);
		free(dom->bufNorth);
		free(dom->bufSouth);
		free(dom);
	    }
	    SelfDestroy();
	}
};
class JacobiInit :public Task{
    private:
	domain *dom;
	int rankx, ranky, rankz;
    public:
	JacobiInit(){
	}
	void Job(){
	    dom= new domain();
	    TaskName taskID= MyName();
	    rankx= taskID[0];
	    ranky= taskID[1];
	    rankz= taskID[2];
	    dom->c1=1.0/6.0;
	    dom->c2=1.0;
	    dom->nx= Nx/tx+2;
	    dom->ny= Ny/ty+2;
	    dom->nz= Nz/tz+2;
	    dom->U = Array3DMalloc(dom->nx, dom->ny, dom->nz);
	    dom->Un = Array3DMalloc(dom->nx, dom->ny, dom->nz);
	    dom->b = Array3DMalloc(dom->nx, dom->ny, dom->nz);
	    for (int k=0; k<dom->nz; k++)
		for (int j=0; j<dom->ny; j++)
		    for(int i=0; i<dom->nx; i++){
			if(((rankz==0)&&(k==0)) ||((rankz==tz-1)&&(k==dom->nz-1)) ||((ranky==0)&&(j==0)) ||\
				((ranky==ty-1)&&(j==dom->ny-1))||((rankx==0)&&(i==0)) ||((rankx==tx-1)&&(i==dom->nx-1))){
			    dom->U[k][j][i] = 0;
			    dom->Un[k][j][i] = 0;
			}else{
			    dom->U[k][j][i] = 1;
			    dom->Un[k][j][i] = 1;
			}
			dom->b[k][j][i] = 0;
		    }
	    dom->bufEast= (double*)malloc(dom->ny*dom->nz*sizeof(double));
	    dom->bufWest= (double*)malloc(dom->ny*dom->nz*sizeof(double));
	    dom->bufNorth= (double*)malloc(dom->nx*dom->nz*sizeof(double));
	    dom->bufSouth= (double*)malloc(dom->nx*dom->nz*sizeof(double));
	}
	bool Dependency(){return true;}
	void PostCompletion(){
	    Jacobi *jac= new Jacobi(dom);
	    TaskName newName= MyName();
	    newName[3]++;
	    jac->SetName(newName);
	    RegisterTask(jac);
	    SelfDestroy();
	}
};

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
	if(!strcmp(argv[argCount], "-tx")) tx= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-ty")) ty= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-tz")) tz= atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-Nx")) Nx = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-Ny")) Ny = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-Nz")) Nz = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-i"))  nIters = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-v"))  verbose = true;
    }
    global_err=0.;
    RTS rts;
    rts.Init();
    rank= rts.MyProc();
    nProcs= rts.ProcCount();
    string graphName= "3DJacobi";
    if(verbose && rank==0){
	cout<< "Creating a 3DJacobi Graph containing ( "<< tx << ", " << ty <<", " << tz << ") tasks" << "for iteration 0"<< endl;
	cout<< "Running the graph with "<< nProcs << " processes" <<endl;
	cout<< "The graph evolves over time"<<endl; 
    }
    double time= -rts.Time();
    rts.Barrier();
    ArrayGraph<JacobiInit, 4, DynamicTaskAssociate> *JacobiGraph= new ArrayGraph<JacobiInit, 4, DynamicTaskAssociate>(graphName, PointVect<4>(tx, ty, tz, 1), rank, nProcs);
    rts.Iterate(JacobiGraph);
    double res= global_err;
    double finalErr;
    rts.ReductionSum(&res, &finalErr, 1, 0); //reduce to process 0
    rts.Barrier();
    time +=rts.Time();
    if(rank==0) {
	cout<<"Residual: " << sqrt(finalErr/((double)(Nx+1)*(double)(Ny+1)*(double)(Nz+1))) <<endl;
	cout<<"Graph execution takes "<< time << " seconds"<<endl;
	double gflops = nIters*(double)Nx*(double)Ny*(double)Nz*8/(1.0e9);
	cout<<"GFLOP/S " << gflops/time <<endl;
    }
    delete JacobiGraph;
    rts.Finalize();
};
