#include <map> 
#include <cmath> 
#include <iostream> 
#include <assert.h> 
#include "Solver.H" 

using namespace std; 

void Solver::InitializeDomainGeometry()
{

    IdxMap.clear();
    CoordMap.clear();

    // Build a index and coordinate map
    int idx = 0;
    for(int i=0; i < nx; i++)  
       for(int j=0; j < ny; j++)  
       {
#if (BL_SPACEDIM == 3)
           for(int k=0; k < nz; k++)  
#endif
           if(isInside(D_DECL(i,j,k))) 
           {
               IdxMap[toCoordIdx(D_DECL(i,j,k))] = idx++;
           }
       }

    idx = 0;
    for(int i=0; i < nx; i++)  
       for(int j=0; j < ny; j++)  
       {
#if (BL_SPACEDIM == 3)
          for (int k=0; k < nz; k++)  
#endif
          if(isInside(D_DECL(i,j,k))) 
               CoordMap[idx++] = toCoordIdx(D_DECL(i,j,k));
       }
}

#if (BL_SPACEDIM == 2)
void Solver::getBoundaryStencil(int i, int j, 
                                       double& W, double& E, double& S, double& N, 
                                       double& C)
{
    ConstantInterpolation(i,j,W,E,S,N,C); 

    // stencil center value has to be positive!
    assert(C > 0);
}

void Solver::getBoundaryStencil(int idx, double& W, double& E, double& S, double& N, 
                                       double& C)
{
    int i = -100000 ,j = -100000;
    getCoord(idx,i,j);
    getBoundaryStencil(i,j,W,E,S,N,C);
}

void Solver::getNeighbours(int idx, double& W, double& E, double& S, double& N)
{
    int i = 0, j = 0;
    getCoord(idx,i,j);
    getNeighbours(i,j,W,E,S,N);
}

#elif (BL_SPACEDIM == 3)
void Solver::getBoundaryStencil(int i, int j, int k, 
                                       double& W, double& E, double& S, double& N, 
                                       double& F, double& B, double& C) 
{
    ConstantInterpolation(i,j,k,W,E,S,N,F,B,C); 

    // stencil center value has to be positive!
    assert(C > 0);
}

void Solver::getBoundaryStencil(int idx, double& W, double& E, double& S, double& N, 
                                       double& F, double& B, double& C) 
{
    int i = 0, j = 0, k = 0;
    getCoord(idx,i,j,k);
    getBoundaryStencil(i,j,k,W,E,S,N,F,B,C);
}

void Solver::getNeighbours(int idx, double& W, double& E, double& S, double& N,
                                  double& F, double& B) 
{
    int i = 0, j = 0, k = 0;
    getCoord(idx,i,j,k);
    getNeighbours(i,j,k,W,E,S,N,F,B);
}
#endif

void Solver::getNeighbours(int i, int j, 
#if (BL_SPACEDIM == 3)
                                  int k, 
#endif
                                  double& W, double& E, double& S, double& N  
#if (BL_SPACEDIM == 3)
                                  ,double& F, double& B
#endif
                                  )
{
    if(i > 0)
        W = getIdx(D_DECL(i-1,j,k));
    else
        W = -1;
    if(i < nx-1)
        E = getIdx(D_DECL(i+1,j,k));
    else
        E = -1;

    if(j < ny-1)
        N = getIdx(D_DECL(i,j+1,k));
    else
        N = -1;
    if(j > 0)
        S = getIdx(D_DECL(i,j-1,k));
    else
        S = -1;

#if (BL_SPACEDIM == 3)
    if(k > 0)
        F = getIdx(i,j,k-1);
    else
        F = -1;
    if(k < nz-1)
        B = getIdx(i,j,k+1);
    else
        B = -1;
#endif
}

void Solver::ConstantInterpolation(int x, int y, 
#if (BL_SPACEDIM == 3)
                                          int z, 
#endif
                                          double& W, double& E, double& S, double& N, 
#if (BL_SPACEDIM == 3)
                                          double& F, double& B, 
#endif
                                          double& C)
{
    W = -1./(hr[0]*hr[0]);
    E = -1./(hr[0]*hr[0]);
    N = -1./(hr[1]*hr[1]);
    S = -1./(hr[1]*hr[1]);

    // we are a right boundary point
    if(!isInside(D_DECL(x+1,y,z))) {
        E = 0.0;
    } 

    // we are a left boundary point
    if(!isInside(D_DECL(x-1,y,z))) {
        W = 0.0;
    }

    // we are a upper boundary point
    if(!isInside(D_DECL(x,y+1,z))) {
        N = 0.0;
    }

    // we are a lower boundary point
    if(!isInside(D_DECL(x,y-1,z))) {
        S = 0.0;
    } 
    C = 2./(hr[0]*hr[0]) + 2./(hr[1]*hr[1]);

#if (BL_SPACEDIM == 3)
    F = -1./(hr[2]*hr[2]);
    B = -1./(hr[2]*hr[2]);
    // we are a front boundary point
    if(!isInside(x,y,z+1)) {
        F = 0.0;
    } 

    // we are a back boundary point
    if(!isInside(x,y,z-1)) {
        B = 0.0;
    } 
    C = C + 2./(hr[2]*hr[2]);
#endif
}
