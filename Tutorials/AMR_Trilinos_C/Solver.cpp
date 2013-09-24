#include <EpetraExt_MatrixMatrix.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <ParallelDescriptor.H>

#include "Solver.H"

void Solver::SetupProblem(const Box& domain, MultiFab& rhs, MultiFab& soln) 
{

    // Now we create the map
    // Note: 1) it is important to set  numMyGridPoints = 0 outside the MFIter loop
    //       2) "bx" here is the cell-centered non-overlapping box
    int numMyGridPoints = 0;
    std::vector<int> MyGlobalElements;

    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.validbox();
       const int* bx_lo = bx.loVect();
       const int* bx_hi = bx.hiVect();

       for(int i = bx_lo[0]; i<= bx_hi[0]; i++)
          for(int j = bx_lo[1]; j<= bx_hi[1]; j++)
          {
#if (BL_SPACEDIM == 3)
             for(int k = bx_lo[2]; k<= bx_hi[2]; k++)
#endif
                if (isInside(D_DECL(i,j,k)))
                {
                    MyGlobalElements.push_back(getIdx(D_DECL(i,j,k)));
                    numMyGridPoints++;
                }
          }
    }

    // Define the Map
    Map = new Epetra_Map(-1, numMyGridPoints, &MyGlobalElements[0], 0, Comm_m);

    // Allocate the RHS and LHS with the new Epetra Map
    RHS = rcp(new Epetra_Vector(*Map));
    LHS = rcp(new Epetra_Vector(*Map));
    
    // Note: it is important to set idx = 0 outside the MFIter loop
    int idx = 0;

    // Copy the values from rhs into RHS->Values()
    // Copy the values from soln into LHS->Values()
    // Note: "bx" here is the cell-centered non-overlapping box
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.validbox();
       const int* bx_lo = bx.loVect();
       const int* bx_hi = bx.hiVect();

       FArrayBox& fab_rhs =  rhs[mfi];
       FArrayBox& fab_lhs = soln[mfi];

       for(int i = bx_lo[0]; i <= bx_hi[0]; i++) 
          for(int j = bx_lo[1]; j <= bx_hi[1]; j++) 
          {
#if (BL_SPACEDIM == 3)
             for(int k = bx_lo[2]; k<= bx_hi[2]; k++)
#endif
                if (isInside(D_DECL(i,j,k)))
                {
                    IntVect cells(D_DECL(i,j,k));
                    RHS->Values()[idx  ] = fab_rhs(cells,0);
                    LHS->Values()[idx++] = fab_lhs(cells,0);
                }
          }
    }

    if(verbose_m)
        this->printLoadBalanceStats();

    A = rcp(new Epetra_CrsMatrix(Copy, *Map, 7));
    ComputeStencil();
}

void Solver::extrapolateLHS() 
{
// Aitken-Neville
// Pi0 (x) := yi , i = 0 : n
// Pik (x) := (x − xi ) Pi+1,k−1(x) − (x − xi+k ) Pi,k−1(x) /(xi+k − xi )
// k = 1, . . . , n, i = 0, . . . , n − k.

    std::deque< Epetra_Vector >::iterator it = OldLHS.begin();

    if(nLHS_m == 0)
        LHS->PutScalar(1.0);
    else if(OldLHS.size() == 1)
        *LHS = *it;
    else if(OldLHS.size() == 2){
        LHS->Update (2.0, *it++, -1.0, *it, 0.0);
    }
    else if(OldLHS.size() > 0){
        int n = OldLHS.size();
        for(int i=0; i<n; ++i){
            *(*P)(i) = *it++;
        }
        for(int k = 1; k < n; ++k){// x==0 & xi==i+1
            for(int i = 0; i < n-k; ++i){
                (*P)(i)->Update(-(i+1)/(float)k, *(*P)(i+1), (i+k+1)/(float)k);//TODO test
            }
        }
        *LHS = *(*P)(0);
    }
    else
        std::cout << "Invalid number of old LHS: " + OldLHS.size() << std::endl;
}

void Solver::CopySolution(Box& domain, MultiFab& soln) 
{
    int idx = 0;

    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.validbox();
       const int* bx_lo = bx.loVect();
       const int* bx_hi = bx.hiVect();

       FArrayBox& fab_lhs = soln[mfi];

       for(int i = bx_lo[0]; i <= bx_hi[0]; i++) 
          for(int j = bx_lo[1]; j <= bx_hi[1]; j++) 
#if (BL_SPACEDIM == 3)
             for(int k = bx_lo[2]; k<= bx_hi[2]; k++)
#endif
                if (isInside(D_DECL(i,j,k)))
                {
                    IntVect cells(D_DECL(i,j,k));
                    fab_lhs(cells,0) = LHS->Values()[idx++];
                }
    }
}

int Solver::getNumIters()
{
    return solver->getNumIters();
}

void Solver::Compute()
{
    //LHS->Random();
    //LHS->PutScalar(1.0);
    extrapolateLHS();

    // create the preconditioner object and compute hierarchy
    // true -> create the multilevel hirarchy
    // ML allows the user to cheaply recompute the preconditioner. You can
    // simply uncomment the following line:
    //
    // MLPrec->ReComputePreconditioner();
    //
    // It is supposed that the linear system matrix has different values, but
    // **exactly** the same structure and layout. The code re-built the
    // hierarchy and re-setup the smoothers and the coarse solver using
    // already available information on the hierarchy. A particular
    // care is required to use ReComputePreconditioner() with nonzero
    // threshold.

    if(MLPrec == Teuchos::null) // first repetition we need to create a new preconditioner
        MLPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
    else if(isReusingHierarchy_m)
        MLPrec->ReComputePreconditioner();
    else if(isReusingPreconditioner_m) {
        // do nothing since we are reusing old preconditioner
    } else { // create a new preconditioner in every repetition
        delete MLPrec.get();//MLPrec now RCP => delete??? TODO
        MLPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, MLList_m));
    }

    // Setup problem
    problem.setOperator(A);
    problem.setLHS(LHS);
    problem.setRHS(RHS);
    prec = rcp(new Belos::EpetraPrecOp(MLPrec));
    problem.setLeftPrec(prec);
    solver->setParameters(rcp(&belosList, false));
    solver->setProblem(rcp(&problem,false));
    if(!problem.isProblemSet()){
        if (problem.setProblem() == false) {
            std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        }
    }


    // Solve problem
    // Timer
    MPI_Barrier(MPI_COMM_WORLD);
    Real dRunTime1 = ParallelDescriptor::second();

    solver->solve();

    // Timer
    MPI_Barrier(MPI_COMM_WORLD);
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    Real dRunTime2 = ParallelDescriptor::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {   
        std::cout << "Run time = " << dRunTime2 << std::endl;
    }   

    // Store new LHS in OldLHS
    OldLHS.push_front(*(LHS.get()));
    if(OldLHS.size() > nLHS_m) OldLHS.pop_back();
    std::cout<<"#OldLHS: "<<OldLHS.size()<<std::endl;
}

void Solver::ComputeStencil() {

    int NumMyElements = Map->NumMyElements();
    int* MyGlobalElements = Map->MyGlobalElements();

    std::vector<double> Values(6);
    std::vector<int> Indices(6);

    printf("NumMyElements = %d\n",NumMyElements);
    for (int i = 0; i < NumMyElements; i++)
    {
        int NumEntries = 0;

#if (BL_SPACEDIM == 2)
        double WV, EV, SV, NV, CV;
        // use boundary object to get stencil values
        getBoundaryStencil(MyGlobalElements[i], WV, EV, SV, NV, CV);

        double W, E, S, N;
        getNeighbours(MyGlobalElements[i], W, E, S, N);
#elif (BL_SPACEDIM == 3)
        double WV, EV, SV, NV, FV, BV, CV;
        // use boundary object to get stencil values
        getBoundaryStencil(MyGlobalElements[i], WV, EV, SV, NV, FV, BV, CV);

        double W, E, S, N, F, B;
        getNeighbours(MyGlobalElements[i], W, E, S, N, F, B);
#endif

        if(E != -1) {
            Indices[NumEntries] = E;
            Values[NumEntries++] = EV;
        }
        if(W != -1) {
            Indices[NumEntries] = W;
            Values[NumEntries++] = WV;
        }
        if(S != -1) {
            Indices[NumEntries] = S;
            Values[NumEntries++] = SV;
        }
        if(N != -1) {
            Indices[NumEntries] = N;
            Values[NumEntries++] = NV;
        }
#if (BL_SPACEDIM == 3)
        if(F != -1) {
            Indices[NumEntries] = F;
            Values[NumEntries++] = FV;
        }
        if(B != -1) {
            Indices[NumEntries] = B;
            Values[NumEntries++] = BV;
        }
#endif

        // put the off-diagonal entries
        A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);

        // put in the diagonal entry
        A->InsertGlobalValues(MyGlobalElements[i], 1, &CV, MyGlobalElements + i);
    }

    A->FillComplete();
    A->OptimizeStorage();
}

void Solver::printLoadBalanceStats() {

    //compute some load balance statistics
    size_t myNumPart = Map->NumMyElements();
    size_t NumPart = Map->NumGlobalElements() * 1.0/Comm_m.NumProc();
    double imbalance = 1.0;
    if(myNumPart >= NumPart)
        imbalance += (myNumPart-NumPart)/NumPart;
    else
        imbalance += (NumPart-myNumPart)/NumPart;

    double max=0.0, min=0.0, avg=0.0;
    int minn=0, maxn=0;
    MPI_Reduce(&imbalance, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&imbalance, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&imbalance, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myNumPart, &minn, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&myNumPart, &maxn, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    avg /= Comm_m.NumProc();
    if(Comm_m.MyPID() == 0) cout << "LBAL min = " << min << ", max = " << max << ", avg = " << avg << endl;
    if(Comm_m.MyPID() == 0) cout << "min nr gridpoints = " << minn << ", max nr gridpoints = " << maxn << endl;

}
