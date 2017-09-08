#ifndef COLLECTIVE_IMPL
#define COLLECTIVE_IMPL

//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 07-21-2017

#include <iostream>
#include <queue>
using namespace std;
#include <mpi.h>
#include <cassert>
using std::is_same;

namespace amrex{

    template<typename T>
    void ReductionSum_impl(T *local, T *global, int length, int root){
        MPI_Datatype datatype;
        if(is_same<T, double>::value) datatype= MPI_DOUBLE;
        MPI_Reduce(local, global, length, datatype, MPI_SUM, root, MPI_COMM_WORLD);
    }

}//end namespace

#endif
