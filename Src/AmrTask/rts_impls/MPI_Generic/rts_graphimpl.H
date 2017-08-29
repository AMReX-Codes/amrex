//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 07-21-2017

#include <iostream>
#include <queue>
using namespace std;
#include <cassert>

namespace amrex{

    template<typename T>
    void ReductionSum_impl(T *local, T *global, int length, int root){
        for(int i=0; i<length; i++) global[i]=local[i];
    }

}//end namespace

