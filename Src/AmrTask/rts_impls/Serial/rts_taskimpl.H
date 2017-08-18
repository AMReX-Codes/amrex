//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 07-21-2017

#include <iostream>
#include <queue>
using namespace std;
#include <cassert>

namespace amrex{

    template<typename T> void LocalAtomicAdd_impl(T *addr, T val){
        *addr+= val;
    }
    template<typename T> void GlobalAtomicAdd_impl(T *addr, T val){
        *addr+= val;
    }

}//end namespace

