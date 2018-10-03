#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MultiFab.H>

class MyTest
{
public:

    MyTest ();

    void runTest ();

private:

    void readParameters ();
    void initData ();

    int n_tasks;
    amrex::MultiFab data_split, data_single, data_all;
    std::string outdir;
    bool fj_verbose;
};

#endif
