#ifndef AMREX_ASYNCOUT_H_
#define AMREX_ASYNCOUT_H_

#include <functional>

namespace amrex {
namespace AsyncOut {

struct WriteInfo {
    int ifile;
    int ispot;
    int nspots;
};

void Initialize ();
void Finalize ();

bool UseAsyncOut ();

WriteInfo GetWriteInfo (int rank);

void Submit (std::function<void()>&& a_f);
void Submit (std::function<void()> const& a_f);

void Finish (); // If you want to wait for jobs submitted to finish

//
// These functions are used inside user's job funciton.
//
void Wait ();   // Wait for my turn to write file.  This is not for waiting for job to finish.
void Notify (); // Notify next MPI process in the same file.

}}

#endif
