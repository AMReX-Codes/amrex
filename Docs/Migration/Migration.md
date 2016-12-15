
# Migration From BoxLib

To help C++ BoxLib users migrate to AMReX, we have provided a set of
tools in `Tools/Migration`.  We recommend you go through the following
steps in the order listed.  You probably should perform the migration
in a branch.  For each step, it is recommended that you start with a
fresh clone of your latest migration branch.  It is assume that the
environment variable `AMREX_HOME` is set to the AMReX directory. 

**The scripts used in this process typically perform search and
  replace in the directory tree rooted at the current directory.  So
  do not run them outside your application code.**

## Step 0

Make sure your code work with the [latest version of BoxLib on
github](https://github.com/BoxLib-Codes/BoxLib). 

## Step 1

AMReX `migration/1-amrex_home` branch should be used in this step.  In
this step, we replace `BOXLIB_HOME` with `AMREX_HOME` in the GNU Make
system by running the following command in your application directory.

    $AMREX_HOME/Tools/Migration/step-1-amrex_home/amrex_home.sh

## Step 2

AMRex `migration/2-parray` branch should be used in this step.  In
this step, we need to remove `PArray`, which is a BoxLib class that
has been removed from AMReX.  This step has to be done manually.

* Instead of `PArray`, one should use `Array<std::unique_ptr<T> >` for
  owning pointers and `Array<T*>` for non-owning pointers.  Note that
  `Array` is an AMReX class derived from `std::vector`.

* The `PArray` class does pointer dereferencing automatically, whereas
  now one must use either `*` or `->` like working with a raw pointer.

* Another difference is how to change a pointer stored.  For `PArray`,
  one often uses the `set` function.  The `Array` class does not have
  that function.  Instead one need to use the subscript operator `[]`
  to access the element.

* The `clear()` function in `PArray` class does not change the size of
  the `PArray` container.  But the `clear` function in `Array` reduces
  the size to zero.  If you do not want to resize the container, you
  can call `AMReX::FillNull(myArray)` to fill it with null pointers.
  `PArray` also has a `clear(int i)` function.  For that, you can do
  `a[i].reset()` for `unique_ptr` and `a[i] = nullptr` for `T*`.

* BoxLib codes often use `PArray<T>&` as function arguments.  In
  AMReX, depending on the situation, one should use
  `Array<std::unique_ptr<T> >&` for functions that need to modify the
  pointers (e.g., allocate memory and store the pointers in the
  function parameter), and `const Array<T*>&` for cases where the
  function does not change the pointers even though they may modify
  the data pointed by the pointers.  AMReX provides a number of
  functions that can convert from `unique_ptr<T>` to `T*`.   For
  example, 
  `Array<T*> GetArrOfPtrs (const Array<std::unique_ptr<T> >& a)`.
  These functions are in `Src/C_BaseLib/Array.H` as of writing,
  and they will be moved to `Src/Base/AMReX_Array.H`.

