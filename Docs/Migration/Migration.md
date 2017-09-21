
# Migration From BoxLib

To help C++ BoxLib users migrate to AMReX, we have provided a set of
tools in `Tools/Migration`.  We recommend you go through the following
steps in the order listed.  You probably should perform the migration
in a branch.  For each step, it is recommended that you start with a
fresh clone of your latest migration branch.  It is assumed that the
environment variable `AMREX_HOME` is set to the AMReX directory.

**The scripts used in this process typically perform search and
  replace in the directory tree rooted at the current directory.  So
  do not run them outside your application code.  After you run a
  script, you should do `git diff` to check the results.**

## Step 0

Make sure your code works with the [latest version of BoxLib on
github](https://github.com/BoxLib-Codes/BoxLib). 

## Step 1

AMReX `migration/1-amrex_home` branch should be used in this step.  In
this step, we replace `BOXLIB_HOME` with `AMREX_HOME` in the GNU Make
system by running the following command in your application directory.

    $AMREX_HOME/Tools/Migration/step-1-amrex_home/amrex_home.sh

## Step 2

AMReX `migration/2-parray` branch should be used in this step.  In
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
  can call `amrex::FillNull(myArray)` to fill it with null pointers.
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

### Step 3

AMReX `migration/3-amrex-prefix` branch should be used in this step.
We have added `AMReX_` to source file names.  Thus we must update the
header file names in the source codes.  For example, `#include
<Box.H>` needs to be changed to `#include <AMReX_Box.H>`.  A script,
`Tools/Migration/step-3-amrex-prefix/amrexprefix.sh`, can be used to
do the search and replace.

### Step 4

AMReX `migration/4-dirname` branch should be used in this step.  In
AMRex, we have renamed several directories:

* Tools/C_mk --> Tools/GNUMake
* C_AMRLib --> Amr
* C_AmrCoreLib --> AmrCore
* C_BaseLib --> Base
* C_BoundaryLib --> Boundary
* C_ParticleLib --> Particle

A script `Tools/Migration/step-4-dirname/dirname.sh` can be used to
perform the search and replace.

### Step 5

AMReX `migration/5-namespace` branch should be used in this step.
BoxLib has a `BoxLib` namspace, but most of BoxLib classes are not in
the namespace.  In AMReX, all classes and functions have been moved
into the `amrex` namespace.  In this step, you can use
`Tools/Migration/step-5-amrex-namespace/amrex-namespace.sh` to replace
`BoxLib::` with `amrex::` for those already in `BoxLib` namespace.
However, the rest of work is expected to be performed manually,
because C++ is too a complicated language for shell scripting.  

For most `.cpp` files, you can put a `using namespace amrex;` line
after the last `include` line, or `using amrex::MultiFab` etc., or you
can add `amrex::` to wherever needed.  Note that having both `using
namespace amrex` and `using namespace std` in one file may cause
conflicts because some names like `min` and `max` exist in both
namespace. 

For header files, it is considered bad practice to have `using
namespace amrex` because it pollutes the namespace.  So you need to
manually add `amrex::` in front of AMReX names likes `MultiFab` and
`BoxArray`.

### Step 6

AMReX `migration/6-distributionmap` branch should be used in this step. 

In BoxLib, there is a `DistributionMapping` cache implemented with
`std::map` with the number of `Box`es as the key.  Utilizing the
cache, `MultiFab`s can be built with shared `DistributionMapping`
without explicitly requiring a `DistributionMapping`.  In AMReX, the
`DistributionMapping` cache is removed for more flexibility.  Because
of the change, `DistributionMapping` is now a required argument to the
non-default constructors of `MultiFab`, `FluxRegister`, `BndryData`,
`AmrLevel`, `LevelBld`, and some other objects.  Many classes
including `MultiFab`, `AmrLevel`, `AmrCore` have functions returning
`DistributionMapping`.  One may also explicitly construct
`DistributionMapping` from `BoxArray`.  For example,

    DistributionMapping dm {a_BoxArray};

It should be emphasized that the result of the code above does **not**
solely depend on `BoxArray`.  Thus, different `DistributionMapping`s
may be built when this is called multiple times with the same
`BoxArray`.  If two `MultiFab`s need to share the same
`DistributionMapping`, only one `DistributionMapping` should be built.

For extensibility, the `MultiFab` constructor now also takes an `MFInfo`
object.  To construct a `MultiFab` without allocating memory for the
data,

    MFInfo info;
    info.SetAlloc(false);   // The default MFInfo is true.

`VisMF::Read` function used to only take an empty `MultiFab` built
with the default constructor.  The function reads `BoxArray` from
disk, builds the `MultiFab` with the `BoxArray` and a
`DistributionMapping` possibly from the `DistributionMapping` cache if
there is one cached for the same number of boxes, and then read the
data from the disk into the `MultiFab`.  Because the cache is removed,
the `VisMF::Read` function has been modified to also take a pre-built
`MultiFab`.  In that case, it is an error if the `MultiFab`'s
`BoxArray` does not match the `BoxArray` on the disk.  If an empty
`MultiFab` is passed into `VisMF::Read`, a **new**
`DistributionMapping` will be built, and this maybe not be desired
behavior.  `VisMF::Read` function is rarely called directly by an
application code.  But if it is, care must be taken to ensure that
`MultiFab`s read from the disk are distributed correctly among
processes.

### Step 7

AMReX `migration/7-bindc` branch should be used in this step.  In
BoxLib, the `Amr` class calls a Fortran subroutine `probinit` to perform
problem-specific initialization on Fortran side.  This is a subroutine
provided by application codes.  In AMReX, this subroutine is renamed
to `amrex_probint`, and it is expected to have `bind(c)` making the
interface of calling it from C++ clean.  A script
`Tools/Migration/step-7-bindc/bindc.sh` can be used to for the
migration.

### Step 8

AMReX `migration/8-deboxlib` branch should be used in this step.  In
the branch, `AMReX_winstd.H` is removed because Windows is not
supported.  `AMReX_BoxLib.H` is renamed `AMReX.H`.  There are also
some changes to `buildInfo` related to `AMReX_` prefix for file names
and `amrex` namespace.  A script
`Tools/Migration/step-8-deboxlib/deboxlib.sh` can be used to for the
migration.

In BoxLib, there are some runtime parameters with the `boxlib.` prefix
(e.g., `boxlib.fpe_trap_invalid` and `boxlib.verbose`).  The prefix is
now changed to `amrex.`.

### Step 9

AMReX `migration/9-pointers-n-sentinel` branch should be used in this
step.

`Amr` and `AmrLevel` used to return `MultiFab *` in the `derive`
function; they now return `std::unique_ptr<MultiFab>` for safety.
Application codes using `Amr` and `AmrLevel` classes need to update by
changing `MultiFab* mf = derive(...)` to `auto mf = derive(...)` and
removing the `delete` call.

AMReX's own smart pointer classes have been removed.  For application
codes that do not directly use them (which is usually the case), no
changes are needed.

We have removed the sentinel from the `DistributionMapping` class.
`DistributionMapping` has an `Array<int>` storing MPI ranks for all
grids in a `BoxArray`.  It used to be that the length of the `Array`
is the number of boxes plus one.  The last element was used to store
the MPI rank of the process as a sentinel.  But that information is no
longer needed.  Thus we have remove the sentinel, and the length of
the array equals the number of boxes.  If an application code builds
its own array of MPI ranks and passes it to `DistributionMapping`,
update is needed.  One might also search for `ProcessorMap` function
that returns the array to see if the array is being used and needs
update.  Most codes probably do not need any changes.


### Step 10

AMReX `migration/10-amrex-real` branch should be used in this step.

AMReX can be compiled with `PRECISION = DOUBLE` (default) or `FALSE`.
It used to be that `Real` was either `double` or `float`.  This has
now been put into the `amrex` namespace and becomes `amrex::Real`.
Fortran code can do `use amrex_fort_module, only : amrex_real` and
define variables as `real(amrex_real) :: x`.  C code can do `#include
<AMReX_REAL.H>` and use `amrex_real`.

