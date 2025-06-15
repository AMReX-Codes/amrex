# GDB xmethod support

GDB doesn't have access to inlined methods, so it allows providing replacement
implementations for these methods in Python using the Xmethod API.
`gdb_amrex_xmethods.py` contains implementations for the non-trivial
subscription operators for the various `amrex::Array*` classes.

## Usage

Run `source <path to amrex>/Tools/GDB/gdb_amrex_xmethods.py` during an
interactive session, or add that line to `~/.gdbinit` to be run in every
session.

Note: built-in support for user-defined function call operators (`operator()`)
was only added in GDB 16 (included in GNU Binutils 2.44). If you're using an
older version, you must explicitly call the `operator()` function, like this:

```
(gdb) p array4.operator()(1, 2, 3, 4)
$1 = 1847849.7974222905
```
