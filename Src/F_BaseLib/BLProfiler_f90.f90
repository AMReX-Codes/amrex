! do nothing; we need these so that F90 BoxLib code can compile.

subroutine bl_proffortfuncstart(str)
  character*(*) str
end subroutine bl_proffortfuncstart

subroutine bl_proffortfuncstop(str)
  character*(*) str
end subroutine bl_proffortfuncstop

subroutine bl_proffortfuncstart_int(i)
  integer i
end subroutine bl_proffortfuncstart_int

subroutine bl_proffortfuncstop_int(i)
  integer i
end subroutine bl_proffortfuncstop_int
