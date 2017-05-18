module amrex_fabio_c_module

  implicit none

  interface

     subroutine fabio_close(fd) bind(c, name='fabio_close')
       integer, intent(out) :: fd
     end subroutine fabio_close

     subroutine fabio_read_skip_d(fd, offset, skip, d, count) bind(c, name='fabio_read_skip_d')
       use iso_c_binding
       integer(kind=c_int), intent(in) :: fd
       integer(kind=c_long), intent(in) :: offset, skip, count
       real(kind=c_double), intent(out) :: d(count)
     end subroutine fabio_read_skip_d

     subroutine fabio_read_skip_s(fd, offset, skip, s, count) bind(c, name='fabio_read_skip_s')
       use iso_c_binding
       integer(kind=c_int), intent(in) :: fd
       integer(kind=c_long), intent(in) :: offset, skip, count
       real(kind=c_float), intent(out) :: s(count)
     end subroutine fabio_read_skip_s

     subroutine fabio_read_d(fd, offset, d, count) bind(c, name='fabio_read_d')
       use iso_c_binding
       integer(kind=c_int), intent(in) :: fd
       integer(kind=c_long), intent(in) :: offset, count
       real(kind=c_double), intent(out) :: d(count)
     end subroutine fabio_read_d

     subroutine fabio_read_s(fd, offset, s, count) bind(c, name='fabio_read_s')
       use iso_c_binding
       integer(kind=c_int), intent(in) :: fd
       integer(kind=c_long), intent(in) :: offset, count
       real(kind=c_float), intent(out) :: s(count)
     end subroutine fabio_read_s

     subroutine fabio_write_raw_d(fd, offset, d, count, dm, lo, hi, nd, nc) bind(c, name='fabio_write_raw_d')
       use iso_c_binding
       integer(kind=c_int), intent(in) :: fd, dm, lo(dm), hi(dm), nd(dm), nc
       integer(kind=c_long), intent(in) :: count
       real(kind=c_double), intent(in) :: d(count)
       integer(kind=c_long), intent(out) :: offset
     end subroutine fabio_write_raw_d

     subroutine fabio_write_raw_s(fd, offset, s, count, dm, lo, hi, nd, nc) bind(c, name='fabio_write_raw_s')
       use iso_c_binding
       integer(kind=c_int), intent(in) :: fd, dm, lo(dm), hi(dm), nd(dm), nc
       integer(kind=c_long), intent(in) :: count
       real(kind=c_float), intent(in) :: s(count)
       integer(kind=c_long), intent(out) :: offset
     end subroutine fabio_write_raw_s

     subroutine fabio_open_str(fd, ifilename, mode) bind(c, name='fabio_open_str')
       integer, intent(out) :: fd
       integer, intent(in) :: ifilename(*)
       integer, intent(in) :: mode
     end subroutine fabio_open_str

     subroutine fabio_unlink_if_empty_str(ifilename) bind(c, name='fabio_unlink_if_empty_str')
       integer, intent(in) :: ifilename(*)
     end subroutine fabio_unlink_if_empty_str

     subroutine fabio_mkdir_str(ifilename, stat) bind(c, name='fabio_mkdir_str')
       integer, intent(in) :: ifilename(*)
       integer, intent(out) :: stat
     end subroutine fabio_mkdir_str

     !
     ! these are used by the particle code.
     !
     subroutine fabio_write_raw_array_i(fd, iv, count) bind(c, name='fabio_write_raw_array_i')
       use iso_c_binding
       integer(kind=c_int), intent(in) :: fd, count
       integer(kind=c_int), intent(in) :: iv(count)
     end subroutine fabio_write_raw_array_i

     subroutine fabio_write_raw_array_d(fd, rv, count) bind(c, name='fabio_write_raw_array_d')
       use iso_c_binding
       integer(kind=c_int), intent(in) :: fd, count
       real(kind=c_double), intent(in) :: rv(count)
     end subroutine fabio_write_raw_array_d

     subroutine fabio_read_raw_array_i(fd, iv, count) bind(c, name='fabio_read_raw_array_i')
       use iso_c_binding
       integer(kind=c_int), intent(in)    :: fd, count
       integer(kind=c_int), intent(inout) :: iv(count)
     end subroutine fabio_read_raw_array_i

     subroutine fabio_read_raw_array_d(fd, rv, count) bind(c, name='fabio_read_raw_array_d')
       use iso_c_binding
       integer(kind=c_int), intent(in)    :: fd, count
       real(kind=c_double), intent(inout) :: rv(count)
     end subroutine fabio_read_raw_array_d

     subroutine val_is_inf(v, res) bind(c, name='val_is_inf')
       use bl_types
       real(dp_t), intent(in)  :: v
       integer,    intent(out) :: res
     end subroutine val_is_inf

     subroutine val_is_nan(v, res) bind(c, name='val_is_nan')
       use bl_types
       real(dp_t), intent(in)  :: v
       integer,    intent(out) :: res
     end subroutine val_is_nan

     subroutine fab_contains_nan(dptr, count, res) bind(c, name='fab_contains_nan')
       use bl_types
       integer,    intent(in)  :: count
       real(dp_t), intent(in)  :: dptr(count)
       integer,    intent(out) :: res
     end subroutine fab_contains_nan

     subroutine fab_contains_inf(dptr, count, res) bind(c, name='fab_contains_inf')
       use bl_types
       integer,    intent(in)  :: count
       real(dp_t), intent(in)  :: dptr(count)
       integer,    intent(out) :: res
     end subroutine fab_contains_inf

  end interface

end module amrex_fabio_c_module
