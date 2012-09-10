!! The FABIO_MODULE manages doing input and output of fabs
!! and multifabs in a 'legacy' fashion, that is so that they
!! can be procesed by amrvis?d and the amrderive code.
module fabio_module

  use bl_types
  use fab_module
  use multifab_module
  use ml_boxarray_module
  use ml_multifab_module

  implicit none

  integer, parameter :: FABIO_DOUBLE = 1
  integer, parameter :: FABIO_SINGLE = 2
  

  interface
     subroutine fabio_close(fd)
       integer, intent(out) :: fd
     end subroutine fabio_close

     subroutine fabio_read_d(fd, offset, d, count)
       use bl_types
       integer, intent(in) :: offset, fd, count
       real(kind=dp_t), intent(out) :: d(count)
     end subroutine fabio_read_d
     subroutine fabio_read_s(fd, offset, s, count)
       use bl_types
       integer, intent(in) :: offset, fd, count
       real(kind=sp_t), intent(out) :: s(count)
     end subroutine fabio_read_s

     subroutine fabio_read_comp_d(fd, offset, skip, d, count)
       use bl_types
       integer, intent(in) :: offset, fd, count, skip
       real(kind=dp_t), intent(out) :: d(count)
     end subroutine fabio_read_comp_d
     subroutine fabio_read_comp_s(fd, offset, skip, s, count)
       use bl_types
       integer, intent(in) :: offset, fd, count, skip
       real(kind=sp_t), intent(out) :: s(count)
     end subroutine fabio_read_comp_s

     subroutine fabio_write_raw_d(fd, offset, d, count, dm, lo, hi, nd, nc)
       use bl_types
       integer, intent(in) :: fd, count, dm, lo(dm), hi(dm), nd(dm), nc
       real(kind=dp_t), intent(in) :: d(count)
       integer, intent(out) :: offset
     end subroutine fabio_write_raw_d
     subroutine fabio_write_raw_s(fd, offset, s, count, dm, lo, hi, nd, nc)
       use bl_types
       integer, intent(in) :: fd, count, dm, lo(dm), hi(dm), nd(dm), nc
       real(kind=sp_t), intent(in) :: s(count)
       integer, intent(out) :: offset
     end subroutine fabio_write_raw_s
     !
     ! These are used by the particle code.
     !
     subroutine fabio_write_raw_array_i(fd, iv, count)
       integer, intent(in) :: fd, count
       integer, intent(in) :: iv(count)
     end subroutine fabio_write_raw_array_i
     subroutine fabio_write_raw_array_d(fd, rv, count)
       integer, intent(in)          :: fd, count
       double precision, intent(in) :: rv(count)
     end subroutine fabio_write_raw_array_d
     subroutine fabio_read_raw_array_i(fd, iv, count)
       integer, intent(in)    :: fd, count
       integer, intent(inout) :: iv(count)
     end subroutine fabio_read_raw_array_i
     subroutine fabio_read_raw_array_d(fd, rv, count)
       integer, intent(in)             :: fd, count
       double precision, intent(inout) :: rv(count)
     end subroutine fabio_read_raw_array_d

  end interface

  private :: fabio_write_raw_d
  private :: fabio_write_raw_s

  interface fabio_write
     module procedure fabio_fab_write_d
     module procedure fabio_multifab_write_d
  end interface

  interface fabio_ml_write
     module procedure fabio_ml_multifab_write_d
     module procedure fabio_ml_mf_write
  end interface

  integer, parameter :: FABIO_RDONLY = 0
  integer, parameter :: FABIO_WRONLY = 1
  integer, parameter :: FABIO_RDWR   = 2
  integer, parameter :: FABIO_APPEND = 3

  integer, parameter :: FABIO_MAX_VAR_NAME = 20
  integer, parameter :: FABIO_MAX_PATH_NAME = 128

contains

  subroutine fabio_mkdir(dirname, stat)
    use bl_string_module
    character(len=*), intent(in) :: dirname
    integer, intent(out), optional :: stat
    interface
       subroutine fabio_mkdir_str(ifilename, stat)
         integer, intent(in) :: ifilename(*)
         integer, intent(out) :: stat
       end subroutine fabio_mkdir_str
    end interface
    integer :: istr(128)
    integer :: lstat

    ! octal conversion 0755
    lstat = 0; if ( present(stat) ) lstat = 1
    call str2int(istr, 128, dirname)
    call fabio_mkdir_str(istr, lstat)
    if ( present(stat) ) stat = lstat

  end subroutine fabio_mkdir

  subroutine fabio_open(fd, filename, mode)
    use bl_string_module
    character(len=*), intent(in):: filename
    integer, intent(out) :: fd
    integer, intent(in), optional :: mode
    interface
       subroutine fabio_open_str(fd, ifilename, mode)
         integer, intent(out) :: fd
         integer, intent(in) :: ifilename(*)
         integer, intent(in) :: mode
       end subroutine fabio_open_str
    end interface
    integer :: istr(128)
    integer :: lmode

    lmode = FABIO_RDONLY
    if ( present(mode) ) then
       if ( mode /= 0 ) lmode = mode
    end if

    call str2int(istr, 128, filename)
    call fabio_open_str(fd, istr, lmode)

  end subroutine fabio_open

  subroutine fabio_unlink_if_empty(filename)
    use bl_string_module
    character(len=*), intent(in):: filename
    interface
       subroutine fabio_unlink_if_empty_str(ifilename)
         integer, intent(in) :: ifilename(*)
       end subroutine fabio_unlink_if_empty_str
    end interface
    integer :: istr(128)

    call str2int(istr, 128, filename)
    call fabio_unlink_if_empty_str(istr)

  end subroutine fabio_unlink_if_empty

  subroutine fabio_fab_write_d(fd, offset, fb, idx, nodal, all, prec)
    use bl_error_module
    integer, intent(in) :: fd, idx
    integer, intent(out) :: offset
    type(multifab), intent(in) :: fb
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: all
    integer, intent(in), optional :: prec
    type(box) :: bx
    logical :: lall
    real(kind=dp_t), pointer :: fbp(:,:,:,:)

    real(kind=sp_t), allocatable :: sbp(:)

    integer :: count, lo(get_dim(fb)), hi(get_dim(fb)), nd(get_dim(fb)), nc, lprec

    integer :: i,j,k,l,m

    lprec = FABIO_DOUBLE; if ( present(prec) ) lprec = prec
    if ( lprec /= FABIO_DOUBLE .and. lprec /= FABIO_SINGLE ) then
       call bl_error("FABIO_WRITE: prec is wrong ", lprec)
    end if
    lall = .false.; if ( present(all) ) lall = all
    if ( lall ) then
       bx = get_pbox(fb,idx)
    else
       bx = get_ibox(fb,idx)
    end if
    count = volume(bx)
    fbp => dataptr(fb, idx, bx)
    nc = ncomp(fb)
    lo = lwb(bx)
    hi = upb(bx)
    nd = 0
    if ( present(nodal) ) then
       where ( nodal ) nd = 1
    end if
    if ( lprec == FABIO_DOUBLE ) then
       call fabio_write_raw_d(fd, offset, fbp, count, get_dim(fb), lo, hi, nd, nc)
    else
       allocate(sbp(nc*count))
       m = 1
       do l = 1, nc
          do k = 1, size(fbp,3)
             do j = 1, size(fbp,2)
                do i = 1, size(fbp,1)
                   sbp(m) = fbp(i,j,k,l)
                   m = m + 1
                end do
             end do
          end do
       end do
       call fabio_write_raw_s(fd, offset, sbp, count, get_dim(fb), lo, hi, nd, nc)
    end if

  end subroutine fabio_fab_write_d

  subroutine fabio_multifab_write_d(mf, dirname, header, all, prec, nOutFiles, lUsingNFiles)
    use parallel
    use bl_IO_module
    use bl_error_module

    type(multifab),   intent(in)           :: mf
    character(len=*), intent(in)           :: dirname, header
    logical,          intent(in), optional :: all
    integer,          intent(in), optional :: prec
    integer,          intent(in), optional :: nOutFiles
    logical,          intent(in), optional :: lUsingNFiles

    type(layout) :: mf_la
    character(len=128) :: fname
    integer :: un
    integer :: nc, nb, i, fd, j, ng, ii
    integer, allocatable :: offset(:), loffset(:)
    type(box) :: bx
    real(kind=dp_t), allocatable :: mx(:,:), mn(:,:)
    real(kind=dp_t), allocatable :: mxl(:),  mnl(:)

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    integer, parameter :: MSG_TAG = 1010

    integer :: nSets, mySet, iSet, wakeUpPID, waitForPID, tag, iBuff(2)
    integer :: nOutFilesLoc
    logical :: lUsingNFilesLoc

    logical :: lall, nodalflags(get_dim(mf))

    lall = .false. ; if ( present(all) ) lall = all

    nodalflags = nodal_flags(mf)

    mf_la = get_layout(mf)

    nOutFilesLoc = 64        ; if (present(nOutFiles))    nOutFilesLoc    = nOutFiles
    lUsingNFilesLoc = .true. ; if (present(lUsingNFiles)) lUsingNFilesLoc = lUsingNFiles

    nOutFilesLoc = max(1, min(nOutFilesLoc, parallel_nprocs()))

    nc = multifab_ncomp(mf)
    nb = nboxes(mf%la)
    ng = 0
    if ( lall ) ng = nghost(mf)
    allocate(offset(nb),loffset(nb))
    allocate(mnl(nc), mxl(nc))
    if ( parallel_IOProcessor() ) then
       allocate(mx(nc,nb), mn(nc,nb))
       un = unit_new()
       call fabio_mkdir(dirname)
       open(unit=un, &
            file = trim(dirname) // "/" // trim(header) // "_H", &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
       write(unit=un, fmt='(i0/i0/i0/i0)') 1, 0, nc, ng
       write(unit=un, fmt='("(",i0," 0")') nb
       do i = 1, nb
          bx = get_box(mf%la, i)
          call box_print(bx, unit = un, legacy = .True., nodal = nodalflags)
       end do
       write(unit=un, fmt='(")")')
    end if
    call parallel_barrier()

    offset = -Huge(offset)
    !
    ! Each processor writes his own FABS.
    !
    if ( lUsingNFilesLoc ) then
      write(unit=fname, fmt='(a,"_D_",i5.5)') trim(header), mod(parallel_myproc(), nOutFilesLoc)

      nSets = (parallel_nprocs() + (nOutFilesLoc - 1)) / nOutFilesLoc
      mySet = parallel_myproc() / nOutFilesLoc

      do iSet = 0, nSets - 1
        if (mySet == iSet) then
          if (iSet == 0) then
            call fabio_open(fd, trim(dirname) // "/" // trim(fname), FABIO_WRONLY)
          else
            call fabio_open(fd, trim(dirname) // "/" // trim(fname), FABIO_APPEND)
          end if

          do i = 1, nfabs(mf)
             call fabio_write(fd, offset(global_index(mf,i)), mf, i, nodal = nodalflags, all = all, prec = prec)
          end do

          call fabio_close(fd)

          wakeUpPID = parallel_myproc() + nOutFilesLoc
          tag       = mod(parallel_myproc(), nOutFilesLoc)
          iBuff(1)  = tag
          iBuff(2)  = wakeUpPID
          if (wakeUpPID < parallel_nprocs()) &
              call parallel_send(iBuff, wakeUpPID, tag)
        end if

        if (mySet == (iSet + 1)) then
          waitForPID = parallel_myproc() - nOutFilesLoc
          tag        = mod(parallel_myproc(), nOutFilesLoc)
          iBuff(1)   = tag
          iBuff(2)   = waitForPID
          call parallel_recv(iBuff, waitForPID, tag)
        end if
      end do
    else
      write(unit=fname, fmt='(a,"_D_",i5.5)') trim(header), parallel_myproc()
      call fabio_open(fd, trim(dirname) // "/" // trim(fname), FABIO_WRONLY)
      do i = 1, nfabs(mf)
         call fabio_write(fd, offset(global_index(mf,i)), mf, i, nodal = nodalflags, all = all, prec = prec)
      end do
      call fabio_close(fd)
    end if

    call parallel_reduce(loffset, offset, MPI_MAX, parallel_IOProcessorNode())

    do i = 1, nb
       if ( local(mf%la,i) ) then
          ii = local_index(mf,i)
          do j = 1, nc
             if ( lall ) then
                dp => dataptr(mf, ii, get_pbox(mf,ii), j, 1)
             else
                dp => dataptr(mf, ii, get_ibox(mf,ii), j, 1)
             end if
             mnl(j) = minval(dp)
             mxl(j) = maxval(dp)
          end do
       end if
       if ( parallel_IOProcessor() ) then
          if ( remote(mf%la,i) ) then
             call parallel_recv(mn(:,i), get_proc(mf_la,i), MSG_TAG)
             call parallel_recv(mx(:,i), get_proc(mf_la,i), MSG_TAG + 1)
          else
             mx(:,i) = mxl
             mn(:,i) = mnl
          end if
       else
          if ( local(mf%la,i) ) then
            call parallel_send(mnl, parallel_IOProcessorNode(), MSG_TAG)
            call parallel_send(mxl, parallel_IOProcessorNode(), MSG_TAG + 1)
          end if
       end if
    end do

    if ( parallel_IOProcessor() ) then
       write(unit=un, fmt='(i0)') nb
       if ( lUsingNFilesLoc ) then
         do i = 1, nb
            write(unit=fname, fmt='(a,"_D_",i5.5)') trim(header), mod(get_proc(mf_la, i), nOutFilesLoc)
            write(unit=un, fmt='("FabOnDisk: ", a, " ", i0)') trim(fname), loffset(i)
         end do
       else
         do i = 1, nb
            write(unit=fname, fmt='(a,"_D_",i5.5)') trim(header), get_proc(mf_la, i)
            write(unit=un, fmt='("FabOnDisk: ", a, " ", i0)') trim(fname), loffset(i)
         end do
       end if
       write(unit=un, fmt='()')
       write(unit=un, fmt='(i0,",",i0)') nb, nc
       do i = 1, nb
          do j = 1, nc
             write(unit=un, fmt='(es27.17e3,1x,",")') mn(j,i)
          end do
       end do
       write(unit=un, fmt='()')
       write(unit=un, fmt='(i0,",",i0)') nb, nc
       do i = 1, nb
          do j = 1, nc
             write(unit=un, fmt='(es27.17e3,1x,",")') mx(j,i)
          end do
       end do
       close(unit=un)
    end if
  end subroutine fabio_multifab_write_d

  subroutine fabio_multifab_read_d(mf, dirname, header)
    use bl_stream_module
    use bl_IO_module
    use bl_error_module
    type(multifab), intent(out) :: mf
    character(len=*), intent(in) :: dirname, header
    integer :: lun
    type(bl_stream) :: strm
    character(len=256) :: str

    call build(strm)
    lun = bl_stream_the_unit(strm)

    open(unit=lun, &
         file = trim(dirname) // "/" // trim(header) // "_H" , &
         status = 'old', action = 'read')
    read(unit=lun,fmt='(a)') str
    if ( str == '&MULTIFAB' ) then
       call bl_error("PLOTFILE_BUILD: not implemented")
    else
       !! if it is not a namelist, we assume it is a VISMF_MULTIFAB
       call build_vismf_multifab
    end if
    call destroy(strm)

  contains

    !       Open Header of sub-directory
    !       : iiii: dummy, dummy, ncomponents, ng
    !       : i ; '(', nboxes dummy
    !       For each box, j
    !       [
    !         : b : bx[j]
    !       ]
    !       :  : ')'
    !       For each box, j
    !       [
    !         : ci : 'FabOnDisk: ' Filename[j], Offset[j]
    !       ]
    !       : i : nboxes, ncomponents
    !       For each box, j
    !       [
    !         : r : min[j]
    !       ]
    !       : i : nboxes, ncomponents
    !       For each box, j
    !       [
    !         : r : man[j]
    !       ]

    subroutine build_vismf_multifab()
      use bl_error_module
      integer :: j, nc
      character(len=FABIO_MAX_PATH_NAME) :: cdummy, filename
      integer :: offset
      integer :: idummy, sz, fd
      integer :: dm
      type(box), allocatable :: bxs(:)
      type(box) :: bx
      type(boxarray) :: ba
      type(layout) :: la
      integer :: nboxes, ng, itype
      real(kind=dp_t), pointer :: pp(:,:,:,:)
      logical :: nodal(MAX_SPACEDIM)

      rewind(lun)
      read(unit=lun, fmt=*) idummy, itype, nc, ng
      if ( itype /= 0 ) then
         call bl_error("BUILD_VISMF_MULTIFAB: can't read this kind of file", itype)
      end if
      call bl_stream_expect(strm, '(')
      nboxes = bl_stream_scan_int(strm)
      allocate(bxs(nboxes))
      idummy = bl_stream_scan_int(strm)
      do j = 1, nboxes
         call box_read(bx, unit = lun, nodal = nodal)
         bxs(j) = box_denodalize(bx, nodal = nodal)
      end do
      call bl_stream_expect(strm, ')')
      call build(ba, bxs)
      call build(la, ba, boxarray_bbox(ba))
      dm = get_dim(ba)
      call build(mf, la, nc = nc, ng = ng, nodal = nodal(1:dm))
      read(unit=lun, fmt=*) idummy
      do j = 1, nboxes
         read(unit=lun, fmt=*) cdummy, &
              filename, offset
         call fabio_open(fd,                         &
              trim(dirname) // "/" // &
              trim(filename))
         bx = get_pbox(mf, j)
         pp => dataptr(mf, j, bx)
         sz = volume(bxs(j))
         call fabio_read_d(fd, offset, pp(:,:,:,:), sz*nc)
         call fabio_close(fd)
      end do
      deallocate(bxs)
      call destroy(ba)
      close(unit=lun)
    end subroutine build_vismf_multifab
  end subroutine fabio_multifab_read_d

  subroutine fabio_ml_multifab_write_d(mfs, rrs, dirname, names, bounding_box, &
       prob_lo, prob_hi, time, dx, nOutFiles, lUsingNFiles, prec)
    use parallel
    use bl_IO_module
    use bl_error_module
    type(multifab), intent(in) :: mfs(:)
    integer, intent(in) :: rrs(:)
    character(len=*), intent(in) :: dirname
    character(len=FABIO_MAX_VAR_NAME), intent(in), optional :: names(:)
    type(box), intent(in), optional :: bounding_box
    real(kind=dp_t), intent(in), optional :: time
    real(kind=dp_t), intent(in), optional :: dx(:)
    real(kind=dp_t), intent(in), optional :: prob_lo(:), prob_hi(:)
    integer,          intent(in), optional :: nOutFiles
    logical,          intent(in), optional :: lUsingNFiles
    integer,          intent(in), optional :: prec

    integer :: i, j, k
    character(len=128) :: header, sd_name
    integer :: nc, un, nl, dm
    real(kind=dp_t), allocatable :: plo(:), phi(:), ldx(:), ldxlev(:)
    real(kind=dp_t), allocatable :: gridlo(:), gridhi(:)
    integer, allocatable ::  lo(:),  hi(:)
    integer :: idummy, rdummy
    type(box) :: lbbox
    real(kind=dp_t) :: ltime
    
    if ( size(mfs) < 1 ) then
       call bl_error("FABIO_ML_MULTIFAB_WRITE_D: write a zero length mlmf")
    end if
!    if ( size(mfs) /= size(rrs) + 1 ) then
!       call bl_error("FABIO_ML_MULTIFAB_WRITE_D: size of mfs /= size (rrs,dim=1)+1")
!    end if

!    nl = size(mfs)
    nl = size(rrs)+1
    nc = ncomp(mfs(1))
    if ( nc == 0 ) then
       if ( parallel_IOProcessor() ) then
          call bl_warn("FABIO_ML_MULTIFAB_WRITE_D: no components in mfs")
       end if
       return
    end if
    dm = get_dim(mfs(1))
    allocate(plo(dm),phi(dm),ldx(dm),ldxlev(dm),lo(dm),hi(dm),gridlo(dm),gridhi(dm))
    if ( present(bounding_box) ) then
       lbbox = bounding_box
    else
       lbbox = bbox(get_boxarray(mfs(1)))
    end if
    ltime = 0.0_dp_t; if ( present(time) ) ltime = time

    idummy = 0
    rdummy = 0.0_dp_t
    lo = lwb(lbbox); hi = upb(lbbox)
    ldx = 0
    if ( present(dx) ) then
       ldx = dx
    else
       ldx  = 1.0_dp_t/(maxval(hi-lo+1))
    end if

    if ( present(prob_lo) ) then
       plo = prob_lo(1:dm)
    else
       plo = lwb(lbbox)*ldx
    endif
    if ( present(prob_hi) ) then
       phi = prob_hi(1:dm)
    else
       phi = (upb(lbbox)+1)*ldx
    endif

    if ( parallel_IOProcessor() ) then
       call fabio_mkdir(dirname)
    end if
    call parallel_barrier()

    do i = 1, nl
       write(unit=sd_name, fmt='(a,"/Level_",i2.2)') trim(dirname), i-1
       call fabio_multifab_write_d(mfs(i), sd_name, "Cell", nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles, prec = prec)
    end do

    if ( parallel_IOProcessor() ) then
       header = "Header"
       un = unit_new()
       open(unit=un, &
            file = trim(dirname) // "/" // trim(header), &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
       write(unit=un, fmt='("NavierStokes-V1.1")')
       write(unit=un, fmt='(i0)') nc
       if ( present(names) ) then
          do i = 1, nc
             write(unit=un, fmt='(A)') trim(names(i))
          end do
       else
          do i = 1, nc
             write(unit=un, fmt='("Var-",i3.3)') i
          end do
       end if
       write(unit=un, fmt='(i1)') dm
       write(unit=un, fmt='(es27.17e3)') ltime
       write(unit=un, fmt='(i0)') nl - 1
       write(unit=un, fmt='(3es27.17e3)') plo
       write(unit=un, fmt='(3es27.17e3)') phi
       do i = 1, nl - 1
          write(unit=un, fmt='(i0,1x)', advance='no') rrs(i)
       end do
       write(unit=un, fmt='()')
       do i = 1, nl
          call box_print(lbbox, unit=un, legacy = .True., advance = 'no')
          write(unit=un, fmt='(" ")', advance = 'no')
          if ( i < nl ) lbbox = refine(lbbox, rrs(i))
       end do
       write(unit=un, fmt='()')
       do i = 1, nl
          write(unit=un, fmt='(i0,1x)', advance = 'no') idummy
       end do
       write(unit=un, fmt='()')
       ldxlev = ldx
       do i = 1, nl
          write(unit=un, fmt='(3es27.17e3)') ldx
          if ( i < nl ) ldx = ldx/rrs(i)
       end do
       write(unit=un, fmt='(i0)') idummy
       write(unit=un, fmt='(i0)') idummy
       ! SOME STUFF
       do i = 1, nl
          write(unit=un, fmt='(i1)', advance = 'no') i-1
          write(unit=un, fmt='(i6)', advance = 'no') nboxes(mfs(i)%la)
          write(unit=un, fmt='(i6)') rdummy
          write(unit=un, fmt='(i1)') idummy
          do j = 1, nboxes(mfs(i)%la)
             gridlo = plo + ldxlev*lwb(get_box(mfs(i)%la,j))
             gridhi = plo + ldxlev*(upb(get_box(mfs(i)%la,j))+1)
             do k = 1, dm
                write(unit=un, fmt='(2es27.17e3)') gridlo(k), gridhi(k)
             end do
          end do
          if (i < nl) ldxlev = ldxlev/rrs(i)
          write(unit=un, fmt='("Level_",i2.2,"/Cell")') i-1
       end do
       close(unit=un)
    end if
  end subroutine fabio_ml_multifab_write_d

  subroutine fabio_ml_multifab_read_d(mmf, root, unit, ng)
    use bl_stream_module
    use bl_IO_module
    use bl_error_module
    type(multifab), pointer :: mmf(:)
    character(len=*), intent(in) :: root
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: ng
    integer :: lun
    type(bl_stream) :: strm
    character(len=256) :: str
    integer :: lng
    logical :: useoldplotreader = .false.

    lng = 0; if ( present(ng) ) lng = ng

    call build(strm, unit)
    lun = bl_stream_the_unit(strm)
    open(unit=lun, &
         file = trim(root) // "/" // "Header", &
         status = 'old', action = 'read')
    read(unit=lun,fmt='(a)') str
    if ( str == '&ML_MULTIFAB' ) then
       call bl_error("PLOTFILE_BUILD: not implemented")
    else if ( str == 'NavierStokes-V1.1' .or. str == 'HyperCLaw-V1.1' ) then 
       if(useoldplotreader) then
         call build_ns_plotfile_old
       else
         call build_ns_plotfile
       endif
    else
       call bl_error('FABIO_ML_MULTIFAB_WRITE_D: Header has improper magic string', str)
    end if
    call destroy(strm)

  contains

    ! NavierStokes-V1.1 Plotfile Formats
    ! Record
    !     : c : NavierStokes-V1.1/HyperClaw-V1.1
    !     : c : Numbers of fields = n
    !    n: i : Field Names
    !     : i : Dimension = dm
    !     : r : Time
    !     : i : Number of Levels - 1 : nl
    !     : r : Physical domain lo end [1:dm]
    !     : r : Physical domain hi end [1:dm]
    !     : i : Refinement Ratios [1:nl-1]
    !     : b : Prob domains per level [1:nl]
    !     : i : unused [1:nl]
    !   nl: r : grid spacing, per level, [1:dm]
    !     : i : unused  :
    !     : i : unused
    !     For each level
    !     [
    !       : iiri : dummy, nboxes, dummy, dummy
    !       For each box, j
    !       [
    !         : r :  plo[1:dm,j], phi[1:dm, j]
    !       ]
    !       : c : level directory
    !     ]
    !     Close Header File
    !     For each level
    !     [
    !       Open Header of sub-directory
    !       : iiii: dummy, dummy, ncomponents, dummy
    !       : i ; '(', nboxes dummy
    !       For each box, j
    !       [
    !         : b : bx[j]
    !       ]
    !       :  : ')'
    !       For each box, j
    !       [
    !         : ci : 'FabOnDisk: ' Filename[j], Offset[j]
    !       ]
    !       : i : nboxes, ncomponents
    !       For each box, j
    !       [
    !         : r : min[j]
    !       ]
    !       : i : nboxes, ncomponents
    !       For each box, j
    !       [
    !         : r : man[j]
    !       ]
    !       Close subgrid file
    !     ]



    subroutine build_ns_plotfile_old()
      use bl_error_module
      integer :: i, n, k, j, nc, jj
      character(len=FABIO_MAX_PATH_NAME) :: str, str1, cdummy, filename
      integer :: offset, idummy, sz, fd, llng
      real(kind=dp_t) :: rdummy, tm
      integer :: nvars, dm, flevel
      integer, allocatable :: refrat(:,:), nboxes(:)
      real(kind=dp_t), allocatable :: dxlev(:,:)
      type(box), allocatable :: bxs(:)
      type(box) :: bx_dummy, bx
      type(boxarray) :: ba
      type(layout) :: la
      character(len=256), allocatable :: fileprefix(:)
      character(len=256), allocatable :: header(:)
      real(kind=dp_t), pointer :: pp(:,:,:,:)
      logical :: nodal(MAX_SPACEDIM)

      read(unit=lun,fmt=*) nvars
      do i = 1, nvars
         read(unit=lun,fmt='(a)') cdummy
      end do
      read(unit=lun, fmt=*) dm
      read(unit=lun, fmt=*) tm
      read(unit=lun, fmt=*) flevel
      flevel = flevel + 1

      allocate(mmf(flevel))

      allocate(nboxes(flevel))
      allocate(fileprefix(flevel))
      allocate(header(flevel))

      read(unit=lun, fmt=*) (rdummy, k=1, 2*dm)
      !! Not make this really work correctly, I need to see if these are
      !! IntVects here.  I have no examples of this.
      allocate(refrat(flevel-1,1:dm))
      read(unit=lun, fmt=*) refrat(:,1)
      refrat(:,2:dm) = spread(refrat(:,1), dim=2, ncopies=dm-1)

      do i = 1, flevel
         call box_read(bx_dummy, unit = lun, nodal = nodal(1:dm))
      end do
      read(unit=lun, fmt=*) (idummy, i=1, flevel)
      allocate(dxlev(flevel,1:dm))
      do i = 1, flevel
         read(unit=lun, fmt=*) dxlev(i,:)
      end do

      read(unit=lun, fmt=*) idummy, idummy
      do i = 1, flevel
         read(unit=lun, fmt=*) idummy, nboxes(i), rdummy, idummy
         do j = 1, nboxes(i)
            read(unit=lun, fmt=*) (rdummy, k=1, 2*dm)
         end do
         read(unit=lun, fmt='(a)') str
         str1 = str(:index(str, "/")-1)
         fileprefix(i) = str1
         str1 = trim(str(index(str, "/")+1:)) // "_H"
         header(i) = trim(str1)
      end do
      close(unit=lun)
      do i = 1, flevel
         open(unit=lun, &
              action = 'read', &
              status = 'old', file = trim(trim(root) // "/" //  &
              trim(fileprefix(i)) // "/" // &
              trim(header(i))) )
         read(unit=lun, fmt=*) idummy, idummy, nc, llng
         if ( llng > lng ) then
            call bl_error("BUILD_PLOTFILE: confused lng", lng)
         end if
         allocate(bxs(nboxes(i)))
         if ( nc /= nvars ) &
              call bl_error("BUILD_PLOTFILE: unexpected nc", nc)
         call bl_stream_expect(strm, '(')
         n = bl_stream_scan_int(strm)
         if ( n /= nboxes(i) ) &
              call bl_error("BUILD_PLOTFILE: unexpected n", n)
         idummy = bl_stream_scan_int(strm)
         do j = 1, nboxes(i)
            call box_read(bx, unit = lun, nodal = nodal(1:dm))
            bxs(j) = box_denodalize(bx, nodal = nodal(1:dm))
         end do
         call bl_stream_expect(strm, ')')
         call build(ba, bxs)
         call build(la, ba, boxarray_bbox(ba))
         call build(mmf(i), la, nc = nvars, ng = ng, nodal = nodal(1:dm))
         read(unit=lun, fmt=*) idummy
         do j = 1, nboxes(i)
            read(unit=lun, fmt=*) cdummy, &
                 filename, offset
            if (remote(mmf(i)%la,j)) cycle
            call fabio_open(fd,                         &
               trim(root) // "/" //                &
               trim(fileprefix(i)) // "/" // &
               trim(filename))
            jj = local_index(mmf(i),j)
            bx = grow(get_ibox(mmf(i), jj), lng)
            pp => dataptr(mmf(i), jj, bx)
            sz = volume(get_ibox(mmf(i),jj))
            call fabio_read_d(fd, offset, pp(:,:,:,:), sz*nvars)
            call fabio_close(fd)
         end do
         deallocate(bxs)
         call destroy(ba)
         close(unit=lun)
      end do

      deallocate(nboxes)
      deallocate(fileprefix)
      deallocate(header)
      deallocate(refrat)
      deallocate(dxlev)

    end subroutine build_ns_plotfile_old

    subroutine build_ns_plotfile()
      use bl_error_module
      integer :: i, n, k, j, nc, offset, idummy, sz, fd, llng
      character(len=FABIO_MAX_PATH_NAME) :: str, str1, cdummy, filename
      real(kind=dp_t) :: rdummy, tm
      integer :: nvars, dm, flevel, jj
      integer, allocatable :: refrat(:,:), nboxes(:)
      real(kind=dp_t), allocatable :: dxlev(:,:)
      type(box), allocatable :: bxs(:)
      type(box) :: bx_dummy, bx
      type(boxarray) :: ba
      type(boxarray), allocatable :: balevs(:)
      type(layout) :: la
      character(len=256), allocatable :: fileprefix(:)
      character(len=256), allocatable :: header(:)
      real(kind=dp_t), pointer :: pp(:,:,:,:)
      logical :: nodal(MAX_SPACEDIM)
      integer :: nSets, mySet, iSet, wakeUpPID, waitForPID, tag, nAtOnce, iBuff(2)

      nAtOnce = min(parallel_nprocs(), 64)
      nSets = (parallel_nprocs() + (nAtOnce - 1)) / nAtOnce
      mySet = parallel_myproc() / nAtOnce

      do iSet = 0, nSets - 1
        if (mySet == iSet) then
          read(unit=lun,fmt=*) nvars
          do i = 1, nvars
             read(unit=lun,fmt='(a)') cdummy
          end do
          read(unit=lun, fmt=*) dm
          read(unit=lun, fmt=*) tm
          read(unit=lun, fmt=*) flevel
          flevel = flevel + 1

          allocate(mmf(flevel))

          allocate(nboxes(flevel))
          allocate(fileprefix(flevel))
          allocate(header(flevel))

          read(unit=lun, fmt=*) (rdummy, k=1, 2*dm)
          !! Not make this really work correctly, I need to see if these are
          !! IntVects here.  I have no examples of this.
          allocate(refrat(flevel-1,1:dm))
          read(unit=lun, fmt=*) refrat(:,1)
          refrat(:,2:dm) = spread(refrat(:,1), dim=2, ncopies=dm-1)

          do i = 1, flevel
             call box_read(bx_dummy, unit = lun, nodal = nodal(1:dm))
          end do
          read(unit=lun, fmt=*) (idummy, i=1, flevel)
          allocate(dxlev(flevel,1:dm))
          do i = 1, flevel
             read(unit=lun, fmt=*) dxlev(i,:)
          end do

          read(unit=lun, fmt=*) idummy, idummy
          do i = 1, flevel
             read(unit=lun, fmt=*) idummy, nboxes(i), rdummy, idummy
             do j = 1, nboxes(i)
                read(unit=lun, fmt=*) (rdummy, k=1, 2*dm)
             end do
             read(unit=lun, fmt='(a)') str
             str1 = str(:index(str, "/")-1)
             fileprefix(i) = str1
             str1 = trim(str(index(str, "/")+1:)) // "_H"
             header(i) = trim(str1)
          end do
          close(unit=lun)
          allocate(balevs(flevel))
          do i = 1, flevel
             open(unit=lun, &
                  action = 'read', &
                  status = 'old', file = trim(trim(root) // "/" //  &
                  trim(fileprefix(i)) // "/" // &
                  trim(header(i))) )
             read(unit=lun, fmt=*) idummy, idummy, nc, llng
             if ( llng > lng ) then
                call bl_error("BUILD_PLOTFILE: confused lng", lng)
             end if
             allocate(bxs(nboxes(i)))
             if ( nc /= nvars ) &
                  call bl_error("BUILD_PLOTFILE: unexpected nc", nc)
             call bl_stream_expect(strm, '(')
             n = bl_stream_scan_int(strm)
             if ( n /= nboxes(i) ) &
                  call bl_error("BUILD_PLOTFILE: unexpected n", n)
             idummy = bl_stream_scan_int(strm)
             do j = 1, nboxes(i)
                call box_read(bx, unit = lun, nodal = nodal(1:dm))
                bxs(j) = box_denodalize(bx, nodal = nodal(1:dm))
             end do
             call bl_stream_expect(strm, ')')
             call build(ba, bxs)
             call boxarray_build_copy(balevs(i), ba)

             close(unit=lun)
             deallocate(bxs)
             call destroy(ba)
          end do

          wakeUpPID = parallel_myproc() + nAtOnce
          tag       = mod(parallel_myproc(), nAtOnce)
          iBuff(1)  = tag
          iBuff(2)  = wakeUpPID
          if (wakeUpPID < parallel_nprocs()) then
            call parallel_send(iBuff, wakeUpPID, tag)
          endif

        end if    !  mySet

        if (mySet == (iSet + 1)) then
          waitForPID = parallel_myproc() - nAtOnce
          tag        = mod(parallel_myproc(), nAtOnce)
          iBuff(1)   = tag
          iBuff(2)   = waitForPID
          call parallel_recv(iBuff, waitForPID, tag)
        end if
      enddo     !  iSet

      do i = 1, flevel
         call build(la, balevs(i), boxarray_bbox(balevs(i)))
         call build(mmf(i), la, nc = nvars, ng = ng, nodal = nodal(1:dm))
      end do


      do iSet = 0, nSets - 1
        if (mySet == iSet) then
          do i = 1, flevel
             open(unit=lun, &
                  action = 'read', &
                  status = 'old', file = trim(trim(root) // "/" //  &
                  trim(fileprefix(i)) // "/" // &
                  trim(header(i))) )
             read(unit=lun, fmt=*) idummy, idummy, nc, llng
             if ( llng > lng ) then
                call bl_error("BUILD_PLOTFILE: confused lng", lng)
             end if
             call bl_stream_expect(strm, '(')
             n = bl_stream_scan_int(strm)
             if ( n /= nboxes(i) ) &
                  call bl_error("BUILD_PLOTFILE: unexpected n", n)
             idummy = bl_stream_scan_int(strm)
             do j = 1, nboxes(i)
                call box_read(bx, unit = lun, nodal = nodal(1:dm))
             end do
             call bl_stream_expect(strm, ')')

             read(unit=lun, fmt=*) idummy
             do j = 1, nboxes(i)
                read(unit=lun, fmt=*) cdummy, &
                     filename, offset
                if (remote(mmf(i)%la,j)) cycle
                call fabio_open(fd,                         &
                   trim(root) // "/" //                &
                   trim(fileprefix(i)) // "/" // &
                   trim(filename))
                jj = local_index(mmf(i),j)
                bx = grow(get_ibox(mmf(i), jj), lng)
                pp => dataptr(mmf(i), jj, bx)
                sz = volume(get_ibox(mmf(i),jj))
                call fabio_read_d(fd, offset, pp(:,:,:,:), sz*nvars)
                call fabio_close(fd)
             end do
             close(unit=lun)
          end do

          deallocate(nboxes)
          deallocate(fileprefix)
          deallocate(header)
          deallocate(refrat)
          deallocate(dxlev)
          do i = 1, size(balevs)
             call destroy(balevs(i))
          end do
          deallocate(balevs)

          wakeUpPID = parallel_myproc() + nAtOnce
          tag       = mod(parallel_myproc(), nAtOnce)
          iBuff(1)  = tag
          iBuff(2)  = wakeUpPID
          if (wakeUpPID < parallel_nprocs()) then
            call parallel_send(iBuff, wakeUpPID, tag)
          endif

        end if    !  mySet

        if (mySet == (iSet + 1)) then
          waitForPID = parallel_myproc() - nAtOnce
          tag        = mod(parallel_myproc(), nAtOnce)
          iBuff(1)   = tag
          iBuff(2)   = waitForPID
          call parallel_recv(iBuff, waitForPID, tag)
        end if

      enddo     !  iSet

    end subroutine build_ns_plotfile

  end subroutine fabio_ml_multifab_read_d


  subroutine fabio_ml_boxarray_read(mba, root)
    use bl_stream_module
    use bl_IO_module
    use bl_error_module
    type(ml_boxarray), intent(out) :: mba
    character(len=*), intent(in) :: root
    integer :: lun
    type(bl_stream) :: strm, strm1
    character(len=256) :: str

    call build(strm)
    lun = bl_stream_the_unit(strm)
    open(unit=lun, &
         file = trim(root) // "/" // "Header", &
         status = 'old', action = 'read')
    read(unit=lun,fmt='(a)') str
    if ( str == '&ML_MULTIFAB' ) then
       call bl_error("PLOTFILE_BUILD: not implemented")
    else if ( str == 'NavierStokes-V1.1' .or. str == 'HyperCLaw-V1.1' ) then 
       call build_ns_plotfile
    else
       call bl_error('FABIO_ML_MULTIFAB_WRITE_D: Header has improper magic string', str)
    end if
    call destroy(strm)

  contains

    ! NavierStokes-V1.1 Plotfile Formats
    ! Record
    !     : c : NavierStokes-V1.1/HyperClaw-V1.1
    !     : c : Numbers of fields = n
    !    n: i : Field Names
    !     : i : Dimension = dm
    !     : r : Time
    !     : i : Number of Levels - 1 : nl
    !     : r : Physical domain lo end [1:dm]
    !     : r : Physical domain hi end [1:dm]
    !     : i : Refinement Ratios [1:nl-1]
    !     : b : Prob domains per level [1:nl]
    !     : i : unused [1:nl]
    !   nl: r : grid spacing, per level, [1:dm]
    !     : i : unused  :
    !     : i : unused
    !     For each level
    !     [
    !       : iiri : dummy, nboxes, dummy, dummy
    !       For each box, j
    !       [
    !         : r :  plo[1:dm,j], phi[1:dm, j]
    !       ]
    !       : c : level directory
    !     ]
    !     Close Header File
    !     For each level
    !     [
    !       Open Header of sub-directory
    !       : iiii: dummy, dummy, ncomponents, dummy
    !       : i ; '(', nboxes dummy
    !       For each box, j
    !       [
    !         : b : bx[j]
    !       ]
    !       :  : ')'
    !       For each box, j
    !       [
    !         : ci : 'FabOnDisk: ' Filename[j], Offset[j]
    !       ]
    !       : i : nboxes, ncomponents
    !       For each box, j
    !       [
    !         : r : min[j]
    !       ]
    !       : i : nboxes, ncomponents
    !       For each box, j
    !       [
    !         : r : man[j]
    !       ]
    !       Close subgrid file
    !     ]

    subroutine build_ns_plotfile()
      use bl_error_module
      integer :: i, n, k
      integer :: j, nc
      character(len=256) :: str, str1, cdummy
      integer :: idummy, lun1
      real(kind=dp_t) :: rdummy
      integer :: nvars, dm, flevel
      integer :: nboxes
      type(box), allocatable :: bxs(:)
      type(box) :: bdummy
      character(len=256) :: fileprefix
      character(len=256) :: header

      read(unit=lun,fmt=*) nvars
      do i = 1, nvars
         read(unit=lun,fmt='(a)') cdummy
      end do
      read(unit=lun, fmt=*) dm
      read(unit=lun, fmt=*) rdummy
      read(unit=lun, fmt=*) flevel
      flevel = flevel + 1

      call build(mba, flevel, dm)

      read(unit=lun, fmt=*) (rdummy, k=1, 2*dm)
      !! Not make this really work correctly, I need to see if these are
      !! IntVects here.  I have no examples of this.
      read(unit=lun, fmt=*) mba%rr(:,1)
      mba%rr(:,2:dm) = spread(mba%rr(:,1), dim=2, ncopies=dm-1)

      do i = 1, flevel
         call box_read(bdummy, unit = lun)
      end do
      read(unit=lun, fmt=*) (idummy, i=1, flevel)
      do i = 1, flevel
         read(unit=lun, fmt=*) (rdummy,j=1,dm)
      end do

      read(unit=lun, fmt=*) idummy, idummy
      do i = 1, flevel
         read(unit=lun, fmt=*) idummy, nboxes, rdummy, idummy
         do j = 1, nboxes
            read(unit=lun, fmt=*) (rdummy, k=1, 2*dm)
         end do
         read(unit=lun, fmt='(a)') str
         str1 = str(:index(str, "/")-1)
         fileprefix = str1
         str1 = trim(str(index(str, "/")+1:)) // "_H"
         header = trim(str1)
         call bl_stream_build(strm1)
         lun1 = bl_stream_the_unit(strm1)
         open(unit=lun1, &
              action = 'read', &
              status = 'old', file = trim(trim(root) // "/" //  &
              trim(fileprefix) // "/" // &
              trim(header)) )
         read(unit=lun1, fmt=*) idummy, idummy, nc, idummy
         allocate(bxs(nboxes))
         if ( nc /= nvars ) &
              call bl_error("BUILD_PLOTFILE: unexpected nc", nc)
         call bl_stream_expect(strm1, '(')
         n = bl_stream_scan_int(strm1)
         if ( n /= nboxes ) &
              call bl_error("BUILD_PLOTFILE: unexpected n", n)
         idummy = bl_stream_scan_int(strm1)
         do j = 1, nboxes
            call box_read(bxs(j), unit = lun1)
         end do
         call bl_stream_expect(strm1, ')')
         call build(mba%bas(i), bxs)
         deallocate(bxs)
         close(unit=lun1)
         call destroy(strm1)
      end do
      close(unit=lun)
    end subroutine build_ns_plotfile
  end subroutine fabio_ml_boxarray_read

  ! FABIO
  subroutine fabio_ml_mf_write(xx, fname, names, time, dx)
    type(ml_multifab), intent(in) :: xx
    character(len=FABIO_MAX_VAR_NAME), intent(in), optional :: names(:)
    character(len=*), intent(in) :: fname
    real(kind=dp_t), intent(in), optional :: dx(:)
    real(kind=dp_t), intent(in), optional :: time

    call fabio_ml_multifab_write_d(xx%mf, xx%mla%mba%rr(:,1), fname, &
         names = names, time = time, dx = dx)

  end subroutine fabio_ml_mf_write

end module fabio_module
