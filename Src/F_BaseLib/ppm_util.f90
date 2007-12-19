module ppm_util_module

  interface
     subroutine store_pgm_str(ifilename, width, height, iimage)
       integer, intent(in) :: ifilename(*)
       integer, intent(in) :: width, height
       integer, intent(in) :: iimage(*)
     end subroutine store_pgm_str
     subroutine load_palette_str(ifilename, r, g, b, a)
       integer, intent(in) :: ifilename(*)
       integer, intent(out) :: r(*), g(*), b(*), a(*)
     end subroutine load_palette_str
     subroutine store_ppm_str(ifilename, width, height, iimage, r, g, b)
       integer, intent(in) :: ifilename(*), r(*), g(*), b(*)
       integer, intent(in) :: width, height
       integer, intent(in) :: iimage(*)
     end subroutine store_ppm_str
  end interface

  private :: store_pgm_str, load_palette_str, store_ppm_str

contains

  subroutine load_palette(fname, r, g, b, a)
    use bl_string_module
    integer, intent(out), dimension(:) :: r, g, b, a
    character(len=*), intent(in) :: fname
    integer :: istr(128)
    call str2int(istr, 128, fname)
    call load_palette_str(istr, r, g, b, a)
  end subroutine load_palette

  subroutine store_pgm(fname, image)
    use bl_string_module
    integer, intent(in) :: image(:,:)
    character(len=*), intent(in) :: fname
    integer :: width, height
    integer, allocatable :: iimage(:)
    integer :: istr(128)
    integer :: i, j, n

    call str2int(istr, 128, fname)
    width = size(image,dim=1); height = size(image,dim=2)
    allocate(iimage(width*height))

    n = 1
    do j = height, 1, -1
       do i = 1, width
          iimage(n) = image(i,j)
          n = n + 1
       end do
    end do
    call store_pgm_str(istr, width, height, iimage)

  end subroutine store_pgm

  subroutine store_ppm(fname, image, r, g, b)
    use bl_string_module
    integer, intent(in) :: image(:,:)
    integer, intent(in) :: r(:), g(:), b(:)
    character(len=*), intent(in) :: fname
    integer :: width, height
    integer, allocatable :: iimage(:)
    integer :: istr(128)
    integer :: i, j, n

    call str2int(istr, 128, fname)
    width = size(image,dim=1); height = size(image,dim=2)
    allocate(iimage(width*height))

    n = 1
    do j = height, 1, -1
       do i = 1, width
          iimage(n) = image(i,j)
          n = n + 1
       end do
    end do
    call store_ppm_str(istr, width, height, iimage, r, g, b)

  end subroutine store_ppm

end module ppm_util_module
