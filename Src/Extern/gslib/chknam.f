      subroutine chknam(str,len)
c-----------------------------------------------------------------------
c
c                   Check for a Valid File Name
c                   ***************************
c
c This subroutine takes the character string "str" of length "len" and
c removes all leading blanks and blanks out all characters after the
c first blank found in the string (leading blanks are removed first).
c
c
c
c-----------------------------------------------------------------------
      parameter (MAXLEN=512)
      character str(MAXLEN)*1
c
c Find first two blanks and blank out remaining characters:
c
      do i=1,len-1
            if(str(i)  .eq.' '.and.
     +         str(i+1).eq.' ') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 2
            end if
      end do
 2    continue
c
c Look for "-fi" for file
c
      do i=1,len-2
            if(str(i)  .eq.'-'.and.
     +         str(i+1).eq.'f'.and.
     +         str(i+2).eq.'i') then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 3
            end if
      end do
 3    continue
c
c Look for "\fi" for file
c
      do i=1,len-2
            if((str(i)  .eq. '\\').and.
     +         (str(i+1).eq. 'f').and.
     +         (str(i+2).eq. 'i')) then
                  do j=i+1,len
                        str(j) = ' '
                  end do
                  go to 4
            end if
      end do
 4    continue
c
c Return with modified file name:
c
      return
      end
