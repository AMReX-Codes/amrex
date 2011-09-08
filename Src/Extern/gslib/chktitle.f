      subroutine chktitle(str,len)
c-----------------------------------------------------------------------
c
c                     Check for a Valid Title
c                     ***********************
c
c This subroutine takes the character string "str" of length "len" and
c blanks out all characters after the first back slash
c
c
c
c-----------------------------------------------------------------------
      parameter (MAXLEN=132)
      character str(MAXLEN)*1,strl*132
c
c Remove leading blanks:
c
      do i=1,len-1
            if(str(i).ne.' ') then
                  if(i.eq.1) go to 1
                  do j=1,len-i+1
                        k = j + i - 1
                        str(j) = str(k)
                  end do
                  do j=len,len-i+2,-1
                        str(j) = ' '
                  end do
                  go to 1
            end if
      end do
 1    continue
c
c Find first back slash and blank out the remaining characters:
c
      do i=1,len-1
            if(str(i).eq.'\\') then
                  do j=i,len
                        str(j) = ' '
                  end do
                  go to 2
            end if
      end do
 2    continue
c
c Although inefficient, copy the string:
c
      do i=1,len
            strl(i:i) = str(i)
      end do
c
c Look for a five character pattern with -Titl...
c
      do i=1,len-5
            if(strl(i:i+4).eq.'-Titl'.or.
     +         strl(i:i+4).eq.'-titl'.or.
     +         strl(i:i+4).eq.'-TITL'.or.
     +         strl(i:i+4).eq.'-X la'.or.
     +         strl(i:i+4).eq.'-Y la') then
                  do j=i,len
                        str(j) = ' '
                  end do
                  go to 3
            end if
      end do
 3    continue
c
c Return with modified character string:
c
      return
      end
