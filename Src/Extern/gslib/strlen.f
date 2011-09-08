      subroutine strlen(str,MAXLEN,lostr)
c-----------------------------------------------------------------------
c
c      Determine the length of the string minus trailing blanks
c
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      character str*512
      lostr = MAXLEN
      do i=1,MAXLEN
            j = MAXLEN - i + 1
            if(str(j:j).ne.' ') return
            lostr = lostr - 1
      end do
      return
      end
