c     This function convert a string to the corresponding real number
c     str2real and str2int do not check the correctness of the input

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function str2real(str,str_len)
      implicit none
      character*(*) str
      integer str_len
      integer point_pos
      real*8 decimal
      integer i,num_start,num_end,isminus 

      num_end = len_trim(str)
      if (num_end .gt. str_len) num_end = str_len
      num_start = 1
      do while(str(num_start:num_start) .eq. ' ')
        num_start = num_start + 1
      enddo 
      isminus = 0
      if (str(num_start:num_start) .eq. '-') then
        isminus = 1
        num_start = num_start + 1
      endif 
      point_pos = index(str(num_start:num_end),'.')
      if (point_pos .eq. 0) then
        point_pos = num_end + 1
      else
        point_pos = point_pos + num_start - 1
      endif
      decimal = 0.
      str2real = 0.
      do i = num_start, point_pos - 1
        str2real = 10 * str2real + ichar(str(i:i)) - 48
      enddo
      if ((point_pos > 0) .and. (point_pos .lt. num_end)) then
        do i = point_pos + 1, num_end 
          decimal =  10 * decimal + ichar(str(i:i)) - 48   
        enddo
        decimal = decimal / (10. ** (num_end - point_pos))
      endif
      str2real = str2real + decimal
      if (isminus .eq. 1) str2real = 0. - str2real
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function str2int(str,str_len)
      implicit none 
      character*(*) str
      integer str_len
      integer i,num_start,num_end,isminus

      num_end = len_trim(str)
      if (num_end .gt. str_len) num_end = str_len
      num_start = 1
      do while(str(num_start:num_start) .eq. ' ')
        num_start = num_start + 1
      enddo 
      isminus = 0
      if (str(num_start:num_start) .eq. '-') then
        isminus = 1
        num_start = num_start + 1
      endif 
      str2int = 0
      do i = num_start, num_end 
        str2int = 10 * str2int + ichar(str(i:i)) - 48
      enddo
      if (isminus .eq. 1) str2int = 0 - str2int
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function isequal(a,b)
      implicit none
      real*8 eps
      parameter (eps=1e-5)
      real*8 a,b
      real*8 c

      c = a - b
      isequal = 0
      if (abs(a-b) .lt. eps) isequal = 1
c      print *,c,isequal 
      return
      end
