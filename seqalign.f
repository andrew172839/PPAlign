cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program seqalign
      implicit none 
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      character*2000 seq1, seq2,aseq1,aseq2
      character*255 fn1,fn2
      real*8 Ogap,Egap,Wss,score
      real*8 R(MAX_seq_len,MAX_seq_len)
      real*8 H(0:MAX_seq_len,0:MAX_seq_len)
      real*8 E(0:MAX_seq_len,0:MAX_seq_len)
      real*8 F(0:MAX_seq_len,0:MAX_seq_len)
      integer B(2,MAX_seq_len),i,j,trackI,l,m
      integer align_type
      
      align_type = 2
      aseq1='';
      aseq2='';
      Egap = 0.
      Ogap = 1. 
      call GetArg(1,fn1)
      call GetArg(2,fn2)
      call readseq(seq1,fn1)
      call readseq(seq2,fn2)
      i = len_trim(fn1)
      write (*,*) fn1(1:i)
c      i = len_trim(seq2)
c      write (*,*) seq2(1:i)
c      print *,seq1(1:len_trim(seq1));
c      print *,seq2(1:len_trim(seq2));
c      print *,"========="
      call simple_scoring_matrix(seq1,seq2,
     &           len_trim(seq1),len_trim(seq2),R)
      call align(len_trim(seq1),len_trim(seq2),
     &           R,align_type,Ogap,Egap,H,E,F)
      call trackback(align_type,len_trim(seq1),len_trim(seq2),
     &     R,H,E,F,Ogap,Egap,B,trackI)
        j = 0
        do i=trackI,1,-1
          j = j + 1
          l = B(1,i)
          m = B(2,i)
          write (*,*) l,"\t",m 
          if (l .gt. 0) then
            aseq1(j:j) = seq1(l:l)
          elseif (l .eq. 0) then
            aseq1(j:j) = '-'
          endif 

          if (m .gt. 0) then
            aseq2(j:j) = seq2(m:m)
          elseif (m .eq. 0) then
            aseq2(j:j) = '-'
          endif 
        enddo
c       write (*,*) aseq1(1:len_trim(aseq1))
c       write (*,*) aseq2(1:len_trim(aseq2))
c       write (*,*) "hello"
      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readseq(seq,fn)
      implicit none
      character*2000 seq
      character*255 fn
      integer i,IO_status
      character*200 ts

      seq ="" 
      i =len_trim(fn)
      open (10,file=fn(1:i),status="unknown",IOSTAT=IO_status)
      do while (IO_status .eq. 0)
        read (10,'(A200)',IOSTAT=IO_status) ts
        if ((ts(1:1) .ne. '>') .and. (IO_status .eq.0)) then
          i = len_trim(ts)
          seq = seq(1:len_trim(seq))//ts(1:i)
        endif
      enddo
c     print *,seq(1:len_trim(seq))
      close (10)
      end
