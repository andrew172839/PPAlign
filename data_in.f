c      file format explaination:
c      PSSM and FREQ file are actually the same,derived from psiBLAST:
c      file has 44 columns:
c      1st column is the residue index
c      2nd column is the residue name
c      3rd - 22nd columns is the PSSM matrix, similar to the 
c         BLOSUM matrix, which is read in by read_PSSM subroutine
c      23rd - 42nd columns is the frequece matrix, read in by read_FREQ
c      43rd - 44th columns: ???
      
c      secondary structure file format: 
c      this file only has 1 column
c      the first line is the number of residues
c      the n+1th line is the secodary structure of nth residue
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_esa(str_esa,seq_len,filename)
      implicit none
      character*255 filename
      integer seq_len,IO_status,i
      real*8 str_esa(*)

      i=len_trim(filename)
      write (*,*) filename(1:i)
      open (10,file=filename(1:i),status="old",IOSTAT=IO_status)
      write (*,*) filename(1:i),IO_status
      seq_len=0
      do while (IO_status .eq. 0)
        seq_len = seq_len + 1
        read (10,'(10X,F11.4)',IOSTAT=IO_status) str_esa(seq_len)
      enddo
      seq_len=seq_len-1
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_PSSM(PSSM,seq,seq_len,filename)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      integer PSSM(20,MAX_seq_len)
      character*1 seq(MAX_seq_len)
      character*255 filename
      integer seq_len,IO_status,i

      seq_len = 0
      i = len_trim(filename)
      open (10,file=filename(1:i),status="old",IOSTAT=IO_status)
      do while (IO_status .eq. 0)
        seq_len = seq_len + 1
        read (10,'(6X,A1,2X,20I3)',IOSTAT=IO_status)
     &       seq(seq_len),(PSSM(i,seq_len),i=1,20)
      enddo
      close (10)
      seq_len = seq_len - 1
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trans_index(wu2ncbi)
      implicit none
      integer i
      integer wu2ncbi(20)
      character ncbilist*20,wulist*20,res*1

      ncbilist='ARNDCQEGHILKMFPSTWYV'
      wulist  ='GAVILSTDNEQKRCMFYWHP' 
      do i=1,20
        wu2ncbi(i)=index(ncbilist,wulist(i:i))
      enddo
      return
      end subroutine trans_index

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_str_profile(PROF,seq,seq_len,filename)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      character*1 seq(MAX_seq_len)
      character*3 res
      real*8 PROF(20,MAX_seq_len),temp(20),test(20),Psum
      integer IO_status,i,seq_len,wu2ncbi(20),j
      character*255 filename
      character*20 res_list 
      integer res2index

      res_list ='GAVILSTDNEQKRCMFYWHP' 
      call trans_index(wu2ncbi)
      i=len_trim(filename)
      seq_len=0
      open (10,file=filename(1:i),status="old",IOSTAT=IO_status)
      do while (IO_status .eq. 0)
        seq_len = seq_len + 1
        read (10,'(5X,A3,20F7.3)',IOSTAT=IO_status) res,(temp(i),i=1,20)
        i=res2index(res)
        seq(seq_len) = res_list(i:i)
        Psum=0.
        do i=1,20
          Psum=Psum+temp(i)
        enddo
        if (Psum .eq. 0) then
          do i=1,20
             PROF(i,seq_len)=0.05
          enddo
        else 
          do i=1,20
            j=wu2ncbi(i)
            PROF(j,seq_len)=temp(i)/Psum
          enddo
        endif
      enddo
      close(10)
      seq_len=seq_len-1
      return
      end subroutine read_str_profile

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_FREQ(FREQ,seq,seq_len,filename)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      integer FREQ(20,MAX_seq_len)
      character*1 seq(MAX_seq_len)
      character*255 filename,BLOSUMfile
      integer seq_len,IO_status,i,j
      integer BLOSUM_matrix(20,20),freq_sum,resindex
      character*20 res_list

      res_list='ARNDCQEGHILKMFPSTWYV'
      BLOSUMfile='~/PPAlign'
      call read_BLOSUM(BLOSUM_matrix,BLOSUMfile)
      seq_len = 0
      i = len_trim(filename)
      open (10,file=filename(1:i),status="old",IOSTAT=IO_status)
      do while (IO_status .eq. 0)
        seq_len = seq_len + 1
        read (10,'(6X,A1,63X,20I4)',IOSTAT=IO_status)
     &       seq(seq_len),(FREQ(i,seq_len),i=1,20)
        freq_sum = 0
        do i = 1, 20
          freq_sum = freq_sum + FREQ(i,seq_len)
        enddo
        if (freq_sum .eq. 0) then
          resindex = index(res_list,seq(seq_len))
          if (resindex .ne. 0) then
            do i = 1, 20
              if (BLOSUM_matrix(resindex,i) .ge. 0) then
                FREQ(i,seq_len) = 2 ** BLOSUM_matrix(resindex,i)
              endif
            enddo 
          elseif (resindex .eq. 0) then
            do i = 1, 20
                FREQ(i,seq_len) = 0.05
            enddo 
          endif
        endif
      enddo
      close (10)
      seq_len = seq_len - 1
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_ss(SS,seq_len,filename)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      character*1 SS(MAX_seq_len)
      integer seq_len,i
      character*255 filename
      i = len_trim(filename)
      open (10,file=filename(1:i),status='old')
      read (10,*) seq_len
      do i=1,seq_len
        read (10,*) SS(i)
      enddo
      close (10) 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*3 function index2res(i)
      integer i
      character*80 res_list
      character*40 list1,list2

c               1   2   3   4   5   6   7   8   9   10
      list1 = 'GLY ALA VAL ILE LEU SER THR ASP ASN GLU ' 
c               11  12  13  14  15  16  17  18  19  20
      list2 = 'GLN LYS ARG CYS MET PHE TYR TRP HIS PRO '
      res_list = list1//list2
      index2res = res_list(i*4-3:i*4-1) 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function res2index(res)
      implicit none
      character*3 res
      character*80 res_list
      character*40 list1,list2

c               1   2   3   4   5   6   7   8   9   10
      list1 = 'GLY ALA VAL ILE LEU SER THR ASP ASN GLU ' 
c               11  12  13  14  15  16  17  18  19  20
      list2 = 'GLN LYS ARG CYS MET PHE TYR TRP HIS PRO '
      res_list = list1//list2
      res2index = (index(res_list,res(1:3))+3)/4
      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function r2index(res)
      implicit none
      character*1 res
      character*20 res_list

c               1   2   3   4   5   6   7   8   9   10
c      list1 = 'GLY ALA VAL ILE LEU SER THR ASP ASN GLU ' 
c               11  12  13  14  15  16  17  18  19  20
c      list2 = 'GLN LYS ARG CYS MET PHE TYR TRP HIS PRO '
      res_list ='GAVILSTDNEQKRCMFYWHP' 
      r2index = index(res_list,res(1:1))
      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_BLOSUM(BLOSUM_matrix,filename)
      implicit none
      character*255 filename
      integer BLOSUM_matrix(20,20)
      integer i,j,IO_status

      i = len_trim(filename)
      open (10,file=filename(1:i),status="old",IOSTAT=IO_status)
      do i = 1,20
      read (10,'(1x,20I3)',IOSTAT=IO_status)(BLOSUM_matrix(i,j),j=1,20)
      enddo
      close (10)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_list(filename,str_list,fN)
      implicit none
      character*255 filename,str_list(*)
      integer fN,IO_status,i

      i=len_trim(filename)
      fN=0
      open (10,file=filename(1:i),status="old",IOSTAT=IO_status)
      do while(IO_status.eq.0)
        fN=fN+1 
        read(10,*,IOSTAT=IO_status),str_list(fN)
      enddo
      fN=fN-1
      close(10)
      return
      end
