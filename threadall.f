cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program threadall
      implicit none
      integer MAX_seq_len,MAX_F
      parameter (MAX_seq_len=2000,MAX_F=20000)
      character*255 data_path,flist(MAX_F),flistname
      integer PSSM(20,MAX_seq_len),fN
      real*8  PROF(20,MAX_seq_len)
      character*1 struct_SS(MAX_seq_len)
      character*1 struct_seq(MAX_seq_len)
      integer struct_len,seq_len
      integer FREQ(20,MAX_seq_len),seqPSSM(20,MAX_seq_len)
      character*1 seq_SS(MAX_seq_len)
      character*1 seq(MAX_seq_len)
      character*255 DP 
      character*255 str_prefix,seq_prefix
      character*255 file_PSSM,file_str_SS,file_FREQ,file_seq_SS
      character*255 file_str_PROF,file_seq_strPSSM
      real*8 R(MAX_seq_len,MAX_seq_len)
      real*8 Ogap,Egap,Wss,Wprof,Wshift,score
      real*8 H(0:MAX_seq_len,0:MAX_seq_len)
      real*8 E(0:MAX_seq_len,0:MAX_seq_len)
      real*8 F(0:MAX_seq_len,0:MAX_seq_len)
      integer B(2,MAX_seq_len)
      character*2000 align1,align2
      integer align_type,istrack
      integer i,j,ii,jj,trackI,k,ali_start,ali_end
      real*8 str2real
      integer str2int
      character*10 arg_Ogap,arg_Egap,arg_Wss,arg_align,
     &             arg_istrack,arg_Wprof,arg_Wshift

      DP='~/PPAlign'
      call GetArg(1,flistname)
      call GetArg(2,seq_prefix)
      call GetArg(3,arg_Ogap)
      call GetArg(4,arg_Egap)
      call GetArg(5,arg_Wss)
      call GetArg(6,arg_align)
      call GetArg(7,arg_Wprof)
      call GetArg(8,arg_Wshift)
      Ogap = str2real(arg_Ogap,10)
      Egap = str2real(arg_Egap,10)
      Wss = str2real(arg_Wss,10)
      align_type = str2int(arg_align,10)
      Wprof = str2real(arg_Wprof,10)
      Wshift = str2real(arg_Wshift,10)
      j = Len_Trim(seq_prefix)
      file_FREQ = seq_prefix(1:j)//'.freq'
      file_seq_SS = seq_prefix(1:j)//'.ss'
      file_seq_strPSSM = seq_prefix(1:j)//'.strPSSM'
      call read_FREQ(FREQ,seq,seq_len,file_FREQ)
      call read_ss(seq_SS,seq_len,file_seq_SS)
c      call read_PSSM(seqPSSM,seq,seq_len,file_FREQ)
      call read_PSSM(seqPSSM,seq,seq_len,file_seq_strPSSM)
      call read_list(flistname,flist,fN);
      do k=1,fN
        i=len_trim(DP)
        j=Len_trim(flist(k))
        file_PSSM = DP(1:i)//'pssm/'//flist(k)(1:j)//'.pssm'
        file_str_SS = DP(1:i)//'ss/'//flist(k)(1:j)//'.ss'
        file_str_PROF = DP(1:i)//'sol/'//flist(k)(1:j)//'.sol'
        call read_str_profile(PROF,struct_seq,struct_len,file_str_PROF)
        call read_PSSM(PSSM,struct_seq,struct_len,file_PSSM)
        call read_ss(struct_SS,struct_len,file_str_SS)
c        call Scoring_Matrix(R,FREQ,seq_SS,seq_len,
c     &       PSSM,struct_SS,struct_len,Wss)
        call New_Scoring_Matrix(R,seqPSSM,FREQ,seq_SS,seq_len,
     &       PSSM,PROF,struct_SS,struct_len,Wss,Wprof,Wshift)
        call align(struct_len,seq_len,R,align_type,Ogap,Egap,H,E,F)
        call trackback(align_type,struct_len,seq_len,
     &     R,H,E,F,Ogap,Egap,B,trackI)
        do i=trackI,1,-1
          if (B(1,i)*B(2,i).ne.0) then
            ali_start = B(2,i)
            goto 101
          endif
        enddo
 101    continue;
        do i=1,trackI
          if (B(1,i)*B(2,i).ne.0) then
            ali_end = B(2,i)
            goto 102
          endif
        enddo
 102    continue;
        call end_select2(align_type,ii,jj,H,struct_len,
     &                   seq_len,Ogap,Egap)
        score = H(ii,jj)
        write(*,'(A10,1F10.3,I5)'),flist(k)(1:j),score,ali_end-ali_start
      enddo
      stop
      end
