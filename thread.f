cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program thread
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      character*255 data_path
      integer PSSM(20,MAX_seq_len)
      real*8 PROF(20,MAX_seq_len)
      character*1 struct_SS(MAX_seq_len)
      character*1 struct_seq(MAX_seq_len)
      integer struct_len
      integer FREQ(20,MAX_seq_len),seqPSSM(20,MAX_seq_len)
      character*1 seq_SS(MAX_seq_len)
      character*1 seq(MAX_seq_len)
      character*255 DP 
      integer seq_len
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
      integer i,j,ii,jj,trackI
      real*8 str2real
      integer str2int
      character*10 arg_Ogap,arg_Egap,arg_Wss,arg_align,
     &             arg_istrack,arg_Wprof,arg_Wshift

      DP='~/PPAlign'
      call GetArg(1,str_prefix)
      call GetArg(2,seq_prefix)
      call GetArg(3,arg_Ogap)
      call GetArg(4,arg_Egap)
      call GetArg(5,arg_Wss)
      call GetArg(6,arg_align)
      call GetArg(7,arg_istrack)
      call GetArg(8,arg_Wprof)
      call GetArg(9,arg_Wshift)
      Ogap = str2real(arg_Ogap,10)
      Egap = str2real(arg_Egap,10)
      Wss = str2real(arg_Wss,10)
      align_type = str2int(arg_align,10)
      istrack = str2int(arg_istrack,10)
      Wprof = str2real(arg_Wprof,10)
      Wshift = str2real(arg_Wshift,10)
      i=len_trim(DP)
      j=Len_trim(str_prefix)
      file_PSSM = DP(1:i)//'pssm/'//str_prefix(1:j)//'.pssm'
      file_str_SS = DP(1:i)//'ss/'//str_prefix(1:j)//'.ss'
      file_str_PROF = DP(1:i)//'sol/'//str_prefix(1:j)//'.sol'
      j = Len_Trim(seq_prefix)
      file_FREQ = seq_prefix(1:j)//'.freq'
      file_seq_SS = seq_prefix(1:j)//'.ss'
      file_seq_strPSSM = seq_prefix(1:j)//'.strPSSM'
      call read_str_profile(PROF,struct_seq,struct_len,file_str_PROF)
      call read_PSSM(PSSM,struct_seq,struct_len,file_PSSM)
c      call read_FREQ(PROF,struct_seq,struct_len,file_PSSM)
      call read_ss(struct_SS,struct_len,file_str_SS)
c      call read_PSSM(seqPSSM,seq,seq_len,file_FREQ)
      call read_PSSM(seqPSSM,seq,seq_len,file_seq_strPSSM)
      call read_FREQ(FREQ,seq,seq_len,file_FREQ)
      call read_ss(seq_SS,seq_len,file_seq_SS)
c      call Scoring_Matrix(R,FREQ,seq_SS,seq_len,
c     &       PSSM,struct_SS,struct_len,Wss)
      call New_Scoring_Matrix(R,seqPSSM,FREQ,seq_SS,seq_len,
     &       PSSM,PROF,struct_SS,struct_len,Wss,Wprof,Wshift)
      call align(struct_len,seq_len,R,align_type,Ogap,Egap,H,E,F)
      if (istrack .gt. 0) then
        call trackback(align_type,struct_len,seq_len,
     &     R,H,E,F,Ogap,Egap,B,trackI)
        do i=trackI,1,-1
          write (*,*) B(1,i),"   ",B(2,i) 
        enddo
      else 
        call end_select2(align_type,ii,jj,H,struct_len,
     &   seq_len,Ogap,Egap)
        score = H(ii,jj)
        print *,score
      endif
      stop
      end
