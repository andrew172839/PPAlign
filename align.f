c     align_type = 1, 2, 3, 4
c     1:  local alignment
c     2:  global alignment
c     3:  local-global alignment
c         SEQqi will be penalized at both ends
c         SEQdj will not be penalized at either ends
c     4:  semi-global alignment
c         both ends will not be penalized 

c     trackback find a route from lower-right of matrix to 
c     upleft of matrix
c     The trackback info is stored in the 2XtrackI maxtrix B
c     Alignment of qi and dj
c     The first row contains the info of qi
c     The second row contains the info of dj
c     for the alignmen or gap:
c     qi  B(1,trackI) = i    qi B(1,trackI) = i  -  B(1,trackI) = 0
c     |                      |                   |
c     dj  B(2,trackI) = j    -  B(2,trackI) = 0  dj B(2,trackI) = j
c     alignment              gap                 gap

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trackback(align_type,si_len,sj_len,
     &                      R,H,E,F,Ogap,Egap,B,trackI)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      integer si_len,sj_len,align_type
      real*8 R(MAX_seq_len,MAX_seq_len)
      real*8 H(0:MAX_seq_len,0:MAX_seq_len)
      real*8 E(0:MAX_seq_len,0:MAX_seq_len)
      real*8 F(0:MAX_seq_len,0:MAX_seq_len)
      real*8 Ogap,Egap
      integer isequal
      integer B(2,MAX_seq_len)
      integer trackI 
      integer i,j,ii,jj

      trackI = 0
      call end_select2(align_type,i,j,H,si_len,sj_len,Ogap,Egap)
c      if ( align_type .eq. 3) then
c        jj = sj_len
c        do while (j < jj)
c          trackI = trackI + 1
c          B(1,trackI) = 0
c          B(2,trackI) = jj
c          jj = jj - 1
c        enddo
c      endif
      do ii=si_len,i+1,-1
         trackI=trackI+1
         B(1,trackI) = ii 
         B(2,trackI) = 0
      enddo
      do jj=sj_len,j+1,-1
         trackI=trackI+1
         B(1,trackI)=jj
         B(2,trackI)=0
      enddo
      do while ((i > 0) .and. (j > 0))
        trackI = trackI + 1
        if (isequal(H(i,j),H(i-1,j-1)+R(i,j)) .eq. 1) then
          B(1,trackI) = i
          B(2,trackI) = j
          i = i - 1
          j = j - 1
        elseif (isequal(H(i,j), E(i,j)) .eq. 1) then
          B(1,trackI) = i
          B(2,trackI) = 0
          i = i - 1
        elseif (isequal(H(i,j), F(i,j)) .eq. 1) then
          B(1,trackI) = 0
          B(2,trackI) = j
          j = j - 1
        endif
      enddo
      if (i .eq. 0) then
        do while (j .gt. 0)
          trackI = trackI + 1
          B(1,trackI) = 0
          B(2,trackI) = j
          j = j - 1
        enddo
      elseif (j .eq. 0) then
        do while (i .gt. 0)
          trackI = trackI + 1
          B(1,trackI) = i
          B(2,trackI) = 0
          i = i - 1
        enddo
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine align(si_len,sj_len,R,align_type,Ogap,Egap,H,E,F)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      integer si_len,sj_len,align_type
      real*8 R(MAX_seq_len,MAX_seq_len)
      real*8 Ogap,Egap
      real*8 H(0:MAX_seq_len,0:MAX_seq_len)
      real*8 E(0:MAX_seq_len,0:MAX_seq_len)
      real*8 F(0:MAX_seq_len,0:MAX_seq_len)
      real*8 HE,HG,HF
      integer i,j,k,ii,jj
      real*8 tt

      call init_EFH(E,F,H,si_len,sj_len,Ogap,Egap,align_type)
      do j = 1,sj_len
        do i = 1,si_len
          HE = E(i-1,j) + Egap 
          HG = H(i-1,j) + Ogap + Egap
          HF = F(i-1,j) + Ogap + Egap
          E(i,j) = DMin1(HE,HG,HF)
          HE = E(i,j-1) + Ogap + Egap
          HG = H(i,j-1) + Ogap + Egap
          HF = F(i,j-1) + Egap
          F(i,j) = DMin1(HE,HG,HF)
          H(i,j) = DMin1(E(i,j),F(i,j),H(i-1,j-1)+R(i,j))
        enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_EFH(E,F,H,si_len,sj_len,Ogap,Egap,align_type)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      real*8 INF
      parameter (INF=20000.)
      real*8 E(0:MAX_seq_len,0:MAX_seq_len)
      real*8 F(0:MAX_seq_len,0:MAX_seq_len)
      real*8 H(0:MAX_seq_len,0:MAX_seq_len)
      integer si_len,sj_len,align_type
      real*8 Ogap,Egap
      integer i,ii,jj

      if ((align_type .eq. 1) .or. (align_type .eq. 4)) then
        H(0,0) = 0.
        do i = 1,si_len
          E(i,0) = INF 
          F(i,0) = INF 
          H(i,0) = 0. 
        enddo
        do i = 1,sj_len
          E(0,i) = INF 
          F(0,i) = INF
          H(0,i) = 0. 
        enddo
      elseif (align_type .eq. 2) then
        H(0,0) = 0.
        do i = 1,si_len
          E(i,0) = INF 
          F(i,0) = INF 
          H(i,0) = Ogap + i * Egap 
        enddo
        do i = 1,sj_len
          E(0,i) = INF
          F(0,i) = INF
          H(0,i) = Ogap + i * Egap 
        enddo
      elseif (align_type .eq. 3) then
        H(0,0) = 0.
        do i = 1,si_len
          E(i,0) = INF 
          F(i,0) = INF 
          H(i,0) = 0 
        enddo
        do i = 1,sj_len
          E(0,i) = INF
          F(0,i) = INF
          H(0,i) = Ogap + i * Egap 
        enddo
      endif 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Scoring_Matrix(R,seq_FREQ,seq_SS,seq_len,
     &     struct_PSSM,struct_SS,struct_len,Wss)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      real*8 R(MAX_seq_len,MAX_seq_len)
      integer seq_FREQ(20,MAX_seq_len)
      character*1 seq_SS(MAX_seq_len)
      integer seq_len
      integer struct_PSSM(20,MAX_seq_len)
      character*1 struct_SS(MAX_seq_len)
      integer struct_len
      integer i,j,k,FREQ_sum
      real*8 Wss,Umutation,U2ndary
      character*1 seq(MAX_seq_len)

      do j = 1,seq_len
        do i = 1,struct_len
          Umutation = 0.
          U2ndary = 0.
          FREQ_sum = 0
          do k = 1,20
            FREQ_sum = FREQ_sum + seq_FREQ(k,j)
          enddo
          if (FREQ_sum .eq. 0) then
            do k = 1,20
              seq_FREQ(k,j) = 0.05
            enddo
            FREQ_sum = 1.
          endif
          do k = 1,20
            Umutation = Umutation - seq_FREQ(k,j) * struct_PSSM(k,i)
          enddo
          Umutation = Umutation / FREQ_sum
          if (seq_SS(j) .eq. struct_SS(i)) then
            U2ndary = -1.
          elseif (seq_SS(j) .ne. struct_SS(i)) then
            U2ndary = 1.
          endif
          R(i,j) = Umutation + Wss * U2ndary
        enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine New_Scoring_Matrix(R,seq_PSSM,seq_FREQ,seq_SS,seq_len,
     &   struct_PSSM,struct_PROF,struct_SS,struct_len,Wss,Wprof,Wshift)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      real*8 R(MAX_seq_len,MAX_seq_len)
      integer seq_FREQ(20,MAX_seq_len),seq_PSSM(20,MAX_seq_len)
      character*1 seq_SS(MAX_seq_len)
      integer seq_len
      integer struct_PSSM(20,MAX_seq_len)
      real*8 struct_PROF(20,MAX_seq_len)
      character*1 struct_SS(MAX_seq_len)
      integer struct_len
      integer i,j,k,FREQ_sum
      real*8 PROF_sum
      real*8 Wss,Wprof,Wshift,Umutation,U2ndary,Ustrmut
      character*1 seq(MAX_seq_len)

      do j = 1,seq_len
        do i = 1,struct_len
          Umutation = 0.
          Ustrmut = 0.
          U2ndary = 0.
          FREQ_sum = 0
          PROF_sum = 0.
          do k = 1,20
            FREQ_sum = FREQ_sum + seq_FREQ(k,j)
          enddo
          if (FREQ_sum .eq. 0) then
            do k = 1,20
              seq_FREQ(k,j) = 0.05
            enddo
            FREQ_sum = 1.
          endif
          do k = 1,20
            PROF_sum = PROF_sum + struct_PROF(k,i)
          enddo
          if (PROF_sum .eq. 0) then
            do k = 1,20
              struct_PROF(k,i) = 0.05
            enddo
            PROF_sum = 1.
          endif
          do k = 1,20
            Umutation=Umutation-seq_FREQ(k,j)*struct_PSSM(k,i)
            Ustrmut=Ustrmut-seq_PSSM(k,j)*struct_PROF(k,i)
c            write (*,*) 'prof',seq_PSSM(k,j),struct_PROF(k,i)
          enddo
          Umutation = Umutation / FREQ_sum
          Ustrmut = Ustrmut / PROF_sum
          if (seq_SS(j) .eq. struct_SS(i)) then
            U2ndary = -1.
          elseif (seq_SS(j) .ne. struct_SS(i)) then
            U2ndary = 1.
          endif
c          write (*,*) Umutation,Ustrmut,U2ndary
          R(i,j) = (1.-Wprof)*Umutation+Wprof*Ustrmut+Wss*U2ndary-Wshift
        enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      subroutine simple_scoring_matrix(q,d,q_len,d_len,R)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      character*256 q,d
      integer q_len,d_len
      real*8 R(MAX_seq_len,MAX_seq_len)
      integer i,j

      do j = 1,d_len
        do i = 1,q_len 
          if (d(j:j) .eq. q(i:i)) then
            R(i,j) = -1.
          elseif (d(j:j) .ne. q(i:i)) then
c            if (d(j:j).eq.'M'.and.q(i:i).eq.'X') then
c              R(i,j)=-1.
c           # elseif (d(j:j).eq.'X'.and.q(i:i).eq.'M') then
c              R(i,j)=-1.
c            else
              R(i,j) = 1000.
c            endif
          endif
        enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine end_select(align_type,ii,jj,H,si_len,sj_len)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      integer align_type,ii,jj,si_len,sj_len,ii2,jj2
      real*8 H(0:MAX_seq_len,0:MAX_seq_len),minv
      integer i,j

      if (align_type .eq. 1) then
        ii = si_len
        jj = sj_len
        minv = H(ii,jj)
        do j = sj_len,1,-1
          do i = si_len,1,-1
            if (minv .gt. H(i,j)) then
              ii = i
              jj = j
              minv = H(ii,jj)
            endif
          enddo
        enddo 
      elseif (align_type .eq. 2) then
        ii = si_len
        jj = sj_len
      elseif (align_type .eq. 3) then
        ii = si_len
        jj = sj_len 
        minv = H(ii,jj)
        do i = si_len,1,-1
          if (minv .gt. H(i,jj)) then
            ii = i
            minv = H(ii,jj)
          endif
        enddo
      elseif (align_type .eq. 4) then
        ii = si_len
        jj = sj_len
        minv = H(ii,jj)
        do i = si_len,1,-1
          if (minv .gt. H(i,sj_len)) then
            ii = i
            minv = H(ii,jj)
          endif
        enddo
        ii2 = si_len
        jj2 = sj_len
        do j = sj_len,1,-1
          if (minv .gt. H(si_len,j)) then
            jj2 = j
            minv = H(ii2,jj2)
          endif
        enddo
        if (H(ii2,jj2) .lt. H(ii,jj)) then
          ii = ii2
          jj = jj2
        endif
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine end_select2(align_type,ii,jj,H,si_len,sj_len,Ogap,Egap)
      implicit none
      integer MAX_seq_len
      parameter (MAX_seq_len=2000)
      integer align_type,ii,jj,si_len,sj_len,ii2,jj2
      real*8 H(0:MAX_seq_len,0:MAX_seq_len),minv,Ogap,Egap,minv0
      integer i,j

      if (align_type .eq. 1) then
        ii = si_len
        jj = sj_len
        minv = H(ii,jj)
        do j = sj_len,1,-1
          do i = si_len,1,-1
            if (minv .gt. H(i,j)) then
              ii = i
              jj = j
              minv = H(ii,jj)
            endif
          enddo
        enddo 
      elseif (align_type .eq. 2) then
        ii = si_len
        jj = sj_len
      elseif (align_type .eq. 3) then
        ii = si_len
        jj = sj_len 
        minv = H(ii,jj)
        do i = si_len,1,-1
          do j=sj_len,1,-1
            if (j.eq.sj_len) then
              minv0=H(i,j)
            elseif(j.ne.sj_len) then
              minv0=H(i,j)+Ogap+(sj_len-j)*Egap
            endif
          if (minv .gt. minv0) then
            ii = i
            jj = j
            minv = minv0 
          endif
        enddo
      enddo
      elseif (align_type .eq. 4) then
        ii = si_len
        jj = sj_len
        minv = H(ii,jj)
        do i = si_len,1,-1
          if (minv .gt. H(i,sj_len)) then
            ii = i
            minv = H(ii,jj)
          endif
        enddo
        ii2 = si_len
        jj2 = sj_len
        do j = sj_len,1,-1
          if (minv .gt. H(si_len,j)) then
            jj2 = j
            minv = H(ii2,jj2)
          endif
        enddo
        if (H(ii2,jj2) .lt. H(ii,jj)) then
          ii = ii2
          jj = jj2
        endif
      endif
      return
      end
