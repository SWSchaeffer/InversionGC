      program ngrdrp
c
c-----NG Reg Don Rec Permut
c-----Uses a permutation test to determine whether Donor Recipient
c-----numbers of gene conversion events are similar among different
c-----regions
c
c-----Written by Stephen W. Schaeffer
c-----Date: 29 June 2022
c
c-----Source: NG Reg Don Rec Permute.f
c      
      integer nsec,iseed,ndar(100000),nrar(100000),nreg(100000)
      integer notm(15,6,6),nptm(15,6,6),nper,iflag
      integer imn(15,6,6),imx(15,6,6)
      real apl(15,6,6),apg(15,6,6)
      character*2 dar,rar,da(6),ra(6)
c
c-----nsec  - random seed
c-----iseed - variable to save original random seed to print with the results
c-----Initialize nsec, seconds past midnight
c
     	nsec=-1*int(secnds(0.0)*10)
     	iseed=nsec
c
c-----Initialize variables
c-----da(i) - ith donor arrangement
c-----ra(i) - ith recipient arrangement
c
      da(1)='AR'
      da(2)='CH'
      da(3)='CU'
      da(4)='PP'
      da(5)='ST'
      da(6)='TL'
      ra(1)='AR'
      ra(2)='CH'
      ra(3)='CU'
      ra(4)='PP'
      ra(5)='ST'
      ra(6)='TL'
c-----Input Region, Donor Arrangement, and Recipient Arrangement for each tract
c-----nreg(i) - number of the ith region for each tract
c-----dar     - donor arrangement
c-----rar     - recipient arrangement
c-----ndar(i) - integer value of the     donor arrangement for the ith tract
c-----nrar(i) - integer value of the recipient arrangement for the ith tract
c
     	open(unit=1,file='GC_Region_Donor_Recipient.tsv')
     	nc=0
50    nc=nc+1
      read(1,*,end=100) nreg(nc),dar,rar
      if(dar.eq.'AR') ndar(nc)=1
      if(dar.eq.'CH') ndar(nc)=2
      if(dar.eq.'CU') ndar(nc)=3
      if(dar.eq.'PP') ndar(nc)=4
      if(dar.eq.'ST') ndar(nc)=5
      if(dar.eq.'TL') ndar(nc)=6
      if(rar.eq.'AR') nrar(nc)=1
      if(rar.eq.'CH') nrar(nc)=2
      if(rar.eq.'CU') nrar(nc)=3
      if(rar.eq.'PP') nrar(nc)=4
      if(rar.eq.'ST') nrar(nc)=5
      if(rar.eq.'TL') nrar(nc)=6
     	goto 50
100  	close(unit=1)
      nc=nc-1
c
c-----Count the number of tracts of the ith region, jth donor, and kth recipient
c-----notm(i,j,k) - observed matrix of ith region, jth donor, and kth recipient
c
      notm=0
      do i=1,nc
      notm(nreg(i),ndar(i),nrar(i))=notm(nreg(i),ndar(i),nrar(i))+1
      notm(15,ndar(i),nrar(i))=notm(15,ndar(i),nrar(i))+1
      end do
c
c-----Begin Random permutation
c-----imn(i,j,k) - minimum value of permutations
c-----imx(i,j,k) - maximum value of permutations
c-----apl(i,j,k) - probability of permutations less than observed
c-----apg(i,j,k) - probability of permutations greater than observed
c
      apl=0.
      apg=0.
      imn=100000
      imx=0
      permute: do nper=1,10000
c      write(*,'("Permutation = ",i6)') nper
      if(mod(nper,100).eq.0) write(*,'("Permutation = ",i6)') nper
c
c-----Permute the donor array 
c-----nptm(i,j,k) - observed matrix of ith region, jth donor, and kth recipient
c
      nptm=0
      do i=nc,1,-1
      ia=int(ran1(nsec)*i)+1
      itmp=ndar(ia)
      ndar(ia)=ndar(i)
      ndar(i)=itmp
      end do
c
c-----Permute the recipient array
c-----Constrain to reject self conversion events
c
200   iflag=0
      do i=nc,1,-1
c      ia=int(ran1(nsec)*i)+1
150   ia=int(ran1(nsec)*i)+1
      if(nrar(ia).eq.ndar(i)) then
      iflag=iflag+1
      if(iflag.eq.100) then
c      write(*,'("Restart...")')
      goto 200
      end if
      goto 150
      end if
      iflag=0
      itmp=nrar(ia)
      nrar(ia)=nrar(i)
      nrar(i)=itmp
      end do
c
c-----Re-estimate the matrix
c
      do i=1,nc
      nptm(nreg(i),ndar(i),nrar(i))=nptm(nreg(i),ndar(i),nrar(i))+1
      nptm(15,ndar(i),nrar(i))=nptm(15,ndar(i),nrar(i))+1
      end do
c
c-----Check observed versus permuted matrices
c
      do k=1,15
      do i=1,6
      do j=1,6
      if(nptm(k,i,j).lt.notm(k,i,j)) apl(k,i,j)=apl(k,i,j)+1.
      if(nptm(k,i,j).ge.notm(k,i,j)) apg(k,i,j)=apg(k,i,j)+1.
      if(nptm(k,i,j).lt.imn(k,i,j)) imn(k,i,j)=nptm(k,i,j)
      if(nptm(k,i,j).gt.imx(k,i,j)) imx(k,i,j)=nptm(k,i,j)
      end do
      end do
      end do
      end do permute
      do k=1,15
      do i=1,6
      do j=1,6
      apl(k,i,j)=apl(k,i,j)/10000.
      apg(k,i,j)=apg(k,i,j)/10000.
      end do
      end do
      end do
      open(unit=1,file='GC_Region_Donor_Recipient_Summary.txt')
      do k=1,15
      write(1,'("Region",i3)') k
      write(1,'(4x,6(6x,a2))') (ra(i),i=1,6)
      do i=1,6
      write(1,'(a2,2x,6i8)') da(i),(notm(k,i,j),j=1,6)
      write(1,'("Min",1x,6i8)') (imn(k,i,j),j=1,6)
      write(1,'("Max",1x,6i8)') (imx(k,i,j),j=1,6)
      write(1,'("P<",2x,6f8.4)') (apl(k,i,j),j=1,6)
      write(1,'("P>",2x,6f8.4)') (apg(k,i,j),j=1,6)
      write(1,'(1x)')
      end do
      write(1,'(1x)')
      end do
      write(1,'(1x)')
      write(1,'("Random seed = ",i15)') iseed
      close(unit=1)
      stop
      end
c
c-----Random number generator
c	
	FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
