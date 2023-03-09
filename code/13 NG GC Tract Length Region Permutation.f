      program nggctl
c
c-----NG GC Tract Length Region Permutation
c-----Uses a permutation test to determine whether different regions
c-----have different gene conversion tract lengths
c
c-----Written by Stephen W. Schaeffer
c-----Date: 22 June 2022
c
c-----Source: NG GC Tract Length Region Permutation.f
c
      integer ia,ib,ic,id,ie,ig,ngc,itmp
      integer nbeg(35000),nend(35000),nlen(35000),nreg(35000)
      integer ireg(14),nb(14),ne(14)
      real sgc(35000),amdo(14),amdp(14),aprob(2,14),smin(14),smax(14)
      integer nsec
      character*2 arr
      character*3 a
c
c-----nsec  - random seed
c-----iseed - variable to save original random seed to print with the results
c-----Initialize nsec, seconds past midnight
c
     	nsec=-1*int(secnds(0.0)*10)
     	iseed=nsec
c
c-----chr() - array of nucleotide positions on the third chromosome map
c-----Initialize chr()
c
      write(*,'("Input the inversion.")')
      read(*,*) arr
c
c-----Open file with gene conversion tracts
c
      open(unit=1,file='GC_'//arr//'_Tracts.tsv')
      read(1,'(1x)')
c
c-----ngc      - number of gene conversion tracts
c-----a(i)     - subregion of the ith gene conversion tract
c-----nbeg(i)  - coordinate of the first nucleotide of the ith gene conversion tract
c-----nend(i)  - coordinate of the last  nucleotide in the ith gene conversion tract
c-----nlen(i)  - length of the ith gene conversion tract
c
c-----Read in the ngc coordinate pairs nbeg, nend until End Of File
c
      ngc=0      
50    ngc=ngc+1
      read(1,*,end=100) a,ib,ic,id,ie,ig,nbeg(ngc),nend(ngc)
      nlen(ngc)=nend(ngc)-nbeg(ngc)+1
      read(a(1:2),*) nreg(ngc)
      goto 50      
100   close(unit=1)
      ngc=ngc-1
c
c-----Estimate the median for each region
c-----ireg(i) - count of gene conversion tracts in the ith region
c
      ireg=0
      do i=1,ngc
      ireg(nreg(i))=ireg(nreg(i))+1
      end do
      nb(1)=1
      ne(1)=ireg(1)
      do i=2,14
      nb(i)=ne(i-1)+1
      ne(i)=ne(i-1)+ireg(i)
      end do
c
c-----Estimate gene conversion tract median length
c-----sgc(i)  - length of the ith gene conversion tract, real variable
c-----amdo(i) - median gene conversion tract length in the ith region
c
      do i=1,14
      na=0
      do j=nb(i),ne(i)
      na=na+1
      sgc(na)=float(nlen(j))
      end do
      call sort(na,sgc)
      if(mod(na,2).eq.0) then
      amdo(i)=sgc(na/2)
      else
      amdo(i)=(sgc((na-1)/2)+(sgc((na+1)/2)))/2.
      end if
      end do      
c
c-----Random permutation test 10,000 replicates
c-----Begin random permutation
c-----amdp(i)    - median gene conversion tract length in the ith region from permuted data
c-----aprob(i,j) - number of permutations that are less than observed i=1
c-----             or greater than observed i=2 for the j regions
c    	
      aprob=0
      smin=1000000.
      smax=0.
      do isim=1,10000
c
c-----Shuffle the gene conversion tracts in the list without replacement
c
      do i=ngc,1,-1
      ia=int(ran1(nsec)*i)+1
      itmp=nlen(ia)
      nlen(ia)=nlen(i)
      nlen(i)=itmp
      end do
c
c-----Estimate the median in the 14 regions
c
      do i=1,14
      na=0
      do j=nb(i),ne(i)
      na=na+1
      sgc(na)=float(nlen(j))
      end do
      call sort(na,sgc)
      if(mod(na,2).eq.0) then
      amdp(i)=sgc(na/2)
      else
      amdp(i)=(sgc((na-1)/2)+(sgc((na+1)/2)))/2.
      end if
      end do
      do i=1,14
      if(amdp(i).lt.amdo(i)) aprob(1,i)=aprob(1,i)+1.
      if(amdp(i).ge.amdo(i)) aprob(2,i)=aprob(2,i)+1.
      if(amdp(i).gt.smax(i)) smax(i)=amdp(i)
      if(amdp(i).lt.smin(i)) smin(i)=amdp(i)
      end do
      end do
      do i=1,14
      aprob(1,i)=aprob(1,i)/10000.
      aprob(2,i)=aprob(2,i)/10000.
      end do      
c
c-----End random permutation
c-----Output results
c
      open(unit=1,file='GC_'//arr//'_Tract_Length_Permut.txt')
      write(1,'("Analysis of ",a2," Gene Conversion Tracts")') arr
      write(1,'("Testing homoegeneity of tract lengths among regions")')
      write(1,'("Method: Permute Tract Assignment to Region")')
      write(1,'(1x)')
      write(1,'(1x)')
      write(1,'("Region  Observed P(Permuted<Obs) P(Permuted>Obs)")')
      do i=1,14
      write(1,'(i2,6x,f8.1,2f16.5)') i,amdo(i),aprob(1,i),aprob(2,i)
      end do 
      write(1,'(1x)')
      write(1,'(1x)')
      write(1,'("Region  Minimum Maximum")')
      do i=1,14
      write(1,'(i2,5x,2f8.1)') i,smin(i),smax(i)
      end do
      write(1,'(1x)')
      write(1,'("Random seed = ",i15)') iseed      
      close(unit=1)
      stop
      end
c
c-----Sort
c
      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
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
