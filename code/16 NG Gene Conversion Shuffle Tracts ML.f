      program ngcst
c
c-----NG Gene Conversion Shuffle Tracts ML
c-----Uses shuffled shuffled tracts test to determine whether outlier genes
c-----have different amounts of gene conversion then non-outlier genes
c-----The tract length distribution is based on the ML estimates from
c-----Betran et al. (1997) paper
c
c-----Written by Stephen W. Schaeffer
c-----Date: 6 July 2022
c
c-----Source: NG Gene Conversion Gene Shuffle ML.f
c
      integer ib,ic,id,ie,ig,iseed,iflg
      integer chr(19787792),nbeg(50000),nend(50000),nlen(50000),nl
      integer ntb(11212),nte(11212),ntex(11212),ts(11212)
      integer ngc,ncb,nce,nc,nsec,ia,nprob,nr,nrb(14),nre(14),nev(5,14)
      real anc,b,gcba(11212),agc(2),gcx(2),gcx2(2),xb(2),s2(2)
      real sagc(2),sgcx(2),sxb(2),aprob,amin,amax,amlphi(5,14)
      character*2 arr
      character*3 a
      character*11 gname(11212)
c
c-----nsec  - random seed
c-----iseed - variable to save original random seed to print with the results
c-----Initialize nsec, seconds past midnight
c
     	nsec=-1*int(secnds(0.0)*10)
     	iseed=nsec
c
c-----chr() - array of nucleotide positions on the third chromosome map
c-----nbeg  - coordinate of the first nucleotide in the gene conversion tract
c-----nend  - coordinate of the last  nucleotide in the gene conversion tract
c-----Initialize chr()
c
      chr=0
      write(*,'("Input the inversion.")')
      read(*,*) arr
c
c-----nr    - Number of regions
c-----nrb(i) - First nucleotide of the ith region
c-----nre(i) - Last  nucleotide of the ith region
c
      open(unit=1,file=arr//'_MinMax.tsv')
      read(1,'(1x)')
      nr=14
      do i=1,nr
      read(1,*) nrb(i),nre(i)
      end do
      close(unit=1)
c
c-----Open file with gene conversion tracts
c
      open(unit=1,file='GC_'//arr//'_Tracts.tsv')
      read(1,'(1x)')
c
c-----Read in coordinate pair nbeg, nend until End Of File
c
      ngc=0      
50    ngc=ngc+1
      read(1,*,end=100) a,ib,ic,id,ie,ig,nbeg(ngc),nend(ngc)
      nlen(ngc)=nend(ngc)-nbeg(ngc)+1
c
c-----Increment coordinates in chr()
c      
      do i=nbeg(ngc),nend(ngc)
      chr(i)=chr(i)+1
      end do
      goto 50      
100   close(unit=1)
      ngc=ngc-1
c
c-----Extract regions with gene conversions
c-----ncb - site of first nucleotide
c-----nce - site of last  nucleotide
c-----nc  - counts of gene conversion tracts
c
      open(unit=1,file='GC_'//arr//'_Tracts_Output_MLpara.csv')
      write(1,'("Region,Beg,End,Count")')
      nc=0
      nr=0
      iflag=0
      do i=1,19787792
      if(iflag.eq.0.and.chr(i).eq.0) then
      cycle
      elseif(iflag.eq.0.and.chr(i).gt.0) then
      nr=nr+1
      ncb=i
      nc=nc+chr(i)
      iflag=1
      cycle
      elseif(iflag.eq.1.and.chr(i).gt.0) then
      nc=nc+chr(i)
      cycle
      elseif(iflag.eq.1.and.chr(i).eq.0) then
      nce=i-1
      anc=float(nc)/float(nce-ncb+1)
      nc=0
      iflag=0
      write(1,'(i5,",",i10,",",i10,",",f6.3)') nr,ncb,nce,anc
      endif
      end do
      close(unit=1)
c
c-----Estimate gene conversion tract depth in outlier and non-outlier genes
c
c-----ntr      - Number of exons across all genes
c-----gname(i) - Name of the ith transcript exon
c-----ntex(i)  - Number of the ith transcript exon
c-----ntb(i)   - First nucleotide of the ith transcript exon
c-----nte(i)   - Last nucleotide of the ith transcript exon
c-----ts(i)  - Status of the ith transcript exon [1=non-outlier; 2=outlier]
c-----b        - Number of nucleotides in the transcript exon
c-----gcba(i)  - Sum of gene conversion tracts that cover each nucleotide in the ith transcript exon
c-----agc(i)   - Number of nucleotides [i, 1=non-outlier; 2=outlier]
c-----gcx(i)   - Sum of gene conversion tracts that cover each nucleotide [i, 1=non-outlier; 2=outlier]
c-----gcx2(i)  - Sum of squares of gene conversion tracts that cover each nucleotide [i, 1=non-outlier; 2=outlier]
c-----xb(i)    - Mean of gene conversion tract coverage [i, 1=non-outlier; 2=outlier]
c-----s2(i)    - Variance of gene conversion tract coverage [i, 1=non-outlier; 2=outlier]
c      
      open(unit=1,file='GC_'//arr//'_Transcripts_List.tsv')
      ntr=0
150   ntr=ntr+1
      read(1,*,end=200) gname(ntr),ntex(ntr),ntb(ntr),nte(ntr),ts(ntr)
      goto 150
200   close(unit=1)
      ntr=ntr-1
      open(unit=1,file='GC_'//arr//'_Transcript_Tract_Depth_MLpara.csv')
      write(1,'("Transcript,Exon,Beg,End,Outlier,Mean_GC")')
      gcba=0.
      agc=0.
      gcx=0.
      gcx2=0.
      do i=1,ntr
      b=0.
      do j=ntb(i),nte(i)
      b=b+1
      agc(ts(i))=agc(ts(i))+1.
      gcba(i)=gcba(i)+chr(j)
      gcx(ts(i))=gcx(ts(i))+chr(j)
      gcx2(ts(i))=gcx2(ts(i))+(chr(j)**2.)
      end do
      gcba(i)=gcba(i)/b
      write(1,'(\a11)') gname(i)
      write(1,'(\3(",",i10))') ntex(i),ntb(i),nte(i)
      write(1,'(",",i1,",",f10.3)') ts(i),gcba(i)
      end do
      close(unit=1)
      do i=1,2
      xb(i)=gcx(i)/agc(i)
      s2(i)=(gcx2(i)-((gcx(i)**2.)/agc(i)))/(agc(i)-1)
      end do
c
c-----Random generation of tracts
c
c-----Begin random generation of gene conversion tracts
c-----nev(i,j)    - ML estimate of the number of true gene conversion tracts
c-----              in the ith arrangement and jth region
c-----amlphi(i,j) - ML estimate of mean phi, the probability of extending
c-----              the gene conversion tract to an additional nucleotide 
c-----              in the ith arrangement and jth region 
c     	
      open(unit=2,file='GC_'//arr//'_ML_Exp.tsv')
      do i=1,5
      do j=1,14
      read(2,*) amlphi(i,j),nev(i,j)
      end do
      end do
      close(unit=2)
      open(unit=2,file='GC_'//arr//'_Random_Shuffle_ML_Data.csv')
      amin=xb(1)
      amax=xb(1)
      do isim=1,10000
c      print *,arr,' Simulation =',isim
      if(mod(isim,100).eq.0) print *,arr,' Simulation =',isim
c
c-----Shuffle the gene conversion tracts and map the histogram
c-----Generate nev(i) random gene conversion tracts in the ith region
c-----using parameter amlphi(i)
c-----nl   - geometric random deviate for gene conversion tract length
c-----ia   - position within the region for the beginning of the tract.
c-----       If the tract extends past the end of the region, then
c-----       a new beginning point is chosen.
c-----iflg - keeps track of how many times a new beginning point is
c-----       sampled.  If iflg reaches 100, a new nl is chosen.  This
c-----       occurs when the random tract length exceeds the size of the
c-----       region.  
c
      chr=0
      do i=1,5
      do j=1,14
      do k=1,nev(i,j)
      iflg=0
300   nl=int(log(ran1(nsec))/log(amlphi(i,j)))
250   ia=int(ran1(nsec)*(nre(j)-nrb(j)))+nrb(j)
      if(ia+nl.gt.nre(j)) then
      iflg=iflg+1
      if(iflg.eq.100) then
      iflg=0
      goto 300
      end if
      goto 250
      endif
      do ij=ia,ia+nl+1
      chr(ij)=chr(ij)+1
      end do
      end do
      end do
      end do
c
c-----Re-estimate mean coverage 
c      
      nprob=0
      sagc=0.
      sgcx=0.
      do i=1,ntr
      if(ts(i).eq.2) then
      do j=ntb(i),nte(i)
      sagc(2)=sagc(2)+1.
      sgcx(2)=sgcx(2)+chr(j)
      end do
      end if
      end do
      sxb(2)=sgcx(2)/sagc(2)
      write(2,'(f10.3)') sxb(2)
      if(sxb(2).le.xb(2)) nprob=nprob+1
      if(sxb(2).gt.amax) then
      amax=sxb(2)
      cycle
      end if
      if(sxb(2).lt.amin) then
      amin=sxb(2)
      cycle
      end if
      end do
      close(unit=2)
      aprob=float(nprob)/10000.
c
c-----End random shuffle analysis
c-----Output results
c
      open(unit=1,file='GC_'//arr//'_Outlier_Tract_Depth_MLpara.txt')
      write(1,'("Analysis of ",a2," Gene Conversion Tracts")') arr
      write(1,'(\"Method: Use ML parameters to generate ")')
      write(1,'("Gene Conversion Tracts")')
      write(1,'(1x)')
      write(1,'(\12x,"[Mean GC Bases]")')
      write(1,'(\" [Variance GC Bases]")')
      write(1,'(" [Nucleotides]")')
      write(1,'(\"Non-Outlier ",f15.3,)') xb(1)
      write(1,'(f20.3,f14.0)') s2(1),agc(1)
      write(1,'(\"    Outlier ",f15.3,)') xb(2)
      write(1,'(f20.3,f14.0)') s2(2),agc(2)
      write(1,'(1x)')
      write(1,'("P(x<=",f7.3,")=",f10.6)') xb(2),aprob
      write(1,'("Minimum = ",f7.3,"  Maximum = ",f7.3)') amin,amax
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
