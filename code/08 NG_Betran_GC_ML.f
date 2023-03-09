      program ngbgcml
c
c-----NG Betran GC ML
c-----Maximum likelihood estimation of gene conversion parameters
c-----Reference: Betran, E., et al. (1997). "The estimation of the number
c-----           and the length distribution of gene conversion tracts
c-----           from population DNA sequence data." Genetics 146(1): 89-99.
c
c-----Written by Stephen W. Schaeffer
c-----Date: 20 March 2019
c
c-----Source: NG Betran GC ML.f
c
	integer nct(20000)
	real apsi(20000),alike(20000),bphi(200000),apmle
	character*2 arr1,arr2,reg
	character*16 fname
c-----Enter conversion tract data
c-----nc      -Number of conversion tracts
c-----nct(i)  -Length of the ith conversion tract
c-----apsi(i) -Probability of a site being informative of a conversion tract
c-----aphi    -Probability ofa converting tract elongating to an additional nucleotide
c-----ans     -Number of sequences
c-----anb     -Number of base pairs
c
      open(unit=3,file='GC.fof')
500   read(3,'(a16)',end=1000) fname
      arr1=fname(5:6)
      arr2=fname(8:9)
      reg=fname(11:12)
      print *,fname
      open(unit=1,file=fname)
	read(1,*) nc
	do i=1,nc
	read(1,*) nct(i),apsi(i)
	end do
	read(1,*) ans
	read(1,*) anb
	close(unit=1)
c
c-----Estimate mean conversion tract length
c-----el -Mean conversion tract length
c
	el=0.
	do i=1,nc
      el=el+float(nct(i))
	end do
	el=el/float(nc)
c
c-----Maximum likelihood parameters
c-----apm  -Initial value of phi
c-----ainc -Incremental increase for phi
c-----nmin -Minimum value for interative loop
c-----nsim -Number of interations to estimate phi
c
c
	apm=980000
	nmin=1
	nsim=20000
	alike=0.
c
c-----Set lower and upper limits on phi
c-----apl=initial lower phi value
c-----apu=initial upper phi value
c
	do isim=nmin,nsim
	if(mod(isim,1000).eq.0) print *,'isim=',isim
c
c-----Increment phi (bphi)
c
      bphi(isim)=(apm+float(isim))/1000000.
      plk=0.
	do i=1,nc
      call problen(apsi(i),bphi(isim),nct(i),p)
	if(p.gt.0) then
	p=log(p)
	else
	p=0.
	end if
      plk=plk+p
	end do
	alike(isim)=plk
	end do
	open(unit=1,file='GC_LH_'//arr1//'_'//arr2//'_'//reg//'.csv')
	write(1,'("Iteration,phi,lnL")')
	do i=nmin,nsim
	write(1,'(i6,100(",",f20.7))') i,bphi(i),alike(i)
	end do
	close(unit=1)
	checkamin: do i=nsim,nmin,-1
	if(alike(i).lt.0.) then
	ni=i
	exit checkamin
	end if
	end do checkamin
	nsim=ni
c
c-----Find the maximum value of phi
c-----amax  -
c-----imle  -index of the maximum phi value
c-----apmle -bphi(imle) maximum likelihood estimate of phi
c-----nl    -index of the lower bound of phi
c-----apt   -sum of the likelihood values 
c
	amax=alike(nsim)
c      
c-----Find the maximum likelihood value
c
	findmax: do i=nsim-1,nmin,-1
	if(alike(i).gt.amax) then
	amax=alike(i)
	else
	exit findmax
	endif
	end do findmax
c
c-----Find the index associated with the maximum likelihood value
c	
	findpos: do i=nsim,nmin,-1
	if(alike(i).eq.amax) then
      imle=i
	apmle=bphi(imle)
	exit findpos
	end if
	end do findpos
c
c-----Estimate a lower and upper confidence limits
c
      avar1=(alike(imle+1)-alike(imle))/(bphi(imle+1)-bphi(imle))
	avar2=(alike(imle+2)-alike(imle+1))/(bphi(imle+2)-bphi(imle+1))
	avar=abs(((avar2-avar1)/(bphi(imle+1)-bphi(imle))))
	asd=(1./avar)**0.5
      aphi=apmle
      apll=aphi-(2.*asd)
	apul=aphi+(2.*asd)
      apll=1./(1.-apll)
	apmle=1./(1.-apmle)
	apul=1./(1.-apul)
	open(unit=1,file='A_Gene_Conv_Estimates.txt',status='new',err=300)
	write(1,'(\"Don Rec Reg      Phi    Obs(N)")')
      write(1,'(\" LowLim(N) Expect(N)  UpLim(N)")')
      write(1,'("  Prob(UE)  Obs#(CT)  Exp#(CT)  Prob(TS)")')
      close(unit=1)
300	open(unit=1,file='A_Gene_Conv_Estimates.txt',access='append')
c
c-----apl2 - probability of an undetected gene conversion event
c-----ek   - expected true number of conversion tracts
c-----apt  - probability of a transferred site
c
	apm=1.-(1./apmle)
c	print *,apm,apsi(1)
c	pause
	apl2=(1.-apm)*(1.-apm+(2.*apm*apsi(1))-(apm*(apsi(1)**2.)))
      apl2=apl2/((1.-apm+(apm*apsi(1)))**2.)
	apl2=1.-apl2
	ek=float(nc)/apl2
      apt=(ek*apmle)/(ans*anb)
	write(1,'(1x,a2,2x,a2,2x,a2,f9.6,4(f10.2),f10.3,I10,f10.2,f10.3)')
     +arr1,arr2,reg,aphi,el,apll,apmle,apul,1.-apl2,nc,ek,apt
	close(unit=1)
	goto 500
1000	close(unit=3)
	stop
	end
c
c-----Subroutine problen
c-----Estimates P(L*=l) equation (6) from Betran et al. (1997)
c-----This equation is iterated until the difference between current (eln)
c-----and prior (el1) estimate is less than 0.0000000001   
c 
c-----aps -Estimate of psi
c-----aph -Tested value of phi
c-----nl  -Track length
c-----pr  -Probability
c-----el1 -Prior   estimate of P(L*=l)
c-----eln -Current estimate of P(L*=l)
c
	subroutine problen(aps,aph,nl,pr)
	integer nl
      real aps,aph,el1,eln,plnn,plnd,pln,pn,pr
	el1=0.
	eln=0.
	do il=nl,100000
c
c-----Estimates P(L*|N=n) equation (4) from Betran et al. (1997)
c-----an   -True     tract length
c-----al   -Observed tract length
c-----plnn -P(L*=l|N=n) numerator   equation (4)
c-----plnd -P(L*=l|N=n) denominator equation (4)
c-----pln  -P(L*=l|N=n) equation (4)
c-----pn   -P(N=n) 
c	
	an=float(il)
	al=float(nl)
	plnn=((an-al+1)*(aps**2)*((1-aps)**(an-al)))
	plnd=1-((1-aps)**an)-(an*((1-aps)**(an-1))*aps)
	pln=plnn/plnd
	pn=(aph**(an-2.))*(1.-aph)
      eln=eln+(pln*pn)
	if(abs(eln-el1).lt.0.0000000001) then
	pr=eln
	return
	else
	el1=eln
	end if
	end do
	end subroutine problen
