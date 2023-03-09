      program ngcsm
c
c-----NG Convert SNP Table to Fasta from Reference
c-----The program takes the SNP table “Chr3_FB_SnpTable.tsv” 
c-----and the reference sequence for Chr3 and generates individual
c-----sequences for each strain.
c
c-----Stephen W. Schaeffer
c-----version: 1.0
c-----18 September 2020
c
      parameter(ns=55)
      integer nb,npos(1500000)
      character*1 rseq(20000000),seq(ns,1500000),oseq(20000000)
      character*14 spos,sname(ns)
c
c-----Enter the polymorphism data
c-----nb       - number of variable sites
c-----spos     - column header
c-----sname(i) - name of the ith strain
c-----npos(i)  - positon of the ith SNP
c-----rseq     - reference sequence
c-----oseq     - output sequence with SNP differences incorporated
c
      nb=0
      open(unit=1,file='Chr3_FB_SNPTable.tsv')
      read(1,*) spos,(sname(i),i=1,ns)
200   nb=nb+1
	 read(1,*,end=100) npos(nb),(seq(i,nb),i=1,ns)
      goto 200
100	 nb=nb-1
      close(unit=1)
      open(unit=1,file='AR_REF_Chr3.fas')
      read(1,'(1x)')
      read(1,'(19797792a1)') (rseq(i),i=1,19797792)
      close(unit=1)
      do i=2,ns
      print *,sname(i)
      ic=scan(sname(i),' ')
      open(unit=1,file=sname(i)(:ic-1)//'_Chr3.fas')
      if(ic-1.eq.6) then
      write(1,'(">",a6)') sname(i)(:ic-1)
      elseif(ic-1.eq.7) then
      write(1,'(">",a7)') sname(i)(:ic-1)
      elseif(ic-1.eq.8) then
      write(1,'(">",a8)') sname(i)(:ic-1)
      elseif(ic-1.eq.9) then
      write(1,'(">",a9)') sname(i)(:ic-1)
      elseif(ic-1.eq.10) then
      write(1,'(">",a10)') sname(i)(:ic-1)
      elseif(ic-1.eq.11) then
      write(1,'(">",a11)') sname(i)(:ic-1)
      elseif(ic-1.eq.12) then
      write(1,'(">",a12)') sname(i)(:ic-1)
      elseif(ic-1.eq.13) then
      write(1,'(">",a13)') sname(i)(:ic-1)
      elseif(ic-1.eq.14) then
      write(1,'(">",a14)') sname(i)(:ic-1)
      endif
      oseq=rseq
      do j=1,nb
      oseq(npos(j))=seq(i,j)
      end do
      write(1,'(19787792a1)') (oseq(j),j=1,19797792)
      close(unit=1)
      end do
      stop
      end
