      program ngcftn
c
c-----Converts a Fasta set of aligned sequences to Nexus format
c---
c-----Written by: Stephen W. Schaeffer
c-----The Pennsylvania State University
c-----Dept. of Biology
c-----Version: 1.0
c-----Date: 19 June 2018
c
      integer nsrb(75),nsre(75)
      character*1 seq(55,20000000)
      character*3 csr(75),cytb(75),cyte(75)
      character*15 seqname,sname(55)
      nc=19787792
      nseq=0
c
c-----Input Syntenic Region Information
c-----nsr - Number of syntenic regions
c-----cytb(i) - Beginning cytogenetic section of the ith syntenic region
c-----cyte(i) - End cytogenetic section of the ith syntenic region
c-----nsrb(i) - First nucleotide of the ith syntenic region
c-----nsre(i) - Last nucleotide of the ith syntenic region
c-----nc      - Nucleotide counter
c-----nidx(i) - Index for the syntenic region of the ith nucleotide
c
c-----Other Information
c-----seq(i,j) - nucleotide of the ith sequence at the jth position
c-----seqname  - sequence name read from 'snp_Dpse_All.fof'
c-----sname(i) - sequence name of the ith strain
c
      open(unit=1,file='Syntenic_Regions.txt')
      read(1,*) nsr
      do i=1,nsr
      read(1,*) csr(i),cytb(i),cyte(i),nsrb(i),nsre(i)
      end do
      close(unit=1)
      open(unit=1,file='snp_Dpse_All.fof')
50    nseq=nseq+1
      read(1,*,end=100)  seqname
      print *,seqname
      sname(nseq)=seqname
      ip=scan(seqname,' ')
      open(unit=2,file=seqname(:ip-1)//'_Chr3.fas')
      read(2,'(1x)')
      read(2,'(19787792a1)') (seq(nseq,i),i=1,nc)
      close(unit=2)
      goto 50
100   close(unit=1)
      nseq=nseq-1
c
c-----Output Nexus File
c
      do i=1,nsr
      open(unit=1,file='Chr3_Syn_Reg_'//csr(i)//'.nex')
      open(unit=2,file='Chr3_Syn_Reg_'//csr(i)//'.fas')
      write(1,'("#nexus")')
      write(1,'("BEGIN TAXA;")')
      write(1,'("DIMENSIONS NTAX=",i2)') nseq
      write(1,'("TAXLABELS")')
      do j=1,nseq-1
      write(1,'(a15)') sname(j)
      end do
      write(1,'(a15,";")') seqname
      write(1,'("END;")')
      write(1,'(1x)')
      write(1,'("BEGIN CHARACTERS;")')
      write(1,'("DIMENSIONS NCHAR=",i8,";")') nsre(i)-nsrb(i)+1
      write(1,'("FORMAT DATATYPE=DNA MISSING=N Gap=-;")')
      write(1,'("MATRIX")')
      do j=1,nseq
c
c-----Output sequence to nexus file
c
      write(1,'(\a15,1x)') sname(j)
      write(1,'(2500000a1)') (seq(j,k),k=nsrb(i),nsre(i))
c
c-----Output sequence to fasta file
c
      write(2,'(">",a15)') sname(j)
      write(2,'(2500000a1)') (seq(j,k),k=nsrb(i),nsre(i))
      end do
      write(1,'(";")')
      write(1,'("END;")')
      write(1,'(1x)')
      write(1,'("BEGIN SETS;")')
      write(1,'("TaxSet AR = 1-15;")')
      write(1,'("TaxSet ST = 16-23;")')
      write(1,'("TaxSet PP = 24-33;")')
      write(1,'("TaxSet CH = 34-42;")')
      write(1,'("TaxSet TL = 43-51;")')
      write(1,'("TaxSet CU = 52-54;")')
      write(1,'("TaxSet Dmir = 55;")')
      write(1,'(1x)')
      write(1,'("BEGIN DNASP;")')
      write(1,'("CHROMOSOMALLOCATION= Autosome;")')
      write(1,'("GENOME= Haploid;")')
      write(1,'("END;")')
      close(unit=1)
      end do
      stop
      end
