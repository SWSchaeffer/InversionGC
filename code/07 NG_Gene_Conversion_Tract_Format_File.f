      program ncgctff
c
c-----NG Gene Conversion Tract Format File
c-----Setup files for NG Betran GC ML 
c
c-----Written by Stephen W. Schaeffer
c-----Date: 12 October 2020
c
c-----Source: NG_Gene_Conversion_Tract_Format_File.f
c
      integer lng(20000),ni(11),ncb(14),nce(14),nb(14),reg(20000)
      real psi(20000)
      character*1 t,exc
      character*2 ar1(20000),ar2(20000),regn(14)
      character*8 ori,qry
c
c-----Initiate variables
c-----t      - tab character
c-----ni(i)  - number of the ith inversion
c-----ncb(i) - beginning nucleotide coordinate for the ith region
c-----nce(i) - end nucleotide coordinate for the ith region
c-----nb(i)  - number of nucleotides in the ith syntenic block
c      
      t=char(9)
      ni(1)=15
      ni(2)=9
      ni(3)=3
      ni(4)=10
      ni(5)=8
      ni(6)=9
      regn(1)='01'
      regn(2)='02'
      regn(3)='03'
      regn(4)='04'
      regn(5)='05'
      regn(6)='06'
      regn(7)='07'
      regn(8)='08'
      regn(9)='09'
      regn(10)='10'
      regn(11)='11'
      regn(12)='12'
      regn(13)='13'
      regn(14)='14'
c
c-----Input Minimum and Maximum nucleotide coordinates
c-----for the 14 Arrowhead syntenic blocks
c      
      open(unit=1,file='AR_MinMax.tsv')
      read(1,'(1x)')
      do i=1,14
      read(1,*) ncb(i),nce(i)
      nb(i)=nce(i)-ncb(i)+1
      end do
      close(unit=1)
c
c-----Open the gene conversion data file
c-----nr     - counter for the number of gene conversion tracts in the ar1(i)_ar2(i)_reg(i) data set
c-----ar1(i) - name of the first  arrangement for the ith data point
c-----ar2(i) - name of the second arrangement for the ith data point
c-----reg(i) - number of the syntenic block   for the ith data point
c-----lng(i) - length of the gene conversion tract for the ith data point
c-----psi(i) - psi value for the gene conversion tract for the ith data point
c-----ori    - concatenated ar1_ar2_reg for the original data set
c-----ia1    - converts the ar1 inversion to a number [AR=1,CH=3,CU=5,PP=7,ST=9,TL=10]
c-----ia2    - converts the ar2 inversion to a number [AR=1,CH=3,CU=5,PP=7,ST=9,TL=10]
c-----ir     - converts the syntenic block character to an integer [1-14]
c-----qry    - concatenated ar1_ar2_reg for the new data point
c      
      open(unit=4,file='GC_Tracts_Summary.tsv')
      write(4,'("Comparison",a1,"Tracts",a1,"Bases",a1,"n")') t,t,t
      open(unit=3,file='GC.fof')
      open(unit=1,file='GC_List.txt')
      read(1,'(1x)')
      nr=0
50    nr=nr+1
      if(nr.eq.1) then
      read(1,*,end=100) ar1(nr),ar2(nr),reg(nr),lng(nr),psi(nr),exc
      ori=ar1(nr)//'_'//ar2(nr)//'_'//regn(reg(nr))
      if(ar1(nr).eq.'AR') ia1=1
      if(ar1(nr).eq.'CH') ia1=2
      if(ar1(nr).eq.'CU') ia1=3
      if(ar1(nr).eq.'PP') ia1=4
      if(ar1(nr).eq.'ST') ia1=5
      if(ar1(nr).eq.'TL') ia1=6
      if(ar2(nr).eq.'AR') ia2=1
      if(ar2(nr).eq.'CH') ia2=2
      if(ar2(nr).eq.'CU') ia2=3
      if(ar2(nr).eq.'PP') ia2=4
      if(ar2(nr).eq.'ST') ia2=5
      if(ar2(nr).eq.'TL') ia2=6
      else
      read(1,*) ar1(nr),ar2(nr),reg(nr),lng(nr),psi(nr),exc
      qry=ar1(nr)//'_'//ar2(nr)//'_'//regn(reg(nr))
c
c-----Check the new data point against the original data point
c      
      if(qry.ne.ori) then
c
c-----If different, output the data
c
      write(3,'("GCT_",a8,".txt")') ori
      write(4,'(\"GCT_",a8,a1,i6,a1,i8,a1)') ori,t,nr-1,t,nb(reg(1)),t
      write(4,'(i3)') ni(ia1)+ni(ia2)
      open(unit=2,file='GCT_'//ori//'.txt')
      write(2,'(i4)') nr-1
      do i=1,nr-1
      write(2,'(i6,a1,f7.5)') lng(i),t,psi(i)
      end do
      write(2,'(i10)') nb(reg(1))
      write(2,'(i5,a1,f7.5)') ni(ia1)+ni(ia2)
      close(unit=2)
      if(exc.eq.'Y') goto 100
c
c-----Move the new data point into the first position of a new data set
c-----Reset ori, ia1, ia2, and ir 
c      
      ar1(1)=ar1(nr)
      ar2(1)=ar2(nr)
      reg(1)=reg(nr)
      lng(1)=lng(nr)
      psi(1)=psi(nr)
      nr=1
      ori=ar1(nr)//'_'//ar2(nr)//'_'//regn(reg(nr))
      if(ar1(nr).eq.'AR') ia1=1
      if(ar1(nr).eq.'CH') ia1=2
      if(ar1(nr).eq.'CU') ia1=3
      if(ar1(nr).eq.'PP') ia1=4
      if(ar1(nr).eq.'ST') ia1=5
      if(ar1(nr).eq.'TL') ia1=6
      if(ar2(nr).eq.'AR') ia2=1
      if(ar2(nr).eq.'CH') ia2=2
      if(ar2(nr).eq.'CU') ia2=3
      if(ar2(nr).eq.'PP') ia2=4
      if(ar2(nr).eq.'ST') ia2=5
      if(ar2(nr).eq.'TL') ia2=6
      end if
      end if
      goto 50     
100   close(unit=1)
      close(unit=3)
      close(unit=4)
      stop
      end
