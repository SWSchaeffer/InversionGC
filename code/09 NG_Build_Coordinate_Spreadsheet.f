      program ngbcs
c
c-----NG Build Coordinate Spreadsheet
c-----Input a coordinate file and reorganize into a new spreadsheet
c
c-----Written by Stephen W. Schaeffer
c-----Date: 1 February 2019
c
c-----Source: NG Build Coordinate Spreadsheet.f
c
      character*1 c
      character*11 gname
      integer ori,nexon,nbeg,nend
c
c-----gname - name of the gene
c-----ori   - gene orientation
c-----nbeg  - first nucleotide of the exon
c-----nend  - last  nucleotide of the exon
c      
      c=','
      open(unit=1,file='Coordinates.txt')
      open(unit=2,file='Coordinates_List.csv')
      write(2,'("Transcript_ID,Ori,Exon,Pos_Beg,Pos_End")')
50    read(1,*,end=100) gname,ori
      read(1,*) nexon
      do i=1,nexon
      read(1,*) nbeg,nend
      write(2,'(a11,4(a1,i10))') gname,c,ori,c,i,c,nbeg,c,nend
      end do
      goto 50
100   continue 
      stop
      end
