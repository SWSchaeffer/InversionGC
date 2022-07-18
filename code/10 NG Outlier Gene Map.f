      program ngogm
c
c-----NG Outlier Gene Map
c-----Maps outlier genes to a chromosomal map shown as a box with
c-----dots.  The map goes from the top (5') to bottom (3') and the
c-----outliers represented as dots on the map.
c
c-----Written by Stephen W. Schaeffer
c-----Date: 8 March 2019
c
c-----Source: NG Outlier Gene Map.f
c
      integer ntr, ntm(200),ntb(200),nte(200)
      character*2 arr
      character*11 gname(200)
      write(*,'("Input the gene arrangement.")')
      read(*,*) arr
c
c-----Input outlier genes
c
c-----ntr      - Number of outlier genes
c-----gname(i) - Name of the ith transcript
c-----ntm(i)   - Midpoint coordinate of the ith transcript
c-----ntb(i)   - Beginning coordinate of the ith transcript
c-----nte(i)   - End coordinate of the ith transcript
c      
      open(unit=1,file='GC_'//arr//'_Outlier_List.txt')
      ntr=0
150   ntr=ntr+1
      read(1,*,end=200) gname(ntr),ntm(ntr),ntb(ntr),nte(ntr)
      goto 150
200   close(unit=1)
      ntr=ntr-1
      open(unit=1,file='GC_'//arr//'_Outlier_Plot.ps')
      write(1,'("%!PS-Adobe-3.0")')
c
c-----Add the horizontal scale bar
c      
      write(1,'("newpath")')
      write(1,'("0 0 moveto")')
      write(1,'("0 422.731 rlineto")')
      write(1,'("22.776 0 rlineto")')
      write(1,'("0 -422.731 rlineto")')
      write(1,'("-22.776 0 rlineto")')
      write(1,'("closepath")') 
      write(1,'("stroke")')
c
c-----Add the outliers to the plot
c      
      write(1,'("1.0 1.0 1.0 1.0 setcmykcolor")')
      do i=1,ntr
      y=422.731-((float(ntm(i))/19797792.)*422.731)
      write(1,'("newpath")')
      write(1,'("11 ",f8.3," 2 0 360 arc")') y
      write(1,'("fill")') 
      end do
      stop
      end
