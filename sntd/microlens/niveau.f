**********************************************************************  
**********************************************************************  
      integer function niveau(cellnum)
*
* calculates cell level for cell number  cellnum:
*
*      0+1 =     1 <=  cellnum  <=    0+   4=  4  --------> niveau 1;
*      4+1 =     5 <=  cellnum  <=    4+  16= 20  --------> niveau 2;
*     20+1 =    21 <=  cellnum  <=   20+  64= 84  --------> niveau 3;
*     84+1 =    85 <=  cellnum  <=   84+ 256=340  --------> niveau 4;
*    340+1 =   341 <=  cellnum  <=  340+1024=1364 --------> niveau 5;
*   1364+1 =  1365 <=  cellnum  <= 1364+4096=5460 --------> niveau 6;
*   5460+1 =  5461 <=  cellnum  <= .......        --------> niveau 7;
*
* input:
*       cellnum  - (actual) number of cell, whose niveau has
*               to be determined (example: cellnum = 3766; niveau => 6)
*
* output:
*       niveau   - function value   (1  <=  niveau  <= levelmax )
*
* called in:
*           main
*
* use of subroutines:
*           none
*
* use of functions:
*           none
*					October 31, 1991; J. Wambsganss
*
      implicit none
      integer i,dum,cellnum
*
      dum = cellnum
      do i = 1,200
           dum =  dum - 4**i
           if(dum.le.0)then
                niveau = i
                return
           endif
      enddo
*
      write(44,*)' cellnum,i,dum: ',cellnum,i,dum
      write(44,*)' (possible reason: levelmax = 30, but it should be '
      write(44,*)' levelmax = 15, if integer*8 cannot be specified!)'
      write(44,*)' ------------>   STOP 271 in NIVEAU !!!'
*
      stop 271
      end
