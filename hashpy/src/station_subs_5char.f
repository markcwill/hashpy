
c GETSTAT_SCSN finds station locations for SCSN stations
c
c   inputs:
c     stlfile - file with the stations and locations, in alphabetical order!
c     snam - name of the station of interest, 4 characters
c     scom - station component, 3 characters
c     snet - network, 2 characters
c   outputs:
c     flat,flon,felev - station lat, lon, and elevation
c
c   input file "stlfile" format = SCEDC format
c     columns  format   value
c     -------------------------
c     1-2        a2     network code
c     5-9        a5     station name
c     11-13      a3     station component
c     61-69      f9.5   station latitude (degrees, signed)
c     71-80      f10.5  station longitude (degrees, signed)
c     82-86      i5     station elevation (meters)
c
c
      subroutine GETSTAT_SCSN(stlfile,snam,scom,snet,flat,
     &                       flon,felev)
      parameter(nsta0=20000)
      character stlfile*100
      character*5 snam,stname(nsta0)
      character*3 scom,scompt(nsta0),scom2
      character*2 snet,snetwk(nsta0)
      real slat(nsta0),slon(nsta0),selev(nsta0)
      logical firstcall
      save firstcall,stname,slat,slon,selev,nsta,scompt,snetwk
      data firstcall/.true./
      
c read in station list - in alphabetical order!
      if (firstcall) then
         firstcall=.false.
         open (19,file=stlfile)
         do i=1,nsta0
           read (19,11,end=12) snetwk(i),stname(i),scompt(i),
     &                      slat(i),slon(i),ntemp
           selev(i)=real(ntemp)/1000.0
         end do
11       format (a2,2x,a5,1x,a3,37x,f9.5,1x,f10.5,1x,i5)
12       nsta=i-1
         close (19)
      end if  
       
      scom2=scom                             ! short-period stations are
      if (scom(1:1).eq."V") scom2(1:1)="E"   ! interchangably called V and E   
      if (scom(1:1).eq."E") scom2(1:1)="V"           
       
c binary search for station name
      i1=1
      i2=nsta
      do 30 it=1,30
         i=(i1+i2)/2
         if (snam.eq.stname(i)) then
           goto 40
         else if (i1.eq.i2) then
           goto 999
         else if (snam.lt.stname(i)) then
            i2=i-1
         else 
            i1=i+1
         end if
30    continue
      goto 999
      
c search for proper component/network
40    i1=i
45    continue
        if (i1.gt.nsta) goto 50
        if (snam.ne.stname(i1)) goto 50
        if ((scom(1:2).eq.scompt(i1)(1:2)).and.
     &             snet.eq.snetwk(i1)) goto 900
        if ((scom2(1:2).eq.scompt(i1)(1:2)).and.
     &             snet.eq.snetwk(i1)) goto 900
        if ((scom(1:2).eq.'XX').and.
     &                 (snet.eq.'XX')) goto 900
        i1=i1+1
      goto 45
50    i1=i-1
55    continue
        if (i1.lt.1) goto 999
        if (snam.ne.stname(i1)) goto 999
        if ((scom(1:2).eq.scompt(i1)(1:2)).and.
     &             snet.eq.snetwk(i1)) goto 900
        if ((scom2(1:2).eq.scompt(i1)(1:2)).and.
     &             snet.eq.snetwk(i1)) goto 900
        if ((scom(1:2).eq.'XX').and.
     &                 (snet.eq.'XX')) goto 900
        i1=i1-1
      goto 55
      
900   flat=slat(i1)
      flon=slon(i1)
      felev=selev(i1)
      return
999   print *,'***station not found ',snam,' ',scom,' ',snet
      flat=999.
      flon=999.
      felev=999.
      return
      end

c --------------------------------------------------------------- c


c CHECK_POL determines whether the polarity of a given station
c     was reversed at a given time.
c
c   inputs:
c     polfile - file with the stations and times of reversal;
c               must be in alphabetical order; if a station is
c               reversed during multiple time periods, it is 
c               listed on multiple lines (same as FPFIT)
c     snam - name of the station of interest, 4 characters
c     evyr  - year of interest, 4 digits
c     evmon - month of interest
c     evdy  - day of interest
c     evhr  - hour of interest (not implemented)
c   output:
c     stpol - station polarity: 1=correct, -1=reversed 
c
c   input file "polfile" format:
c     columns  format   value
c     -------------------------
c     1-5        a5     station name
c     6-9        i4     begining of reversal: year
c     10-11      i2                           month
c     12-13      i2                           day
c     15-18      i4     end of reversal: year
c     19-20      i2                      month
c     21-22      i2                      day

      subroutine CHECK_POL(polfile,snam,evyr,evmon,evdy,evhr,stpol)
      parameter(nsta0=300)
      character polfile*100,polfileold*100
      character*5 snam,statname(nsta0)
      integer begtime(nsta0),endtime(nsta0)
      integer evyr,evmon,evdy,evhr,stpol,nstat(nsta0),nth(nsta0)
      integer i,i1,i2,count,itemp,nrev,evtime
      save polfileold,statname,begtime,endtime,nrev,nstat,nth

c read in polarity reversal file - in alphabetical order
      if (polfile.ne.polfileold) then
         print *,'Reading polarity reversal file ',polfile
         open (19,file=polfile)
         do i=1,300
           read (19,11,end=12) statname(i),begtime(i),endtime(i)  
           if (endtime(i).eq.0) then
             endtime(i)=99999999
           end if
         end do
11       format (a5,i8,1x,i8)
12       nrev=i-1
         close (19)
         polfileold=polfile
         nth(1)=1
         nstat(1)=1
         count=1
         do 17 i=2,nrev
           if (statname(i).eq.statname(i-1)) then
             count=count+1
           else
             count=1
           end if
           nth(i)=count
           do itemp=i-count+1,i
             nstat(itemp)=count
           end do
17       continue
      end if
      
      evtime=evyr*10000+evmon*100+evdy
      stpol=1

c binary search for a station
      i1=1
      i2=nrev
20    continue
         i=(i1+i2)/2 
         if (snam.eq.statname(i)) then
            do 25 itemp=i-nth(i)+1,i+nstat(i)-nth(i)
              if ((evtime.ge.begtime(itemp)).and.
     &                (evtime.le.endtime(itemp))) then
                stpol=-1
                goto 30
              end if
25          continue
            goto 30
         else if (i1.ge.i2) then
            goto 30
         else if (snam.lt.statname(i)) then
            i2=i-1
         else
            i1=i+1
         end if
      goto 20
30    continue
      return
      end

c ------------------------------------------------------------------- c
      
c GET_COR reads a file of station amplitude corrections
c
c   inputs:
c     stlfile - file with the stations and locations
c     snam - name of the station of interest, 4 characters
c     scom - station component, 3 characters
c     snet - network, 2 characters
c   outputs:
c     qcor - corrections to be subtracted from log(S/P)
c
c   input file format:
c     columns  format   variable
c     -------------------------
c     1-5        a5     station name
c     7-9        a3     station component (vertical - Z optional)
c     11-12      a2     network code
c     14-20      f7.4   correction to be subtracted from log(S/P)
c
c
      subroutine GET_COR(stlfile,snam,scom,snet,qcor)
      parameter(nsta0=10000)
      character stlfile*100
      character*5 snam,stname(nsta0)
      character*3 scom,scom2,scompt(nsta0)
      character*2 snet,snetwk(nsta0)
      real corr_val(nsta0)
      logical firstcall
      save firstcall,stname,corr_val,nsta,scompt,snetwk
      data firstcall/.true./
      
c read in station list - in alphabetical order!
      if (firstcall) then
         firstcall=.false.
         open (19,file=stlfile)
         do 10 i=1,nsta0
           read (19,11,end=12) stname(i),scompt(i),snetwk(i),
     &                              corr_val(i)
10       continue
11       format (a5,1x,a3,a2,1x,f7.4)
12       nsta=i-1
         close (19)
      end if  
      
      scom2=scom                             ! short-period stations are
      if (scom(1:1).eq."V") scom2(1:1)="E"   ! called both V and E     
      if (scom(1:1).eq."E") scom2(1:1)="V"           

c binary search for station name
      i1=1
      i2=nsta
      do 30 it=1,30
         i=(i1+i2)/2
         if (snam.eq.stname(i)) then
           goto 40
         else if (i1.eq.i2) then
           goto 999
         else if (snam.lt.stname(i)) then
            i2=i-1
         else 
            i1=i+1
         end if
30    continue
      print *,'station not found'
      goto 999
      
c search for proper component/network
40    i1=i
45    continue
        if (i1.gt.nsta) goto 50
        if (snam.ne.stname(i1)) goto 50
        if (scom(1:2).eq.scompt(i1)(1:2)) goto 900
        if (scom2(1:2).eq.scompt(i1)(1:2)) goto 900
        i1=i1+1
      goto 45
50    i1=i-1
55    continue
        if (i1.lt.1) goto 999
        if (snam.ne.stname(i1)) goto 999
        if (scom(1:2).eq.scompt(i1)(1:2)) goto 900
        if (scom2(1:2).eq.scompt(i1)(1:2)) goto 900
        i1=i1-1
      goto 55

900   qcor=corr_val(i1)
      return
999   print *,'GET_COR ***station not found ',snam,' ',scom,' ',snet,
     &  ' in file ',stlfile
      qcor=-999.
      return
      end

