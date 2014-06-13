
! GETSTAT_TRI finds station locations for TriNet stations
!
!   inputs:
!     stlfile - file with the stations and locations, in alphabetical order!
!     snam - name of the station of interest, 4 characters
!     scom - station component, 3 characters
!     snet - network, 2 characters
!   outputs:
!     flat,flon,felev - station lat, lon, and elevation
!
!   input file "stlfile" format:
!     columns  format   value
!     -------------------------
!     1-4        a4     station name
!     6-8        a3     station component
!     42-50      f9.5   station latitude (degrees, signed)
!     52-61      f10.5  station longitude (degrees, signed)
!     63-67      i5     station elevation (meters)
!     91-92      a2     network code
!
!
subroutine GETSTAT_TRI(stlfile, snam, scom, snet, flat, flon, felev)
!f2py intent(in)  stlfile
!f2py intent(in)  snam
!f2py intent(in)  scom 
!f2py intent(in)  snet 
!f2py intent(out) flat  
!f2py intent(out) flon  
!f2py intent(out) felev  
    parameter(nsta0=20000)

    character(len=100) :: stlfile
    character(len=4) :: snam, stname(nsta0)
    character(len=3) :: scom, scompt(nsta0), scom2
    character(len=2) :: snet, snetwk(nsta0)
    character(len=*), parameter :: fmtsta = '(a4,1x,a3,33x,f9.5,1x,f10.5,1x,i5,23x,a2)'
    real :: slat(nsta0), slon(nsta0), selev(nsta0)
    logical :: firstcall
    save firstcall, stname, slat, slon, selev, nsta, scompt, snetwk
    data firstcall /.true./


    ! read in station list - in alphabetical order!
    if (firstcall) then
        firstcall = .false.
        open (19,file=stlfile)
        do i=1, nsta0
            read (19,fmtsta,end=12) stname(i), scompt(i), slat(i), slon(i), ntemp, snetwk(i)
            selev(i) = real(ntemp)/1000.0
        end do
!11      format (a4,1x,a3,33x,f9.5,1x,f10.5,1x,i5,23x,a2)
12      nsta = i-1 !12
        close(19)
    end if  
       
    scom2 = scom                             ! short-period stations are
    if (scom(1:1) == "V") scom2(1:1) = "E"   ! interchangably called V and E   
    if (scom(1:1) == "E") scom2(1:1) = "V"           

    ! binary search for station name
    i1 = 1
    i2 = nsta
    d30: do it=1, 30
        i = (i1+i2)/2
        if (snam == stname(i)) then
            goto 40
        else if (i1 == i2) then
            goto 999
        else if (snam < stname(i)) then
            i2 = i-1
        else 
            i1 = i+1
        end if
    end do d30
    go to 999
      
    ! search for proper component/network
40  i1 = i
45  continue
    if (i1 > nsta) go to 50
    if (snam /= stname(i1)) go to 50
    if ((scom(1:2) == scompt(i1)(1:2)) .and. snet == snetwk(i1)) go to 900
    if ((scom2(1:2) == scompt(i1)(1:2)) .and. snet == snetwk(i1)) go to 900
    if ((scom(1:2) == 'XX') .and. (snet == 'XX')) go to 900
    i1 = i1+1
    go to 45
50  i1 = i-1
55  continue
    if (i1 < 1) go to 999
    if (snam /= stname(i1)) go to 999
    if ((scom(1:2) == scompt(i1)(1:2)) .and. snet == snetwk(i1)) go to 900
    if ((scom2(1:2) == scompt(i1)(1:2)) .and. snet == snetwk(i1)) go to 900
    if ((scom(1:2) == 'XX') .and. (snet == 'XX')) go to 900
    i1 = i1-1
    go to 55
      
900 flat = slat(i1)
    flon = slon(i1)
    felev = selev(i1)
    return
999 print *,'***station not found ',snam,' ',scom,' ',snet
    flat = 999.
    flon = 999.
    felev = 999.
    return
end subroutine

! --------------------------------------------------------------- !


! CHECK_POL determines whether the polarity of a given station
!     was reversed at a given time.
!
!   inputs:
!     polfile - file with the stations and times of reversal;
!               must be in alphabetical order; if a station is
!               reversed during multiple time periods, it is 
!               listed on multiple lines (same as FPFIT)
!     snam - name of the station of interest, 4 characters
!     evyr  - year of interest, 4 digits
!     evmon - month of interest
!     evdy  - day of interest
!     evhr  - hour of interest (not implemented)
!   output:
!     stpol - station polarity: 1=correct, -1=reversed 
!
!   input file "polfile" format:
!     columns  format   value
!     -------------------------
!     1-4        a4     station name
!     6-9        i4     begining of reversal: year
!     10-11      i2                           month
!     12-13      i2                           day
!     15-18      i4     end of reversal: year
!     19-20      i2                      month
!     21-22      i2                      day

subroutine CHECK_POL(polfile,snam,evyr,evmon,evdy,evhr,stpol)
!f2py intent(in) polfile
!f2py intent(in) snam
!f2py intent(in) evyr  
!f2py intent(in) evmon  
!f2py intent(in) evday 
!f2py intent(in) evhr 
!f2py intent(out) stpol   
    parameter(nsta0=300)
    character(len=100) :: polfile, polfileold
    character(len=4) ::  snam, statname(nsta0)
    character(len=*), parameter :: fmtpol = '(a4,1x,i8,1x,i8)'
    integer, dimension(nsta0) :: begtime, endtime, nstat, nth
    integer :: evyr, evmon, evdy, evhr, stpol
    integer :: i, i1, i2, count, itemp, nrev, evtime
    save polfileold, statname, begtime, endtime, nrev, nstat, nth


    ! read in polarity reversal file - in alphabetical order
    if (polfile /= polfileold) then
        print *,'Reading polarity reversal file ',polfile
        open (19,file=polfile)
        do i=1, 300
            read (19,fmtpol,end=12) statname(i),begtime(i),endtime(i)  
            if (endtime(i) == 0) then
                endtime(i) = 99999999
            end if
        end do
12      nrevi = i-1 !12 ^--11 format etc
        close (19)
        polfileold = polfile
        nth(1) = 1
        nstat(1) = 1
        count = 1
        d17: do i=2, nrev
            if (statname(i) == statname(i-1)) then
                count = count+1
            else
                count = 1
            end if
            nth(i) = count
            do itemp=i-count+1, i
                nstat(itemp) = count
            end do
        end do d17
    end if
      
    evtime =evyr*10000 + evmon*100 + evdy
    stpol = 1

    ! binary search for a station
    i1 = 1
    i2 = nrev
20  continue
    i = (i1 + i2)/2 
    if (snam == statname(i)) then
        d25: do itemp=i-nth(i)+1, i+nstat(i)-nth(i)
            if ((evtime >= begtime(itemp)) .and. (evtime <= endtime(itemp))) then
                stpol = -1
                go to 30
            end if
        end do d25
        go to 30
    else if (i1 >= i2) then
        goto 30
    else if (snam < statname(i)) then
        i2 = i-1
    else
        i1 = i+1
    end if
    goto 20
30  continue
    return
end subroutine

! ------------------------------------------------------------------- !
      
! GET_COR reads a file of station amplitude corrections
!
!   inputs:
!     stlfile - file with the stations and locations
!     snam - name of the station of interest, 4 characters
!     scom - station component, 3 characters
!     snet - network, 2 characters
!   outputs:
!     qcor - corrections to be subtracted from log(S/P)
!
!   input file format:
!     columns  format   variable
!     -------------------------
!     1-4        a4     station name
!     7-9        a3     station component (vertical - Z optional)
!     11-12      a2     network code
!     14-20      f7.4   correction to be subtracted from log(S/P)
!
!
subroutine GET_COR(stlfile,snam,scom,snet,qcor)
!f2py intent(in) stlfile  
!f2py intent(in) snam
!f2py intent(in) scom  
!f2py intent(in) snet  
!f2py intent(out) qcor       
    parameter(nsta0=10000)
    character(len=100) ::  stlfile
    character(len=4) :: snam,stname(nsta0)
    character(len=3) :: scom,scom2,scompt(nsta0)
    character(len=2) :: snet,snetwk(nsta0)
    character(len=*), parameter :: fmtcor = '(a4,2x,a3,a2,1x,f7.4)'
    real :: corr_val(nsta0)
    logical :: firstcall
    save firstcall, stname, corr_val, nsta, scompt, snetwk
    data firstcall /.true./
      
! read in station list - in alphabetical order!
    if (firstcall) then
        firstcall=.false.
        open (19,file=stlfile)
        do i=1, nsta0
            read (19,fmtcor,end=12) stname(i), scompt(i), snetwk(i), corr_val(i)
        end do
12      nsta = i-1
        close(19)
    end if  
      
    scom2 = scom                             ! short-period stations are
    if (scom(1:1) == "V") scom2(1:1) = "E"   ! called both V and E     
    if (scom(1:1) == "E") scom2(1:1) = "V"           

    ! binary search for station name
    i1 = 1
    i2 = nsta
    d30: do it=1, 30
        i = (i1+i2)/2
        if (snam == stname(i)) then
            go to 40
        else if (i1 == i2) then
            go to 999
        else if (snam < stname(i)) then
            i2 = i-1
        else 
            i1 = i+1
        end if
    end do d30 !   continue
    print *,'station not found'
    go to 999
      
    ! search for proper component/network
40  i1 = i
45  continue
    if (i1 > nsta) go to 50
    if (snam /= stname(i1)) go to 50
    if (scom(1:2) == scompt(i1)(1:2)) go to 900
    if (scom2(1:2) == scompt(i1)(1:2)) go to 900
    i1=i1+1
    go to 45
50  i1 = i-1
55  continue
    if (i1 < 1) go to 999
    if (snam /= stname(i1)) go to 999
    if (scom(1:2) == scompt(i1)(1:2)) go to 900
    if (scom2(1:2) == scompt(i1)(1:2)) go to 900
    i1 = i1-1
    go to 55

900 qcor = corr_val(i1)
    return
999 print *,'GET_COR ***station not found ',snam,' ',scom,' ',snet,' in file ',stlfile
    qcor = -999.
    return
end subroutine

