! MK_TABLES creates tables of takeoff angles given 1D velocity models.

!   output:
!     ntab - number of tables (max nindex)
!
!   you are prompted for the names of the 1D velocity model files,
!   velocity file format (free format):
!     depth(km) P_velocity(km/s)
      
subroutine MK_TABLE(ntab)
!f2py intent(in,out) ntab

    include 'vel.inc'
    real :: table(nx0,nd0,nindex), delttab(nx0), deptab(nd0)
    integer :: ndel, ndep

!   common block:
!       table(nx0,nd0,nindex)  =  takeoff angle table
!           delttab(nx0)  =  list of ranges for tables
!           deptab(nd0)  =  list of source depths for tables
!           ndel       =  number of distance points in table
!           ndep       =  number of source depths in table
    common /angtable/ table, delttab, deptab, ndel, ndep

    integer :: itab
    integer, parameter :: nray0 = 10001
    real, parameter :: degrad = 180./3.14159265
    real, dimension(1000) :: z , alpha, slow
    real, dimension(20000) ::  xsave, tsave, psave, usave
    real, dimension(nray0) :: deltab, tttab, ptab
    real, dimension(nray0, nd0) :: depxcor, depucor, deptcor, tt
    character(len=100) :: vmodel

    if (ntab /= 1) then
        print *,'Enter number of velocity models (max ',nindex,')'
        read *, ntab
    end if


d300: do itab=1, ntab
    print *,'Enter file name for velocity model ', itab
    read (*,'(a)') vmodel

! set up table
    qtempdep2 = dep2 + dep3/20.
    ndep = int((qtempdep2-dep1)/dep3) + 1
    do idep=1, ndep
        dep = dep1 + dep3*real(idep-1)
        deptab(idep) = dep
    end do

! read velocity model TODO: fix this Mark
    open (7, file=vmodel, status='old')
    do i=1, 1000
        read (7, *, end=30) z(i), alpha(i)
    end do
    print *,'***1000 point maximum exceeded in model'
30  close (7)

! Start
    z(i) = z(i-1)           
    alpha(i) = alpha(i-1)
    npts = i
    npts_old = npts
    do i=npts_old, 2, -1
        do idep=ndep, 1, -1
            if ((z(i-1) <= (deptab(idep)-0.1)) .and. (z(i) >= (deptab(idep)+0.1))) then
                npts = npts+1
                do j=npts, i+1, -1
                    z(j) = z(j-1)
                    alpha(j) = alpha(j-1)
                end do
                z(i) = deptab(idep)
                frac = (z(i)-z(i-1))/(z(i+1)-z(i-1))
                alpha(i) = alpha(i-1) + frac*(alpha(i+1) - alpha(i-1))
            end if
        end do
    end do
    do i=1, npts
        slow(i) = 1./alpha(i)
    end do
    pmax = slow(1)
    plongcut = slow(npts)
    pstep = (pmax-pmin)/float(nump)


! do P-wave ray tracing
    npmax = int((pmax+pstep/2.-pmin)/pstep) + 1
    d200: do np=1, npmax
        p = pmin + pstep*real(np-1)
        ptab(np) = p
        x = 0.
        t = 0.
        imth = 3
        do idep=1, ndep
            if (deptab(idep) == 0.) then
                depxcor(np,idep) = 0.
                deptcor(np,idep) = 0.
                depucor(np,idep) = slow(1)
            else
                depxcor(np,idep) = -999.
                deptcor(np,idep) = -999.
                depucor(np,idep) = -999.
            end if
        end do
        d100: do i=1, npts-1
            if (z(i) >= 9999) then 
                deltab(np) = -999.
                tttab(np) = -999.
                cycle d200
            end if
            h = z(i+1) - z(i)
            if (h==0.) cycle d100 ! skip if interface
            call LAYERTRACE(p, h, slow(i), slow(i+1), imth, dx, dt, irtr)
            x = x+dx
            t = t+dt
            if (irtr==0 .or. irtr==2) exit d100 ! ray has turned
            xdeg = x                       ! actually in km 
            tmin = t                       ! actually in s
            do idep=1, ndep
                if (abs(z(i+1) - deptab(idep)) < 0.1) then
                    depxcor(np,idep) = xdeg
                    deptcor(np,idep) = tmin
                depucor(np,idep) = slow(i+1)            
                end if
            end do
        end do d100
        xdeg = 2.*x  ! 105 ! actually in km 
        tmin = 2.*t        ! actually in s
        deltab(np) = xdeg  ! 110
        tttab(np) = tmin
    end do d200 !  end loop on ray parameter p

! create table
    d250: do idep=1, ndep
        icount = 0
        xold = -999.
        if (deptab(idep)==0.) then
            i2 = np
            go to 223
        end if
        d220: do i=1, np  ! upgoing rays from source
            x2 = depxcor(i, idep)
            if (x2 == -999.) exit d220
            if (x2 <= xold) exit d220     ! stop when heads inward
            t2 = deptcor(i, idep)
            icount = icount+1
            xsave(icount) = x2
            tsave(icount) = t2
            psave(icount) = -ptab(i)
            usave(icount) = depucor(i, idep)
            xold = x2
        end do d220
        i2 = i-1
223     continue
        d225: do i=i2, 1, -1  ! downgoing rays from source
            if (depxcor(i,idep) == -999.) cycle d225
            if (deltab(i) == -999.) cycle d225
            x2 = deltab(i) - depxcor(i, idep)
            t2 = tttab(i) - deptcor(i, idep)
            icount = icount+1
            xsave(icount) = x2
            tsave(icount) = t2
            psave(icount) = ptab(i)
            usave(icount) = depucor(i, idep)
            xold = x2
        end do d225
        ncount = icount

        ndel = int((del2-del1)/del3) + 1
        d240: do idel=1, ndel ! do 240
            del = del1 + del3*real(idel-1)
            delttab(idel) = del
            tt(idel,idep) = 999.
            d230: do i=2, ncount ! do 230
                x1 = xsave(i-1)
                x2 = xsave(i)
                if (x1 > del .or. x2 < del) cycle d230
                if (psave(i) > 0. .and. psave(i) < plongcut) cycle d230
                frac = (del-x1)/(x2-x1)
                t1 = tsave(i-1) + frac*(tsave(i) - tsave(i-1))
                if (t1 < tt(idel, idep)) then
                    tt(idel, idep) = t1
                    scr1 = psave(i)/usave(i)
                    angle = asin(scr1)*degrad
                    if (angle < 0.) then
                        angle = -angle
                    else
                        angle = 180. - angle
                    end if
                    table(idel, idep, itab) = angle
                end if
            end do d230
        end do d240
    end do d250
end do d300

    if (delttab(1)==0.) then
        do idep=1, ndep
            table(1, idep, itab) = 0. ! straight up at zero range
        end do
    end if
      
    return  ! 999

end subroutine

! ------------------------------------------------------------ !

! subroutine GET_TTS obtains the takeoff angle for a velocity model
! at a specified range and earthquake depth by interpolating
! from a table of takeoff angles.  
!    Inputs:    ip     =  index number for model (up to nindex)
!               del    =  range
!               qdep   =  earthquake depth
!    Returns:   tt     =  takeoff angle (degrees)
!               iflag  = -1 if outside depth range
!                      =  0 for interpolation
!                      =  1 for extrapolation in range
!
subroutine GET_TTS(ip, del, qdep, tt, iflag)
!f2py intent(in) ip
!f2py intent(in) del
!f2py intent(in) qdep
!f2py intent(out) tt
!f2py intent(out) iflag

    include 'vel.inc'
    real :: t(nx0, nd0, nindex), x(nx0), d(nd0)
    integer :: nx, nd

!  common block:
!    t(nx0,nd0,nindex)  =  takeoff angle tables
!          x(nx0)  =  list of ranges for tables
!          d(nd0)  =  list of source depths for tables
!           nx     =  number of distance points in table
!           nd     =  number of source depths in table
    common /angtable/ t, x, d, nx, nd

    !
    ! check if outside depth range
    if (qdep < d(1) .or. qdep > d(nd0)) then
        iflag = -1
        tt = 999
        print *,'*** event outside of velocity table depth range', &
            ' event depth=',qdep,' table range=',d(1),d(nd0)
        return
    end if
    ! first check to see if interpolation alone will work
    d30: do id=2, nd
        if (d(id) < qdep) cycle d30
        id1 = id-1
        id2 = id
        go to 35
    end do d30
    id1 = nd-1
    id2 = nd
35  continue
    d35: do ix=2, nx
        if (x(ix) < del) cycle d35
        ix1 = ix-1
        ix2 = ix
        go to 37
    end do d35
    ix1 = nx-1
    ix2 = nx
37  if (t(ix1,id1,ip)==0.) go to 50
    if (t(ix1,id2,ip)==0.) go to 50
    if (t(ix2,id1,ip)==0.) go to 50
    if (t(ix2,id2,ip)==0.) go to 50
    if (x(ix2) < del) go to 50
    iflag = 0
    xfrac = (del-x(ix1))/(x(ix2)-x(ix1))
    t1 = t(ix1,id1,ip) + xfrac*(t(ix2,id1,ip) - t(ix1,id1,ip))
    t2 = t(ix1,id2,ip) + xfrac*(t(ix2,id2,ip) - t(ix1,id2,ip))
    dfrac = (qdep-d(id1))/(d(id2)-d(id1))
    tt = t1 + dfrac*(t2-t1)
    return
    ! extrapolate to get tt
50  iflag = 1
    xoffmin1 = 999.
    xoffmin2 = 999.
    ixbest1 = 999
    ixbest2 = 999
    d60: do ix=2, nx
        if (t(ix-1,id1,ip)==0) go to 55
        if (t(ix,id1,ip)==0) go to 55
        xoff = abs((x(ix-1) + x(ix))/2. - del)
        if (xoff < xoffmin1) then
            xoffmin1 = xoff
            ixbest1 = ix
        end if
55      if (t(ix-1,id2,ip)==0) cycle d60
        if (t(ix,id2,ip)==0) cycle d60
        xoff = abs((x(ix-1) + x(ix))/2. - del)
        if (xoff < xoffmin2) then
            xoffmin2 = xoff
            ixbest2 = ix
        end if
    end do d60
    if (ixbest1==999 .or. ixbest2==999) then
        iflag = -1
        tt = 999
        return
    end if
    xfrac1 = (del - x(ixbest1 - 1)) / (x(ixbest1) - x(ixbest1 - 1))
    t1 = t(ixbest1-1,id1,ip)
    t2 = t(ixbest1,id1,ip)
    tt1 = t1 + xfrac1*(t2-t1)
    xfrac2 = (del - x(ixbest2 - 1)) / (x(ixbest2) - x(ixbest2 - 1))
    t1 = t(ixbest2-1,id2,ip)
    t2 = t(ixbest2,id2,ip)
    tt2 = t1 + xfrac2*(t2-t1)
    dfrac = (qdep-d(id1))/(d(id2)-d(id1))
    tt = tt1 + dfrac*(tt2-tt1)
    return
end subroutine

! ------------------------------------------------------------ !

! LAYERTRACE calculates the travel time and range offset
! for ray tracing through a single layer.
!
! Input:    p     =  horizontal slowness
!           h     =  layer thickness
!           utop  =  slowness at top of layer
!           ubot  =  slowness at bottom of layer
!           imth  =  interpolation method
!                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)
!                         = 2,  v(z) = a - b*z
!                         = 3,  v(z) = a*exp(-b*z)
!
! Returns:  dx    =  range offset
!           dt    =  travel time
!           irtr  =  return code
!                 = -1, zero thickness layer
!                 =  0,  ray turned above layer
!                 =  1,  ray passed through layer
!                 =  2,  ray turned within layer, 1 segment counted
!
! Note:  This version does calculation in double precision,
!        but all i/o is still single precision
!
subroutine LAYERTRACE(p1,h1,utop1,ubot1,imth,dx1,dt1,irtr)
    implicit real(kind=8) (a-h,o-z)
    real(kind=4) :: p1, h1, utop1, ubot1, dx1, dt1
    p = dble(p1)
    h = dble(h1)
    utop = dble(utop1)
    ubot = dble(ubot1)
!
    if (h == 0.) then      !check for zero thickness layer
        dx1 = 0.
        dt1 = 0.
        irtr = -1
        return         
    end if
!
    u = utop
    y = u - p
    if (y <= 0.) then   !complex vertical slowness
        dx1 = 0.
        dt1 = 0.
        irtr = 0
        return
    end if
!
    q = y*(u + p)
    qs = dsqrt(q)
!
    ! special function needed for integral at top of layer
    if (imth == 2) then
        y = u + qs
        if (p /= 0.) y = y/p
        qr = dlog(y)
    else if (imth == 3) then
         qr = atan2(qs, p)
    end if      
!
    if (imth == 1) then
        b = -(utop**2 - ubot**2)/(2.*h)
    else if (imth == 2) then
        vtop = 1./utop
        vbot = 1./ubot
        b = -(vtop - vbot)/h
    else
        b = -dlog(ubot/utop)/h
    end if
!
    if (b == 0.) then     !constant velocity layer
        b = 1./h
        etau = qs
        ex = p/qs
        irtr = 1
        go to 160
    end if
!
    ! integral at upper limit, 1/b factor omitted until end
    if (imth == 1) then
        etau = -q*qs/3.
        ex = -qs*p
    else if (imth == 2) then
        ex = qs/u                 !*** - in some versions (wrongly)
        etau = qr-ex
        if (p /= 0.) ex = ex/p
    else
        etau = qs - p*qr
        ex = qr
    end if
!
    ! check lower limit to see if we have turning point
    u = ubot
    if (u <= p) then   !if turning point,
        irtr = 2          !then no contribution
        go to 160       !from bottom point
    end if 
    irtr = 1
    q = (u-p)*(u+p)
    qs = dsqrt(q)
!
    if (imth == 1) then
        etau = etau + q*qs/3.
        ex = ex + qs*p
    else if (imth == 2) then
        y = u + qs
        z = qs/u
        etau = etau + z
        if (p /= 0.) then
            y = y/p
            z = z/p
        end if
        qr = dlog(y)
        etau = etau - qr
        ex = ex - z
    else
        qr = atan2(qs, p)
        etau = etau - qs + p*qr
        ex = ex - qr
    end if      
!
160 dx = ex/b
    dtau = etau/b
    dt = dtau + p*dx     !convert tau to t
!
    dx1 = sngl(dx)
    dt1 = sngl(dt)
    return
end subroutine
