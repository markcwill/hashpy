! MK_TABLES creates tables of takeoff angles given 1D velocity models.
!
!   input: 
!     ind  -  index of the table to add
!   vmodel -  filename of the velocity model file
!   output:
!     ntab - number of tables (max nindex)
!
!   add a velocity model to the table,
!   velocity file format (free format):
!     depth(km) P_velocity(km/s)
      
subroutine MK_TABLE_ADD(itab, vmodel)

    include 'vel.inc'
    integer, parameter :: nray0 = 10001
    real, parameter :: degrad = 180./3.14159265

    integer, intent(in) :: itab
    character(len=100), intent(in) :: vmodel
    real, dimension(1000) :: z , alpha, slow
    real, dimension(20000) ::  xsave, tsave, psave, usave
    real, dimension(nray0) :: deltab, tttab, ptab
    real, dimension(nray0, nd0) :: depxcor, depucor, deptcor, tt
    
    !   common block:
    !       table(nx0,nd0,nindex)  =  takeoff angle table
    !           delttab(nx0)  =  list of ranges for tables
    !           deptab(nd0)  =  list of source depths for tables
    !           ndel       =  number of distance points in table
    !           ndep       =  number of source depths in table
    real :: table(nx0,nd0,nindex), delttab(nx0), deptab(nd0)
    integer :: ndel, ndep
    common /angtable/ table, delttab, deptab, ndel, ndep
    

    ! set up table
    qtempdep2 = dep2 + dep3/20.
    ndep = int((qtempdep2-dep1)/dep3) + 1
    do idep=1, ndep
        dep = dep1 + dep3*real(idep-1)
        deptab(idep) = dep
    end do

    ! read velocity model TODO: fix this -MCW
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
    two: do np=1, npmax
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
        one: do i=1, npts-1
            if (z(i) >= 9999) then 
                deltab(np) = -999.
                tttab(np) = -999.
                cycle two
            end if
            h = z(i+1) - z(i)
            if (h==0.) cycle one ! skip if interface
            call LAYERTRACE(p, h, slow(i), slow(i+1), imth, dx, dt, irtr)
            x = x+dx
            t = t+dt
            if (irtr==0 .or. irtr==2) exit one ! ray has turned
            xdeg = x                       ! actually in km 
            tmin = t                       ! actually in s
            do idep=1, ndep
                if (abs(z(i+1) - deptab(idep)) < 0.1) then
                    depxcor(np,idep) = xdeg
                    deptcor(np,idep) = tmin
                depucor(np,idep) = slow(i+1)            
                end if
            end do
        end do one
        xdeg = 2.*x  ! 105 ! actually in km 
        tmin = 2.*t        ! actually in s
        deltab(np) = xdeg  ! 110
        tttab(np) = tmin
    end do two !  end loop on ray parameter p

    ! create table
    twofifty: do idep=1, ndep
        icount = 0
        xold = -999.
        if (deptab(idep)==0.) then
            i2 = np
            !go to 223
        else !end if
        twotwenty: do i=1, np  ! upgoing rays from source
            x2 = depxcor(i, idep)
            if (x2 == -999.) exit twotwenty
            if (x2 <= xold) exit twotwenty     ! stop when heads inward
            t2 = deptcor(i, idep)
            icount = icount+1
            xsave(icount) = x2
            tsave(icount) = t2
            psave(icount) = -ptab(i)
            usave(icount) = depucor(i, idep)
            xold = x2
        end do twotwenty
        i2 = i-1
        end if  !223     continue
        
        twotwofive: do i=i2, 1, -1  ! downgoing rays from source
            if (depxcor(i,idep) == -999.) cycle twotwofive
            if (deltab(i) == -999.) cycle twotwofive
            x2 = deltab(i) - depxcor(i, idep)
            t2 = tttab(i) - deptcor(i, idep)
            icount = icount+1
            xsave(icount) = x2
            tsave(icount) = t2
            psave(icount) = ptab(i)
            usave(icount) = depucor(i, idep)
            xold = x2
        end do twotwofive
        ncount = icount

        ndel = int((del2-del1)/del3) + 1
        twoforty: do idel=1, ndel ! do 240
            del = del1 + del3*real(idel-1)
            delttab(idel) = del
            tt(idel,idep) = 999.
            twothirty: do i=2, ncount ! do 230
                x1 = xsave(i-1)
                x2 = xsave(i)
                if (x1 > del .or. x2 < del) cycle twothirty
                if (psave(i) > 0. .and. psave(i) < plongcut) cycle twothirty
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
            end do twothirty
        end do twoforty
    end do twofifty

    if (delttab(1)==0.) then
        do idep=1, ndep
            table(1, idep, itab) = 0. ! straight up at zero range
        end do
    end if
      
    return  ! 999
end subroutine

