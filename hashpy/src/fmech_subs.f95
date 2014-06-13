! subroutine FOCALMC performs grid search to find acceptable focal mechanisms,
!                    for multiple trials of ray azimuths and takeoff angles.
!                    Acceptable mechanisms are those with less than "ntotal"
!                    misfit polarities, or the minimum plus "nextra" if this
!                    is greater.
!
!  Inputs:  
!           p_azi_mc(npsta,nmc)  =  azimuth to station from event (deg. E of N)
!           p_the_mc(npsta,nmc)  =  takeoff angle (from vert, up=0, <90 upgoing, >90 downgoing)
!           p_pol(npsta)  =  first motion, 1=up, -1=down
!           p_qual(npsta) =  quality, 0=impulsive, 1=emergent
!           npsta  =  number of first motions
!           nmc    =  number of trials
!           dang   =  desired angle spacing for grid search
!           maxout =  maximum number of fault planes to return:
!                     if more are found, a random selection will be returned
!           nextra =  number of additional misfits allowed above minimum
!           ntotal =  total number of allowed misfits
!  Outputs: 
!           nf     =  number of fault planes found
!           strike(min(maxout,nf)) = strike
!           dip(min(maxout,nf))    = dip
!           rake(min(maxout,nf))   = rake
!           faults(3,min(maxout,nf)) = fault normal vector
!           slips(3,min(maxout,nf))  = slip vector
!
!
subroutine FOCALMC(p_azi_mc, p_the_mc, p_pol, p_qual, npsta, nmc, dang, &
                   maxout, nextra, ntotal, nf, strike, dip, rake, &
                   faults, slips)
!f2py intent(in)  p_azi_mc
!f2py intent(in)  p_the_mc
!f2py intent(in)  p_pol
!f2py intent(in)  p_qual
!f2py intent(in)  npsta
!f2py intent(in)  nmc
!f2py intent(in)  dang
!f2py intent(in)  maxout
!f2py intent(in)  nextra
!f2py intent(in)  ntotal
!f2py intent(out) nf
!f2py intent(out) strike
!f2py intent(out) dip
!f2py intent(out) rake
!f2py intent(out) faults
!f2py intent(out) slips

    include 'param.inc'
    include 'rot.inc'
     
    ! input and output arrays
    real, dimension(npick0,nmc0) ::  p_azi_mc, p_the_mc
    real, dimension(npick0) ::  p_a1, p_a2, p_a3
    real, dimension(3) :: faultnorm, slip
    real, dimension(3,nmax0) :: faults, slips
    real, dimension(nmax0) ::  strike, dip, rake
    integer, dimension(npsta) :: p_pol, p_qual
    save dangold, nrot, b1, b2, b3

    ! coordinate transformation arrays
    real, dimension(3,ncoor) :: b1, b2, b3
    real, dimension(3) :: bb1, bb2, bb3

    ! fit arrays
    integer :: fit(2,ncoor), nmiss01min(0:npick0)
    integer, dimension(ncoor) :: irotgood, irotgood2
    real :: fran

    real, parameter :: pi = 3.1415927
    real, parameter :: degrad = 180./pi

    if (maxout > nmax0) then
        maxout = nmax0
    end if

    ! Set up array with direction cosines for all coordinate transformations
    if (dang == dangold) go to 8  ! TODO: replace w/ if(dang/=dangold)then, endif @ 8? -MCW
    irot = 0
    d5: do ithe=0, int(90.1/dang)
        the = real(ithe)*dang
        rthe = the/degrad
        costhe = cos(rthe)
        sinthe = sin(rthe)
        fnumang = 360./dang
        numphi = nint(fnumang*sin(rthe))
        if (numphi /= 0) then
            dphi = 360./float(numphi)
        else
            dphi = 10000.
        end if
        d4: do iphi=0, int(359.9/dphi)
            phi = real(iphi)*dphi
            rphi = phi/degrad
            cosphi = cos(rphi)
            sinphi = sin(rphi)
            bb3(3) = costhe
            bb3(1) = sinthe*cosphi
            bb3(2) = sinthe*sinphi
            bb1(3) = -sinthe
            bb1(1) = costhe*cosphi
            bb1(2) = costhe*sinphi
            call CROSS(bb3, bb1, bb2)
            d3: do izeta=0, int(179.9/dang)
                zeta = real(izeta)*dang
                rzeta = zeta/degrad
                coszeta = cos(rzeta)
                sinzeta = sin(rzeta)
                irot = irot+1
                if (irot > ncoor) then
                    print *,'***FOCAL error: # of rotations too big'
                    return
                end if
                b3(3,irot) = bb3(3)
                b3(1,irot) = bb3(1)
                b3(2,irot) = bb3(2)
                b1(1,irot) = bb1(1)*coszeta + bb2(1)*sinzeta
                b1(2,irot) = bb1(2)*coszeta + bb2(2)*sinzeta                
                b1(3,irot) = bb1(3)*coszeta + bb2(3)*sinzeta
                b2(1,irot) = bb2(1)*coszeta - bb1(1)*sinzeta
                b2(2,irot) = bb2(2)*coszeta - bb1(2)*sinzeta                
                b2(3,irot) = bb2(3)*coszeta - bb1(3)*sinzeta
            end do d3
        end do d4
    end do d5
    nrot = irot
    dangold = dang
8   continue

    do irot=1, nrot
        irotgood(irot) = 0
    end do

    ! loop over multiple trials
    d430: do im=1, nmc 

        ! Convert data to Cartesian coordinates
        d40: do i=1, npsta
            call TO_CAR(p_the_mc(i,im), p_azi_mc(i,im), 1., p_a1(i), p_a2(i), p_a3(i))
        end do d40

        ! find misfit for each solution and minimum misfit
        nmiss0min = 999
        nmissmin = 999
        do i=0, npick0
           nmiss01min(i) = 999
        end do
        d420: do irot=1, nrot  
            nmiss = 0
            nmiss0 = 0           
            d400: do ista=1,npsta
                p_b1 = b1(1,irot)*p_a1(ista) + b1(2,irot)*p_a2(ista) + b1(3,irot)*p_a3(ista) 
                p_b3 = b3(1,irot)*p_a1(ista) + b3(2,irot)*p_a2(ista) + b3(3,irot)*p_a3(ista) 
                prod = p_b1*p_b3
                ipol = -1
                if (prod > 0.) ipol = 1    ! predicted polarization
                if (ipol /= p_pol(ista)) then
                    nmiss = nmiss+1                       
                    if (p_qual(ista) == 0) nmiss0 = nmiss0+1
                end if
            end do d400
            fit(1,irot) = nmiss0           ! misfit impulsive polarities
            fit(2,irot) = nmiss            ! total misfit polarities
            if (nmiss0 < nmiss0min) nmiss0min = nmiss0
            if (nmiss < nmissmin) nmissmin = nmiss
            if (nmiss < nmiss01min(nmiss0)) then
                nmiss01min(nmiss0) = nmiss
            end if
        end do d420
        
        ! choose fit criteria
        if (nmiss0min == 0) then 
            nmiss0max = ntotal
            nmissmax = ntotal
        else
            nmiss0max = ntotal
            nmissmax = npsta
        end if
        if (nmiss0max < nmiss0min + nextra) then
            nmiss0max = nmiss0min + nextra
        end if
        if (nmissmax < nmiss01min(nmiss0min) + nextra) then
            nmissmax = nmiss01min(nmiss0min) + nextra
        end if

        ! loop over rotations - find those meeting fit criteria
        d440: do irot=1, nrot        
            nmiss0 = fit(1,irot)
            nmiss = fit(2,irot)
            if ((nmiss0 <= nmiss0max) .and. (nmiss <= nmissmax)) then
                irotgood(irot) = 1
            end if
        end do d440

    end do d430

    nfault = 0
    do irot=1, nrot
        if (irotgood(irot) > 0) then
            nfault = nfault+1
            irotgood2(nfault) = irot
        end if
    end do
 
    !  Select output solutions  
    nf = 0      
    if (nfault <= maxout) then
        do i=1, nfault
            irot = irotgood2(i)
            nf = nf+1
            faultnorm(1) = b3(1,irot)
            faultnorm(2) = b3(2,irot)
            faultnorm(3) = b3(3,irot)
            slip(1) = b1(1,irot)
            slip(2) = b1(2,irot)
            slip(3) = b1(3,irot)
            do m=1, 3
                faults(m,nf) = faultnorm(m)
                slips(m,nf) = slip(m)
            end do
            call FPCOOR(s1, d1, r1, faultnorm, slip, 2)
            strike(nf) = s1
            dip(nf) = d1
            rake(nf) = r1
        end do
    else
        d441: do i=1, 99999
            fran = rand(0)
            iscr = nint(fran*float(nfault) + 0.5)
            if (iscr < 1) iscr = 1
            if (iscr > nfault) iscr = nfault
            if (irotgood2(iscr) <= 0) cycle d441
            irot = irotgood2(iscr)
            irotgood2(iscr) = -1
            nf = nf+1
            faultnorm(1) = b3(1,irot)
            faultnorm(2) = b3(2,irot)
            faultnorm(3) = b3(3,irot)
            slip(1) = b1(1,irot)
            slip(2) = b1(2,irot)
            slip(3) = b1(3,irot)
            do m=1, 3
              faults(m,nf) = faultnorm(m)
              slips(m,nf) = slip(m)
            end do
            call FPCOOR(s1, d1, r1, faultnorm, slip, 2)
            strike(nf) = s1
            dip(nf) = d1
            rake(nf) = r1
            if (nf == maxout) exit d441
        end do d441
    end if

    return 
end subroutine
