! subroutine GET_MISF finds the percent of misfit polarities for a given mechanism  
!    Inputs:    npol   = number of polarity observations
!               p_azi_mc(npol) = azimuths
!               p_the_mc(npol) = takeoff angles
!               p_pol(npol)  = polarity observations
!               p_qual(npol) = quality of polarity observations
!               str_avg,dip_avg,rak_avg = mechanism
!    Outputs:   mfrac = weighted fraction misfit polarities (like FPFIT)
!               stdr = station distribution ratio (like FPFIT)

subroutine GET_MISF(npol, p_azi_mc, p_the_mc, p_pol, p_qual, str_avg, &
                    dip_avg, rak_avg, mfrac, stdr)
!f2py intent(in) npol
!f2py intent(in) p_azi_mc      
!f2py intent(in) p_the_mc  
!f2py intent(in) p_pol  
!f2py intent(in) p_qual
!f2py intent(in) str_avg
!f2py intent(in) dip_avg
!f2py intent(in) rak_avg
!f2py intent(out) mfrac      
!f2py intent(out) stdr      

    integer :: k, npol, p_pol(npol), p_qual(npol)
    real :: str_avg, dip_avg, rak_avg
    real :: strike, dip, rake, mfrac, qcount, azi, toff, pol, wt, wo
    real, dimension(3) :: a, b, bb1, bb2, bb3
    real, dimension(3, 3) :: M
    real, dimension(npol) :: p_azi_mc, p_the_mc
    real, parameter :: rad = 3.14159265/180.

    strike = str_avg*rad
    dip = dip_avg*rad
    rake = rak_avg*rad
      
    M(1,1) = -sin(dip)*cos(rake)*sin(2*strike) - sin(2*dip)*sin(rake)* & 
        sin(strike)*sin(strike)
    M(2,2) = sin(dip)*cos(rake)*sin(2*strike) - sin(2*dip)*sin(rake)* &
        cos(strike)*cos(strike)
    M(3,3) = sin(2*dip) * sin(rake)
    M(1,2) = sin(dip)*cos(rake)*cos(2*strike) + 0.5*sin(2*dip)*sin(rake)* & 
        sin(2*strike)
    M(2,1) = M(1,2)
    M(1,3) = -cos(dip)*cos(rake)*cos(strike) - cos(2*dip)*sin(rake)* & 
        sin(strike)
    M(3,1) = M(1,3)
    M(2,3) = -cos(dip)*cos(rake)*sin(strike) + cos(2*dip)*sin(rake)* & 
        cos(strike)
    M(3,2) = M(2,3)
    
    call FPCOOR(strike, dip, rake, bb3, bb1, 1)
    call CROSS(bb3, bb1, bb2)
    
    mfrac = 0.
    qcount = 0.
    scount = 0.
      
    do k=1, npol
        call TO_CAR(p_the_mc(k), p_azi_mc(k), 1., p_a1, p_a2, p_a3)
        p_b1 = bb1(1)*p_a1 + bb1(2)*p_a2 + bb1(3)*p_a3 
        p_b3 = bb3(1)*p_a1 + bb3(2)*p_a2 + bb3(3)*p_a3
        p_proj1 = p_a1 - p_b3*bb3(1)
        p_proj2 = p_a2 - p_b3*bb3(2)
        p_proj3 = p_a3 - p_b3*bb3(3)
        plen = sqrt(p_proj1*p_proj1+p_proj2*p_proj2+p_proj3*p_proj3)
        
        pp_b1 = bb1(1)*p_proj1 + bb1(2)*p_proj2 + bb1(3)*p_proj3
        pp_b2 = bb2(1)*p_proj1 + bb2(2)*p_proj2 + bb2(3)*p_proj3
        phi = atan2(pp_b2, pp_b1)
        theta = acos(p_b3)
        p_amp = abs(sin(2*theta) * cos(phi))     
        wt = sqrt(p_amp)
        azi = rad*p_azi_mc(k)
        toff = rad*p_the_mc(k)        
        a(1) = sin(toff) * cos(azi)
        a(2) = sin(toff) * sin(azi)
        a(3) = -cos(toff)
        do in=1, 3
            b(in) = 0
            do jn=1, 3
                b(in) = b(in) + M(in,jn)*a(jn)
            end do
        end do
        if ((a(1)*b(1) + a(2)*b(2) + a(3)*b(3)) < 0) then
            pol = -1
        else
            pol = 1
        end if
        if (p_qual(k) == 0) then
            wo = 1
        else
            wo = 0.5
        end if
        if ((pol*p_pol(k)) < 0) then
            mfrac = mfrac + wt*wo
        end if
        qcount = qcount + wt*wo
        scount = scount + wo
    end do
    mfrac = mfrac/qcount
    stdr = qcount/scount
      
    return 
end subroutine

! --------------------------------------------------------------- !


! subroutine GET_GAP finds the maximum azimuthal and takeoff angle gaps  
!    Inputs:    npol   = number of polarity observations
!               p_azi_mc(npol) = azimuths
!               p_the_mc(npol) = takeoff angles
!    Outputs:   magap  = maximum azimuthal gap
!               mpgap  = maximum takeoff angle gap

subroutine GET_GAP(npol, p_azi_mc, p_the_mc, magap, mpgap)
!f2py intent(in) npol
!f2py intent(in) p_azi_mc      
!f2py intent(in) p_the_mc  
!f2py intent(out) magap  
!f2py intent(out) mpgap  

    include 'param.inc'
    real, dimension(npol) :: p_azi_mc, p_the_mc
    real, dimension(npick0) :: p2_azi, p2_the

    do k=1, npol
        if (p_the_mc(k) > 90) then
            p2_the(k) = 180. - p_the_mc(k)
            p2_azi(k) = p_azi_mc(k) - 180.
            if (p2_azi(k) < 0) then
                p2_azi(k) = p2_azi(k) + 360.
            end if
        else
            p2_the(k)=p_the_mc(k)
            p2_azi(k)=p_azi_mc(k)
        end if
    end do
    call sort(npol, p2_azi)
    call sort(npol, p2_the)
    magap = 0
    mpgap = 0
    do k=2, npol
        if (p2_azi(k) - p2_azi(k-1) > magap) then
            magap = p2_azi(k) - p2_azi(k-1)
        end if
        if (p2_the(k) - p2_the(k-1) > mpgap) then
            mpgap = p2_the(k) - p2_the(k-1)
        end if
    end do
    if (p2_azi(1) - p2_azi(npol) + 360 > magap) then
        magap = p2_azi(1) - p2_azi(npol) + 360
    end if
    if (90. - p2_the(npol) > mpgap) then
        mpgap = 90. - p2_the(npol)
    end if
    if (p2_the(1) > mpgap) then
        mpgap = p2_the(1)
    end if
      
    return 
end subroutine

! --------------------------------------------------------------- !

subroutine sort(N, RA)
      
    real, dimension(N) :: RA
      
    if (n==0) then
        print *,'***n=0 in SORT'
        return
    end if
    L =N/2 + 1
    IR = N
10  continue
    if(L > 1) then
        L = L-1
        RRA = RA(L)
    else
        RRA = RA(IR)
        RA(IR) = RA(1)
        IR = IR - 1
        if (IR==1) then
            RA(1) = RRA
            return
        endif
    endif
    I = L
    J = L + L
20  if (J <= IR) then
        if (J < IR) then
            if (RA(J) < RA(J+1)) J = J+1
        endif
        if (RRA < RA(J)) then
            RA(I) = RA(J)
            I = J
            J = J+J
        else
            J = IR+1
        endif
        goto 20
    endif
    RA(I) = RRA
    goto 10
end subroutine
