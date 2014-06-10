!
!  cross product of two vectors, sets v3 = v1 x v2
!
subroutine CROSS(v1, v2, v3)
!f2py intent(in) v1
!f2py intent(in) v2
!f2py intent(out) v3
    real, dimension(3) :: v1, v2, v3

    v3(1) = v1(2)*v2(3) - v1(3)*v2(2)   
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
    return
end subroutine  
!--------------------------------------------------------------!


!
!  transforms spherical co-ordinates to cartesian
!
subroutine TO_CAR(the, phi, r, x, y, z)
!f2py intent(in) the
!f2py intent(in) phi
!f2py intent(in) r
!f2py intent(out) x
!f2py intent(out) y
!f2py intent(out) z
      degrad = 3.1415927/180.
      z= -r * cos(the*degrad)
      x= r * sin(the*degrad) * cos(phi*degrad)
      y= r * sin(the*degrad) * sin(phi*degrad)
      return
end
!--------------------------------------------------------------!


! subroutine FPCOOR gets fault normal vector,fnorm, and slip 
! vector, slip, from (strike,dip,rake) or vice versa.
!   idir = 1 compute fnorm,slip
!   idir = c compute strike,dip,rake
! Reference:  Aki and Richards, p. 115
!   uses (x,y,z) coordinate system with x=north, y=east, z=down
subroutine FPCOOR(strike, dip, rake, fnorm, slip, idir)
!f2py intent(in) idir
!f2py intent(in, out) strike
!f2py intent(in, out) dip
!f2py intent(in, out) rake
!f2py intent(in, out) fnorm
!f2py intent(in, out) slip
    real ::  fnorm(3), slip(3), phi, del, lam, a, clam, slam
    
    degrad = 180./3.1415927
    pi = 3.1415927
    phi = strike/degrad
    del = dip/degrad
    lam = rake/degrad
    if (idir == 1) then
        fnorm(1) = -sin(del) * sin(phi)
        fnorm(2) = sin(del) * cos(phi)
        fnorm(3)= -cos(del)
        slip(1)= cos(lam)*cos(phi) + cos(del)*sin(lam)*sin(phi)
        slip(2)= cos(lam)*sin(phi) - cos(del)*sin(lam)*cos(phi)
        slip(3)= -sin(lam) * sin(del)
    else
        if ((1.-abs(fnorm(3))) <= 1e-7) then
!           print *,'***FPCOOR warning, horz fault, strike undefined'
            del = 0.
            phi = atan2(-slip(1), slip(2))
            clam = cos(phi)*slip(1) + sin(phi)*slip(2)
            slam = sin(phi)*slip(1) - cos(phi)*slip(2)
            lam = atan2(slam, clam)
        else
            phi = atan2(-fnorm(1), fnorm(2))
            a = sqrt(fnorm(1)*fnorm(1) + fnorm(2)*fnorm(2))
            del = atan2(a, -fnorm(3))
            clam = cos(phi)*slip(1) + sin(phi)*slip(2)
            slam= -slip(3)/sin(del)
            lam = atan2(slam, clam)
            if (del > (0.5*pi)) then
                del = pi - del
                phi = phi + pi
                lam = -lam
            end if
         end if
         strike = phi*degrad
         if (strike < 0.) then 
            strike = strike+360.
         end if
         dip = del*degrad
         rake = lam*degrad
         if (rake <= -180.) then
            rake = rake+360.
         end if
         if (rake >= 180.) then 
            rake = rake-360.
         end if
    end if
    return
end subroutine
!------------------------------------------------------------!

!
! normally-distributed random numbers, from numerical recipes      
!
subroutine RAN_NORM(fran)
!f2py intent(out) fran
    save jran, ifirst
    im = 120050  ! overflow at 2**28
    ia = 2311
    ic = 25367
    if (ifirst /= 12345) then
        jran = 314159
        ifirst = 12345
    end if
    fran = 0
    do i = 1, 12
        jran = mod(jran*ia+ic,im)
        fran = fran + (float(jran)/float(im))
    end do
    fran=fran-6.
    return
end subroutine

