c
c  cross product of two vectors, sets v3 = v1 x v2
c
      subroutine CROSS(v1,v2,v3)
      real v1(3),v2(3),v3(3)

cf2py intent(in) v1
cf2py intent(in) v2
cf2py intent(out) v3

      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)   
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      return
      end   

c ------------------------------------------------------------ c

c
c  transforms spherical co-ordinates to cartesian
c
      subroutine TO_CAR(the,phi,r,x,y,z)
      degrad=3.1415927/180.
      z=-r*cos(the*degrad)
      x=r*sin(the*degrad)*cos(phi*degrad)
      y=r*sin(the*degrad)*sin(phi*degrad)
      return
      end


c ------------------------------------------------------------ c

c subroutine FPCOOR gets fault normal vector,fnorm, and slip 
c vector, slip, from (strike,dip,rake) or vice versa.
c   idir = 1 compute fnorm,slip
c   idir = c compute strike,dip,rake
c Reference:  Aki and Richards, p. 115
c   uses (x,y,z) coordinate system with x=north, y=east, z=down
      subroutine FPCOOR(strike,dip,rake,fnorm,slip,idir)
      real fnorm(3),slip(3),phi,del,lam,a,clam,slam
      degrad=180./3.1415927
      pi=3.1415927
      phi=strike/degrad
      del=dip/degrad
      lam=rake/degrad
      if (idir.eq.1) then
         fnorm(1)=-sin(del)*sin(phi)
         fnorm(2)= sin(del)*cos(phi)
         fnorm(3)=-cos(del)
         slip(1)= cos(lam)*cos(phi)+cos(del)*sin(lam)*sin(phi)
         slip(2)= cos(lam)*sin(phi)-cos(del)*sin(lam)*cos(phi)
         slip(3)=-sin(lam)*sin(del)
      else
         if ((1.-abs(fnorm(3))).le.1e-7) then
c          print *,'***FPCOOR warning, horz fault, strike undefined'
           del=0.
           phi=atan2(-slip(1),slip(2))
           clam=cos(phi)*slip(1)+sin(phi)*slip(2)
           slam=sin(phi)*slip(1)-cos(phi)*slip(2)
           lam=atan2(slam,clam)
         else
           phi=atan2(-fnorm(1),fnorm(2))
           a=sqrt(fnorm(1)*fnorm(1)+fnorm(2)*fnorm(2))
           del=atan2(a,-fnorm(3))
           clam=cos(phi)*slip(1)+sin(phi)*slip(2)
           slam=-slip(3)/sin(del)
           lam=atan2(slam,clam)
           if (del.gt.(0.5*pi)) then
             del=pi-del
             phi=phi+pi
             lam=-lam
           end if
         end if
         strike=phi*degrad
         if (strike.lt.0.) strike=strike+360.
         dip=del*degrad
         rake=lam*degrad
         if (rake.le.-180.) rake=rake+360.
         if (rake.gt.180.) rake=rake-360.
      end if
      return
      end


c ------------------------------------------------------------ c

c normally-distributed random numbers, from numerical recipes      
      subroutine RAN_NORM(fran)
      save jran,ifirst
cf2py intent(out) fran
      im=120050                          !overflow at 2**28
      ia=2311
      ic=25367
      if (ifirst.ne.12345) then
         jran=314159
         ifirst=12345
      end if
      fran=0
      do 10 i=1,12
        jran=mod(jran*ia+ic,im)
        fran=fran+(float(jran)/float(im))
10    continue
      fran=fran-6.
      return
      end

