c subroutine GET_MISF finds the percent of misfit polarities for a given mechanism  
c    Inputs:    npol   = number of polarity observations
c               p_azi_mc(npol) = azimuths
c               p_the_mc(npol) = takeoff angles
c               p_pol(npol)  = polarity observations
c               p_qual(npol) = quality of polarity observations
c               str_avg,dip_avg,rak_avg = mechanism
c    Outputs:   mfrac = weighted fraction misfit polarities (like FPFIT)
c               stdr = station distribution ratio (like FPFIT)

      subroutine GET_MISF(npol,p_azi_mc,p_the_mc,p_pol,p_qual,str_avg,
     &                  dip_avg,rak_avg,mfrac,stdr)

      dimension p_azi_mc(npol),p_the_mc(npol)
      real str_avg,dip_avg,rak_avg,M(3,3),a(3),b(3)
      real strike,dip,rake,mfrac,qcount,azi,toff,pol,wt,wo
      integer k,npol,p_pol(npol),p_qual(npol)
      real bb1(3),bb2(3),bb3(3)

cf2py intent(in) npol
cf2py intent(in) p_azi_mc      
cf2py intent(in) p_the_mc  
cf2py intent(in) p_pol  
cf2py intent(in) p_qual
cf2py intent(in) str_avg
cf2py intent(in) dip_avg
cf2py intent(in) rak_avg
cf2py intent(out) mfrac      
cf2py intent(out) stdr      

      rad=3.14159265/180.

      strike=str_avg*rad
      dip=dip_avg*rad
      rake=rak_avg*rad
      
      M(1,1)=-sin(dip)*cos(rake)*sin(2*strike)-sin(2*dip)*sin(rake)*
     & sin(strike)*sin(strike)
      M(2,2)=sin(dip)*cos(rake)*sin(2*strike)-sin(2*dip)*sin(rake)*
     & cos(strike)*cos(strike)
      M(3,3)=sin(2*dip)*sin(rake)
      M(1,2)=sin(dip)*cos(rake)*cos(2*strike)+0.5*sin(2*dip)*sin(rake)*
     & sin(2*strike)
      M(2,1)=M(1,2)
      M(1,3)=-cos(dip)*cos(rake)*cos(strike)-cos(2*dip)*sin(rake)*
     & sin(strike)
      M(3,1)=M(1,3)
      M(2,3)=-cos(dip)*cos(rake)*sin(strike)+cos(2*dip)*sin(rake)*
     & cos(strike)
      M(3,2)=M(2,3)
      call FPCOOR(strike,dip,rake,bb3,bb1,1)
      call CROSS(bb3,bb1,bb2)
      mfrac=0.
      qcount=0.
      scount=0.
      
      do 600 k=1,npol
          call TO_CAR(p_the_mc(k),p_azi_mc(k),1.,p_a1,
     &                p_a2,p_a3)
          p_b1= bb1(1)*p_a1
     &              +bb1(2)*p_a2
     &              +bb1(3)*p_a3 
          p_b3= bb3(1)*p_a1
     &              +bb3(2)*p_a2
     &              +bb3(3)*p_a3
          p_proj1=p_a1-p_b3*bb3(1)
          p_proj2=p_a2-p_b3*bb3(2)
          p_proj3=p_a3-p_b3*bb3(3)
          plen=sqrt(p_proj1*p_proj1+p_proj2*p_proj2+
     &                    p_proj3*p_proj3)
          pp_b1=bb1(1)*p_proj1+bb1(2)*p_proj2
     &              +bb1(3)*p_proj3
          pp_b2=bb2(1)*p_proj1+bb2(2)*p_proj2
     &              +bb2(3)*p_proj3
          phi=atan2(pp_b2,pp_b1)
          theta=acos(p_b3)
          p_amp=abs(sin(2*theta)*cos(phi))     
          wt=sqrt(p_amp)
         azi=rad*p_azi_mc(k)
         toff=rad*p_the_mc(k)        
         a(1)=sin(toff)*cos(azi)
         a(2)=sin(toff)*sin(azi)
         a(3)=-cos(toff)
         do in=1,3
           b(in)=0
           do jn=1,3
              b(in)=b(in)+M(in,jn)*a(jn)
           end do
         end do
         if ((a(1)*b(1)+a(2)*b(2)+a(3)*b(3)).lt.0) then
           pol=-1
         else
           pol=1
         end if
         if (p_qual(k).eq.0) then
           wo=1
         else
           wo=0.5
         end if
         if ((pol*p_pol(k)).lt.0) then
           mfrac=mfrac+wt*wo
         end if
         qcount=qcount+wt*wo
         scount=scount+wo
600   continue
      mfrac=mfrac/qcount
      stdr=qcount/scount
      
      return 
      end

c --------------------------------------------------------------- c


c subroutine GET_GAP finds the maximum azimuthal and takeoff angle gaps  
c    Inputs:    npol   = number of polarity observations
c               p_azi_mc(npol) = azimuths
c               p_the_mc(npol) = takeoff angles
c    Outputs:   magap  = maximum azimuthal gap
c               mpgap  = maximum takeoff angle gap

      subroutine GET_GAP(npol,p_azi_mc,p_the_mc,magap,mpgap)

      include 'param.inc'
      dimension p_azi_mc(npol),p_the_mc(npol)
      real p2_azi(npick0),p2_the(npick0)

cf2py intent(in) npol
cf2py intent(in) p_azi_mc      
cf2py intent(in) p_the_mc  
cf2py intent(out) magap  
cf2py intent(out) mpgap  

      do 403 k=1,npol
        if (p_the_mc(k).gt.90) then
          p2_the(k)=180.-p_the_mc(k)
          p2_azi(k)=p_azi_mc(k)-180.
          if (p2_azi(k).lt.0) p2_azi(k)=p2_azi(k)+360.
        else
          p2_the(k)=p_the_mc(k)
          p2_azi(k)=p_azi_mc(k)
        end if
403   continue
      call sort(npol,p2_azi)
      call sort(npol,p2_the)
      magap=0
      mpgap=0
      do 405 k=2,npol
        if (p2_azi(k)-p2_azi(k-1).gt.magap) then
          magap=p2_azi(k)-p2_azi(k-1)
        end if
        if (p2_the(k)-p2_the(k-1).gt.mpgap) then
          mpgap=p2_the(k)-p2_the(k-1)
        end if
405   continue
      if (p2_azi(1)-p2_azi(npol)+360.gt.magap) then
        magap=p2_azi(1)-p2_azi(npol)+360
      end if
      if (90.-p2_the(npol).gt.mpgap) then
        mpgap=90.-p2_the(npol)
      end if
      if (p2_the(1).gt.mpgap) then
        mpgap=p2_the(1)
      end if
      
      return 
      end

c --------------------------------------------------------------- c

      subroutine sort(N,RA)
      
      dimension RA(N)
      
      if (n.eq.0) then
         print *,'***n=0 in SORT'
         return
      end if
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
