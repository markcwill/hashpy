c subroutine FOCALAMP_MC performs grid search to find focal mechanisms, using both
c            P-polarity and S/P amplitude ratio information
c  Inputs:  
c           p_azi_mc(npsta,nmc)  =  azimuth to station from event (deg. E of N)
c           p_the_mc(npsta,nmc)  =  takeoff angle (from vert, <90 upgoing, >90 downgoing)
c           sp_amp(npsta)  =  amplitude ratios
c           p_pol(nspta)   =  P polarities
c           npsta  =  number of observations
c           nmc    =  number of trials
c           dang   =  desired angle spacing for grid search
c           maxout =  maximum number of fault planes to return:
c                     if more are found, a random selection will be returned
c           nextra =  number of polarity additional misfits allowed above minimum
c           ntotal =  total number of allowed polarity misfits
c           qextra =  additional amplitude misfit allowed above minimum
c           qtotal =  total allowed amplitude misfit
c  Outputs: 
c           nf     =  number of fault planes found
c           strike(min(maxout,nf)) = strike
c           dip(min(maxout,nf))    = dip
c           rake(min(maxout,nf))   = rake
c           faults(3,min(maxout,nf)) = fault normal vector
c           slips(3,min(maxout,nf))  = slip vector
c

      subroutine FOCALAMP_MC(p_azi_mc,p_the_mc,sp_amp,p_pol,npsta,nmc,
     &    dang,maxout,nextra,ntotal,qextra,qtotal,nf,strike,dip,rake,
     &    faults,slips)
     
      include 'param.inc'
      include 'rot.inc'
      parameter (ntab=180) 

c input and output arrays
      dimension p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0)
      real sp_amp(npsta)
      real p_a1(npick0),p_a2(npick0),p_a3(npick0)
      real faultnorm(3),slip(3),faults(3,nmax0),slips(3,nmax0)
      real strike(nmax0),dip(nmax0),rake(nmax0)
      integer p_pol(npsta)
      save dangold,nrot,b1,b2,b3    
      save amptable,phitable,thetable
cf2py intent(in)  p_azi_mc
cf2py intent(in)  p_the_mc
cf2py intent(in)  sp_amp
cf2py intent(in)  p_pol
cf2py intent(in)  npsta
cf2py intent(in)  nmc
cf2py intent(in)  dang
cf2py intent(in)  maxout
cf2py intent(in)  nextra
cf2py intent(in)  ntotal
cf2py intent(in)  qextra
cf2py intent(in)  qtotal
cf2py intent(out)  nf
cf2py intent(out)  strike
cf2py intent(out)  dip
cf2py intent(out)  rake
cf2py intent(out)  faults
cf2py intent(out)  slips

c coordinate transformation arrays
      real b1(3,ncoor),bb1(3)
      real b2(3,ncoor),bb2(3)
      real b3(3,ncoor),bb3(3)
c P and S amplitude arrays
      real amptable(2,ntab,2*ntab)
      real phitable(2*ntab+1,2*ntab+1)
      real thetable(2*ntab+1)
c misfit arrays
      real qmis(ncoor)
      integer nmis(ncoor)
      integer irotgood(ncoor),irotgood2(ncoor)
      
      pi=3.1415927
      degrad=180./pi

      if (maxout.gt.nmax0) then
        maxout=nmax0
      end if

c Set up array with direction cosines for all coordinate transformations
      if (dang.eq.dangold) go to 8
      irot=0
       do 5 ithe=0,int(90.1/dang)
          the=real(ithe)*dang
         rthe=the/degrad
         costhe=cos(rthe)
         sinthe=sin(rthe)
         fnumang=360./dang
         numphi=nint(fnumang*sin(rthe))
         if (numphi.ne.0) then
            dphi=360./float(numphi)
         else
            dphi=10000.
         end if
         do 4 iphi=0,int(359.9/dphi)
            phi=real(iphi)*dphi
            rphi=phi/degrad
            cosphi=cos(rphi)
            sinphi=sin(rphi)
            bb3(3)=costhe
            bb3(1)=sinthe*cosphi
            bb3(2)=sinthe*sinphi
            bb1(3)=-sinthe
            bb1(1)=costhe*cosphi
            bb1(2)=costhe*sinphi
            call CROSS(bb3,bb1,bb2)
            do 3 izeta=0,int(179.9/dang)
               zeta=real(izeta)*dang
               rzeta=zeta/degrad
               coszeta=cos(rzeta)
               sinzeta=sin(rzeta)
               irot=irot+1
               if (irot.gt.ncoor) then
                  print *,'***FOCAL error: # of rotations too big'
                  return
               end if
               b3(3,irot)=bb3(3)
               b3(1,irot)=bb3(1)
               b3(2,irot)=bb3(2)
               b1(1,irot)=bb1(1)*coszeta+bb2(1)*sinzeta
               b1(2,irot)=bb1(2)*coszeta+bb2(2)*sinzeta                
               b1(3,irot)=bb1(3)*coszeta+bb2(3)*sinzeta
               b2(1,irot)=bb2(1)*coszeta-bb1(1)*sinzeta
               b2(2,irot)=bb2(2)*coszeta-bb1(2)*sinzeta                
               b2(3,irot)=bb2(3)*coszeta-bb1(3)*sinzeta
3           continue
4        continue
5     continue
      nrot=irot
      dangold=dang
      
      astep=1./real(ntab)
      do 150 i=1,2*ntab+1
        bbb3=-1.+real(i-1)*astep
        thetable(i)=acos(bbb3)
        do 140 j=1,2*ntab+1
          bbb1=-1.+real(j-1)*astep
          phitable(i,j)=atan2(bbb3,bbb1)
          if (phitable(i,j).lt.0.) then
            phitable(i,j)=phitable(i,j)+2.*pi
          end if
140     continue
150   continue

      do 250 i=1,2*ntab
        phi=real(i-1)*pi*astep
        do 240 j=1,ntab
          theta=real(j-1)*pi*astep
          amptable(1,j,i)=abs(sin(2*theta)*cos(phi))                
          s1=cos(2*theta)*cos(phi)  
          s2=-cos(theta)*sin(phi)
          amptable(2,j,i)=sqrt(s1*s1+s2*s2)
240     continue
250   continue

8     continue

      do irot=1,nrot
        irotgood(irot)=0
      end do

c loop over multiple trials
      do 430 im=1,nmc 

c  Convert data to Cartesian coordinates
      do 40 i=1,npsta
        call TO_CAR(p_the_mc(i,im),p_azi_mc(i,im),1.,
     &               p_a1(i),p_a2(i),p_a3(i))
40    continue

c  find misfit for each solution and minimum misfit
         nmis0min=1e5
         qmis0min=1.0e5
         do 420 irot=1,nrot  
           qmis(irot)=0.
           nmis(irot)=0
           do 400 ista=1,npsta
             p_b1= b1(1,irot)*p_a1(ista)
     &              +b1(2,irot)*p_a2(ista)
     &              +b1(3,irot)*p_a3(ista) 
             p_b3= b3(1,irot)*p_a1(ista)
     &              +b3(2,irot)*p_a2(ista)
     &              +b3(3,irot)*p_a3(ista) 
             if (sp_amp(ista).ne.0.) then
               p_proj1=p_a1(ista)-p_b3*b3(1,irot)
               p_proj2=p_a2(ista)-p_b3*b3(2,irot)
               p_proj3=p_a3(ista)-p_b3*b3(3,irot)
               plen=sqrt(p_proj1*p_proj1+p_proj2*p_proj2+
     &                    p_proj3*p_proj3)
               p_proj1=p_proj1/plen
               p_proj2=p_proj2/plen
               p_proj3=p_proj3/plen
               pp_b1=b1(1,irot)*p_proj1+b1(2,irot)*p_proj2
     &                +b1(3,irot)*p_proj3
               pp_b2=b2(1,irot)*p_proj1+b2(2,irot)*p_proj2
     &              +b2(3,irot)*p_proj3
               i=nint((p_b3+1.)/astep)+1
               theta=thetable(i)
               i=nint((pp_b2+1.)/astep)+1
               j=nint((pp_b1+1.)/astep)+1
               phi=phitable(i,j)
               i=nint(phi/(pi*astep))+1
               if (i.gt.2*ntab) i=1
               j=nint(theta/(pi*astep))+1
               if (j.gt.ntab) j=1
               p_amp=amptable(1,j,i)
               s_amp=amptable(2,j,i)
               if (p_amp.eq.0.0) then
                 sp_ratio=4.0
               else if (s_amp.eq.0.0) then
                 sp_ratio=-2.0
               else
                 sp_ratio=log10(4.9*s_amp/p_amp)
               end if
               qmis(irot)=qmis(irot)+abs(sp_amp(ista)-sp_ratio)

             end if
             if (p_pol(ista).ne.0) then
               prod=p_b1*p_b3
               ipol=-1
               if (prod.gt.0.) ipol=1 
               if (ipol.ne.p_pol(ista)) then
                  nmis(irot)=nmis(irot)+1                       
               end if
             end if
400         continue
            if (nmis(irot).lt.nmis0min) nmis0min=nmis(irot)
            if (qmis(irot).lt.qmis0min) qmis0min=qmis(irot)
420      continue
 
         nmis0max=ntotal
         if (nmis0max.lt.nmis0min+nextra) then
            nmis0max=nmis0min+nextra
         end if
         qmis0max=qtotal
         if (qmis0max.lt.qmis0min+qextra) then
            qmis0max=qmis0min+qextra
         end if

c loop over rotations - find those meeting fit criteria
425      nadd=0
         do irot=1,nrot        
            if ((nmis(irot).le.nmis0max).and. 
     &            (qmis(irot).le.qmis0max)) then
              irotgood(irot)=1
              nadd=nadd+1
            end if
         end do
       
c         print *,im,nmis0max,qmis0max,nadd

         if (nadd.eq.0) then  ! if there are none that meet criteria
           qmis0min=1.0e5     ! loosen the amplitude criteria
           do irot=1,nrot        
             if ((nmis(irot).le.nmis0max).and.
     &           (qmis(irot).lt.qmis0min)) then
               qmis0min=qmis(irot)
             end if
           end do
           qmis0max=qtotal
           if (qmis0max.lt.qmis0min+qextra) then
              qmis0max=qmis0min+qextra
           end if
           goto 425
         end if

430     continue

        nfault=0
        do irot=1,nrot
          if (irotgood(irot).gt.0) then
            nfault=nfault+1
            irotgood2(nfault)=irot
          end if
        end do

c  Select output solutions  
        nf=0      
        if (nfault.le.maxout) then
          do i=1,nfault
            irot=irotgood2(i)
            nf=nf+1
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
            end do
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
          end do
        else
          do 441 i=1,99999
            fran=rand(0)
            iscr=nint(fran*float(nfault)+0.5)
            if (iscr.lt.1) iscr=1
            if (iscr.gt.nfault) iscr=nfault
            if (irotgood2(iscr).le.0) goto 441
            irot=irotgood2(iscr)
            irotgood2(iscr)=-1
            nf=nf+1
            faultnorm(1)=b3(1,irot)
            faultnorm(2)=b3(2,irot)
            faultnorm(3)=b3(3,irot)
            slip(1)=b1(1,irot)
            slip(2)=b1(2,irot)
            slip(3)=b1(3,irot)
            do m=1,3
              faults(m,nf)=faultnorm(m)
              slips(m,nf)=slip(m)
            end do
            call FPCOOR(s1,d1,r1,faultnorm,slip,2)
            strike(nf)=s1
            dip(nf)=d1
            rake(nf)=r1
            if (nf.eq.maxout) go to 445
441       continue
445       continue
        end if
         
      return 
      end

c ------------------------------------------------------------------- c
      

c subroutine GET_MISF_AMP finds the percent of misfit polarities and the
c                         average S/P ratio misfit for a given mechanism  
c    Inputs:    npol   = number of polarity observations
c               p_azi_mc(npol) = azimuths
c               p_the_mc(npol) = takeoff angles
c               sp_ratio(npol) = S/P ratio
c               p_pol(npol)  = polarity observations
c               str_avg,dip_avg,rak_avg = mechanism
c    Outputs:   mfrac = weighted fraction misfit polarities
c               mavg = average S/P misfit (log10)
c               stdr = station distribution ratio

      subroutine GET_MISF_AMP(npol,p_azi_mc,p_the_mc,sp_ratio,p_pol,
     &          str_avg,dip_avg,rak_avg,mfrac,mavg,stdr)
cf2py intent(in)  npol
cf2py intent(in)  p_azi_mc
cf2py intent(in)  p_the_mc
cf2py intent(in)  sp_ratio
cf2py intent(in)  p_pol
cf2py intent(in)  str_avg
cf2py intent(in)  dip_avg
cf2py intent(in)  rak_avg
cf2py intent(out)  mfrac
cf2py intent(out)  mavg
cf2py intent(out)  stdr

      dimension p_azi_mc(npol),p_the_mc(npol)
      real str_avg,dip_avg,rak_avg,M(3,3),a(3),b(3),sp_ratio(npol)
      real strike,dip,rake,mfrac,mavg,qcount,azi,toff,pol,wt,wo
      integer k,npol,p_pol(npol)
      real bb1(3),bb2(3),bb3(3)
      
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
      stdr=0.
      scount=0.
      mavg=0.
      acount=0.
      
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
          p_proj1=p_proj1/plen
          p_proj2=p_proj2/plen
          p_proj3=p_proj3/plen
          pp_b1=bb1(1)*p_proj1+bb1(2)*p_proj2
     &              +bb1(3)*p_proj3
          pp_b2=bb2(1)*p_proj1+bb2(2)*p_proj2
     &              +bb2(3)*p_proj3
          phi=atan2(pp_b2,pp_b1)
          theta=acos(p_b3)
          p_amp=abs(sin(2*theta)*cos(phi))     
          wt=sqrt(p_amp)
          if (p_pol(k).ne.0) then
            azi=rad*p_azi_mc(k)
            toff=rad*p_the_mc(k)        
            a(1)=sin(toff)*cos(azi)
            a(2)=sin(toff)*sin(azi)
            a(3)=-cos(toff)
            do 615 in=1,3
              b(in)=0
              do 610 jn=1,3
                 b(in)=b(in)+M(in,jn)*a(jn)
610           continue
615         continue
            if ((a(1)*b(1)+a(2)*b(2)+a(3)*b(3)).lt.0) then
              pol=-1
            else
             pol=1
            end if
            if ((pol*p_pol(k)).lt.0) then
              mfrac=mfrac+wt
            end if
            qcount=qcount+wt
            stdr=stdr+wt
            scount=scount+1.0
          end if
          if (sp_ratio(k).ne.0.) then
            s1=cos(2*theta)*cos(phi)  
            s2=-cos(theta)*sin(phi)
            s_amp=sqrt(s1*s1+s2*s2)
            sp_rat=log10(4.9*s_amp/p_amp)
            mavg=mavg+abs(sp_ratio(k)-sp_rat)
            acount=acount+1.0
            stdr=stdr+wt
            scount=scount+1.0
          end if
600    continue
       mfrac=mfrac/qcount
       if (qcount.eq.0.0) mfrac=0.0
       mavg=mavg/acount
       if (acount.eq.0.0) mavg=0.0
       stdr=stdr/scount
       if (scount.eq.0.0) stdr=0.0
       
       return 
       end
