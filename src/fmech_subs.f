c subroutine FOCALMC performs grid search to find acceptable focal mechanisms,
c                    for multiple trials of ray azimuths and takeoff angles.
c                    Acceptable mechanisms are those with less than "ntotal"
c                    misfit polarities, or the minimum plus "nextra" if this
c                    is greater.
c
c  Inputs:  
c           p_azi_mc(npsta,nmc)  =  azimuth to station from event (deg. E of N)
c           p_the_mc(npsta,nmc)  =  takeoff angle (from vert, up=0, <90 upgoing, >90 downgoing)
c           p_pol(npsta)  =  first motion, 1=up, -1=down
c           p_qual(npsta) =  quality, 0=impulsive, 1=emergent
c           npsta  =  number of first motions
c           nmc    =  number of trials
c           dang   =  desired angle spacing for grid search
c           maxout =  maximum number of fault planes to return:
c                     if more are found, a random selection will be returned
c           nextra =  number of additional misfits allowed above minimum
c           ntotal =  total number of allowed misfits
c  Outputs: 
c           nf     =  number of fault planes found
c           strike(min(maxout,nf)) = strike
c           dip(min(maxout,nf))    = dip
c           rake(min(maxout,nf))   = rake
c           faults(3,min(maxout,nf)) = fault normal vector
c           slips(3,min(maxout,nf))  = slip vector
c
c
      subroutine FOCALMC(p_azi_mc,p_the_mc,p_pol,p_qual,npsta,nmc,
     &    dang,maxout,nextra,ntotal,nf,strike,dip,
     &    rake,faults,slips)
     
      include 'param.inc'
      include 'rot.inc'
     
c input and output arrays
      dimension p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0)
      integer p_pol(npsta),p_qual(npsta)
      real p_a1(npick0),p_a2(npick0),p_a3(npick0)
      real faultnorm(3),slip(3),faults(3,nmax0),slips(3,nmax0)
      real strike(nmax0),dip(nmax0),rake(nmax0)
      save dangold,nrot,b1,b2,b3
cf2py intent(in)  p_azi_mc
cf2py intent(in)  p_the_mc
cf2py intent(in)  p_pol
cf2py intent(in)  p_qual
cf2py intent(in)  npsta
cf2py intent(in)  nmc
cf2py intent(in)  dang
cf2py intent(in)  maxout
cf2py intent(in)  nextra
cf2py intent(in)  ntotal
cf2py intent(out) nf
cf2py intent(out) strike
cf2py intent(out) dip
cf2py intent(out) rake
cf2py intent(out) faults
cf2py intent(out) slips

c coordinate transformation arrays
      real b1(3,ncoor),bb1(3)
      real b2(3,ncoor),bb2(3)
      real b3(3,ncoor),bb3(3)

c fit arrays
      integer fit(2,ncoor),nmiss01min(0:npick0)
      integer irotgood(ncoor),irotgood2(ncoor)
      real fran

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
         nmiss0min=999
         nmissmin=999
         do i=0,npick0
           nmiss01min(i)=999
         end do
         do 420 irot=1,nrot  
            nmiss=0
            nmiss0=0           
            do 400 ista=1,npsta
               p_b1= b1(1,irot)*p_a1(ista)
     &              +b1(2,irot)*p_a2(ista)
     &              +b1(3,irot)*p_a3(ista) 
               p_b3= b3(1,irot)*p_a1(ista)
     &              +b3(2,irot)*p_a2(ista)
     &              +b3(3,irot)*p_a3(ista) 
               prod=p_b1*p_b3
               ipol=-1
               if (prod.gt.0.) ipol=1    ! predicted polarization
               if (ipol.ne.p_pol(ista)) then
                  nmiss=nmiss+1                       
                  if (p_qual(ista).eq.0) nmiss0=nmiss0+1
               end if
400         continue
            fit(1,irot)=nmiss0           ! misfit impulsive polarities
            fit(2,irot)=nmiss            ! total misfit polarities
            if (nmiss0.lt.nmiss0min) nmiss0min=nmiss0
            if (nmiss.lt.nmissmin) nmissmin=nmiss
            if (nmiss.lt.nmiss01min(nmiss0)) then
               nmiss01min(nmiss0)=nmiss
            end if
420      continue

c choose fit criteria
         if (nmiss0min.eq.0) then 
            nmiss0max=ntotal
            nmissmax=ntotal
         else
            nmiss0max=ntotal
            nmissmax=npsta
         end if
         if (nmiss0max.lt.nmiss0min+nextra) then
            nmiss0max=nmiss0min+nextra
         end if
         if (nmissmax.lt.nmiss01min(nmiss0min)+nextra) then
            nmissmax=nmiss01min(nmiss0min)+nextra
         end if

c loop over rotations - find those meeting fit criteria
         do 440 irot=1,nrot        
            nmiss0=fit(1,irot)
            nmiss=fit(2,irot)
            if ((nmiss0.le.nmiss0max).and.(nmiss.le.nmissmax)) then
              irotgood(irot)=1
            end if
440     continue

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
