c Sample main driver program for using the focal mechanism subroutines.  
c Takeoff angle uncertainty is given in the polarity input file, otherwise similar
c input format to FPFIT.
c
c INPUT FILE: "fpfile" FPFIT-like format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
c
c event line:
c     columns  format   value
c     -------------------------
c     1-10       5i2    origin time, year, month, day, hour, minute
c     11-14      f4.2   origin time, seconds
c     15-16      i2     latitude, degrees
c     17         a1     latitude, 'S'=south
c     18-21      f4.2   latitude, minutes
c     22-24      i3     longitude, degrees
c     25         a1     longitude, 'E'=east
c     26-29      f4.2   longitude, minutes
c     30-34      f5.2   depth, km
c     35-36      f2.1   magnitude
c     81-88      2f4.2  horizontal and vertical uncertainty, km
c     123-138    a16    event ID
c
c polarity line:
c     columns  format   value
c     -------------------------
c     1-4        a4     station name
c     7          a1     polarity : U, u, +, D, d, or -
c     8          i1     quality: 0=good quality, 1=lower quality, etc
c     59-62      f4.1   source-station distance (km)
c     66-68      i3     takeoff angle
c     79-81      i3     azimuth
c     83-85      i3     takeoff angle uncertainty  ** NOT in standard FPFIT files
c     87-89      i3     azimuth uncertainty        ** NOT in standard FPFIT files
c

      include 'param.inc'
      include 'rot.inc'

c variables for storing earthquake input information  
      character*16 icusp            ! event ID
      real qlat,qlon,qdep           ! location
      real qmag                     ! magnitude
      integer iyr,imon,idy,ihr,imn  ! origin time, year, month, day, hour, minute
      real qsec                     ! origin time, seconds
      real seh, sez                 ! location error, horizontal & vertical 
      real rms                      ! residual travel time error 
      real terr                     ! origin time error 
      character*1 evtype            ! event type
      character*1 magtype           ! magnitude type
      character*1 locqual           ! location quality
      character*1 cns,cew           ! north/south and east/west codes
c
c variables for polarity information, input to HASH subroutines
      character*4 sname(npick0)                        ! station name
      character*1 pickpol(npick0)                      ! polarity pick : U, u, +, D, d, or -
      integer p_pol(npick0),p_qual(npick0),spol        ! polarity pick (+/-1),  quality (0/1), and reversal (+/-1)
      real p_azi_mc(npick0,nmc0),p_the_mc(npick0,nmc0) ! azimuth and takeoff angle for each trial
      integer nmc                                      ! number of trails with different azimuth and take-off angles
      integer npol                                     ! number of polarity readings
c
c variables for set of acceptable mechanisms, output of HASH subroutines
      integer nmult                                    ! number of solutions (1 if no mulitples)
      real str_avg(5),dip_avg(5),rak_avg(5)            ! solution(s)
      integer nout2                                    ! number of acceptable mechanisms returned
      real f1norm(3,nmax0),f2norm(3,nmax0)             ! normal vectors to the fault planes
      real strike2(nmax0),dip2(nmax0),rake2(nmax0)     ! strike, dip and rake
      real var_est(2,5),var_avg(5)                     ! variance of each plane, average
      real mfrac(5),stdr(5)                            ! fraction misfit polarities, station distribution ratio
      real prob(5)                                     ! probability true mechanism is "close" to preferred solution(s)
      character*1 qual(5),mflag                        ! solution quality rating, multiple flag
c
c control parameters
      integer npolmin,max_agap,max_pgap                ! minimum polarities, maximum azimuthal & "plungal" gap
      real delmax                                      ! maximum station distance
      real dang,dang2                                  ! grid angle to sample focal sphere
      integer maxout                                   ! max number of acceptable mechanisms output
      real badfrac                                     ! assumed rate of polarity error (fraction)
      real cangle                                      ! definition of "close" == 45 degrees
c
c file names
      character*100 outfile1,outfile2,plfile,fpfile

      degrad=180./3.1415927
      rad=1./degrad
      
      print *,'Enter station polarity reversal file'       
      read (*,'(a)') plfile

      print *,'Enter name of input file (FPFIT-like format)'       
      read (*,'(a)') fpfile
      open (12,file=fpfile,status='old')

      print *,'Enter output file name for focal mechanisms'
      read (*,'(a)') outfile1
      open (13,file=outfile1)

      print *,'Enter output file name for acceptable planes'
      read (*,'(a)') outfile2
      open (11,file=outfile2)

      print *,'Enter mininum number of polarities (e.g., 8)'
      read *,npolmin

      print *,'Enter maximum azimuthal gap (e.g., 90)'
      read *,max_agap

      print *,'Enter maximum takeoff angle gap (e.g., 60)'
      read *,max_pgap

      print *,'Enter grid angle for focal mech search, in degrees 
     &  (min ',dang0,')'
      read *,dang
      dang2=max(dang0,dang) ! don't do finer than dang0
      
      print *,'Enter number of trials (e.g., 30)'
      read *,nmc

      print *,'Enter maxout for focal mech. output (e.g., 500)'
      read *,maxout

      print *,'Enter fraction of picks presumed bad (e.g., 0.10)'
      read *,badfrac

      print *,'Enter maximum allowed source-station distance, 
     &         in km (e.g., 120)'
      read *,delmax

      print *,'Enter angle for computing mechanisms probability, 
     &         in degrees (e.g., 45)'
      read *,cangle
      
      print *,'Enter probability threshold for multiples (e.g., 0.1)'
      read *,prob_max
      
20    continue

c read in earthquake location, etc      ! FPFIT-like format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
      read (12,25,end=999) iyr,imon,idy,ihr,imn,qsec,ilatd,cns,qlatm,
     &                ilond,cew,qlonm,qdep,qmag,seh,sez,icusp
25    format (5i2,f4.2,i2,a1,f4.2,i3,a1,f4.2,f5.2,f2.1,44x,2f4.2,
     &                34x,a16)
      if (iyr.lt.50) then     ! fix Y2K problem
        iyr=iyr+2000
      else
        iyr=iyr+1900
      end if        
      qlat=real(ilatd)+(qlatm/60.0)
      if (cns.eq.'S') then
        qlat=-qlat
      end if
      qlon=-(real(ilond)+(qlonm/60.0))
      if (cew.eq.'E') then
        qlon=-qlon
      end if
      if (sez.eq.0.) sez=1.   ! set parameters not given in input file
      terr=-9               
      rms=-9
      nppick=-9
      nspick=-9
      evtype='L'
      magtype='X'
      locqual='X'
            
c read in polarities                 ! FPFIT-like format - ** YOU MAY NEED TO CHANGE THE INPUT FORMAT **
      k=1
30    continue
        read (12,35) sname(k),pickpol(k),p_qual(k),
     &           qdist,ith,iaz,isthe,isazi
35      format (a4,2x,a1,i1,50x,f4.1,i3,10x,i3,1x,i3,1x,i3)
        if (sname(k).eq.'    ')  goto 40 ! end of data for this event
        if (qdist.gt.delmax) goto 30
        if (p_qual(k).gt.1) goto 30
        if (pickpol(k).eq.'U'.or.
     &                    pickpol(k).eq.'u'.or.
     &                    pickpol(k).eq.'+') then
          p_pol(k)=1
        else if (pickpol(k).eq.'D'.or.
     &                    pickpol(k).eq.'d'.or.
     &                    pickpol(k).eq.'-') then
          p_pol(k)=-1
        else
          goto 30
        end if
        call CHECK_POL(plfile,sname(k),iyr,imon,idy,ihr,spol)  
                          ! SCSN station polarity reversal information - ** YOU MAY NEED TO USE YOUR OWN SUBROUTINE **        
        p_pol(k)=p_pol(k)*spol
        qazi=real(iaz)
        qthe=180.-real(ith)
        sazi=real(isazi)
        sthe=real(isthe)
        if (qazi.lt.0.) qazi=qazi+360.
        p_azi_mc(k,1)=qazi      ! set azimuth and takeoff angle for each trial, given uncertainty
        p_the_mc(k,1)=qthe
        do 105 nm=2,nmc
          call RAN_NORM(val)
          p_azi_mc(k,nm)=qazi+sazi*val
          call RAN_NORM(val)
          p_the_mc(k,nm)=qthe+sthe*val   
105     continue
        k=k+1
      goto 30
40    continue
      npol=k-1

cc view polarity data
c      do 250 k=1,npol 
c          print *,k,' ',sname(k),p_azi_mc(k,1),p_the_mc(k,1),
c     &         p_azi_mc(k,2),p_the_mc(k,2),p_pol(k)
c250   continue

c stop if there aren't enough polarities
      print *,'cid = ',icusp,'  npol = ',npol
      if (npol.lt.npolmin) then
        str_avg(1)=999
        dip_avg(1)=99
        rak_avg(1)=999
        var_est(1,1)=99
        var_est(2,1)=99
        mfrac(1)=0.99
        qual(1)='F'
        prob(1)=0.0
        nout1=0
        nout2=0
        nmult=0
        goto 400
      end if

c determine maximum azimuthal and takeoff gap in polarity observations
c and stop if either gap is too big
      call GET_GAP(npol,p_azi_mc,p_the_mc,magap,mpgap)
      if ((magap.gt.max_agap).or.(mpgap.gt.max_pgap)) then
        str_avg(1)=999
        dip_avg(1)=99
        rak_avg(1)=999
        var_est(1,1)=99
        var_est(2,1)=99
        mfrac(1)=0.99
        qual(1)='E'
        prob(1)=0.0
        nout1=0
        nout2=0
        nmult=0
        goto 400
      end if

c determine maximum acceptable number misfit polarities
      nmismax=max(nint(npol*badfrac),2)                    
      nextra=max(nint(npol*badfrac*0.5),2)  

c find the set of acceptable focal mechanisms for all trials            
      call FOCALMC(p_azi_mc,p_the_mc,p_pol,p_qual,npol,nmc,
     &       dang2,nmax0,nextra,nmismax,nf2,strike2,dip2,
     &       rake2,f1norm,f2norm)
      nout2=min(nmax0,nf2)  ! number mechs returned from sub
      nout1=min(maxout,nf2)  ! number mechs to return
 
c find the probable mechanism from the set of acceptable solutions          
      call MECH_PROB(nout2,f1norm,f2norm,cangle,prob_max,nmult,
     &          str_avg,dip_avg,rak_avg,prob,var_est)           

      do 390 imult=1,nmult
      
      var_avg(imult)=(var_est(1,imult)+var_est(2,imult))/2.
      print *,'cid = ',icusp,imult,'  mech = ',
     &          str_avg(imult),dip_avg(imult),rak_avg(imult)

c find misfit for prefered solution
      call GET_MISF(npol,p_azi_mc,p_the_mc,p_pol,p_qual,str_avg(imult),
     &       dip_avg(imult),rak_avg(imult),mfrac(imult),stdr(imult))
      
c solution quality rating  ** YOU MAY WISH TO DEVELOP YOUR OWN QUALITY RATING SYSTEM **
      if ((prob(imult).gt.0.8).and.(var_avg(imult).le.25).and.
     &     (mfrac(imult).le.0.15).and.(stdr(imult).ge.0.5)) then
        qual(imult)='A'
      else if ((prob(imult).gt.0.6).and.(var_avg(imult).le.35).and.
     &  (mfrac(imult).le.0.2).and.(stdr(imult).ge.0.4)) then
        qual(imult)='B'
      else if ((prob(imult).gt.0.5).and.(var_avg(imult).le.45).and.
     &  (mfrac(imult).le.0.3).and.(stdr(imult).ge.0.3)) then
        qual(imult)='C'
      else
        qual(imult)='D'
      end if

390   continue

400   continue
       
c output prefered mechanism  ** YOU MAY WISH TO CHANGE THE OUTPUT FORMAT **
      if (nmult.gt.1) then
        mflag='*'
      else
        mflag=' '
      end if
      do i=1,nmult
      write (13,411) icusp,iyr,imon,idy,ihr,imn,qsec,evtype,
     &   qmag,magtype,qlat,qlon,qdep,locqual,rms,seh,sez,terr,
     &   nppick+nspick,nppick,nspick,
     &   nint(str_avg(i)),nint(dip_avg(i)),nint(rak_avg(i)),
     &   nint(var_est(1,i)),nint(var_est(2,i)),npol,nint(mfrac(i)*100.),
     &   qual(i),nint(100*prob(i)),nint(100*stdr(i)),mflag
      end do
411   format(a16,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,1x,a1,1x,
     &  f5.3,1x,a1,1x,f9.5,1x,f10.5,1x,f7.3,1x,a1,1x,f7.3,1x,f7.3,
     &  1x,f7.3,1x,f7.3,3x,i4,1x,i4,1x,i4,1x,i4,1x,i3,1x,i4,3x,i2,
     &  1x,i2,1x,i3,1x,i2,1x,a1,1x,i3,1x,i2,1x,a1)

c output set of acceptable mechanisms  ** YOU MAY WISH TO CHANGE THE OUTPUT FORMAT **
      write (11,412) iyr,imon,idy,ihr,imn,qsec,qmag,
     &   qlat,qlon,qdep,sez,seh,npol,nout2,icusp,
     &   str_avg(1),dip_avg(1),rak_avg(1),var_est(1,1),var_est(2,1),
     &   mfrac(1),qual(1),prob(1),stdr(1)
412   format (i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,2x,f3.1,1x,f9.4,1x,
     &   f10.4,1x,f6.2,1x,f8.4,1x,f8.4,1x,i5,1x,i5,1x,a16,1x,f7.1,1x,
     &   f6.1,1x,f7.1,1x,f6.1,1x,f6.1,1x,f7.3,2x,a1,1x,f7.3,1x,f4.2)
      do 500 ic=1,nout1
        write (11,550) strike2(ic),dip2(ic),rake2(ic),f1norm(1,ic),
     &      f1norm(2,ic),f1norm(3,ic),f2norm(1,ic),f2norm(2,ic),
     &      f2norm(3,ic)
500   continue
550   format (5x,3f9.2,6f9.4)

800   continue
      goto 20

999   close (11)
      close (12)
      close (13)
      stop
      end
