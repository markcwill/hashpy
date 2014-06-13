c MK_TABLES creates tables of takeoff angles given 1D velocity models.

c   output:
c     ntab - number of tables (max nindex)
c
c   you are prompted for the names of the 1D velocity model files,
c   velocity file format (free format):
c     depth(km) P_velocity(km/s)
      
      subroutine MK_TABLE(ntab)

      include 'vel.inc'
      real table(nx0,nd0,nindex),delttab(nx0),deptab(nd0)
      integer ndel,ndep

cf2py intent(in,out) ntab

      common /angtable/ table,delttab,deptab,ndel,ndep
c  common block:
c    table(nx0,nd0,nindex)  =  takeoff angle table
c        delttab(nx0)  =  list of ranges for tables
c         deptab(nd0)  =  list of source depths for tables
c           ndel       =  number of distance points in table
c           ndep       =  number of source depths in table

      parameter (nray0=10001)
      real z(1000),alpha(1000),slow(1000)
      character vmodel*100
      real deltab(nray0),tttab(nray0),ptab(nray0),tt(nray0,nd0)
      real depxcor(nray0,nd0),depucor(nray0,nd0),deptcor(nray0,nd0)
      real xsave(20000),tsave(20000),psave(20000),usave(20000)
c
      degrad=180./3.14159265

      if (ntab.ne.1) then
        print *,'Enter number of velocity models (max ',nindex,')'
        read *,ntab
      end if
      
      do 300 itab=1,ntab
      
      print *,'Enter file name for velocity model ',itab
      read (*,'(a)') vmodel

c set up table
      qtempdep2=dep2+dep3/20.
      ndep=int((qtempdep2-dep1)/dep3)+1
      do idep=1,ndep
         dep=dep1+dep3*real(idep-1)
         deptab(idep)=dep
      end do

c read velocity model
      open (7,file=vmodel,status='old')
      do i=1,1000
         read (7,*,end=30) z(i),alpha(i)
      end do
      print *,'***1000 point maximum exceeded in model'
30    close (7)
38    z(i)=z(i-1)           
      alpha(i)=alpha(i-1)
      npts=i
      npts_old=npts
      do i=npts_old,2,-1
        do idep=ndep,1,-1
          if ((z(i-1).le.(deptab(idep)-0.1)).and.
     &        (z(i).ge.(deptab(idep)+0.1))) then
            npts=npts+1
            do j=npts,i+1,-1
              z(j)=z(j-1)
              alpha(j)=alpha(j-1)
            end do
            z(i)=deptab(idep)
            frac=(z(i)-z(i-1))/(z(i+1)-z(i-1))
            alpha(i)=alpha(i-1)+frac*(alpha(i+1)-alpha(i-1))
          end if
        end do
      end do
      do i=1,npts
         slow(i)=1./alpha(i)
      end do
      pmax=slow(1)
      plongcut=slow(npts)
      pstep=(pmax-pmin)/float(nump)


c do P-wave ray tracing
       npmax=int((pmax+pstep/2.-pmin)/pstep)+1
       do 200 np=1,npmax
         p=pmin+pstep*real(np-1)
         ptab(np)=p
         x=0.
         t=0.
         imth=3
         do 70 idep=1,ndep
            if (deptab(idep).eq.0.) then
               depxcor(np,idep)=0.
               deptcor(np,idep)=0.
               depucor(np,idep)=slow(1)
            else
               depxcor(np,idep)=-999.
               deptcor(np,idep)=-999.
               depucor(np,idep)=-999.
            end if
70       continue
         do 100 i=1,npts-1
           if (z(i).ge.9999) then 
              deltab(np)=-999.
              tttab(np)=-999.
              go to 200
           end if
           h=z(i+1)-z(i)
           if (h.eq.0.) go to 100          !skip if interface
           call LAYERTRACE(p,h,slow(i),slow(i+1),imth,dx,dt,irtr)
           x=x+dx
           t=t+dt
           if (irtr.eq.0.or.irtr.eq.2) go to 105   !ray has turned
           xdeg=x                       ! actually in km 
           tmin=t                       ! actually in s
           do 80 idep=1,ndep
            if (abs(z(i+1)-deptab(idep)).lt.0.1) then
               depxcor(np,idep)=xdeg
               deptcor(np,idep)=tmin
               depucor(np,idep)=slow(i+1)            
            end if
80         continue
100      continue
105      xdeg=2.*x                     ! actually in km 
         tmin=2.*t                     ! actually in s
110      deltab(np)=xdeg
         tttab(np)=tmin
200   continue    !  end loop on ray parameter p

c create table
      do 250 idep=1,ndep
         icount=0
         xold=-999.
         if (deptab(idep).eq.0.) then
            i2=np
            go to 223
         end if
         do 220 i=1,np                    !upgoing rays from source
            x2=depxcor(i,idep)
            if (x2.eq.-999.) go to 221
            if (x2.le.xold) go to 221     !stop when heads inward
            t2=deptcor(i,idep)
            icount=icount+1
            xsave(icount)=x2
            tsave(icount)=t2
            psave(icount)=-ptab(i)
            usave(icount)=depucor(i,idep)
            xold=x2
220      continue
221      continue
         i2=i-1
223      do 225 i=i2,1,-1                 !downgoing rays from source
            if (depxcor(i,idep).eq.-999.) go to 225
            if (deltab(i).eq.-999.) go to 225
            x2=deltab(i)-depxcor(i,idep)
            t2=tttab(i)-deptcor(i,idep)
            icount=icount+1
            xsave(icount)=x2
            tsave(icount)=t2
            psave(icount)=ptab(i)
            usave(icount)=depucor(i,idep)
            xold=x2
225      continue
226      ncount=icount

         ndel=int((del2-del1)/del3)+1
         do 240 idel=1,ndel
            del=del1+del3*real(idel-1)
            delttab(idel)=del
            tt(idel,idep)=999.
            do 230 i=2,ncount
               x1=xsave(i-1)
               x2=xsave(i)
               if (x1.gt.del.or.x2.lt.del) go to 230
               if (psave(i).gt.0..and.psave(i).lt.plongcut) go to 230
               frac=(del-x1)/(x2-x1)
               t1=tsave(i-1)+frac*(tsave(i)-tsave(i-1))
               if (t1.lt.tt(idel,idep)) then
                  tt(idel,idep)=t1
                  scr1=psave(i)/usave(i)
                  angle=asin(scr1)*degrad
                  if (angle.lt.0.) then
                     angle=-angle
                  else
                     angle=180.-angle
                  end if
                  table(idel,idep,itab)=angle
               end if
230         continue
240      continue


250   continue
      if (delttab(1).eq.0.) then
         do idep=1,ndep
            table(1,idep,itab)=0.        !straight up at zero range
         end do
      end if
300   continue

c      do idel=1,ndel
c        do idep=1,ndep
c          print *,idel,delttab(idel),idep,deptab(idep),
c     &     table(idel,idep,1)
c        end do
c      end do
      
999   return
      end


c ------------------------------------------------------------ c


c subroutine GET_TTS obtains the takeoff angle for a velocity model
c at a specified range and earthquake depth by interpolating
c from a table of takeoff angles.  
c    Inputs:    ip     =  index number for model (up to nindex)
c               del    =  range
c               qdep   =  earthquake depth
c    Returns:   tt     =  takeoff angle (degrees)
c               iflag  = -1 if outside depth range
c                      =  0 for interpolation
c                      =  1 for extrapolation in range
c
      subroutine GET_TTS(ip,del,qdep,tt,iflag)

      include 'vel.inc'
      real t(nx0,nd0,nindex),x(nx0),d(nd0)
      integer nx,nd

cf2py intent(in) ip
cf2py intent(in) del
cf2py intent(in) qdep
cf2py intent(out) tt
cf2py intent(out) iflag

      common /angtable/ t,x,d,nx,nd
c  common block:
c    t(nx0,nd0,nindex)  =  takeoff angle tables
c          x(nx0)  =  list of ranges for tables
c          d(nd0)  =  list of source depths for tables
c           nx     =  number of distance points in table
c           nd     =  number of source depths in table

c
c check if outside depth range
      if (qdep.lt.d(1).or.qdep.gt.d(nd0)) then
         iflag=-1
         tt=999
         print *,'*** event outside of velocity table depth range, 
     &           event depth=',qdep,' table range=',d(1),d(nd0)
         return
      end if
c first check to see if interpolation alone will work
      do 30 id=2,nd
         if (d(id).lt.qdep) go to 30
         id1=id-1
         id2=id
         go to 32
30    continue
      id1=nd-1
      id2=nd
32    do 35 ix=2,nx
         if (x(ix).lt.del) go to 35
         ix1=ix-1
         ix2=ix
         go to 37
35    continue
      ix1=nx-1
      ix2=nx
37    if (t(ix1,id1,ip).eq.0.) go to 50
      if (t(ix1,id2,ip).eq.0.) go to 50
      if (t(ix2,id1,ip).eq.0.) go to 50
      if (t(ix2,id2,ip).eq.0.) go to 50
      if (x(ix2).lt.del) go to 50
      iflag=0
      xfrac=(del-x(ix1))/(x(ix2)-x(ix1))
      t1=t(ix1,id1,ip)+xfrac*(t(ix2,id1,ip)-t(ix1,id1,ip))
      t2=t(ix1,id2,ip)+xfrac*(t(ix2,id2,ip)-t(ix1,id2,ip))
      dfrac=(qdep-d(id1))/(d(id2)-d(id1))
      tt=t1+dfrac*(t2-t1)
      return
c extrapolate to get tt
50    iflag=1
      xoffmin1=999.
      xoffmin2=999.
      ixbest1=999
      ixbest2=999
      do 60 ix=2,nx
         if (t(ix-1,id1,ip).eq.0) go to 55
         if (t(ix,id1,ip).eq.0) go to 55
         xoff=abs((x(ix-1)+x(ix))/2.-del)
         if (xoff.lt.xoffmin1) then
            xoffmin1=xoff
            ixbest1=ix
         end if
55       if (t(ix-1,id2,ip).eq.0) go to 60
         if (t(ix,id2,ip).eq.0) go to 60
         xoff=abs((x(ix-1)+x(ix))/2.-del)
         if (xoff.lt.xoffmin2) then
            xoffmin2=xoff
            ixbest2=ix
         end if
60    continue
      if (ixbest1.eq.999.or.ixbest2.eq.999) then
         iflag=-1
         tt=999
         return
      end if
      xfrac1=(del-x(ixbest1-1))/(x(ixbest1)-x(ixbest1-1))
      t1=t(ixbest1-1,id1,ip)
      t2=t(ixbest1,id1,ip)
      tt1=t1+xfrac1*(t2-t1)
      xfrac2=(del-x(ixbest2-1))/(x(ixbest2)-x(ixbest2-1))
      t1=t(ixbest2-1,id2,ip)
      t2=t(ixbest2,id2,ip)
      tt2=t1+xfrac2*(t2-t1)
      dfrac=(qdep-d(id1))/(d(id2)-d(id1))
      tt=tt1+dfrac*(tt2-tt1)
999   return
      end


c ------------------------------------------------------------ c

c LAYERTRACE calculates the travel time and range offset
c for ray tracing through a single layer.
c
c Input:    p     =  horizontal slowness
c           h     =  layer thickness
c           utop  =  slowness at top of layer
c           ubot  =  slowness at bottom of layer
c           imth  =  interpolation method
c                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)
c                         = 2,  v(z) = a - b*z
c                         = 3,  v(z) = a*exp(-b*z)
c
c Returns:  dx    =  range offset
c           dt    =  travel time
c           irtr  =  return code
c                 = -1, zero thickness layer
c                 =  0,  ray turned above layer
c                 =  1,  ray passed through layer
c                 =  2,  ray turned within layer, 1 segment counted
c
c Note:  This version does calculation in double precision,
c        but all i/o is still single precision
c
      subroutine LAYERTRACE(p1,h1,utop1,ubot1,imth,dx1,dt1,irtr)
      implicit real*8 (a-h,o-z)
      real*4 p1,h1,utop1,ubot1,dx1,dt1
      p=dble(p1)
      h=dble(h1)
      utop=dble(utop1)
      ubot=dble(ubot1)
c
      if (h.eq.0.) then      !check for zero thickness layer
         dx1=0.
         dt1=0.
         irtr=-1
         return         
      end if
c
      u=utop
      y=u-p
      if (y.le.0.) then   !complex vertical slowness
         dx1=0.
         dt1=0.
         irtr=0
         return
      end if
c
      q=y*(u+p)
      qs=dsqrt(q)
c
c special function needed for integral at top of layer
      if (imth.eq.2) then
         y=u+qs
         if (p.ne.0.) y=y/p
         qr=dlog(y)
      else if (imth.eq.3) then
         qr=atan2(qs,p)
      end if      
c
      if (imth.eq.1) then
          b=-(utop**2-ubot**2)/(2.*h)
      else if (imth.eq.2) then
          vtop=1./utop
          vbot=1./ubot
          b=-(vtop-vbot)/h
      else
          b=-dlog(ubot/utop)/h
      end if  
c
      if (b.eq.0.) then     !constant velocity layer
         b=1./h
         etau=qs
         ex=p/qs
         irtr=1
         go to 160
      end if
c
c integral at upper limit, 1/b factor omitted until end
      if (imth.eq.1) then
         etau=-q*qs/3.
         ex=-qs*p
      else if (imth.eq.2) then
         ex=qs/u                 !*** - in some versions (wrongly)
         etau=qr-ex
         if (p.ne.0.) ex=ex/p
      else
         etau=qs-p*qr
         ex=qr
      end if
c
c check lower limit to see if we have turning point
      u=ubot
      if (u.le.p) then   !if turning point,
         irtr=2          !then no contribution
         go to 160       !from bottom point
      end if 
      irtr=1
      q=(u-p)*(u+p)
      qs=dsqrt(q)
c
      if (imth.eq.1) then
         etau=etau+q*qs/3.
         ex=ex+qs*p
      else if (imth.eq.2) then
         y=u+qs
         z=qs/u
         etau=etau+z
         if (p.ne.0.) then
            y=y/p
            z=z/p
         end if
         qr=dlog(y)
         etau=etau-qr
         ex=ex-z
      else
         qr=atan2(qs,p)
         etau=etau-qs+p*qr
         ex=ex-qr
      end if      
c
160   dx=ex/b
      dtau=etau/b
      dt=dtau+p*dx     !convert tau to t
c
      dx1=sngl(dx)
      dt1=sngl(dt)
      return
      end
