c MK_TABLES creates tables of takeoff angles given 1D velocity models.
c
c   input: 
c     ind  -  index of the table to add
c   vmodel -  filename of the velocity model file
c   output:
c     ntab - number of tables (max nindex)
c
c   add a velocity model to the table,
c   velocity file format (free format):
c     depth(km) P_velocity(km/s)
      
      subroutine MK_TABLE_ADD(ind,ntab,vmodel)

      include 'vel.inc'
      real table(nx0,nd0,nindex),delttab(nx0),deptab(nd0)
      integer ndel,ndep

cf2py intent(in) ind
cf2py intent(in) vmodel
cf2py intent(out) ntab

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

c      if (ntab.ne.1) then
c        print *,'Enter number of velocity models (max ',nindex,')'
c        read *,ntab
c      end if
      ntab=ind
      
      do 300 itab=ind,ntab
      
c      print *,'Enter file name for velocity model ',itab
c      read (*,'(a)') vmodel

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
