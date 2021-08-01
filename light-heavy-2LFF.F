#ifdef SMH
      subroutine interpolate_F1F2_LH_HL(t,tminin,tmaxin,m,mh,f1lh,f2lh,f1hl,f2hl)
      implicit none
      double precision t,m,mh,f1lh,f2lh,tminin,tmaxin,f1hl,f2hl
      double complex f1syed,f2syed
      real startt,endt
      integer nint,k,ifirst,i,ir
      parameter (nint=2000)
      double precision tmincut,tmaxcut,tcur,tmin,tmax,dt,tr
      common/tmintmax/tmincut,tmaxcut
      double precision tv(0:nint+1)
      double precision f1lhv(0:nint+1),f2lhv(0:nint+1)
      double precision f1hlv(0:nint+1),f2hlv(0:nint+1)
      common/f1f2v/f1lhv,f2lhv,f1hlv,f2hlv,tv,tmin,tmax,dt,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         tmin = 1.1d0*tminin
         tmax = 0.9d0*tmaxin
         dt = (tmax-tmin)/nint
         tr = tmin
         print*,'-- Populating arrays for light-to-heavy bubbles '//
     .        '(in quad precision) --'
         call cpu_time(startt)
         do k=1,nint
            ir = 0
            if (k.eq.1) ir = 1
            call twoloop_light_heavy_Syed_quad(tr,m,mh,ir,f1syed,f2syed)
            tv(k)    = tr
            f1lhv(k) = dreal(f1syed)
            f2lhv(k) = dreal(f2syed)            
            tr       = tmin + k*dt
         enddo
         call cpu_time(endt)
         print*,'-- finished in ',endt-startt,'seconds --'

         print*,'-- Populating arrays for heavy-to-light bubbles '//
     .        '(in quad precision) --'
         call cpu_time(startt)
         tr = tmin
         do k=1,nint
            ir = 0
            if (k.eq.1) ir = 1
            call twoloop_light_heavy_Syed_quad(tr,mh,m,ir,f1syed,f2syed)
            f1hlv(k) = dreal(f1syed)
            f2hlv(k) = dreal(f2syed)            
            tr       = tmin + k*dt
         enddo
         call cpu_time(endt)
         print*,'-- finished in ',endt-startt,'seconds --'
         
         tv(0)  = tv(1) - dt

         f1lhv(0) = f1lhv(1)
         f2lhv(0) = f2lhv(1)
         tv(nint+1)  = tv(nint) + dt
         f1lhv(nint+1) = (f1lhv(nint-1)-f1lhv(nint-2)) +f1lhv(nint)
         f2lhv(nint+1) = (f2lhv(nint-1)-f2lhv(nint-2)) +f2lhv(nint)

         f1hlv(0) = f1hlv(1)
         f2hlv(0) = f2hlv(1)
         f1hlv(nint+1) = (f1hlv(nint-1)-f1hlv(nint-2)) +f1hlv(nint)
         f2hlv(nint+1) = (f2hlv(nint-1)-f2hlv(nint-2)) +f2hlv(nint)
         
         ifirst = 1
      endif

      k = int((t-tmin)/dt) + 1

c      tr = tmin
c      i = 0
c      do while (t.gt.tr)
c         i = i+1
c         tr = tmin + i*dt
c      enddo
      
      f1lh = (f1lhv(k+1)-f1lhv(k))/dt * (t-tv(k)) + f1lhv(k)
      f2lh = (f2lhv(k+1)-f2lhv(k))/dt * (t-tv(k)) + f2lhv(k)

      f1hl = (f1hlv(k+1)-f1hlv(k))/dt * (t-tv(k)) + f1hlv(k)
      f2hl = (f2hlv(k+1)-f2hlv(k))/dt * (t-tv(k)) + f2hlv(k)
      
      return
      end
************************************************************************************
      subroutine burgersF1timelike(burg,s,m2,mh2)
      implicit none
      double precision burg,m2,mh2,pi,z3,L,s,C
      parameter (z3=1.202056903159594285399738161511449990764986292d0)
      parameter (pi=4.d0*atan(1.d0))

      L = log(s/mh2)

*** m = mh      
      C = 383.d0/108.d0 - 11.d0/36.d0*pi*pi
      
*** mh diverso da m      
      C = -1.d0/3.d0*z3 - 19.d0/18.d0*pi*pi/6.d0 + 3355d0/1296d0

      burg = -1.d0/36.d0 * L**3 + 19.d0/72.d0*L**2
     .     + (pi*pi/18.d0-265d0/216d0)*L + C
      
      burg = burg * 4.d0

      return
      end
      
************************************************************************************      
      subroutine twoloop_light_heavy_Syed(t,mo,mho,ires,f1syed,f2syed)
      use handyg
*      implicit none      
      implicit double complex (Z)
      double complex f1syed,f2syed,f11,f12,f13,f14
      double complex f21, f22,f11ex,f12ex,f13ex,f14ex
      double complex f1comb,f2comb
      double precision z
      double precision m,mh,t,del,mo,mho
      double complex cone,x,s,xs,xvar,carg,xh,umxs
      parameter(cone=(1.d0,0.d0))
      double complex i_
      parameter(i_=(0.d0,1.d0))
      double precision log2,sq2
      double precision zeta2,zeta3
      double precision alpha,convfac,pi,ame,amu,egmin
      
      common/parameters/ame,amu,convfac,alpha,pi      
*
      double precision b5a1,b5a2,b5a3,b5a4,b6a3,b6a4
      double precision b5a1b,b5a2b,b5a3b,b5a4b,b6a3b,b6a4b
      integer ifirst,ico,i,ires
      double precision PL4unm,logmomh
      double precision pi2,pi4
      double complex Li22,Li32,logupz,logumz,li2momh,li2momh2
      double complex logA,logB,logC,atanhz,eps
      common/twolcmnlightheavy/log2,sq2,pi2,pi4,zeta2,zeta3,
     .     b5a1,b5a2,b5a3,b5a4,b6a3,b6a4,Li22,Li32,logmomh,
     .     li2momh,li2momh2,logupz,logumz,atanhz,PL4unm,ifirst,ico
      double complex GPL(1:882)
      common/GPLscommon/GPL,f11,f12,f13,f14,f21,f22,f1comb,f2comb
      data ifirst /0/
      double complex c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
      double complex c11,c12,c13,c14,c15,c16,c17,c18,c19,c20
      double complex c21,c22,c23,c24,c25,c26,c27

*********** for HPLOG
      integer n1,n2,nw,i1,i2,k1,k2,k3
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
      parameter(n2=1)
      parameter(nw=3)
      double complex Hc1,Hc2,Hc3,Hc4,li3A,li3B,li2A,li2B
      double precision  Hr1,Hr2,Hr3,Hr4,Hupx,Humx
      double precision  Hi1,Hi2,Hi3,Hi4
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2),
     .     Hc4(n1:n2,n1:n2,n1:n2,n1:n1)
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),
     .     Hr4(n1:n2,n1:n2,n1:n2,n1:n1)
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),
     .     Hi4(n1:n2,n1:n2,n1:n2,n1:n1)
***** 
*     quads
      complex*32 xsq,sq
      real*16 mq
*****      
      call clearcache
*************************      
      s  = t!/mo/mo
      m  = mo!1.d0!mo
      mh = mho!/mo
      z = m/mh

      xs = sqrt(1.d0-4.d0*m*m/s)
      umxs = 1.d0 - xs
      
c      xs = (1.d0+i_*2.d0*m/sqrt(-s))*(1.d0-i_*2.d0*m/sqrt(-s))
c      xs = sqrt(xs)
**************************      

c       print*,
c     &  (+ Sqrt(s**2 - 4*m**2*s) * ( s - 2*m**2 ))/m**4
c     &,(+ s**2 - 4*m**2*s)/m**4

c      print*,
c     & (+ Sqrt(s**2 - 4*m**2*s) * ( s - 2*m**2 ))/m**4
c     &,(+ s**2 - 4*m**2*s)/m**4
      
      if (ifirst.eq.0.or.ires.gt.0) then
         ico = 0
         log2= log(2.d0)
         logmomh = log(m/mh)
         
         sq2 = sqrt(2.d0)
         zeta2=pi*pi/6.d0
         zeta3=1.202056903159594285399738161511449990764986292d0         
         pi2= pi*pi
         pi4= pi2*pi2

         b5a1 =  2.77258872223978109212d0
         b5a2 =  1.36805607802364681993647982810813818502d0
         b5a3 =  1.13556770870411119047364593456472384651d0
         b5a4 =  1.0583097143302935d0
         b6a3 =  8.41439832211715975835164054594208875911d0
         b6a4 = -8.30035652771381970800065987017540898127d0
* more precise
         b5a1 = 2.77258872223978123766892848583270669d0
         b5a2 = 1.36805607802364717427642022798538871d0
         b5a3 = 1.13556770870411139095090943038957078d0
         b5a4 = 1.05830971433029308016425544951004003d0
         b6a3 = 8.41439832211715999779816713065087404d0
         b6a4 =-8.30035652771381950065162465679101891d0

         PL4unm = 0.517479061673899386330758161898862945622377475141d0
         
c         Li22 = myli2(cone*2.d0)
c         call hplog(2.d0,nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c         Li22 = Hc2(0,1)         
c         call hplog(2.d0,nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c         Li32 = Hc3(0,0,1)
         
         Li22=-G([inum(0.d0,0),inum(1.d0,0)],inum(2.d0,di0))
         Li32=-G([inum(0.d0,0),inum(0.d0,0),inum(1.,0)],inum(2.d0,di0))  

c     call hplog(z,nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c     li2momh  = Hc2(0,1) ! not ok for mh<m!
         li2momh=-G([inum(0.d0,0),inum(1.d0,0)],inum(z,di0))
c     call hplog(z*z,nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c     li2momh2 = Hc2(0,1)  ! not ok for mh<m!      
         li2momh2=-G([inum(0.d0,0),inum(1.d0,0)],inum(z*z,di0))
         
c         logupz  = log(cone*(1.d0+z))
c         logumz  = log(cone*(1.d0-z))
         logupz=G([inum(cone*0.d0,0)],inum(cone*(1.d0+z),di0))
         logumz=G([inum(cone*0.d0,0)],inum(cone*(1.d0-z),di0))

         atanhz = 0.5d0*(logupz-logumz) ! sostituito nel form direttamente
         
         include 'constgpl_defs.f'
         
         ifirst = 1
      endif

      ico = ico + 1

      x = sqrt(cone*(4.d0*m*m-s))*sqrt(cone*(4.d0*m*m-4.d0*mh*mh-s))
      x = -2.d0*m*m -x+s + sq2*sqrt(cone*(-s*(4.d0*m*m-2.d0*mh*mh+x-s)))
      x = -x/m/mh*0.5d0
      
      xvar =
     &  + m**(-1) * (  - 1.D0/2.D0*s*mh**(-1) + 1.D0/2.D0*sqrt( - s)*
     &    sqrt( - 4*mh**2 - s*xs**2)*mh**(-1)*xs - 1.D0/2.D0*sqrt(2*s*
     &    mh**2 + s**2*xs**2 - sqrt( - s)*sqrt( - 4*mh**2 - s*xs**2)*s*
     &    xs)*mh**(-1)*sq2 )
      xvar = xvar + z

cc      x = xvar

      eps = x - 1.d0
      
c      print*,eps
      
      include 'vargpl_defs.f'
c      include 'vargpl_defs_eqmasses.f'

c      call hplog(dreal( -1d0 + 2d0/(1d0 - xs**2) - 2d0/(1d0 - xs**2)*xs),
c     +     nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c      li2A = Hc2(0,1)
c      li3A = Hc3(0,0,1)

      carg = -1d0 + 2d0/(1d0 - xs**2) - 2d0/(1d0 - xs**2)*xs
      li2A = -G([inum(0.d0,0),inum(1.d0,0)],inum(carg,di0))
      li3A = -G([inum(0.d0,0),inum(0.d0,0),inum(1.d0,0)],inum(carg,di0))

c      mq = m
c      sq = s
c      xsq = sqrt(1.q0-4.q0*mq*mq/sq)
      
c$$$      print*,'--'
c$$$      print*,xs,s/m/m,xsq,sq/mq/mq
c$$$      print*,-0.125d0/( - 1.d0/(1 - xs**2) + 1.d0/2d0/(1 - xs**2)*xs
c$$$     .     + 1.d0/(1 - xs**2)/(1 - xs**2)
c$$$     .     - 1.d0/(1 - xs**2)/(1 - xs**2)*xs),
c$$$     .     -0.125q0/( - 1.q0/(1 - xsq**2) + 1.q0/2q0/(1 - xsq**2)*xsq
c$$$     .     + 1.q0/(1 - xsq**2)/(1 - xsq**2) - 1.q0/(1 - xsq**2)/(1 - xsq**2)*xsq)
c$$$      print*,-m*m*3.d0*1.d0/( - 6*s + 3*s*xs + 6d0/(1 - xs**2)*s
c$$$     .     - 6d0/(1 - xs**2)*s*xs),
c$$$     .     -mq*mq*3.q0/( - 6*sq + 3*sq*xsq + 6q0/(1 - xsq**2)*sq
c$$$     .     - 6q0/(1 - xsq**2)*sq*xsq)
c$$$
c$$$      print*,-m*m*3.d0/(3.d0/2.d0*s*( (s/m**2 - 2.d0)*(1.d0-xs) - 2.d0)
c$$$     .     ),-mq*mq*3.q0/(3.q0/2.q0*s*( (sq/mq**2 - 2.q0)*(1.q0-xsq) - 2.q0)
c$$$     .     ),-mq*mq*3.q0/(( - 6*sq
c$$$     .     + 3*sq*xsq + 6q0/(1 - xsq**2)*sq
c$$$     .     - 6q0/(1 - xsq**2)*sq*xsq))
      
c      call hplog(dreal(2d0 - 2d0/(1d0 - xs**2) + 2d0/(1d0 - xs**2)*xs),
c     +     nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c      li2B = Hc2(0,1)
c      li3B = Hc3(0,0,1)

!      carg = 2d0 - 2d0/(1d0 - xs**2) + 2d0/(1d0 - xs**2)*xs, i.e.
      carg = 1.d0 - carg
      li2B = -G([inum(0.d0,0),inum(1.d0,0)],inum(carg,di0))
      li3B = -G([inum(0.d0,0),inum(0.d0,0),inum(1.d0,0)],inum(carg,di0))

      logA = Log(cone*( - 1 + 2/(1 - xs**2) - 2/(1 - xs**2)*xs))
      logB = Log(cone*(2 - 2/(1 - xs**2) + 2/(1 - xs**2)*xs))
c not neede logA - logC = i_*pi      logC = Log(cone*(1 - 2/(1 - xs**2) + 2/(1 - xs**2)*xs))
c      ie logA = G([inum(0.d0,0)],inum(cone*(( - 1 + 2/(1 - xs**2) - 2/(1-xs**2)*xs)),di0))
c      ie logB = G([inum(0.d0,0)],inum(cone*(( + 1 - 2/(1 - xs**2) + 2/(1-xs**2)*xs)),di0))
ccc      print*,'!!!!REINCLUDE!!!!'
c      include 'f1comb_light_heavy.f'
c     include 'f1_onlyf11_expanded_in_z.f'
      
      include 'f1_light_heavy.f'
c      include 'Syed/gpl/F1eqmasses_eps.f'
c      f11 = 0.d0
c      print*,''
c      print*,f11
c      print*,f12
c      print*,f13
c      print*,f14
c      print*,f11 - f12 - f13 + f14
      f1syed = f11 - f12 - f13 + f14
      include 'f2_light_heavy.f'
c      include 'Syed/gpl/F2eqmasses_eps.f'
      f2syed = f21 - f22
      
c      print*,f2syed
c      print*,'** F1 **'
c      write(*,'(A,f25.18," + ",f25.18," * i")')' f1_1  = ',f11
c      write(*,'(A,f25.18," + ",f25.18," * i")')' f2_1  = ',f12
c      write(*,'(A,f25.18," + ",f25.18," * i")')' f3_1  = ',f13
c      write(*,'(A,f25.18," + ",f25.18," * i")')' f4_1  = ',f14
c      write(*,'(A,f25.18," + ",f25.18," * i")')
c     .     ' sum = f1_1 - f2_1 - f3_1  + f4_1  = ',f1syed
c      print*,'** F2 **'
c      write(*,'(A,f25.18," + ",f25.18," * i")')' f1_2  = ',f21
c      write(*,'(A,f25.18," + ",f25.18," * i")')' f2_2  = ',f22
c      write(*,'(A,f25.18," + ",f25.18," * i")')
c     .     ' sum = f2_1 - f2_2  = ',f2syed            
c     if (ico.eq.10) stop
      
      f1syed = f1syed*0.25d0
      f2syed = f2syed*0.25d0
      return
      end
**************************************************      
      subroutine twoloop_light_heavy_Syed_quad(t,mo,mho,ires,f1syed,f2syed)
      use handyg
*      implicit none      
      implicit complex*32 (Z)
      double complex f1syed,f2syed
      complex*32 f11,f12,f13,f14
      complex*32 f21, f22
      real*16 z
      double precision t,del,mo,mho
      real*16 m,mh
      complex*32 cone,x,s,xs,xvar,carg,xh,umxs
      parameter(cone=(1.q0,0.q0))
      complex*32 i_
      parameter(i_=(0.q0,1.q0))
      real*16 log2,sq2
      real*16 zeta2,zeta3
      double precision alpha,convfac,pis,ame,amu      
      common/parameters/ame,amu,convfac,alpha,pis
*
      real*16 b5a1,b5a2,b5a3,b5a4,b6a3,b6a4
      complex*32 b5a1b,b5a2b,b5a3b,b5a4b,b6a3b,b6a4b
      integer ifirst,ico,i,ires
      real*16 PL4unm,logmomh
      real*16 pi2,pi4,pi
      parameter (pi=4q0*atan(1q0))
      complex*32 Li22,Li32,logupz,logumz,li2momh,li2momh2
      complex*32 logA,logB,logC,atanhz
      common/twolcmnlightheavyquad/log2,sq2,pi2,pi4,zeta2,zeta3,
     .     b5a1,b5a2,b5a3,b5a4,b6a3,b6a4,Li22,Li32,logmomh,
     .     li2momh,li2momh2,logupz,logumz,atanhz,PL4unm,ifirst,ico
      complex*32 GPL(1:882)
      common/GPLscommonquad/GPL,f11,f12,f13,f14,f21,f22,f1comb,f2comb
      data ifirst /0/

*********** for HPLOG
      integer n1,n2,nw,i1,i2,k1,k2,k3
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
      parameter(n2=1)
      parameter(nw=3)
      double complex Hc1,Hc2,Hc3,Hc4,li3A,li3B,li2A,li2B
      double precision  Hr1,Hr2,Hr3,Hr4,Hupx,Humx
      double precision  Hi1,Hi2,Hi3,Hi4
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2),
     .     Hc4(n1:n2,n1:n2,n1:n2,n1:n1)
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),
     .     Hr4(n1:n2,n1:n2,n1:n2,n1:n1)
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),
     .     Hi4(n1:n2,n1:n2,n1:n2,n1:n1)
***** 
      call clearcache
*************************      
      s  = t!/mo/mo
      m  = mo!1.q0!mo
      mh = mho!/mo
      z = m/mh

      xs = sqrt(1.q0-4.q0*m*m/s)
      umxs = 1q0 - xs
      
      if (ifirst.eq.0.or.ires.gt.0) then

c         call GPLopts(1q-20, 10000q0,  2.q0)
         
         ico = 0
         log2= log(2.q0)
         logmomh = log(m/mh)
         
         sq2 = sqrt(2.q0)
         zeta2=pi*pi/6.q0
         zeta3=1.202056903159594285399738161511449990764986292q0         
         pi2= pi*pi
         pi4= pi2*pi2

         b5a1 =  2.77258872223978109212q0
         b5a2 =  1.36805607802364681993647982810813818502q0
         b5a3 =  1.13556770870411119047364593456472384651q0
         b5a4 =  1.0583097143302935q0
         b6a3 =  8.41439832211715975835164054594208875911q0
         b6a4 = -8.30035652771381970800065987017540898127q0
* more precise
         b5a1 = 2.77258872223978123766892848583270669q0
         b5a2 = 1.36805607802364717427642022798538871q0
         b5a3 = 1.13556770870411139095090943038957078q0
         b5a4 = 1.05830971433029308016425544951004003q0
         b6a3 = 8.41439832211715999779816713065087404q0
         b6a4 =-8.30035652771381950065162465679101891q0

         
         PL4unm = 0.517479061673899386330758161898862945622377475141q0

c         Li22 = myli2(cone*2.q0)
c         call hplog(2.q0,nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c         Li22 = Hc2(0,1)         
c         call hplog(2.q0,nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c         Li32 = Hc3(0,0,1)
         
         Li22=-G([inum(0.q0,0),inum(1.q0,0)],inum(2.q0,di0))
         Li32=-G([inum(0.q0,0),inum(0.q0,0),inum(1.,0)],inum(2.q0,di0))  
         
c     call hplog(z,nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c     li2momh  = Hc2(0,1) ! not ok for mh<m!
         li2momh=-G([inum(0.q0,0),inum(1.q0,0)],inum(z,di0))
c     call hplog(z*z,nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c     li2momh2 = Hc2(0,1)  ! not ok for mh<m!      
         li2momh2=-G([inum(0.q0,0),inum(1.q0,0)],inum(z*z,di0))
         
c         logupz  = log(cone*(1.q0+z))
c         logumz  = log(cone*(1.q0-z))
         logupz=G([inum(cone*0.q0,0)],inum(cone*(1.q0+z),di0))
         logumz=G([inum(cone*0.q0,0)],inum(cone*(1.q0-z),di0))

         atanhz = 0.5q0*(logupz-logumz) ! sostituito nel form direttamente
         
         include 'constgpl_defs.f'

c         print*,b5a1
c         print*,b5a2
c         print*,b5a3
c         print*,b5a4
c         print*,b6a3
c         print*,b6a4
c         include 'b56_withgpls.f'
c         print*,b5a1b
c         print*,b5a2b
c         print*,b5a3b
c         print*,b5a4b
c         print*,b6a3b
c         print*,b6a4b
c         print*,b5a1b-b5a1
c         print*,b5a2b-b5a2
c         print*,b5a3b-b5a3
c         print*,b5a4b-b5a4
c         print*,b6a3b-b6a3
c         print*,b6a4b-b6a4
         
         ifirst = 1
      endif

      ico = ico + 1
      
      x = sqrt(cone*(4.q0*m*m-s))*sqrt(cone*(4.q0*m*m-4.q0*mh*mh-s))
      x = -2.q0*m*m -x+s + sq2*sqrt(cone*(-s*(4.q0*m*m-2.q0*mh*mh+x-s)))
      x = -x/m/mh*0.5q0
      
      xvar =
     &  + m**(-1) * (  - 1.Q0/2.Q0*s*mh**(-1) + 1.Q0/2.Q0*sqrt( - s)*
     &    sqrt( - 4*mh**2 - s*xs**2)*mh**(-1)*xs - 1.Q0/2.Q0*sqrt(2*s*
     &    mh**2 + s**2*xs**2 - sqrt( - s)*sqrt( - 4*mh**2 - s*xs**2)*s*
     &    xs)*mh**(-1)*sq2 )
      xvar = xvar + z

c      x = xvar
            
      include 'vargpl_defs.f'

c      call hplog(dreal( -1q0 + 2q0/(1q0 - xs**2) - 2q0/(1q0 - xs**2)*xs),
c     +     nw,Hc1,Hc2,Hc3,Hc4,
c     +        Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)         
c      li2A = Hc2(0,1)
c      li3A = Hc3(0,0,1)

      carg = -1q0 + 2q0/(1q0 - xs**2) - 2q0/(1q0 - xs**2)*xs
      li2A = -G([inum(0.q0,0),inum(1.q0,0)],inum(carg,di0))
      li3A = -G([inum(0.q0,0),inum(0.q0,0),inum(1.q0,0)],inum(carg,di0))
      
!      carg = 2q0 - 2q0/(1q0 - xs**2) + 2q0/(1q0 - xs**2)*xs, i.e.
      carg = 1.q0 - carg
      li2B = -G([inum(0.q0,0),inum(1.q0,0)],inum(carg,di0))
      li3B = -G([inum(0.q0,0),inum(0.q0,0),inum(1.q0,0)],inum(carg,di0))

      logA = Log(cone*( - 1 + 2/(1 - xs**2) - 2/(1 - xs**2)*xs))
      logB = Log(cone*(2 - 2/(1 - xs**2) + 2/(1 - xs**2)*xs))
c not neede logA - logC = i_*pi      logC = Log(cone*(1 - 2/(1 - xs**2) + 2/(1 - xs**2)*xs))
c      ie logA = G([inum(0.q0,0)],inum(cone*(( - 1 + 2/(1 - xs**2) - 2/(1-xs**2)*xs)),di0))
c      ie logB = G([inum(0.q0,0)],inum(cone*(( + 1 - 2/(1 - xs**2) + 2/(1-xs**2)*xs)),di0))
ccc      print*,'!!!!REINCLUDE!!!!'
c      include 'f1comb_light_heavy.f'
c     include 'f1_onlyf11_expanded_in_z.f'

c      print*,'REINCLUDE IN QUAD! stopping'
c      stop
      
      include 'f1_light_heavy_quad.f'      
      f1syed = f11 - f12 - f13 + f14
      include 'f2_light_heavy_quad.f'
      f2syed = f21 - f22
      
      f1syed = f1syed*0.25q0
      f2syed = f2syed*0.25q0
      return
      end
**************************************************
      subroutine getGPLbfp(G0,Gmu0,Gmumu0,x,xb)
      use handyG
      implicit none      
      double complex G0,Gmu0,Gmumu0,i,c1,c0
      double precision zeta2,pi,zeta3
      double complex xb
      double precision x
      parameter (pi=4.d0*atan(1.d0))
      parameter (zeta2= pi*pi/6.d0)
      parameter (i=(0.d0,1.d0))
      parameter (c1=(1.d0,0.d0))
      parameter (c0=(0.d0,0.d0))
      parameter (zeta3=
     .     1.202056903159594285399738161511449990764986292d0)
* for hplog ------------
      integer n1,n2,nw,i1,i2
      parameter(n1=-1)
      parameter(n2=1)
      parameter(nw=3)
      double complex Hc1,Hc2,Hc3,Hc4
      double precision     Hr1,Hr2,Hr3,Hr4
      double precision     Hi1,Hi2,Hi3,Hi4
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $     Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
*  FOR CHAPLIN
      double complex HPL1,HPL2,HPL3,HPL4
      external HPL1,HPL2,HPL3,HPL4
      double complex z
      double precision epsz
      double precision rez,imz
      parameter(epsz=1.d-16)
      integer j1,j2,j3,j4
************
      z   =  xb
c     call hplog(dble(xb),nw,Hc1,Hc2,Hc3,Hc4,
c     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

******THIS IS FOR CHAPLIN **********************************
c      do j1=-1,1
c         Hc1(j1)= HPL1(j1,z)
c         do j2= -1,1
c            Hc2(j1,j2)= HPL2(j1,j2,z)
c            do j3= -1,1
c               Hc3(j1,j2,j3)= HPL3(j1,j2,j3,z)
c            enddo
c         enddo
c      enddo
     
c      print*,HC1
      
c      G0     =  2.d0*HC1(-1)-HC1(0)
c      Gmu0   = -i*zeta2+2.d0*i*HC2(0,-1)-i*HC2(0,0)
c     Gmumu0 = -2.d0*zeta3+zeta2*HC1(0)-2.d0*HC3(0,0,-1)+HC3(0,0,0)

      G0     =  2.d0*HPL1(-1,z)-HPL1(0,z)
      Gmu0   = -i*zeta2+2.d0*i*HPL2(0,-1,z)-i*HPL2(0,0,z)
      Gmumu0 = -2.d0*zeta3+zeta2*HPL1(0,z)
     .     -2.d0*HPL3(0,0,-1,z)+HPL3(0,0,0,z)
      
c      call clearcache      
c      G0     =  2.d0*G([-c1,z])-G([c0,z])
c      Gmu0   = -i*zeta2+2.d0*i*G([c0,-c1,z])-i*G([c0,c0,z])
c      Gmumu0 = -2.d0*zeta3+zeta2*G([c0,z])
c     .     -2.d0*G([c0,c0,-c1,z])+G([c0,c0,c0,z])
      
c      print*,G0,Gmu0,Gmumu0
c      print*,G([0.d0,x])
c      print*,G([-1.d0/3.d0,0.d0,x])
      
c      print
      
      return
      end
****************************************************
      subroutine boncianiferrogliapeninLH(t,mh2,f1)
* https://arxiv.org/pdf/0802.2215.pdf eq. 179      
      implicit none
      double precision t,mh,mh2
      double precision x,arg
      double complex xb,sq
      double complex f1,c1,i
      parameter (c1 = (1.d0,0.d0))
      parameter (i  = (0.d0,1.d0))
      double complex Gmu0,G0,Gmumu0
      if (t.lt.0.d0) then
         x = -t/mh2      
!     yb eq. 80
         arg = -t-4.d0*mh2 
         if (arg.ge.0.d0) then
            sq = sqrt(arg)
         else
            sq = -i*sqrt(-arg)
         endif
         xb = sqrt(-t) - sq
         xb = xb/(sqrt(-t) + sq)
      else
         x = -t/mh2
!     xb eq. 79
         arg = -t-4.d0*mh2 
         if (arg.ge.0.d0) then
            sq = sqrt(arg)
            xb = sqrt(t) - sq
            xb = xb/(sqrt(t) + sq)
         else
            sq = -i*sqrt(-arg)
            xb = sqrt(-t*c1) - sq
            xb = xb/(sqrt(-t*c1) + sq)
! ovvero            
c            xb = -i*sqrt(t)+i*sqrt(t+4.d0*mh2)
c            xb = xb/(-i*sqrt(t)-i*sqrt(t+4.d0*mh2))
         endif
      endif
      
      call getGPLbfp(G0,Gmu0,Gmumu0,x,xb)

      arg = x * (4.d0-x)
      if (arg.ge.0.d0) then
         sq = sqrt(arg)
      else
         if (x.lt.0.d0) then
            sq = -i*sqrt(-arg)
         else
            sq = i*sqrt(-arg)
         endif
      endif
      
      if (abs(sq).eq.0.d0) sq = 1.d0 ! to catch the case t == -4.d0*amu*amu, when sq = 0 and Gmu0 = 0
      
      f1 = (x-4.d0)*(46.d0-19.d0*x)*Gmu0/36.d0/x/sq
      f1 = f1 
     .     -1.d0/1296.d0/x/x*(
     .     8568*x-3355d0*x*x - x*(3960d0-1590d0*x)*G0
     .     + (1296.d0-216.d0*x*x)*Gmumu0
     .     )

      f1 = 4.d0*f1
      return
      end
***

#endif
