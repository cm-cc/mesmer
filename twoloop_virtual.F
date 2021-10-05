      subroutine integrate547(me,del,aint,aints,aintex)
      implicit double precision (a-h,k,m,o-z)
      double precision aints(10)
      parameter (pi    = 4.d0*atan(1.d0))
      parameter (zeta2 = pi*pi/6.d0)

      k2min = 4.d0*me*me
      k2max = del*del

      np  = 1000000000
      dk2 = (k2max - k2min)/np

      k2      = k2min - dk2
      aint    = 0.d0
      aints   = 0.d0
      aintex  = 0.d0
      do i = 1,np
         k2 = k2min + (i-1)*dk2

         dbar = del + sqrt(abs(del*del-k2))
         dbar = dbar/(del - sqrt(del*del-k2))         
         
         ovall = sqrt(1.d0-4.d0*me*me/k2)/k2

         fex = (1.d0+2.d0*me*me/k2)*(
     .        0.25d0*log(dbar)**2-dilogmy(-1.d0/dbar)
     .        +log(dbar)*log(1.d0+1.d0/dbar)-zeta2+pi*pi/12.d0
     .        +sqrt(1.d0-k2/del/del)-0.5d0*log(dbar)
     .        )
         fex = fex * ovall
         
         f1 = 0.25d0*log(4.d0*del*del/k2)**2
         f2 = - 0.25d0*k2/del/del
         f3 = 0.25d0*k2/del/del*log(4.d0*del*del/k2)
         f4 = - zeta2 + pi*pi/12.d0 + 1.d0
         f5 = - 0.5d0 * log(4.d0*del*del/k2)

         fa = (f1+f2+f3+f4+f5) * ovall
         fb = (f1+f2+f3+f4+f5) * ovall *2.d0*me*me/k2                  
*     fa: f1  = c9*F9, f2 = c8*F8, f3 = wrong, f4 = c7*F7, f5 = c6*F6
*     fb: f1  = c1*F1, f2 = c3*F3, f3 = c2*F2, f4 = c4*F4, f5 = c5*F5
         f = fa + fb
         aint = aint + f
         aints(9)  = aints(9)  + f1 * ovall
         aints(8)  = aints(8)  + f2 * ovall
         aints(10) = aints(10) + f3 * ovall
         aints(7)  = aints(7)  + f4 * ovall
         aints(6)  = aints(6)  + f5 * ovall
         aints(1)  = aints(1)  + f1 * ovall *2.d0*me*me/k2
         aints(3)  = aints(3)  + f2 * ovall *2.d0*me*me/k2
         aints(2)  = aints(2)  + f3 * ovall *2.d0*me*me/k2
         aints(4)  = aints(4)  + f4 * ovall *2.d0*me*me/k2
         aints(5)  = aints(5)  + f5 * ovall *2.d0*me*me/k2
         aintex    = aintex + fex
      enddo
      aint   = aint   * dk2
      aints  = aints  * dk2
      aintex = aintex * dk2
      return
      end
*************************************************************      
      double precision function Fsoftrealpairs(del, me)
! written by Ettore Budassi      
      implicit none
      double precision del, me,Fsr,zeta2,zeta3,aint
      double precision x, upx,umx,dilogmy,pi,dnli3,aintex
      double precision c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
      double precision f1,f2,f3,f4,f5,f6,f7,f8,f9,f10
      double precision aints(10)
      external dilogmy
      external dnli3
      integer k
      parameter (pi=4d0*atan(1d0))
      double precision ame,ammu,convfac,alpha,pip
      common/parameters/ame,ammu,convfac,alpha,pip

*********** for HPLOG
      integer n1,n2,nw,i1,i2,k1,k2,k3
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
      parameter(n2=1)
      parameter(nw=3)
      double complex Hc1,Hc2,Hc3,Hc4,Hupx,Humx
      double precision     Hr1,Hr2,Hr3,Hr4
      double precision     Hi1,Hi2,Hi3,Hi4
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $     Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
*****

      x=sqrt(1.d0 - 4.d0*me*me/del/del)
      upx=1.d0+x
      umx=1.d0-x

c      call hplog(upx*0.5d0,nw,Hc1,Hc2,Hc3,Hc4,
c     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
c      Hupx = Hc2(0,1)
c      call hplog(umx*0.5d0,nw,Hc1,Hc2,Hc3,Hc4,
c     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
c      Humx = Hc2(0,1)
c      print*,del,dilogmy(upx*0.5d0)/Hupx,dilogmy(umx*0.5d0)/Humx

*NEWNEW************************************************************************
      c1=1d0/8d0
      f1=log(4d0)*log(4d0)*2d0/3d0*x*x*x-4d0/3d0*(log(me*me/del/del)
     + *(log(upx/umx)-2d0/3d0*x*x*x-2d0*x)+2d0/3d0*x*x*x*log(1d0-x*x)
     + -4d0/9d0*x*x*x-4d0/3d0*x+2d0/3d0*log(upx/umx)-log(2d0)*log(upx)
     + +dilogmy(upx/2d0)-pi*pi/12d0+log(2d0)*log(2d0)/2d0+log(umx)
     + *log(umx)/2d0-2d0*umx*log(umx)-2d0*x-log(upx)*log(upx)/2d0
     + +log(2d0)*log(umx)-dilogmy(umx/2d0)+pi*pi/12d0-log(2d0)*log(2d0)
     + /2d0 + 2d0*(upx)*log(upx)-2d0*x)

      Fsr =c1*f1
      
      c2=me*me/del/del      
      f2=log(del*del/me/me)*(log(upx/umx)/2d0-x)
     + +log(upx)*log(2d0)/2d0
     + -dilogmy(upx/2d0)/2d0+pi*pi/24d0-log(2d0)*log(2d0)/4d0
     + -log(umx)*log(umx)/4d0+umx*log(umx)+x+log(upx)*log(upx)/4d0
     + -log(2d0)*log(umx)/2d0 +dilogmy(umx/2d0)/2d0 -pi*pi/24d0
     + +log(2d0)*log(2d0)/4d0-upx*log(upx)+x
      
      Fsr = Fsr+c2*f2
     
      c3=-me*me/del/del
      f3=log(upx/umx)/2d0-x

      Fsr = Fsr+c3*f3

      c4=1d0-pi*pi/6d0+pi*pi/12d0
      f4=x*x*x/3d0

      Fsr = Fsr+c4*f4

      c5= 1d0
      f5= x/3d0+x*x*x/9d0-log(upx/umx)/6d0-x*x*x*log(4d0)/6d0

      Fsr = Fsr +c5*f5

      c6=-1d0
      f6 = f2

      Fsr = Fsr +c6*f6

      c7=1d0-pi*pi/6d0+pi*pi/12d0
      f7=-2d0*x+log((upx)/(umx))

      Fsr = Fsr +c7*f7

      c8=-me*me/del/del
      f8=del*del/4d0/me/me*x-log((upx)/(umx))/2d0

      Fsr = Fsr +c8*f8

      c9=1d0
      f9=log(me*me/del/del)*log(me*me/del/del)/2d0*(log(upx/umx)/2d0
     + -x)-log(me*me/del/del)*(log(2d0)*log(upx)/2d0-dilogmy(upx/2d0)
     + /2d0+pi*pi/24d0-log(2d0)*log(2d0)/4d0-log(umx)*log(umx)/4d0
     + +umx*log(umx)+x+log(upx)*log(upx)/4d0-log(2d0)*log(umx)/2d0
     + +dilogmy(umx/2d0)/2d0-pi*pi/24d0+log(2d0)*log(2d0)/4d0+x
     + -upx*log(upx))-log(umx)*log(umx)*log(umx)/12d0+log(upx)*log(upx)
     + *log(upx)/12d0+(log(umx)*log(umx)*log(upx/2d0)-2d0*dnli3(umx/2d0)
     + +2d0*log(umx)*dilogmy(umx/2d0))/4d0+(2d0*dnli3(upx/2d0)-log(upx)
     + *log(upx)*log(umx/2d0)-2d0*log(upx)*dilogmy(upx/2d0))/4d0
     + +(log(umx)*dilogmy(umx/2d0)-dnli3(umx/2d0)-log(2d0)*log(umx)
     + *log(umx)/2d0)/2d0+(dnli3(upx/2d0)-dilogmy(upx/2d0)*log(upx)
     + +log(2d0)*log(upx)*log(upx)/2d0)/2d0-log(umx)*x*log(upx)+x
     + *log(upx)+log(upx)*log(umx)+log(upx)-2d0*log(umx)*log(upx/2d0)
     + +log(umx)*x-2d0*x-2d0*dilogmy(umx/2d0)-log(umx)-log(2d0)*log(2d0)
     + +pi*pi/6d0-(upx*log(upx)*log(upx)-2d0*upx*log(upx)+2d0*x)/2d0
     + -(2d0*x+2d0*umx*log(umx)-umx*log(umx)*log(umx))/2d0

      Fsr = Fsr +c9*f9
**** nuovo EB ottimo
c      print*,'10-',c10*f10
c      c10=me*me/del/del/2d0
c      f10=-2d0*log(me*me/del/del)*(del*del/4d0/me/me*x-log(upx/umx)
c     + /2d0)-log(upx)*log(upx)/2d0-log(upx)*log(2d0)+dilogmy(upx/2d0)
c     + -pi*pi/12d0+log(2d0)*log(2d0)/2d0+log(umx)*log(2d0)
c     + -dilogmy(umx/2d0)+pi*pi/12d0-log(2d0)*log(2d0)/2d0+log(umx)
c     + *log(umx)/2d0+log(upx)/umx-log(upx/umx)/2d0-1+(log(umx)+1)/umx
c     + +1-(log(upx)+1)/upx-log(umx)/upx-log(upx/umx)/2
c      print*,'10-',c10*f10
     
***** mio !! OK
      c10 = me*me/del/del
      f10 = 1d0/24d0*
     .     (
     .     12*x/upx
     .     + pi*pi
     .     - 6*log(upx/umx)-6*log(2d0)**2
     .     + 12*log(2d0)*log(umx)
     .     + 24d0*x*log(upx)/upx/umx
     .     -6*log(upx)**2
     .     -12d0*dilogmy(umx/2.d0)
     .     ) 
     .     - 1d0/24d0/upx/umx*
     .     (
     .     pi*pi
     .     - x * (12d0 + x*(12d0+pi*pi) )
     .     - 6d0*upx*umx*log(2d0)**2
     .     + 6.d0*( -1d0+(x-4d0)*x )*log(umx)
     .     - 6d0*upx*umx*log(umx)**2
     .     + 6d0*upx*umx*( 1d0+2.d0*log(2d0) )*log(upx)
     .     -12d0*upx*umx*dilogmy(upx/2.d0) 
     .     )
      f10 = f10
     .     + log(del*del/me/me)*(x/umx/upx-1.d0/2.d0*log(upx/umx))

      Fsr= Fsr + c10*f10
**********
cc      print*,'EB ',Fsr
c      call integrate547(me,del,aint,aints,aintex)
c      print*,Fsr,aint,aintex,'ratio app/ex=',aint/aintex,
c     .     'ex - approx',aintex-aint
c      do k = 1,10
c         print*,k,aints(k)
c      enddo
******************************************************************************
      Fsoftrealpairs = Fsr/3.d0*(alpha/pi)**2 * 0.25d0 ! need to be divided by 4, to match the a^2/(2*pi)^2 normalization!

      return
      end
***************************************************************************      
*     2L vv and vrsoft with form factors (according to Mastrolia-Remiddi
*     NPB664 (2003) and hep-ph/0406046 and Burgers 1985) 
      FUNCTION TWOLOOP_TWOBODY(S,T,U)
      implicit none
      double precision twoloop_twobody,z,xx,carlolo
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin
      integer ihard,k
      double precision tmin,tmax,dt
      
      double complex F1NF1, F2NF1, F1NF1ex, F2NF1ex
      double complex f1lh,f2lh,f1lhq,f2lhq

      common/carlolocmn/carlolo
      
      common/parameters/ame,amu,convfac,alpha,pi      

cc      COMMON/CONST/EMUBEAM,EEBEAM,ALPHA,CONVFAC,PI,AME,AMU
      COMMON/CUTOFF/EGMIN,IHARD
*
      double precision s,t,u
*
      double precision mm2,me2,aa,bb,cc,dd,burg
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision deltauv,mu2dim,lambda2,m2min
      double precision rsoft,rsoft2
      double precision rflux
      double precision kallen
      external kallen
      double precision x,omx,opx,omx2,omx3,omx4,omx5,opx2,opx3,opx4,opx5
      double precision    zeta2,zeta3
      double complex f1cur_1L_epsm1,f1cur_1L_eps0,f1cur_1L_eps1
      double complex f2cur_1L_eps0


      double complex f1e_1L_epsm1,f1e_1L_eps0,f1e_1L_eps1
      double complex f2e_1L_eps0
      double complex f1mu_1L_epsm1,f1mu_1L_eps0,f1mu_1L_eps1
      double complex f2mu_1L_eps0
      double complex F1e,F1mu,F2e,F2mu
      double complex F1ed,F1mud,F2ed,F2mud
      
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0
      double complex f1cur_2L_log2_expanded_sh,f1cur_2L_log1_expanded_sh
     +           ,f1cur_2L_log0_expanded_sh
      double complex f2cur_2L_log1_expanded_sh,f2cur_2L_log0_expanded_sh
      double complex f1cur_2L_log2_b,f1cur_2L_log1_b,f1cur_2L_log0_b
      double complex f1cur_2L_log0_pairs_b
      double complex f2cur_2L_log1_b,f2cur_2L_log0_b
      double complex F1,F1d,F2,F2d
      double complex F1oal,F1doal,F2oal,F2doal
      double complex F1oal2,F1doal2,F2oal2,F2doal2
      double complex f1f1d_2L,f1f2d_2L
      double complex matrixel,matrixelvv,matrixelrv
      double complex matrixelvve,matrixelvvmu
      double complex arg
      double complex soft
      double complex cone
      parameter(cone=(1.d0,0.d0))
      double precision log2
      double precision f1timesf1d,f2timesf2d
      double precision virt
      integer itmu,ite
      common/subset/itmu,ite
      common/gammamass/lambda2

      integer ifirst,ic
      double precision ao2pi,ao2pi2,pi2,pi4
      common/twolcmn1/zeta2,zeta3,mm2,me2,ao2pi,ao2pi2,pi2,pi4,ic,ifirst
      data ifirst,ic /0,0/

      double complex virtvupvdwn,F12L,F22L,F12Ld,F22Ld
      double precision carlo2L
      common/testvirtvevmu/virtvupvdwn
      
      double precision phmass
      common/photonmasslambda/phmass
      double precision tmincut,tmaxcut
      common/tmintmax/tmincut,tmaxcut
**** charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu      
*     
#include "looptools.h"
c#ifdef QUAD
c      real*16 uquad
c      parameter (uquad=1.q0)
c#else
c      real*16 uquad
c      parameter (uquad=1.d0)
c#endif
*
c      deltauv= getdelta()
c      mu2dim= getmudim()
c      lambda2= getlambda()
c      m2min= getminmass()
*
      if (ifirst.eq.0) then
         lambda2= phmass*phmass

         pi2 = pi*pi
         pi4 = pi2*pi2
         
         zeta2=pi2/6.d0
         zeta3=1.202056903159594285399738161511449990764986292d0
*     
         mm2= amu*amu
         me2= ame*ame

         ao2pi  = alpha*0.5d0/pi
         ao2pi2 = ao2pi*ao2pi
         
         itmu = 0
         ite  = 1
         
         ifirst = 1
      endif
      ic = ic + 1
*
c      p1p2 = (  s - mm2 - me2 ) * 0.5d0 
c      p1p3 = ( -t + 2.d0*mm2  ) * 0.5d0 
c      p1p4 = ( -u + mm2 + me2 ) * 0.5d0 
c      p2p3 = ( -u + mm2 + me2 ) * 0.5d0 
c      p2p4 = ( -t + 2.d0*me2  ) * 0.5d0 
c      p3p4 = (  s - mm2 - me2 ) * 0.5d0 
*
*     F1, F2       {F1= 1 + (alpha/2/pi) F11 + (alpha/2/pi)^2 F12}
*                  {F2=     (alpha/2/pi) F21 + (alpha/2/pi)^2 F22}
      matrixelvve  = (0.d0,0.d0)
      matrixelvvmu = (0.d0,0.d0)

** E LINE ********************************************************************************
      IF (QRe.eq.1) THEN         
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
         x = (sqrt(-t+4.d0*me2)-sqrt(-t))/(sqrt(-t+4.d0*me2) + sqrt(-t))
      
c      call twoloopff_NPB664_2003_pp_sh(x,me2,t,
c     +     f1cur_2L_log2,
c     +     f1cur_2L_log1,
c     +     f1cur_2L_log0,
c     +     f2cur_2L_log1,
c     +     f2cur_2L_log0)
*************************************************************************
      call twoloopff_Pp_Carlo(x,me2,t,
     +     f1cur_2L_log2,
     +     f1cur_2L_log1,
     +     f1cur_2L_log0,
     +     f2cur_2L_log1,
     +     f2cur_2L_log0)
*************************************************************************
      
carlo      call twoloopff_NPB664_2003_expanded_sh(me2,t,
carlo     +     f1cur_2L_log2,
carlo     +     f1cur_2L_log1,
carlo     +     f1cur_2L_log0,
carlo     +     f2cur_2L_log1,
carlo     +        f2cur_2L_log0)
c      call twoloopff_burgers(me2,t,
c     +         f1cur_2L_log2_b,f1cur_2L_log1_b,f1cur_2L_log0_b,
c     +         f2cur_2L_log1_b,f2cur_2L_log0_b,f1cur_2L_log0_pairs_b)
!! subtracting pairs !!!      
      call twoloop_NF1(x,me2,t,F1NF1,F2NF1)
      f1cur_2L_log0 =  f1cur_2L_log0 - F1NF1
      f2cur_2L_log0 =  f2cur_2L_log0 - F2NF1

****      
c      call burgersF1timelike(burg,1000d0,me2,mm2)
c      print*,burg
c      call twoloop_light_heavy_Syed(-200.d0*amu*amu,ame,amu,0,f1lh,f2lh)
c      print*,f1lh
c      call twoloop_NF1_TL(10.d0,mm2,t,F1NF1,F2NF1)
c      call boncianiferrogliapeninLH(-200.d0*amu*amu,mm2,f1lh)
c      print*,f1lh
c      print*,'stop twoloop_virtual.F line 403'
c      stop
***      
      
c      call twoloop_light_heavy_Syed(t,ame,amu,0,f1lh,f2lh)
c      print*,F1NF1,F2NF1,f1lh,f2lh
c      f1cur_2L_log0 =  f1cur_2L_log0 + f1lh
c      f2cur_2L_log0 =  f2cur_2L_log0 + f2lh

c      call interpolate_F1F2_LH_HL(t,tmincut,tmaxcut,ame,amu,aa,bb,cc,dd)
c      f1cur_2L_log0 =  f1cur_2L_log0 + aa
c      f2cur_2L_log0 =  f2cur_2L_log0 + bb
      
c      tmin = -0.14d0
c      tmax = -5d-4
c      dt = (tmax - tmin)/10000
c      t = tmin - dt
c      do k = 1,10000
c         t = tmin + k*dt
c        call twoloop_light_heavy_Syed(t,ame,amu,0,f1lh,f2lh)
cc         call interpolate_F1F2_LH_HL(t,tmin,tmax,ame,amu,aa,bb,cc,dd)
c         x = (sqrt(-t+4.d0*me2)-sqrt(-t))/(sqrt(-t+4.d0*me2) + sqrt(-t))
c         call twoloop_NF1(x,me2,t,F1NF1,F2NF1)
cc         print*,t,dble(F1NF1),dble(F2NF1)
cc         call twoloop_light_heavy_Syed_quad(t,ame,amu,f1lhq,f2lhq)
c     print*,t,dble(f1lh),dble(f2lh),dble(F1NF1),dble(F2NF1)
c        call boncianiferrogliapeninLH(t,mm2,f2lh)
c        print*,t,dble(f1lh),dble(f2lh)
c      enddo
c      stop
      
cc ONLY PAIRS !      
c      f1cur_2L_log2 = 0.d0
c      f1cur_2L_log1 = 0.d0
c      f1cur_2L_log0 = f1lh ! F1NF1
c      f2cur_2L_log1 = 0.d0
c      f2cur_2L_log0 = f2lh ! F2NF1
* contributions to two loops from pure two loops (with lambda regularization)
*
      f1f1d_2L= 2.d0*f1cur_2L_log2*(log(lambda2/me2))**2
     +     +2.d0*f1cur_2L_log1*log(lambda2/me2)
     +     +2.d0*f1cur_2L_log0
      f1f1d_2L= ao2pi2*f1f1d_2L
      F12L = f1f1d_2L

      F22L = f2cur_2L_log1 * log(lambda2/me2)
     +     + f2cur_2L_log0
      F22L = 2.d0*ao2pi2 * F22L
*
*     this is the two-loop virtual contribution
      carlo2L = 0.d0
     &  + F22L * ( 8*t**2 + 16*mm2*t )
      carlo2L = carlo2L + F12L * ( 8*t**2 + 16*s*t + 16*s**2 - 32*me2*s
     &     + 16*me2**2 - 32*mm2*s + 32*mm2*me2 + 16*mm2**2 )
      carlo2L = carlo2L * 64.d0*pi4*ao2pi2/t/t !overall factor for born

!! this is ok      carloLO = ( 8*t**2 + 16*s*t + 16*s**2 - 32*me2*s
c     &     + 16*me2**2 - 32*mm2*s + 32*mm2*me2 + 16*mm2**2 )
c      carloLO = carloLO * 64.d0*pi4*ao2pi2/t/t !overall factor for born

      matrixelvve = carlo2L
      ENDIF
********************************************************************************
** MU LINE *********************************************************************
      IF (QRmu.eq.1) THEN
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
      x= (sqrt(-t+4.d0*mm2) - sqrt(-t))/(sqrt(-t+4.d0*mm2) + sqrt(-t))
c      call twoloopff_NPB664_2003_pp_sh(x,mm2,t,
c     +     f1cur_2L_log2,
c     +     f1cur_2L_log1,
c     +     f1cur_2L_log0,
c     +     f2cur_2L_log1,
c     +     f2cur_2L_log0)
      call twoloopff_Pp_Carlo(x,mm2,t,
     +     f1cur_2L_log2,
     +     f1cur_2L_log1,
     +     f1cur_2L_log0,
     +     f2cur_2L_log1,
     +     f2cur_2L_log0)
      
carlo      call twoloopff_NPB664_2003_expanded_sh(mm2,t,
carlo     +     f1cur_2L_log2,
carlo     +     f1cur_2L_log1,
carlo     +     f1cur_2L_log0,
carlo     +     f2cur_2L_log1,
carlo     +     f2cur_2L_log0)      

!! subtracting pairs !!!
      call twoloop_NF1(x,mm2,t,F1NF1,F2NF1)
c      call twoloop_NF1_expanded(x,mm2,t,F1NF1ex,F2NF1ex)      
      f1cur_2L_log0 = f1cur_2L_log0 - F1NF1
      f2cur_2L_log0 = f2cur_2L_log0 - F2NF1

c      call twoloop_light_heavy_Syed_quad(t,amu,ame,0,f1lh,f2lh)
c      print*,t/amu/amu,F1NF1,F2NF1,f1lh,f2lh
c      stop
      
c      call twoloop_light_heavy_Syed(t,amu,ame,f1lh,f2lh)      
c      f1cur_2L_log0 =  f1cur_2L_log0 + f1lh
c      f2cur_2L_log0 =  f2cur_2L_log0 + f2lh

c      call interpolate_F1F2_LH_HL(t,tmincut,tmaxcut,ame,amu,aa,bb,cc,dd)
c      f1cur_2L_log0 =  f1cur_2L_log0 + cc
c      f2cur_2L_log0 =  f2cur_2L_log0 + dd
      
*
* contributions to two loops from pure two loops (with lambda regularization)
*
      f1f1d_2L= 2.d0*f1cur_2L_log2*(log(lambda2/mm2))**2
     +     +2.d0*f1cur_2L_log1*log(lambda2/mm2)
     +     +2.d0*f1cur_2L_log0
      f1f1d_2L= ao2pi2*f1f1d_2L

      F12L = f1f1d_2L
      F22L = f2cur_2L_log1 * log(lambda2/mm2)
     +     + f2cur_2L_log0
      F22L = 2.d0*ao2pi2 * F22L
*
*     this is the two-loop virtual contribution
      carlo2L = 0.d0
     &  + F22L * ( 8*t**2 + 16*me2*t )
      carlo2L = carlo2L + F12L * ( 8*t**2 + 16*s*t + 16*s**2 - 32*mm2*s
     &     + 16*mm2**2 - 32*me2*s + 32*me2*mm2 + 16*me2**2 )
      carlo2L = carlo2L * 64.d0*pi4*ao2pi2/t/t !overall factor for born

      matrixelvvmu = carlo2L      
      ENDIF
***********************************************************************************
* this line below is for two loop
      twoloop_twobody = dble(matrixelvve) + dble(matrixelvvmu)

      if (QRmu.eq.1.and.QRe.eq.1) then
         x = (sqrt(-t+4.d0*me2)-sqrt(-t))/(sqrt(-t+4.d0*me2) + sqrt(-t))
         call twoloopff_0406046(x,f1e_1L_epsm1,f1e_1L_eps0,
     +        f2e_1L_eps0)
         x= (sqrt(-t+4.d0*mm2)-sqrt(-t))/(sqrt(-t+4.d0*mm2) + sqrt(-t))
         call twoloopff_0406046(x,f1mu_1L_epsm1,f1mu_1L_eps0,
     +        f2mu_1L_eps0)
         
         F1e  = ao2pi * f1e_1L_epsm1 * log(lambda2/me2)
     +        + ao2pi *f1e_1L_eps0
         F2e  = ao2pi *f2e_1L_eps0

         F1mu  = ao2pi * f1mu_1L_epsm1 * log(lambda2/mm2)
     +        + ao2pi *f1mu_1L_eps0
         F2mu  = ao2pi *f2mu_1L_eps0

         virtvupvdwn =
     &        + F2e*F2mu * ( 12*t**2 )
         virtvupvdwn = virtvupvdwn + F2e*F1mu * ( 8*t**2 + 16*mm2*t )
         virtvupvdwn = virtvupvdwn + F1e*F2mu * ( 8*t**2 + 16*me2*t )
         virtvupvdwn = virtvupvdwn + F1e*F1mu * ( 8*t**2 + 16*s*t + 16*
     &        s**2 - 32*me2*s + 16*me2**2 - 32*mm2*s + 32*mm2*me2 + 16*
     &        mm2**2 )

         virtvupvdwn = virtvupvdwn * 64.d0*pi4*ao2pi2/t/t

         twoloop_twobody = twoloop_twobody + 2.d0*dreal(virtvupvdwn)
      endif
      return
      END
*************************************************************************************
* PLB 164 (1985), Burgers two-loop form factors
      subroutine twoloopff_burgers(m2,t,
     +                   f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0,
     +         f2cur_2L_log1,f2cur_2L_log0,f1cur_2L_log0_pairs)
      implicit none
      double precision m2,t
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0,f1cur_2L_log0_pairs
      double complex cone
      double precision q2,u
      double precision log2
      double precision zeta2,zeta3
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin

      common/parameters/ame,amu,convfac,alpha,pi
c      COMMON/CONST/EMUBEAM,EEBEAM,ALPHA,CONVFAC,PI,AME,AMU
      double precision pi2,pi4
*
      integer ifirst
      common/twolcmn2/cone,log2,pi2,pi4,zeta2,zeta3,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         cone= (1.d0,0.d0)
         log2= log(2.d0)
         zeta2=pi*pi/6.d0
         zeta3=1.202056903159594d0
         zeta3=1.202056903159594285399738161511449990764986292d0
         pi2= pi*pi
         pi4= pi2*pi2
         ifirst = 1
      endif
      
      q2 = -t
      u  = -log(q2/m2)
*
* coefficient of log^2(lambda/m) for F1
*
* coefficient of log(lambda/m) for F1
*
* coefficient not IR (no power of log(lambda/m) for F1 
*
      f1cur_2L_log2= (0.d0,0.d0)
      f1cur_2L_log1= (0.d0,0.d0)
      f1cur_2L_log0= (0.d0,0.d0)
      
      f1cur_2L_log0= u**4/32.d0*cone + 31.d0/144.d0*u**3*cone
     +     +(229.d0/288.d0-1.d0/8.d0*pi2/6.d0)*u**2*cone
     +     +(1627.d0/864.d0+13.d0/24.d0*pi2/6.d0-3.d0/2.d0*zeta3)*u*cone
     +     -59.d0/40.d0*pi4/36.d0*cone-9.d0/4.d0*zeta3*cone
     +     -3.d0*(pi2/6.d0)*log2*cone+127.d0/48.d0*pi2/6.d0*cone
     +     +1171.d0/216.d0*cone
      f1cur_2L_log0= f1cur_2L_log0 * 4.d0 ! to match the (alpha/2pi)^2
*                    expansion of hep-ph/0406046 instead of (alpha/pi)^2
*
*     coefficient of log(lambda/m) for F2
*
      f2cur_2L_log1= (0.d0,0.d0)
      f2cur_2L_log0= (0.d0,0.d0)
*
* electron pair virtual contribution
*
      f1cur_2L_log0_pairs= cone*(- 1.d0/36.d0*(log(q2/m2))**3
     +                           +19.d0/72.d0*(log(q2/m2))**2
     +     +(1.d0/3.d0*pi2/6.d0-265.d0/216.d0)*log(q2/m2)
     +                           +383.d0/108.d0-11.d0/6.d0*pi2/6.d0)
      f1cur_2L_log0_pairs= f1cur_2L_log0_pairs * 4.d0 ! to match the (alpha/2pi)^2
*                    expansion of hep-ph/0406046 instead of (alpha/pi)^2
      return
      end
*******************************************************************************
* Equations 4.1 and 4.2 provided by Pierpaolo Mastrolia
      subroutine twoloopff_NPB664_2003_pp_sh(x,m2,t,
     +                   f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0,
     +                                 f2cur_2L_log1,f2cur_2L_log0)
c     implicit none
      implicit double complex (Z)
      double precision x,m2,t
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0
      double precision q2
* for hplog ------------
      double precision omx,opx,omx2,omx3,omx4,omx5,opx2,opx3,opx4,opx5
      integer n1,n2,nw,i1,i2
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
      parameter(n2=1)
      parameter(nw=4)
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
      double complex HPL1,HPL2,HPL3,HPL4,z
      double precision epsz
      double precision rez,imz
      parameter(epsz=1.d-15)
      integer j1,j2,j3,j4
*------------------------------------------------
      double complex cone
      parameter(cone=(1.d0,0.d0))
      double precision log2
      double precision zeta2,zeta3
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin

      common/parameters/ame,amu,convfac,alpha,pi      
c      COMMON/CONST/EMUBEAM,EEBEAM,ALPHA,CONVFAC,PI,AME,AMU
      double precision pi2,pi4
*
      integer ifirst
      common/twolcmn3/log2,pi2,pi4,zeta2,zeta3,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         log2= log(2.d0)
         zeta2=pi*pi/6.d0
         zeta3=1.202056903159594d0
         zeta3=1.202056903159594285399738161511449990764986292d0         
         pi2= pi*pi
         pi4= pi2*pi2
         ifirst = 1
      endif
*     
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
*
****** THIS IS FOR HPLOG ************************************
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,
     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
******THIS IS FOR CHAPLIN **********************************
c      rez =  x
c      imz = -epsz
c      z   =  dcmplx(rez,imz)
c      do j1=-1,1,1
c         Hc1(j1)= HPL1(j1,z)
c         do j2= -1,1,1
c            Hc2(j1,j2)= HPL2(j1,j2,z)
c            do j3= -1,1,1
c               Hc3(j1,j2,j3)= HPL3(j1,j2,j3,z)
c               do j4= -1,1,1
c                  Hc4(j1,j2,j3,j4)= HPL4(j1,j2,j3,j4,z)
c               enddo
c            enddo
c         enddo
c      enddo
*************************************************************
      omx= 1.d0-x
      opx= 1.d0+x
c      omx2= omx*omx
c      omx3= omx2*omx
c      omx4= omx3*omx
c      omx5= omx4*omx
c      opx2= opx*opx
c      opx3= opx2*opx
c      opx4= opx3*opx
c      opx5= opx4*opx
*
      q2= -t
      
      Z1_=opx**(-1)
      Z2_=omx**(-1)
      Z3_=HC1(-1)
      Z4_=HC1(0)
      Z5_=HC1(1)
      Z6_=HC2(0,-1)
      Z7_=HC2(0,1)
      Z8_=HC3(0,0,-1)
      Z9_=HC3(0,0,1)
      Z10_=HC4(0,0,0,-1)
      Z11_=HC4(0,0,0,1)
      Z12_=HC3(0,-1,1)
      Z13_=HC3(0,1,-1)
      Z14_=HC3(0,1,1)
      Z15_=HC4(0,-1,0,1)
      Z16_=HC4(0,0,1,1)
      Z17_=2*Z1_
      Z18_=Z17_ - 5
      Z19_=6*Z1_
      Z20_=Z18_*Z19_
      Z20_=Z20_ + 25
      Z20_=Z20_*Z17_
      Z20_=Z20_ - 15
      Z20_=Z20_*Z1_
      Z20_=Z20_ + 1.D0/2.D0
      Z20_=Z20_*Z1_
      Z21_=1.D0/2.D0*Z2_
      Z20_=Z20_ + Z21_
      Z22_=Z20_*Z5_
      Z23_=1.D0/4.D0*Z2_
      Z24_=21.D0/2.D0 - 13.D0/3.D0*Z2_
      Z24_=Z24_*Z23_
      Z25_=7*Z1_
      Z26_=Z25_ - 23
      Z27_=Z26_*Z17_
      Z27_=293.D0/6.D0 + Z27_
      Z27_=Z1_*Z27_
      Z27_= - 44.D0/3.D0 + Z27_
      Z27_=Z1_*Z27_
      Z27_= - 89.D0/24.D0 + Z27_
      Z27_=Z1_*Z27_
      Z28_=Z2_ - 3.D0/2.D0
      Z28_=Z28_*Z2_
      Z29_= - 11.D0/12.D0 - Z28_
      Z29_=Z2_*Z29_
      Z30_=Z1_ - 5.D0/2.D0
      Z30_=Z30_*Z1_
      Z31_=1.D0/3.D0 + 1.D0/4.D0*Z30_
      Z31_=Z1_*Z31_
      Z31_=1.D0/8.D0 + Z31_
      Z31_=Z1_*Z31_
      Z31_= - 11.D0/192.D0 + Z31_
      Z31_=Z1_*Z31_
      Z29_=1.D0/16.D0*Z29_ + Z31_
      Z29_=Z4_*Z29_
      Z24_=Z29_ + Z22_ + Z24_ + Z27_
      Z24_=Z4_*Z24_
      Z27_=3*Z2_
      Z29_=5.D0/4.D0 - Z2_
      Z29_=Z29_*Z27_
      Z31_=4*Z1_
      Z32_=211 - 83*Z1_
      Z32_=Z32_*Z31_
      Z32_= - 697 + Z32_
      Z32_=Z1_*Z32_
      Z32_=731.D0/4.D0 + Z32_
      Z33_=1.D0/3.D0*Z1_
      Z32_=Z32_*Z33_
      Z34_=Z2_ - 1
      Z35_=Z34_*Z2_
      Z36_=Z1_ - 2
      Z37_=Z36_*Z1_
      Z38_= - 89 - 90*Z37_
      Z38_=Z1_*Z38_
      Z38_= - 1 + Z38_
      Z38_=Z1_*Z38_
      Z38_= - Z35_ + Z38_
      Z38_=Z5_*Z38_
      Z39_=253.D0/2.D0 + 126*Z37_
      Z39_=Z1_*Z39_
      Z39_= - 1.D0/2.D0 + Z39_
      Z39_=Z1_*Z39_
      Z39_=5.D0/2.D0*Z35_ + Z39_
      Z39_=Z3_*Z39_
      Z40_= - 7.D0/4.D0 + Z28_
      Z40_=Z40_*Z21_
      Z41_=Z18_*Z1_
      Z42_=63*Z41_
      Z43_=248 + Z42_
      Z43_=Z1_*Z43_
      Z43_= - 57 + Z43_
      Z43_=Z1_*Z43_
      Z43_= - 7.D0/8.D0 + Z43_
      Z43_=Z1_*Z43_
      Z40_=Z40_ + Z43_
      Z40_=Z6_*Z40_
      Z24_=Z24_ + Z40_ + Z39_ + 2*Z38_ + Z29_ + Z32_
      Z24_=Z4_*Z24_
      Z29_=Z17_ - 3
      Z19_=Z29_*Z19_
      Z19_=Z19_ + 7
      Z19_=Z19_*Z1_
      Z19_=Z19_ - Z2_
      Z32_=4*Z5_
      Z38_=Z19_*Z32_
      Z39_=Z29_*Z1_
      Z40_=89 + 96*Z39_
      Z40_=Z40_*Z1_
      Z43_=7*Z2_
      Z40_=Z40_ + Z43_
      Z44_= - Z3_*Z40_
      Z45_=113 + 30*Z41_
      Z45_=Z45_*Z17_
      Z45_= - 39 + Z45_
      Z45_=Z45_*Z17_
      Z45_= - 7 + Z45_
      Z45_=Z1_*Z45_
      Z45_= - Z43_ + Z45_
      Z46_=4*Z9_
      Z45_=Z45_*Z46_
      Z47_= - 853 - 222*Z41_
      Z47_=Z47_*Z17_
      Z47_=339 + Z47_
      Z47_=Z1_*Z47_
      Z47_=35.D0/2.D0 + Z47_
      Z47_=Z1_*Z47_
      Z47_=35.D0/2.D0*Z2_ + Z47_
      Z47_=Z8_*Z47_
      Z48_=5*Z35_
      Z49_= - 649 - 636*Z37_
      Z49_=Z1_*Z49_
      Z49_=13 + Z49_
      Z49_=Z1_*Z49_
      Z49_= - Z48_ + Z49_
      Z49_=Z6_*Z49_
      Z50_=2.D0/3.D0*Z1_
      Z51_=Z50_ - 1
      Z51_=Z51_*Z1_
      Z52_= - 607.D0/36.D0 - 62*Z51_
      Z52_=Z1_*Z52_
      Z24_=Z24_ + Z49_ + Z47_ + Z45_ + Z44_ + Z38_ - 137.D0/36.D0*Z2_
     &  + Z52_
      Z24_=Z4_*Z24_
      Z38_=Z25_ - 18
      Z44_=Z38_*Z17_
      Z44_=173.D0/6.D0 + Z44_
      Z44_=Z1_*Z44_
      Z44_= - 7.D0/3.D0 + Z44_
      Z44_=Z1_*Z44_
      Z44_= - 145.D0/24.D0 + Z44_
      Z44_=Z1_*Z44_
      Z45_=1.D0/3.D0*Z2_
      Z47_=Z45_ - 1.D0/2.D0
      Z49_= - Z47_*Z43_
      Z49_= - 1.D0/12.D0 + Z49_
      Z23_=Z49_*Z23_
      Z49_=11*Z30_
      Z52_=45.D0/2.D0 + Z49_
      Z52_=Z1_*Z52_
      Z52_= - 25.D0/4.D0 + Z52_
      Z52_=Z1_*Z52_
      Z52_= - 1.D0/48.D0 + Z52_
      Z52_=Z1_*Z52_
      Z23_=Z23_ + Z52_
      Z52_=1.D0/2.D0*Z4_
      Z23_=Z23_*Z52_
      Z53_=71.D0/2.D0 - 17*Z2_
      Z53_=Z2_*Z53_
      Z22_=Z23_ + Z22_ + 1.D0/12.D0*Z53_ + Z44_
      Z22_=Z4_*Z22_
      Z23_= - 32 + 13*Z1_
      Z23_=Z23_*Z50_
      Z23_=Z23_ + 43.D0/3.D0 + 4*log2
      Z23_=Z23_*Z17_
      Z23_=Z23_ - 25.D0/12.D0 - 8*log2
      Z23_=Z1_*Z23_
      Z44_= - 77.D0/2.D0 - 30*Z37_
      Z44_=Z1_*Z44_
      Z44_=17.D0/2.D0 + Z44_
      Z44_=Z1_*Z44_
      Z44_=3.D0/2.D0*Z35_ + Z44_
      Z44_=Z3_*Z44_
      Z50_=7.D0/4.D0 + 3*Z28_
      Z21_=Z50_*Z21_
      Z50_=15*Z41_
      Z53_= - 58 - Z50_
      Z53_=Z1_*Z53_
      Z53_=12 + Z53_
      Z53_=Z1_*Z53_
      Z53_=7.D0/8.D0 + Z53_
      Z53_=Z1_*Z53_
      Z21_=Z21_ + Z53_
      Z21_=Z6_*Z21_
      Z21_=Z22_ + Z21_ + Z44_ - 5.D0/4.D0*Z2_ + Z23_
      Z21_=pi2*Z21_
      Z22_=1 + Z37_
      Z22_=Z22_*Z1_**2
      Z23_= - 485 - 126*Z41_
      Z23_=Z1_*Z23_
      Z23_=195.D0/2.D0 + Z23_
      Z23_=Z1_*Z23_
      Z23_=19.D0/4.D0 + Z23_
      Z23_=Z1_*Z23_
      Z23_=19.D0/4.D0*Z2_ + Z23_
      Z23_=Z4_*Z23_
      Z22_=456*Z22_ + Z23_
      Z22_=Z4_*Z22_
      Z23_= - pi2*Z20_
      Z19_=Z23_ - 4*Z19_ + Z22_
      Z19_=Z7_*Z19_
      Z22_= - 13 - 12*Z37_
      Z22_=Z1_*Z22_
      Z22_=1 + Z22_
      Z22_=Z22_*Z25_
      Z22_= - Z48_ + Z22_
      Z20_=Z4_*Z20_
      Z20_=2*Z22_ - 7*Z20_
      Z20_=zeta3*Z20_
      Z22_= - 139 - 138*Z37_
      Z22_=Z1_*Z22_
      Z22_=1 + Z22_
      Z22_=Z1_*Z22_
      Z22_=Z35_ + Z22_
      Z22_=Z9_*Z22_
      Z23_=Z1_ - 1
      Z23_=Z23_*Z1_
      Z22_= - 17.D0/3.D0*Z23_ + Z22_
      Z44_=209 + 204*Z37_
      Z44_=Z1_*Z44_
      Z44_= - 5 + Z44_
      Z44_=Z1_*Z44_
      Z44_=Z35_ + Z44_
      Z44_=Z8_*Z44_
      Z40_=Z6_*Z40_
      Z48_=605 + 159*Z41_
      Z48_=Z48_*Z17_
      Z48_= - 225 + Z48_
      Z48_=Z1_*Z48_
      Z48_= - 63.D0/4.D0 + Z48_
      Z48_=Z1_*Z48_
      Z28_= - 63.D0/4.D0 - Z28_
      Z28_=Z2_*Z28_
      Z28_=Z28_ + Z48_
      Z48_=3*Z10_
      Z28_=Z28_*Z48_
      Z47_=Z2_*Z47_
      Z47_= - 109.D0/12.D0 - 23*Z47_
      Z47_=Z2_*Z47_
      Z53_=19*Z30_
      Z54_=1097.D0/30.D0 + Z53_
      Z54_=Z1_*Z54_
      Z54_= - 147.D0/20.D0 + Z54_
      Z54_=Z1_*Z54_
      Z54_= - 109.D0/240.D0 + Z54_
      Z54_=Z1_*Z54_
      Z47_=1.D0/20.D0*Z47_ + Z54_
      Z54_=1.D0/4.D0*pi4
      Z47_=Z47_*Z54_
      Z55_= - 67 - 18*Z41_
      Z55_=Z55_*Z17_
      Z55_=21 + Z55_
      Z55_=Z1_*Z55_
      Z55_=5.D0/2.D0 + Z55_
      Z55_=Z1_*Z55_
      Z55_=5.D0/2.D0*Z2_ + Z55_
      Z55_=Z11_*Z55_
      Z19_=21*Z55_ + Z47_ + Z20_ + Z28_ + Z19_ + Z21_ + Z24_ + Z40_ + 4
     & *Z22_ + 5*Z44_
      Z20_=Z35_ - Z23_
      Z20_=Z4_*Z20_
      Z20_=Z20_ + Z2_ - Z1_
      Z20_=Z4_*Z20_
      Z21_=Z35_ + Z23_ + 1
      Z22_=Z21_*Z3_
      Z24_= - Z1_*Z26_
      Z24_= - 25 + Z24_
      Z24_=Z1_*Z24_
      Z24_=23.D0/3.D0 + Z24_
      Z24_=Z24_*Z17_
      Z24_=235.D0/36.D0 + Z24_
      Z24_=Z1_*Z24_
      Z26_=3*Z1_
      Z18_=Z18_*Z26_
      Z18_=Z18_ + 13
      Z28_= - Z18_*Z17_
      Z28_=29.D0/3.D0 + Z28_
      Z28_=Z28_*Z17_
      Z28_= - 3 + Z28_
      Z28_=Z1_*Z28_
      Z40_=4*Z2_
      Z44_= - 5 + Z40_
      Z44_=Z2_*Z44_
      Z44_=4 + Z44_
      Z28_=1.D0/3.D0*Z44_ + Z28_
      Z28_=Z5_*Z28_
      Z44_= - 185.D0/8.D0 + 13*Z2_
      Z44_=Z2_*Z44_
      Z44_=23 + Z44_
      Z30_= - 17.D0/12.D0 - Z30_
      Z30_=Z1_*Z30_
      Z30_=43.D0/24.D0 + Z30_
      Z30_=Z1_*Z30_
      Z30_= - 169.D0/48.D0 + Z30_
      Z30_=Z1_*Z30_
      Z30_=1.D0/6.D0*Z44_ + Z30_
      Z30_=Z4_*Z30_
      Z44_= - 11.D0/12.D0 - Z2_
      Z44_=Z2_*Z44_
      Z44_= - 85.D0/12.D0 + Z44_
      Z24_=1.D0/4.D0*Z30_ - Z22_ + Z28_ + 1.D0/3.D0*Z44_ + Z24_
      Z24_=Z4_*Z24_
      Z28_=Z40_*Z34_
      Z30_= - 19 - 18*Z37_
      Z30_=Z1_*Z30_
      Z30_=1 + Z30_
      Z30_=Z30_*Z25_
      Z44_=2*Z3_
      Z47_=Z21_*Z44_
      Z30_=Z47_ + Z30_ - 1 - Z28_
      Z30_=Z3_*Z30_
      Z47_=2*Z2_
      Z55_=Z47_ - 9
      Z56_=Z55_*Z47_
      Z57_= - 629 + 247*Z1_
      Z57_=Z1_*Z57_
      Z57_=556 + Z57_
      Z57_=Z1_*Z57_
      Z57_= - 1039.D0/4.D0 + 4.D0/3.D0*Z57_
      Z33_=Z57_*Z33_
      Z57_=19.D0/8.D0 + Z2_
      Z57_=Z2_*Z57_
      Z58_= - 517.D0/2.D0 - Z42_
      Z58_=Z1_*Z58_
      Z58_=295.D0/4.D0 + Z58_
      Z58_=Z1_*Z58_
      Z58_= - 29.D0/8.D0 + Z58_
      Z58_=Z1_*Z58_
      Z57_=Z58_ - 1 + Z57_
      Z57_=Z6_*Z57_
      Z58_=197 + 180*Z37_
      Z58_=Z1_*Z58_
      Z58_= - 17 + Z58_
      Z58_=Z1_*Z58_
      Z58_=8 + Z58_
      Z58_=Z5_*Z58_
      Z24_=Z24_ + Z57_ + Z30_ + Z58_ + Z33_ + 1903.D0/72.D0 + Z56_
      Z24_=Z4_*Z24_
      Z30_=8*Z5_
      Z33_=Z21_*Z30_
      Z56_=355 + 318*Z37_
      Z56_=Z1_*Z56_
      Z56_= - 37 + Z56_
      Z56_=Z56_*Z17_
      Z35_= - 4*Z22_ - Z33_ + Z56_ + 19 + 12*Z35_
      Z35_=Z6_*Z35_
      Z56_= - 1627.D0/4.D0 + 553*Z2_
      Z51_=2657.D0/36.D0 + 178*Z51_
      Z51_=Z1_*Z51_
      Z51_=1.D0/18.D0*Z56_ + Z51_
      Z29_=Z29_*Z26_
      Z29_=Z29_ + 5
      Z29_=Z29_*Z1_
      Z29_=Z29_ + 2*Z34_
      Z30_= - Z29_*Z30_
      Z39_=39 + 32*Z39_
      Z39_=Z39_*Z26_
      Z39_=Z39_ - 61.D0/2.D0 + 40*Z2_
      Z56_=Z3_*Z39_
      Z27_= - 1 + Z27_
      Z27_=Z2_*Z27_
      Z57_= - 59 - Z50_
      Z57_=Z1_*Z57_
      Z57_=111 + 8*Z57_
      Z57_=Z1_*Z57_
      Z57_= - 9 + Z57_
      Z57_=Z1_*Z57_
      Z27_=Z57_ + 7 + Z27_
      Z27_=Z27_*Z46_
      Z46_=Z13_*Z21_
      Z57_=445 + 111*Z41_
      Z57_=Z57_*Z17_
      Z57_= - 227 + Z57_
      Z57_=Z57_*Z17_
      Z57_=27 + Z57_
      Z57_=Z1_*Z57_
      Z58_= - 9 - Z40_
      Z58_=Z2_*Z58_
      Z57_=Z57_ - 8 + Z58_
      Z57_=Z8_*Z57_
      Z24_=Z24_ + Z35_ + Z57_ + 8*Z46_ + Z27_ + Z56_ + 1.D0/3.D0*Z51_
     &  + Z30_
      Z24_=Z4_*Z24_
      Z27_= - Z1_*Z38_
      Z27_= - 15 + Z27_
      Z27_=Z27_*Z17_
      Z27_=13.D0/6.D0 + Z27_
      Z27_=Z1_*Z27_
      Z27_=349.D0/36.D0 + Z27_
      Z27_=Z1_*Z27_
      Z30_=Z18_*Z1_
      Z30_=Z30_ - 14.D0/3.D0
      Z30_=Z30_*Z31_
      Z30_=Z30_ + 7.D0/3.D0
      Z30_=Z30_*Z1_
      Z35_=2.D0/3.D0*Z2_
      Z38_=Z35_ - 1
      Z38_=Z38_*Z2_
      Z30_=Z30_ - Z38_ - 2.D0/3.D0
      Z38_= - Z5_*Z30_
      Z43_= - 121.D0/8.D0 + Z43_
      Z43_=Z43_*Z45_
      Z43_=5 + Z43_
      Z45_= - 281.D0/12.D0 - Z49_
      Z45_=Z1_*Z45_
      Z45_=211.D0/24.D0 + Z45_
      Z45_=Z1_*Z45_
      Z45_= - 51.D0/16.D0 + Z45_
      Z45_=Z1_*Z45_
      Z43_=1.D0/2.D0*Z43_ + Z45_
      Z43_=Z43_*Z52_
      Z45_=13.D0/36.D0 - Z2_
      Z45_=Z2_*Z45_
      Z22_=Z43_ + 1.D0/3.D0*Z22_ + Z38_ + Z27_ - 103.D0/36.D0 + Z45_
      Z22_=Z4_*Z22_
      Z27_=12*log2
      Z38_=100 - 41*Z1_
      Z38_=Z1_*Z38_
      Z38_=4.D0/9.D0*Z38_ - 221.D0/9.D0 - Z27_
      Z38_=Z1_*Z38_
      Z27_=Z38_ - 31.D0/12.D0 + Z27_
      Z27_=Z1_*Z27_
      Z36_=Z36_*Z17_
      Z36_=3 + Z36_
      Z36_=Z1_*Z36_
      Z36_= - 1 + Z36_
      Z38_=15*Z1_
      Z36_=Z36_*Z38_
      Z36_=2 + Z36_
      Z36_=Z3_*Z36_
      Z38_=35.D0/8.D0 - Z2_
      Z38_=Z2_*Z38_
      Z38_= - 7 + Z38_
      Z43_=121.D0/2.D0 + Z50_
      Z43_=Z1_*Z43_
      Z43_= - 193.D0/12.D0 + Z43_
      Z43_=Z1_*Z43_
      Z43_=83.D0/24.D0 + Z43_
      Z43_=Z1_*Z43_
      Z38_=1.D0/3.D0*Z38_ + Z43_
      Z38_=Z6_*Z38_
      Z22_=Z22_ + Z38_ + Z36_ + Z27_ - 2*log2 + 15.D0/8.D0 + Z35_
      Z22_=pi2*Z22_
      Z27_=Z32_*Z21_
      Z32_= - 85 - 76*Z37_
      Z32_=Z1_*Z32_
      Z32_=9 + Z32_
      Z32_=Z32_*Z26_
      Z32_=Z27_ + Z32_ - 13 - Z28_
      Z35_=253 + Z42_
      Z35_=Z35_*Z17_
      Z35_= - 133 + Z35_
      Z35_=Z1_*Z35_
      Z35_=23.D0/2.D0 + Z35_
      Z35_=Z1_*Z35_
      Z36_=3.D0/2.D0 - Z40_
      Z36_=Z2_*Z36_
      Z35_=Z35_ - 6 + Z36_
      Z35_=Z4_*Z35_
      Z32_=2*Z32_ + Z35_
      Z32_=Z4_*Z32_
      Z35_=Z21_*Z6_
      Z29_= - Z35_ + Z29_
      Z30_=pi2*Z30_
      Z36_=Z7_*Z21_
      Z29_=4*Z36_ + Z30_ + 8*Z29_ + Z32_
      Z29_=Z7_*Z29_
      Z28_=Z28_ + 9
      Z30_=95 + 84*Z37_
      Z30_=Z1_*Z30_
      Z30_= - 11 + Z30_
      Z30_=Z30_*Z17_
      Z18_=Z18_*Z25_
      Z18_= - 31 + Z18_
      Z18_=Z18_*Z31_
      Z18_=3 + Z18_
      Z18_=Z1_*Z18_
      Z25_=Z2_*Z55_
      Z18_=Z18_ + 10 + Z25_
      Z18_=Z4_*Z18_
      Z18_=Z18_ + Z27_ + Z30_ - Z28_
      Z18_=zeta3*Z18_
      Z25_= - Z33_ + 2*Z28_
      Z27_=313 + 276*Z37_
      Z27_=Z1_*Z27_
      Z27_= - 37 + Z27_
      Z27_=Z1_*Z27_
      Z27_=Z27_ + Z25_
      Z27_=Z9_*Z27_
      Z28_= - 577 - 510*Z37_
      Z28_=Z1_*Z28_
      Z28_=67 + Z28_
      Z28_=Z1_*Z28_
      Z25_=Z28_ - Z25_
      Z25_=Z8_*Z25_
      Z25_=Z27_ + Z25_
      Z27_=2*Z35_
      Z28_=Z27_ - Z39_
      Z28_=Z6_*Z28_
      Z30_=8*Z12_ - 16*Z14_
      Z30_=Z4_*Z30_
      Z30_=8*Z15_ + 16*Z16_ + Z30_
      Z30_=Z21_*Z30_
      Z31_= - 421 - 106*Z41_
      Z31_=Z31_*Z26_
      Z31_=617.D0/2.D0 + Z31_
      Z31_=Z1_*Z31_
      Z31_= - 87.D0/4.D0 + Z31_
      Z31_=Z1_*Z31_
      Z32_=9.D0/4.D0 + Z40_
      Z32_=Z2_*Z32_
      Z31_=Z31_ + 12 + Z32_
      Z31_=Z31_*Z48_
      Z32_= - 1219.D0/8.D0 + 61*Z2_
      Z32_=Z2_*Z32_
      Z32_=181 + Z32_
      Z33_= - 763.D0/20.D0 - Z53_
      Z33_=Z1_*Z33_
      Z33_=749.D0/72.D0 + Z33_
      Z33_=Z1_*Z33_
      Z33_= - 577.D0/240.D0 + Z33_
      Z33_=Z1_*Z33_
      Z32_=1.D0/90.D0*Z32_ + Z33_
      Z32_=Z32_*Z54_
      Z33_=35 + 9*Z41_
      Z33_=Z1_*Z33_
      Z33_= - 325 + 42*Z33_
      Z17_=Z33_*Z17_
      Z17_=57 + Z17_
      Z17_=Z1_*Z17_
      Z33_=13 - 20*Z2_
      Z33_=Z2_*Z33_
      Z17_=Z17_ - 56 + Z33_
      Z17_=Z11_*Z17_
      Z23_=1171.D0/6.D0 + 196*Z23_
      Z17_=Z17_ + Z32_ + Z18_ + Z31_ + Z29_ + Z22_ + Z24_ + Z28_ + 1.D0/
     & 9.D0*Z23_ + Z30_ + 2*Z25_
      Z18_=Z52_ - Z44_
      Z18_=Z21_*Z18_
      Z22_= - 3.D0/2.D0 + Z47_
      Z22_=Z2_*Z22_
      Z23_= - 1.D0/2.D0 + Z1_
      Z23_=Z1_*Z23_
      Z18_=Z23_ + 1 + Z22_ + Z18_
      Z18_=Z4_*Z18_
      Z22_=Z34_ + Z1_
      Z23_= - Z22_*Z44_
      Z18_=Z18_ + Z27_ + Z23_ + Z26_ - 7.D0/2.D0 + Z40_
      Z18_=Z4_*Z18_
      Z23_= - Z4_*Z21_
      Z23_=Z23_ - Z22_
      Z23_=pi2*Z23_
      Z24_=Z6_*Z22_
      Z24_=1 + Z24_
      Z18_=1.D0/6.D0*Z23_ + 2*Z24_ + Z18_
      Z21_=Z21_*Z52_
      Z21_=Z21_ + Z22_
      Z21_=Z4_*Z21_
      Z21_=1.D0/2.D0 + Z21_
      
*      H=u*Z19_ + u**2*Z20_ + u**3*Z17_ + u**4*Z18_ + u**5*Z21_      
*** STATS: original	207P 1646M 547A : 2639
*** STATS: optimized 0P 385M 349A : 734

      f1cur_2L_log2 = Z21_
      f1cur_2L_log1 = Z18_
      f1cur_2L_log0 = Z17_
      f2cur_2L_log1 = Z20_
      f2cur_2L_log0 = Z19_

 999  return
      end
*
*     Eqs. (4.1) and (4.2) of NPB 664 (2003) valid for Q2 --> infinity (Q2=-t)
*     as obtained by Syed starting from the exact mathematica version from Pp
*
      subroutine twoloopff_NPB664_2003_expanded_sh(m2,t,
     +                   f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0,
     +         f2cur_2L_log1,f2cur_2L_log0)
      implicit none
      double precision m2,t
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0,f1cur_2L_log0_pairs
      double complex cone
      double precision q2,L,L2,L3,L4,y
      double precision log2
      double precision zeta2,zeta3
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin

      common/parameters/ame,amu,convfac,alpha,pi      
c      COMMON/CONST/EMUBEAM,EEBEAM,ALPHA,CONVFAC,PI,AME,AMU
      double precision pi2,pi4
*
      integer ifirst
      common/twolcmn4/cone,log2,pi2,pi4,zeta2,zeta3,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         cone= (1.d0,0.d0)
         log2= log(2.d0)
         zeta2=pi*pi/6.d0
         zeta3=1.202056903159594d0
         zeta3=1.202056903159594285399738161511449990764986292d0         
         pi2= pi*pi
         pi4= pi2*pi2
         ifirst = 1
      endif
*****************
      q2 = -t
      y  = m2/q2
      L  = log(q2/m2)
      L2 = L  * L
      L3 = L2 * L
      L4 = L3 * L
*
* coefficient of log^2(lambda/m) for F1
*
      f1cur_2L_log2= (0.d0,0.d0)
     +  -2.d0*y*(1.d0 + (-L)) + (1.d0 + 2.d0*(-L) + L2)/2.d0 + 
     -  y**2*(5.d0 + 5.d0*(-L) + 2.d0*L2) - 
     -  (2.d0*y**3*(25.d0 + 34.d0*(-L) + 12.d0*L2))/3.d0 + 
     -  (y**4*(392.d0 + 549.d0*(-L) + 192.d0*L2))/6.d0
*
*     coefficient of log(lambda/m) for F1
*
      f1cur_2L_log1= (0.d0,0.d0)
     + +(y*(-15.d0 + pi2 - 21.d0*(-L) - 12.d0*L2))/6.d0 + 
     -  (y**3*(-212.d0 + 17.d0*pi2 - 500.d0*(-L) + 12.d0*pi2*(-L) - 
     -       291.d0*L2 - 36.d0*(-L3)))/9.d0 + 
     -  (12.d0 - pi2 + 21.d0*(-L) - pi2*(-L) + 12.d0*L2 + 
     -     3.d0*(-L3))/12.d0 + 
     -  (y**2*(72.d0 - 5.d0*pi2 + 165.d0*(-L) - 4.d0*pi2*(-L) + 
     -       87.d0*L2 + 12.d0*(-L3)))/12.d0 + 
     -  (y**4*(13813.d0 - 1098.d0*pi2 + 33303.d0*(-L)-768.d0*pi2*(-L) + 
     -       19422.d0*L2 + 2304.d0*(-L3)))/144.d0
      f1cur_2L_log1= f1cur_2L_log1 * 2.d0 ! to match the (alpha/2pi)^2
*     
* coefficient not IR (no power of log(lambda/m) for F1 
*
      f1cur_2L_log0= (0.d0,0.d0)
     + + (y**3*(-271469.d0 + 19858*pi2 + 2724.d0*pi4 - 213264.d0*zeta3 + 
     -       6912.d0*pi2*3.d0*log2 + 6620.d0*(-L) + 12732.d0*pi2*(-L) - 
     -       131328.d0*zeta3*(-L) - 61200.d0*L2 + 
     -       3672.d0*pi2*L2 - 6384.d0*(-L3) - 684.d0*L4))/
     -   432. + (y*(-25565.d0 - 3135.d0*pi2 + 96.d0*pi4 - 5400.d0*zeta3+ 
     -       2160.d0*pi2*3.d0*log2 - 12735.d0*(-L) + 900.d0*pi2*(-L) - 
     -       6480.d0*zeta3*(-L) - 8235.d0*L2 + 180.d0*pi2*L2 - 
     -       1530.d0*(-L3) - 45.d0*L4))/2160.d0 + 
     -  (23420.d0 + 1755.d0*pi2 - 2160.d0*log2*pi2 - 177.d0*pi4 - 
     -     9720.d0*zeta3 + 8135.d0*(-L) + 390.d0*pi2*(-L) - 
     -     6480.d0*zeta3*(-L) + 3435.d0*L2 - 90.d0*pi2*L2 + 
     -     930.d0*(-L3) + 135.d0*L4)/4320.d0 + 
     -  (y**2*(177405.d0 - 1515.d0*pi2 - 1509.d0*pi4 + 138240.d0*zeta3 - 
     -       8640.d0*pi2*3.d0*log2 + 39770.d0*(-L) - 10500.d0*pi2*(-L) + 
     -       84240.d0*zeta3*(-L) + 53115.d0*L2 - 
     -       2070.d0*pi2*L2 + 5160.d0*(-L3) + 900.d0*L4))/
     -   2160.d0 + (y**4*(412256945.d0-41373420.d0*pi2-4437792.d0*pi4 + 
     -     316846080.d0*zeta3 - 6635520.d0*pi2*3.d0*log2
     -     - 53549160.d0*(-L) - 
     -       17258400.d0*pi2*(-L) + 199272960.d0*zeta3*(-L) + 
     -       78322500.d0*L2 - 5814720.d0*pi2*L2 + 
     -       8035200.d0*(-L3) + 617760.d0*L4))/103680.d0
      f1cur_2L_log0= f1cur_2L_log0 * 4.d0 ! to match the (alpha/2pi)^2
*     expansion of hep-ph/0406046 instead of (alpha/pi)^2
*
*     coefficient of log(lambda/m) for F2 (expansion by Syed)
*
      f2cur_2L_log1= (0.d0,0.d0)
     + + y*(-L)*(1.d0 + (-L)) - 2.d0*y**2*(1 + 3.d0*(-L) + L2) + 
     -  y**3*(11.d0 + 20.d0*(-L) + 8.d0*L2) - 
     -  (2.d0*y**4*(67.d0 + 116.d0*(-L) + 48.d0*L2))/3.d0
      f2cur_2L_log1= 2.d0*f2cur_2L_log1 !to account for the expansion a/2/pi
*
*     coefficient not IR for F2 (expansion by Syed)
*
      f2cur_2L_log0= (0.d0,0.d0)
     +     + (y*(408.d0 + 39.d0*pi2 - 144.d0*log2*pi2 + 72.d0*zeta3
     +     + 235.d0*(-L) - 30.d0*pi2*(-L) + 3.d0*L2
     +     + 18.d0*(-L3)))/72.d0 + 
     -     (y**4*(-879190.d0 + 94620.d0*pi2 + 34560.d0*log2*pi2
     +     + 9693.d0*pi4 - 763560.d0*zeta3 + 98835.d0*(-L)
     +     + 34740.d0*pi2*(-L) - 453600.d0*zeta3*(-L) - 183825.d0*L2 + 
     -       13770.d0*pi2*L2 - 17550.d0*(-L3) - 135.d0*L4))/
     -     270.d0 + (y**2*(-8045.d0 + 405.d0*pi2 + 1440.d0*log2*pi2
     +     + 66.d0*pi4 - 10080.d0*zeta3 - 2315.d0*(-L)
     +     + 540.d0*pi2*(-L) - 5040.d0*zeta3*(-L)
     +     - 3075.d0*L2 + 120.d0*pi2*L2 - 330.d0*(-L3)-15.d0*L4))/180.d0
     +     + (y**3*(165605.d0 - 14850.d0*pi2 - 11520.d0*log2*pi2
     +     - 1716.d0*pi4 + 161280.d0*zeta3 - 570.d0*(-L)
     +     - 7500.d0*pi2*(-L) + 90720.d0*zeta3*(-L) + 43110.d0*L2 - 
     -       2640.d0*pi2*L2 + 4260.d0*(-L3) + 90.d0*L4))/360.d0
      f2cur_2L_log0= 4.d0*f2cur_2L_log0  !to account for the expansion a/2/pi
c     ALREADY MULTIPLIED BY 4 bY SYED
      return
      end
*******************************************************************************
      subroutine twoloopff_0406046(x,f1cur_1L_epsm1,f1cur_1L_eps0,
     +                         f2cur_1L_eps0)
      implicit none
      double precision x
      double complex f1cur_1L_epsm1,f1cur_1L_eps0
      double complex f2cur_1L_eps0

* for hplog ------------
      double precision omx,opx,omx2,omx3,omx4,omx5,opx2,opx3,opx4,opx5
      integer n1,n2,nw,i1,i2
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
      parameter(n2=1)
      parameter(nw=4)
      double complex Hc1,Hc2,Hc3,Hc4
      double precision     Hr1,Hr2,Hr3,Hr4
      double precision     Hi1,Hi2,Hi3,Hi4
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
*  FOR CHAPLIN
      double complex HPL1,HPL2,HPL3,HPL4,z
      double precision epsz
      double precision rez,imz
      parameter(epsz=1.d-15)
      integer j1,j2,j3,j4
*------------------------------------------------
      double complex cone
      parameter(cone=(1.d0,0.d0))
      double precision log2
      double precision zeta2,zeta3
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin

      common/parameters/ame,amu,convfac,alpha,pi      
c      COMMON/CONST/EMUBEAM,EEBEAM,ALPHA,CONVFAC,PI,AME,AMU
*

      integer ifirst
      common/twolcmn5/zeta2,zeta3,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         zeta2=pi*pi/6.d0
         zeta3=1.202056903159594d0
         zeta3=1.202056903159594285399738161511449990764986292d0
         ifirst = 1
      endif
*     
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
*
****** THIS IS FOR HPLOG ************************************
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,
     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
******THIS IS FOR CHAPLIN **********************************
c      rez= x
c      imz= -epsz
c      z= dcmplx(rez,imz)
c      do j1=-1,1,1
c         Hc1(j1)= HPL1(j1,z)
c         do j2= -1,1,1
c            Hc2(j1,j2)= HPL2(j1,j2,z)
c         enddo
c      enddo
*************************************************************
      omx= 1.d0-x
      opx= 1.d0+x
      omx2= omx*omx
      omx3= omx2*omx
      omx4= omx3*omx
      omx5= omx4*omx
      opx2= opx*opx
      opx3= opx2*opx
      opx4= opx3*opx
      opx5= opx4*opx
*
* RENORMALIZED FORMFACTORS at one loop, Eq. (51) of hep-ph/0406046
*
      f1cur_1L_epsm1= -cone
     +                + Hc1(0)*(cone-cone/(omx)-cone/(opx) )
      f1cur_1L_eps0= 2.d0*cone+Hc1(0)*(cone/2.d0-cone/(opx))
     +             +(cone-cone/(omx)-cone/(opx))
     +      *(zeta2*cone-2.d0*Hc1(0)-Hc2(0,0)+2.d0*Hc2(-1,0))
      f1cur_1L_eps0= -f1cur_1L_eps0
*
* Eq.(52) of hep-ph/0406046
*
      f2cur_1L_eps0= -(cone/(omx)-cone/(opx))*Hc1(0)
*
      return
      end
*
      function kallen(x,y,z)
      double precision kallen
      double precision x,y,z
*
      kallen= x**2+y**2+z**2-2.d0*(x*y+y*z+z*x)
*
      return
      end
*
***************************************************************
**  DILOGARITHM STUFF
** if |z| > 1 calculate li2(1/z) and convert
      double complex  function myli2(zz)
      implicit none
*
      double complex cone,zero,ii
      parameter ( cone=(1d0,0d0),
     1            zero=(0d0,0d0),
     2            ii=  (0d0,1d0) )
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi*pi)
*
      double complex zz
*
      double complex myli22
      external myli22
*
c      if(abs(zz-cone).lt.1.d-8) then
      if(abs(zz-cone).lt.1.d-12) then
         myli2=cone*pi2/6.d0
         return
      endif

      if (abs(zz).le.1.d0) then
         myli2=myli22(zz)
      else  
         myli2=-log(-zz)**2/2d0-cone*pi2/6d0-myli22(cone/zz)
      endif

      end function myli2
*
** if |z|>1/2 calculate li2(1-z) and convert
*
      double complex function myli22(zz)
      implicit none
*
      double complex cone,zero,ii
      parameter ( cone=(1d0,0d0),
     1            zero=(0d0,0d0),
     2            ii=  (0d0,1d0) )
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi*pi)
*
      double complex zz
*
      double complex myli23
      external myli23
*
      if (abs(zz).lt.0.5d0) then
         myli22=myli23(zz)
      else
         myli22=cone*pi2/6d0-log(zz)*log(cone-zz)-myli23(cone-zz)
      endif
*
      end function myli22
*
**  main function for dilogarithms
*
      double complex function myli23(zz)
      implicit none
*
      double complex cone,zero,ii
      parameter ( cone=(1d0,0d0),
     1            zero=(0d0,0d0),
     2            ii=  (0d0,1d0) )
      real * 8 pi,pi2
      parameter (pi=3.141592653589793238462643383279502884197D0,
     1           pi2=pi*pi)
*
      double precision epsilon
      parameter (epsilon=1.d-20)
      double complex zz
*
      integer dim
      parameter (dim = 25)
      double precision b(dim)
      integer i,n
      double complex z,add
*
* even n Bernoulli numbers, already divided by (n + 1)!
*
      data b /
     $  2.777777777777777777777777777777777777777777777777777777778d-2,
     $ -2.777777777777777777777777777777777777777777777777777777778d-4,
     $  4.724111866969009826152683295540438397581254724111866969009d-6,
     $ -9.185773074661963550852439741328630217519106407995296884185d-8,
     $  1.897886998897099907200917301927402937503947604957705967806d-9,
     $ -4.064761645144225526805909386291966674547057127439707822288d-11,
     $  8.921691020456452555217987316752748851514283613049045147810d-13,
     $ -1.993929586072107568723644347793789705630694749653880147036d-14,
     $  4.518980029619918191650476552855593228396819014466618405199d-16,
     $ -1.035651761218124701448341154221865666596091238168650515964d-17,
     $  2.395218621026186745740283743000980381678949001942974256251d-19,
     $ -5.581785874325009336283074505625419905567054667644398095136d-21,
     $  1.309150755418321285812307399186592301749849838783303836854d-22,
     $ -3.087419802426740293242279764866462431595565256132745695326d-24,
     $  7.315975652702203420357905609252148591033401063690875035693d-26,
     $ -1.740845657234000740989055147759702545340841421754271264171d-27,
     $  4.157635644613899719617899620775226673488254159511563860825d-29,
     $ -9.962148488284622103194006702455838849854860017394488768062d-31,
     $  2.394034424896165300521167987893749562934279156932915750221d-32,
     $ -5.768347355367390084291793161877654244072332317926275110062d-34,
     $  1.393179479647007977827886603911548331732411625673399565806d-35,
     $ -3.372121965485089470468473635254930958979742891656539304386d-37,
     $  8.178208777562102621764777214872834267876189462495503276198d-39,
     $ -1.987010831152385925564820669234786567541858995824743201790d-40,
     $  4.835778518040550896287059373115378207694465369420827842732d-42
     $  /
*
      z=-log(cone-zz)

      myli23=z-z**2/4d0

      if (abs(myli23).gt.epsilon) then
          do i=1,dim
             n    = 2*i
             add  = z**(n+1)*b(i)
             if (abs(add/myli23).gt.epsilon) then
                 myli23 = myli23 + add
             else
                 exit
             endif
          enddo
      endif

      end function myli23
********************** CARLO *******************
      subroutine twoloop_NF1(x,m2,t,F12L,F22L)
!      https://arxiv.org/pdf/hep-ph/0405275.pdf
      implicit none
c     implicit double complex (Z)
      double complex F12L,F22L,f1syed,f2syed
      double precision x,m2,t
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0
      double precision q2
* for hplog ------------
      double precision omx,opx,omx2,omx3,omx4,omx5,opx2,opx3,opx4,opx5
      integer n1,n2,nw,i1,i2
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
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
      double complex HPL1,HPL2,HPL3,HPL4,z
      double precision epsz
      double precision rez,imz
      parameter(epsz=1.d-15)
      integer j1,j2,j3,j4
*------------------------------------------------
      double complex cone
      parameter(cone=(1.d0,0.d0))
      double precision log2
      double precision zeta2,zeta3
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin

      common/parameters/ame,amu,convfac,alpha,pi      
      double precision pi2,pi4
*
      integer ifirst
      common/twolcmnNF1/log2,pi2,pi4,zeta2,zeta3,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         log2= log(2.d0)
         zeta2=pi*pi/6.d0
         zeta3=1.202056903159594285399738161511449990764986292d0         
         pi2= pi*pi
         pi4= pi2*pi2
         ifirst = 1
      endif
*     
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
*
****** THIS IS FOR HPLOG ************************************
c      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,
c     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
!     no need to call HPLOGS HERE!
! https://arxiv.org/pdf/hep-ph/0507152.pdf eq.3
!     hc_n(n0) = 1/n! * log^n(x)      
      HC1(0)     = log(x)
      HC2(0,0)   = 0.5d0*HC1(0)*HC1(0)
      HC3(0,0,0) = HC2(0,0)*HC1(0)/3.d0
******THIS IS FOR CHAPLIN **********************************
c      rez =  x
c      imz = -epsz
c      z   =  dcmplx(rez,imz)
c      do j1=-1,1,1
c         Hc1(j1)= HPL1(j1,z)
c         do j2= -1,1,1
c            Hc2(j1,j2)= HPL2(j1,j2,z)
c            do j3= -1,1,1
c               Hc3(j1,j2,j3)= HPL3(j1,j2,j3,z)
c               do j4= -1,1,1
c                  Hc4(j1,j2,j3,j4)= HPL4(j1,j2,j3,j4,z)
c               enddo
c            enddo
c         enddo
c      enddo
*************************************************************
      omx  = 1.d0-x
      opx  = 1.d0+x
      opx2 = opx*opx
      opx3 = opx2*opx
      opx4 = opx3*opx
      opx5 = opx4*opx
*
c      print*,'<<',x
c      print*,HC1(0)-log(x)
c      print*,HC2(0,0)-0.5d0*log(x)**2
c      print*,HC3(0,0,0)-log(x)**3/6.d0
      
      F12L = 383.d0/108.d0 + 49.d0/9.d0/opx2 - 49.d0/9.d0/opx
     .     + (
     .     -0.25d0 + 98.d0/3.d0/opx4 - 196.d0/3.d0/opx3
     .     +229d0/6.d0/opx2 - 5.5d0/opx 
     .     )*zeta2
     .     + (
     .     -265d0/216d0 + 89d0/9.d0/opx3 - 89d0/6.d0/opx2
     .     + 563d0/108d0/opx + 59.d0/27.d0/omx
     .     )*HC1(0)
     .     + (
     .     19.d0/36.d0 + 62.d0/9.d0/opx4 - 124.d0/9.d0/opx3
     .     +163d0/18d0/opx2 - 13.d0/6.d0/opx
     .     )*HC2(0,0)
     .     + (
     .     -1.d0/6.d0 - 6.d0/opx5 + 15.d0/opx4
     .     -11.d0/opx3 + 1.5d0/opx2 + 5.d0/12.d0/opx +5.d0/12.d0/omx
     .     )*( HC3(0,0,0) + zeta2 * HC1(0) )

      F22L = -17.d0/3.d0/opx2 + 17.d0/3.d0/opx
     .     + (
     .     -34.d0/opx4 + 68d0/opx3 - 33.d0/opx2 - 1.d0/opx
     .     ) * zeta2
     .     + (
     .     -31.d0/3.d0/opx3 + 15.5d0/opx2 -59.d0/18.d0/opx
     .      - 17.d0/9.d0/omx
     .     ) * HC1(0)
     .     + (
     .     -22.d0/3.d0/opx4 + 44.d0/3.d0/opx3 - 23.d0/3.d0/opx2
     .     + 1.d0/3.d0/opx
     .     ) * HC2(0,0)
     .     + (
     .     6.d0/opx5 - 15.d0/opx4 + 10.5d0/opx3
     .     - 3.d0/4.d0/opx2 - 3.d0/8.d0/opx - 3.d0/8.d0/omx
     .     ) * (zeta2* HC1(0) + HC3(0,0,0))

      F12L = 4.d0*F12L  ! 4 to agree with PSI and to be consistent with alpha^2/(2*pi)^2
      F22L = 4.d0*F22L

c$$$      f1syed = (947.d0/18.d0 - (16d0*pi**2)/3.d0 - 
c$$$     -    (223d0 + 252d0*pi**2 - 3812d0*x + 216d0*pi**2*x - 8293d0*x**2+ 
c$$$     -       4380d0*pi**2*x**2 + 8293d0*x**4 - 4380d0*pi**2*x**4 + 
c$$$     -  3812d0*x**5 - 216d0*pi**2*x**5 - 223d0*x**6 - 252d0*pi**2*x**6 + 
c$$$     -   4d0*((1d0 + x)**2*(265d0 - 124d0*x + 3494d0*x**2 - 124d0*x**3 + 
c$$$     -             265d0*x**4) + 
c$$$     -      6*pi**2*(1d0 + 4d0*x + x**2 + 68d0*x**3 + x**4 + 4d0*x**5 + 
c$$$     -             x**6))*HC1(0) - 
c$$$     -     24d0*(-19d0 + 2d0*x - 187d0*x**2 + 187d0*x**4 - 2d0*x**5
c$$$     - + 19d0*x**6)*
c$$$     -        HC2(0,0) + 144d0*HC3(0,0,0) + 576d0*x*HC3(0,0,0) + 
c$$$     -       144d0*x**2*HC3(0,0,0) + 9792d0*x**3*HC3(0,0,0) + 
c$$$     -       144d0*x**4*HC3(0,0,0) + 576d0*x**5*HC3(0,0,0) + 
c$$$     -       144d0*x**6*HC3(0,0,0))/(54d0*(-1d0 + x)*(1d0 + x)**5))/16d0
c$$$
c$$$      f2syed =
c$$$     -     -(x*((25d0 - 272d0*x - 18d0*(33d0 + 2d0*pi**2)*x**2 
c$$$     - - 272d0*x**3 + 25d0*x**4)*
c$$$     -        HC1(0) + 3d0*((-1d0 + x**2)*
c$$$     -           (-34d0*(1d0 + x)**2 + pi**2*(1d0 + 36d0*x + x**2)) - 
c$$$     -          2d0*(-1d0 + 20d0*x - 20d0*x**3 + x**4)*HC2(0,0) - 
c$$$     -          72d0*x**2*HC3(0,0,0))))/(18.d0*(-1d0 + x)*(1d0 + x)**5)
c$$$      f1syed = 4.d0*f1syed
c$$$      f2syed = 4.d0*f2syed
c$$$      print*,'x,F1, F2',x,(F12L-f1syed)/F12L,(F22L-f2syed)/F22L
      
      return
      end
***********************************************************************
********************** TEST TIMELIKE *******************
      subroutine twoloop_NF1_TL(s,m2,t,F12L,F22L)
!      https://arxiv.org/pdf/hep-ph/0405275.pdf
      implicit none
c     implicit double complex (Z)
      double complex F12L,F22L,f1syed,f2syed
      double precision s
      double complex x
      double precision m2,t
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0
      double precision q2
* for hplog ------------
      double precision omx,opx,omx2,omx3,omx4,omx5,opx2,opx3,opx4,opx5
      integer n1,n2,nw,i1,i2
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
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
      double complex HPL1,HPL2,HPL3,HPL4,z
      double precision epsz
      double precision rez,imz
      parameter(epsz=1.d-15)
      integer j1,j2,j3,j4
*------------------------------------------------
      double complex cone
      parameter(cone=(1.d0,0.d0))
      double precision log2
      double precision zeta2,zeta3
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin

      common/parameters/ame,amu,convfac,alpha,pi      
      double precision pi2,pi4
*
      integer ifirst
      common/twolcmnNF1/log2,pi2,pi4,zeta2,zeta3,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         log2= log(2.d0)
         zeta2=pi*pi/6.d0
         zeta3=1.202056903159594285399738161511449990764986292d0         
         pi2= pi*pi
         pi4= pi2*pi2
         ifirst = 1
      endif
*     
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046

      x = (sqrt(-s*cone+4.d0*m2)-sqrt(-s*cone))/
     .     (sqrt(-s*cone+4.d0*m2) + sqrt(-s*cone))

      
****** THIS IS FOR HPLOG ************************************
c      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,
c     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
!     no need to call HPLOGS HERE!
! https://arxiv.org/pdf/hep-ph/0507152.pdf eq.3
!     hc_n(n0) = 1/n! * log^n(x)      
      HC1(0)     = log(x)
      HC2(0,0)   = 0.5d0*HC1(0)*HC1(0)
      HC3(0,0,0) = HC2(0,0)*HC1(0)/3.d0
******THIS IS FOR CHAPLIN **********************************
c      rez =  x
c      imz = -epsz
c      z   =  dcmplx(rez,imz)
c      do j1=-1,1,1
c         Hc1(j1)= HPL1(j1,z)
c         do j2= -1,1,1
c            Hc2(j1,j2)= HPL2(j1,j2,z)
c            do j3= -1,1,1
c               Hc3(j1,j2,j3)= HPL3(j1,j2,j3,z)
c               do j4= -1,1,1
c                  Hc4(j1,j2,j3,j4)= HPL4(j1,j2,j3,j4,z)
c               enddo
c            enddo
c         enddo
c      enddo
*************************************************************
      omx  = 1.d0-x
      opx  = 1.d0+x
      opx2 = opx*opx
      opx3 = opx2*opx
      opx4 = opx3*opx
      opx5 = opx4*opx
*
c      print*,'<<',x
c      print*,HC1(0)-log(x)
c      print*,HC2(0,0)-0.5d0*log(x)**2
c      print*,HC3(0,0,0)-log(x)**3/6.d0
      
      F12L = 383.d0/108.d0 + 49.d0/9.d0/opx2 - 49.d0/9.d0/opx
     .     + (
     .     -0.25d0 + 98.d0/3.d0/opx4 - 196.d0/3.d0/opx3
     .     +229d0/6.d0/opx2 - 5.5d0/opx 
     .     )*zeta2
     .     + (
     .     -265d0/216d0 + 89d0/9.d0/opx3 - 89d0/6.d0/opx2
     .     + 563d0/108d0/opx + 59.d0/27.d0/omx
     .     )*HC1(0)
     .     + (
     .     19.d0/36.d0 + 62.d0/9.d0/opx4 - 124.d0/9.d0/opx3
     .     +163d0/18d0/opx2 - 13.d0/6.d0/opx
     .     )*HC2(0,0)
     .     + (
     .     -1.d0/6.d0 - 6.d0/opx5 + 15.d0/opx4
     .     -11.d0/opx3 + 1.5d0/opx2 + 5.d0/12.d0/opx +5.d0/12.d0/omx
     .     )*( HC3(0,0,0) + zeta2 * HC1(0) )

      F22L = -17.d0/3.d0/opx2 + 17.d0/3.d0/opx
     .     + (
     .     -34.d0/opx4 + 68d0/opx3 - 33.d0/opx2 - 1.d0/opx
     .     ) * zeta2
     .     + (
     .     -31.d0/3.d0/opx3 + 15.5d0/opx2 -59.d0/18.d0/opx
     .      - 17.d0/9.d0/omx
     .     ) * HC1(0)
     .     + (
     .     -22.d0/3.d0/opx4 + 44.d0/3.d0/opx3 - 23.d0/3.d0/opx2
     .     + 1.d0/3.d0/opx
     .     ) * HC2(0,0)
     .     + (
     .     6.d0/opx5 - 15.d0/opx4 + 10.5d0/opx3
     .     - 3.d0/4.d0/opx2 - 3.d0/8.d0/opx - 3.d0/8.d0/omx
     .     ) * (zeta2* HC1(0) + HC3(0,0,0))

      F12L = 4.d0*F12L  ! 4 to agree with PSI and to be consistent with alpha^2/(2*pi)^2
      F22L = 4.d0*F22L      
      return
      end      
*************************************************************************
      subroutine twoloop_NF1_expanded(x,m2,t,F12L,F22L)
!      https://arxiv.org/pdf/hep-ph/0405275.pdf
      implicit none
c     implicit double complex (Z)
      double complex F12L,F22L
      double precision x,m2,t
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0
      double precision q2,xapprox
* for hplog ------------
      double precision omx,opx,omx2,omx3,omx4,omx5,opx2,opx3,opx4,opx5
      integer n1,n2,nw,i1,i2
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
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
      double complex HPL1,HPL2,HPL3,HPL4,z
      double precision epsz
      double precision rez,imz
      parameter(epsz=1.d-15)
      integer j1,j2,j3,j4
*------------------------------------------------
      double complex cone
      parameter(cone=(1.d0,0.d0))
      double precision log2
      double precision zeta2,zeta3
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin

      common/parameters/ame,amu,convfac,alpha,pi      
      double precision pi2,pi4
*
      integer ifirst
      common/twolcmnNF1/log2,pi2,pi4,zeta2,zeta3,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         log2= log(2.d0)
         zeta2=pi*pi/6.d0
         zeta3=1.202056903159594285399738161511449990764986292d0         
         pi2= pi*pi
         pi4= pi2*pi2
         ifirst = 1
      endif
*     
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
*
****** THIS IS FOR HPLOG ************************************
c      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,
c     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
!     no need to call HPLOGS HERE!

!     hc_n(n0) = 1/n! * log^n(x)

      xapprox = -m2/t
      
      HC1(0)     = log(xapprox)
      HC2(0,0)   = 0.5d0*HC1(0)*HC1(0)
      HC3(0,0,0) = HC2(0,0)*HC1(0)/3.d0
     
******THIS IS FOR CHAPLIN **********************************
c      rez =  x
c      imz = -epsz
c      z   =  dcmplx(rez,imz)
c      do j1=-1,1,1
c         Hc1(j1)= HPL1(j1,z)
c         do j2= -1,1,1
c            Hc2(j1,j2)= HPL2(j1,j2,z)
c            do j3= -1,1,1
c               Hc3(j1,j2,j3)= HPL3(j1,j2,j3,z)
c               do j4= -1,1,1
c                  Hc4(j1,j2,j3,j4)= HPL4(j1,j2,j3,j4,z)
c               enddo
c            enddo
c         enddo
c      enddo
*************************************************************
      omx  = 1.d0
      opx  = 1.d0
      opx2 = opx*opx
      opx3 = opx2*opx
      opx4 = opx3*opx
      opx5 = opx4*opx
*      
      F12L = 383.d0/108.d0 + 49.d0/9.d0/opx2 - 49.d0/9.d0/opx
     .     + (
     .     -0.25d0 + 98.d0/3.d0/opx4 - 196.d0/3.d0/opx3
     .     +229d0/6.d0/opx2 - 5.5d0/opx 
     .     )*zeta2
     .     + (
     .     -265d0/216d0 + 89d0/9.d0/opx3 - 89d0/6.d0/opx2
     .     + 563d0/108d0/opx + 59.d0/27.d0/omx
     .     )*HC1(0)
     .     + (
     .     19.d0/36.d0 + 62.d0/9.d0/opx4 - 124.d0/9.d0/opx3
     .     +163d0/18d0/opx2 - 13.d0/6.d0/opx
     .     )*HC2(0,0)
     .     + (
     .     -1.d0/6.d0 - 6.d0/opx5 + 15.d0/opx4
     .     -11.d0/opx3 + 1.5d0/opx2 + 5.d0/12.d0/opx +5.d0/12.d0/omx
     .     )*( HC3(0,0,0) + zeta2 * HC1(0) )

      F22L = -17.d0/3.d0/opx2 + 17.d0/3.d0/opx
     .     + (
     .     -34.d0/opx4 + 68d0/opx3 - 33.d0/opx2 - 1.d0/opx
     .     ) * zeta2
     .     + (
     .     -31.d0/3.d0/opx3 + 15.5d0/opx2 -59.d0/18.d0/opx
     .      - 17.d0/9.d0/omx
     .     ) * HC1(0)
     .     + (
     .     -22.d0/3.d0/opx4 + 44.d0/3.d0/opx3 - 23.d0/3.d0/opx2
     .     + 1.d0/3.d0/opx
     .     ) * HC2(0,0)
     .     + (
     .     6.d0/opx5 - 15.d0/opx4 + 10.5d0/opx3
     .     - 3.d0/4.d0/opx2 - 3.d0/8.d0/opx - 3.d0/8.d0/omx
     .     ) * (zeta2* HC1(0) + HC3(0,0,0))

      F12L = 4.d0*F12L ! 4 to agree with PSI and to be consistent with alpha^2/(2*pi)^2
      F22L = 4.d0*F22L 
      return
      end
****************************************************************************      
*******************************************************************************
* Equations 4.1 and 4.2 provided by Pierpaolo Mastrolia
      subroutine twoloopff_Pp_Carlo(x,m2,t,
     +                   f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0,
     +                                 f2cur_2L_log1,f2cur_2L_log0)
c     implicit none
      implicit double precision (Z)
      double precision x,m2,t
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0
      double precision q2
* for hplog ------------
      double precision omx,opx,omx2,omx3,omx4,omx5,opx2,opx3,opx4,opx5
      double precision omxm1,opxm1
      integer n1,n2,nw,i1,i2
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
      parameter(n2=1)
      parameter(nw=4)
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
      double complex HPL1,HPL2,HPL3,HPL4,z
      double precision epsz
      double precision rez,imz
      parameter(epsz=1.d-15)
      integer j1,j2,j3,j4
*------------------------------------------------
      double complex cone
      parameter(cone=(1.d0,0.d0))
      double precision ln2
      double precision z2,z3
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin

      common/parameters/ame,amu,convfac,alpha,pi      
      double precision pi2,pi4
      integer ifirst
      double precision cover
      common/twolcmn3Carlo/ln2,pi2,pi4,z2,z3,cover,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         ln2 = log(2.d0)
         z2  = pi*pi/6.d0
         z3  = 1.202056903159594d0
         z3  = 1.202056903159594285399738161511449990764986292d0
         pi2 = pi*pi
         pi4 = pi2*pi2
         cover = 1.d0 / 34560d0
         ifirst = 1
      endif
*     
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
*
****** THIS IS FOR HPLOG ************************************
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,
     +     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
******THIS IS FOR CHAPLIN **********************************
c      rez =  x
c      imz = -epsz
c      z   =  dcmplx(rez,imz)
c      do j1=-1,1,1
c         Hc1(j1)= HPL1(j1,z)
c         do j2= -1,1,1
c            Hc2(j1,j2)= HPL2(j1,j2,z)
c            do j3= -1,1,1
c               Hc3(j1,j2,j3)= HPL3(j1,j2,j3,z)
c               do j4= -1,1,1
c                  Hc4(j1,j2,j3,j4)= HPL4(j1,j2,j3,j4,z)
c               enddo
c            enddo
c         enddo
c      enddo
*************************************************************
c      omx= 1.d0-x
c      opx= 1.d0+x

      omxm1 = 1.d0/(1.d0-x)
      opxm1 = 1.d0/(1.d0+x)

      Z1_=HC1(-1)
      Z2_=HC1(0)
      Z3_=HC2(0,-1)
      Z4_=HC1(1)
      Z5_=HC2(0,1)
      Z6_=HC3(0,0,-1)
      Z7_=HC3(0,0,1)
      Z8_=HC4(0,0,0,-1)
      Z9_=HC4(0,0,0,1)
      Z10_=HC3(0,-1,1)
      Z11_=HC3(0,1,-1)
      Z12_=HC3(0,1,1)
      Z13_=HC4(0,-1,0,1)
      Z14_=HC4(0,0,1,1)
      Z15_=opxm1 - 1
      Z16_=Z15_*opxm1
      Z17_=omxm1 - 1
      Z18_=Z17_*omxm1
      Z19_=Z18_ + Z16_ + 1
      Z20_=4*Z1_
      Z21_=Z2_ - Z20_
      Z21_=Z19_*Z21_
      Z22_=4*omxm1
      Z23_= - 3 + Z22_
      Z23_=omxm1*Z23_
      Z24_=2*opxm1
      Z25_= - 1 + Z24_
      Z25_=opxm1*Z25_
      Z21_=Z25_ + 2 + Z23_ + Z21_
      Z21_=Z2_*Z21_
      Z23_=Z17_ + opxm1
      Z25_= - Z23_*Z20_
      Z26_=8*omxm1
      Z27_=6*opxm1
      Z21_=Z21_ + Z25_ + Z27_ - 7 + Z26_
      Z21_=Z2_*Z21_
      Z21_=4 + Z21_
      Z25_=Z19_*Z2_
      Z28_=Z25_ + Z23_
      Z29_=12*Z3_ - pi2
      Z28_=Z28_*Z29_
      Z21_=3*Z21_ + Z28_
      Z21_=2880*Z21_
      Z23_=2*Z23_ + Z25_
      Z23_=Z2_*Z23_
      Z23_=1 + Z23_
      Z23_=17280*Z23_
      Z25_=4*opxm1
      Z28_=211 - 83*opxm1
      Z28_=Z28_*Z25_
      Z28_= - 697 + Z28_
      Z28_=Z28_*Z25_
      Z28_=731 + Z28_
      Z28_=opxm1*Z28_
      Z29_=5*Z18_
      Z30_=opxm1 - 2
      Z31_=Z30_*opxm1
      Z32_=253 + 252*Z31_
      Z32_=opxm1*Z32_
      Z32_= - 1 + Z32_
      Z32_=opxm1*Z32_
      Z32_=Z29_ + Z32_
      Z33_=6*Z1_
      Z32_=Z32_*Z33_
      Z34_= - 89 - 90*Z31_
      Z34_=opxm1*Z34_
      Z34_= - 1 + Z34_
      Z34_=opxm1*Z34_
      Z34_= - Z18_ + Z34_
      Z34_=Z4_*Z34_
      Z35_=Z22_ - 5
      Z35_=Z35_*omxm1
      Z28_=24*Z34_ + Z32_ - 9*Z35_ + Z28_
      Z32_=Z24_ - 5
      Z34_=Z32_*Z27_
      Z36_=Z34_ + 25
      Z36_=Z36_*Z24_
      Z36_=Z36_ - 15
      Z36_=Z36_*Z24_
      Z36_=Z36_ + 1
      Z36_=Z36_*opxm1
      Z36_=Z36_ + omxm1
      Z37_=12*Z4_
      Z38_=Z36_*Z37_
      Z39_=7*opxm1
      Z40_=Z39_ - 23
      Z41_=12*opxm1
      Z42_=Z40_*Z41_
      Z42_=293 + Z42_
      Z42_=opxm1*Z42_
      Z42_= - 88 + Z42_
      Z42_=Z42_*Z25_
      Z42_= - 89 + Z42_
      Z42_=opxm1*Z42_
      Z43_=63 - 26*omxm1
      Z43_=omxm1*Z43_
      Z42_=Z38_ + Z43_ + Z42_
      Z43_=3*opxm1
      Z44_=Z32_*Z43_
      Z45_=8 + Z44_
      Z45_=opxm1*Z45_
      Z45_=3 + Z45_
      Z46_=8*opxm1
      Z45_=Z45_*Z46_
      Z45_= - 11 + Z45_
      Z45_=opxm1*Z45_
      Z47_=2*omxm1
      Z48_=Z47_ - 3
      Z49_=Z48_*omxm1
      Z50_=6*Z49_
      Z51_= - 11 - Z50_
      Z51_=omxm1*Z51_
      Z45_=Z51_ + Z45_
      Z45_=Z2_*Z45_
      Z42_=8*Z42_ + Z45_
      Z42_=Z2_*Z42_
      Z28_=16*Z28_ + Z42_
      Z42_=3*Z2_
      Z28_=Z28_*Z42_
      Z45_=Z24_ - 3
      Z51_=Z45_*opxm1
      Z52_=89 + 96*Z51_
      Z52_=Z52_*opxm1
      Z53_=7*omxm1
      Z52_=Z52_ + Z53_
      Z54_=Z1_*Z52_
      Z55_=Z45_*Z27_
      Z55_=Z55_ + 7
      Z55_=Z55_*opxm1
      Z55_=Z55_ - omxm1
      Z56_=Z4_*Z55_
      Z57_= - 607 - 744*Z51_
      Z57_=opxm1*Z57_
      Z54_=144*Z56_ - 36*Z54_ - 137*omxm1 + Z57_
      Z28_=16*Z54_ + Z28_
      Z28_=Z2_*Z28_
      Z48_=Z48_*Z47_
      Z54_= - 63 - Z48_
      Z54_=omxm1*Z54_
      Z32_=Z32_*opxm1
      Z56_=605 + 159*Z32_
      Z56_=Z56_*Z24_
      Z56_= - 225 + Z56_
      Z56_=Z56_*Z25_
      Z56_= - 63 + Z56_
      Z56_=opxm1*Z56_
      Z54_=Z54_ + Z56_
      Z54_=Z8_*Z54_
      Z56_=126*Z32_
      Z57_= - 485 - Z56_
      Z57_=Z57_*Z24_
      Z57_=195 + Z57_
      Z57_=Z57_*Z24_
      Z57_=19 + Z57_
      Z57_=opxm1*Z57_
      Z57_=19*omxm1 + Z57_
      Z57_=Z2_*Z57_
      Z58_=1 + Z31_
      Z58_=Z58_*opxm1**2
      Z57_=1824*Z58_ + Z57_
      Z57_=Z2_*Z57_
      Z55_= - 16*Z55_ + Z57_
      Z55_=Z5_*Z55_
      Z28_=144*Z55_ + 432*Z54_ - 13056*Z16_ + Z28_
      Z54_=209 + 204*Z31_
      Z54_=opxm1*Z54_
      Z54_= - 5 + Z54_
      Z54_=opxm1*Z54_
      Z54_=Z18_ + Z54_
      Z55_= - 853 - 222*Z32_
      Z55_=Z55_*Z24_
      Z55_=339 + Z55_
      Z55_=Z55_*Z24_
      Z55_=35 + Z55_
      Z55_=opxm1*Z55_
      Z55_=35*omxm1 + Z55_
      Z55_=Z2_*Z55_
      Z54_=10*Z54_ + Z55_
      Z54_=Z6_*Z54_
      Z55_= - 13 - 12*Z31_
      Z55_=opxm1*Z55_
      Z55_=1 + Z55_
      Z55_=Z55_*Z39_
      Z55_= - Z29_ + Z55_
      Z57_=Z2_*Z36_
      Z55_=4*Z55_ - 7*Z57_
      Z55_=z3*Z55_
      Z54_=Z54_ + Z55_
      Z55_= - 32 + 13*opxm1
      Z55_=Z55_*Z24_
      Z55_=43 + Z55_
      Z55_=Z55_*Z46_
      Z55_= - 25 + Z55_
      Z55_=opxm1*Z55_
      Z57_= - 77 - 60*Z31_
      Z57_=opxm1*Z57_
      Z57_=17 + Z57_
      Z57_=opxm1*Z57_
      Z57_=3*Z18_ + Z57_
      Z33_=Z57_*Z33_
      Z57_=ln2*Z16_
      Z33_=96*Z57_ + Z33_ - 15*omxm1 + Z55_
      Z55_=Z39_ - 18
      Z57_=Z55_*Z41_
      Z57_=173 + Z57_
      Z57_=opxm1*Z57_
      Z57_= - 14 + Z57_
      Z57_=Z57_*Z25_
      Z57_= - 145 + Z57_
      Z57_=opxm1*Z57_
      Z58_=71 - 34*omxm1
      Z58_=omxm1*Z58_
      Z38_=Z38_ + Z58_ + Z57_
      Z57_=45 + 11*Z32_
      Z57_=Z57_*Z24_
      Z57_= - 25 + Z57_
      Z57_=Z57_*Z41_
      Z57_= - 1 + Z57_
      Z57_=opxm1*Z57_
      Z58_= - 1 - 14*Z49_
      Z58_=omxm1*Z58_
      Z57_=Z58_ + Z57_
      Z57_=Z2_*Z57_
      Z38_=4*Z38_ + Z57_
      Z38_=Z2_*Z38_
      Z36_=Z5_*Z36_
      Z33_= - 48*Z36_ + 8*Z33_ + Z38_
      Z36_=1097 + 285*Z32_
      Z36_=Z36_*Z24_
      Z36_= - 441 + Z36_
      Z36_=Z36_*Z25_
      Z36_= - 109 + Z36_
      Z36_=opxm1*Z36_
      Z38_= - 109 - 46*Z49_
      Z38_=omxm1*Z38_
      Z36_=Z38_ + Z36_
      Z36_=pi2*Z36_
      Z33_=10*Z33_ + Z36_
      Z36_=3*pi2
      Z33_=Z33_*Z36_
      Z38_= - 7 + Z48_
      Z38_=omxm1*Z38_
      Z48_=63*Z32_
      Z57_=248 + Z48_
      Z57_=opxm1*Z57_
      Z57_= - 57 + Z57_
      Z57_=Z57_*Z46_
      Z57_= - 7 + Z57_
      Z57_=opxm1*Z57_
      Z38_=Z38_ + Z57_
      Z38_=Z2_*Z38_
      Z57_= - 649 - 636*Z31_
      Z57_=opxm1*Z57_
      Z57_=13 + Z57_
      Z57_=opxm1*Z57_
      Z29_= - Z29_ + Z57_
      Z29_=8*Z29_ + Z38_
      Z29_=Z2_*Z29_
      Z38_=15*Z32_
      Z57_= - 58 - Z38_
      Z57_=opxm1*Z57_
      Z57_=12 + Z57_
      Z57_=Z57_*Z46_
      Z57_=7 + Z57_
      Z57_=opxm1*Z57_
      Z50_=7 + Z50_
      Z50_=omxm1*Z50_
      Z50_=Z50_ + Z57_
      Z50_=pi2*Z50_
      Z29_=Z50_ + 8*Z52_ + Z29_
      Z50_=360*Z3_
      Z29_=Z29_*Z50_
      Z52_=30*Z32_
      Z57_=113 + Z52_
      Z57_=Z57_*Z24_
      Z57_= - 39 + Z57_
      Z57_=Z57_*Z24_
      Z57_= - 7 + Z57_
      Z57_=opxm1*Z57_
      Z53_= - Z53_ + Z57_
      Z53_=Z2_*Z53_
      Z57_= - 139 - 138*Z31_
      Z57_=opxm1*Z57_
      Z57_=1 + Z57_
      Z57_=opxm1*Z57_
      Z53_=Z53_ + Z18_ + Z57_
      Z53_=Z7_*Z53_
      Z57_= - 67 - 18*Z32_
      Z57_=Z57_*Z24_
      Z57_=21 + Z57_
      Z57_=Z57_*Z24_
      Z57_=5 + Z57_
      Z57_=opxm1*Z57_
      Z57_=5*omxm1 + Z57_
      Z57_=Z9_*Z57_
      Z28_=11520*Z53_ + 30240*Z57_ + Z29_ + 5*Z28_ + Z33_ + 1440*Z54_
      Z28_=3*Z28_
      Z29_=Z18_ - Z16_
      Z29_=Z2_*Z29_
      Z29_=Z29_ + omxm1 - opxm1
      Z29_=17280*Z2_*Z29_
      Z33_=Z19_*Z4_
      Z33_=4*Z33_
      Z53_=Z22_*Z17_
      Z54_=Z53_ + 9
      Z57_=95 + 84*Z31_
      Z57_=opxm1*Z57_
      Z57_= - 11 + Z57_
      Z57_=Z57_*Z24_
      Z44_=Z44_ + 13
      Z58_=Z44_*Z39_
      Z58_= - 31 + Z58_
      Z58_=Z58_*Z25_
      Z58_=3 + Z58_
      Z58_=opxm1*Z58_
      Z47_=Z47_ - 9
      Z47_=Z47_*omxm1
      Z58_=Z58_ + 10 + Z47_
      Z58_=Z2_*Z58_
      Z57_=Z58_ + Z33_ + Z57_ - Z54_
      Z57_=z3*Z57_
      Z58_=8*Z4_
      Z59_=Z58_*Z19_
      Z54_= - Z59_ + 2*Z54_
      Z59_= - 577 - 510*Z31_
      Z59_=opxm1*Z59_
      Z59_=67 + Z59_
      Z59_=opxm1*Z59_
      Z59_=Z59_ - Z54_
      Z60_=445 + 111*Z32_
      Z60_=Z60_*Z24_
      Z60_= - 227 + Z60_
      Z60_=Z60_*Z24_
      Z60_=27 + Z60_
      Z60_=opxm1*Z60_
      Z22_= - 9 - Z22_
      Z22_=omxm1*Z22_
      Z22_=Z60_ - 8 + Z22_
      Z22_=Z2_*Z22_
      Z22_=2*Z59_ + Z22_
      Z22_=Z6_*Z22_
      Z59_=35 + 9*Z32_
      Z59_=opxm1*Z59_
      Z59_= - 325 + 42*Z59_
      Z59_=Z59_*Z24_
      Z59_=57 + Z59_
      Z59_=opxm1*Z59_
      Z60_=13 - 20*omxm1
      Z60_=omxm1*Z60_
      Z59_=Z59_ - 56 + Z60_
      Z59_=Z9_*Z59_
      Z22_=Z59_ + Z22_ + Z57_
      Z57_=Z19_*Z1_
      Z59_= - 19 - 18*Z31_
      Z59_=opxm1*Z59_
      Z59_=1 + Z59_
      Z39_=Z59_*Z39_
      Z39_=2*Z57_ + Z39_ - 1 - Z53_
      Z59_=72*Z1_
      Z39_=Z39_*Z59_
      Z60_= - 629 + 247*opxm1
      Z60_=opxm1*Z60_
      Z60_=556 + Z60_
      Z60_=opxm1*Z60_
      Z60_= - 3117 + 16*Z60_
      Z60_=Z60_*Z24_
      Z61_=197 + 180*Z31_
      Z61_=opxm1*Z61_
      Z61_= - 17 + Z61_
      Z61_=opxm1*Z61_
      Z61_=8 + Z61_
      Z61_=Z4_*Z61_
      Z39_=72*Z61_ + Z39_ + Z60_ + 1903 + 144*Z47_
      Z40_= - opxm1*Z40_
      Z40_= - 25 + Z40_
      Z40_=Z40_*Z43_
      Z40_=23 + Z40_
      Z40_=opxm1*Z40_
      Z40_=235 + 24*Z40_
      Z40_=opxm1*Z40_
      Z47_= - Z44_*Z27_
      Z47_=29 + Z47_
      Z47_=Z47_*Z24_
      Z47_= - 9 + Z47_
      Z47_=opxm1*Z47_
      Z35_=Z47_ + 4 + Z35_
      Z35_=Z35_*Z37_
      Z47_= - 11 - 12*omxm1
      Z47_=omxm1*Z47_
      Z35_=Z35_ - 36*Z57_ + Z40_ - 85 + Z47_
      Z34_= - 17 - Z34_
      Z34_=Z34_*Z24_
      Z34_=43 + Z34_
      Z34_=Z34_*Z24_
      Z34_= - 169 + Z34_
      Z34_=opxm1*Z34_
      Z40_= - 185 + 104*omxm1
      Z40_=omxm1*Z40_
      Z34_=Z34_ + 184 + Z40_
      Z34_=Z34_*Z42_
      Z34_=16*Z35_ + Z34_
      Z34_=Z2_*Z34_
      Z34_=8*Z39_ + Z34_
      Z34_=Z34_*Z42_
      Z35_=2657 + 2136*Z51_
      Z35_=Z35_*Z24_
      Z39_=39 + 32*Z51_
      Z39_=Z39_*Z27_
      Z39_=Z39_ - 61 + 80*omxm1
      Z40_=Z1_*Z39_
      Z45_=Z45_*Z43_
      Z45_=Z45_ + 5
      Z45_=Z45_*opxm1
      Z17_=Z45_ + 2*Z17_
      Z45_=Z4_*Z17_
      Z35_= - 1728*Z45_ + 108*Z40_ + Z35_ - 1627 + 2212*omxm1
      Z40_=Z12_*Z19_
      Z34_= - 27648*Z40_ + 8*Z35_ + Z34_
      Z34_=Z2_*Z34_
      Z35_= - 85 - 76*Z31_
      Z35_=opxm1*Z35_
      Z35_=9 + Z35_
      Z35_=Z35_*Z43_
      Z33_=Z33_ + Z35_ - 13 - Z53_
      Z35_=253 + Z48_
      Z35_=Z35_*Z24_
      Z35_= - 133 + Z35_
      Z35_=Z35_*Z24_
      Z35_=23 + Z35_
      Z35_=opxm1*Z35_
      Z40_=3 - Z26_
      Z40_=omxm1*Z40_
      Z35_=Z35_ - 12 + Z40_
      Z35_=Z2_*Z35_
      Z33_=4*Z33_ + Z35_
      Z33_=Z2_*Z33_
      Z35_=Z19_*Z5_
      Z17_=8*Z35_ + 16*Z17_ + Z33_
      Z17_=Z5_*Z17_
      Z33_= - 421 - 106*Z32_
      Z33_=Z33_*Z27_
      Z33_=617 + Z33_
      Z33_=Z33_*Z24_
      Z33_= - 87 + Z33_
      Z33_=opxm1*Z33_
      Z40_=16*omxm1
      Z45_=9 + Z40_
      Z45_=omxm1*Z45_
      Z33_=Z33_ + 48 + Z45_
      Z33_=Z8_*Z33_
      Z16_=1171 + 1176*Z16_
      Z16_=864*Z17_ + 1296*Z33_ + 32*Z16_ + Z34_
      Z17_=100 - 41*opxm1
      Z17_=Z17_*Z25_
      Z17_= - 221 + Z17_
      Z17_=Z17_*Z25_
      Z17_= - 93 + Z17_
      Z17_=Z17_*Z24_
      Z30_=Z30_*Z24_
      Z30_=3 + Z30_
      Z30_=opxm1*Z30_
      Z30_= - 1 + Z30_
      Z30_=opxm1*Z30_
      Z30_=2 + 15*Z30_
      Z30_=Z30_*Z59_
      Z15_= - Z15_*Z27_
      Z15_= - 1 + Z15_
      Z15_=ln2*Z15_
      Z33_=45 + Z40_
      Z15_=144*Z15_ + Z30_ + 3*Z33_ + Z17_
      Z17_= - opxm1*Z55_
      Z17_= - 15 + Z17_
      Z17_=Z17_*Z41_
      Z17_=13 + Z17_
      Z17_=Z17_*Z27_
      Z17_=349 + Z17_
      Z17_=opxm1*Z17_
      Z30_=Z44_*Z43_
      Z30_=Z30_ - 14
      Z25_=Z30_*Z25_
      Z25_=Z25_ + 7
      Z25_=Z25_*opxm1
      Z25_=Z25_ - Z49_ - 2
      Z30_= - Z25_*Z37_
      Z33_=13 - 36*omxm1
      Z33_=omxm1*Z33_
      Z17_=Z30_ + 12*Z57_ + Z17_ - 103 + Z33_
      Z30_= - 281 - 66*Z32_
      Z30_=Z30_*Z24_
      Z30_=211 + Z30_
      Z30_=Z30_*Z24_
      Z30_= - 153 + Z30_
      Z30_=opxm1*Z30_
      Z33_= - 121 + 56*omxm1
      Z33_=omxm1*Z33_
      Z30_=Z30_ + 120 + Z33_
      Z30_=Z30_*Z42_
      Z17_=8*Z17_ + Z30_
      Z17_=Z2_*Z17_
      Z25_=Z5_*Z25_
      Z15_=96*Z25_ + 4*Z15_ + Z17_
      Z17_= - 763 - 190*Z32_
      Z17_=opxm1*Z17_
      Z17_=3745 + 18*Z17_
      Z17_=Z17_*Z24_
      Z17_= - 1731 + Z17_
      Z17_=opxm1*Z17_
      Z25_= - 1219 + 488*omxm1
      Z25_=omxm1*Z25_
      Z17_=Z17_ + 1448 + Z25_
      Z17_=pi2*Z17_
      Z15_=10*Z15_ + Z17_
      Z15_=Z15_*Z36_
      Z17_= - Z58_ - Z20_
      Z17_=Z19_*Z17_
      Z20_=355 + 318*Z31_
      Z20_=opxm1*Z20_
      Z20_= - 37 + Z20_
      Z20_=Z20_*Z24_
      Z17_=Z20_ + 19 + 12*Z18_ + Z17_
      Z18_= - 517 - Z56_
      Z18_=Z18_*Z24_
      Z18_=295 + Z18_
      Z18_=Z18_*Z24_
      Z18_= - 29 + Z18_
      Z18_=opxm1*Z18_
      Z20_=19 + Z26_
      Z20_=omxm1*Z20_
      Z18_=Z18_ - 8 + Z20_
      Z18_=Z2_*Z18_
      Z17_=8*Z17_ + Z18_
      Z17_=Z2_*Z17_
      Z17_= - 64*Z35_ - 4*Z39_ + Z17_
      Z18_=121 + Z52_
      Z18_=Z18_*Z27_
      Z18_= - 193 + Z18_
      Z18_=Z18_*Z24_
      Z18_=83 + Z18_
      Z18_=opxm1*Z18_
      Z20_=35 - Z26_
      Z20_=omxm1*Z20_
      Z18_=Z18_ - 56 + Z20_
      Z18_=pi2*Z18_
      Z20_=Z3_*Z19_
      Z17_=48*Z20_ + 3*Z17_ + Z18_
      Z17_=Z17_*Z50_
      Z18_= - 59 - Z38_
      Z18_=Z18_*Z46_
      Z18_=111 + Z18_
      Z18_=opxm1*Z18_
      Z18_= - 9 + Z18_
      Z18_=opxm1*Z18_
      Z20_= - 1 + 3*omxm1
      Z20_=omxm1*Z20_
      Z18_=Z18_ + 7 + Z20_
      Z18_=Z2_*Z18_
      Z20_=313 + 276*Z31_
      Z20_=opxm1*Z20_
      Z20_= - 37 + Z20_
      Z20_=opxm1*Z20_
      Z18_=2*Z18_ + Z20_ + Z54_
      Z18_=Z7_*Z18_
      Z20_=Z10_ + Z11_
      Z20_=69120*Z20_
      Z20_=Z2_*Z20_
      Z20_=138240*Z14_ + 69120*Z13_ + Z20_
      Z19_=Z19_*Z20_
      Z15_=17280*Z18_ + Z17_ + 5*Z16_ + Z15_ + Z19_ + 8640*Z22_

*      F1AF2=loglamom*Z21_ + loglamom**2*Z23_ + A*Z28_ + A*loglamom*Z29_
*     &  + Z15_
*** STATS: original	207P 1684M 547A : 2677
*** STATS: optimized 0P 413M 341A : 754

      f1cur_2L_log2 = Z23_ * cover
      f1cur_2L_log1 = Z21_ * cover * 2.d0 ! the factors 2 and 4 are to translate log(lambda/m) <-> log(lambda^2/m^2)
      f1cur_2L_log0 = Z15_ * cover * 4.d0

      f2cur_2L_log1 = Z29_ * cover * 2.d0
      f2cur_2L_log0 = Z28_ * cover * 4.d0

      return
      end
***************************
      function Pi02L_Boncianietal(t,m2)
! hep-ph/0405275  -- eq. 93
      implicit double precision (a-h,m,o-z)
      double precision Pi02L_Boncianietal,F
      double precision t,m2
      double precision ame,amu,convfac,alpha,picommon
      common/parameters/ame,amu,convfac,alpha,picommon
      common/twolcmnPi02L/zeta2,zeta3,pi,ao2pi,ao2pi2,pi2,pi4,ic,
     .     ifirst
      data ifirst /1/
      double precision dalo(3),danlo(3),dannlo(3)
      common/dalonlonnlo/dalo,danlo,dannlo

      common/FsPi02Ldebug/F1,F2,F3,F4,F5,F6,F7,F8
      
      integer n1,n2,nw,i1,i2
      parameter(n1=-1)   !(-1,1) (0,1) (-1,0)
      parameter(n2=1)
      parameter(nw=3)
      double complex Hc1,Hc2,Hc3,Hc4
      double precision     Hr1,Hr2,Hr3,Hr4
      double precision     Hi1,Hi2,Hi3,Hi4
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     .          Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     .          Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     .     Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      if (ifirst.eq.1) then
         pi     = 4.d0*atan(1.d0)
         pi2    = pi * pi
         pi4    = pi2*pi2
         zeta2  = pi2/6.d0
         ao2pi  = alpha/pi*0.5d0
         ao2pi2 = ao2pi*ao2pi
         zeta3  = 1.202056903159594285399738161511449990764986292d0
         ifirst = 0
      endif
      
      x = (sqrt(-t+4.d0*m2)-sqrt(-t))/(sqrt(-t+4.d0*m2) + sqrt(-t))

      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4,
     .     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)

      umx   = 1.d0 - x
      umxm1 = 1.d0/umx
      umxm2 = umxm1 * umxm1
      umxm3 = umxm2 * umxm1
      umxm4 = umxm3 * umxm1

c      umxm1 = 1
c      umxm2 = 1
c      umxm3 = 1
c      umxm4 = 1

      
c      F = 0.d0      
c      F1 = -5.d0/24.d0 + 13.d0/6.d0*umxm2 - 13.d0/6.d0*umxm1
c      F2 = F + (
c     .     1.d0 - 4.d0*umxm4 + 8.d0*umxm3 - 4.d0*umxm2
c     .     ) * zeta3
c      F3 = F + (
c     .     0.25d0  + 3.d0*umxm3 - 4.5d0*umxm2 + umxm1
c     .     )*Hr1(0)
c      F4 = F + (
c     .     -4.d0/3.d0 - 16.d0/3.d0*umxm3 + 8.d0*umxm2
c     .     ) * Hr2(-1,0)
c      F5 = F + (
c     .     2.d0 - 7.d0/3.d0*umxm4 + 26.d0/3.d0*umxm3 - 23.d0/3.d0*umxm2
c     .     -2.d0/3.d0*umxm1
c     .     )
c      F6 = F + (
c     .     2.d0/3.d0 + 8.d0/3.d0*umxm3 - 4.d0*umxm2
c     .     ) * Hr2(1,0)
c      F7 = F + (
c     .     4.d0/3.d0 - 16.d0/3.d0*umxm4 + 32.d0/3.d0*umxm3
c     .     -16.d0/3.d0 * umxm2
c     .     ) * ( Hr3(0,-1,0) - Hr3(-1,0,0) )
c      F8 = F + (
c     .     2.d0/3.d0 - 8.d0/3.d0*umxm4 + 16.d0/3.d0*umxm3
c     .     -8.d0/3.d0 * umxm2
c     .     ) * ( Hr3(1,0,0) - Hr3(0,1,0))

      F  = 0.d0
      F1 = -5 + 52 * umxm2 - 52*umxm1
      F2 = F + (
     .     24 - 96*umxm4 + 192*umxm3 - 96*umxm2
     .     ) * zeta3
      F3 = F + (
     .     6  + 72*umxm3 - 108*umxm2 + 24*umxm1
     .     ) * Hr1(0)
      F4 = F + (
     .     -32 - 128*umxm3 + 192*umxm2
     .     ) * Hr2(-1,0)
      F5 = F + (
     .     48 - 56*umxm4 + 208*umxm3 - 184*umxm2 -16*umxm1
     .     )
      F6 = F + (
     .     16 + 64*umxm3 - 96*umxm2
     .     ) * Hr2(1,0)
      F7 = F + (
     .     32 - 128*umxm4 + 256*umxm3 -128 * umxm2
     .     ) * ( Hr3(0,-1,0) - Hr3(-1,0,0) )
      F8 = F + (
     .     16 - 64*umxm4 + 128*umxm3 -64 * umxm2
     .     ) * ( Hr3(1,0,0) - Hr3(0,1,0))
      F1 = F1/24
      F2 = F2/24
      F3 = F3/24
      F4 = F4/24
      F5 = F5/24
      F6 = F6/24
      F7 = F7/24
      F8 = F8/24      
      F  = F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8
      
c      print*,umxm1,umxm2,umxm3,umxm4     
c      z = danlo(2) / (ao2pi2*4.d0 )
c      print*,F1,F2,F3,F4,F5,F6,F7,F8
c      stop
      
***   expansion -- eq. 115
cc      F = -5.d0/24.d0 + zeta3 -0.25d0* log(m2/(-t))
      
!! 1-Loop -- eq. 83 -- this is exactly the one I have (and Teubner et al.'s)
! this works for e and mu bubble     
c      F1L = 5.d0/9.d0 - 4.d0/3.d0*umxm2 + 4.d0/3.d0*umxm1
c      F1L = F1L - (
c     .     1.d0/3.d0 + 4.d0/3.d0*umxm3 - 2.d0*umxm2
c     .     ) * Hr1(0)
c      F1L = -F1L / alpha *pi
!!!      
      Pi02L_Boncianietal = F*ao2pi2*4.d0 


c      F1 = F1 *ao2pi2*4.d0 
c      F2 = F2 *ao2pi2*4.d0 
c      F3 = F3 *ao2pi2*4.d0 
c      F4 = F4 *ao2pi2*4.d0 
c      F5 = F5 *ao2pi2*4.d0 
c      F6 = F6 *ao2pi2*4.d0 
c      F7 = F7 *ao2pi2*4.d0 
c      F8 = F8 *ao2pi2*4.d0       
      return
      end
**********
***** boxes with bubble NF=1 for *** bhabha (mmu must be = me!!)***
      double complex function mlog(x)
      double precision x,pi
      parameter (pi = 4.d0*atan(1.d0))
      if (x.gt.0.d0) then
         mlog = log(x)
      else
         mlog = log(abs(x)) + (0.d0,1.d0)*pi
      endif
      return
      end
***!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
****** for dalpha_lep_NLO with dispersive integral
      double precision function Rnlobit(s,m2)
      implicit double precision (a-h,m,o-z)
      common/parameters/ame,amu,convfac,alpha,pi
      common/Rnlobitcommon/ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         do k = 1,100
c            print*,'only LO in Rnlobit!'
         enddo
         ifirst = 1
      endif
      Rnlobit = 0.d0      
      if (s.le.4.d0*m2) return
      v = sqrt(1.d0-4.d0*m2/s)
      Rnlobit = v*(3.d0-v*v)*0.5d0 ! this is LO! works perfectly
      
c      Rnlobit = 3.d0/4.d0*alpha/pi * ! this works, except that I need the exact mass dependency...
c     .     (
c     .     1.d0 + 12.d0*m2/s 
c     .    + m2*m2/s/s * (10.d0 - 24.d0*log(m2/s))
c     .     )
!!! https://cds.cern.ch/record/598118/files/0212294.pdf
!!!   https://arxiv.org/pdf/hep-ph/0212294.pdf - eq. 18
!!! WORKS GREAT FOR BOTH MU AND E!!!!
      p = (1.d0-v)/(1.d0+v)
      Rnlobit = (1.d0+v*v)* (
     .     dli2(p*p) + 2.d0*dli2(p)
     .     - log(p)*log((1.d0+v)**3/8.d0/v/v)
     .     )
     . + 3.d0*v*log((1.d0-v*v)/4.d0/v) - v*log(v)
     . + (33.d0+22.d0*v*v-7.d0*v**4)/8.d0/(3.d0-v*v)*log(1.d0/p)
     . + (15.d0*v-9.d0*v**3)/4.d0/(3.d0-v*v)
      Rnlobit = Rnlobit*(3.d0-v*v)*0.5d0
      Rnlobit = Rnlobit*alpha/pi
      return
      end
***!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* for NF bubbles with dispersive integral
**************************
***************************************************************************      
*     2L vv and vrsoft with form factors (according to Mastrolia-Remiddi
*     NPB664 (2003) and hep-ph/0406046 and Burgers 1985) 
      FUNCTION TWOLOOP_ONLYNF1(S,T,U)
      implicit none
      double precision twoloop_onlyNF1,z,xx,carlolo
      double precision emubeam,EEBEAM,alpha,convfac,pi,ame,amu,egmin
      integer ihard,k
      double precision tmin,tmax,dt
      
      double complex F1NF1, F2NF1, F1NF1ex, F2NF1ex
      double complex f1lh,f2lh,f1lhq,f2lhq

      common/carlolocmnNF1/carlolo
      
      common/parameters/ame,amu,convfac,alpha,pi      

cc      COMMON/CONST/EMUBEAM,EEBEAM,ALPHA,CONVFAC,PI,AME,AMU
      COMMON/CUTOFF/EGMIN,IHARD
*
      double precision s,t,u
*
      double precision mm2,me2,aa,bb,cc,dd,burg
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
      double precision deltauv,mu2dim,lambda2,m2min
      double precision rsoft,rsoft2
      double precision rflux
      double precision kallen
      external kallen
      double precision x,omx,opx,omx2,omx3,omx4,omx5,opx2,opx3,opx4,opx5
      double precision    zeta2,zeta3
      double complex f1cur_1L_epsm1,f1cur_1L_eps0,f1cur_1L_eps1
      double complex f2cur_1L_eps0


      double complex f1e_1L_epsm1,f1e_1L_eps0,f1e_1L_eps1
      double complex f2e_1L_eps0
      double complex f1mu_1L_epsm1,f1mu_1L_eps0,f1mu_1L_eps1
      double complex f2mu_1L_eps0
      double complex F1e,F1mu,F2e,F2mu
      double complex F1ed,F1mud,F2ed,F2mud
      
      double complex f1cur_2L_log2,f1cur_2L_log1,f1cur_2L_log0
      double complex f2cur_2L_log1,f2cur_2L_log0
      double complex f1cur_2L_log2_expanded_sh,f1cur_2L_log1_expanded_sh
     +           ,f1cur_2L_log0_expanded_sh
      double complex f2cur_2L_log1_expanded_sh,f2cur_2L_log0_expanded_sh
      double complex f1cur_2L_log2_b,f1cur_2L_log1_b,f1cur_2L_log0_b
      double complex f1cur_2L_log0_pairs_b
      double complex f2cur_2L_log1_b,f2cur_2L_log0_b
      double complex F1,F1d,F2,F2d
      double complex F1oal,F1doal,F2oal,F2doal
      double complex F1oal2,F1doal2,F2oal2,F2doal2
      double complex f1f1d_2L,f1f2d_2L
      double complex matrixel,matrixelvv,matrixelrv
      double complex matrixelvve,matrixelvvmu
      double complex arg
      double complex soft
      double complex cone
      parameter(cone=(1.d0,0.d0))
      double precision log2
      double precision f1timesf1d,f2timesf2d
      double precision virt
      integer itmu,ite
      common/subsetNF1/itmu,ite
      common/gammamassNF1/lambda2

      integer ifirst,ic
      double precision ao2pi,ao2pi2,pi2,pi4
      common/twolcmn1NF1/zeta2,zeta3,mm2,me2,ao2pi,ao2pi2,pi2,pi4,ic,
     .     ifirst
      data ifirst,ic /0,0/

      double complex virtvupvdwn,F12L,F22L,F12Ld,F22Ld
      double precision carlo2L
      common/testvirtvevmu/virtvupvdwn
      
      double precision phmass
      common/photonmasslambda/phmass
      double precision tmincut,tmaxcut
      common/tmintmax/tmincut,tmaxcut
**** charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
      integer ionlybubbles,ibubbleonv,ibubbleonboxsub,iexnnlobubbles
      common/ionlybubblescmn/iexnnlobubbles,ionlybubbles,ibubbleonv,
     .     ibubbleonboxsub      
      integer bubblemask(5)
      double precision colch2(5)
      common/colorandch2factors/colch2,bubblemask
      
*     
#include "looptools.h"
*
      if (ifirst.eq.0) then
         lambda2= phmass*phmass

         pi2 = pi*pi
         pi4 = pi2*pi2
         
         zeta2=pi2/6.d0
         zeta3=1.202056903159594285399738161511449990764986292d0
*     
         mm2= amu*amu
         me2= ame*ame

         ao2pi  = alpha*0.5d0/pi
         ao2pi2 = ao2pi*ao2pi
         
         itmu = 0
         ite  = 1
         
         ifirst = 1

         if (sum(bubblemask(3:5)).gt.0) then
            print*,'twoloop_onlynf1'
            print*,'should not reach here. Stopping.'
            stop
         endif
         if (sum(bubblemask(1:2)).eq.0) then
            print*,'twoloop_onlynf1'
            print*,'useless reaching here. Stopping.'
            stop
         endif

         if (QRe.eq.1.and.QRmu.eq.1) then
            print*,'twoloop_onlynf1'
            print*,'Should reach here only for QRe Qrmu = 1 0 or 1 0'
            print*,'Stopping'
            stop
         endif         
      endif
      ic = ic + 1
      matrixelvve  = (0.d0,0.d0)
      matrixelvvmu = (0.d0,0.d0)

** E LINE ********************************************************************************
      IF (QRe.eq.1) THEN         
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
         x = (sqrt(-t+4.d0*me2)-sqrt(-t))/(sqrt(-t+4.d0*me2) + sqrt(-t))
*************************************************************************
      f1cur_2L_log2 = 0.d0
      f1cur_2L_log1 = 0.d0
      f1cur_2L_log0 = 0.d0
      f2cur_2L_log1 = 0.d0
      f2cur_2L_log0 = 0.d0

      if (bubblemask(1).eq.1) then
         call twoloop_NF1(x,me2,t,F1NF1,F2NF1)
cc ONLY PAIRS !      
         f1cur_2L_log0 = F1NF1
         f2cur_2L_log0 = F2NF1
      endif
#ifdef SMH
      if (bubblemask(2).eq.1) then
c         call twoloop_light_heavy_Syed(t,ame,amu,0,f1lh,f2lh)      
         call interpolate_F1F2_LH_HL(t,tmincut,tmaxcut,ame,amu,aa,bb
     .        ,cc,dd)
         f1lh = aa
         f2lh = bb
         f1cur_2L_log0 = f1cur_2L_log0 + f1lh
         f2cur_2L_log0 = f2cur_2L_log0 + f2lh
      endif
#endif      
* contributions to two loops from pure two loops (with lambda regularization)
*
      f1f1d_2L= 2.d0*f1cur_2L_log2*(log(lambda2/me2))**2
     +     +2.d0*f1cur_2L_log1*log(lambda2/me2)
     +     +2.d0*f1cur_2L_log0
      f1f1d_2L= ao2pi2*f1f1d_2L
      F12L = f1f1d_2L

      F22L = f2cur_2L_log1 * log(lambda2/me2)
     +     + f2cur_2L_log0
      F22L = 2.d0*ao2pi2 * F22L
*
*     this is the two-loop virtual contribution
      carlo2L = 0.d0
     &  + F22L * ( 8*t**2 + 16*mm2*t )
      carlo2L = carlo2L + F12L * ( 8*t**2 + 16*s*t + 16*s**2 - 32*me2*s
     &     + 16*me2**2 - 32*mm2*s + 32*mm2*me2 + 16*mm2**2 )
      carlo2L = carlo2L * 64.d0*pi4*ao2pi2/t/t !overall factor for born

!! this is ok      carloLO = ( 8*t**2 + 16*s*t + 16*s**2 - 32*me2*s
c     &     + 16*me2**2 - 32*mm2*s + 32*mm2*me2 + 16*mm2**2 )
c      carloLO = carloLO * 64.d0*pi4*ao2pi2/t/t !overall factor for born

      matrixelvve = carlo2L
      ENDIF
********************************************************************************
** MU LINE *********************************************************************
      IF (QRmu.eq.1) THEN
* x argument of harmonic polylog, eq. (14) of hep-ph/0406046
      x= (sqrt(-t+4.d0*mm2) - sqrt(-t))/(sqrt(-t+4.d0*mm2) + sqrt(-t))

      f1cur_2L_log2 = 0.d0
      f1cur_2L_log1 = 0.d0
      f1cur_2L_log0 = 0.d0
      f2cur_2L_log1 = 0.d0
      f2cur_2L_log0 = 0.d0

      if (bubblemask(2).eq.1) then
         call twoloop_NF1(x,mm2,t,F1NF1,F2NF1)
         f1cur_2L_log0 = F1NF1
         f2cur_2L_log0 = F2NF1
      endif
#ifdef SMH
      if (bubblemask(1).eq.1) then
c         call twoloop_light_heavy_Syed(t,amu,ame,0,f1lh,f2lh)      
         call interpolate_F1F2_LH_HL(t,tmincut,tmaxcut,ame,amu,aa,bb,
     .        cc,dd)
         f1lh = cc
         f2lh = dd
         f1cur_2L_log0 = f1cur_2L_log0 + f1lh
         f2cur_2L_log0 = f2cur_2L_log0 + f2lh
      endif
#endif
*
* contributions to two loops from pure two loops (with lambda regularization)
*
      f1f1d_2L= 2.d0*f1cur_2L_log2*(log(lambda2/mm2))**2
     +     +2.d0*f1cur_2L_log1*log(lambda2/mm2)
     +     +2.d0*f1cur_2L_log0
      f1f1d_2L= ao2pi2*f1f1d_2L

      F12L = f1f1d_2L
      F22L = f2cur_2L_log1 * log(lambda2/mm2)
     +     + f2cur_2L_log0
      F22L = 2.d0*ao2pi2 * F22L
*
*     this is the two-loop virtual contribution
      carlo2L = 0.d0
     &  + F22L * ( 8*t**2 + 16*me2*t )
      carlo2L = carlo2L + F12L * ( 8*t**2 + 16*s*t + 16*s**2 - 32*mm2*s
     &     + 16*mm2**2 - 32*me2*s + 32*me2*mm2 + 16*me2**2 )
      carlo2L = carlo2L * 64.d0*pi4*ao2pi2/t/t !overall factor for born

      matrixelvvmu = carlo2L      
      ENDIF
***********************************************************************************
* this line below is for two loop
      twoloop_onlyNF1 = dble(matrixelvve) + dble(matrixelvvmu)
      return
      END
*************************************************************************************
      subroutine leptonicbubbles(p1,p2,p3,p4,ibubblev,ibubbleb,
     .     bubbleonv,bubbleonboxsub)
#ifdef COLLIER
      use collier
#endif      
      implicit double precision (a-h,l,m,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision kernf,xi(1),kernfcv,kernfb,leptonicbubble
      real*4 xir4(1)
      double precision mph2cv,mph2b,mph2fact(5)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/tautopmasses/taumass,topmass
      integer iarun
      common/ialpharunning/iarun
      double precision p1c(0:3),p2c(0:3),p3c(0:3),p4c(0:3)
      double precision probs(5)
      common/leptbubbleDR/ovall,me,me2,mm,mm2,p1c,p2c,p3c,p4c,mph2fact,
     .     probs,ifirst
      data ifirst /0/
**** charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
      integer ioneloop
      common/oneloopcalc/ioneloop
      common/forbubbles/ctandv0,box0
      
      common/tmintmax/tmin,tmax
********************************************************
***   Carlo -- common block for energies
      parameter(imaxknt=1586)
      dimension wKNT(imaxknt),wk(0:imaxknt)
      common/energiesRhadKNT/wKNT
******
      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
**** from invariants.h      
      double precision vpdal
      double complex vpdah
      double precision vpvect(5)
      common/vp_dalphas/vpdal,vpdah,vpvect
***
      if (ifirst.eq.0) then
         xx = wrapRhad(1.d0) ! just to fill the common enegiesRhadKNT
         
         me = ame
         mm = ammu
         me2 = me*me
         mm2 = mm*mm
**********************
         mph2fact(1) = me2
         mph2fact(2) = mm2
         mph2fact(3) = taumass*taumass
         mph2fact(4) = topmass*topmass
         mph2fact(5) = wKNT(1)*wKNT(1)*0.25d0

         print*,'Including',sum(bubblemask),' blob(s)'
         print*,'blob ''mask''',bubblemask
         print*,'blob masses',sqrt(mph2fact)
**********************
         ovall = alpha/3.d0/pi

         probs(1) = 0.999d0    * bubblemask(1)
         probs(2) = 0.001d0    * bubblemask(2)
         probs(3) = 0.00001d0  * bubblemask(3)
         probs(4) = 0.000001d0 * bubblemask(4)
         probs(5) = 0.001d0    * bubblemask(5)

c         probs(1) = 0.999d0    * bubblemask(1)
c         probs(2) = 0.001d0    * bubblemask(2)
c         probs(3) = 0.00001d0  * bubblemask(3)
c         probs(4) = 0.000001d0 * bubblemask(4)
c         probs(5) = 0.1d0    * bubblemask(5)
       
         probs = probs/sum(probs)

         ifirst = 1
      endif
      
      bubbleonv      = 0.d0
      bubbleonboxsub = 0.d0
      if (ibubblev.eq.0.and.ibubbleb.eq.0) return
      if (iarun.eq.0) return

      if (ioneloop.eq.0) then
         mph2cv = 0.d0
         mph2b  = 0.d0
         box0   = 0.d0
         if (ibubbleb.eq.1) then
#ifdef COLLIER
            call InitEvent_cll(3)
#else
            call clearcache
#endif
            call getoneloopforbubble(p1,p2,p3,p4,0,1,mph2b,
     .           born,ctandvuu,box0)
            box0 = box0 * 0.5d0 ! because for 0 ph mass box0 is twice the actual box that I need
         endif
      endif

#ifdef COLLIER
      call InitEvent_cll(4)
#else
      call clearcache
#endif
      
      wcv = 1.d0
      wb  = 1.d0
      
      p1c = p1
      p2c = p2
      p3c = p3
      p4c = p4
      
      t = dot(p1-p3,p1-p3)

      xmin = 0.5d0
      xmax = 1.d0

!      doing a MC sum over flavours
      wfl = 1.d0
      iwhichblob = 1
      call anotherrng(xi,1)
      if (xi(1).lt.probs(1)) then
         iwhichblob = 1
         wfl        = 1.d0/probs(1)
      elseif(xi(1).lt.(probs(1)+probs(2))) then
         iwhichblob = 2
         wfl        = 1.d0/probs(2)
      elseif(xi(1).lt.(probs(1)+probs(2)+probs(3))) then
         iwhichblob = 3
         wfl        = 1.d0/probs(3)
      elseif(xi(1).lt.(probs(1)+probs(2)+probs(3)+probs(4))) then
         iwhichblob = 4
         wfl        = 1.d0/probs(4)
      else
         iwhichblob = 5
         wfl        = 1.d0/probs(5)
      endif
      
      kernfcv = 0.d0
      kernfb  = 0.d0
      if (ibubblev.gt.0) then
         if (iwhichblob.le.4) then
c      an = log(xmax/xmin)
c      call getrnd(xi,1)
c      xcv = xmin*exp(xi(1)*an)
c      wcv = wcv * an * xcv
** flat sampling
c      an = xmax - xmin
c      call getrnd(xi,1)
c      x = an * xi(1) + xmin
c      w = w * an
** optimized sampling, as 1/(1-x)^1/2, for vertex bubble...
            a   = 0.5d0
            uma = 1.d0-a
            an  = (1.d0-xmin)**uma - (1.d0-xmax)**uma
            an  = an/uma
!     I use another prng in order to not break the sequence and lose correlations, just in case
            call anotherrng(xi,1)            
            xcv = 1.d0 - ((1.d0-xmin)**uma-uma*an*xi(1))**(1.d0/uma)            
            wcv = an *(1.d0-xcv)**a

            ycv     = xcv / (1.d0-xcv)
            mph2cv  = mph2fact(iwhichblob) * 4.d0*ycv
            kernfcv = sqrt(abs(2.d0*xcv-1.d0)/xcv)*
     .           0.5d0*(1.d0 + xcv)/xcv  * colch2(iwhichblob)
         else
!! this for Rhad!!
cc            call getxsampleddiscrete(xmin,xmax,xcv,wcv,mph2cv)
            call getxsampleddiscretev2(1,xmin,xmax,xcv,wcv,mph2cv)
            kernfcv = wrapRhad(mph2cv) * colch2(iwhichblob)
         endif
         kernfcv = kernfcv/xcv/(1.d0-xcv)
         kernfcv = kernfcv*ovall       !* 0.d0
      endif

      if (ibubbleb.gt.0) then
         if (iwhichblob.le.4) then
            th = t*0.25d0/mph2fact(iwhichblob)
            Q  = (th + 1.d0) / th ! this is always < 1...
** optimized sampling, as 1/(1-Q*x), for boxes, where Q = (th+1) / th and th = t/4/mp2
            an = 1.d0/Q*log((1.d0-Q*xmin)/(1.d0-Q*xmax))
            call anotherrng(xi,1)

            xb = 1.d0 - (1.d0-xmin*Q)*exp(-an*Q*xi(1))
            xb = xb/Q
            wb = an * (1.d0-Q*xb)
         
            yb    = xb / (1.d0-xb)
            mph2b = mph2fact(iwhichblob)*4.d0*yb

            kernfb = sqrt(abs(2.d0*xb-1.d0)/xb)*
     .           0.5d0*(1.d0 + xb)/xb * colch2(iwhichblob)
         else            
!! this for Rhad!!
cc            call getxsampleddiscrete(xmin,xmax,xb,wb,mph2b)
            call getxsampleddiscretev2(2,xmin,xmax,xb,wb,mph2b)
            kernfb = wrapRhad(mph2b) * colch2(iwhichblob)
         endif
         kernfb = kernfb/xb/(1.d0-xb)
         kernfb = kernfb*ovall      !* 0.d0
      endif

      ctandv = 0.d0
      box    = 0.d0
      call getoneloopforbubble(p1,p2,p3,p4,ibubblev,0,
     .     mph2cv,born,ctandv,boxuu)
      call getoneloopforbubble(p1,p2,p3,p4,0,ibubbleb,
     .     mph2b,born,ctandvuu,box)

**********
      bubbleonv      = kernfcv* ctandv * wcv
      bubbleonboxsub = kernfb * box    * wb
*********

      if (iarun.ne.0) then
         vpsubbox = t/(t-mph2b) * box0 * kernfb*wb ! kernfb has already colch2
         bubbleonboxsub = bubbleonboxsub - vpsubbox
!     if using the integrated vp subtraction. For fixed order calculation, vpsubint is
!     added back in the calling routine
c         vpsubint       = vpvect(iwhichblob)
c         bubbleonboxsub = bubbleonboxsub - vpsubint*box0
      endif


c      if (iwhichblob.eq.5) wfl = 0.d0
      
      bubbleonv      = bubbleonv      * wfl
      bubbleonboxsub = bubbleonboxsub * wfl
      
c      if ((.not.(abs(bubbleonv).lt.1d50).or.isnan(bubbleonv)).or.
c     .  (.not.(abs(bubbleonboxsub).lt.1d50).or.isnan(bubbleonboxsub)))
c     .     then
c         print*,'print in twoloopvirtual'
c         print*,p1
c         print*,p2
c         print*,p3
c         print*,p4
c         print*,iwhichblob,colch2(iwhichblob),mph2fact(iwhichblob)
c         print*,'vert',t,bubbleonv,xcv,ycv,wcv,kernfcv,mph2cv,ctandv,
c     .        xiv,anv,wfl,probs,probs(iwhichblob)
c         print*,'box',t,bubbleonboxsub,xb,yb,wb,kernfb,mph2b,vpboxsub,
c     .        box,xib,anb,Q,th,wfl,probs,probs(iwhichblob)
c         print*,'end print in twoloopvirtual'
c      endif
     
      return
      end
**********************************
      function leptonicbubblequadpack_fun(p1,p2,p3,p4,ictandv,iboxsub)
#ifdef COLLIER
      use collier
#endif      
      implicit double precision (a-h,l,m,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision kernf,xi(1),mph2b
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision p1c(0:3),p2c(0:3),p3c(0:3),p4c(0:3)
      common/leptbubblequad/ovall,me,me2,mm,mm2,p1c,p2c,p3c,p4c,
     .     t,sffbox,born,mp2,ifirst
      data ifirst /0/

      double precision ABSERR,EPSABS,EPSREL,RESULT,WORK
      INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL
      DIMENSION IWORK(100),WORK(400)

****  charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu      
      integer ioneloop
      common/oneloopcalc/ioneloop
***      
      double precision wrapintegrand
      EXTERNAL wrapintegrand

      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask

      common/whichtocalculatewrapper/ictv,ibsub

      common/forbubbles/ctandv0,box0
      
      if (ifirst.eq.0) then
         me = ame
         mm = ammu
         me2 = me*me
         mm2 = mm*mm
         ovall = alpha/3.d0/pi
         ifirst = 1
      endif

c      t      = dot(p1-p3,p1-p3)      
c      sff    = soft_integral_mue(sqrt(-t),p1,p2,p3,p4)      
c      sffbox = sffij(1,2)+sffij(3,4)+sffij(1,4)+sffij(2,3)
      if (ioneloop.eq.0) then
         mph2b  = 0.d0
         box0   = 0.d0
         if (iboxsub.eq.1) then
#ifdef COLLIER
         call InitEvent_cll(3)
#else
         call clearcache
#endif         
            call getoneloopforbubble(p1c,p2c,p3c,p4c,0,1,mph2b,
     .           born,ctandvuu,box0)
            box0 = box0 * 0.5d0 ! because for 0 ph mass box0 is twice the actual box that I need
         endif
      endif      
#ifdef COLLIER
      call InitEvent_cll(4)
#else
      call clearcache
#endif         

      print*,'mh... not to be used, because it sums box and vertexes...'
      print*,'stopping'
      stop
      
      ictv  = ictandv
      ibsub = iboxsub
      
      p1c = p1
      p2c = p2
      p3c = p3
      p4c = p4      

      xmin = 0.5d0
      xmax = 1.d0

      epsabs = 0d0
      epsrel = 1d-5
      key    = 4
      limit  = 100
      lenw   = 4 * limit
      
c!!!      call dqag(wrapintegrand,xmin,xmax,EPSABS,EPSREL,KEY,RESULT,
c!!!     .     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

      call dqags(wrapintegrand,xmin,xmax,EPSABS,EPSREL,RESULT,
     .     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
cc      call dqags(wrapintegrandwithRhad,xmin,xmax,EPSABS,EPSREL,RESULT,
cc    .     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

ccc      print*,dot(p1c-p3c,p1c-p3c),neval,result,abserr
      
      leptonicbubblequadpack_fun = result

      return
      end
**********************************
      subroutine leptonicbubblequadpack(p1,p2,p3,p4,ictandv,iboxsub,
     .     ctandv,boxsub)
#ifdef COLLIER
      use collier
#endif      
      implicit double precision (a-h,l,m,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision kernf,xi(1),mph2b
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision p1c(0:3),p2c(0:3),p3c(0:3),p4c(0:3)
      common/leptbubblequad/ovall,me,me2,mm,mm2,p1c,p2c,p3c,p4c,
     .     t,sffbox,born,mp2,ifirst
      data ifirst /0/

      double precision ABSERR,EPSABS,EPSREL,RESULT,WORK
      INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL
      DIMENSION IWORK(100),WORK(400)

****  charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu      
      integer ioneloop
      common/oneloopcalc/ioneloop
***      
      double precision wrapintegrandctandv
      EXTERNAL wrapintegrandctandv
      double precision wrapintegrandboxsub
      EXTERNAL wrapintegrandboxsub

      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask

      common/whichtocalculatewrapper/ictv,ibsub

      common/forbubbles/ctandv0,box0
      
      if (ifirst.eq.0) then
         me = ame
         mm = ammu
         me2 = me*me
         mm2 = mm*mm
         ovall = alpha/3.d0/pi
         ifirst = 1
      endif

c      t      = dot(p1-p3,p1-p3)      
c      sff    = soft_integral_mue(sqrt(-t),p1,p2,p3,p4)      
c      sffbox = sffij(1,2)+sffij(3,4)+sffij(1,4)+sffij(2,3)
      if (ioneloop.eq.0) then
         mph2b  = 0.d0
         box0   = 0.d0
         if (iboxsub.eq.1) then
#ifdef COLLIER
            call InitEvent_cll(3)
#else
            call clearcache
#endif         
            call getoneloopforbubble(p1c,p2c,p3c,p4c,0,1,mph2b,
     .           born,ctandvuu,box0)
            box0 = box0 * 0.5d0 ! because for 0 ph mass box0 is twice the actual box that I need
         endif
      endif

      ctandv = 0.d0
      boxsub = 0.d0
      if (ictandv.eq.0.and.iboxsub.eq.0) return
      
#ifdef COLLIER
      call InitEvent_cll(3)
#else
      call clearcache
#endif         

      ictv  = ictandv
      ibsub = iboxsub
      
      p1c = p1
      p2c = p2
      p3c = p3
      p4c = p4

      xmin = 0.5d0
      xmax = 1.d0

      epsabs = 0d0
      epsrel = 1d-5
      key    = 4
      limit  = 100
      lenw   = 4 * limit
      
c!!!      call dqag(wrapintegrand,xmin,xmax,EPSABS,EPSREL,KEY,RESULT,
c!!!     .     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

      if (ictv.eq.1) then
         call dqags(wrapintegrandctandv,xmin,xmax,EPSABS,EPSREL,RESULT,
     .        ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         ctandv = result
      endif
      if (ibsub.eq.1) then
         call dqags(wrapintegrandboxsub,xmin,xmax,EPSABS,EPSREL,RESULT,
     .        ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
         boxsub = result
      endif

      return
      end
**********************************
      function dalnloquadpack(t,m2)
      implicit double precision (a-h,l,m,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi

      common/forRnlobit/m2c,alphac,q2c,pic
      
      double precision ABSERR,EPSABS,EPSREL,RESULT,WORK
      INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL
      DIMENSION IWORK(100),WORK(400)

      double precision wrapdispersivedanlo
      EXTERNAL wrapdispersivedanlo

      m2c    = m2
      alphac = alpha
      q2c    = t
      pic    = pi

      xmin = 0.5d0
      xmax = 1.d0

      epsabs = 0d0
      epsrel = 1d-10
      key    = 4
      limit  = 100
      lenw   = 4 * limit
      
      call dqags(wrapdispersivedanlo,xmin,xmax,EPSABS,EPSREL,RESULT,
     .     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

      dalnloquadpack = result

      return
      end
*********************************
      function wrapdispersivedanlo(x)
      implicit double precision (a-h,m,o-z)
      common/forRnlobit/m2,alpha,q2,pi
      s = 4.d0*m2*x/(1.d0-x)
      wrapdispersivedanlo =
     .     1.d0/(1.d0-x)/x*Rnlobit(s,m2)*1.d0/(q2-s)
      wrapdispersivedanlo = wrapdispersivedanlo * q2*alpha/3.d0/pi
      return
      end
*********************************      
**********************************
      function dahadquadpack(t)
      implicit double precision (a-h,l,m,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi

      common/forRhadquadpack/m2c,alphac,q2c,pic
      
      double precision ABSERR,EPSABS,EPSREL,RESULT,WORK
      INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL
      DIMENSION IWORK(100),WORK(400)

      double precision wrapdispersivedahad
      EXTERNAL wrapdispersivedahad

      common/dahadquadpackcmn/ifirst
      data ifirst /0/
***   Carlo -- common block for energies
      parameter(imaxknt=1586)
      dimension wKNT(imaxknt),wk(0:imaxknt)
      common/energiesRhadKNT/wKNT
******

      if (ifirst.eq.0) then
         xx = wrapRhad(1.d0)    ! just to fill the common enegiesRhadKNT

         ifirst = 1
      endif
         
      
      m2c    = wKNT(1)*wKNT(1)*0.25d0
      alphac = alpha
      q2c = t
      pic = pi

      xmin = 0.5d0
      xmax = 1.d0

      epsabs = 0d0
      epsrel = 1d-5
      key    = 4
      limit  = 100
      lenw   = 4 * limit
      
      call dqags(wrapdispersivedahad,xmin,xmax,EPSABS,EPSREL,RESULT,
     .     ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)

      dahadquadpack = result

      return
      end
*********************************

**********************************
      function dahadquadpackv2(t)
      implicit double precision (a-h,l,m,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi

      common/forRhadquadpack/m2c,alphac,q2c,pic
      
      double precision ABSERR,EPSABS,EPSREL,RESULT,WORK
      INTEGER IER,IWORK,KEY,LAST,LENW,LIMIT,NEVAL,LENIW

      parameter (nbreak=2)

      dimension break(nbreak)
      DIMENSION IWORK(100+nbreak),WORK(400+nbreak)

      double precision wrapdispersivedahad
      EXTERNAL wrapdispersivedahad

      common/dahadquadpackcmn/ifirst
      data ifirst /0/
***   Carlo -- common block for energies
      parameter(imaxknt=1586)
      dimension wKNT(imaxknt),wk(0:imaxknt)
      common/energiesRhadKNT/wKNT
******

      parameter (np = 3885)     ! number of points in Rhad-scan.dat
      double precision rhens(np),rhad(np,2)
      common/rhadscans/rhens,rhad

      dimension x(0:np)
      
      if (ifirst.eq.0) then
         xx = wrapRhad(1.d0)    ! just to fill the common enegiesRhadKNT
         ifirst = 1
      endif
      
      m2c    = wKNT(1)*wKNT(1)*0.25d0
      alphac = alpha
      q2c = t
      pic = pi

      xmin = 0.5d0
      xmax = 1.d0

      epsabs = 0d0
      epsrel = 1d-12
      key    = 4
      limit  = 100
      lenw   = 4 * limit + nbreak
      leniw  = 2 * limit + nbreak

      break = 0.d0
      
      sred = rhens(1180)**2/rhens(1)**2
c      break(1) = 0.5d0!sred/(1.d0+sred)

      sred = rhens(2892)**2/rhens(1)**2
c      break(2) = 1.d0!sred/(1.d0+sred)

      
c      call dqagp(wrapdispersivedahad,xmin,xmax,nbreak,break,EPSABS,
c     .     EPSREL,RESULT,ABSERR,NEVAL,IER,leniw,LENW,LAST,IWORK,WORK)

      do k = 0,np-1
         sred = rhens(k+1)**2/rhens(1)**2
         x(k) = sred/(1.d0+sred)
      enddo
      x(np) = 1.d0         
      result = 0.d0
      do k = 0,np-2
         dx = x(k+1)-x(k)
         fx = wrapdispersivedahad(x(k))
         dfx = (wrapdispersivedahad(x(k+1))-fx)/dx
         result = result + dx * (fx + 0.5d0*dfx*dx)
      enddo
      
      dahadquadpackv2 = result

      return
      end
*********************************
      function wrapdispersivedahad(x)
      implicit double precision (a-h,m,o-z)
      common/forRhadquadpack/m2,alpha,q2,pi
      s = 4.d0*m2*x/(1.d0-x)
      rh = wrapRhad(s)
      wrapdispersivedahad =
     .     1.d0/(1.d0-x)/x*rh*1.d0/(q2-s)
      wrapdispersivedahad = wrapdispersivedahad * q2*alpha/3.d0/pi
      return
      end
*********************************      
      function wrapintegrand(x) ! this sums boxes and vertexes
      implicit double precision (a-h,m,o-z)
      double precision p1c(0:3),p2c(0:3),p3c(0:3),p4c(0:3)
      double precision kernf
      double precision mph2,mph2b
      common/leptbubblequad/ovall,me,me2,mm,mm2,p1c,p2c,p3c,p4c,
     .     t,sffbox,born,mp2,ifirstaa
      data ifirst /0/      
****  charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu      
***      
      common/whichtocalculatewrapper/ictv,ibsub

      double precision mph2fact(5)
      integer ifirst
      common/wrapmasses/mph2fact,ifirst
c      data ifirst /0/

      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
      
      common/parameters/ame,ammu,convfac,alpha,pi
      common/tautopmasses/taumass,topmass
      common/forbubbles/ctandv0,box0
      integer ioneloop
      common/oneloopcalc/ioneloop

***   Carlo -- common block for energies
      parameter(imaxknt=1586)
      dimension wKNT(imaxknt),wk(0:imaxknt)
      common/energiesRhadKNT/wKNT
******
      print*,'mh... not to be used, because it sums box and vertexes...'
      print*,'stopping'
      stop
      
      if (ifirst.eq.0) then
         xx = wrapRhad(1.d0) ! just to fill the common enegiesRhadKNT
         me = ame
         mm = ammu
         me2 = me*me
         mm2 = mm*mm
         mph2fact(1) = me2
         mph2fact(2) = mm2
         mph2fact(3) = taumass*taumass
         mph2fact(4) = topmass*topmass
         mph2fact(5) = wKNT(1)*wKNT(1)*0.25d0
         ifirst = 1
      endif

      wrapintegrand = 0.d0
      
      do iwhichblob = 1,5
         if (bubblemask(iwhichblob).gt.0) then
            
            y      = x / (1.d0-x)
            mph2   = 4.d0*mph2fact(iwhichblob)*y

            if(iwhichblob.eq.5) then
               kernf = wrapRhad(mph2)
            else
               kernf  = sqrt((2.d0*x-1.d0)/x)*0.5d0*(1.d0 + x)/x ! this is Rlep
            endif
            kernf  = kernf/x/(1.d0-x)*ovall*colch2(iwhichblob)
            
            ctandv = 0.d0
            box    = 0.d0
            ivert  = ictv
            ibox   = ibsub
            t      = dot(p1c-p3c,p1c-p3c)
            call getoneloopforbubble(p1c,p2c,p3c,p4c,ivert,ibox,
     .           mph2,born,ctandv,box)
            
            if (ibox.eq.1) then
               box = box - t/(t-mph2)*box0
            endif
            
            ampm2 = ivert*ctandv + ibox*box

            wrapintegrand = wrapintegrand + kernf*ampm2

         endif
      enddo
      
      return
      end
***************************
      function wrapintegrandctandv(x) ! this is only vertexes and ct
      implicit double precision (a-h,m,o-z)
      double precision p1c(0:3),p2c(0:3),p3c(0:3),p4c(0:3)
      double precision kernf
      double precision mph2,mph2b
      common/leptbubblequad/ovall,me,me2,mm,mm2,p1c,p2c,p3c,p4c,
     .     t,sffbox,born,mp2,ifirstaa
      data ifirst /0/      
****  charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu      
***      
      common/whichtocalculatewrapper/ictv,ibsub

      double precision mph2fact(5)
      integer ifirst
      common/wrapmasses1/mph2fact,ifirst
      data ifirst /0/

      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
      
      common/parameters/ame,ammu,convfac,alpha,pi
      common/tautopmasses/taumass,topmass
      common/forbubbles/ctandv0,box0
      integer ioneloop
      common/oneloopcalc/ioneloop

***   Carlo -- common block for energies
      parameter(imaxknt=1586)
      dimension wKNT(imaxknt),wk(0:imaxknt)
      common/energiesRhadKNT/wKNT
******
      if (ifirst.eq.0) then
         xx = wrapRhad(1.d0) ! just to fill the common enegiesRhadKNT
         me = ame
         mm = ammu
         me2 = me*me
         mm2 = mm*mm
         mph2fact(1) = me2
         mph2fact(2) = mm2
         mph2fact(3) = taumass*taumass
         mph2fact(4) = topmass*topmass
         mph2fact(5) = wKNT(1)*wKNT(1)*0.25d0
         ifirst = 1
      endif

      wrapintegrandctandv = 0.d0      
      do iwhichblob = 1,5
         if (bubblemask(iwhichblob).gt.0) then
            
            y      = x / (1.d0-x)
            mph2   = 4.d0*mph2fact(iwhichblob)*y

            if(iwhichblob.eq.5) then
               kernf = wrapRhad(mph2)
            else
               kernf  = sqrt((2.d0*x-1.d0)/x)*0.5d0*(1.d0 + x)/x ! this is Rlep
            endif
            kernf  = kernf/x/(1.d0-x)*ovall*colch2(iwhichblob)
            
            ctandv = 0.d0
            box    = 0.d0
            ivert  = ictv
            ibox   = 0
            t      = dot(p1c-p3c,p1c-p3c)
            call getoneloopforbubble(p1c,p2c,p3c,p4c,ivert,ibox,
     .           mph2,born,ctandv,box)

            wrapintegrandctandv = wrapintegrandctandv + kernf*ctandv
         endif
      enddo
      
      return
      end
***************************
*********************************      
      function wrapintegrandboxsub(x) ! this has only subtracted boxes
      implicit double precision (a-h,m,o-z)
      double precision p1c(0:3),p2c(0:3),p3c(0:3),p4c(0:3)
      double precision kernf
      double precision mph2,mph2b
      common/leptbubblequad/ovall,me,me2,mm,mm2,p1c,p2c,p3c,p4c,
     .     t,sffbox,born,mp2,ifirstaa
      data ifirst /0/      
****  charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu      
***      
      common/whichtocalculatewrapper/ictv,ibsub

      double precision mph2fact(5)
      integer ifirst
      common/wrapmasses2/mph2fact,ifirst
      data ifirst /0/

      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
      
      common/parameters/ame,ammu,convfac,alpha,pi
      common/tautopmasses/taumass,topmass
      common/forbubbles/ctandv0,box0
      integer ioneloop
      common/oneloopcalc/ioneloop

***   Carlo -- common block for energies
      parameter(imaxknt=1586)
      dimension wKNT(imaxknt),wk(0:imaxknt)
      common/energiesRhadKNT/wKNT
******
      if (ifirst.eq.0) then
         xx = wrapRhad(1.d0) ! just to fill the common enegiesRhadKNT
         me = ame
         mm = ammu
         me2 = me*me
         mm2 = mm*mm
         mph2fact(1) = me2
         mph2fact(2) = mm2
         mph2fact(3) = taumass*taumass
         mph2fact(4) = topmass*topmass
         mph2fact(5) = wKNT(1)*wKNT(1)*0.25d0
         ifirst = 1
      endif

      wrapintegrandboxsub = 0.d0
      
      do iwhichblob = 1,5
         if (bubblemask(iwhichblob).gt.0) then
            y      = x / (1.d0-x)
            mph2   = 4.d0*mph2fact(iwhichblob)*y

            if(iwhichblob.eq.5) then
               kernf = wrapRhad(mph2)
            else
               kernf  = sqrt((2.d0*x-1.d0)/x)*0.5d0*(1.d0 + x)/x ! this is Rlep
            endif
            kernf  = kernf/x/(1.d0-x)*ovall*colch2(iwhichblob)
            
            ctandv = 0.d0
            box    = 0.d0
            ivert  = 0
            ibox   = ibsub
            t      = dot(p1c-p3c,p1c-p3c)
            call getoneloopforbubble(p1c,p2c,p3c,p4c,ivert,ibox,
     .           mph2,born,ctandv,box)
            
            if (ibox.eq.1) then
               box = box - t/(t-mph2)*box0
            endif

            wrapintegrandboxsub = wrapintegrandboxsub + kernf*ibox*box
         endif
      enddo
      return
      end
***************************      
      subroutine getoneloopforbubble(p1,p2,p3,p4,icv,ib,mph2,born,
     .     ctandv,box)
#ifdef COLLIER
      use collier
#endif      
      implicit double precision (a-h,m,o-z)
**** wraps the one loop with a given "photon mass"
      double precision mph2
#include "invariants.h"
      integer hh,icv,ib,n
      double complex amplborn
      double complex boxamp,vertamp

      born = 0.d0
      ctandv = 0.d0
      box = 0.d0
      if (icv.eq.0.and.ib.eq.0) return
      
      call loadinvariants(p1,p2,p3,p4)
      cts = 0.d0
      
      if (icv.eq.1) call getQEDcts_withphmass(cts,mph2)
      if (icv.eq.1) call loadQEDvertexfuns_withphmass(mph2)
      if (ib.eq.1 ) call loadQEDboxesfuns_withphmass(mph2)
      
      t   = dot(p1-p3,p1-p3)
      tm1 = 1.d0/t

      ampm2    = 0.d0
      born     = 0.d0
      oneloop  = 0.d0
      vertamp  = czero
      boxamp   = czero
      ctandv   = 0.d0
      box      = 0.d0
      vert     = 0.d0
      
      do hh = 0, 7
         h1  = mod(hh,16)/8 + 1
         h2  = mod(hh, 8)/4 + 1
         h3  = mod(hh, 4)/2 + 1
         h4  = mod(hh, 2)   + 1

         if (Qmu.lt.0.d0) then
            call getcurrentnew(Qmu,h3,h1,p3,p1,mm,mm,cmu,cImu,c5mu)
         else
            call getcurrentnew(Qmu,h3,h1,p1,p3,mm,mm,cmu,cImu,c5mu)
         endif            
         call getcurrentnew(-1.d0,h4,h2,p4,p2,me,me,ce,cIe,c5e)

         cmuce    = dotcmom(cmu,ce)
         amplborn = cmuce*im*ec2*tm1 
                  
         bbar    = amplborn*conjg(amplborn)
         born    = born + bbar

         call loadspinvariants(p1,p2,p3,p4,cmu,ce)
         if (icv.eq.1) call getQEDvertexes_withphmass(vertamp)
         if (ib.eq.1)  call getQEDboxes_withphmass(boxamp)

         ctandv = ctandv + bbar * (1.d0*cts )
     .        + 1.d0 * 2.d0*dreal(amplborn*conjg(vertamp))
         box    = box  + 2.d0*dreal(amplborn*conjg(boxamp))
         vert   = vert + 2.d0*dreal(amplborn*conjg(vertamp))
      enddo
      
      born   = 2.d0*born
      ctandv = 2.d0*ctandv
      box    = 2.d0*box
      vert   = 2.d0*vert


c      if (ib.gt.0) print*,dot(p1+p2,p1+p2),t,mph2,box
      
      return
      end      
*********************************************************
      subroutine getQEDvertexes_withphmass(vertexes)
      implicit none
      double complex C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu
      double complex C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte
      double complex C00ltmu,C00lte,vertexes,cover,amplvertexes
      double complex amplvertexesmu,amplvertexese,vertexmu,vertexe
      integer ifirst,n
      common/ifirstvertwphmass/cover,ifirst
      data ifirst /0/
#include "invariants.h"
      common/vertexfunswphmass/
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte      

      double complex amplve,amplvmu
      common/vertexeslegswphmass/amplve,amplvmu

      double precision colch2(5),an
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask      
c      common/colorandch2factors/colch2
      
      if (ifirst.eq.0) then
         cover = ec4*im*pi216m1*4.d0
         ifirst = 1
      endif

      an = 1.d0
      
      if (Qmu.lt.0.d0) then
      amplvmu =
     .  + QRmu * (  - an*1.D0/2.D0*cmuce + cmuce*C00ltmu + cmuce*p1p3*
     .    C12ltmu + 2*cmuce*p1p3*C1ltmu + cmuce*p1p3*C0ltmu + mm2*cmuce
     .    *C11ltmu + 2*mm2*cmuce*C1ltmu - cImu*mm*p3ce*C12ltmu - cImu*
     .    mm*p3ce*C11ltmu - cImu*mm*p3ce*C1ltmu - cImu*mm*p1ce*C12ltmu
     .     - cImu*mm*p1ce*C11ltmu - cImu*mm*p1ce*C1ltmu )
      amplve =
     . + QRe * (  - an*1.D0/2.D0*cmuce + cmuce*
     .    C00lte + cmuce*p2p4*C12lte + 2*cmuce*p2p4*C1lte + cmuce*p2p4*
     .    C0lte + me2*cmuce*C11lte + 2*me2*cmuce*C1lte - cIe*me*p4cmu*
     .    C12lte - cIe*me*p4cmu*C11lte - cIe*me*p4cmu*C1lte - cIe*me*
     .    p2cmu*C12lte - cIe*me*p2cmu*C11lte - cIe*me*p2cmu*C1lte )         
      else
         amplvmu =
     .  + QRmu * (  - an*1.D0/2.D0*cmuce + cmuce*C00ltmu + cmuce*p1p3*
     .    C12ltmu + 2*cmuce*p1p3*C1ltmu + cmuce*p1p3*C0ltmu + mm2*cmuce
     .    *C11ltmu + 2*mm2*cmuce*C1ltmu + cImu*mm*p3ce*C12ltmu + cImu*
     .    mm*p3ce*C11ltmu + cImu*mm*p3ce*C1ltmu + cImu*mm*p1ce*C12ltmu
     .     + cImu*mm*p1ce*C11ltmu + cImu*mm*p1ce*C1ltmu )
         amplve =
     .    + QRe * (  - an*1.D0/2.D0*cmuce + cmuce*
     .    C00lte + cmuce*p2p4*C12lte + 2*cmuce*p2p4*C1lte + cmuce*p2p4*
     .    C0lte + me2*cmuce*C11lte + 2*me2*cmuce*C1lte - cIe*me*p4cmu*
     .    C12lte - cIe*me*p4cmu*C11lte - cIe*me*p4cmu*C1lte - cIe*me*
     .    p2cmu*C12lte - cIe*me*p2cmu*C11lte - cIe*me*p2cmu*C1lte )
      endif
      amplve  = amplve  * tm1 * cover
      amplvmu = amplvmu * tm1 * cover

      vertexes = amplve+amplvmu      
      return
      end
**************************************************************************
**************************************************************************
      subroutine getQEDboxes_withphmass(boxes)
      implicit none
      double complex Z1_,Z2_,Z3_,Z4_,Z5_,Z6_,Z7_
      double complex Z8_,Z9_,Z10_,Z11_,Z12_,Z13_
      double complex Z14_,Z15_,Z16_,Z17_,Z18_,Z19_,Z20_
      double complex Z21_,Z22_,Z23_,Z24_,Z25_,Z26_
      double complex Z27_
      double complex D0ltd,D1ltd,D2ltd,D3ltd
      double complex D0ltc,D1ltc,D2ltc,D3ltc
      double complex D00ltd,D11ltd,D12ltd,D13ltd
      double complex D22ltd,D23ltd,D33ltd
      double complex D00ltc,D11ltc,D12ltc,D13ltc
      double complex D22ltc,D23ltc,D33ltc
      double complex boxes,cover
      integer ifirst
      data ifirst /0/
      common/ifirstboxwphmass/cover,ifirst
#include "invariants.h"
      common/boxesfunswphmass/
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,     
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc

      if (QRmu.eq.0.or.QRe.eq.0) then
         boxes = (0.d0,0.d0)
         return
      endif
      
      if (ifirst.eq.0) then
         cover  = ec4*im*pi216m1*2.d0
         ifirst = 1
      endif

      if (Qmu.lt.0.d0) then
      Z1_= - D3ltd - D33ltd
      Z1_=cmuce*Z1_
      Z2_=D3ltc*cmuce
      Z3_= - D33ltd*c5muc5e
      Z4_=cmuce - c5muc5e
      Z5_=D33ltc*Z4_
      Z1_=Z2_ + Z5_ + Z3_ + Z1_
      Z1_=me2*Z1_
      Z3_=p3ce*p4cmu
      Z5_=p1ce*p4cmu
      Z3_=Z3_ - Z5_
      Z6_=p3c5mu - p1c5mu
      Z7_= - p4c5e*Z6_
      Z8_=p1c5e - p3c5e
      Z9_=p4c5mu*Z8_
      Z10_=p3p4*Z4_
      Z7_= - 2*Z10_ + Z9_ + Z7_ - Z3_
      Z7_=D23ltc*Z7_
      Z9_=p3ce*p2cmu
      Z10_=p1ce*p2cmu
      Z9_=Z9_ - Z10_
      Z11_=p1p2*cmuce
      Z12_=Z11_ + e12cec5mu
      Z13_=mm2 - p1p3
      Z13_=Z13_*cmuce
      Z14_= - p2p3*cmuce
      Z14_=e23cec5mu + Z14_ - e13cmuc5e - Z13_ - Z9_ + Z12_
      Z14_=D2ltd*Z14_
      Z15_=mm*cImu
      Z16_=p1ce - p3ce
      Z17_= - Z16_*Z15_
      Z18_=p1p4 - p3p4
      Z18_=cmuce*Z18_
      Z3_=Z17_ - e13cmuc5e - e34cec5mu + Z13_ + Z18_ - Z3_
      Z3_=D2ltc*Z3_
      Z17_=p3c5mu*p1c5e
      Z18_=p3c5e*p1c5mu
      Z19_=p1c5e*p1c5mu
      Z17_= - 2*Z19_ + Z17_ + Z18_
      Z20_=p1p3*c5muc5e
      Z21_=mm2*c5muc5e
      Z20_=Z20_ - Z21_
      Z22_=Z20_ - Z13_
      Z22_=2*Z22_
      Z23_=Z22_ - Z17_
      Z23_=D12ltd*Z23_
      Z13_=Z20_ + Z13_
      Z13_=2*Z13_
      Z17_=Z13_ - Z17_
      Z17_=D12ltc*Z17_
      Z20_=Z8_*p3c5mu
      Z18_= - Z19_ + Z20_ + Z18_
      Z13_=Z13_ - Z18_
      Z13_=D22ltc*Z13_
      Z20_=cmuce*mm2
      Z12_= - Z20_ + Z10_ + Z12_
      Z12_=D1ltd*Z12_
      Z24_=p4c5e*p1c5mu
      Z25_=p4c5mu*p1c5e
      Z24_=Z5_ + Z24_ + Z25_
      Z24_=D13ltc*Z24_
      Z25_=p1p2 - p2p3
      Z26_=cmuce + c5muc5e
      Z25_=Z26_*Z25_
      Z6_=p2c5e*Z6_
      Z6_=Z6_ + 2*Z25_ - Z9_
      Z6_=D23ltd*Z6_
      Z9_=Z22_ - Z18_
      Z9_=D22ltd*Z9_
      Z18_=p1p2*Z26_
      Z22_= - p2c5e*p1c5mu
      Z10_=Z22_ + 2*Z18_ + Z10_
      Z10_=D13ltd*Z10_
      Z18_=p2c5e*D33ltd
      Z8_= - D23ltd*Z8_
      Z22_= - D13ltd*p1c5e
      Z8_=Z22_ + Z18_ + Z8_
      Z8_=p2c5mu*Z8_
      Z18_=D23ltc + D13ltc
      Z18_=2*Z18_
      Z4_=Z4_*Z18_
      Z18_=D1ltc + 2*D0ltc
      Z18_=cmuce*Z18_
      Z2_=Z2_ + Z4_ + Z18_
      Z2_=p1p4*Z2_
      Z4_=D12ltc - D12ltd
      Z18_= - p3ce*Z4_
      Z4_= - D11ltd + Z4_
      Z4_=p1ce*Z4_
      Z16_=D2ltd*Z16_
      Z4_=Z16_ + Z4_ + Z18_
      Z4_=Z4_*Z15_
      Z16_=Z21_ - Z19_
      Z18_=p1ce*Z15_
      Z18_=Z18_ + Z20_ - Z16_
      Z18_=D11ltc*Z18_
      Z19_= - D1ltc + D13ltd + D13ltc - D1ltd - D3ltd - D3ltc
      Z15_=Z19_*Z15_
      Z19_= - D33ltd*p2cmu
      Z21_=p4cmu*D33ltc
      Z15_=Z15_ + Z19_ + Z21_
      Z15_=cIe*me*Z15_
      Z19_=2*D0ltd + D3ltd
      Z11_=Z11_*Z19_
      Z16_= - Z20_ - Z16_
      Z16_=D11ltd*Z16_
      Z19_=D3ltd*p2cmu
      Z21_=D3ltc*p4cmu
      Z19_=Z19_ + Z21_
      Z19_=p1ce*Z19_
      Z5_=Z20_ + Z5_
      Z5_=D1ltc*Z5_
      Z20_=3*c5muc5e
      Z21_=5*cmuce
      Z22_= - Z20_ - Z21_
      Z22_=D00ltd*Z22_
      Z20_= - Z20_ + Z21_
      Z20_=D00ltc*Z20_
      Z21_=p4c5mu*D33ltc*p4c5e
      Z25_=D3ltc*e14cmuc5e
      Z26_= - e12cmuc5e*D3ltd
      Z27_=D1ltc + D2ltc
      Z27_=e14cec5mu*Z27_
      boxes=Z27_ + Z15_ + Z3_ + Z20_ + Z26_ + Z18_ + Z4_ + Z2_ + Z22_
     &  + Z14_ + Z7_ + Z1_ + Z5_ + Z8_ + Z10_ + Z9_ + Z6_ + Z24_ + Z12_
     &  + Z13_ + Z19_ + Z25_ + Z16_ + Z17_ + Z23_ + Z21_ + Z11_        
      else
      Z1_=D11ltd - D11ltc
      Z2_=D2ltd + D1ltd
      Z3_=D2ltc + D1ltc
      Z4_= - Z3_ + Z2_ + Z1_
      Z4_=mm2*Z4_
      Z5_=D13ltd + D23ltd
      Z6_=Z5_ + D0ltd
      Z2_= - D3ltd - Z2_ - 2*Z6_
      Z2_=p1p2*Z2_
      Z6_=D23ltc + D13ltc
      Z7_=2*Z6_
      Z8_= - 2*D0ltc - D3ltc - Z7_ - Z3_
      Z8_=p1p4*Z8_
      Z9_=p2p3*D23ltd
      Z10_=p3p4*D23ltc
      Z11_=Z9_ + Z10_
      Z12_=D2ltc*p3p4
      Z13_=D2ltd*p2p3
      Z14_=D3ltd - D3ltc - D33ltc + D33ltd
      Z14_=me2*Z14_
      Z15_=D2ltc - D2ltd
      Z15_=p1p3*Z15_
      Z2_=Z4_ + 5*D00ltd + Z15_ + Z14_ + Z13_ + 2*Z11_ + Z12_ + Z8_ + 
     & Z2_
      Z2_=cmuce*Z2_
      Z4_=D3ltd + D1ltd
      Z8_= - Z5_ - Z4_
      Z8_=p1ce*Z8_
      Z11_=p3ce - p1ce
      Z12_=Z11_*D2ltd
      Z13_=p3ce*D23ltd
      Z14_=cIe*me
      Z15_=D33ltd*Z14_
      Z8_=Z15_ + Z12_ + Z13_ + Z8_
      Z8_=p2cmu*Z8_
      Z13_=D3ltc + D1ltc
      Z15_= - Z6_ - Z13_
      Z15_=p1ce*Z15_
      Z16_=Z11_*D2ltc
      Z17_=p3ce*D23ltc
      Z18_= - D33ltc*Z14_
      Z15_=Z18_ + Z16_ + Z15_ + Z17_
      Z15_=p4cmu*Z15_
      Z17_=p1p2*Z5_
      Z9_=Z17_ - Z9_ + Z10_
      Z10_= - D33ltc - D33ltd
      Z10_=me2*Z10_
      Z7_= - p1p4*Z7_
      Z17_=D11ltd + D11ltc
      Z18_= - mm2*Z17_
      Z7_=Z18_ - 3*D00ltd + Z7_ + Z10_ + 2*Z9_
      Z7_=c5muc5e*Z7_
      Z4_=D13ltd + D13ltc - Z13_ - Z4_
      Z4_=Z4_*Z14_
      Z1_= - p1ce*Z1_
      Z1_=Z4_ - Z12_ + Z1_ + Z16_
      Z1_=cImu*Z1_
      Z4_=Z11_*cImu
      Z9_=D12ltd*Z4_
      Z1_=Z1_ + Z9_
      Z1_=mm*Z1_
      Z9_= - p1c5e*Z5_
      Z10_=p3c5e*D23ltd
      Z11_=p2c5e*D33ltd
      Z9_=Z11_ + Z9_ + Z10_
      Z9_=p2c5mu*Z9_
      Z10_= - p3c5mu*D23ltc
      Z11_=p1c5mu*Z6_
      Z12_=p4c5mu*D33ltc
      Z10_=Z12_ + Z10_ + Z11_
      Z10_=p4c5e*Z10_
      Z11_=2*c5muc5e
      Z12_=mm2 - p1p3
      Z11_=Z11_*Z12_
      Z13_=p1c5e*p3c5mu
      Z11_=Z11_ + Z13_
      Z13_=p1c5mu - p3c5mu
      Z13_=Z13_*p3c5e
      Z14_=p1c5mu*p1c5e
      Z13_=Z11_ + Z13_ - Z14_
      Z16_=2*cmuce
      Z12_=Z16_*Z12_
      Z16_= - Z12_ - Z13_
      Z16_=D22ltc*Z16_
      Z18_=p3c5e*p1c5mu
      Z11_= - Z11_ - Z18_ + 2*Z14_
      Z18_=Z12_ + Z11_
      Z18_=D12ltd*Z18_
      Z13_=Z12_ - Z13_
      Z13_=D22ltd*Z13_
      Z4_= - mm*Z4_
      Z4_=Z4_ - Z12_ + Z11_
      Z4_=D12ltc*Z4_
      Z11_=Z17_*Z14_
      Z12_=p3c5mu*D23ltd
      Z5_= - p1c5mu*Z5_
      Z5_=Z12_ + Z5_
      Z5_=p2c5e*Z5_
      Z6_=p1c5e*Z6_
      Z12_= - p3c5e*D23ltc
      Z6_=Z6_ + Z12_
      Z6_=p4c5mu*Z6_
      Z12_=e12cec5mu*D1ltd
      Z14_=e12cmuc5e*D3ltd
      Z3_=e14cec5mu*Z3_
      Z17_=e23cec5mu + e12cec5mu
      Z17_=D2ltd*Z17_
      Z19_=D2ltc + D2ltd
      Z19_=e13cmuc5e*Z19_
      Z20_= - 3*c5muc5e - 5*cmuce
      Z20_=D00ltc*Z20_
      Z21_= - e34cec5mu*D2ltc
      Z22_= - e14cmuc5e*D3ltc
      boxes=Z4_ + Z22_ + Z1_ + Z13_ + Z21_ + Z20_ + Z18_ + Z16_ + Z19_
     &  + Z10_ + Z2_ + Z7_ + Z9_ + Z15_ + Z8_ + Z17_ + Z6_ + Z3_ + Z5_
     &  + Z14_ + Z12_ + Z11_         
      endif
      boxes = boxes * cover
      return
      end
***************************************************************************
      subroutine getxsampleddiscrete(xmin,xmax,x,wx,s) ! this is for Rhad, works very well
!     from Hto4l
#ifdef COLLIER
      use collier
#endif
      implicit double precision (a-h,m,o-z)

      double precision pin1(0:3),pin2(0:3),p3(0:3),p4(0:3)
      double precision p1(0:3),p2(0:3)
      common/momentainitial/pin1,pin2
      common/tmintmax/tmin,tmax
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
c      common/colorandch2factors/colch2

***   Carlo -- common block for energies
      parameter(imaxknt=1586)
      dimension wKNT(imaxknt),wk(0:imaxknt)
      common/energiesRhadKNT/wKNT
*****
!     parameter (np = 5000)
c      parameter (np = 20000)
      parameter (np = imaxknt)
      dimension xk(0:np),fk(0:np),ck(0:np),tgk(0:np)
      dimension ckcumul(0:np),xi(1)
      double precision mph2
      real*4 xi4(1)
      integer ifirst,imidguess(1)
      common/xsamplingdiscrete/dx,totint,xk,fk,ck,tgk,ckcumul,
     .     wk,imidguess,ifirst
      data ifirst /0/

      print*,'this is only for Rhad from KNT. v2 is also for Fred'
      print*,'stopping'
      stop
      
      if (ifirst.eq.0) then

c         print*,'It''s now specialized for Rhad!!'
c         print*,'It''s now specialized for Rhad!!'
c         print*,'It''s now specialized for Rhad!!'
         
         xx         = wrapRhad(1.d0) ! just to fill the common
         wk(0:1585) = wKNT(1:1586)
         wk(np)     = wk(np-1)
         
c         do k = 0,imaxknt-1
c            sred = wk(k)**2/wk(0)**2
c            xk(k) = sred/(1.d0+sred)
c            print*,xk(k)
c         enddo
c         stop
                     
c in userinterface         colch2(1:3) = 1.d0
c         colch2(4)   = 4.d0/3.d0 ! i.e. for top 3 colors * (2/3 charge)^2 = 12/9 = 4/3
         
         pm    = sqrt(tridot(pin1,pin1))
         tmean = (tmin -tmax)*0.5d0
         cth   = 1.d0+tmean*0.5d0/pm/pm
         
         p1 = pin1
         p2 = pin2
         p3 = 0d0
         p4 = 0d0
         p3(0) = sqrt(pm*pm+ammu*ammu)
         p4(0) = sqrt(pm*pm+ame*ame)
         p3(1) = pm*sqrt(1.d0-cth*cth)
         p3(3) = pm*cth
         p4(1:3) = -p3(1:3)         
         
         dx    = (xmax - xmin)/np
         xk(0) = xmin
         do k = 1,np
            xk(k) = xk(0) + k*dx
         enddo

         do k = 0,np-1
            sred = wk(k)**2/wk(0)**2
            xk(k) = sred/(1.d0+sred)
         enddo
         xk(np) = 1.d0         
         
c      kernfcv = sqrt(abs(2.d0*xcv-1.d0)/xcv)*0.5d0*(1.d0 + xcv)/xcv
c      kernfcv = kernfcv/xcv/(1.d0-xcv)

#ifdef COLLIER         
         call SwitchOffCacheSystem_cll
#else
         call clearcache
#endif 
         do k = 0,np-1
            xcv    = xk(k)
            ycv    = xcv / (1.d0-xcv)
            mph2 = 4.d0*ame**2 * ycv

            mph2 = wk(0)**2 * ycv
            rhad = wrapRhad(mph2)

            call getoneloopforbubble(p1,p2,p3,p4,1,0,mph2,tree,v,b)
            fk(k) = v * rhad
!     .           sqrt(abs(2.d0*xcv-1.d0)/xcv)*0.5d0*(1.d0 + xcv)/xcv
            fk(k) = fk(k)/xcv/(1.d0-xcv)
            fk(k) = abs(fk(k))
         enddo
         fk(np) = fk(np-1)
#ifdef COLLIER         
         call SwitchOnCacheSystem_cll
#else
         call clearcache
#endif         

         do k = 0,np-1
            tgk(k) = (fk(k+1)-fk(k))/(xk(k+1)-xk(k))
         enddo
         tgk(np) = 0.d0
         
         ck(0) = 0.d0
         do k = 1,np
            dx    = xk(k)-xk(k-1)
            ck(k) = dx * (fk(k-1)+0.5d0*tgk(k-1)*dx)
         enddo
         
         ckcumul(0) = 0.d0
         do k = 1,np
            ckcumul(k) = ckcumul(k-1) + ck(k)
         enddo
         totint = ckcumul(np)
         
         imidguess = maxloc(ck) - 1

c         do k = 0,np
c            write(99,*)k,xk(k),fk(k),tgk(k),
c     .           ck(k)/totint,ckcumul(k)/totint
c         enddo
c         stop
         ifirst = 1
      endif

      call anotherrng(xi,1)      
      csiN = xi(1)*totint

** not binary search...      
c      k = 0
c      do while(csiN.gt.ckcumul(k))
c         k = k + 1
c      enddo
c      k = k - 1
c      cumul = ckcumul(k)
**      
***   binary search
      istart  = 0
      ifinish = np - 1
      imid = imidguess(1)       !(istart + ifinish)/2
      do while(.not.(ckcumul(imid+1).gt.csiN.and.ckcumul(imid).le.csiN))
         if (csiN.gt.ckcumul(imid)) then
            istart = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
      cumul = ckcumul(k)
****      
      a = tgk(k)
      b = fk(k)
      c = 2.d0*(cumul-csiN)
      if (abs(a).gt.0.d0) then
         xt = (-b + sqrt(b*b-a*c))/a
      else
         xt = -c/b*0.5d0
      endif
      x  = xt + xk(k)
      fx = fk(k) + tgk(k)*(x-xk(k))

      wx = totint/fx

      s = wk(0)**2 * x / (1.d0-x)

      return
      end
******************************
      subroutine getxsampleddiscretev2(i,xmin,xmax,x,wx,s) ! this is for Rhad, works very well
! i = 1 vertexes + cts, i = 2 boxes
!     from Hto4l
#ifdef COLLIER
      use collier
#endif
      implicit double precision (a-h,m,o-z)
      integer i
      double precision pin1(0:3),pin2(0:3),p3(0:3),p4(0:3)
      double precision p1(0:3),p2(0:3)
      common/momentainitial/pin1,pin2
      common/tmintmax/tmin,tmax
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
c      common/colorandch2factors/colch2

***   Carlo -- common block for energies
      parameter(imaxknt=1586)
      dimension wKNT(imaxknt)
      common/energiesRhadKNT/wKNT

      parameter (np = 3885)     ! number of points in Rhad-scan.dat
      double precision rhens(np),rhad(np,2)
      common/rhadscans/rhens,rhad      
*****
      dimension xk(0:np,2),fk(0:np,2),ck(0:np,2),tgk(0:np,2)
      dimension ckcumul(0:np,2),xi(1),wk(0:np,2)
      double precision mph2
      real*4 xi4(1)
      integer ifirst,imidguess(1)
      double precision totint(2)
      common/xsamplingdiscretev2/dx,totint,xk,fk,ck,tgk,ckcumul,
     .     wk,imidguess,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         
         xx           = wrapRhad(1.d0) ! just to fill the commons
         wk(0:np-1,1) = rhens(1:np)
         wk(0:np-1,2) = rhens(1:np)
         wk(np,1)     = wk(np-1,1)
         wk(np,2)     = wk(np-1,2)
         
         pm    = sqrt(tridot(pin1,pin1))
         tmean = (tmin -tmax)*0.5d0
         cth   = 1.d0+tmean*0.5d0/pm/pm
         
         p1 = pin1
         p2 = pin2
         p3 = 0d0
         p4 = 0d0
         p3(0) = sqrt(pm*pm+ammu*ammu)
         p4(0) = sqrt(pm*pm+ame*ame)
         p3(1) = pm*sqrt(1.d0-cth*cth)
         p3(3) = pm*cth
         p4(1:3) = -p3(1:3)         
         
         dx    = (xmax - xmin)/np
         xk(0,1) = xmin
         xk(0,2) = xmin
         do k = 1,np
            xk(k,1) = xk(0,1) + k*dx
            xk(k,2) = xk(0,2) + k*dx
         enddo

         do k = 0,np-1
            sred    = wk(k,1)**2/wk(0,1)**2
            xk(k,1) = sred/(1.d0+sred)
            sred    = wk(k,2)**2/wk(0,2)**2
            xk(k,2) = sred/(1.d0+sred)
         enddo
         xk(np,1) = 1.d0
         xk(np,2) = 1.d0
         
#ifdef COLLIER         
         call SwitchOffCacheSystem_cll
#else
         call clearcache
#endif

         call getoneloopforbubble(p1,p2,p3,p4,0,1,0.d0,
     .        born,ctandvuu,box0)
         box0 = box0 * 0.5d0    ! because for 0 ph mass box0 is twice the actual box that I need

         do k = 0,np-1
            do j = 1,2
               xcv    = xk(k,j)
               ycv    = xcv / (1.d0-xcv)

               mph2 = wk(0,j)**2 * ycv
               
               rhadcur = wrapRhad(mph2)

               if (j.eq.1) then
                 call getoneloopforbubble(p1,p2,p3,p4,1,0,mph2,tree,u,b)
               else
                 call getoneloopforbubble(p1,p2,p3,p4,0,1,mph2,tree,b,u)
                 u = u - tmean / (tmean - mph2)*box0 
               endif
               fk(k,j) = u * rhadcur
               fk(k,j) = fk(k,j)/xcv/(1.d0-xcv)
               fk(k,j) = abs(fk(k,j))
            enddo
         enddo
         fk(np,1) = fk(np-1,1)
         fk(np,2) = fk(np-1,2)
         
#ifdef COLLIER         
         call SwitchOnCacheSystem_cll
#else
         call clearcache
#endif         

         do k = 0,np-1
            tgk(k,1) = (fk(k+1,1)-fk(k,1))/(xk(k+1,1)-xk(k,1))
            tgk(k,2) = (fk(k+1,2)-fk(k,2))/(xk(k+1,2)-xk(k,2))
         enddo
         tgk(np,1) = 0.d0
         tgk(np,2) = 0.d0
         
         ck(0,1) = 0.d0
         ck(0,2) = 0.d0
         do k = 1,np
            dx      = xk(k,1)-xk(k-1,1)
            ck(k,1) = dx * (fk(k-1,1)+0.5d0*tgk(k-1,1)*dx)
            dx      = xk(k,2)-xk(k-1,2)
            ck(k,2) = dx * (fk(k-1,2)+0.5d0*tgk(k-1,2)*dx)
         enddo

         ckcumul(0,1) = 0.d0
         ckcumul(0,2) = 0.d0
         do k = 1,np
            ckcumul(k,1) = ckcumul(k-1,1) + ck(k,1)
            ckcumul(k,2) = ckcumul(k-1,2) + ck(k,2)
         enddo
         totint(1) = ckcumul(np,1)
         totint(2) = ckcumul(np,2)

         imidguess = maxloc(ck(:,1)) - 1

         ifirst = 1
      endif
      
      call anotherrng(xi,1)
      csiN = xi(1)*totint(i)

** not binary search...      
c      k = 0
c      do while(csiN.gt.ckcumul(k))
c         k = k + 1
c      enddo
c      k = k - 1
c      cumul = ckcumul(k)
**      
***   binary search
      istart  = 0
      ifinish = np - 1
      
      imid    = imidguess(1)    !(istart + ifinish)/2
      do while(.not.(ckcumul(imid+1,i).gt.csiN
     .     .and.ckcumul(imid,i).le.csiN))
         if (csiN.gt.ckcumul(imid,i)) then
            istart = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k     = imid
      cumul = ckcumul(k,i)
****      
      a = tgk(k,i)
      b = fk(k,i)
      c = 2.d0*(cumul-csiN)
      if (abs(a).gt.0.d0) then
         xt = (-b + sqrt(b*b-a*c))/a
      else
         xt = -c/b*0.5d0
      endif
      x  = xt + xk(k,i)
      fx = fk(k,i) + tgk(k,i)*(x-xk(k,i))

      wx = totint(i)/fx

      s = wk(0,i)**2 * x / (1.d0-x)

      return
      end
******************************

      subroutine getxsampleddiscretelept(xmin,xmax,x,wx)
!     from Hto4l
#ifdef COLLIER      
      use collier
#endif
      implicit double precision (a-h,m,o-z)

      double precision pin1(0:3),pin2(0:3),p3(0:3),p4(0:3)
      double precision p1(0:3),p2(0:3)
      common/momentainitial/pin1,pin2
      common/tmintmax/tmin,tmax
      common/parameters/ame,ammu,convfac,alpha,pi
      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
c      common/colorandch2factors/colch2

      parameter(imaxknt=2000)
*****
      parameter (np = imaxknt)
      dimension xk(0:np),fk(0:np),ck(0:np),tgk(0:np)
      dimension ckcumul(0:np),xi(1)
      double precision mph2
      real*4 xi4(1)
      integer ifirst,imidguess(1)
      common/xsamplingdiscretelept/dx,totint,xk,fk,ck,tgk,ckcumul,
     .     imidguess,ifirst
      data ifirst /0/
      
      if (ifirst.eq.0) then

         do k = 1,10
            print*,'the ''standard'' sampling seems to be better...'
         enddo
         pm    = sqrt(tridot(pin1,pin1))
         tmean = (tmin -tmax)*0.5d0
         cth   = 1.d0+tmean*0.5d0/pm/pm
         
         p1 = pin1
         p2 = pin2
         p3 = 0d0
         p4 = 0d0
         p3(0) = sqrt(pm*pm+ammu*ammu)
         p4(0) = sqrt(pm*pm+ame*ame)
         p3(1) = pm*sqrt(1.d0-cth*cth)
         p3(3) = pm*cth
         p4(1:3) = -p3(1:3)         
         
         dx    = (xmax - xmin)/np
         xk(0) = xmin
         do k = 1,np
            xk(k) = xk(0) + k*dx
         enddo
         
#ifdef COLLIER         
         call SwitchOffCacheSystem_cll
#else
         call clearcache
#endif 
         do k = 0,np-1
            xcv    = xk(k)
            ycv    = xcv / (1.d0-xcv)
            mph2 = 4.d0*ame**2 * ycv

            call getoneloopforbubble(p1,p2,p3,p4,1,0,mph2,tree,v,b)
            fk(k) = v * 
     .           sqrt(abs(2.d0*xcv-1.d0)/xcv)*0.5d0*(1.d0 + xcv)/xcv
            fk(k) = fk(k)/xcv/(1.d0-xcv)
            fk(k) = abs(fk(k))
         enddo
         fk(np) = fk(np-1)
#ifdef COLLIER         
         call SwitchOnCacheSystem_cll
#else
         call clearcache
#endif
         do k = 0,np-1
            tgk(k) = (fk(k+1)-fk(k))/(xk(k+1)-xk(k))
         enddo
         tgk(np) = 0.d0
         
         ck(0) = 0.d0
         do k = 1,np
            dx    = xk(k)-xk(k-1)
            ck(k) = dx * (fk(k-1)+0.5d0*tgk(k-1)*dx)
         enddo
         
         ckcumul(0) = 0.d0
         do k = 1,np
            ckcumul(k) = ckcumul(k-1) + ck(k)
         enddo
         totint = ckcumul(np)
         
         imidguess = maxloc(ck) - 1

c         do k = 0,np
c            write(99,*)xk(k),fk(k),tgk(k),ck(k),ckcumul(k)
c         enddo
         
         ifirst = 1
      endif

      call anotherrng(xi,1)
      
      csiN = xi(1)*totint

** not binary search...      
c      k = 0
c      do while(csiN.gt.ckcumul(k))
c         k = k + 1
c      enddo
c      k = k - 1
c      cumul = ckcumul(k)
**      
***   binary search
      istart  = 0
      ifinish = np - 1
      imid = imidguess(1)       !(istart + ifinish)/2
      do while(.not.(ckcumul(imid+1).gt.csiN.and.ckcumul(imid).le.csiN))
         if (csiN.gt.ckcumul(imid)) then
            istart = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
      cumul = ckcumul(k)
****      
      a = tgk(k)
      b = fk(k)
      c = 2.d0*(cumul-csiN)
      if (abs(a).gt.0.d0) then
         xt = (-b + sqrt(b*b-a*c))/a
      else
         xt = -c/b*0.5d0
      endif
      x  = xt + xk(k)
      fx = fk(k) + tgk(k)*(x-xk(k))

      wx = totint/fx

      return
      end      
*****************************************************
      subroutine getQEDcts_withphmass(cts,mph2)
      implicit none
      double precision cts,ctstmp,mph2,mph2tmp
      double complex tmp
      double complex B1lt,dB0lt,dB1lt
#include "invariants.h"
      integer ifirst,n,k
      common/qedctswphmass/ctstmp,ifirst
      data ifirst /0/

      double precision ctstmpmu,ctstmpe      
      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask

      ctstmp = 0.d0
      do k = 1,1
         mph2tmp = mph2
#ifdef COLLIER
        call getBsforctscollier_withphmass(mm2,mph2tmp,B1lt,dB0lt,dB1lt)
#else
        call getBsforcts_withphmass(mm2,mph2tmp,B1lt,dB0lt,dB1lt)
#endif         
      tmp   = B1lt + 0.5d0 + 2.d0*mm2*dB1lt + 4.d0*mm2*dB0lt
      
      cts = 2.d0*tmp *pi216m1 * ec2 * QRmu

#ifdef COLLIER
      call getBsforctscollier_withphmass(me2,mph2tmp,B1lt,dB0lt,dB1lt)
#else
      call getBsforcts_withphmass(me2,mph2tmp,B1lt,dB0lt,dB1lt)         
#endif
      tmp   = B1lt + 0.5d0 + 2.d0*me2*dB1lt + 4.d0*me2*dB0lt
         
      cts = cts + 2.d0*tmp *pi216m1 * ec2 * QRe
      
      ctstmp   = ctstmp + 2.d0 * cts

      enddo
      
      cts = ctstmp*1.d0
      return
      end
******************************************************************
      subroutine getBsforcts_withphmass(m2,mph2,B1o,dB0o,dB1o)
      implicit none
      double complex B1o,dB0o,dB1o
      double precision m2,mph2
      integer*8 i
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      i    = Bget(m2*uquad,m2*uquad,mph2*uquad)
      B1o  = Bval(bb1,i)
      dB0o = Bval(dbb0,i)
      dB1o = Bval(dbb1,i)
      return
      end
***********************************************
      subroutine loadQEDvertexfuns_withphmass(mph2)
      implicit none
      double precision mph2,mph2tmp
      integer n,k
      double complex C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu
      double complex C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte
      double complex C00ltmu,C00lte      
#include "invariants.h"
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      integer*8 imu,ie
      common/vertexfunswphmass/
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte      
      data 
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte
     .     /14*(0.d0,0.d0)/

      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
c      common/colorandch2factors/colch2
      
#ifdef COLLIER
      call Cfunscollier_withphmass(mph2)
      return
#else
      C0ltmu  = czero
      C1ltmu  = czero
      C00ltmu = czero
      C11ltmu = czero
      C12ltmu = czero
      C0lte   = czero
      C1lte   = czero
      C00lte  = czero
      C11lte  = czero
      C12lte  = czero
      
      if (QRmu.gt.0) then
         do k = 1,1
            mph2tmp  = mph2
            imu = Cget(mm2*uquad,
     .           2.d0*mm2*uquad-2.d0*p1p3*uquad,
     .           mm2*uquad,mph2tmp*uquad,mm2*uquad,mm2*uquad)
            C0ltmu  = Cval(cc0,imu)
            C1ltmu  = Cval(cc1,imu)
            C00ltmu = Cval(cc00,imu)
            C11ltmu = Cval(cc11,imu)
            C12ltmu = Cval(cc12,imu)
         enddo
         
      endif

      if (QRe.gt.0) then
         do k = 1,1
            mph2tmp  = mph2
            ie  = Cget(me2*uquad,
     .           2.d0*me2*uquad-2.d0*p2p4*uquad,
     .           me2*uquad,mph2tmp*uquad,me2*uquad,me2*uquad)
            C0lte   = Cval(cc0,ie)
            C1lte   = Cval(cc1,ie)
            C00lte  = Cval(cc00,ie)
            C11lte  = Cval(cc11,ie)
            C12lte  = Cval(cc12,ie)
         enddo
      endif
#endif
      return
      end
*********************************************************
      subroutine loadQEDboxesfuns_withphmass(mph2)
      implicit none
      double precision mph2,mph2tmp
      integer n,k
      double complex D0ltd,D1ltd,D2ltd,D3ltd
      double complex D0ltc,D1ltc,D2ltc,D3ltc
      double complex D00ltd,D11ltd,D12ltd,D13ltd
      double complex D22ltd,D23ltd,D33ltd
      double complex D00ltc,D11ltc,D12ltc,D13ltc
      double complex D22ltc,D23ltc,D33ltc
      double precision pard(10),parc(10)
#include "invariants.h"
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      
      integer*8 id1,ic1,id2,ic2
      common/boxesfunswphmass/
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,     
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc
      integer ifirst
      common/loadboxeswphmass/pard,parc,ifirst
      data ifirst /0/

      data 
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,     
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc
     .     /22*(0.d0,0.d0)/

      double precision colch2(5),an
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
c      common/colorandch2factors/colch2
*******************
      if (QRmu.eq.0.or.QRe.eq.0) return
******************      
#ifdef COLLIER
      call Dfunscollier_withphmass(mph2)
      return
#else
      if (ifirst.eq.0) then
         pard(1)  = mm2
         pard(2)  = mm2
         pard(3)  = me2
         pard(4)  = me2
         pard(8)  = mm2
         pard(10) = me2

         parc = pard
         
         ifirst = 1
      endif

      pard(5) = 2.d0*mm2-2.d0*p1p3
      pard(6) = mm2+me2+2.d0*p3p4
      pard(9)  = 0.d0
      parc(5) = 2.d0*mm2-2.d0*p1p3
      parc(6) = mm2+me2-2.d0*p2p3
      parc(9)  = 0.d0

      D0ltd = czero
      D1ltd = czero
      D2ltd = czero
      D3ltd = czero
      D0ltc = czero
      D1ltc = czero
      D2ltc = czero
      D3ltc = czero
      D00ltd = czero
      D11ltd = czero
      D12ltd = czero
      D13ltd = czero
      D22ltd = czero
      D23ltd = czero
      D33ltd = czero
      D00ltc = czero
      D11ltc = czero
      D12ltc = czero
      D13ltc = czero
      D22ltc = czero
      D23ltc = czero
      D33ltc = czero
      
      do k = 1,1

         mph2tmp = mph2
         
         pard(7)  = mph2tmp
         parc(7)  = mph2tmp
               
         id1 = Dget(pard(1)*uquad,
     .        pard(2)*uquad,
     .        pard(3)*uquad,
     .        pard(4)*uquad,
     .        pard(5)*uquad,
     .        pard(6)*uquad,
     .        pard(7)*uquad,
     .        pard(8)*uquad,
     .        pard(9)*uquad,
     .        pard(10)*uquad)

!! the second box is not needed, it's enough to multiply by two the first!      
         D0ltd  = Dval(dd0,id1)  * 2.d0 ! + Dval(dd0,id2)
         D1ltd  = Dval(dd1,id1)  * 2.d0 ! + Dval(dd1,id2)
         D2ltd  = Dval(dd2,id1)  * 2.d0 ! + Dval(dd2,id2)
         D3ltd  = Dval(dd3,id1)  * 2.d0 ! + Dval(dd3,id2)
         D00ltd = Dval(dd00,id1) * 2.d0 ! + Dval(dd00,id2)
         D11ltd = Dval(dd11,id1) * 2.d0 ! + Dval(dd11,id2)
         D12ltd = Dval(dd12,id1) * 2.d0 ! + Dval(dd12,id2)
         D13ltd = Dval(dd13,id1) * 2.d0 ! + Dval(dd13,id2)
         D22ltd = Dval(dd22,id1) * 2.d0 ! + Dval(dd22,id2)
         D23ltd = Dval(dd23,id1) * 2.d0 ! + Dval(dd23,id2)
         D33ltd = Dval(dd33,id1) * 2.d0 ! + Dval(dd33,id2)
         
         ic1 = Dget(parc(1)*uquad,
     .        parc(2)*uquad,
     .        parc(3)*uquad,
     .        parc(4)*uquad,
     .        parc(5)*uquad,
     .        parc(6)*uquad,
     .        parc(7)*uquad,
     .        parc(8)*uquad,
     .        parc(9)*uquad,
     .        parc(10)*uquad)
      
!! the second box is not needed, it's enough to multiply by two the first!      
         D0ltc  = Dval(dd0,ic1)  * 2.d0 ! + Dval(dd0,ic2)
         D1ltc  = Dval(dd1,ic1)  * 2.d0 ! + Dval(dd1,ic2)
         D2ltc  = Dval(dd2,ic1)  * 2.d0 ! + Dval(dd2,ic2)
         D3ltc  = Dval(dd3,ic1)  * 2.d0 ! + Dval(dd3,ic2)
         D00ltc = Dval(dd00,ic1) * 2.d0 ! + Dval(dd00,ic2)
         D11ltc = Dval(dd11,ic1) * 2.d0 ! + Dval(dd11,ic2)
         D12ltc = Dval(dd12,ic1) * 2.d0 ! + Dval(dd12,ic2)
         D13ltc = Dval(dd13,ic1) * 2.d0 ! + Dval(dd13,ic2)
         D22ltc = Dval(dd22,ic1) * 2.d0 ! + Dval(dd22,ic2)
         D23ltc = Dval(dd23,ic1) * 2.d0 ! + Dval(dd23,ic2)
         D33ltc = Dval(dd33,ic1) * 2.d0 ! + Dval(dd33,ic2)
      enddo
#endif
      return
      end
******************************************************************
**** FOR COLLIER
#ifdef COLLIER
      subroutine Cfunscollier_withphmass(mph2)
      use collier
      double complex C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu
      double complex C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte
      double complex C00ltmu,C00lte
      double precision mph2,mph2tmp
      integer n,k
#include "invariants.h"
      common/vertexfunswphmass/
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte      
      data 
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte
     .     /14*(0.d0,0.d0)/

      integer rank ! max rank
      parameter (rank = 2)
      double complex cnmu(0:rank/2,0:rank,0:rank)
      double complex cnmuuv(0:rank/2,0:rank,0:rank)
      double complex cne(0:rank/2,0:rank,0:rank)
      double complex cneuv(0:rank/2,0:rank,0:rank)

      double complex cparmu(6),cpare(6)
      common/ccparwphmass/cparmu,cpare,ifirst
      data ifirst /0/

      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
c      common/colorandch2factors/colch2

      
      if (ifirst.eq.0) then
         cparmu(1)  = cuno*mm2
         cparmu(3)  = cuno*mm2
         cparmu(5)  = cuno*mm2
         cparmu(6)  = cuno*mm2

         cpare(1)  = cuno*me2
         cpare(3)  = cuno*me2
         cpare(5)  = cuno*me2
         cpare(6)  = cuno*me2
         ifirst    = 1
      endif
      cparmu(2) = cuno*(2.d0*mm2-2.d0*p1p3)
      cpare(2)  = cuno*(2.d0*me2-2.d0*p2p4)      
      
      C0ltmu  = czero
      C1ltmu  = czero
      C00ltmu = czero
      C11ltmu = czero
      C12ltmu = czero
      C0lte  = czero
      C1lte  = czero
      C00lte = czero
      C11lte = czero
      C12lte = czero
      
      do k = 1,1

         mph2tmp   = mph2
         cparmu(4) = cuno*mph2tmp
         cpare(4)  = cuno*mph2tmp
         
         if (QRmu.gt.0) then
! muon      
            call C_cll(cnmu,cnmuuv,
     .           cparmu(1),
     .           cparmu(2),
     .           cparmu(3),
     .           cparmu(4),
     .           cparmu(5),
     .           cparmu(6),
     .           rank)
            C0ltmu  = cnmu(0,0,0)
            C1ltmu  = cnmu(0,1,0)
            C00ltmu = cnmu(1,0,0)
            C11ltmu = cnmu(0,2,0)
            C12ltmu = cnmu(0,1,1)
         endif
! electron
         if (QRe.gt.0) then
            call C_cll(cne,cneuv,
     .           cpare(1),
     .           cpare(2),
     .           cpare(3),
     .           cpare(4),
     .           cpare(5),
     .           cpare(6),
     .           rank)
            C0lte  = cne(0,0,0)
            C1lte  = cne(0,1,0)
            C00lte = cne(1,0,0)
            C11lte = cne(0,2,0)
            C12lte = cne(0,1,1)
         endif
      enddo

      return
      end
**************************************************************************
      subroutine Dfunscollier_withphmass(mph2)
      use collier
      double precision mph2,mph2tmp
      integer n,k
      double complex D0ltd,D1ltd,D2ltd,D3ltd
      double complex D0ltc,D1ltc,D2ltc,D3ltc
      double complex D00ltd,D11ltd,D12ltd,D13ltd
      double complex D22ltd,D23ltd,D33ltd
      double complex D00ltc,D11ltc,D12ltc,D13ltc
      double complex D22ltc,D23ltc,D33ltc
      common/boxesfunswphmass/
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,     
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc

      integer rank ! max rank
      parameter (rank = 2)
      double complex dnd1(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnduv1(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnc1(0:rank/2,0:rank,0:rank,0:rank)
      double complex dncuv1(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnd2(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnduv2(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnc2(0:rank/2,0:rank,0:rank,0:rank)
      double complex dncuv2(0:rank/2,0:rank,0:rank,0:rank)

#include "invariants.h"
      double complex d0cll,cpard(10),cparc(10)
      common/dcparwphmass/cpard,cparc,ifirst
      data ifirst /0/
      double precision colch2(5),an
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
      
      if (ifirst.eq.0) then
         cpard(1)  = cuno*mm2
         cpard(2)  = cuno*mm2
         cpard(3)  = cuno*me2
         cpard(4)  = cuno*me2
         cpard(8)  = cuno*mm2
         cpard(10) = cuno*me2

         cparc = cpard
         
         ifirst   = 1
      endif

      cpard(5) = cuno*(2.d0*mm2-2.d0*p1p3)
      cpard(6) = cuno*(mm2+me2+2.d0*p3p4)
      cpard(9)  = cuno*0.d0
      cparc(5) = cuno*(2.d0*mm2-2.d0*p1p3)
      cparc(6) = cuno*(mm2+me2-2.d0*p2p3)
      cparc(9)  = cuno*0.d0

      D0ltd = czero
      D1ltd = czero
      D2ltd = czero
      D3ltd = czero
      D0ltc = czero
      D1ltc = czero
      D2ltc = czero
      D3ltc = czero
      D00ltd = czero
      D11ltd = czero
      D12ltd = czero
      D13ltd = czero
      D22ltd = czero
      D23ltd = czero
      D33ltd = czero
      D00ltc = czero
      D11ltc = czero
      D12ltc = czero
      D13ltc = czero
      D22ltc = czero
      D23ltc = czero
      D33ltc = czero

      do k = 1,1
         mph2tmp = mph2
      
         cpard(7)  = cuno*mph2tmp
         cparc(7)  = cuno*mph2tmp
      
! direct box      
         call D_cll(dnd1,dnduv1,
     .        cpard(1),
     .        cpard(2),
     .        cpard(3),
     .        cpard(4),
     .        cpard(5),
     .        cpard(6),
     .        cpard(7),
     .        cpard(8),
     .        cpard(9),
     .        cpard(10),rank)
     
! crossed box      
         call D_cll(dnc1,dncuv1,
     .        cparc(1),
     .        cparc(2),
     .        cparc(3),
     .        cparc(4),
     .        cparc(5),
     .        cparc(6),
     .        cparc(7),
     .        cparc(8),
     .        cparc(9),
     .        cparc(10),rank)

!! the second box is not needed, it's enough to multiply by two the first!      
         D0ltd  = dnd1(0,0,0,0) * 2.d0 ! + dnd2(0,0,0,0)
         D1ltd  = dnd1(0,1,0,0) * 2.d0 ! + dnd2(0,1,0,0)
         D2ltd  = dnd1(0,0,1,0) * 2.d0 ! + dnd2(0,0,1,0)
         D3ltd  = dnd1(0,0,0,1) * 2.d0 ! + dnd2(0,0,0,1)
         D00ltd = dnd1(1,0,0,0) * 2.d0 ! + dnd2(1,0,0,0)
         D11ltd = dnd1(0,2,0,0) * 2.d0 ! + dnd2(0,2,0,0) 
         D12ltd = dnd1(0,1,1,0) * 2.d0 ! + dnd2(0,1,1,0)
         D13ltd = dnd1(0,1,0,1) * 2.d0 ! + dnd2(0,1,0,1)
         D22ltd = dnd1(0,0,2,0) * 2.d0 ! + dnd2(0,0,2,0)
         D23ltd = dnd1(0,0,1,1) * 2.d0 ! + dnd2(0,0,1,1)
         D33ltd = dnd1(0,0,0,2) * 2.d0 ! + dnd2(0,0,0,2)

         D0ltc  = dnc1(0,0,0,0) * 2.d0 ! + dnc2(0,0,0,0)
         D1ltc  = dnc1(0,1,0,0) * 2.d0 ! + dnc2(0,1,0,0)
         D2ltc  = dnc1(0,0,1,0) * 2.d0 ! + dnc2(0,0,1,0)
         D3ltc  = dnc1(0,0,0,1) * 2.d0 ! + dnc2(0,0,0,1)
         D00ltc = dnc1(1,0,0,0) * 2.d0 ! + dnc2(1,0,0,0)
         D11ltc = dnc1(0,2,0,0) * 2.d0 ! + dnc2(0,2,0,0)
         D12ltc = dnc1(0,1,1,0) * 2.d0 ! + dnc2(0,1,1,0)
         D13ltc = dnc1(0,1,0,1) * 2.d0 ! + dnc2(0,1,0,1)
         D22ltc = dnc1(0,0,2,0) * 2.d0 ! + dnc2(0,0,2,0)
         D23ltc = dnc1(0,0,1,1) * 2.d0 ! + dnc2(0,0,1,1)
         D33ltc = dnc1(0,0,0,2) * 2.d0 ! + dnc2(0,0,0,2)
      enddo
      return
      end
*******************************
      subroutine getBsforctscollier_withphmass(m2,mph2,B1o,dB0o,dB1o)
      use collier
      implicit none
      double complex B1o,dB0o,dB1o
      double precision m2,mph2
      integer rank
      parameter (rank=1)
      double complex b(0:rank/2,0:rank),buv(0:rank/2,0:rank)
      
      call B_cll(b,buv,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,mph2*(1.d0,0.d0),1)
      B1o = b(0,1)
      
      call DB0_cll(dB0o,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,mph2*(1.d0,0.d0))
      call DB1_cll(dB1o,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,mph2*(1.d0,0.d0))

      return
      end
#endif


