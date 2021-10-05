      subroutine initmuescattering(s)
      implicit double precision (a-h,o-z)
#include "invariants.h"
      common/parameters/ame,ammu,convfac,alpha,pi
      character*6 ord
      integer irecola
      common/recolamode/irecola
      common/qedORDER/ord
      common/nphot_mode/nphotmode
      common/photoncutoff/egmin,egmin2
      double precision lambda
      external lambda
      common/epssoft/eps
      
      ifixedorder = 0! if VP is at fixed-order alpha (1), or resummed (also in NLO diagrams) (0)

*** Qmu,QR[mu,e] read in userinterface
      
      if (ord.eq.'born')   ioneloop = 0
      if (ord.eq.'alpha')  ioneloop = 1
      if (ord.eq.'alpha2') ioneloop = 1
         
      ec2 = 4.d0*pi*alpha       ! electron charge^2
      ec4 = ec2*ec2
      ec  = sqrt(ec2)

      mm  = ammu
      mm2 = mm*mm
      me  = ame
      me2 = me*me
      
      mm4 = mm2*mm2
      mm6 = mm2*mm4
      me4 = me2*me2
      me6 = me2*me4

      pi216   = pi*pi*16.d0 ! needed in the cts...
      pi216m1 = 1.d0/pi216

*** default      
cc      print*,' egmin/egmin2 governed by eps in input! '
      egmin = eps * 0.5d0*sqrt(s)
c      egmin = 1d-6 ! cfr fulvio
      
      egmin2 = egmin ! it must be egmin2 >= egmin
      
c      if (ioneloop.eq.1) then
      ! initialize anyway looptools/collier...
      call initloops(s)
c      endif

      irecola = 1
      call recola_init
      
      return
      end
**********************************************************************
      subroutine initloops(s)
      implicit double precision (a-h,m,o-z)
#include "looptools.h"
      integer ifirst
      data ifirst /0/
      double precision egmin
      common/photoncutoff/egmin,egmin2
      common/initloop/ifirst
      RealType deltauv,mudim
      double precision phmass
      common/photonmasslambda/phmass
      double precision duvcoll,mucoll
      common/uvparams/duvcoll,mucoll

      common/parameters/ame,ammu,convfac,alpha,pi

      RealType uquad      
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif

*** photon mass read in userinterface
      
      deltauv = 0.d0
      mudim   = 1.d0
      
      if (ifirst.eq.0) then
         call initlooplibrary
         ifirst = 1      
      endif

* both initialized always, because I uset getlambda from LT in soft photon factor
      duvcoll = deltauv
      mucoll  = mudim
      phmcoll = phmass
#ifdef COLLIER
#ifndef RECOLA
!     Collier
      call setregcollier(duvcoll,mucoll,phmcoll)
#endif
#else
!     Looptools
      call setdelta(deltauv)
      call setmudim(mudim)

c      call setversionkey(0*KeyBget + 0*KeyC0   + 0*KeyD0
c     .                 + 0*KeyD0C  + 0*KeyEget + 0*KeyEgetC)
c      call setversionkey(1*KeyAll)
c      call setversionkey(3*KeyEget)
! default 9
c      print*,'WARNING: changed warning digits from 9 to 12 in LT'
c       call setwarndigits(14) ! seems to be safe to not see may "rotation" prints from LT

c     call setdebugkey(DebugAll)
c      call setdebugkey(DebugE)      
c      call setmaxdev(1d-5)
c      print*,getmaxdev()

      call setlambda(uquad*phmass*phmass)
#endif      
      return
      end
****************************************************
c      subroutine clearthecache ! only looptools...
c#include "looptools.h"
c#ifndef COLLIER      
c      call clearcache ! I clear anyway the LT cache...
c#endif
c      return
c      end
*****************************************************
      subroutine initlooplibrary
#include "looptools.h"
#ifdef COLLIER
#ifndef RECOLA
      call initcollier
#endif
#else      
      call ltini
#endif
      return
      end
*****************************************************
      subroutine exitlooplibraries
#include "looptools.h"
!     if (ioneloop.eq.1)
      ! exiting anyway
#ifndef COLLIER      
      call ltexi
#endif      
      return
      end
******************************************************
***   COLLIER ROUTINES/FUNCTIONS
******************************************************
#ifdef COLLIER
      subroutine initcollier
      use collier
      double precision acc
      character*200  path
      common/cmnpath/path
      
c     call Init_cll(4,2,path(1:len_trim(path))//'cll_out')      
      call Init_cll(5,3,path(1:len_trim(path))//'cll_out')
      call SetMode_cll(1)       ! 1 is COLI, 2 DD, 3 both. 1 is faster

      call InitCacheSystem_cll(5,5) ! init 5 caches, 1 for 2->2, 2 for 2->3, 3,4,5 for leptonic bubbles
      
c      call SetCritAcc_cll(1d-16)
c      call SetReqAcc_cll(1d-16)
c      call GetReqAcc_cll(acc)
c      print*,acc
c      stop
      return
      end
******************************************************
      subroutine setregcollier(duv,mu,m)
      use collier
      double precision duv,mu,m,ame
      double precision m2reg(1)
      call SetDeltaUV_cll(duv)
      call SetMuUV2_cll(mu)
      call SetmuIR2_cll(m*m)
      return
      end
****************************************************
c      subroutine clearthecachecollier
c      use collier
c      integer ifirst,ncache,n
c      data ifirst /0/
c      common/cachecollier/ifirst,ncache
c      if (ifirst.eq.0) then
c         ncache = 1
c         call InitCacheSystem_cll(ncache,5)
c         ifirst = 1
c      endif
c      call InitEvent_cll(ncache)
c      return
c      end
*****************************************************
      subroutine getBsforctscollier(m2,B1o,dB0o,dB1o)
      use collier
      implicit none
      double complex B1o,dB0o,dB1o
      double precision m2
      integer rank
      parameter (rank=1)
      double complex b(0:rank/2,0:rank),buv(0:rank/2,0:rank)

      call SwitchOffCacheSystem_cll
      
      call B_cll(b,buv,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,(0.d0,0.d0),1)
      B1o = b(0,1)
      
      call DB0_cll(dB0o,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,(0.d0,0.d0))
      call DB1_cll(dB1o,(1.d0,0.d0)*m2,(1.d0,0.d0)*m2,(0.d0,0.d0))

      call SwitchOnCacheSystem_cll
      
      return
      end
**************************************************************************
      subroutine Dfunscollier
      use collier
      double complex D0ltd,D1ltd,D2ltd,D3ltd
      double complex D0ltc,D1ltc,D2ltc,D3ltc
      double complex D00ltd,D11ltd,D12ltd,D13ltd
      double complex D22ltd,D23ltd,D33ltd
      double complex D00ltc,D11ltc,D12ltc,D13ltc
      double complex D22ltc,D23ltc,D33ltc
      common/boxesfuns/
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,     
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc

      integer rank ! max rank
      parameter (rank = 2)
      double complex dnd(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnduv(0:rank/2,0:rank,0:rank,0:rank)
      double complex dnc(0:rank/2,0:rank,0:rank,0:rank)
      double complex dncuv(0:rank/2,0:rank,0:rank,0:rank)

#include "invariants.h"
      double complex d0cll,cpard(10),cparc(10)
      common/dcpar/cpard,cparc,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         cpard(1)  = cuno*mm2
         cpard(2)  = cuno*mm2
         cpard(3)  = cuno*me2
         cpard(4)  = cuno*me2
         cpard(7)  = cuno*0.d0
         cpard(8)  = cuno*mm2
         cpard(9)  = cuno*0.d0
         cpard(10) = cuno*me2

         cparc = cpard
         
         ifirst   = 1
      endif
      cpard(5) = cuno*(2.d0*mm2-2.d0*p1p3)
      cpard(6) = cuno*(mm2+me2+2.d0*p3p4)
      cparc(5) = cuno*(2.d0*mm2-2.d0*p1p3)
      cparc(6) = cuno*(mm2+me2-2.d0*p2p3)

! direct box      
      call D_cll(dnd,dnduv,
     .     cpard(1),
     .     cpard(2),
     .     cpard(3),
     .     cpard(4),
     .     cpard(5),
     .     cpard(6),
     .     cpard(7),
     .     cpard(8),
     .     cpard(9),
     .     cpard(10),rank)

! crossed box      
      call D_cll(dnc,dncuv,
     .     cparc(1),
     .     cparc(2),
     .     cparc(3),
     .     cparc(4),
     .     cparc(5),
     .     cparc(6),
     .     cparc(7),
     .     cparc(8),
     .     cparc(9),
     .     cparc(10),rank)

      D0ltd  = dnd(0,0,0,0)
      D1ltd  = dnd(0,1,0,0)
      D2ltd  = dnd(0,0,1,0)
      D3ltd  = dnd(0,0,0,1)
      D00ltd = dnd(1,0,0,0)
      D11ltd = dnd(0,2,0,0)
      D12ltd = dnd(0,1,1,0)
      D13ltd = dnd(0,1,0,1)
      D22ltd = dnd(0,0,2,0)
      D23ltd = dnd(0,0,1,1)
      D33ltd = dnd(0,0,0,2)

      D0ltc  = dnc(0,0,0,0)
      D1ltc  = dnc(0,1,0,0)
      D2ltc  = dnc(0,0,1,0)
      D3ltc  = dnc(0,0,0,1)
      D00ltc = dnc(1,0,0,0)
      D11ltc = dnc(0,2,0,0)
      D12ltc = dnc(0,1,1,0)
      D13ltc = dnc(0,1,0,1)
      D22ltc = dnc(0,0,2,0)
      D23ltc = dnc(0,0,1,1)
      D33ltc = dnc(0,0,0,2)

      return
      end
**************************************************************************
******************************************************************
      subroutine Cfunscollier
      use collier
      double complex C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu
      double complex C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte
      double complex C00ltmu,C00lte
#include "invariants.h"
      common/vertexfuns/
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
      common/ccpar/cparmu,cpare,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         cparmu(1)  = cuno*mm2
         cparmu(3)  = cuno*mm2
         cparmu(4)  = cuno*0.d0
         cparmu(5)  = cuno*mm2
         cparmu(6)  = cuno*mm2

         cpare(1)  = cuno*me2
         cpare(3)  = cuno*me2
         cpare(4)  = cuno*0.d0
         cpare(5)  = cuno*me2
         cpare(6)  = cuno*me2
         ifirst    = 1
      endif
      cparmu(2) = cuno*(2.d0*mm2-2.d0*p1p3)
      cpare(2)  = cuno*(2.d0*me2-2.d0*p2p4)      

      if (QRmu.gt.0) then
! muon      
         call C_cll(cnmu,cnmuuv,
     .        cparmu(1),
     .        cparmu(2),
     .        cparmu(3),
     .        cparmu(4),
     .        cparmu(5),
     .        cparmu(6),
     .        rank)
         C0ltmu  = cnmu(0,0,0)
         C1ltmu  = cnmu(0,1,0)
         C00ltmu = cnmu(1,0,0)
         C11ltmu = cnmu(0,2,0)
         C12ltmu = cnmu(0,1,1)
      endif
! electron
      if (QRe.gt.0) then
         call C_cll(cne,cneuv,
     .        cpare(1),
     .        cpare(2),
     .        cpare(3),
     .        cpare(4),
     .        cpare(5),
     .        cpare(6),
     .        rank)
         C0lte  = cne(0,0,0)
         C1lte  = cne(0,1,0)
         C00lte = cne(1,0,0)
         C11lte = cne(0,2,0)
         C12lte = cne(0,1,1)
      endif
      return
      end
*** END  COLLIER ROUTINES/FUNCTIONS
#endif
******************************************************
      double precision function radmue(p1,p2,p3,p4,k) !,ce,cmu)
      implicit none !double precision (a-h,o-z)
#include "invariants.h"
      double complex sp1(2,4),sp2(2,4),sp3(2,4),sp4(2,4)
      double complex s1l(4),s2l(4),s3l(4),s4l(4)
      double precision t13m1,t24m1,dp1km1,dp2km1,dp3km1,dp4km1
      double precision k(0:3),e1(0:3),e2(0:3),eps(0:3)
      double complex evkcec5mu(0:3),evkcmuc5e(0:3)
      double complex rad(0:3)
      double complex kce,kcmu,cover
      double complex Z6_(0:3),Z7_(0:3),Z8_(0:3),Z9_(0:3)
      double precision Z1_(0:3),Z2_(0:3),Z3_(0:3),Z4_(0:3),Z5_(0:3)
      double precision Z14_(0:3),Z15_(0:3),Z16_(0:3)
      double precision Z10_,Z11_,Z12_,Z13_
      double precision coverm2,csi(1)
      double precision sqm,pr1,pr2,w
      double complex ct(0:3),c5t(0:3),cIt
      double complex cmum(0:3,2,2)
      double complex c5mum(0:3,2,2)
      double complex cImum(2,2)
      double complex cem(0:3,2,2)
      double complex c5em(0:3,2,2)
      double complex cIem(2,2)
      integer hh,ifirst,h1b,h2b,h3b,h4b,icmu,ice
      common/ifirstrad/cover,coverm2,ifirst
      double complex vpolc
      external vpolc
      double precision radmuesum,dah13,dah24
      external radmuesum
      double precision radmuesumsoft,radmuesumonlye,radmuesumonlyeSD
      external radmuesumsoft,radmuesumonlye,radmuesumonlyeSD
      data ifirst /0/
      integer kkk
      
      if (ifirst.eq.0) then
         cover   = ec2*ec*im
         coverm2 = cover*conjg(cover)
         ifirst  = 1
      endif
      
      radmue = radmuesum(p1,p2,p3,p4,k)
      
      return
*********************************************************
c     goto 666
      
      t13m1 = 1.d0/dot(p1-p3,p1-p3)
      t24m1 = 1.d0/dot(p2-p4,p2-p4)

      if (ifixedorder.eq.0) then    
         t13m1 = t13m1 * vpolc(1.d0/t13m1)
         t24m1 = t24m1 * vpolc(1.d0/t24m1)
      endif

      dp1km1 = 0.5d0/dot(p1,k)
      dp2km1 = 0.5d0/dot(p2,k)
      dp3km1 = 0.5d0/dot(p3,k)
      dp4km1 = 0.5d0/dot(p4,k)

      call get_us(p2,me,sp2) ! for electron-
      call get_us(p4,me,sp4) ! for electron-

      radmue = 0.d0

      Z1_=p1
      Z2_=p3
      Z3_=p2
      Z4_=p4
      Z5_=k

      Z10_=dp3km1 - dp1km1
      Z11_=dp1km1 + dp3km1
      Z12_=dp4km1 - dp2km1
      Z13_=dp2km1 + dp4km1
      
      if (Qmu.lt.0.d0) then     

         call get_us(p1,mm,sp1) ! for mu-
         call get_us(p3,mm,sp3) ! for mu-
         
         Z14_=dp1km1*Z1_
         Z15_= - Z2_*dp3km1
         Z14_=Z14_ + Z15_
         Z14_=2*Z14_
         Z15_=dp2km1*Z3_
         Z16_= - Z4_*dp4km1
         Z15_=Z15_ + Z16_
         Z15_=2*Z15_

c$$$         do h1 = 1,1
c$$$            do h3 = 1,2
c$$$               call getcurrentnew(h3,h1,p3,p1,mm,mm,cmu,cImu,c5mu)               
c$$$               cmum(:,h3,h1) = cmu
c$$$               c5mum(:,h3,h1) = c5mu               
c$$$               cImum(h3,h1) = cImu
c$$$            enddo
c$$$         enddo
c$$$         do h4 = 1,2
c$$$            do h2 = 1,2
c$$$               call getcurrentnew(h4,h2,p4,p2,me,me,ce,cIe,c5e)
c$$$               cem(:,h4,h2) = ce
c$$$               c5em(:,h4,h2) = c5e
c$$$               cIem(h4,h2) = cIe
c$$$            enddo
c$$$         enddo
c            print*,

         h1b = -1
         h2b = -1
         h3b = -1
         h4b = -1         
         do hh = 0,7            !! the 2 below is because I'm summing over half helicities
c            h1  = mod(hh,16)/8 + 1
c            h2  = mod(hh, 8)/4 + 1
c            h3  = mod(hh, 4)/2 + 1
c            h4  = mod(hh, 2)   + 1

            h1  = mod(hh,16)/8 + 1
            h3  = mod(hh, 8)/4 + 1
            h2  = mod(hh, 4)/2 + 1
            h4  = mod(hh, 2)   + 1

            if (h1.ne.h1b) then
               s1l = sp1(h1,:)
               h1b = h1
            endif
            if (h2.ne.h2b) then
               s2l = sp2(h2,:)
               h2b = h2
            endif
            if (h3.ne.h3b) then
               s3l = sp3(h3,:)
               h3b = h3
            endif
            if (h4.ne.h4b) then
               s4l = sp4(h4,:)           
               h4b = h4
            endif

            if (hh.eq.0.or.hh.eq.4) then
               call getcurrent(s3l,s1l,cmu,cImu,c5mu)
            endif
            call getcurrent(s4l,s2l,ce,cIe,c5e)
c           call getcurrentnew(h3,h1,p3,p1,mm,mm,cmu,cImu,c5mu)
c oppure
c            cmu  = cmum(:,h3,h1)
c            c5mu = c5mum(:,h3,h1)
c            cImu = cImum(h3,h1)
            
c            call getcurrentnew(h4,h2,p4,p2,me,me,ce,cIe,c5e)
c oppure
c            ce  = cem(:,h4,h2)
c            c5e = c5em(:,h4,h2)
c            cIe = cIem(h4,h2)

c     print*,cImu,cIe,hh,h3,h1,h4,h2
            
c            print*,abs(c5mu-c5t),h3,h1,hh
c            stop
           
            cmuce = dotcmom(cmu,ce)

            kce   = dot1cmom(k,ce)
            kcmu  = dot1cmom(k,cmu)      
            
            call getev(evkcec5mu,k,ce,c5mu)
            call getev(evkcmuc5e,k,cmu,c5e)
            
** this because in the form I have e_(k,ce,c5mu,sig) but in my routine sig is in the first place
            evkcec5mu = -evkcec5mu 
            evkcmuc5e = -evkcmuc5e

            rad=t24m1*evkcec5mu*Z10_ + t24m1*kcmu*ce*Z11_
     &           - t24m1*kce*cmu*Z11_ + 
     &           t13m1*evkcmuc5e*Z12_ - t13m1*kcmu*ce*Z13_
     &           + t13m1*kce*cmu*Z13_+ 
     &           cmuce*t24m1*Z14_ - cmuce*t24m1*Z5_*Z11_
     &           + cmuce*t13m1*Z15_- cmuce*t13m1*Z5_*Z13_         
            
            radmue = radmue - dotcmom(rad,conjg(rad))         
         enddo

      else ! if Qmu > 0

         call get_vs(p3,mm,sp1) ! these are for mu+!!
         call get_vs(p1,mm,sp3) ! these are for mu+!!

         Z14_= - dp1km1*Z1_
         Z15_=Z2_*dp3km1
         Z14_=Z14_ + Z15_
         Z14_=2*Z14_
         Z15_=dp2km1*Z3_
         Z16_= - Z4_*dp4km1
         Z15_=Z15_ + Z16_
         Z15_=2*Z15_
         
         do hh = 0,7            !! the 2 below is because I'm summing over half helicities
            h1  = mod(hh,16)/8 + 1
            h2  = mod(hh, 8)/4 + 1
            h3  = mod(hh, 4)/2 + 1
            h4  = mod(hh, 2)   + 1
            s1l = sp1(h1,:)
            s2l = sp2(h2,:)
            s3l = sp3(h3,:)
            s4l = sp4(h4,:)
            
            call getcurrent(s4l,s2l,ce,cIe,c5e)
            call getcurrent(s3l,s1l,cmu,cImu,c5mu)
            
            cmuce = dotcmom(cmu,ce)
            kce   = dot1cmom(k,ce)
            kcmu  = dot1cmom(k,cmu)      
            
            call getev(evkcec5mu,k,ce,c5mu)
            call getev(evkcmuc5e,k,cmu,c5e)
            
** this because in the form I have e_(k,ce,c5mu,sig) but in my routine sig is in the first place
            evkcec5mu = -evkcec5mu 
            evkcmuc5e = -evkcmuc5e

            rad=t24m1*evkcec5mu*Z10_ - t24m1*kcmu*ce*Z11_
     &           + t24m1*kce*cmu*Z11_ + 
     &           t13m1*evkcmuc5e*Z12_ - t13m1*kcmu*ce*Z13_
     &         + t13m1*kce*cmu*Z13_ + 
     &           cmuce*t24m1*Z14_ + cmuce*t24m1*Z5_*Z11_
     &         + cmuce*t13m1*Z15_ - cmuce*t13m1*Z5_*Z13_
            
            radmue = radmue - dotcmom(rad,conjg(rad))         
         enddo

      endif
      radmue = 2.d0*coverm2*radmue
      return
      end
******************************************************
      double precision function radmuesum(p1,p2,p3,p4,k)
** summed over spin/elicities. it's much faster (~0.35)      
*      implicit none !double precision (a-h,o-z)
      implicit double precision (Z)
#include "invariants.h"
      double precision t13m1,t24m1,dp1km1,dp2km1,dp3km1,dp4km1
      double precision dp1km2,dp2km2,dp3km2,dp4km2
      double precision p1km1,p2km1,p3km1,p4km1,muemueg1L
      double precision k(0:3),t13m2,t24m2
      double complex cover
      double precision vpleps13(1:4)
      double precision vpleps24(1:4)
      double precision coverm2,dah13,dah24,pi
      double precision dalep13,dalep24
      double precision t13m1v(0:2),t24m1v(0:2)
      double precision t13m2v(0:2),t24m2v(0:2)
      double precision p1p22,p1p32,p1p42,p2p32,p2p42,p3p42
      common/ifirstradsum/cover,coverm2,ifirst
      double complex vpolc
      external vpolc
      integer ifirst
      data ifirst /0/
      double precision radmue_mu,radmue_e,radmue_int
*** for 1g1L
      double precision dalp24,dalp13,radmuefermionloop
      common/dalphasfor1g1Lcfrrecola/dalp24,dalp13,radmuefermionloop
***
      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask

      
      if (ifirst.eq.0) then
         cover   = ec2*ec*im
         coverm2 = cover*conjg(cover)
         ifirst = 1
         radmuefermionloop = 0.d0
         dalp24 = 0.d0
         dalp13 = 0.d0
      endif

      reweightLO = 0.d0

      call loadinvariants1g(p1,p2,p3,p4,k)
      
      t13m1v(0) = 1.d0/t13
      t24m1v(0) = 1.d0/t24
      
      if (ifixedorder.eq.0) then
         call getVP1loop(1.d0/t13m1v(0),dal,dah)
         dah13   = dah

         dalep13 = dot_product(vpvect,bubblemask)
         t13m1v(1) = t13m1v(0) * 1.d0/(1.d0 - dal - dah)
         t13m1v(2) = t13m1v(0) * 1.d0/(1.d0 - dal)

         call getVP1loop(1.d0/t24m1v(0),dal,dah)
         dah24   = dah
         dalep24 = dot_product(vpvect,bubblemask)
         t24m1v(1) = t24m1v(0) * 1.d0/(1.d0 - dal - dah)
         t24m1v(2) = t24m1v(0) * 1.d0/(1.d0 - dal)
      else
         t13m1v(1:2) = t13m1v(0)
         t24m1v(1:2) = t24m1v(0)
      endif
      
      t13m2v = t13m1v*t13m1v
      t24m2v = t24m1v*t24m1v
      
      
      dp1km1 = 0.5d0/p1k
      dp2km1 = 0.5d0/p2k
      dp3km1 = 0.5d0/p3k
      dp4km1 = 0.5d0/p4k

      dp1km2 = dp1km1*dp1km1
      dp2km2 = dp2km1*dp2km1
      dp3km2 = dp3km1*dp3km1
      dp4km2 = dp4km1*dp4km1
      
      p1p22 = p1p2*p1p2
      p1p32 = p1p3*p1p3
      p1p42 = p1p4*p1p4
      p2p32 = p2p3*p2p3
      p2p42 = p2p4*p2p4
      p3p42 = p3p4*p3p4
      
      if (Qmu.lt.0.d0) then     
cc         include 'oneloop/radmuesumminus.f'

      Z1_=p1p4 - p4k
      Z1_=Z1_*p2p3
      Z2_=p2k*p3p4
      Z3_=p1p2*p3p4
      Z4_= - Z1_ + Z2_ - Z3_
      Z4_=Z4_*mm2*dp1km2
      Z5_=p3p4 + p4k
      Z5_=Z5_*p1p2
      Z6_=p2k*p1p4
      Z5_=Z5_ + Z6_
      Z7_=p2p3*p1p4
      Z8_= - Z7_ - Z5_
      Z8_=mm2*Z8_
      Z9_=p1p3*mm2
      Z9_= - Z9_ + 2*mm4
      Z10_=mm2*p1k
      Z10_=Z10_ - Z9_
      Z10_=me2*Z10_
      Z11_=p2p4*mm4
      Z8_=Z10_ + Z11_ + Z8_
      Z8_=dp3km2*Z8_
      Z10_= - mm2*p3k
      Z9_=Z10_ - Z9_
      Z9_=me2*Z9_
      Z9_=Z9_ + Z11_
      Z9_=dp1km2*Z9_
      Z4_=Z8_ + Z4_ + Z9_
      Z8_=2*p1p4
      Z9_= - p3p4 - Z8_
      Z9_=p1p2*Z9_
      Z10_=2*p1p3
      Z11_=Z10_ + mm2
      Z12_= - p3k + Z11_
      Z12_=me2*Z12_
      Z1_=Z12_ - Z1_ + Z2_ + Z9_
      Z1_=dp1km1*Z1_
      Z2_=2*p3p4
      Z9_=Z2_ + p4k
      Z9_=Z9_*p1p2
      Z12_=2*mm2
      Z13_=Z12_*p2p4
      Z9_=Z9_ - Z13_
      Z14_=Z8_ - p4k
      Z14_=Z14_*p2p3
      Z15_= - p3p4 + p1p4
      Z15_=p2k*Z15_
      Z15_=Z14_ + Z15_ + Z9_
      Z15_=p1p3*Z15_
      Z16_=2*me2
      Z17_=p1p3*Z12_
      Z17_= - p1p32 + Z17_
      Z17_=Z17_*Z16_
      Z18_= - p2k*p4k*Z12_
      Z15_=Z17_ + Z18_ + Z15_
      Z15_=dp1km1*Z15_
      Z17_=Z2_ + p1p4
      Z17_=p2p3*Z17_
      Z11_= - p1k - Z11_
      Z11_=me2*Z11_
      Z5_=2*Z15_ + Z11_ + Z17_ + Z5_
      Z5_=dp3km1*Z5_
      Z1_=Z5_ + 2*Z4_ + Z1_
      Z1_=64*Z1_
      Z4_=p3p4*p1k
      Z5_=p1p2 - p1p4
      Z11_=p3k*Z5_
      Z15_=p1k + Z8_
      Z15_=p2p3*Z15_
      Z11_=Z15_ + Z11_ - Z4_ + 2*Z3_
      Z11_=p2p4*Z11_
      Z15_= - p2p42*Z12_
      Z11_=Z15_ + Z11_
      Z15_=2*dp4km1
      Z11_=Z11_*Z15_
      Z17_=mm2*p4k
      Z18_=Z17_ - Z4_
      Z19_=p3k*p1p4
      Z20_=Z3_ + Z18_ - Z19_
      Z21_= - p3k*p1k
      Z22_= - p1p3*p2p4
      Z21_=Z22_ + Z21_ + Z13_
      Z21_=dp4km1*Z21_
      Z21_=mm2 + 4*Z21_
      Z21_=me2*Z21_
      Z22_=2*p1p2
      Z23_= - p1p4 - Z22_
      Z23_=p2p3*Z23_
      Z11_=Z21_ + Z11_ + Z23_ + Z13_ - Z20_
      Z11_=dp2km1*Z11_
      Z21_=p1p4 + p1k
      Z21_=Z21_*p2p3
      Z23_=p3k*p1p2
      Z21_=Z21_ + Z23_ + Z3_
      Z24_=mm2*p2k
      Z25_=Z13_ + Z24_
      Z26_=Z8_*p3p4
      Z27_=Z26_ - Z25_ + Z21_
      Z27_=dp4km1*Z27_
      Z28_=p2p4*mm2
      Z21_=Z28_ + Z24_ - Z21_
      Z21_=dp4km2*Z21_
      Z24_= - dp4km1*mm2
      Z21_=2*Z21_ + Z24_
      Z21_=me2*Z21_
      Z20_= - Z7_ + Z28_ - Z20_
      Z16_=dp2km2*Z20_*Z16_
      Z20_= - dp4km2 - dp2km2
      Z12_=Z12_ - p1p3
      Z12_=me4*Z12_*Z20_
      Z11_=Z11_ + 2*Z12_ + Z16_ + Z27_ + Z21_
      Z11_=64*Z11_
      Z12_=p3p4*Z25_
      Z16_= - Z4_ - Z26_
      Z16_=p2p3*Z16_
      Z6_=Z6_*p3p4
      Z20_= - p3p42*Z22_
      Z12_=Z16_ + Z20_ - Z6_ + Z12_
      Z12_=Z12_*Z15_
      Z16_=4*mm2
      Z10_=Z16_ - Z10_
      Z20_= - p1k + Z10_
      Z20_=me2*Z20_
      Z9_=Z20_ + Z18_ + Z9_
      Z9_=p2p3*Z9_
      Z18_=p2p32*Z8_
      Z9_=Z18_ + Z9_
      Z9_=dp2km1*Z9_
      Z18_= - p3p4*Z16_
      Z2_=p1p3*Z2_
      Z2_=Z2_ + Z4_ + Z18_
      Z4_=Z15_*me2
      Z2_=Z2_*Z4_
      Z18_=p1p2 + p1p4
      Z20_= - Z18_ - p2k + p3p4 - p4k
      Z20_=Z20_*mm2
      Z21_=p1p3*p3p4
      Z22_=p1p3 + p1p4 - mm2
      Z22_=p2p3*Z22_
      Z2_=2*Z9_ + Z2_ + Z12_ + Z22_ + Z21_ - Z3_ - Z20_
      Z2_=dp3km1*Z2_
      Z9_= - p1p4*Z25_
      Z12_=Z3_*Z8_
      Z21_=Z23_*p1p4
      Z22_=p2p3*p1p42
      Z6_=2*Z22_ + Z21_ - Z6_ + Z12_ + Z9_
      Z6_=Z6_*Z15_
      Z9_=p1p4*Z16_
      Z8_= - p1p3*Z8_
      Z8_=Z8_ + Z19_ + Z9_
      Z4_=Z8_*Z4_
      Z8_= - p1p3*Z18_
      Z9_= - p1p4 - mm2
      Z9_=p2p3*Z9_
      Z4_=Z4_ + Z6_ + Z9_ + Z8_ + Z3_ - Z20_
      Z4_=dp1km1*Z4_
      Z6_= - Z17_ - Z14_ + Z13_
      Z6_=p1p2*Z6_
      Z8_= - p1p2*Z10_
      Z8_= - Z23_ + Z8_
      Z8_=me2*Z8_
      Z6_=Z8_ + Z21_ + Z6_
      Z8_=p1p22*p3p4
      Z6_= - 4*Z8_ + 2*Z6_
      Z6_=dp1km1*Z6_
      Z5_= - p3k - Z5_ - p2p3 + p3p4 - p1k
      Z5_=Z5_*me2
      Z8_= - p2p4*p1p2
      Z9_= - p1p4 - p2p4
      Z9_=p2p3*Z9_
      Z6_= - Z5_ + Z9_ + Z3_ + Z8_ + Z6_
      Z6_=dp2km1*Z6_
      Z8_=p1p4 + p3p4
      Z8_=p2p4*Z8_
      Z3_= - Z5_ + Z7_ - Z3_ + Z8_
      Z3_=dp4km1*Z3_
      Z2_=Z2_ + Z6_ + Z4_ - p2p3 - p3p4 - Z18_ + Z3_
      Z2_=64*Z2_

      radmuesum =  QRmu * t24m2v(1)*Z1_
     .     + QRe  * t13m2v(1)*Z11_
     .     + QRmu*QRe * t13m1v(1)*t24m1v(1)*Z2_

      if (ionlybubbles.gt.0) then
         radmuesum =  QRmu * t24m2v(0)*(2.d0*dalep24)*Z1_
     .        + QRe  * t13m2v(0)*(2.d0*dalep13)*Z11_
     .        + QRmu*QRe * t13m1v(0)*t24m1v(0)*(dalep13+dalep24)*Z2_
      endif
      
cfrpasserafael
c      radmuesum=t24m2v(0)*Z1_*2.d0*dah24 + t13m2v(0)*Z11_*2.d0*dah13
c     .     + t13m1v(0)*t24m1v(0)*Z2_*(dah13+dah24)
cfrpasserafael      

      extraweights(1)=
     .     QRmu*t24m2v(0)*Z1_ + QRe*t13m2v(0)*Z11_
     .     + QRmu*QRe * t13m1v(0)*t24m1v(0)*Z2_
      extraweights(2)=
     .     QRmu * t24m2v(2)*Z1_ + QRe * t13m2v(2)*Z11_
     .     + QRmu*QRe * t13m1v(2)*t24m1v(2)*Z2_      
      extraweights(1:2) = extraweights(1:2)/radmuesum
      
      else  ! if Qmu > 0
cc         include 'oneloop/radmuesumplus.f'

      Z1_=dp3km2 + dp1km2
      Z2_=mm2*Z1_
      Z3_=dp1km1*mm2
      Z4_= - 1 + 4*Z3_
      Z5_=dp3km1*Z4_
      Z2_=Z5_ + dp1km1 + Z2_
      Z2_=p1p3*Z2_
      Z1_=Z1_*mm4
      Z5_=dp3km1*mm2
      Z6_=2*mm2
      Z7_=Z6_*dp1km2
      Z8_=Z7_ + dp1km1
      Z9_= - p3k*Z8_
      Z10_=Z6_*dp3km2
      Z11_=Z10_ - dp3km1
      Z12_=p1k*Z11_
      Z13_=4*dp3km1
      Z14_= - p1p32*dp1km1*Z13_
      Z2_= - 4*Z1_ + Z14_ + Z12_ + 2*Z2_ - Z5_ + Z3_ + Z9_
      Z2_=me2*Z2_
      Z9_=4*dp1km1
      Z12_=p1p3*dp3km1
      Z14_=Z9_*Z12_
      Z15_=dp3km1 - dp1km1
      Z7_= - Z14_ - Z15_ + Z7_ + Z10_
      Z10_= - p2p3*Z7_
      Z14_=2*dp1km1
      Z16_= - p1p2*Z14_
      Z10_=Z10_ + Z16_
      Z10_=p1p4*Z10_
      Z7_= - p1p2*Z7_
      Z16_=2*dp3km1
      Z17_=Z16_*p2p3
      Z7_=Z17_ + Z7_
      Z7_=p3p4*Z7_
      Z18_=Z14_*Z12_
      Z11_=Z11_ - Z18_
      Z19_= - p1p4*Z11_
      Z8_=Z8_ - Z18_
      Z18_=p3p4*Z8_
      Z20_= - p4k*Z3_*Z13_
      Z18_=Z20_ + Z19_ + Z18_
      Z18_=p2k*Z18_
      Z8_=p2p3*Z8_
      Z11_= - p1p2*Z11_
      Z8_=Z8_ + Z11_
      Z8_=p4k*Z8_
      Z11_=Z3_*Z12_
      Z1_= - 2*Z11_ + Z1_
      Z11_=2*p2p4
      Z1_=Z1_*Z11_
      Z1_=Z1_ + Z18_ + Z2_ + Z8_ + Z10_ + Z7_
      Z1_=64*Z1_
      Z2_= - me2*p1k*p3k
      Z7_= - p2p42*mm2
      Z2_=Z2_ + Z7_
      Z2_=dp4km1*Z2_
      Z7_=p4k*mm2
      Z8_=me2*mm2
      Z10_=2*p1p2
      Z18_= - p2p3*Z10_
      Z19_=p3k - p2p3
      Z19_=p1p4*Z19_
      Z20_=p1k - p1p2
      Z20_=p3p4*Z20_
      Z2_=4*Z2_ + Z8_ - Z7_ + Z20_ + Z18_ + Z19_
      Z2_=dp2km1*Z2_
      Z7_=p1p2*p3k
      Z18_=p1k*p2p3
      Z19_=Z7_ + Z18_
      Z20_=p1p4*p2p3
      Z21_=p1p2 + 2*p1p4
      Z21_=p3p4*Z21_
      Z22_= - p2k*mm2
      Z20_=Z22_ - Z8_ + Z21_ + Z20_ + Z19_
      Z20_=dp4km1*Z20_
      Z21_= - dp4km2*Z19_
      Z22_=p3k*dp2km2
      Z23_=dp4km2 + dp2km2
      Z24_= - p2p3*Z23_
      Z22_=Z22_ + Z24_
      Z22_=p1p4*Z22_
      Z24_=p1k*dp2km2
      Z25_= - p1p2*Z23_
      Z24_=Z24_ + Z25_
      Z24_=p3p4*Z24_
      Z25_=dp2km2*mm2
      Z26_= - p4k*Z25_
      Z21_=Z26_ + Z24_ + Z22_ + Z21_
      Z21_=me2*Z21_
      Z22_=dp4km2*mm2
      Z22_=Z22_ + Z25_
      Z23_=p1p3*Z23_
      Z23_= - 2*Z22_ + Z23_
      Z23_=me4*Z23_
      Z8_=p2k*dp4km2*Z8_
      Z8_=Z8_ + Z21_ + Z23_
      Z21_= - p3k + 2*p2p3
      Z21_=p1p4*Z21_
      Z23_= - p1k + Z10_
      Z23_=p3p4*Z23_
      Z6_=Z6_ - p1p3
      Z6_=me2*Z6_
      Z6_=2*Z6_ + Z23_ + Z21_ + Z19_
      Z6_=dp4km1*Z6_
      Z6_=mm2 + Z6_
      Z6_=dp2km1*Z6_
      Z19_=me2*Z22_
      Z21_= - dp4km1*mm2
      Z6_=Z6_ + Z19_ + Z21_
      Z6_=Z6_*Z11_
      Z2_=Z6_ + Z2_ + 2*Z8_ + Z20_
      Z2_=64*Z2_
      Z6_=Z16_*Z18_
      Z8_=p2p3*dp3km1
      Z11_= - p1p2*dp1km1
      Z8_=Z8_ + Z11_
      Z8_=p1p4*Z8_
      Z6_=4*Z8_ + Z6_ + p1p2
      Z6_=p3p4*Z6_
      Z8_=Z14_*p3k
      Z11_=p1p3*dp1km1
      Z8_= - Z8_ - 8*Z3_ + 1 + 4*Z11_
      Z19_=p1p4*Z8_
      Z20_= - 4*Z12_ + 1 + 8*Z5_
      Z16_= - p1k*Z16_
      Z16_=Z16_ + Z20_
      Z16_=p3p4*Z16_
      Z16_=Z16_ + Z19_ - p1p2 - p1k - p3k - p2p3
      Z16_=me2*Z16_
      Z7_= - Z14_*Z7_
      Z7_= - p2p3 + Z7_
      Z7_=p1p4*Z7_
      Z19_=dp3km1 + dp1km1
      Z21_=p1p4*Z19_
      Z21_= - Z5_ + Z21_
      Z21_=p3p4*Z21_
      Z22_=p1p4*Z3_
      Z21_=Z22_ + Z21_
      Z21_=p2k*Z21_
      Z9_= - p1p42*p2p3*Z9_
      Z22_=p3p42*p1p2*Z13_
      Z6_=2*Z21_ + Z16_ + Z22_ + Z9_ + Z7_ + Z6_
      Z6_=dp4km1*Z6_
      Z7_= - p1p2*Z8_
      Z8_= - p2p3*Z20_
      Z9_= - 1 + Z17_
      Z9_=p1k*Z9_
      Z7_=p3p4 + p1p4 + Z7_ + Z9_ - p3k + Z8_
      Z7_=me2*Z7_
      Z8_= - p3k*dp1km1
      Z9_=p2p3*Z14_
      Z8_=Z8_ + Z9_
      Z8_=Z8_*Z10_
      Z9_= - p2p32*Z13_
      Z8_=Z8_ + Z9_ + p2p3
      Z8_=p1p4*Z8_
      Z9_=p1p22*Z14_
      Z10_=dp3km1*Z18_
      Z9_=Z9_ + Z10_
      Z10_= - p2p3*Z13_
      Z10_= - 1 + Z10_
      Z10_=p1p2*Z10_
      Z9_=2*Z9_ + Z10_
      Z9_=p3p4*Z9_
      Z10_= - p2p3*Z19_
      Z10_=Z3_ + Z10_
      Z10_=p1p2*Z10_
      Z13_= - p2p3*Z5_
      Z10_=Z13_ + Z10_
      Z10_=p4k*Z10_
      Z7_=Z7_ + 2*Z10_ + Z8_ + Z9_
      Z7_=dp2km1*Z7_
      Z8_=p1p4*Z4_
      Z9_=1 + 4*Z5_
      Z10_= - p3p4*Z9_
      Z8_=Z8_ + Z10_
      Z8_=dp4km1*Z8_
      Z9_=p2p3*Z9_
      Z4_= - p1p2*Z4_
      Z4_=Z9_ + Z4_
      Z4_=dp2km1*Z4_
      Z4_=Z8_ + Z4_
      Z4_=p2p4*Z4_
      Z3_=Z5_ + Z3_
      Z5_= - p2k - p4k
      Z5_=Z3_*Z5_
      Z8_= - Z12_ + Z3_ + 1
      Z9_=p2p3*Z8_
      Z3_= - Z11_ + Z3_ - 1
      Z10_= - p1p2*Z3_
      Z11_= - p2p3*Z15_
      Z3_=Z11_ - Z3_
      Z3_=p1p4*Z3_
      Z11_=p1p2*Z15_
      Z8_=Z11_ + Z8_
      Z8_=p3p4*Z8_
      Z3_=Z4_ + Z7_ + Z6_ + Z8_ + Z3_ + Z9_ + Z10_ + Z5_
      Z3_=64*Z3_

      radmuesum= QRmu * t24m2v(1)*Z1_
     .     + QRe * t13m2v(1)*Z2_
     .     + QRmu*QRe * t13m1v(1)*t24m1v(1)*Z3_

      if (ionlybubbles.gt.0) then
         radmuesum= QRmu * t24m2v(0)*(2.d0*dalep24)*Z1_
     .        + QRe * t13m2v(0)*(2.d0*dalep13)*Z2_
     .        + QRmu*QRe * t13m1v(0)*t24m1v(0)*(dalep13+dalep24)*Z3_
      endif
      
cfrpasserafael
c      radmuesum=t24m2v(0)*Z1_*2.d0*dah24 + t13m2v(0)*Z2_*2.d0*dah13 +
c     .     t13m1v(1)*t24m1v(0)*Z3_ * (dah13+dah24)
cfrpasserafael      
      
      extraweights(1)=
     .     QRmu * t24m2v(0)*Z1_ +QRe * t13m2v(0)*Z2_
     .     +QRmu*QRe* t13m1v(0)*t24m1v(0)*Z3_
      extraweights(2)=
     .     QRmu * t24m2v(2)*Z1_+QRe*t13m2v(2)*Z2_
     . + QRmu * QRe * t13m1v(2)*t24m1v(2)*Z3_
      extraweights(1:2) = extraweights(1:2)/radmuesum
      
      endif

      extraweights(0) = 1.d0
      
      radmuesum = coverm2*radmuesum
      
      return
      end
****************************************************
      double precision function apvertex(x)
      implicit double precision (a-h,m,o-z)
      apvertex = (1.d0+x*x)/(1.d0-x)
      return
      end
****************************************************      
      double precision function m2born(s,t,me2,mm2,alpha)
      implicit double precision (a-h,m,o-z)
      pi = 4.d0*atan(1.d0)
      t2 = t*t
      s2 = s*s
      m2born = t2+2.d0*s*t+2.d0*s2+2.d0*mm2*mm2
      m2born = m2born - 4.d0*mm2*s+2.d0*me2*me2 -4.d0*me2*s
      m2born = m2born + 4.d0*me2*mm2
      m2born = 4.d0*pi*alpha*alpha/t2
      return
      end
*************************************      
      subroutine getVP1loop(q2,dal,dah)
      implicit none
      double precision q2,dal,dah
      double complex vpolc,c
      external vpolc
** from invariants.h
      double complex vpdah
      double precision vpdal
      double precision vpvect(5)
      common/vp_dalphas/vpdal,vpdah,vpvect
      data vpvect /5*0.d0/
**      
      integer iarun
      common/ialpharunning/iarun

      dal = 0.d0
      dah = 0.d0
      if (iarun.eq.0) return
      
      c   = vpolc(q2)
      dal = vpdal
      dah = dble(vpdah)
      return
      end
*****************************************************
      subroutine getBsforcts(m2,B1o,dB0o,dB1o)
      implicit none
      double complex B1o,dB0o,dB1o
      double precision m2
      integer*8 i
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      i    = Bget(m2*uquad,m2*uquad,0.d0*uquad)
      B1o  = Bval(bb1,i)
      dB0o = Bval(dbb0,i)
      dB1o = Bval(dbb1,i)
      return
      end
***********************************************
      double complex function Sigrlii(q2,m2)
#ifdef COLLIER
      use collier
#endif
      implicit none
      double precision q2,m2
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi      
#include "invariants.h"
#ifdef COLLIER
      double complex B1o,dB0o,dB1o
      integer rank
      parameter (rank=1)
      double complex b(0:rank/2,0:rank),buv(0:rank/2,0:rank)
      call SwitchOffCacheSystem_cll
      call B_cll(b,buv,q2*cuno,m2*cuno,0.d0*cuno,1)
      Sigrlii = 2.d0*b(0,1)+1.d0
      Sigrlii = -alpha*0.25d0/pi* Sigrlii      
      call SwitchOnCacheSystem_cll
#else
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      ComplexType B(Nbb)
      call Bput(B,q2*uquad,m2*uquad,0.d0*uquad)
c     Sigrlii = 2.d0*B1(q2,m2,0.d0)+1.d0
      Sigrlii = 2.d0*B(bb1)+1.d0
      Sigrlii = -alpha*0.25d0/pi* Sigrlii
#endif
      return
      end
***********************************************
      double complex function dSigrlii(q2,m2)
#ifdef COLLIER
      use collier
#endif
      implicit none
      double precision q2,m2
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi      
#include "invariants.h"
#ifdef COLLIER
      double complex B1o,dB0o,dB1o
      call SwitchOffCacheSystem_cll
      call DB1_cll(dB1o,q2*cuno,m2*cuno,0.d0*cuno)
      dSigrlii = 2.d0*dB1o
      dSigrlii = -alpha*0.25d0/pi* dSigrlii
      call SwitchOnCacheSystem_cll
#else
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      ComplexType B(Nbb)

      call Bput(B,q2*uquad,m2*uquad,0.d0*uquad)
c     dSigrlii = 2.d0*DB1(q2,m2,0.d0)
      dSigrlii = 2.d0*B(dbb1)      
      dSigrlii = -alpha*0.25d0/pi* dSigrlii
#endif      
      return
      end
***********************************************
      double complex function SigSii(q2,m2)
#ifdef COLLIER
      use collier
#endif
      implicit none
      double precision q2,m2
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi      
#include "invariants.h"

#ifdef COLLIER
      double complex B1o,dB0o,dB1o
      integer rank
      parameter (rank=1)
      double complex b(0:rank/2,0:rank),buv(0:rank/2,0:rank)
      call SwitchOffCacheSystem_cll
      call B_cll(b,buv,q2*cuno,m2*cuno,0.d0*cuno,1)
      SigSii = 4.d0*b(0,0)-2.d0
      SigSii = -alpha*0.25d0/pi* SigSii
      call SwitchOnCacheSystem_cll      
#else      
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      ComplexType B(Nbb)
      call Bput(B,q2*uquad,m2*uquad,0.d0*uquad)
c     SigSii = 4.d0*B0(q2,m2,0.d0)-2.d0
      SigSii = 4.d0*B(bb0)-2.d0
      SigSii = -alpha*0.25d0/pi* SigSii
#endif
      return
      end
***********************************************
      double complex function dSigSii(q2,m2)
#ifdef COLLIER
      use collier
#endif
      implicit none
      double precision q2,m2
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi      
#include "invariants.h"

#ifdef COLLIER
      double complex B1o,dB0o,dB1o
      integer rank
      call SwitchOffCacheSystem_cll
      call DB0_cll(dB0o,q2*cuno,m2*cuno,0.d0*cuno)
      dSigSii = 4.d0*dB0o
      dSigSii = -alpha*0.25d0/pi* dSigSii
      call SwitchOnCacheSystem_cll      
#else      
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      ComplexType B(Nbb)
      call Bput(B,q2*uquad,m2*uquad,0.d0*uquad)
c      dSigSii = 4.d0*DB0(q2,m2,0.d0)
      dSigSii = 4.d0*B(dbb0)     
      dSigSii = -alpha*0.25d0/pi* dSigSii
#endif
      return
      end
*****************************************************
      subroutine getQEDcts(cts)
      implicit none
      double precision cts,ctstmp
      double complex tmp
      double complex B1lt,dB0lt,dB1lt
#include "invariants.h"
      integer ifirst
      common/qedcts/ctstmp,ifirst
      data ifirst /0/

      double precision ctstmpmu,ctstmpe      
      common/qedcts_e_mu/ctstmpmu,ctstmpe

      double complex Sigrlii
      double complex dSigrlii
      double complex SigSii
      double complex dSigSii

      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi      
      
      if (ifirst.eq.0) then      
c         i     = Bget(mm2,mm2,0.d0)
c         B1lt  = Bval(bb1,i)
c        dB0lt = Bval(dbb0,i)
c         dB1lt = Bval(dbb1,i)
#ifdef COLLIER
         call getBsforctscollier(mm2,B1lt,dB0lt,dB1lt)
#else
         call getBsforcts(mm2,B1lt,dB0lt,dB1lt)
#endif         
         tmp   = B1lt + 0.5d0 + 2.d0*mm2*dB1lt + 4.d0*mm2*dB0lt

         cts = 2.d0*tmp *pi216m1 * ec2 * QRmu

         ctstmpmu = cts *QRmu
        
c         i     = Bget(me2,me2,0.d0)
c         B1lt  = Bval(bb1,i)
c         dB0lt = Bval(dbb0,i)
c         dB1lt = Bval(dbb1,i)

#ifdef COLLIER
         call getBsforctscollier(me2,B1lt,dB0lt,dB1lt)
#else
         call getBsforcts(me2,B1lt,dB0lt,dB1lt)         
#endif
         tmp   = B1lt + 0.5d0 + 2.d0*me2*dB1lt + 4.d0*me2*dB0lt
         
         cts = cts + 2.d0*tmp *pi216m1 * ec2 * QRe
         ctstmpe = 2.d0*tmp *pi216m1 * ec2 * QRe
         
         ctstmp   = 2.d0 * cts
         ctstmpmu = 2.d0*ctstmpmu
         ctstmpe  = 2.d0*ctstmpe

******the following for the "cross" (Denner A.5)
         dmasse  = 0.5d0*me*(2.d0*Sigrlii(me2,me2)+2.d0*SigSii(me2,me2))   
         dmassmu = 0.5d0*mm*(2.d0*Sigrlii(mm2,mm2)+2.d0*SigSii(mm2,mm2))   
         Crlmuoe2  = 0.d0 ! gia' nel form. C'e' dietro tutto il ragionamento di Fulvio
         Crleoe2   = 0.d0
         CSpmmuoe2 = 1.d0*dmassmu/ec2 ! solo la parte di massa, c'e' dietro tutto il ragionamento di Fulvio
         CSpmeoe2  = 1.d0*dmasse/ec2
***********************************************
         ifirst = 1 !!!!!!!!!!!!!!!!!!!!!!
      endif
      cts = ctstmp*1.d0
      return
      end
*****************
      subroutine loadQEDvertexfuns
      implicit none
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
      common/vertexfuns/
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte      

      data 
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte
     .     /14*(0.d0,0.d0)/
      
#ifdef COLLIER
      call Cfunscollier
      return
#else
      if (QRmu.gt.0) then
         imu = Cget(mm2*uquad,
     .        2.d0*mm2*uquad-2.d0*p1p3*uquad,
     .        mm2*uquad,0.d0*uquad,mm2*uquad,mm2*uquad)
         C0ltmu  = Cval(cc0,imu)
         C1ltmu  = Cval(cc1,imu)
c     C2ltmu  = Cval(cc2,imu)
c     C2ltmu  = C1ltmu
         C00ltmu = Cval(cc00,imu)
         C11ltmu = Cval(cc11,imu)
         C12ltmu = Cval(cc12,imu)
c     C22ltmu = Cval(cc11,imu)
c      C22ltmu = C11ltmu
      endif

      if (QRe.gt.0) then      
         ie  = Cget(me2*uquad,
     .        2.d0*me2*uquad-2.d0*p2p4*uquad,
     .        me2*uquad,0.d0*uquad,me2*uquad,me2*uquad)
         
         C0lte   = Cval(cc0,ie)
         C1lte   = Cval(cc1,ie)
c     C2lte   = Cval(cc2,ie)
c     C2lte   = C1lte      
         C00lte  = Cval(cc00,ie)
         C11lte  = Cval(cc11,ie)
         C12lte  = Cval(cc12,ie)
c     C22lte  = Cval(cc11,ie)
c     C22lte  = C11lte      
      endif
      
c      C2 = C1 and C22 = C11, already assumed in form
#endif
      return
      end
*********************************************************
      subroutine loadQEDfuns1g1L
*** load all functions (vertexes, boxes, self and pentagons)
      implicit none
      double precision k(0:3)
      integer ifirst,iprint
      parameter (iprint = 0)
      integer ponoff(160)
      common/loadltfuns1g1L/ifirst,ponoff
      data ifirst,ponoff /0,160*0/
#include "invariants.h"
#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif

#include "funsdeccmn1g1L.h"
      ComplexType BmuSf1(Nbb),BmuSf2(Nbb)
      ComplexType BeSf1(Nbb),BeSf2(Nbb)
      ComplexType CmuV0(Ncc),CmuV1(Ncc),CmuV2(Ncc),CmuV3(Ncc),
     .     CmuV4(Ncc)
      ComplexType CeV0(Ncc),CeV1(Ncc),CeV2(Ncc),CeV3(Ncc),
     .     CeV4(Ncc)
      ComplexType DmuBx1(Ndd),DmuBx2(Ndd)
      ComplexType DeBx1(Ndd),DeBx2(Ndd)
      ComplexType dDBx1(Ndd),dDBx2(Ndd),dDBx3(Ndd),dDBx4(Ndd)
      ComplexType cDBx1(Ndd),cDBx2(Ndd),cDBx3(Ndd),cDBx4(Ndd)
      ComplexType dPt1(Nee),dPt2(Nee),cPt1(Nee),cPt2(Nee)
      ComplexType Ctest(Ncc)

      if (ifirst.eq.0) then
         include 'funssettozero.f'
         ponoff = 1
         
c         ponoff(1) = 0
c         ponoff(4) = 0
c         ponoff(14) = 0
c         ponoff(37) = 0

c         ponoff(41) = 0
c         ponoff(43) = 0
c         ponoff(51) = 0
c         ponoff(72) = 0

c         ponoff(81) = 0
c         ponoff(84) = 0
c         ponoff(94) = 0
c         ponoff(117) = 0

c         ponoff(121) = 0
c         ponoff(123) = 0
c         ponoff(131) = 0
c         ponoff(151) = 0
         
#ifndef QUAD
         do ifirst = 1,50
            print*,'Are u sure you want to use LT for these'//
     .           ' functions? It''s in quad precision??'
         enddo
#endif
      endif
*** the rest are already in invariants common

      if (QRmu.eq.1) then      
      call Bput(BmuSf1,p1mk2*uquad,0.d0*uquad,mm2*uquad)
      call Bput(BmuSf2,p3pk2*uquad,0.d0*uquad,mm2*uquad)
      B0muSf1 = BmuSf1(bb0)
      B1muSf1 = BmuSf1(bb1)
      B0muSf2 = BmuSf2(bb0)
      B1muSf2 = BmuSf2(bb1)
      call Cput(CmuV0,mm2*uquad,t13*uquad,mm2*uquad,0.d0*uquad,
     .     mm2*uquad,mm2*uquad)
      C0muV0  = CmuV0(cc0)
      C1muV0  = CmuV0(cc1)
      C2muV0  = CmuV0(cc2)
      C00muV0 = CmuV0(cc00)
      C11muV0 = CmuV0(cc11)
      C12muV0 = CmuV0(cc12)
      C22muV0 = CmuV0(cc22)      
      call Cput(CmuV1,mm2*uquad,0.d0*uquad,p1mk2*uquad,0.d0*uquad,
     .     mm2*uquad,mm2*uquad)
      C0muV1  = CmuV1(cc0)
      C1muV1  = CmuV1(cc1)
      C2muV1  = CmuV1(cc2)
      C00muV1 = CmuV1(cc00)
      C11muV1 = CmuV1(cc11)
      C12muV1 = CmuV1(cc12)
      C22muV1 = CmuV1(cc22)
      call Cput(CmuV2,p1mk2*uquad,t24*uquad,mm2*uquad,0.d0*uquad,
     .     mm2*uquad,mm2*uquad)
      C0muV2  = CmuV2(cc0)
      C1muV2  = CmuV2(cc1)
      C2muV2  = CmuV2(cc2)
      C00muV2 = CmuV2(cc00)
      C11muV2 = CmuV2(cc11)
      C12muV2 = CmuV2(cc12)
      C22muV2 = CmuV2(cc22)
      call Cput(CmuV3,p3pk2*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,
     .     mm2*uquad,mm2*uquad)
      C0muV3  = CmuV3(cc0)
      C1muV3  = CmuV3(cc1)
      C2muV3  = CmuV3(cc2)
      C00muV3 = CmuV3(cc00)
      C11muV3 = CmuV3(cc11)
      C12muV3 = CmuV3(cc12)
      C22muV3 = CmuV3(cc22)
      call Cput(CmuV4,mm2*uquad,t24*uquad,p3pk2*uquad,0.d0*uquad,
     .     mm2*uquad,mm2*uquad)
      C0muV4  = CmuV4(cc0)
      C1muV4  = CmuV4(cc1)
      C2muV4  = CmuV4(cc2)
      C00muV4 = CmuV4(cc00)
      C11muV4 = CmuV4(cc11)
      C12muV4 = CmuV4(cc12)
      C22muV4 = CmuV4(cc22)

      call Dput(DmuBx1,mm2*uquad,0.d0*uquad,t24*uquad,mm2*uquad,
     .   p1mk2*uquad,t13*uquad,0.d0*uquad,mm2*uquad,mm2*uquad,mm2*uquad)
      call Dput(DmuBx2,mm2*uquad,t24*uquad,0.d0*uquad,mm2*uquad,
     .   p3pk2*uquad,t13*uquad,0.d0*uquad,mm2*uquad,mm2*uquad,mm2*uquad)
      D0muBx1 = DmuBx1(dd0)
      D1muBx1 = DmuBx1(dd1)
      D2muBx1 = DmuBx1(dd2)
      D3muBx1 = DmuBx1(dd3)
      D00muBx1 = DmuBx1(dd00)
      D11muBx1 = DmuBx1(dd11)
      D12muBx1 = DmuBx1(dd12)
      D13muBx1 = DmuBx1(dd13)      
      D22muBx1 = DmuBx1(dd22)
      D23muBx1 = DmuBx1(dd23)
      D33muBx1 = DmuBx1(dd33)
      D001muBx1 = DmuBx1(dd001)
      D002muBx1 = DmuBx1(dd002)
      D003muBx1 = DmuBx1(dd003)
      D111muBx1 = DmuBx1(dd111)
      D112muBx1 = DmuBx1(dd112)
      D113muBx1 = DmuBx1(dd113)
      D122muBx1 = DmuBx1(dd122)
      D123muBx1 = DmuBx1(dd123)
      D133muBx1 = DmuBx1(dd133)
      D222muBx1 = DmuBx1(dd222)
      D223muBx1 = DmuBx1(dd223)
      D233muBx1 = DmuBx1(dd233)
      D333muBx1 = DmuBx1(dd333)

      D0muBx2 = DmuBx2(dd0)
      D1muBx2 = DmuBx2(dd1)
      D2muBx2 = DmuBx2(dd2)
      D3muBx2 = DmuBx2(dd3)
      D00muBx2 = DmuBx2(dd00)
      D11muBx2 = DmuBx2(dd11)
      D12muBx2 = DmuBx2(dd12)
      D13muBx2 = DmuBx2(dd13)      
      D22muBx2 = DmuBx2(dd22)
      D23muBx2 = DmuBx2(dd23)
      D33muBx2 = DmuBx2(dd33)
      D001muBx2 = DmuBx2(dd001)
      D002muBx2 = DmuBx2(dd002)
      D003muBx2 = DmuBx2(dd003)
      D111muBx2 = DmuBx2(dd111)
      D112muBx2 = DmuBx2(dd112)
      D113muBx2 = DmuBx2(dd113)
      D122muBx2 = DmuBx2(dd122)
      D123muBx2 = DmuBx2(dd123)
      D133muBx2 = DmuBx2(dd133)
      D222muBx2 = DmuBx2(dd222)
      D223muBx2 = DmuBx2(dd223)
      D233muBx2 = DmuBx2(dd233)
      D333muBx2 = DmuBx2(dd333)
      endif

      if (QRe.eq.1) then      
      call Bput(BeSf1,p2mk2*uquad,0.d0*uquad,me2*uquad)
      call Bput(BeSf2,p4pk2*uquad,0.d0*uquad,me2*uquad)
      B0eSf1  = BeSf1(bb0)
      B1eSf1  = Besf1(bb1)
      B0eSf2  = BeSf2(bb0)
      B1eSf2  = BeSf2(bb1)
      
      call Cput(CeV0,me2*uquad,t24*uquad,me2*uquad,0.d0*uquad,
     .     me2*uquad,me2*uquad)
      C0eV0   = CeV0(cc0)
      C1eV0   = CeV0(cc1)
      C2eV0   = CeV0(cc2)
      C00eV0  = CeV0(cc00)
      C11eV0  = CeV0(cc11)
      C12eV0  = CeV0(cc12)
      C22eV0  = CeV0(cc22)
      
      call Cput(CeV1,me2*uquad,0.d0*uquad,p2mk2*uquad,0.d0*uquad,
     .     me2*uquad,me2*uquad)
      C0eV1   = CeV1(cc0)
      C1eV1   = CeV1(cc1)
      C2eV1   = CeV1(cc2)
      C00eV1  = CeV1(cc00)
      C11eV1  = CeV1(cc11)
      C12eV1  = CeV1(cc12)
      C22eV1  = CeV1(cc22)

      call Cput(CeV2,p2mk2*uquad,t13*uquad,me2*uquad,0.d0*uquad,
     .     me2*uquad,me2*uquad)      
      C0eV2   = CeV2(cc0)
      C1eV2   = CeV2(cc1)
      C2eV2   = CeV2(cc2)
      C00eV2  = CeV2(cc00)
      C11eV2  = CeV2(cc11)
      C12eV2  = CeV2(cc12)
      C22eV2  = CeV2(cc22)

      call Cput(CeV3,p4pk2*uquad,0.d0*uquad,me2*uquad,0.d0*uquad,
     .     me2*uquad,me2*uquad)
      C0eV3   = CeV3(cc0)
      C1eV3   = CeV3(cc1)
      C2eV3   = CeV3(cc2)
      C00eV3  = CeV3(cc00)
      C11eV3  = CeV3(cc11)
      C12eV3  = CeV3(cc12)
      C22eV3  = CeV3(cc22)

      call Cput(CeV4,me2*uquad,t13*uquad,p4pk2*uquad,0.d0*uquad,
     .     me2*uquad,me2*uquad)
      C0eV4   = CeV4(cc0)
      C1eV4   = CeV4(cc1)
      C2eV4   = CeV4(cc2)
      C00eV4  = CeV4(cc00)
      C11eV4  = CeV4(cc11)
      C12eV4  = CeV4(cc12)
      C22eV4  = CeV4(cc22)
      
      call Dput(DeBx1,me2*uquad,0.d0*uquad,t13*uquad,me2*uquad,
     .   p2mk2*uquad,t24*uquad,0.d0*uquad,me2*uquad,me2*uquad,me2*uquad)
      call Dput(DeBx2,me2*uquad,t13*uquad,0.d0*uquad,me2*uquad,
     .   p4pk2*uquad,t24*uquad,0.d0*uquad,me2*uquad,me2*uquad,me2*uquad)
      D0eBx1 = DeBx1(dd0)
      D1eBx1 = DeBx1(dd1)
      D2eBx1 = DeBx1(dd2)
      D3eBx1 = DeBx1(dd3)
      D00eBx1 = DeBx1(dd00)
      D11eBx1 = DeBx1(dd11)
      D12eBx1 = DeBx1(dd12)
      D13eBx1 = DeBx1(dd13)      
      D22eBx1 = DeBx1(dd22)
      D23eBx1 = DeBx1(dd23)
      D33eBx1 = DeBx1(dd33)
      D001eBx1 = DeBx1(dd001)
      D002eBx1 = DeBx1(dd002)
      D003eBx1 = DeBx1(dd003)
      D111eBx1 = DeBx1(dd111)
      D112eBx1 = DeBx1(dd112)
      D113eBx1 = DeBx1(dd113)
      D122eBx1 = DeBx1(dd122)
      D123eBx1 = DeBx1(dd123)
      D133eBx1 = DeBx1(dd133)
      D222eBx1 = DeBx1(dd222)
      D223eBx1 = DeBx1(dd223)      
      D233eBx1 = DeBx1(dd233)
      D333eBx1 = DeBx1(dd333)

      D0eBx2 = DeBx2(dd0)
      D1eBx2 = DeBx2(dd1)
      D2eBx2 = DeBx2(dd2)
      D3eBx2 = DeBx2(dd3)
      D00eBx2 = DeBx2(dd00)
      D11eBx2 = DeBx2(dd11)
      D12eBx2 = DeBx2(dd12)
      D13eBx2 = DeBx2(dd13)      
      D22eBx2 = DeBx2(dd22)
      D23eBx2 = DeBx2(dd23)
      D33eBx2 = DeBx2(dd33)
      D001eBx2 = DeBx2(dd001)
      D002eBx2 = DeBx2(dd002)
      D003eBx2 = DeBx2(dd003)
      D111eBx2 = DeBx2(dd111)
      D112eBx2 = DeBx2(dd112)
      D113eBx2 = DeBx2(dd113)
      D122eBx2 = DeBx2(dd122)
      D123eBx2 = DeBx2(dd123)
      D133eBx2 = DeBx2(dd133)
      D222eBx2 = DeBx2(dd222)
      D223eBx2 = DeBx2(dd223)      
      D233eBx2 = DeBx2(dd233)
      D333eBx2 = DeBx2(dd333)
      endif
      
********************
** no up-down boxes or pentagons if one of the charges is 0!!      
      if (QRmu.eq.0.or.QRe.eq.0) return
********************      

      call Dput(dDBx1,p1mk2*uquad,mm2*uquad,me2*uquad,me2*uquad,
     .    t24*uquad,s34*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      call Dput(dDBx2,mm2*uquad,p3pk2*uquad,me2*uquad,me2*uquad,
     .    t24*uquad,s12*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      call Dput(dDBx3,mm2*uquad,mm2*uquad,me2*uquad,p2mk2*uquad,
     .    t13*uquad,s34*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      call Dput(dDBx4,mm2*uquad,mm2*uquad,p4pk2*uquad,me2*uquad,
     .    t13*uquad,s12*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      
      call Dput(cDBx1,p1mk2*uquad,mm2*uquad,me2*uquad,me2*uquad,
     .    t24*uquad,u23*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      call Dput(cDBx2,mm2*uquad,p3pk2*uquad,me2*uquad,me2*uquad,
     .    t24*uquad,u14*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      call Dput(cDBx3,mm2*uquad,mm2*uquad,p2mk2*uquad,me2*uquad,
     .    t13*uquad,u14*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      call Dput(cDBx4,mm2*uquad,mm2*uquad,me2*uquad,p4pk2*uquad,
     .    t13*uquad,u23*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      
      d1D0  = dDBx1(dd0)
      d1D1  = dDBx1(dd1)
      d1D2  = dDBx1(dd2)
      d1D3  = dDBx1(dd3)
      d1D00 = dDBx1(dd00) 
      d1D11 = dDBx1(dd11)
      d1D12 = dDBx1(dd12)
      d1D13 = dDBx1(dd13)
      d1D22 = dDBx1(dd22)
      d1D23 = dDBx1(dd23)
      d1D33 = dDBx1(dd33)

      c1D0  = cDBx1(dd0)
      c1D1  = cDBx1(dd1)
      c1D2  = cDBx1(dd2)
      c1D3  = cDBx1(dd3)
      c1D00 = cDBx1(dd00) 
      c1D11 = cDBx1(dd11)
      c1D12 = cDBx1(dd12)
      c1D13 = cDBx1(dd13)
      c1D22 = cDBx1(dd22)
      c1D23 = cDBx1(dd23)
      c1D33 = cDBx1(dd33)

      d2D0  = dDBx2(dd0)
      d2D1  = dDBx2(dd1)
      d2D2  = dDBx2(dd2)
      d2D3  = dDBx2(dd3)
      d2D00 = dDBx2(dd00)
      d2D11 = dDBx2(dd11)
      d2D12 = dDBx2(dd12)
      d2D13 = dDBx2(dd13)
      d2D22 = dDBx2(dd22)
      d2D23 = dDBx2(dd23)
      d2D33 = dDBx2(dd33)

      c2D0  = cDBx2(dd0)
      c2D1  = cDBx2(dd1)
      c2D2  = cDBx2(dd2)
      c2D3  = cDBx2(dd3)
      c2D00 = cDBx2(dd00) 
      c2D11 = cDBx2(dd11)
      c2D12 = cDBx2(dd12)
      c2D13 = cDBx2(dd13)
      c2D22 = cDBx2(dd22)
      c2D23 = cDBx2(dd23)
      c2D33 = cDBx2(dd33)

      d3D0  = dDBx3(dd0)
      d3D1  = dDBx3(dd1)
      d3D2  = dDBx3(dd2)
      d3D3  = dDBx3(dd3)
      d3D00 = dDBx3(dd00) 
      d3D11 = dDBx3(dd11)
      d3D12 = dDBx3(dd12)
      d3D13 = dDBx3(dd13)
      d3D22 = dDBx3(dd22)
      d3D23 = dDBx3(dd23)
      d3D33 = dDBx3(dd33)

      c3D0  = cDBx3(dd0)
      c3D1  = cDBx3(dd1)
      c3D2  = cDBx3(dd2)
      c3D3  = cDBx3(dd3)
      c3D00 = cDBx3(dd00) 
      c3D11 = cDBx3(dd11)
      c3D12 = cDBx3(dd12)
      c3D13 = cDBx3(dd13)
      c3D22 = cDBx3(dd22)
      c3D23 = cDBx3(dd23)
      c3D33 = cDBx3(dd33)
      
      d4D0  = dDBx4(dd0)
      d4D1  = dDBx4(dd1)
      d4D2  = dDBx4(dd2)
      d4D3  = dDBx4(dd3)
      d4D00 = dDBx4(dd00) 
      d4D11 = dDBx4(dd11)
      d4D12 = dDBx4(dd12)
      d4D13 = dDBx4(dd13)
      d4D22 = dDBx4(dd22)
      d4D23 = dDBx4(dd23)
      d4D33 = dDBx4(dd33)

      c4D0  = cDBx4(dd0)
      c4D1  = cDBx4(dd1)
      c4D2  = cDBx4(dd2)
      c4D3  = cDBx4(dd3)
      c4D00 = cDBx4(dd00) 
      c4D11 = cDBx4(dd11)
      c4D12 = cDBx4(dd12)
      c4D13 = cDBx4(dd13)
      c4D22 = cDBx4(dd22)
      c4D23 = cDBx4(dd23)
      c4D33 = cDBx4(dd33)

c      call Eput(dPt1,mm2*uquad,0.d0*uquad,mm2*uquad,me2*uquad,me2*uquad,
c     .     p1mk2*uquad,p3pk2*uquad,s34*uquad,t24*uquad,s12*uquad,
c     .     0.d0*uquad,mm2*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
c      call Eput(dPt2,mm2*uquad,mm2*uquad,me2*uquad,0.d0*uquad,me2*uquad,
c     .     t13*uquad,s34*uquad,p4pk2*uquad,p2mk2*uquad,s12*uquad,
c     .     0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad,me2*uquad)
c      call Eput(cPt1,mm2*uquad,0.d0*uquad,mm2*uquad,me2*uquad,me2*uquad,
c     .     p1mk2*uquad,p3pk2*uquad,u23*uquad,t24*uquad,u14*uquad,
c     .     0.d0*uquad,mm2*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
c      call Eput(cPt2,mm2*uquad,mm2*uquad,me2*uquad,0.d0*uquad,me2*uquad,
c     .     t13*uquad,u23*uquad,p2mk2*uquad,p4pk2*uquad,u14*uquad,
c     .     0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad,me2*uquad)
*** NEW ROUTING!
      call Eput(dPt1,0.d0*uquad,mm2*uquad,me2*uquad,me2*uquad,mm2*uquad,
     .     p3pk2*uquad,s34*uquad,t24*uquad,s12*uquad,p1mk2*uquad,
     .     mm2*uquad,mm2*uquad,0.d0*uquad,me2*uquad,0.d0*uquad)
      call Eput(dPt2,me2*uquad,mm2*uquad,mm2*uquad,me2*uquad,0.d0*uquad,
     .     s12*uquad,t13*uquad,s34*uquad,p4pk2*uquad,p2mk2*uquad,
     .     me2*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      call Eput(cPt1,0.d0*uquad,mm2*uquad,me2*uquad,me2*uquad,mm2*uquad,
     .     p3pk2*uquad,u23*uquad,t24*uquad,u14*uquad,p1mk2*uquad,
     .     mm2*uquad,mm2*uquad,0.d0*uquad,me2*uquad,0.d0*uquad)
      call Eput(cPt2,me2*uquad,mm2*uquad,mm2*uquad,me2*uquad,0.d0*uquad,
     .     u23*uquad,t13*uquad,u14*uquad,p4pk2*uquad,p2mk2*uquad,
     .     me2*uquad,0.d0*uquad,mm2*uquad,0.d0*uquad,me2*uquad)
      
 333  d1E0   = dPt1(ee0)    !* ponoff(1)
      d1E1   = dPt1(ee1)    !* ponoff(2)
      d1E2   = dPt1(ee2)    !* ponoff(3)
      d1E3   = dPt1(ee3)    !* ponoff(4)
      d1E4   = dPt1(ee4)    !* ponoff(5)
      d1E00  = dPt1(ee00)   !* ponoff(6)
      d1E11  = dPt1(ee11)   !* ponoff(7)
      d1E12  = dPt1(ee12)   !* ponoff(8)
      d1E13  = dPt1(ee13)   !* ponoff(9)
      d1E14  = dPt1(ee14)   !* ponoff(10)
      d1E22  = dPt1(ee22)   !* ponoff(11)
      d1E23  = dPt1(ee23)   !* ponoff(12)
      d1E24  = dPt1(ee24)   !* ponoff(13)
      d1E33  = dPt1(ee33)   !* ponoff(14)
      d1E34  = dPt1(ee34)   !* ponoff(15)
      d1E44  = dPt1(ee44)   !* ponoff(16)
      d1E001 = dPt1(ee001)  !* ponoff(17)
      d1E002 = dPt1(ee002)  !* ponoff(18)
      d1E003 = dPt1(ee003)  !* ponoff(19)
      d1E004 = dPt1(ee004)  !* ponoff(20)
      d1E111 = dPt1(ee111)  !* ponoff(21)
      d1E112 = dPt1(ee112)  !* ponoff(22)
      d1E113 = dPt1(ee113)  !* ponoff(23)
      d1E114 = dPt1(ee114)  !* ponoff(24)
      d1E122 = dPt1(ee122)  !* ponoff(25)
      d1E123 = dPt1(ee123)  !* ponoff(26)
      d1E124 = dPt1(ee124)  !* ponoff(27)
      d1E133 = dPt1(ee133)  !* ponoff(28)
      d1E134 = dPt1(ee134)  !* ponoff(29)
      d1E144 = dPt1(ee144)  !* ponoff(30)
      d1E222 = dPt1(ee222)  !* ponoff(31)
      d1E223 = dPt1(ee223)  !* ponoff(32)
      d1E224 = dPt1(ee224)  !* ponoff(33)
      d1E233 = dPt1(ee233)  !* ponoff(34)
      d1E234 = dPt1(ee234)  !* ponoff(35)
      d1E244 = dPt1(ee244)  !* ponoff(36)
      d1E333 = dPt1(ee333)  !* ponoff(37)
      d1E334 = dPt1(ee334)  !* ponoff(38)
      d1E344 = dPt1(ee344)  !* ponoff(39)
      d1E444 = dPt1(ee444)  !* ponoff(40)

      d2E0   = dPt2(ee0)        !* ponoff(41)
      d2E1   = dPt2(ee1)    !* ponoff(42)
      d2E2   = dPt2(ee2)    !* ponoff(43)
      d2E3   = dPt2(ee3)    !* ponoff(44)
      d2E4   = dPt2(ee4)    !* ponoff(45)
      d2E00  = dPt2(ee00)   !* ponoff(46)
      d2E11  = dPt2(ee11)   !* ponoff(47)
      d2E12  = dPt2(ee12)   !* ponoff(48)
      d2E13  = dPt2(ee13)   !* ponoff(49)
      d2E14  = dPt2(ee14)   !* ponoff(50)
      d2E22  = dPt2(ee22)   !* ponoff(51)
      d2E23  = dPt2(ee23)   !* ponoff(52)
      d2E24  = dPt2(ee24)   !* ponoff(53)
      d2E33  = dPt2(ee33)   !* ponoff(54)
      d2E34  = dPt2(ee34)   !* ponoff(55)
      d2E44  = dPt2(ee44)   !* ponoff(56)
      d2E001 = dPt2(ee001)  !* ponoff(57)
      d2E002 = dPt2(ee002)  !* ponoff(58)
      d2E003 = dPt2(ee003)  !* ponoff(59)
      d2E004 = dPt2(ee004)  !* ponoff(60)    
      d2E111 = dPt2(ee111)  !* ponoff(61)
      d2E112 = dPt2(ee112)  !* ponoff(62)
      d2E113 = dPt2(ee113)  !* ponoff(63)
      d2E114 = dPt2(ee114)  !* ponoff(64)
      d2E122 = dPt2(ee122)  !* ponoff(65)
      d2E123 = dPt2(ee123)  !* ponoff(66)
      d2E124 = dPt2(ee124)  !* ponoff(67)
      d2E133 = dPt2(ee133)  !* ponoff(68)
      d2E134 = dPt2(ee134)  !* ponoff(69)
      d2E144 = dPt2(ee144)  !* ponoff(70)
      d2E222 = dPt2(ee222)  !* ponoff(71)
      d2E223 = dPt2(ee223)  !* ponoff(72)
      d2E224 = dPt2(ee224)  !* ponoff(73)
      d2E233 = dPt2(ee233)  !* ponoff(74)
      d2E234 = dPt2(ee234)  !* ponoff(75)
      d2E244 = dPt2(ee244)  !* ponoff(76)
      d2E333 = dPt2(ee333)  !* ponoff(77)
      d2E334 = dPt2(ee334)  !* ponoff(78)
      d2E344 = dPt2(ee344)  !* ponoff(79)
      d2E444 = dPt2(ee444)  !* ponoff(80)

      c1E0   = cPt1(ee0)    !* ponoff(81)
      c1E1   = cPt1(ee1)    !* ponoff(82)
      c1E2   = cPt1(ee2)    !* ponoff(83)
      c1E3   = cPt1(ee3)    !* ponoff(84)
      c1E4   = cPt1(ee4)    !* ponoff(85)
      c1E00  = cPt1(ee00)   !* ponoff(86)
      c1E11  = cPt1(ee11)   !* ponoff(87)
      c1E12  = cPt1(ee12)   !* ponoff(88)
      c1E13  = cPt1(ee13)   !* ponoff(89)
      c1E14  = cPt1(ee14)   !* ponoff(90)
      c1E22  = cPt1(ee22)   !* ponoff(91)
      c1E23  = cPt1(ee23)   !* ponoff(92)
      c1E24  = cPt1(ee24)   !* ponoff(93)
      c1E33  = cPt1(ee33)   !* ponoff(94)
      c1E34  = cPt1(ee34)   !* ponoff(95)
      c1E44  = cPt1(ee44)   !* ponoff(96)
      c1E001 = cPt1(ee001)  !* ponoff(97)
      c1E002 = cPt1(ee002)  !* ponoff(98)
      c1E003 = cPt1(ee003)  !* ponoff(99)
      c1E004 = cPt1(ee004)  !* ponoff(100)
      c1E111 = cPt1(ee111)  !* ponoff(101)
      c1E112 = cPt1(ee112)  !* ponoff(102)
      c1E113 = cPt1(ee113)  !* ponoff(103)
      c1E114 = cPt1(ee114)  !* ponoff(104)
      c1E122 = cPt1(ee122)  !* ponoff(105)
      c1E123 = cPt1(ee123)  !* ponoff(106)
      c1E124 = cPt1(ee124)  !* ponoff(107)
      c1E133 = cPt1(ee133)  !* ponoff(108)
      c1E134 = cPt1(ee134)  !* ponoff(109)
      c1E144 = cPt1(ee144)  !* ponoff(110)
      c1E222 = cPt1(ee222)  !* ponoff(111)
      c1E223 = cPt1(ee223)  !* ponoff(112)
      c1E224 = cPt1(ee224)  !* ponoff(113)
      c1E233 = cPt1(ee233)  !* ponoff(114)
      c1E234 = cPt1(ee234)  !* ponoff(115)
      c1E244 = cPt1(ee244)  !* ponoff(116)
      c1E333 = cPt1(ee333)  !* ponoff(117)
      c1E334 = cPt1(ee334)  !* ponoff(118)
      c1E344 = cPt1(ee344)  !* ponoff(119)
      c1E444 = cPt1(ee444)  !* ponoff(120)

      c2E0   = cPt2(ee0)    !* ponoff(121)
      c2E1   = cPt2(ee1)    !* ponoff(122)
      c2E2   = cPt2(ee2)    !* ponoff(123)
      c2E3   = cPt2(ee3)    !* ponoff(124)
      c2E4   = cPt2(ee4)    !* ponoff(125)
      c2E00  = cPt2(ee00)   !* ponoff(126)
      c2E11  = cPt2(ee11)   !* ponoff(127)
      c2E12  = cPt2(ee12)   !* ponoff(128)
      c2E13  = cPt2(ee13)   !* ponoff(129)
      c2E14  = cPt2(ee14)   !* ponoff(130)
      c2E22  = cPt2(ee22)   !* ponoff(131)
      c2E23  = cPt2(ee23)   !* ponoff(132)
      c2E24  = cPt2(ee24)   !* ponoff(133)
      c2E33  = cPt2(ee33)   !* ponoff(134)
      c2E34  = cPt2(ee34)   !* ponoff(135)
      c2E44  = cPt2(ee44)   !* ponoff(136)
      c2E001 = cPt2(ee001)  !* ponoff(137)
      c2E002 = cPt2(ee002)  !* ponoff(138)
      c2E003 = cPt2(ee003)  !* ponoff(139)
      c2E004 = cPt2(ee004)  !* ponoff(140)
      c2E111 = cPt2(ee111)  !* ponoff(141)
      c2E112 = cPt2(ee112)  !* ponoff(142)
      c2E113 = cPt2(ee113)  !* ponoff(143)
      c2E114 = cPt2(ee114)  !* ponoff(144)
      c2E122 = cPt2(ee122)  !* ponoff(145)
      c2E123 = cPt2(ee123)  !* ponoff(146)
      c2E124 = cPt2(ee124)  !* ponoff(147)
      c2E133 = cPt2(ee133)  !* ponoff(148)
      c2E134 = cPt2(ee134)  !* ponoff(149)
      c2E144 = cPt2(ee144)  !* ponoff(150)
      c2E222 = cPt2(ee222)  !* ponoff(151)
      c2E223 = cPt2(ee223)  !* ponoff(152)
      c2E224 = cPt2(ee224)  !* ponoff(153)
      c2E233 = cPt2(ee233)  !* ponoff(154)
      c2E234 = cPt2(ee234)  !* ponoff(155)
      c2E244 = cPt2(ee244)  !* ponoff(156)
      c2E333 = cPt2(ee333)  !* ponoff(157)
      c2E334 = cPt2(ee334)  !* ponoff(158)
      c2E344 = cPt2(ee344)  !* ponoff(159)
      c2E444 = cPt2(ee444)  !* ponoff(160)
c      if (iprint.eq.1) then
c         include 'printltfun'
c      endif
      return
      end
**********************************************************************
c      subroutine loadQEDfuns1g1LCOLLIERten(iprint,iten,p1,p2,p3,p4,k)
      subroutine loadQEDfuns1g1LCOLLIER
*** load all functions (vertexes, boxes, self and pentagons)
#ifdef COLLIER
      use collier
#endif      
      implicit none
      double precision k(0:3)
c      double precision p1k,p2k,p3k,p4k
#include "invariants.h"
#include "funsdeccmn1g1L.h"
      integer iprint,iten
      parameter (iprint = 0)
      real*16 qdot
      external qdot
#ifdef COLLIER
      integer r1,r2,r3,iii
      parameter (r1 = 1)
      parameter (r2 = 2)
      parameter (r3 = 3)
      double complex CmuV0(0:r2/2,0:r2,0:r2)
      double complex CmuuvV0(0:r2/2,0:r2,0:r2)
      double complex CmuV1(0:r2/2,0:r2,0:r2)
      double complex CmuuvV1(0:r2/2,0:r2,0:r2)
      double complex CmuV2(0:r2/2,0:r2,0:r2)
      double complex CmuuvV2(0:r2/2,0:r2,0:r2)
      double complex CmuV3(0:r2/2,0:r2,0:r2)
      double complex CmuuvV3(0:r2/2,0:r2,0:r2)
      double complex CmuV4(0:r2/2,0:r2,0:r2)
      double complex CmuuvV4(0:r2/2,0:r2,0:r2)
      double complex CeV0(0:r2/2,0:r2,0:r2)
      double complex CeuvV0(0:r2/2,0:r2,0:r2)
      double complex CeV1(0:r2/2,0:r2,0:r2)
      double complex CeuvV1(0:r2/2,0:r2,0:r2)
      double complex CeV2(0:r2/2,0:r2,0:r2)
      double complex CeuvV2(0:r2/2,0:r2,0:r2)
      double complex CeV3(0:r2/2,0:r2,0:r2)
      double complex CeuvV3(0:r2/2,0:r2,0:r2)
      double complex CeV4(0:r2/2,0:r2,0:r2)
      double complex CeuvV4(0:r2/2,0:r2,0:r2)

      double complex DmuBx1(0:r3/2,0:r3,0:r3,0:r3)
      double complex DmuBx2(0:r3/2,0:r3,0:r3,0:r3)
      double complex DeBx1(0:r3/2,0:r3,0:r3,0:r3)
      double complex DeBx2(0:r3/2,0:r3,0:r3,0:r3)
      double complex DmuuvBx1(0:r3/2,0:r3,0:r3,0:r3)
      double complex DmuuvBx2(0:r3/2,0:r3,0:r3,0:r3)
      double complex DeuvBx1(0:r3/2,0:r3,0:r3,0:r3)
      double complex DeuvBx2(0:r3/2,0:r3,0:r3,0:r3)

      double complex dDBx1(0:r2/2,0:r2,0:r2,0:r2)
      double complex dDBx2(0:r2/2,0:r2,0:r2,0:r2)
      double complex dDBx3(0:r2/2,0:r2,0:r2,0:r2)
      double complex dDBx4(0:r2/2,0:r2,0:r2,0:r2)
      double complex dDuvBx1(0:r2/2,0:r2,0:r2,0:r2)
      double complex dDuvBx2(0:r2/2,0:r2,0:r2,0:r2)
      double complex dDuvBx3(0:r2/2,0:r2,0:r2,0:r2)
      double complex dDuvBx4(0:r2/2,0:r2,0:r2,0:r2)
      double complex cDBx1(0:r2/2,0:r2,0:r2,0:r2)
      double complex cDBx2(0:r2/2,0:r2,0:r2,0:r2)
      double complex cDBx3(0:r2/2,0:r2,0:r2,0:r2)
      double complex cDBx4(0:r2/2,0:r2,0:r2,0:r2)
      double complex cDuvBx1(0:r2/2,0:r2,0:r2,0:r2)
      double complex cDuvBx2(0:r2/2,0:r2,0:r2,0:r2)
      double complex cDuvBx3(0:r2/2,0:r2,0:r2,0:r2)
      double complex cDuvBx4(0:r2/2,0:r2,0:r2,0:r2)

      double complex BmuSf1(0:r1/2,0:r1)
      double complex BmuSf2(0:r1/2,0:r1)
      double complex BmuuvSf1(0:r1/2,0:r1)
      double complex BmuuvSf2(0:r1/2,0:r1)
      double complex BeSf1(0:r1/2,0:r1)
      double complex BeSf2(0:r1/2,0:r1)
      double complex BeuvSf1(0:r1/2,0:r1)
      double complex BeuvSf2(0:r1/2,0:r1)

      double complex dEPt1(0:r3/2,0:r3,0:r3,0:r3,0:r3)
      double complex dEPt2(0:r3/2,0:r3,0:r3,0:r3,0:r3)
      double complex cEPt1(0:r3/2,0:r3,0:r3,0:r3,0:r3)
      double complex cEPt2(0:r3/2,0:r3,0:r3,0:r3,0:r3)
      double complex dEuvPt1(0:r3/2,0:r3,0:r3,0:r3,0:r3)
      double complex dEuvPt2(0:r3/2,0:r3,0:r3,0:r3,0:r3)
      double complex cEuvPt1(0:r3/2,0:r3,0:r3,0:r3,0:r3)
      double complex cEuvPt2(0:r3/2,0:r3,0:r3,0:r3,0:r3)

      double precision PtErr(0:r3)
      double precision acc0,acc1,acc2
      integer ierr
      
      integer indexes(0:3,4),i,l,j,m,ii,mask(35*4)
      double complex c,tmp
      integer ifirst
      integer ponoff(160)
      common/funs1g1LCOLLIERtenpt/c,ifirst,indexes,mask,ponoff,iten
      data ifirst,ponoff,iten /0,160*0,0/
      data mask /
     .     0,0,0,0, ! 1
     .     1,0,0,0, ! 2
     .     0,1,0,0, ! 3
     .     0,0,1,0, ! 4
     .     0,0,0,1, ! 5
     .     2,0,0,0, ! 6
     .     1,1,0,0, ! 7
     .     1,0,1,0, ! 8
     .     1,0,0,1, ! 9
     .     0,2,0,0, !10
     .     0,1,1,0, !11
     .     0,1,0,1, !12
     .     0,0,2,0, !13
     .     0,0,1,1, !14
     .     0,0,0,2, !15
     .     3,0,0,0, !16
     .     2,1,0,0, !17
     .     2,0,1,0, !18
     .     2,0,0,1, !19
     .     1,2,0,0, !20
     .     1,1,1,0, !21
     .     1,1,0,1, !22
     .     1,0,2,0, !23
     .     1,0,1,1, !24
     .     1,0,0,2, !25
     .     0,3,0,0, !26
     .     0,2,1,0, !27
     .     0,2,0,1, !28
     .     0,1,2,0, !29
     .     0,1,1,1, !30
     .     0,1,0,2, !31
     .     0,0,3,0, !32
     .     0,0,2,1, !33
     .     0,0,1,2, !34
     .     0,0,0,3  !35
     .     / 
      if (ifirst.eq.0) then
         indexes = 0
         indexes(0,1) = 1
         indexes(1,2) = 1
         indexes(2,3) = 1
         indexes(3,4) = 1

         ponoff = 1
c         ponoff(1) = 0
c         ponoff(4) = 0
c         ponoff(14) = 0
c         ponoff(37) = 0

c         ponoff(41) = 0
c         ponoff(43) = 0
c         ponoff(51) = 0
c         ponoff(72) = 0

c         ponoff(81) = 0
c         ponoff(84) = 0
c         ponoff(94) = 0
c         ponoff(117) = 0

c         ponoff(121) = 0
c         ponoff(123) = 0
c         ponoff(131) = 0
c         ponoff(151) = 0
                 
         include 'funssettozero.f'
         c = cuno
         ifirst = 1
      endif

c      if (iten.eq.1) call InitCacheSystem_cll(3,5)      
      call InitEvent_cll(2)

      
      if (QRmu.eq.1) then
      call B_cll(BmuSf1,BmuuvSf1,c*p1mk2,c*0.d0,c*mm2,r1)
      call B_cll(BmuSf2,BmuuvSf2,c*p3pk2,c*0.d0,c*mm2,r1)
      B0muSf1 = BmuSf1(0,0)
      B1muSf1 = BmuSf1(0,1)
      B0muSf2 = BmuSf2(0,0)
      B1muSf2 = BmuSf2(0,1)
      
      call C_cll(CmuV0,CmuuvV0,c*mm2,c*t13,c*mm2,c*0.d0,
     .     c*mm2,c*mm2,r2)
      C0muV0  = CmuV0(0,0,0)
      C1muV0  = CmuV0(0,1,0)
      C2muV0  = CmuV0(0,0,1)
      C00muV0 = CmuV0(1,0,0)
      C11muV0 = CmuV0(0,2,0)
      C12muV0 = CmuV0(0,1,1)
      C22muV0 = CmuV0(0,0,2)
      
      call C_Cll(CmuV1,CmuuvV1,c*mm2,c*0.d0,c*p1mk2,c*0.d0,c*mm2,c*mm2,
     .     r2)
      C0muV1  = CmuV1(0,0,0)
      C1muV1  = CmuV1(0,1,0)
      C2muV1  = CmuV1(0,0,1)
      C00muV1 = CmuV1(1,0,0)
      C11muV1 = CmuV1(0,2,0)
      C12muV1 = CmuV1(0,1,1)
      C22muV1 = CmuV1(0,0,2)
      
      call C_cll(CmuV2,CmuuvV2,c*p1mk2,c*t24,c*mm2,c*0.d0,c*mm2,c*mm2,
     .     r2)      
      C0muV2  = CmuV2(0,0,0)
      C1muV2  = CmuV2(0,1,0)
      C2muV2  = CmuV2(0,0,1)
      C00muV2 = CmuV2(1,0,0)
      C11muV2 = CmuV2(0,2,0)
      C12muV2 = CmuV2(0,1,1)
      C22muV2 = CmuV2(0,0,2)

      call C_cll(CmuV3,CmuuvV3,c*p3pk2,c*0.d0,c*mm2,c*0.d0,c*mm2,c*mm2,
     .     r2)
      C0muV3  = CmuV3(0,0,0)
      C1muV3  = CmuV3(0,1,0)
      C2muV3  = CmuV3(0,0,1)
      C00muV3 = CmuV3(1,0,0)
      C11muV3 = CmuV3(0,2,0)
      C12muV3 = CmuV3(0,1,1)
      C22muV3 = CmuV3(0,0,2)
      
      call C_cll(CmuV4,CmuuvV4,c*mm2,c*t24,c*p3pk2,c*0.d0,c*mm2,c*mm2,
     .     r2)
      C0muV4  = CmuV4(0,0,0)
      C1muV4  = CmuV4(0,1,0)
      C2muV4  = CmuV4(0,0,1)
      C00muV4 = CmuV4(1,0,0)
      C11muV4 = CmuV4(0,2,0)
      C12muV4 = CmuV4(0,1,1)
      C22muV4 = CmuV4(0,0,2)
      
      call D_cll(DmuBx1,DmuuvBx1,
     .     c*mm2,c*0.d0,c*t24,c*mm2,c*p1mk2,c*t13,
     .     c*0.d0,c*mm2,c*mm2,c*mm2,r3)
      call D_cll(DmuBx2,DmuuvBx2,
     .     c*mm2,c*t24,c*0.d0,c*mm2,c*p3pk2,c*t13,
     .     c*0.d0,c*mm2,c*mm2,c*mm2,r3)            
      D0muBx1   = DmuBx1(0,0,0,0)
      D1muBx1   = DmuBx1(0,1,0,0)
      D2muBx1   = DmuBx1(0,0,1,0)
      D3muBx1   = DmuBx1(0,0,0,1)
      D00muBx1  = DmuBx1(1,0,0,0)
      D11muBx1  = DmuBx1(0,2,0,0)
      D12muBx1  = DmuBx1(0,1,1,0)
      D13muBx1  = DmuBx1(0,1,0,1)
      D22muBx1  = DmuBx1(0,0,2,0)
      D23muBx1  = DmuBx1(0,0,1,1)
      D33muBx1  = DmuBx1(0,0,0,2)
      D001muBx1 = DmuBx1(1,1,0,0)
      D002muBx1 = DmuBx1(1,0,1,0)
      D003muBx1 = DmuBx1(1,0,0,1)
      D111muBx1 = DmuBx1(0,3,0,0)
      D112muBx1 = DmuBx1(0,2,1,0)
      D113muBx1 = DmuBx1(0,2,0,1)
      D122muBx1 = DmuBx1(0,1,2,0)
      D123muBx1 = DmuBx1(0,1,1,1)
      D133muBx1 = DmuBx1(0,1,0,2)
      D222muBx1 = DmuBx1(0,0,3,0)
      D223muBx1 = DmuBx1(0,0,2,1)
      D233muBx1 = DmuBx1(0,0,1,2)
      D333muBx1 = DmuBx1(0,0,0,3)
      D0muBx2   = DmuBx2(0,0,0,0)
      D1muBx2   = DmuBx2(0,1,0,0)
      D2muBx2   = DmuBx2(0,0,1,0)
      D3muBx2   = DmuBx2(0,0,0,1)
      D00muBx2  = DmuBx2(1,0,0,0)
      D11muBx2  = DmuBx2(0,2,0,0)
      D12muBx2  = DmuBx2(0,1,1,0)
      D13muBx2  = DmuBx2(0,1,0,1)
      D22muBx2  = DmuBx2(0,0,2,0)
      D23muBx2  = DmuBx2(0,0,1,1)
      D33muBx2  = DmuBx2(0,0,0,2)
      D001muBx2 = DmuBx2(1,1,0,0)
      D002muBx2 = DmuBx2(1,0,1,0)
      D003muBx2 = DmuBx2(1,0,0,1)
      D111muBx2 = DmuBx2(0,3,0,0)
      D112muBx2 = DmuBx2(0,2,1,0)
      D113muBx2 = DmuBx2(0,2,0,1)
      D122muBx2 = DmuBx2(0,1,2,0)
      D123muBx2 = DmuBx2(0,1,1,1)
      D133muBx2 = DmuBx2(0,1,0,2)
      D222muBx2 = DmuBx2(0,0,3,0)
      D223muBx2 = DmuBx2(0,0,2,1)
      D233muBx2 = DmuBx2(0,0,1,2)
      D333muBx2 = DmuBx2(0,0,0,3)
      endif

      if (QRe.eq.1) then
      call B_cll(BeSf1,BeuvSf1,c*p2mk2,c*0.d0,c*me2,r1)
      call B_cll(BeSf2,BeuvSf2,c*p4pk2,c*0.d0,c*me2,r1)
      B0eSf1 = BeSf1(0,0)
      B1eSf1 = BeSf1(0,1)
      B0eSf2 = BeSf2(0,0)
      B1eSf2 = BeSf2(0,1)

      call C_cll(CeV0,CeuvV0,c*me2,c*t24,c*me2,c*0.d0,
     .     c*me2,c*me2,r2)
      C0eV0  = CeV0(0,0,0)
      C1eV0  = CeV0(0,1,0)
      C2eV0  = CeV0(0,0,1)
      C00eV0 = CeV0(1,0,0)
      C11eV0 = CeV0(0,2,0)
      C12eV0 = CeV0(0,1,1)
      C22eV0 = CeV0(0,0,2)

      call C_Cll(CeV1,CeuvV1,c*me2,c*0.d0,c*p2mk2,c*0.d0,c*me2,c*me2,
     .     r2)
      C0eV1  = CeV1(0,0,0)
      C1eV1  = CeV1(0,1,0)
      C2eV1  = CeV1(0,0,1)
      C00eV1 = CeV1(1,0,0)
      C11eV1 = CeV1(0,2,0)
      C12eV1 = CeV1(0,1,1)
      C22eV1 = CeV1(0,0,2)

      call C_cll(CeV2,CeuvV2,c*p2mk2,c*t13,c*me2,c*0.d0,c*me2,c*me2,
     .     r2)
      C0eV2  = CeV2(0,0,0)
      C1eV2  = CeV2(0,1,0)
      C2eV2  = CeV2(0,0,1)
      C00eV2 = CeV2(1,0,0)
      C11eV2 = CeV2(0,2,0)
      C12eV2 = CeV2(0,1,1)
      C22eV2 = CeV2(0,0,2)

      call C_cll(CeV3,CeuvV3,c*p4pk2,c*0.d0,c*me2,c*0.d0,c*me2,c*me2,
     .     r2)      
      C0eV3  = CeV3(0,0,0)
      C1eV3  = CeV3(0,1,0)
      C2eV3  = CeV3(0,0,1)
      C00eV3 = CeV3(1,0,0)
      C11eV3 = CeV3(0,2,0)
      C12eV3 = CeV3(0,1,1)
      C22eV3 = CeV3(0,0,2)
      
      call C_cll(CeV4,CeuvV4,c*me2,c*t13,c*p4pk2,c*0.d0,c*me2,c*me2,
     .     r2)      
      C0eV4  = CeV4(0,0,0)
      C1eV4  = CeV4(0,1,0)
      C2eV4  = CeV4(0,0,1)
      C00eV4 = CeV4(1,0,0)
      C11eV4 = CeV4(0,2,0)
      C12eV4 = CeV4(0,1,1)
      C22eV4 = CeV4(0,0,2)
      
      call D_cll(DeBx1,DeuvBx1,
     .     c*me2,c*0.d0,c*t13,c*me2,c*p2mk2,c*t24,
     .     c*0.d0,c*me2,c*me2,c*me2,r3)
      call D_cll(DeBx2,DeuvBx2,
     .     c*me2,c*t13,c*0.d0,c*me2,c*p4pk2,c*t24,
     .     c*0.d0,c*me2,c*me2,c*me2,r3)      
      D0eBx1   = DeBx1(0,0,0,0)
      D1eBx1   = DeBx1(0,1,0,0)
      D2eBx1   = DeBx1(0,0,1,0)
      D3eBx1   = DeBx1(0,0,0,1)
      D00eBx1  = DeBx1(1,0,0,0)
      D11eBx1  = DeBx1(0,2,0,0)
      D12eBx1  = DeBx1(0,1,1,0)
      D13eBx1  = DeBx1(0,1,0,1)
      D22eBx1  = DeBx1(0,0,2,0)
      D23eBx1  = DeBx1(0,0,1,1)
      D33eBx1  = DeBx1(0,0,0,2)
      D001eBx1 = DeBx1(1,1,0,0)
      D002eBx1 = DeBx1(1,0,1,0)
      D003eBx1 = DeBx1(1,0,0,1)
      D111eBx1 = DeBx1(0,3,0,0)
      D112eBx1 = DeBx1(0,2,1,0)
      D113eBx1 = DeBx1(0,2,0,1)
      D122eBx1 = DeBx1(0,1,2,0)
      D123eBx1 = DeBx1(0,1,1,1)
      D133eBx1 = DeBx1(0,1,0,2)
      D222eBx1 = DeBx1(0,0,3,0)
      D223eBx1 = DeBx1(0,0,2,1)
      D233eBx1 = DeBx1(0,0,1,2)
      D333eBx1 = DeBx1(0,0,0,3)
      D0eBx2   = DeBx2(0,0,0,0)
      D1eBx2   = DeBx2(0,1,0,0)
      D2eBx2   = DeBx2(0,0,1,0)
      D3eBx2   = DeBx2(0,0,0,1)
      D00eBx2  = DeBx2(1,0,0,0)
      D11eBx2  = DeBx2(0,2,0,0)
      D12eBx2  = DeBx2(0,1,1,0)
      D13eBx2  = DeBx2(0,1,0,1)
      D22eBx2  = DeBx2(0,0,2,0)
      D23eBx2  = DeBx2(0,0,1,1)
      D33eBx2  = DeBx2(0,0,0,2)
      D001eBx2 = DeBx2(1,1,0,0)
      D002eBx2 = DeBx2(1,0,1,0)
      D003eBx2 = DeBx2(1,0,0,1)
      D111eBx2 = DeBx2(0,3,0,0)
      D112eBx2 = DeBx2(0,2,1,0)
      D113eBx2 = DeBx2(0,2,0,1)
      D122eBx2 = DeBx2(0,1,2,0)
      D123eBx2 = DeBx2(0,1,1,1)
      D133eBx2 = DeBx2(0,1,0,2)
      D222eBx2 = DeBx2(0,0,3,0)
      D223eBx2 = DeBx2(0,0,2,1)
      D233eBx2 = DeBx2(0,0,1,2)
      D333eBx2 = DeBx2(0,0,0,3)
      endif
      
**********************************
      if (QRmu.eq.0.or.QRe.eq.0) return
**********************************      
      IF (iten.eq.0) THEN

      call D_cll(dDBx1,dDuvBx1,c*p1mk2,c*mm2,c*me2,c*me2,c*t24,c*s34,
     .     c*0.d0,c*mm2,c*0.d0,c*me2,r2)
      call D_cll(dDBx2,dDuvBx2,c*mm2,c*p3pk2,c*me2,c*me2,c*t24,c*s12,
     .     c*0.d0,c*mm2,c*0.d0,c*me2,r2)
      call D_cll(dDBx3,dDuvBx3,c*mm2,c*mm2,c*me2,c*p2mk2,c*t13,c*s34,
     .     c*0.d0,c*mm2,c*0.d0,c*me2,r2)
      call D_cll(dDBx4,dDuvBx4,c*mm2,c*mm2,c*p4pk2,c*me2,c*t13,c*s12,
     .     c*0.d0,c*mm2,c*0.d0,c*me2,r2)      

      call D_cll(cDBx1,cDuvBx1,c*p1mk2,c*mm2,c*me2,c*me2,c*t24,c*u23,
     .     c*0.d0,c*mm2,c*0.d0,c*me2,r2)
      call D_cll(cDBx2,cDuvBx2,c*mm2,c*p3pk2,c*me2,c*me2,c*t24,c*u14,
     .     c*0.d0,c*mm2,c*0.d0,c*me2,r2)
      call D_cll(cDBx3,cDuvBx3,c*mm2,c*mm2,c*p2mk2,c*me2,c*t13,c*u14,
     .     c*0.d0,c*mm2,c*0.d0,c*me2,r2)
      call D_cll(cDBx4,cDuvBx4,c*mm2,c*mm2,c*me2,c*p4pk2,c*t13,c*u23,
     .     c*0.d0,c*mm2,c*0.d0,c*me2,r2)

      d1D0  = dDBx1(0,0,0,0)
      d1D1  = dDBx1(0,1,0,0)
      d1D2  = dDBx1(0,0,1,0)
      d1D3  = dDBx1(0,0,0,1)
      d1D00 = dDBx1(1,0,0,0) 
      d1D11 = dDBx1(0,2,0,0)
      d1D12 = dDBx1(0,1,1,0)
      d1D13 = dDBx1(0,1,0,1)
      d1D22 = dDBx1(0,0,2,0)
      d1D23 = dDBx1(0,0,1,1)
      d1D33 = dDBx1(0,0,0,2)

      d2D0  = dDBx2(0,0,0,0)
      d2D1  = dDBx2(0,1,0,0)
      d2D2  = dDBx2(0,0,1,0)
      d2D3  = dDBx2(0,0,0,1)
      d2D00 = dDBx2(1,0,0,0) 
      d2D11 = dDBx2(0,2,0,0)
      d2D12 = dDBx2(0,1,1,0)
      d2D13 = dDBx2(0,1,0,1)
      d2D22 = dDBx2(0,0,2,0)
      d2D23 = dDBx2(0,0,1,1)
      d2D33 = dDBx2(0,0,0,2)

      d3D0  = dDBx3(0,0,0,0)
      d3D1  = dDBx3(0,1,0,0)
      d3D2  = dDBx3(0,0,1,0)
      d3D3  = dDBx3(0,0,0,1)
      d3D00 = dDBx3(1,0,0,0) 
      d3D11 = dDBx3(0,2,0,0)
      d3D12 = dDBx3(0,1,1,0)
      d3D13 = dDBx3(0,1,0,1)
      d3D22 = dDBx3(0,0,2,0)
      d3D23 = dDBx3(0,0,1,1)
      d3D33 = dDBx3(0,0,0,2)

      d4D0  = dDBx4(0,0,0,0)
      d4D1  = dDBx4(0,1,0,0)
      d4D2  = dDBx4(0,0,1,0)
      d4D3  = dDBx4(0,0,0,1)
      d4D00 = dDBx4(1,0,0,0) 
      d4D11 = dDBx4(0,2,0,0)
      d4D12 = dDBx4(0,1,1,0)
      d4D13 = dDBx4(0,1,0,1)
      d4D22 = dDBx4(0,0,2,0)
      d4D23 = dDBx4(0,0,1,1)
      d4D33 = dDBx4(0,0,0,2)
      
      c1D0  = cDBx1(0,0,0,0)
      c1D1  = cDBx1(0,1,0,0)
      c1D2  = cDBx1(0,0,1,0)
      c1D3  = cDBx1(0,0,0,1)
      c1D00 = cDBx1(1,0,0,0) 
      c1D11 = cDBx1(0,2,0,0)
      c1D12 = cDBx1(0,1,1,0)
      c1D13 = cDBx1(0,1,0,1)
      c1D22 = cDBx1(0,0,2,0)
      c1D23 = cDBx1(0,0,1,1)
      c1D33 = cDBx1(0,0,0,2)

      c2D0  = cDBx2(0,0,0,0)
      c2D1  = cDBx2(0,1,0,0)
      c2D2  = cDBx2(0,0,1,0)
      c2D3  = cDBx2(0,0,0,1)
      c2D00 = cDBx2(1,0,0,0) 
      c2D11 = cDBx2(0,2,0,0)
      c2D12 = cDBx2(0,1,1,0)
      c2D13 = cDBx2(0,1,0,1)
      c2D22 = cDBx2(0,0,2,0)
      c2D23 = cDBx2(0,0,1,1)
      c2D33 = cDBx2(0,0,0,2)

      c3D0  = cDBx3(0,0,0,0)
      c3D1  = cDBx3(0,1,0,0)
      c3D2  = cDBx3(0,0,1,0)
      c3D3  = cDBx3(0,0,0,1)
      c3D00 = cDBx3(1,0,0,0) 
      c3D11 = cDBx3(0,2,0,0)
      c3D12 = cDBx3(0,1,1,0)
      c3D13 = cDBx3(0,1,0,1)
      c3D22 = cDBx3(0,0,2,0)
      c3D23 = cDBx3(0,0,1,1)
      c3D33 = cDBx3(0,0,0,2)

      c4D0  = cDBx4(0,0,0,0)
      c4D1  = cDBx4(0,1,0,0)
      c4D2  = cDBx4(0,0,1,0)
      c4D3  = cDBx4(0,0,0,1)
      c4D00 = cDBx4(1,0,0,0) 
      c4D11 = cDBx4(0,2,0,0)
      c4D12 = cDBx4(0,1,1,0)
      c4D13 = cDBx4(0,1,0,1)
      c4D22 = cDBx4(0,0,2,0)
      c4D23 = cDBx4(0,0,1,1)
      c4D33 = cDBx4(0,0,0,2)

* OLD ROUTING      
c      call E_cll(dEPt1,dEuvPt1,c*mm2,c*0.d0,c*mm2,c*me2,c*me2,
c     .     c*p1mk2,c*p3pk2,c*s34,c*t24,c*s12,c*0.d0,c*mm2,c*mm2,
c     .     c*0.d0,c*me2,r3,PtErr)
c      call E_cll(dEPt2,dEuvPt2,c*mm2,c*mm2,c*me2,c*0.d0,c*me2,c*t13,
c     .     c*s34,c*p4pk2,c*p2mk2,c*s12,c*0.d0,c*mm2,c*0.d0,c*me2,
c     .     c*me2,r3,PtErr)
c      call E_cll(cEPt1,cEuvPt1,c*mm2,c*0.d0,c*mm2,c*me2,c*me2,c*p1mk2,
c     .     c*p3pk2,c*u23,c*t24,c*u14,c*0.d0,c*mm2,c*mm2,c*0.d0,c*me2,
c     .     r3,PtErr)
c      call E_cll(cEPt2,cEuvPt2,c*mm2,c*mm2,c*me2,c*0.d0,c*me2,c*t13,
c     .     c*u23,c*p2mk2,c*p4pk2,c*u14,c*0.d0,c*mm2,c*0.d0,c*me2,
c     .     c*me2,r3,PtErr)
*** NEW ROUTING!
      call E_cll(dEPt1,dEuvPt1,0.d0*c,mm2*c,me2*c,me2*c,mm2*c,
     .     p3pk2*c,s34*c,t24*c,s12*c,p1mk2*c,
     .     mm2*c,mm2*c,0.d0*c,me2*c,0.d0*c,r3,PtErr)
      call E_cll(dEPt2,dEuvPt2,me2*c,mm2*c,mm2*c,me2*c,0.d0*c,
     .     s12*c,t13*c,s34*c,p4pk2*c,p2mk2*c,
     .     me2*c,0.d0*c,mm2*c,0.d0*c,me2*c,r3,PtErr)
      call E_cll(cEPt1,cEuvPt1,0.d0*c,mm2*c,me2*c,me2*c,mm2*c,
     .     p3pk2*c,u23*c,t24*c,u14*c,p1mk2*c,
     .     mm2*c,mm2*c,0.d0*c,me2*c,0.d0*c,
     .     r3,PtErr)
      call E_cll(cEPt2,cEuvPt2,me2*c,mm2*c,mm2*c,me2*c,0.d0*c,
     .     u23*c,t13*c,u14*c,p4pk2*c,p2mk2*c,
     .     me2*c,0.d0*c,mm2*c,0.d0*c,me2*c,r3,PtErr)

      d1E0   = dEPt1(0,0,0,0,0) !* ponoff(1)  
      d1E1   = dEPt1(0,1,0,0,0) !* ponoff(2)  
      d1E2   = dEPt1(0,0,1,0,0) !* ponoff(3)  
      d1E3   = dEPt1(0,0,0,1,0) !* ponoff(4)  
      d1E4   = dEPt1(0,0,0,0,1) !* ponoff(5)  
      d1E00  = dEPt1(1,0,0,0,0) !* ponoff(6)  
      d1E11  = dEPt1(0,2,0,0,0) !* ponoff(7)  
      d1E12  = dEPt1(0,1,1,0,0) !* ponoff(8)  
      d1E13  = dEPt1(0,1,0,1,0) !* ponoff(9)  
      d1E14  = dEPt1(0,1,0,0,1) !* ponoff(10) 
      d1E22  = dEPt1(0,0,2,0,0) !* ponoff(11) 
      d1E23  = dEPt1(0,0,1,1,0) !* ponoff(12) 
      d1E24  = dEPt1(0,0,1,0,1) !* ponoff(13) 
      d1E33  = dEPt1(0,0,0,2,0) !* ponoff(14) 
      d1E34  = dEPt1(0,0,0,1,1) !* ponoff(15) 
      d1E44  = dEPt1(0,0,0,0,2) !* ponoff(16) 
      d1E001 = dEPt1(1,1,0,0,0) !* ponoff(17) 
      d1E002 = dEPt1(1,0,1,0,0) !* ponoff(18) 
      d1E003 = dEPt1(1,0,0,1,0) !* ponoff(19) 
      d1E004 = dEPt1(1,0,0,0,1) !* ponoff(20) 
      d1E111 = dEPt1(0,3,0,0,0) !* ponoff(21) 
      d1E112 = dEPt1(0,2,1,0,0) !* ponoff(22) 
      d1E113 = dEPt1(0,2,0,1,0) !* ponoff(23) 
      d1E114 = dEPt1(0,2,0,0,1) !* ponoff(24) 
      d1E122 = dEPt1(0,1,2,0,0) !* ponoff(25) 
      d1E123 = dEPt1(0,1,1,1,0) !* ponoff(26) 
      d1E124 = dEPt1(0,1,1,0,1) !* ponoff(27) 
      d1E133 = dEPt1(0,1,0,2,0) !* ponoff(28) 
      d1E134 = dEPt1(0,1,0,1,1) !* ponoff(29) 
      d1E144 = dEPt1(0,1,0,0,2) !* ponoff(30) 
      d1E222 = dEPt1(0,0,3,0,0) !* ponoff(31) 
      d1E223 = dEPt1(0,0,2,1,0) !* ponoff(32) 
      d1E224 = dEPt1(0,0,2,0,1) !* ponoff(33) 
      d1E233 = dEPt1(0,0,1,2,0) !* ponoff(34) 
      d1E234 = dEPt1(0,0,1,1,1) !* ponoff(35) 
      d1E244 = dEPt1(0,0,1,0,2) !* ponoff(36) 
      d1E333 = dEPt1(0,0,0,3,0) !* ponoff(37) 
      d1E334 = dEPt1(0,0,0,2,1) !* ponoff(38) 
      d1E344 = dEPt1(0,0,0,1,2) !* ponoff(39) 
      d1E444 = dEPt1(0,0,0,0,3) !* ponoff(40) 
                                             
      d2E0   = dEPt2(0,0,0,0,0) !* ponoff(41)
      d2E1   = dEPt2(0,1,0,0,0) !* ponoff(42) 
      d2E2   = dEPt2(0,0,1,0,0) !* ponoff(43) 
      d2E3   = dEPt2(0,0,0,1,0) !* ponoff(44) 
      d2E4   = dEPt2(0,0,0,0,1) !* ponoff(45) 
      d2E00  = dEPt2(1,0,0,0,0) !* ponoff(46) 
      d2E11  = dEPt2(0,2,0,0,0) !* ponoff(47) 
      d2E12  = dEPt2(0,1,1,0,0) !* ponoff(48) 
      d2E13  = dEPt2(0,1,0,1,0) !* ponoff(49) 
      d2E14  = dEPt2(0,1,0,0,1) !* ponoff(50) 
      d2E22  = dEPt2(0,0,2,0,0) !* ponoff(51) 
      d2E23  = dEPt2(0,0,1,1,0) !* ponoff(52) 
      d2E24  = dEPt2(0,0,1,0,1) !* ponoff(53) 
      d2E33  = dEPt2(0,0,0,2,0) !* ponoff(54) 
      d2E34  = dEPt2(0,0,0,1,1) !* ponoff(55) 
      d2E44  = dEPt2(0,0,0,0,2) !* ponoff(56) 
      d2E001 = dEPt2(1,1,0,0,0) !* ponoff(57) 
      d2E002 = dEPt2(1,0,1,0,0) !* ponoff(58) 
      d2E003 = dEPt2(1,0,0,1,0) !* ponoff(59) 
      d2E004 = dEPt2(1,0,0,0,1) !* ponoff(60) 
      d2E111 = dEPt2(0,3,0,0,0) !* ponoff(61) 
      d2E112 = dEPt2(0,2,1,0,0) !* ponoff(62) 
      d2E113 = dEPt2(0,2,0,1,0) !* ponoff(63) 
      d2E114 = dEPt2(0,2,0,0,1) !* ponoff(64) 
      d2E122 = dEPt2(0,1,2,0,0) !* ponoff(65) 
      d2E123 = dEPt2(0,1,1,1,0) !* ponoff(66) 
      d2E124 = dEPt2(0,1,1,0,1) !* ponoff(67) 
      d2E133 = dEPt2(0,1,0,2,0) !* ponoff(68) 
      d2E134 = dEPt2(0,1,0,1,1) !* ponoff(69) 
      d2E144 = dEPt2(0,1,0,0,2) !* ponoff(70) 
      d2E222 = dEPt2(0,0,3,0,0) !* ponoff(71) 
      d2E223 = dEPt2(0,0,2,1,0) !* ponoff(72) 
      d2E224 = dEPt2(0,0,2,0,1) !* ponoff(73) 
      d2E233 = dEPt2(0,0,1,2,0) !* ponoff(74) 
      d2E234 = dEPt2(0,0,1,1,1) !* ponoff(75) 
      d2E244 = dEPt2(0,0,1,0,2) !* ponoff(76) 
      d2E333 = dEPt2(0,0,0,3,0) !* ponoff(77) 
      d2E334 = dEPt2(0,0,0,2,1) !* ponoff(78) 
      d2E344 = dEPt2(0,0,0,1,2) !* ponoff(79) 
      d2E444 = dEPt2(0,0,0,0,3) !* ponoff(80) 
                                             
      c1E0   = cEPt1(0,0,0,0,0) !* ponoff(81) 
      c1E1   = cEPt1(0,1,0,0,0) !* ponoff(82) 
      c1E2   = cEPt1(0,0,1,0,0) !* ponoff(83) 
      c1E3   = cEPt1(0,0,0,1,0) !* ponoff(84) 
      c1E4   = cEPt1(0,0,0,0,1) !* ponoff(85) 
      c1E00  = cEPt1(1,0,0,0,0) !* ponoff(86) 
      c1E11  = cEPt1(0,2,0,0,0) !* ponoff(87) 
      c1E12  = cEPt1(0,1,1,0,0) !* ponoff(88) 
      c1E13  = cEPt1(0,1,0,1,0) !* ponoff(89) 
      c1E14  = cEPt1(0,1,0,0,1) !* ponoff(90) 
      c1E22  = cEPt1(0,0,2,0,0) !* ponoff(91) 
      c1E23  = cEPt1(0,0,1,1,0) !* ponoff(92) 
      c1E24  = cEPt1(0,0,1,0,1) !* ponoff(93) 
      c1E33  = cEPt1(0,0,0,2,0) !* ponoff(94) 
      c1E34  = cEPt1(0,0,0,1,1) !* ponoff(95) 
      c1E44  = cEPt1(0,0,0,0,2) !* ponoff(96) 
      c1E001 = cEPt1(1,1,0,0,0) !* ponoff(97) 
      c1E002 = cEPt1(1,0,1,0,0) !* ponoff(98) 
      c1E003 = cEPt1(1,0,0,1,0) !* ponoff(99) 
      c1E004 = cEPt1(1,0,0,0,1) !* ponoff(100)
      c1E111 = cEPt1(0,3,0,0,0) !* ponoff(101)
      c1E112 = cEPt1(0,2,1,0,0) !* ponoff(102)
      c1E113 = cEPt1(0,2,0,1,0) !* ponoff(103)
      c1E114 = cEPt1(0,2,0,0,1) !* ponoff(104)
      c1E122 = cEPt1(0,1,2,0,0) !* ponoff(105)
      c1E123 = cEPt1(0,1,1,1,0) !* ponoff(106)
      c1E124 = cEPt1(0,1,1,0,1) !* ponoff(107)
      c1E133 = cEPt1(0,1,0,2,0) !* ponoff(108)
      c1E134 = cEPt1(0,1,0,1,1) !* ponoff(109)
      c1E144 = cEPt1(0,1,0,0,2) !* ponoff(110)
      c1E222 = cEPt1(0,0,3,0,0) !* ponoff(111)
      c1E223 = cEPt1(0,0,2,1,0) !* ponoff(112)
      c1E224 = cEPt1(0,0,2,0,1) !* ponoff(113)
      c1E233 = cEPt1(0,0,1,2,0) !* ponoff(114)
      c1E234 = cEPt1(0,0,1,1,1) !* ponoff(115)
      c1E244 = cEPt1(0,0,1,0,2) !* ponoff(116)
      c1E333 = cEPt1(0,0,0,3,0) !* ponoff(117)
      c1E334 = cEPt1(0,0,0,2,1) !* ponoff(118)
      c1E344 = cEPt1(0,0,0,1,2) !* ponoff(119)
      c1E444 = cEPt1(0,0,0,0,3) !* ponoff(120)
                                             
      c2E0   = cEPt2(0,0,0,0,0) !* ponoff(121)
      c2E1   = cEPt2(0,1,0,0,0) !* ponoff(122)
      c2E2   = cEPt2(0,0,1,0,0) !* ponoff(123)
      c2E3   = cEPt2(0,0,0,1,0) !* ponoff(124)
      c2E4   = cEPt2(0,0,0,0,1) !* ponoff(125)
      c2E00  = cEPt2(1,0,0,0,0) !* ponoff(126)
      c2E11  = cEPt2(0,2,0,0,0) !* ponoff(127)
      c2E12  = cEPt2(0,1,1,0,0) !* ponoff(128)
      c2E13  = cEPt2(0,1,0,1,0) !* ponoff(129)
      c2E14  = cEPt2(0,1,0,0,1) !* ponoff(130)
      c2E22  = cEPt2(0,0,2,0,0) !* ponoff(131)
      c2E23  = cEPt2(0,0,1,1,0) !* ponoff(132)
      c2E24  = cEPt2(0,0,1,0,1) !* ponoff(133)
      c2E33  = cEPt2(0,0,0,2,0) !* ponoff(134)
      c2E34  = cEPt2(0,0,0,1,1) !* ponoff(135)
      c2E44  = cEPt2(0,0,0,0,2) !* ponoff(136)
      c2E001 = cEPt2(1,1,0,0,0) !* ponoff(137)
      c2E002 = cEPt2(1,0,1,0,0) !* ponoff(138)
      c2E003 = cEPt2(1,0,0,1,0) !* ponoff(139)
      c2E004 = cEPt2(1,0,0,0,1) !* ponoff(140)
      c2E111 = cEPt2(0,3,0,0,0) !* ponoff(141)
      c2E112 = cEPt2(0,2,1,0,0) !* ponoff(142)
      c2E113 = cEPt2(0,2,0,1,0) !* ponoff(143)
      c2E114 = cEPt2(0,2,0,0,1) !* ponoff(144)
      c2E122 = cEPt2(0,1,2,0,0) !* ponoff(145)
      c2E123 = cEPt2(0,1,1,1,0) !* ponoff(146)
      c2E124 = cEPt2(0,1,1,0,1) !* ponoff(147)
      c2E133 = cEPt2(0,1,0,2,0) !* ponoff(148)
      c2E134 = cEPt2(0,1,0,1,1) !* ponoff(149)
      c2E144 = cEPt2(0,1,0,0,2) !* ponoff(150)
      c2E222 = cEPt2(0,0,3,0,0) !* ponoff(151)
      c2E223 = cEPt2(0,0,2,1,0) !* ponoff(152)
      c2E224 = cEPt2(0,0,2,0,1) !* ponoff(153)
      c2E233 = cEPt2(0,0,1,2,0) !* ponoff(154)
      c2E234 = cEPt2(0,0,1,1,1) !* ponoff(155)
      c2E244 = cEPt2(0,0,1,0,2) !* ponoff(156)
      c2E333 = cEPt2(0,0,0,3,0) !* ponoff(157)
      c2E334 = cEPt2(0,0,0,2,1) !* ponoff(158)
      c2E344 = cEPt2(0,0,0,1,2) !* ponoff(159)
      c2E444 = cEPt2(0,0,0,0,3) !* ponoff(160)
      
      ELSE ! if iten = 1

         print*,'the new routing is not implemented here! stopping'
         print*,'TENSORS NOT IMPLEMENTED ANYMORE!'
         print*,'in case you need to pass in input also momenta'
         stop

c$$$         call InitCacheSystem_cll(3,5) ! without this, it doesnt work...
c$$$         call InitEvent_cll(3)
c$$$
c$$$c     call InitCacheSystem_cll(3,5)
c$$$c     call InitEvent_cll(3)
c$$$         
c$$$         call Dten_cll(dDBx1ten,dDBx1tenuv,c*(p1-k),c*(p4-p2),-c*p2,
c$$$     .        c*p1mk2,c*mm2,c*me2,c*me2,c*t24,c*s34,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,r2)
c$$$         call Dten_cll(dDBx2ten,dDBx2tenuv,c*p1,c*(p4-p2),-c*p2,
c$$$     .        c*mm2,c*p3pk2,c*me2,c*me2,c*t24,
c$$$     .        c*s12,c*0.d0,c*mm2,c*0.d0,c*me2,r2)
c$$$         call Dten_cll(dDBx3ten,dDBx3tenuv,c*p1,c*(p1-p3),c*(k-p2),
c$$$     .        c*mm2,c*mm2,c*me2,c*p2mk2,c*t13,c*s34,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,r2)
c$$$         call Dten_cll(dDBx4ten,dDBx4tenuv,c*p1,c*(p1-p3),-c*p2,
c$$$     .        c*mm2,c*mm2,c*p4pk2,c*me2,c*t13,c*s12,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,r2)      
c$$$
c$$$         call Dten_cll(cDBx1ten,cDBx1tenuv,c*(p1-k),c*(p4-p2),c*p4,
c$$$     .        c*p1mk2,c*mm2,c*me2,c*me2,c*t24,c*u23,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,r2)
c$$$         call Dten_cll(cDBx2ten,cDBx2tenuv,c*p1,c*(p4-p2),c*p4,
c$$$     .        c*mm2,c*p3pk2,c*me2,c*me2,c*t24,c*u14,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,r2)
c$$$         call Dten_cll(cDBx3ten,cDBx3tenuv,c*p1,c*(p1-p3),c*p4,
c$$$     .        c*mm2,c*mm2,c*p2mk2,c*me2,c*t13,c*u14,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,r2)
c$$$         call Dten_cll(cDBx4ten,cDBx4tenuv,c*p1,c*(p1-p3),c*(p4+k),
c$$$     .        c*mm2,c*mm2,c*me2,c*p4pk2,c*t13,c*u23,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,r2)
c$$$         
c$$$****  standard routing         
c$$$         call Eten_cll(dEPt1ten,dEPt1tenuv,
c$$$     .        c*p1, c*(p1-k), c*(p4-p2),-c*p2,
c$$$     .        c*mm2,c*0.d0,c*mm2,c*me2,c*me2,
c$$$     .        c*p1mk2,c*p3pk2,c*s34,c*t24,c*s12,
c$$$     .        c*0.d0,c*mm2,c*mm2,c*0.d0,c*me2,r3)
c$$$         call Eten_cll(dEPt2ten,dEPt2tenuv,
c$$$     .        c*p1,c*(p1-p3),c*(k-p2),-c*p2,
c$$$     .        c*mm2,c*mm2,c*me2,c*0.d0,c*me2,
c$$$     .        c*t13,c*s34,c*p4pk2,c*p2mk2,c*s12,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,c*me2,r3)
c$$$         call Eten_cll(cEPt1ten,cEPt1tenuv,
c$$$     .        c*p1,c*(p1-k),c*(p4-p2),c*p4,
c$$$     .        c*mm2,c*0.d0,c*mm2,c*me2,c*me2,
c$$$     .        c*p1mk2,c*p3pk2,c*u23,c*t24,c*u14,
c$$$     .        c*0.d0,c*mm2,c*mm2,c*0.d0,c*me2,r3)
c$$$         call Eten_cll(cEPt2ten,cEPt2tenuv,
c$$$     .        c*p1,c*(p1-p3),c*(k+p4),c*p4,
c$$$     .        c*mm2,c*mm2,c*me2,c*0.d0,c*me2,
c$$$     .        c*t13,c*u23,c*p2mk2,c*p4pk2,c*u14,
c$$$     .        c*0.d0,c*mm2,c*0.d0,c*me2,c*me2,r3)
c$$$         
c$$$c     dEPt1ten = (0.d0,0.d0)
c$$$c     dEPt2ten = (0.d0,0.d0)
c$$$c     cEPt1ten = (0.d0,0.d0)
c$$$c     cEPt2ten = (0.d0,0.d0)
c$$$***** second routing, doesn't change much... uff
c$$$***** trying recola routing         
c$$$c     call Eten_cll(dEPt1ten,dEPt1tenuv,
c$$$c     .        c*p1,c*(p1+p2),c*(p3+k),c*k,
c$$$c     .        c*mm2,c*me2,c*me2,c*mm2,c*0.d0,
c$$$c     .        c*s12,c*t24,c*s34,c*p3pk2,c*p1mk2,
c$$$c     .        c*mm2,c*0.d0,c*me2,c*0.d0,c*mm2,r3)
c$$$c     call Eten_cll(dEPt2ten,dEPt2tenuv,
c$$$c     .        c*p1,c*(p1+p2),c*(p3+p4),c*p3,
c$$$c     .        c*mm2,c*me2,c*0.d0,c*me2,c*mm2,
c$$$c     .        c*s12,c*p2mk2,c*p4pk2,c*s34,c*t13,
c$$$c     .        c*mm2,c*0.d0,c*me2,c*me2,c*0.d0,r3)
c$$$c     call Eten_cll(cEPt1ten,cEPt1tenuv,
c$$$c     .        c*p1,c*(p1-p4),c*(p3+k),c*k,
c$$$c     .        c*mm2,c*me2,c*me2,c*mm2,c*0.d0,
c$$$c     .        c*u14,c*t24,c*u23,c*p3pk2,c*p1mk2,
c$$$c     .        c*mm2,c*0.d0,c*me2,c*0.d0,c*mm2,r3)
c$$$c     call Eten_cll(cEPt2ten,cEPt2tenuv,
c$$$c     .        c*p1,c*(p1-p4),c*(p3-p2),c*p3,
c$$$c     .        c*mm2,c*me2,c*0.d0,c*me2,c*mm2,
c$$$c     .        c*u14,c*p4pk2,c*p2mk2,c*u23,c*t13,
c$$$c     .        c*mm2,c*0.d0,c*me2,c*me2,c*0.d0,r3)
c$$$***********************
      ENDIF

c      if (iprint.eq.1) then 
c         include 'printltfun'
c      endif
#endif
      return
      end
*********************************************************
      double complex function dotten(ten,nact,rank,pa,pb,pc)
      implicit none
      integer i,j,k,nact,rank
      double precision pa(0:3),pb(0:3),pc(0:3)
      double precision pam(0:3),pbm(0:3),pcm(0:3)
      double complex ten(0:rank,0:rank,0:rank,0:rank),tmp
      integer ifirst,indexes(0:3,4)
      double precision st
      data ifirst /0/
      double precision u,uo2,uo3,uo6
      common/dottencmn/u,uo2,uo3,uo6,ifirst,indexes

      if (ifirst.eq.0) then
         indexes      = 0
         indexes(0,1) = 1
         indexes(1,2) = 1
         indexes(2,3) = 1
         indexes(3,4) = 1
         u   = 1.d0
         uo2 = u/2.d0
         uo3 = u/3.d0
         uo6 = u/6.d0
         ifirst = 1
      endif

      dotten = (0.d0,0.d0)
      if (nact.eq.-1) then
         dotten = ten(2,0,0,0)-ten(0,2,0,0)-ten(0,0,2,0)-ten(0,0,0,2)
         return
      endif

      if (nact.eq.0) then
       dotten =
     .+(ten(3,0,0,0) - ten(1,2,0,0) - ten(1,0,2,0) - ten(1,0,0,2))*pa(0)
     .-(ten(2,1,0,0) - ten(0,3,0,0) - ten(0,1,2,0) - ten(0,1,0,2))*pa(1)
     .-(ten(2,0,1,0) - ten(0,2,1,0) - ten(0,0,3,0) - ten(0,0,1,2))*pa(2)
     .-(ten(2,0,0,1) - ten(0,2,0,1) - ten(0,0,2,1) - ten(0,0,0,3))*pa(3)         
       return
      endif

      pam(0)   =  pa(0)
      pam(1:3) = -pa(1:3)      
      if (nact.eq.1) then
         do i = 0,3
            dotten = dotten + pam(i)*ten(
     .           indexes(i,1),
     .           indexes(i,2),
     .           indexes(i,3),
     .           indexes(i,4))
         enddo
         return
      endif

      pbm(0)   =  pb(0)
      pbm(1:3) = -pb(1:3)
      if (nact.eq.2) then
         do i = 0,3
            do j = i,3
               if (j.eq.i) then
                  st = uo2
               else
                  st = u
               endif               
               dotten = dotten + st*(pam(i)*pbm(j)+pam(j)*pbm(i))
     .              *ten(
     .              indexes(i,1)+indexes(j,1),
     .              indexes(i,2)+indexes(j,2),
     .              indexes(i,3)+indexes(j,3),
     .              indexes(i,4)+indexes(j,4)
     .              )
            enddo
         enddo
         return
      endif

      pcm(0)   =  pc(0)
      pcm(1:3) = -pc(1:3)
      if (nact.eq.3) then
         do i = 0,3
            do j = i,3
               do k = j,3
                  if (i.eq.j.or.k.eq.j)  then
                     st = uo2
                     if (j.eq.i.and.k.eq.j) st = uo6
                  else
                     st = u
                  endif
                  dotten = dotten +
     .                 st*(
     .                 pam(i)*pbm(j)*pcm(k)+pam(i)*pbm(k)*pcm(j)+
     .                 pam(j)*pbm(i)*pcm(k)+pam(j)*pbm(k)*pcm(i)+
     .                 pam(k)*pbm(i)*pcm(j)+pam(k)*pbm(j)*pcm(i)
     .                 )
     .                 *ten(
     .                 indexes(i,1)+indexes(j,1)+indexes(k,1),
     .                 indexes(i,2)+indexes(j,2)+indexes(k,2),
     .                 indexes(i,3)+indexes(j,3)+indexes(k,3),
     .                 indexes(i,4)+indexes(j,4)+indexes(k,4)
     .                 )
               enddo
            enddo
         enddo
      endif

c$$$      if (nact.eq.3) then
c$$$         do i = 0,3
c$$$            do j = 0,3
c$$$               do k = 0,3
c$$$                  st = 1.d0
c$$$                  dotten = dotten +
c$$$     .                 st*(
c$$$     .                 pam(i)*pbm(j)*pcm(k)
c$$$     .                 )
c$$$     .                 *ten(
c$$$     .                 indexes(i,1)+indexes(j,1)+indexes(k,1),
c$$$     .                 indexes(i,2)+indexes(j,2)+indexes(k,2),
c$$$     .                 indexes(i,3)+indexes(j,3)+indexes(k,3),
c$$$     .                 indexes(i,4)+indexes(j,4)+indexes(k,4)
c$$$     .                 )
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      endif      

      return
      end
*********************************************************
      subroutine getQEDvertexes(vertexes)
      implicit none
      double complex C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu
      double complex C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte
      double complex C00ltmu,C00lte,vertexes,cover,amplvertexes
      double complex amplvertexesmu,amplvertexese,vertexmu,vertexe
      integer ifirst
      data ifirst /0/
      common/ifirstvert/cover,ifirst
#include "invariants.h"
      common/vertexfuns/
     .     C0ltmu,C1ltmu,C2ltmu,C11ltmu,C12ltmu,C22ltmu,
     .     C0lte ,C1lte ,C2lte ,C11lte ,C12lte ,C22lte,
     .     C00ltmu,C00lte      

      double complex amplve,amplvmu
      common/vertexeslegs/amplve,amplvmu
      
      if (ifirst.eq.0) then
         cover = ec4*im*pi216m1*4.d0
         ifirst = 1
      endif

      if (Qmu.lt.0.d0) then
c         include 'oneloop/vertexesminus.f'
      amplvmu =
     .  + QRmu * (  - 1.D0/2.D0*cmuce + cmuce*C00ltmu + cmuce*p1p3*
     .    C12ltmu + 2*cmuce*p1p3*C1ltmu + cmuce*p1p3*C0ltmu + mm2*cmuce
     .    *C11ltmu + 2*mm2*cmuce*C1ltmu - cImu*mm*p3ce*C12ltmu - cImu*
     .    mm*p3ce*C11ltmu - cImu*mm*p3ce*C1ltmu - cImu*mm*p1ce*C12ltmu
     .     - cImu*mm*p1ce*C11ltmu - cImu*mm*p1ce*C1ltmu )
      amplve =
     . + QRe * (  - 1.D0/2.D0*cmuce + cmuce*
     .    C00lte + cmuce*p2p4*C12lte + 2*cmuce*p2p4*C1lte + cmuce*p2p4*
     .    C0lte + me2*cmuce*C11lte + 2*me2*cmuce*C1lte - cIe*me*p4cmu*
     .    C12lte - cIe*me*p4cmu*C11lte - cIe*me*p4cmu*C1lte - cIe*me*
     .    p2cmu*C12lte - cIe*me*p2cmu*C11lte - cIe*me*p2cmu*C1lte )         
      else
c         include 'oneloop/vertexesplus.f'
         amplvmu =
     .  + QRmu * (  - 1.D0/2.D0*cmuce + cmuce*C00ltmu + cmuce*p1p3*
     .    C12ltmu + 2*cmuce*p1p3*C1ltmu + cmuce*p1p3*C0ltmu + mm2*cmuce
     .    *C11ltmu + 2*mm2*cmuce*C1ltmu + cImu*mm*p3ce*C12ltmu + cImu*
     .    mm*p3ce*C11ltmu + cImu*mm*p3ce*C1ltmu + cImu*mm*p1ce*C12ltmu
     .     + cImu*mm*p1ce*C11ltmu + cImu*mm*p1ce*C1ltmu )
         amplve =
     .    + QRe * (  - 1.D0/2.D0*cmuce + cmuce*
     .    C00lte + cmuce*p2p4*C12lte + 2*cmuce*p2p4*C1lte + cmuce*p2p4*
     .    C0lte + me2*cmuce*C11lte + 2*me2*cmuce*C1lte - cIe*me*p4cmu*
     .    C12lte - cIe*me*p4cmu*C11lte - cIe*me*p4cmu*C1lte - cIe*me*
     .    p2cmu*C12lte - cIe*me*p2cmu*C11lte - cIe*me*p2cmu*C1lte )
      endif
      amplve  = amplve * tm1 * cover
      amplvmu = amplvmu * tm1 * cover

      vertexes = amplve+amplvmu
      
c      vertexes = amplvertexes * tm1 * cover
!     !vertexes = vertexes * vpt ! done outside, it's better
      return
      end
*********************************************************
      real*16 function qdot(p1,p2)
      double precision p1(0:3),p2(0:3)
      qdot = p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)      
      return
      end
********************************************************      
      subroutine loadQEDboxesfuns
      implicit none
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
      
      integer*8 id,ic
      common/boxesfuns/
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,     
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc
      integer ifirst
      common/loadboxes/pard,parc,ifirst
      data ifirst /0/

      data 
     .     D0ltd,D1ltd,D2ltd,D3ltd,
     .     D0ltc,D1ltc,D2ltc,D3ltc,     
     .     D00ltd,D11ltd,D12ltd,D13ltd,
     .     D22ltd,D23ltd,D33ltd,
     .     D00ltc,D11ltc,D12ltc,D13ltc,
     .     D22ltc,D23ltc,D33ltc
     .     /22*(0.d0,0.d0)/

      if (QRmu.eq.0.or.QRe.eq.0) return
      
#ifdef COLLIER
      call Dfunscollier
      return
#else
      if (ifirst.eq.0) then
         pard(1)  = mm2
         pard(2)  = mm2
         pard(3)  = me2
         pard(4)  = me2
         pard(7)  = 0.d0
         pard(8)  = mm2
         pard(9)  = 0.d0
         pard(10) = me2

         parc = pard
         
         ifirst = 1
      endif
      pard(5) = 2.d0*mm2-2.d0*p1p3
      pard(6) = mm2+me2+2.d0*p3p4
      parc(5) = 2.d0*mm2-2.d0*p1p3
      parc(6) = mm2+me2-2.d0*p2p3
      
      id = Dget(pard(1)*uquad,
     .     pard(2)*uquad,
     .     pard(3)*uquad,
     .     pard(4)*uquad,
     .     pard(5)*uquad,
     .     pard(6)*uquad,
     .     pard(7)*uquad,
     .     pard(8)*uquad,
     .     pard(9)*uquad,
     .     pard(10)*uquad)
      D0ltd  = Dval(dd0,id)
      D1ltd  = Dval(dd1,id)
      D2ltd  = Dval(dd2,id)
      D3ltd  = Dval(dd3,id)
      D00ltd = Dval(dd00,id)
      D11ltd = Dval(dd11,id)
      D12ltd = Dval(dd12,id)
      D13ltd = Dval(dd13,id)
      D22ltd = Dval(dd22,id)
      D23ltd = Dval(dd23,id)
      D33ltd = Dval(dd33,id)
      
      ic = Dget(parc(1)*uquad,
     .     parc(2)*uquad,
     .     parc(3)*uquad,
     .     parc(4)*uquad,
     .     parc(5)*uquad,
     .     parc(6)*uquad,
     .     parc(7)*uquad,
     .     parc(8)*uquad,
     .     parc(9)*uquad,
     .     parc(10)*uquad)
      D0ltc  = Dval(dd0,ic)
      D1ltc  = Dval(dd1,ic)
      D2ltc  = Dval(dd2,ic)
      D3ltc  = Dval(dd3,ic)
      D00ltc = Dval(dd00,ic)
      D11ltc = Dval(dd11,ic)
      D12ltc = Dval(dd12,ic)
      D13ltc = Dval(dd13,ic)
      D22ltc = Dval(dd22,ic)
      D23ltc = Dval(dd23,ic)
      D33ltc = Dval(dd33,ic)      
#endif
      return
      end
**************************************************************************
      subroutine getQEDboxes(boxes)
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
      common/ifirstbox/cover,ifirst
#include "invariants.h"
      common/boxesfuns/
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
         include 'oneloop/boxesminus.f'
      else
         include 'oneloop/boxesplus.f'
      endif
      boxes = boxes * cover

!!!      boxes = boxes * vpt ! done outside, it's better
      
      return
      end
***************************************************************************
      subroutine loadinvariants(p1,p2,p3,p4)
#include "invariants.h"
      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p3p4 = dot(p3,p4)      
      return
      end
***************************************************************************
c$$$      subroutine loadvp1g1L
c$$$      implicit none
c$$$#include "invariants.h"
c$$$      integer ifirst
c$$$      common/ifirstvp1g1L/ifirst
c$$$      data ifirst /0/
c$$$      if (ifirst.eq.0) then
c$$$         ifirst = 1
c$$$         print*,'loadvp1g1L must be called after loadinvariants1g!!'
c$$$         print*,'For the moment, loadvp1g1L is a dummy routine'
c$$$      endif
c$$$!     here t13 and t24 are known
c$$$      vp24 = 1.d0
c$$$      vp13 = 1.d0      
c$$$      return
c$$$      end
***************************************************************************
      subroutine loadinvariants1g(p1o,p2o,p3o,p4o,ko)
      implicit none
#include "invariants.h"
      double precision k(0:3),ko(0:3) ! o for original
      double precision p1o(0:3),p2o(0:3),p3o(0:3),p4o(0:3)
      double precision kh(0:3) ! h for hammer
      double precision p1h(0:3),p2h(0:3),p3h(0:3),p4h(0:3)
      double precision pIn(0:3),pOut(0:3)
      logical correct,goon
      double precision zerocut,Eold,Enew,check
      integer j,lin,legs,i1,i2
      parameter (lin=2, legs=5)
c      parameter (zerocut = 5d-16)
c      parameter (zerocut = 1d-18)
      parameter (zerocut = 1d-15)
c      parameter (correct = .TRUE.)
      parameter (correct = .FALSE.)
      double precision tridot
      external tridot
      double precision p(0:3,legs),p0(0:3,legs),p1r(3),p2r(3)
      double precision qmur(0:3),qr(3),d1,d2,d3,delta
      double precision E,qsr,p1Xq,p2Xq,q1q(3),q2q(3),p12(3),p12Xq
      double precision m1,m2,E1,E2,qtm
      double precision p1q(3),p1t(3),p2q(3),p2t(3),p1ts,p2ts,p1s,p2s
      integer ifirst,ii
      common/inv1g/ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         do ii = 1,1
            if (Qmu.gt.0.d0.and.correct) then
               print*,'This routine with correct on'//
     .              ' does not work for crossing symmetry!!'
               print*,'Switch off correct or upgrade the routine!'
            endif
c            print*,'correct variable is',correct
         enddo
         ifirst = 1
      endif

      if (.not.correct) then ! default now, the "hammer" is essentially useless, also for pentagons
         
         p1p2  = dot(p1o,p2o)
         p1p3  = dot(p1o,p3o)
         p1p4  = dot(p1o,p4o)
         p2p3  = dot(p2o,p3o)
         p2p4  = dot(p2o,p4o)
         p3p4  = dot(p3o,p4o)
         p1k   = dot(p1o,ko)
         p2k   = dot(p2o,ko)
         p3k   = dot(p3o,ko)
         p4k   = dot(p4o,ko)
         s12   = mm2 + me2 + 2.d0*p1p2 !dot(p1+p2,p1+p2)
         s34   = mm2 + me2 + 2.d0*p3p4 !dot(p3+p4,p3+p4)
         t13   = 2.d0*mm2  - 2.d0*p1p3 !dot(p1-p3,p1-p3)
         t24   = 2.d0*me2  - 2.d0*p2p4 !dot(p2-p4,p2-p4)
         u14   = mm2 + me2 - 2.d0*p1p4 !dot(p1-p4,p1-p4)
         u23   = mm2 + me2 - 2.d0*p2p3 ! dot(p2-p3,p2-p3)
         p1mk2 = mm2 - 2.d0*p1k ! dot(p1-k,p1-k)
         p2mk2 = me2 - 2.d0*p2k !dot(p2-k,p2-k)
         p3pk2 = mm2 + 2.d0*p3k !dot(p3+k,p3+k)
         p4pk2 = me2 + 2.d0*p4k !dot(p4+k,p4+k)

      else ! if correct

         p1h = p1o
         p2h = p2o
         p3h = p3o
         p4h = p4o
         kh  = ko
         
         Eold = p1h(0)
         Enew = sqrt(tridot(p1h,p1h) + mm2 )
         if (abs(Eold-Enew)/Enew.gt.zerocut) then
c            print*,'p1 corrected'
            p1h(0) = Enew
         endif
         Eold = p2h(0)
         Enew = sqrt(tridot(p2h,p2h) + me2 )
         if (abs(Eold-Enew)/Enew.gt.zerocut) then
c            print*,'p2 corrected'
            p2h(0) = Enew
         endif
         Eold = p3h(0)
         Enew = sqrt(tridot(p3h,p3h) + mm2 )
         if (abs(Eold-Enew)/Enew.gt.zerocut) then
c            print*,'p3 corrected'
            p3h(0) = Enew
         endif
         Eold = p4h(0)
         Enew = sqrt(tridot(p4h,p4h) + me2 )
         if (abs(Eold-Enew)/Enew.gt.zerocut) then
c            print*,'p4 corrected'
            p4h(0) = Enew
         endif
         Eold = kh(0)
         Enew = sqrt(tridot(kh,kh)+0.d0)
         if (abs(Eold-Enew)/Enew.gt.zerocut) then
c            print*,'k corrected'
            kh(0) = Enew
         endif

         pIn  = p1h+p2h
         pOut = p3h+p4h+kh
         
         check = 0.d0
         do j = 0,3
            check = max(check,abs(pIn(j)-pOut(j))/abs(pIn(0)))
         enddo
         if (check.gt.zerocut) then
c            print*,'conservation corrected'

            p(0:3,1) = p1h
            p(0:3,2) = p2h
            p(0:3,3) = p3h
            p(0:3,4) = p4h
            p(0:3,5) = kh

            do i1 = 3,legs
               do  i2 = i1+1,legs
                  p1r = p(1:3,i1)
                  p2r = p(1:3,i2)
                  qmur = p(:,i1) + p(:,i2) - pOut(:) + pIn(:)

                  E = qmur(0)
                  qr = qmur(1:3)
                  qsr = tridot(qmur,qmur)

                  p1Xq = dot_product(p1r,qr)
                  p2Xq = dot_product(p2r,qr)

                  p12 = p1r+p2r
                  p12Xq = dot_product(p12,qr)

                  if (qsr.ne.0.d0) then
                     q1q = p1Xq/p12Xq*qr
                     q2q = p2Xq/p12Xq*qr
                  else
                     q1q = 0.d0
                     q2q = 0.d0
                  endif

                  if (i1.eq.3) then
                     m1 = mm
                  elseif(i1.eq.4) then
                     m1 = me
                  else
                     m1 = 0.d0
                  endif
                  if (i2.eq.4) then
                     m2 = me
                  else
                     m2 = 0.d0
                  endif

                  d1 = m1*m1 + dot_product(q1q,q1q)
                  d2 = m2*m2 + dot_product(q2q,q2q)
                  d3 = E*E
                  delta = d1*d1+d2*d2+d3*d3-2*d1*d2-2*d1*d3-2*d2*d3

                  if (d3+d1-d2.ge.0d0.and.
     .                 d3+d2-d1.ge.0d0.and.
     .                 delta.ge.0d0) then
                     goto 666
                  endif
               enddo
            enddo
 666        continue

cc            p0 = p ! for cross-checks

            E1  = (d3+d1-d2)/(2*E)
            E2  = (d3+d2-d1)/(2*E)
            qtm = sqrt(delta)/(2*E)

            p(0,i1) = E1
            p(0,i2) = E2
            if (qsr.gt.0d0) then
               p1q = p1Xq/qsr * qr
               p2q = p2Xq/qsr * qr
               p1t = p1r - p1q
               p1ts = dot_product(p1t,p1t)
               p2t = p2r - p2q
               p2ts = dot_product(p2t,p2t)
               p(1:3,i1) = q1q + qtm/sqrt(p1ts)*p1t
               p(1:3,i2) = q2q + qtm/sqrt(p2ts)*p2t
            else
               p1s = dot_product(p1r,p1r)
               p2s = dot_product(p2r,p2r)
               if (p1s.eq.0d0.and.p2s.eq.0d0) then
                  p(1:3,i1) = 0d0
                  p(1:3,i2) = 0d0
               elseif (p1s.ge.p2s) then
                  p(1:3,i1) = + qtm/sqrt(p1s)*p1r
                  p(1:3,i2) = - qtm/sqrt(p1s)*p1r
               else
                  p(1:3,i1) = - qtm/sqrt(p2s)*p2r
                  p(1:3,i2) = + qtm/sqrt(p2s)*p2r
               endif
            endif

c            print*,' '
c            print*,p0(0:3,3),dot(p0(:,3),p0(:,3))/mm2
c            print*,p(0:3,3),dot(p(:,3),p(:,3))/mm2
c            print*,p0(0:3,4),dot(p0(:,4),p0(:,4))/me2
c            print*,p(0:3,4),dot(p(:,4),p(:,4))/me2
c            print*,p0(:,3)+p0(:,4)+p0(:,5)-p0(:,1)-p0(:,2)
c            print*,p(:,3)+p(:,4)+p(:,5)-p(:,1)-p(:,2)

            p1h = p(:,1)
            p2h = p(:,2)
            p3h = p(:,3)
            p4h = p(:,4)
            kh  = p(:,5)
            
         endif
      
         p1p2  = dot(p1h,p2h)
         p1p3  = dot(p1h,p3h)
         p1p4  = dot(p1h,p4h)
         p2p3  = dot(p2h,p3h)
         p2p4  = dot(p2h,p4h)
         p3p4  = dot(p3h,p4h)
         p1k   = dot(p1h,kh)
         p2k   = dot(p2h,kh)
         p3k   = dot(p3h,kh)
         p4k   = dot(p4h,kh)
         s12   = mm2 + me2 + 2.d0*p1p2 !dot(p1+p2,p1+p2)
         s34   = mm2 + me2 + 2.d0*p3p4 !dot(p3+p4,p3+p4)
         t13   = 2.d0*mm2  - 2.d0*p1p3 !dot(p1-p3,p1-p3)
         t24   = 2.d0*me2  - 2.d0*p2p4 !dot(p2-p4,p2-p4)
         u14   = mm2 + me2 - 2.d0*p1p4 !dot(p1-p4,p1-p4)
         u23   = mm2 + me2 - 2.d0*p2p3 ! dot(p2-p3,p2-p3)
         p1mk2 = mm2 - 2.d0*p1k ! dot(p1-k,p1-k)
         p2mk2 = me2 - 2.d0*p2k !dot(p2-k,p2-k)
         p3pk2 = mm2 + 2.d0*p3k !dot(p3+k,p3+k)
         p4pk2 = me2 + 2.d0*p4k !dot(p4+k,p4+k)
         
      endif
      
      return
      end
**************************************************************************
      subroutine loadspinvariants(p1,p2,p3,p4,cmu,ce)
#include "invariants.h"
!     double precision f
!     c5X are in common!

      p1cmu = dot1cmom(p1,cmu)
      p2cmu = dot1cmom(p2,cmu)
      p3cmu = dot1cmom(p3,cmu)
      p4cmu = dot1cmom(p4,cmu)
      p1ce  = dot1cmom(p1,ce)
      p2ce  = dot1cmom(p2,ce)
      p3ce  = dot1cmom(p3,ce)
      p4ce  = dot1cmom(p4,ce)

      e13cmuc5e = epsfun2cmom(p1,p3,cmu,c5e)
      e14cmuc5e = epsfun2cmom(p1,p4,cmu,c5e)
      e14cec5mu = epsfun2cmom(p1,p4,ce,c5mu)
      e34cec5mu = epsfun2cmom(p3,p4,ce,c5mu)
      e12cmuc5e = epsfun2cmom(p1,p2,cmu,c5e)
      e12cec5mu = epsfun2cmom(p1,p2,ce,c5mu)
      e23cec5mu = epsfun2cmom(p2,p3,ce,c5mu)

c      print*,'<<<'
c      print*,e13cmuc5e
c      print*,e14cmuc5e
c      print*,e14cec5mu
c      print*,e34cec5mu
c      print*,e12cmuc5e
c      print*,e12cec5mu
c      print*,e23cec5mu

      p1c5mu = dot1cmom(p1,c5mu)
      p2c5mu = dot1cmom(p2,c5mu)
      p3c5mu = dot1cmom(p3,c5mu)
      p4c5mu = dot1cmom(p4,c5mu)
      p1c5e  = dot1cmom(p1,c5e)
      p2c5e  = dot1cmom(p2,c5e)
      p3c5e  = dot1cmom(p3,c5e)
      p4c5e  = dot1cmom(p4,c5e)

      c5muc5e = dotcmom(c5mu,c5e)            
      
c$$$      emy12cmuc5e = epsfuncmom(cuno*p1,cuno*p2,cmu,c5e)*im
c$$$      emy12cec5mu = epsfuncmom(cuno*p1,cuno*p2,ce,c5mu)*im
c$$$      emy13cmuc5e = epsfuncmom(cuno*p1,cuno*p3,cmu,c5e)*im
c$$$      emy14cmuc5e = epsfuncmom(cuno*p1,cuno*p4,cmu,c5e)*im
c$$$      emy14cec5mu = epsfuncmom(cuno*p1,cuno*p4,ce,c5mu)*im
c$$$      emy23cec5mu = epsfuncmom(cuno*p2,cuno*p3,ce,c5mu)*im
c$$$      emy34cec5mu = epsfuncmom(cuno*p3,cuno*p4,ce,c5mu)*im
      return
      end
***
      subroutine getcurrent(s4l,s2l,ce,cIe,c5e)
      implicit none
      double complex s4l(1:4),s2l(1:4),sa(1:4),sb(1:4),ce(0:3),im,cIe
      double complex c5e(0:3),s5a(1:4)
      parameter (im = (0.d0,1.d0))
*      double complex dot_product
*      external dot_product

** c(0) =        u1(1)d * u2(1) + u1(2)d * u2(2) + u1(3)d * u2(3) + u1(4)d * u2(4)
** c(1) =        u1(1)d * u2(4) + u1(2)d * u2(3) + u1(3)d * u2(2) + u1(4)d * u2(1)
** c(2) = im * (-u1(1)d * u2(4) + u1(2)d * u2(3) - u1(3)d * u2(2) + u1(4)d * u2(1))
** c(3) =        u1(1)d * u2(3) - u1(2)d * u2(4) + u1(3)d * u2(1) - u1(4)d * u2(2)      
      sa    = s4l

      s5a(1) = -sa(3)
      s5a(2) = -sa(4)
      s5a(3) = -sa(1)
      s5a(4) = -sa(2)
      
      sb    = s2l
      ce(0)= dot_product(sa,sb)

      c5e(0) = dot_product(s5a,sb)
      
      sb(1) = s2l(4)
      sb(2) = s2l(3)
      sb(3) = s2l(2)
      sb(4) = s2l(1)
      ce(1)= dot_product(sa,sb)

      c5e(1)= dot_product(s5a,sb)

      sb(1) = -im*sb(1)         ! occhio! e' sb definito appena sopra, non s2l!
      sb(2) =  im*sb(2)
      sb(3) = -im*sb(3)
      sb(4) =  im*sb(4)
      ce(2)=  dot_product(sa,sb)

      c5e(2)=  dot_product(s5a,sb)

      sb(1) =  s2l(3)
      sb(2) = -s2l(4)
      sb(3) =  s2l(1)
      sb(4) = -s2l(2)
      ce(3) = dot_product(sa,sb)

      c5e(3) = dot_product(s5a,sb)
      
      sb(1) =  s2l(1)
      sb(2) =  s2l(2)
      sb(3) = -s2l(3)
      sb(4) = -s2l(4)

      cIe = dot_product(sa,sb)
      
      return
      end
***************************************************
      subroutine getcurrentnew(ch,h1,h2,p1,p2,m1,m2,c,cI,c5)
      implicit none
      integer h1,h2
      double precision p1(0:3),p2(0:3),den12,den22,m1,m2,ch
      double complex c0,c1,c2,c3
      double complex c50,c51,c52,c53
      double precision nm1
      double complex c(0:3),im,cI,c5(0:3)
      parameter (im = (0.d0,1.d0))

      den12 = p1(0) + m1
      den22 = p2(0) + m2

      if (ch.lt.0.d0) then
      
      if (h1.eq.1.and.h2.eq.1) then
      c0 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c0 = c0 + den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c1 =
     &  + im * (  - p1(2)*den22 + p2(2)*den12 )
      c1 = c1 + p1(1)*den22 + p2(1)*den12

      c2 =
     &  + im * ( p1(1)*den22 - p2(1)*den12 )
      c2 = c2 + p1(2)*den22 + p2(2)*den12

      c3 =
     &  + p1(3)*den22 + p2(3)*den12

      c50 =
     &  - p1(3)*den22 - p2(3)*den12

      c51 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c51 = c51 - p1(1)*p2(3) - p1(3)*p2(1)

      c52 =
     &  + im * (  - p1(1)*p2(3) + p1(3)*p2(1) )
      c52 = c52 - p1(2)*p2(3) - p1(3)*p2(2)

      c53 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c53 = c53 - den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) - p1(3)*p2(3)

      cI =
     &  + im * (  - p1(1)*p2(2) + p1(2)*p2(1) )
      cI = cI + den12*den22 - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)
      elseif (h1.eq.1.and.h2.eq.2) then
      c0 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c0 = c0 - p1(1)*p2(3) + p1(3)*p2(1)

      c1 =
     &  + p1(3)*den22 - p2(3)*den12

      c2 =
     &  + im * (  - p1(3)*den22 + p2(3)*den12 )

      c3 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c3 = c3 - p1(1)*den22 + p2(1)*den12

      c50 =
     &  + im * ( p1(2)*den22 + p2(2)*den12 )
      c50 = c50 - p1(1)*den22 - p2(1)*den12

      c51 =
     &  + im * ( p1(1)*p2(2) + p1(2)*p2(1) )
      c51 = c51 - den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c52 =
     &  + im * ( den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) - p1(3)*p2(3)
     &     )
      c52 = c52 - p1(1)*p2(2) - p1(2)*p2(1)

      c53 =
     &  + im * ( p1(2)*p2(3) + p1(3)*p2(2) )
      c53 = c53 - p1(1)*p2(3) - p1(3)*p2(1)

      cI =
     &  + im * (  - p1(2)*p2(3) + p1(3)*p2(2) )
      cI = cI + p1(1)*p2(3) - p1(3)*p2(1)
      elseif (h1.eq.2.and.h2.eq.1) then         
      c0 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c0 = c0 + p1(1)*p2(3) - p1(3)*p2(1)

      c1 =
     &  - p1(3)*den22 + p2(3)*den12

      c2 =
     &  + im * (  - p1(3)*den22 + p2(3)*den12 )

      c3 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c3 = c3 + p1(1)*den22 - p2(1)*den12

      c50 =
     &  + im * (  - p1(2)*den22 - p2(2)*den12 )
      c50 = c50 - p1(1)*den22 - p2(1)*den12

      c51 =
     &  + im * (  - p1(1)*p2(2) - p1(2)*p2(1) )
      c51 = c51 - den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c52 =
     &  + im * (  - den12*den22 + p1(1)*p2(1) - p1(2)*p2(2) + p1(3)*
     &    p2(3) )
      c52 = c52 - p1(1)*p2(2) - p1(2)*p2(1)

      c53 =
     &  + im * (  - p1(2)*p2(3) - p1(3)*p2(2) )
      c53 = c53 - p1(1)*p2(3) - p1(3)*p2(1)

      cI =
     &  + im * (  - p1(2)*p2(3) + p1(3)*p2(2) )
      cI = cI - p1(1)*p2(3) + p1(3)*p2(1)
      elseif (h1.eq.2.and.h2.eq.2) then
      c0 =
     &  + im * (  - p1(1)*p2(2) + p1(2)*p2(1) )
      c0 = c0 + den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c1 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c1 = c1 + p1(1)*den22 + p2(1)*den12

      c2 =
     &  + im * (  - p1(1)*den22 + p2(1)*den12 )
      c2 = c2 + p1(2)*den22 + p2(2)*den12

      c3 =
     &  + p1(3)*den22 + p2(3)*den12

      c50 =
     &  + p1(3)*den22 + p2(3)*den12

      c51 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c51 = c51 + p1(1)*p2(3) + p1(3)*p2(1)

      c52 =
     &  + im * (  - p1(1)*p2(3) + p1(3)*p2(1) )
      c52 = c52 + p1(2)*p2(3) + p1(3)*p2(2)

      c53 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c53 = c53 + den12*den22 - p1(1)*p2(1) - p1(2)*p2(2) + p1(3)*p2(3)

      cI =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      cI = cI + den12*den22 - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)

      endif

      else ! if ch > 0

      if (h1.eq.1.and.h2.eq.1) then
      c0 =
     &  + im * (  - p1(1)*p2(2) + p1(2)*p2(1) )
      c0 = c0 + den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c1 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c1 = c1 + p1(1)*den22 + p2(1)*den12

      c2 =
     &  + im * (  - p1(1)*den22 + p2(1)*den12 )
      c2 = c2 + p1(2)*den22 + p2(2)*den12

      c3 =
     &  + p1(3)*den22 + p2(3)*den12

      c50 =
     &  + p1(3)*den22 + p2(3)*den12

      c51 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c51 = c51 + p1(1)*p2(3) + p1(3)*p2(1)

      c52 =
     &  + im * (  - p1(1)*p2(3) + p1(3)*p2(1) )
      c52 = c52 + p1(2)*p2(3) + p1(3)*p2(2)

      c53 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c53 = c53 + den12*den22 - p1(1)*p2(1) - p1(2)*p2(2) + p1(3)*p2(3)

      cI =
     &  + im * (  - p1(1)*p2(2) + p1(2)*p2(1) )
      cI = cI - den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)
         
      elseif (h1.eq.1.and.h2.eq.2) then
      c0 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c0 = c0 + p1(1)*p2(3) - p1(3)*p2(1)

      c1 =
     &  - p1(3)*den22 + p2(3)*den12

      c2 =
     &  + im * (  - p1(3)*den22 + p2(3)*den12 )

      c3 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c3 = c3 + p1(1)*den22 - p2(1)*den12

      c50 =
     &  + im * (  - p1(2)*den22 - p2(2)*den12 )
      c50 = c50 - p1(1)*den22 - p2(1)*den12

      c51 =
     &  + im * (  - p1(1)*p2(2) - p1(2)*p2(1) )
      c51 = c51 - den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c52 =
     &  + im * (  - den12*den22 + p1(1)*p2(1) - p1(2)*p2(2) + p1(3)*
     &    p2(3) )
      c52 = c52 - p1(1)*p2(2) - p1(2)*p2(1)

      c53 =
     &  + im * (  - p1(2)*p2(3) - p1(3)*p2(2) )
      c53 = c53 - p1(1)*p2(3) - p1(3)*p2(1)

      cI =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      cI = cI + p1(1)*p2(3) - p1(3)*p2(1)
         
      elseif (h1.eq.2.and.h2.eq.1) then         
      c0 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c0 = c0 - p1(1)*p2(3) + p1(3)*p2(1)

      c1 =
     &  + p1(3)*den22 - p2(3)*den12

      c2 =
     &  + im * (  - p1(3)*den22 + p2(3)*den12 )

      c3 =
     &  + im * ( p1(2)*den22 - p2(2)*den12 )
      c3 = c3 - p1(1)*den22 + p2(1)*den12

      c50 =
     &  + im * ( p1(2)*den22 + p2(2)*den12 )
      c50 = c50 - p1(1)*den22 - p2(1)*den12

      c51 =
     &  + im * ( p1(1)*p2(2) + p1(2)*p2(1) )
      c51 = c51 - den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c52 =
     &  + im * ( den12*den22 - p1(1)*p2(1) + p1(2)*p2(2) - p1(3)*p2(3)
     &     )
      c52 = c52 - p1(1)*p2(2) - p1(2)*p2(1)

      c53 =
     &  + im * ( p1(2)*p2(3) + p1(3)*p2(2) )
      c53 = c53 - p1(1)*p2(3) - p1(3)*p2(1)

      cI =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      cI = cI - p1(1)*p2(3) + p1(3)*p2(1)
         
      elseif (h1.eq.2.and.h2.eq.2) then
      c0 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c0 = c0 + den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)

      c1 =
     &  + im * (  - p1(2)*den22 + p2(2)*den12 )
      c1 = c1 + p1(1)*den22 + p2(1)*den12

      c2 =
     &  + im * ( p1(1)*den22 - p2(1)*den12 )
      c2 = c2 + p1(2)*den22 + p2(2)*den12

      c3 =
     &  + p1(3)*den22 + p2(3)*den12

      c50 =
     &  - p1(3)*den22 - p2(3)*den12

      c51 =
     &  + im * ( p1(2)*p2(3) - p1(3)*p2(2) )
      c51 = c51 - p1(1)*p2(3) - p1(3)*p2(1)

      c52 =
     &  + im * (  - p1(1)*p2(3) + p1(3)*p2(1) )
      c52 = c52 - p1(2)*p2(3) - p1(3)*p2(2)

      c53 =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      c53 = c53 - den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) - p1(3)*p2(3)

      cI =
     &  + im * ( p1(1)*p2(2) - p1(2)*p2(1) )
      cI = cI - den12*den22 + p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)         
      endif

      endif
         
      c(0) = c0
      c(1) = c1
      c(2) = c2
      c(3) = c3
      c5(0) = c50
      c5(1) = c51
      c5(2) = c52
      c5(3) = c53
      
      nm1 = 1.d0/sqrt(den12*den22)
      
      c  = c*nm1
      c5 = c5*nm1
      cI = cI*nm1
      
      return
      end
***************************************************

      
*** Levi-Civita tensor contracted with 4 COMPLEX momenta
      double complex function epsfuncmom(q1,q2,q3,q4)
      implicit none
      double complex q1(0:3),q2(0:3),q3(0:3),q4(0:3)
      double complex e,e1
      
c      e = q1(1)*q2(2)*q3(3)*q4(0) - q1(1)*q2(2)*q3(0)*q4(3) +
c     .    q1(1)*q2(3)*q3(0)*q4(2) - q1(1)*q2(3)*q3(2)*q4(0) +
c     .    q1(1)*q2(0)*q3(2)*q4(3) - q1(1)*q2(0)*q3(3)*q4(2) +
c     .    q1(2)*q2(1)*q3(0)*q4(3) - q1(2)*q2(1)*q3(3)*q4(0) +
c     .    q1(2)*q2(3)*q3(1)*q4(0) - q1(2)*q2(3)*q3(0)*q4(1) +
c     .    q1(2)*q2(0)*q3(3)*q4(1) - q1(2)*q2(0)*q3(1)*q4(3) +
c     .    q1(3)*q2(1)*q3(2)*q4(0) - q1(3)*q2(1)*q3(0)*q4(2) +
c     .    q1(3)*q2(2)*q3(0)*q4(1) - q1(3)*q2(2)*q3(1)*q4(0) +
c     .    q1(3)*q2(0)*q3(1)*q4(2) - q1(3)*q2(0)*q3(2)*q4(1) +
c     .    q1(0)*q2(1)*q3(3)*q4(2) - q1(0)*q2(1)*q3(2)*q4(3) +
c     .    q1(0)*q2(2)*q3(1)*q4(3) - q1(0)*q2(2)*q3(3)*q4(1) +
c     .    q1(0)*q2(3)*q3(2)*q4(1) - q1(0)*q2(3)*q3(1)*q4(2)

      
      e = q1(1) * (q2(2)*(q3(3)*q4(0) - q3(0)*q4(3))  +
     .             q2(3)*(q3(0)*q4(2) - q3(2)*q4(0))  +
     .             q2(0)*(q3(2)*q4(3) - q3(3)*q4(2))) +
     .    q1(2) * (q2(1)*(q3(0)*q4(3) - q3(3)*q4(0))  +
     .             q2(3)*(q3(1)*q4(0) - q3(0)*q4(1))  +
     .             q2(0)*(q3(3)*q4(1) - q3(1)*q4(3))) +
     .    q1(3) * (q2(1)*(q3(2)*q4(0) - q3(0)*q4(2))  +
     .             q2(2)*(q3(0)*q4(1) - q3(1)*q4(0))  +
     .             q2(0)*(q3(1)*q4(2) - q3(2)*q4(1))) +
     .    q1(0) * (q2(1)*(q3(3)*q4(2) - q3(2)*q4(3))  +
     .             q2(2)*(q3(1)*q4(3) - q3(3)*q4(1))  +
     .             q2(3)*(q3(2)*q4(1) - q3(1)*q4(2)))

      
      epsfuncmom = -(0.d0,1.d0) * e!! e chittemmuorto!!!
      return
      end
**
***************************************************
*** Levi-Civita tensor contracted with 2 real and 2 COMPLEX momenta
      double complex function epsfun2cmom(q1,q2,q3,q4)
      implicit none
      double complex q3(0:3),q4(0:3)
      double precision q1(0:3),q2(0:3)
      double complex e,ev(0:3),dot1cmom
      external dot1cmom
c      e = q1(1)*q2(2)*q3(3)*q4(0) - q1(1)*q2(2)*q3(0)*q4(3) +
c     .    q1(1)*q2(3)*q3(0)*q4(2) - q1(1)*q2(3)*q3(2)*q4(0) +
c     .    q1(1)*q2(0)*q3(2)*q4(3) - q1(1)*q2(0)*q3(3)*q4(2) +
c     .    q1(2)*q2(1)*q3(0)*q4(3) - q1(2)*q2(1)*q3(3)*q4(0) +
c     .    q1(2)*q2(3)*q3(1)*q4(0) - q1(2)*q2(3)*q3(0)*q4(1) +
c     .    q1(2)*q2(0)*q3(3)*q4(1) - q1(2)*q2(0)*q3(1)*q4(3) +
c     .    q1(3)*q2(1)*q3(2)*q4(0) - q1(3)*q2(1)*q3(0)*q4(2) +
c     .    q1(3)*q2(2)*q3(0)*q4(1) - q1(3)*q2(2)*q3(1)*q4(0) +
c     .    q1(3)*q2(0)*q3(1)*q4(2) - q1(3)*q2(0)*q3(2)*q4(1) +
c     .    q1(0)*q2(1)*q3(3)*q4(2) - q1(0)*q2(1)*q3(2)*q4(3) +
c     .    q1(0)*q2(2)*q3(1)*q4(3) - q1(0)*q2(2)*q3(3)*q4(1) +
c     .    q1(0)*q2(3)*q3(2)*q4(1) - q1(0)*q2(3)*q3(1)*q4(2)
            
      e = q1(1) * (q2(2)*(q3(3)*q4(0) - q3(0)*q4(3))  +
     .             q2(3)*(q3(0)*q4(2) - q3(2)*q4(0))  +
     .             q2(0)*(q3(2)*q4(3) - q3(3)*q4(2))) +
     .    q1(2) * (q2(1)*(q3(0)*q4(3) - q3(3)*q4(0))  +
     .             q2(3)*(q3(1)*q4(0) - q3(0)*q4(1))  +
     .             q2(0)*(q3(3)*q4(1) - q3(1)*q4(3))) +
     .    q1(3) * (q2(1)*(q3(2)*q4(0) - q3(0)*q4(2))  +
     .             q2(2)*(q3(0)*q4(1) - q3(1)*q4(0))  +
     .             q2(0)*(q3(1)*q4(2) - q3(2)*q4(1))) +
     .    q1(0) * (q2(1)*(q3(3)*q4(2) - q3(2)*q4(3))  +
     .             q2(2)*(q3(1)*q4(3) - q3(3)*q4(1))  +
     .             q2(3)*(q3(2)*q4(1) - q3(1)*q4(2)))

c      call getev(ev,q2,q3,q4)
c      print*,'here',e-dot1cmom(q1,ev)

      
      epsfun2cmom = -(0.d0,1.d0) * e!! e chittemmuorto!!!
      return
      end
*****************************************
***************************************************
*** Levi-Civita tensor contracted with 2 real and 2 COMPLEX momenta
      subroutine getev(ev,q2,q3,q4)
      implicit none
      double complex q3(0:3),q4(0:3)
      double precision q1(0:3),q2(0:3)
      double complex e,ev(0:3)
c      e = q1(1)*q2(2)*q3(3)*q4(0) - q1(1)*q2(2)*q3(0)*q4(3) +
c     .    q1(1)*q2(3)*q3(0)*q4(2) - q1(1)*q2(3)*q3(2)*q4(0) +
c     .    q1(1)*q2(0)*q3(2)*q4(3) - q1(1)*q2(0)*q3(3)*q4(2) +
c     .    q1(2)*q2(1)*q3(0)*q4(3) - q1(2)*q2(1)*q3(3)*q4(0) +
c     .    q1(2)*q2(3)*q3(1)*q4(0) - q1(2)*q2(3)*q3(0)*q4(1) +
c     .    q1(2)*q2(0)*q3(3)*q4(1) - q1(2)*q2(0)*q3(1)*q4(3) +
c     .    q1(3)*q2(1)*q3(2)*q4(0) - q1(3)*q2(1)*q3(0)*q4(2) +
c     .    q1(3)*q2(2)*q3(0)*q4(1) - q1(3)*q2(2)*q3(1)*q4(0) +
c     .    q1(3)*q2(0)*q3(1)*q4(2) - q1(3)*q2(0)*q3(2)*q4(1) +
c     .    q1(0)*q2(1)*q3(3)*q4(2) - q1(0)*q2(1)*q3(2)*q4(3) +
c     .    q1(0)*q2(2)*q3(1)*q4(3) - q1(0)*q2(2)*q3(3)*q4(1) +
c     .    q1(0)*q2(3)*q3(2)*q4(1) - q1(0)*q2(3)*q3(1)*q4(2)

      ev(0) = 
     .      + q2(1)*(q3(3)*q4(2) - q3(2)*q4(3))
     .      + q2(2)*(q3(1)*q4(3) - q3(3)*q4(1))
     .      + q2(3)*(q3(2)*q4(1) - q3(1)*q4(2))
      
      ev(1) =
     .      - q2(2)*(q3(3)*q4(0) - q3(0)*q4(3))
     .      - q2(3)*(q3(0)*q4(2) - q3(2)*q4(0))
     .      - q2(0)*(q3(2)*q4(3) - q3(3)*q4(2))

      ev(2) =
     .      - q2(1)*(q3(0)*q4(3) - q3(3)*q4(0))
     .      - q2(3)*(q3(1)*q4(0) - q3(0)*q4(1))
     .      - q2(0)*(q3(3)*q4(1) - q3(1)*q4(3))

      ev(3) = 
     .      - q2(1)*(q3(2)*q4(0) - q3(0)*q4(2))
     .      - q2(2)*(q3(0)*q4(1) - q3(1)*q4(0))
     .      - q2(0)*(q3(1)*q4(2) - q3(2)*q4(1))

**--- the minus is already above
cc      ev(1:3) = -ev(1:3)

      ev = -(0.d0,1.d0)*ev
      
      return
      end
*********************************************
      double complex function dotcmom(p,q)
      implicit none
      double complex p(0:3),q(0:3)
      dotcmom=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end
********************************************
      double complex function dot1cmom(p,q)
      implicit none
      double complex q(0:3)
      double precision p(0:3)
      dot1cmom=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end
*************** PASSERA EXP ******************************
      subroutine KMP(t,s,oal,oal2,resum)
      implicit double precision (a-h,l,m,o-z)
      common/parameters/ame,ammu,convfac,alpha,pi
      me2  = ame*ame
      mm2  = ammu*ammu
      tmin = -lambda(s,me2,mm2)/s
      twoaopi = 2.d0*alpha/pi
      alog    = log((tmin-t)*2.d0/tmin)
      sq = sqrt(t*(t-4.d0*mm2))
      cm = log(-t/me2) - 2.d0 + 2.d0*log(1.d0-t/(mm2-s)) +
     .     (2.d0*mm2-t)/sq*log((sq+2d0*mm2-t)*0.5d0/mm2)
      oal   = 1.d0 + twoaopi*cm*alog
      oal2  = oal + 0.5d0*(oal-1.d0)**2
      resum = exp(twoaopi*cm*alog)
      return
      end
***********************************************************
      subroutine YFSvirt_works_perfectly(YFS)
! This works, I try to separate the ii part in the default routine...      
!     https://s3.cern.ch/inspire-prod-files-2/2977b9e588a3dbe4077d3ae5016e31c5
!     https://arxiv.org/pdf/hep-ph/0603034.pdf
!     https://arxiv.org/pdf/1809.10650.pdf
!     https://arxiv.org/pdf/0810.5071.pdf
!     https://arxiv.org/pdf/hep-ph/0006359.pdf      
#ifdef COLLIER
      use collier
#endif
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      dimension pi(0:3),pj(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)      
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirstaaa
      common/pmatrix/pmat
      common/parameters/ame,ammu,convfac,alpha,pig
      double precision z(4),th(4)
      double complex bdummy(4)
      common/yfscmn/overall,z,th,bdummy,ifirst
      data ifirst/0/
      double complex c,C0singular
      double complex B0ii,B0ij,B0jj,C0ij
      external C0singular
      parameter (c = (1.d0,0.d0))
      double precision phmass
      common/photonmasslambda/phmass      

      double precision Qmu
      common/muoncharge/Qmu
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
      
#ifdef QUAD
      real*16 uquad
      parameter (uquad=1.q0)
#else
      double precision uquad
      parameter (uquad=1.d0)
#endif

#ifndef COLLIER           
#include "looptools.h"
      ComplexType Cfuns(Ncc),Bfuns(Nbb)
#else
      integer rank              ! max rank
      parameter (rank = 0)
      double complex Bfuns(0:rank/2,0:rank)
      double complex Bfunsuv(0:rank/2,0:rank)
#endif

      if (ifirst.eq.0) then
         do k = 1,1
c            print*,'YFS virtual must be called after soft_integral_mue!'
         enddo
         overall = 1.d0/pig/pig/16.d0  * alpha/8.d0/pig/pig/pig

         overall = alpha /4.d0/pig

         z(1) =  Qmu*QRmu
         z(2) = -1.d0*QRe
         z(3) =  Qmu*QRmu
         z(4) = -1.d0*QRe

         th(1) = -1.d0
         th(2) = -1.d0
         th(3) =  1.d0
         th(4) =  1.d0

         ifirst = 1
      endif

      npartl = npart

      YFS = 0.d0
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = z(i)*z(j)*th(i)*th(j)
            YFSij = 0.d0
            if (abs(etaij).gt.0.d0) then
               pi  = th(i)*pmat(i,:)
               pj  = th(j)*pmat(j,:)
               mi2 = masses(i)
               mj2 = masses(j)

               pi_pj = dot(pi,pj)
               s = mi2 + mj2 + 2.d0*pi_pj
#ifdef COLLIER
               call B0_cll(B0ii,c*0.d0,c*mi2,c*mi2)
               call B0_cll(B0jj,c*0.d0,c*mj2,c*mj2)
               call B0_cll(B0ij,c*s,c*mi2,c*mj2)
c      call C_cll(Bfuns,Bfunsuv,c*mi2,c*s,c*mj2,c*0.d0,c*mi2,c*mj2,rank) ! bah, non funziona
              call C0_cll(C0ij,c*mi2,c*s,c*mj2,c*0.d0,c*mi2,c*mj2)
#else
               B0ii = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
               B0jj = B0(uquad*0.d0,uquad*mj2,uquad*mj2)
               B0ij = B0(uquad*s,uquad*mi2,uquad*mj2)
               C0ij = C0(uquad*mi2,uquad*s,uquad*mj2,
     .              uquad*0.d0,uquad*mi2,uquad*mj2)               
#endif               
               YFSij = 2.d0*B0ij - B0ii - B0jj
     .         + 4.d0 * mi2 * C0singular(mi2,phmass)
     .         + 4.d0 * mj2 * C0singular(mj2,phmass)
     .         + 8.d0*pi_pj * C0ij
            endif
            YFS = YFS + YFSij*etaij
         enddo
      enddo

      YFS = overall * YFS

      return
      end
***********************************************************
      subroutine YFSvirt(YFSmat)
! I try to separate the ii part here... yep, it's identical to the previous one
!     https://s3.cern.ch/inspire-prod-files-2/2977b9e588a3dbe4077d3ae5016e31c5
!     https://arxiv.org/pdf/hep-ph/0603034.pdf
!     https://arxiv.org/pdf/1809.10650.pdf
!     https://arxiv.org/pdf/0810.5071.pdf
!     https://arxiv.org/pdf/hep-ph/0006359.pdf      
#ifdef COLLIER
      use collier
#endif
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      dimension pi(0:3),pj(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)      
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirstaaa
      common/pmatrix/pmat
      common/parameters/ame,ammu,convfac,alpha,pig
      double precision z(4),th(4)
      double complex B0ii(4),B0ij,B0jj,C0ij
      common/yfscmn/overall,z,th,B0ii,ifirst
      data ifirst/0/
      double complex c,C0singular,cYFSij,cYFSii
      external C0singular
      parameter (c = (1.d0,0.d0))
      double precision phmass
      common/photonmasslambda/phmass

      double precision Qmu
      common/muoncharge/Qmu
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu

      double complex YFSmat(4,4),YFSij,YFSii
      integer maskud(4,4),masknoud(4,4),maske(4,4),maskmu(4,4)
      common/yfsmasks/maskud,masknoud,maske,maskmu
      
      double complex C0sd,csd
      external C0sd
      
#ifdef QUAD
      real*16 uquad
      parameter (uquad=1.q0)
#else
      double precision uquad
      parameter (uquad=1.d0)
#endif

#ifndef COLLIER           
#include "looptools.h"
      ComplexType Cfuns(Ncc),Bfuns(Nbb)
#else
      integer rank              ! max rank
      parameter (rank = 0)
      double complex Bfuns(0:rank/2,0:rank)
      double complex Bfunsuv(0:rank/2,0:rank)
#endif

      if (ifirst.eq.0) then
         do k = 1,1
c            print*,'YFS virtual must be called after soft_integral_mue!'
         enddo

         YFSmat = 0.d0
         
         maskud   = 0
         masknoud = 0
         maske  = 0
         maskmu = 0
         
         maskud(1,2) = 1
         maskud(1,4) = 1
         maskud(2,3) = 1
         maskud(3,4) = 1
         maskud(2,1) = 1
         maskud(4,1) = 1
         maskud(3,2) = 1
         maskud(4,3) = 1
         
         masknoud(1,1) = 1
         masknoud(2,2) = 1
         masknoud(3,3) = 1
         masknoud(4,4) = 1

         masknoud(1,3) = 1
         masknoud(3,1) = 1
         masknoud(2,4) = 1
         masknoud(4,2) = 1

         maskmu(1,1) = 1
         maske(2,2)  = 1
         maskmu(3,3) = 1
         maske(4,4)  = 1

         maskmu(1,3) = 1
         maskmu(3,1) = 1
         maske(2,4)  = 1
         maske(4,2)  = 1
         
         overall = alpha /4.d0/pig

         z(1) =  Qmu*QRmu
         z(2) = -1.d0*QRe
         z(3) =  Qmu*QRmu
         z(4) = -1.d0*QRe

         th(1) = -1.d0
         th(2) = -1.d0
         th(3) =  1.d0
         th(4) =  1.d0

#ifdef COLLIER
         call SwitchOffCacheSystem_cll
         mi2 = masses(1)
         call B0_cll(B0ii(1),c*0.d0,c*mi2,c*mi2)
         mi2 = masses(2)
         call B0_cll(B0ii(2),c*0.d0,c*mi2,c*mi2)
         mi2 = masses(3)
         call B0_cll(B0ii(3),c*0.d0,c*mi2,c*mi2)
         mi2 = masses(4)
         call B0_cll(B0ii(4),c*0.d0,c*mi2,c*mi2)
         call SwitchOnCacheSystem_cll
#else
         call markcache
         mi2 = masses(1)
         B0ii(1) = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
         mi2 = masses(2)
         B0ii(2) = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
         mi2 = masses(3)
         B0ii(3) = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
         mi2 = masses(4)
         B0ii(4) = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
         call restorecache
#endif               
         ifirst = 1
      endif

      npartl = npart

! non diagonal contributions
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = z(i)*z(j)*th(i)*th(j)
            YFSij = 0.d0
            if (abs(etaij).gt.0.d0) then
               pi  = th(i)*pmat(i,:)
               pj  = th(j)*pmat(j,:)
               mi2 = masses(i)
               mj2 = masses(j)

               pi_pj = dot(pi,pj)
               s = mi2 + mj2 + 2.d0*pi_pj
#ifdef COLLIER
               call B0_cll(B0ij,c*s,c*mi2,c*mj2)
               call C0_cll(C0ij,c*mi2,c*s,c*mj2,c*0.d0,c*mi2,c*mj2)
c               csd = C0sd(s,mi2,mj2,phmass*phmass)
c               print*,s,C0ij,csd,C0ij-csd!,realpart(csd)/realpart(C0ij),
c     .              aimag(csd)/aimag(C0ij) ! OK
#else
               B0ij = B0(uquad*s,uquad*mi2,uquad*mj2)
               C0ij = C0(uquad*mi2,uquad*s,uquad*mj2,
     .              uquad*0.d0,uquad*mi2,uquad*mj2)
#endif
               YFSij = 2.d0*B0ij + 8.d0*pi_pj * C0ij
            endif
            YFSmat(i,j) = YFSij*etaij
         enddo
      enddo

! diagonal contributions      
      do i = 1,npartl
         YFSii = 0.d0
         etaii = z(i)*z(i)
         if (abs(etaii).gt.0.d0) then
            mi2 = masses(i)
            YFSii = 4.d0*mi2*C0singular(mi2,phmass) - B0ii(i)
         endif
         YFSmat(i,i)    =  - YFSii*etaii
      enddo
*
      YFSmat = overall*YFSmat
      return
      end      
***********************************************************************************************
      double complex function C0singular(m2,l)
! here I assume l going to 0 and p2 = m2
      implicit double precision (a-h,l,m,o-z)
! https://it.wikipedia.org/wiki/Tavola_degli_integrali_indefiniti_di_funzioni_razionali
      a  = m2
      b  = -2.d0*m2
      c  = m2
!      abc = a + b + c
      abc = l*l
      dueapb = l*l
      
!      disc = 4.d0*a*c - b*b
!      disc = l*l * (4.d0*m2 - l*l)
      
!      sq = sqrt(disc)
!      sq = l*sqrt(4.d0*m2 - l*l) 

! not needed if XX is commented      sq = 2.d0*l*sqrt(m2)      ! I set to 0 the term in l4

      C0singular = 0.5d0/a*( log(abc/c))
!XX     .     + dueapb/a/sq*atan(dueapb/sq) ! questi due sono prop a l, quindi vanno messi a 0!!
!XX     .     - dueapb/a/sq*atan(b/sq)

      return
      end
*********************************************************************************************
      function C0sd(s,m12,m22,l2)
! B.5 from hep-ph/0308246
      implicit none
      double complex C0sd
      double precision s,m1,m2,l
      double precision m12,m22,l2

      double precision pi,pi2o6,eps,xx

      double precision xs
      double precision sq
      double precision omxs2
      double complex lxs,i,li21,li22
      
      parameter (pi    = 4.d0*atan(1.d0))
      parameter (pi2o6 = pi*pi/6.d0)
      parameter (eps = 1.d-24)
c      parameter (eps = 0.d0)
      parameter (i = (0.d0,1.d0))
      
      double precision dilogmy,z
      external dilogmy
      double complex myli2
      external myli2
      
      m1 = sqrt(m12)
      m2 = sqrt(m22)
      l  = sqrt(l2)

      sq = sqrt(1.d0-4.d0*m1*m2/(s-(m1-m2)**2))
      
      xs = (sq - 1.d0)/(sq + 1.d0)

      omxs2 = 1.d0 - xs*xs
      
      lxs = log(dcmplx(xs))

! needed to rewrite dilogmy(1-z)!!      
      z = xs*m1/m2
      li21 = -dilogmy(z)+pi2o6 - log(dcmplx(z))*log(dcmplx(1.d0-z))
      z = xs*m2/m1
      li22 = -dilogmy(z)+pi2o6 - log(dcmplx(z))*log(dcmplx(1.d0-z))
!1      
      C0sd = -log(l2/m1/m2)*lxs
     . -0.5d0*lxs*lxs +2.d0*lxs*log(omxs2) +0.5d0*(log(m1/m2))**2-pi2o6
c     .     +myli2(dcmplx(xs*xs)) + myli2(dcmplx(1.d0-xs*m1/m2))
c     . + myli2(dcmplx(1.d0-xs*m2/m1)) ! complex myli2 is needed!
     .     +dilogmy(xs*xs) + li21 + li22

      C0sd = xs/m1/m2/omxs2*C0sd
      
      return
      end
***********************************************************
      function C0sd_onlyIR(s,m12,m22,l2)
! B.5 from hep-ph/0308246
      implicit none
      double complex C0sd_onlyIR
      double precision s,m1,m2,l
      double precision m12,m22,l2

      double precision pi,pi2o6,eps,xx

      double precision xs
      double precision sq
      double precision omxs2
      double complex lxs,i,li21,li22
      
      parameter (pi    = 4.d0*atan(1.d0))
      parameter (pi2o6 = pi*pi/6.d0)
      parameter (eps = 1.d-24)
c      parameter (eps = 0.d0)
      parameter (i = (0.d0,1.d0))
      
      double precision dilogmy,z
      external dilogmy
      double complex myli2
      external myli2
      
      m1 = sqrt(m12)
      m2 = sqrt(m22)
      l  = sqrt(l2)

      sq = sqrt(1.d0-4.d0*m1*m2/(s-(m1-m2)**2))
      
      xs = (sq - 1.d0)/(sq + 1.d0)

c      print*,s,m1,m2,-m1*m2/(s-m22)/xs

      omxs2 = 1.d0 - xs*xs
      
      lxs = log(dcmplx(xs))

! needed to rewrite dilogmy(1-z)!!      
      z = xs*m1/m2
      li21 = -dilogmy(z)+pi2o6 - log(dcmplx(z))*log(dcmplx(1.d0-z))
      z = xs*m2/m1
      li22 = -dilogmy(z)+pi2o6 - log(dcmplx(z))*log(dcmplx(1.d0-z))
!      print*,s,m1,m2,li21,li22,dilogmy(xs*xs)
!!
      C0sd_onlyIR = -log(l2/m1/m2)*lxs
     .     -0.5d0*lxs*lxs                     *1.d0 ! 1 
     .     +2.d0*lxs*log(omxs2)               *0.d0 ! 2
     .    +0.5d0*(log(m1/m2))**2              *1.d0 ! 3 
     .     +dilogmy(xs*xs) + li21 + li22-pi2o6 *0.d0 ! 4

c      print*,''
c      print*,-0.5d0*lxs*lxs
c      print*,2.d0*lxs*log(omxs2)
c      print*,0.5d0*(log(m1/m2))**2
c      print*,dilogmy(xs*xs) + li21 + li22-pi2o6
      
      
      C0sd_onlyIR = xs/m1/m2/omxs2*C0sd_onlyIR
      
      return
      end
***********************************************************
      subroutine YFSvirt_onlyIR(YFSmat)
! I try to separate the ii part here... yep, it's identical to the previous one
!     https://s3.cern.ch/inspire-prod-files-2/2977b9e588a3dbe4077d3ae5016e31c5
!     https://arxiv.org/pdf/hep-ph/0603034.pdf
!     https://arxiv.org/pdf/1809.10650.pdf
!     https://arxiv.org/pdf/0810.5071.pdf
!     https://arxiv.org/pdf/hep-ph/0006359.pdf      
#ifdef COLLIER
      use collier
#endif
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      dimension pi(0:3),pj(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)      
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirstaaa
      common/pmatrix/pmat
      common/parameters/ame,ammu,convfac,alpha,pig
      double precision z(4),th(4)
      double complex B0ii(4),B0ij,B0jj,C0ij
      common/yfsIRcmn/overall,z,th,B0ii,ifirst
      data ifirst/0/
      double complex c,C0singular,cYFSij,cYFSii
      external C0singular
      parameter (c = (1.d0,0.d0))
      double precision phmass
      common/photonmasslambda/phmass

      double precision Qmu
      common/muoncharge/Qmu
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu

      double complex YFSmat(4,4),YFSij,YFSii
      integer maskud(4,4),masknoud(4,4),maske(4,4),maskmu(4,4)
      common/yfsIRmasks/maskud,masknoud,maske,maskmu
      
      double complex C0sd_onlyIR,csd,C0sd
      external C0sd_onlyIR,C0sd
      
#ifdef QUAD
      real*16 uquad
      parameter (uquad=1.q0)
#else
      double precision uquad
      parameter (uquad=1.d0)
#endif

#ifndef COLLIER           
#include "looptools.h"
      ComplexType Cfuns(Ncc),Bfuns(Nbb)
#else
      integer rank              ! max rank
      parameter (rank = 0)
      double complex Bfuns(0:rank/2,0:rank)
      double complex Bfunsuv(0:rank/2,0:rank)
#endif

      if (ifirst.eq.0) then
         do k = 1,1
c            print*,'YFS virtual must be called after soft_integral_mue!'
         enddo

         YFSmat = 0.d0
         
         maskud   = 0
         masknoud = 0
         maske  = 0
         maskmu = 0
         
         maskud(1,2) = 1
         maskud(1,4) = 1
         maskud(2,3) = 1
         maskud(3,4) = 1
         maskud(2,1) = 1
         maskud(4,1) = 1
         maskud(3,2) = 1
         maskud(4,3) = 1
         
         masknoud(1,1) = 1
         masknoud(2,2) = 1
         masknoud(3,3) = 1
         masknoud(4,4) = 1

         masknoud(1,3) = 1
         masknoud(3,1) = 1
         masknoud(2,4) = 1
         masknoud(4,2) = 1

         maskmu(1,1) = 1
         maske(2,2)  = 1
         maskmu(3,3) = 1
         maske(4,4)  = 1

         maskmu(1,3) = 1
         maskmu(3,1) = 1
         maske(2,4)  = 1
         maske(4,2)  = 1
         
         overall = alpha /4.d0/pig

         z(1) =  Qmu*QRmu
         z(2) = -1.d0*QRe
         z(3) =  Qmu*QRmu
         z(4) = -1.d0*QRe

         th(1) = -1.d0
         th(2) = -1.d0
         th(3) =  1.d0
         th(4) =  1.d0

#ifdef COLLIER
         call SwitchOffCacheSystem_cll
         mi2 = masses(1)
         call B0_cll(B0ii(1),c*0.d0,c*mi2,c*mi2)
         mi2 = masses(2)
         call B0_cll(B0ii(2),c*0.d0,c*mi2,c*mi2)
         mi2 = masses(3)
         call B0_cll(B0ii(3),c*0.d0,c*mi2,c*mi2)
         mi2 = masses(4)
         call B0_cll(B0ii(4),c*0.d0,c*mi2,c*mi2)
         call SwitchOnCacheSystem_cll
#else
         call markcache
         mi2 = masses(1)
         B0ii(1) = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
         mi2 = masses(2)
         B0ii(2) = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
         mi2 = masses(3)
         B0ii(3) = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
         mi2 = masses(4)
         B0ii(4) = B0(uquad*0.d0,uquad*mi2,uquad*mi2)
         call restorecache
#endif
!!!
ccc         B0ii = (0.d0,0.d0)
!!!          
         
         ifirst = 1
      endif
    
      npartl = npart

! non diagonal contributions
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = z(i)*z(j)*th(i)*th(j)
            YFSij = 0.d0
            if (abs(etaij).gt.0.d0) then
               pi  = th(i)*pmat(i,:)
               pj  = th(j)*pmat(j,:)
               mi2 = masses(i)
               mj2 = masses(j)

               pi_pj = dot(pi,pj)
               s = mi2 + mj2 + 2.d0*pi_pj

               B0ij = (0.d0,0.d0)
               C0ij = C0sd_onlyIR(s,mi2,mj2,phmass*phmass)
               
#ifdef COLLIER
               call B0_cll(B0ij,c*s,c*mi2,c*mj2)
cc               call C0_cll(C0ij,c*mi2,c*s,c*mj2,c*0.d0,c*mi2,c*mj2)
#else
               B0ij = B0(uquad*s,uquad*mi2,uquad*mj2)
cc               C0ij = C0(uquad*mi2,uquad*s,uquad*mj2,
cc     .              uquad*0.d0,uquad*mi2,uquad*mj2)
#endif
               YFSij = 2.d0*B0ij + 8.d0*pi_pj * C0ij
            endif
            YFSmat(i,j) = YFSij*etaij            
         enddo
      enddo

! diagonal contributions      
      do i = 1,npartl
         YFSii = 0.d0
         etaii = z(i)*z(i)
         if (abs(etaii).gt.0.d0) then
            mi2 = masses(i)
            YFSii = 4.d0*mi2*C0singular(mi2,phmass) - B0ii(i)
         endif
         YFSmat(i,i)    =  - YFSii*etaii
      enddo
*
      YFSmat = overall*YFSmat
      return
      end      
***********************************************************************************************      
      double complex function C0singular_notassumedlgoesto0(m2,l)
      implicit double precision (a-h,l,m,o-z)
! https://it.wikipedia.org/wiki/Tavola_degli_integrali_indefiniti_di_funzioni_razionali

      p2 = m2

      a  = p2
      b  = -(p2-l*l+m2)
      c  = m2
!      abc = a + b + c
      abc = l*l

      dueapb = p2 + l*l - m2
      
!      disc = 4.d0*a*c - b*b
      disc = l*l * (4.d0*m2 - l*l)
      
!      sq = sqrt(disc)
      sq = l*sqrt(4.d0*m2 - l*l) 

      f1 = 0.5d0/a*log(abs(abc))
      f1 = f1 + (dueapb)/a/sq*atan((2.d0*a+b)/sq)

      f0 = 0.5d0/a*log(abs(c))
      f0 = f0 + (dueapb)/a/sq*atan((b)/sq)

      C0singular_notassumedlgoesto0 = f1-f0

      return
      end
***************************************      
      function soft_integral_mue(egmax,p1,p2,p3,p4)
! this should work also for massless momenta...
      implicit double precision (a-h,l,m,o-z)
      double precision phmass
      common/photonmasslambda/phmass      
      parameter (npart = 4)
      dimension p(0:3),q(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)      
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirst
      common/pmatrix/pmat
      common/parameters/ame,ammu,convfac,alpha,pi
      data ifirst/0/
      
      double precision sffij(4,4)
      common/softFFij/sffij

      integer maskud(4,4),masknoud(4,4)
      common/sffmasks/maskud,masknoud
      common/sffudandnoud/sffud,sffnoud
      
#include "invariants.h"

! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle

      if (ifirst.eq.0) then

         sffdebugIR = 0.d0

         eta(1) = -1.d0 * Qmu     * QRmu
         eta(2) = -1.d0 * (-1.d0) * QRe
         eta(3) =  1.d0 * Qmu     * QRmu
         eta(4) =  1.d0 * (-1.d0) * QRe

         masses(1) = mm*mm
         masses(2) = me*me
         masses(3) = mm*mm
         masses(4) = me*me
         
         ao2pi = alpha/2.d0/pi
         alog2 = log(2.d0)

***************************************
         maskud   = 0
         masknoud = 0

         maskud(1,2) = 1
         maskud(1,4) = 1
         maskud(2,3) = 1
         maskud(3,4) = 1
         maskud(2,1) = 1
         maskud(4,1) = 1
         maskud(3,2) = 1
         maskud(4,3) = 1
         
         masknoud(1,1) = 1
         masknoud(2,2) = 1
         masknoud(3,3) = 1
         masknoud(4,4) = 1

         masknoud(1,3) = 1
         masknoud(3,1) = 1
         masknoud(2,4) = 1
         masknoud(4,2) = 1
***************************************         
         ifirst = 1
      endif
      
      linfra    = log(egmax/phmass)
      linfrapl2 = linfra + alog2

      pmat(1,:) = p1
      pmat(2,:) = p2
      pmat(3,:) = p3
      pmat(4,:) = p4

      npartl = npart

! off-diagonal contributions
      softint = 0.d0
      
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = eta(i) * eta(j)            
            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta_mue(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)
               vl = 0.5d0*(rho*rho*masses(i) - q2)
               v  = vl/(p(0) - q(0))
c               arglog = 1.d0+2.d0*vl/q2
* i.e.    
               arglog=rho*rho*masses(i)/masses(j)
c               if (arglog.gt.0.d0) then 
                 tot=linfrapl2*log(arglog)+termfinite_mue(p,q,v,rho,i,j)
c               else
c                  tot = 0.d0
c               endif
               tot = -tot * 2.d0*dot(p,q)/vl* etaij
               tot =  tot * 2.d0 ! this is the double product when squaring the eikonal
               softint   = softint  + tot

c               sffij(i,j) = tot
c               sffij(j,i) = tot
            endif
         enddo
      enddo      
      soft_integral_mue = softint
      
! diagonal contributions
      softint = 0.d0
      do  i = 1,npartl
         etaii = eta(i)*eta(i)
         if (etaii.gt.1.d-3) then
c            betai=1.d0-masses(i)/pmat(i,0)/pmat(i,0)
c            betai=sqrt(betai)            
            betai = sqrt(tridot(pmat(i,:),pmat(i,:)))/pmat(i,0)
            if (betai.gt.1.d-12) then
               bterm =  log((1.d0-betai)/(1.d0+betai))/betai
            else
               bterm = -2.d0
            endif
            term = 2.d0*linfrapl2  + bterm
            term = -2.d0 * term * etaii
            softint = softint + term
cc              print*, 2.d0*linfrapl2
****            
c            sffij(i,i) = term
         endif
      enddo
! divide by 2 to obtain exactly what is inside {} in Dittmaier & Kraemer
!     paper (hep-ph/0109062)
      
 666  soft_integral_mue = 0.5d0*(soft_integral_mue + softint)*ao2pi

c      sffij = 0.5d0*ao2pi * sffij
c     c      sffdebugIR = 0.5d0*ao2pi*sffdebugIR
      return
      end
***********************************************************************************************
      subroutine rescale_momenta_mue(npart,i,j,pmat,masses,p,q,rho)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3),pmat(npart,0:3),masses(npart)
      dimension p1(0:3),p2(0:3)
      
      p1  = pmat(i,:)
      p2  = pmat(j,:)
      m12 = masses(i)
      m22 = masses(j)
      
      p1p2 = dot(p1,p2)

      rho1 = p1p2 + sqrt(p1p2*p1p2 - m12*m22)
      rho1 = rho1 / m12
c      rho2 = p1p2 - sqrt(p1p2**2 - m12*m22)
c      rho2 = rho2 / m12
!  better numerical solution !
      rho2 = m22/m12/rho1
      if ( (rho1*p1(0)-p2(0)) .gt. 0.d0) then
         rho = rho1
      endif
      if ( (rho2*p1(0)-p2(0)) .gt. 0.d0) then
         rho = rho2
      endif
      p = rho * p1
      q = p2
      return
      end
***********************************************************************************************
      function termfinite_mue(p,q,v,rho,i,j)
      implicit double precision (a-h,m,o-z)
      dimension p(0:3),q(0:3)
      double complex cspen,carg
      parameter (npart=4)
      dimension masses(npart),eta(npart)
      double precision phmass
      common/photonmasslambda/phmass      
      common/sffcommonNEW/ao2pi,masses,eta,alog2,ifirst
      double complex myli2
      external myli2
      
      u0   = p(0)
      umod = sqrt(tridot(p,p))
      pp = u0 + umod
      pm = u0 - umod
      ppopm = pp*pp/rho/rho/masses(i) ! the same as pp/pm

      u0   = q(0)
      umod = sqrt(tridot(q,q))
      qp = u0 + umod
      qm = u0 - umod
      qpoqm = qp*qp/masses(j)  ! the same as qp/qm
      
c      arg1 = (v-pm)/v
c      arg2 = (v-pp)/v
c      arg3 = (v-qm)/v
c      arg4 = (v-qp)/v

      arg1 = 1.d0-pm/v
      arg2 = 1.d0-pp/v
      arg3 = 1.d0-qm/v
      arg4 = 1.d0-qp/v
      
      ddlog1 = dilogmy(arg1)
      ddlog2 = dilogmy(arg2)
      ddlog3 = dilogmy(arg3) !well it works...!
      ddlog4 = dilogmy(arg4)

c      ddlog1 = myli2(dcmplx(arg1)) 
c      ddlog2 = myli2(dcmplx(arg2))
c      ddlog3 = myli2(dcmplx(arg3)) !well it works...!
c      ddlog4 = myli2(dcmplx(arg4))
      
      termfinite_mue =
     .     0.25d0*log(ppopm)**2 + ddlog1 + ddlog2
     .   - 0.25d0*log(qpoqm)**2 - ddlog3 - ddlog4
      return
      end
*****************************************************************      
      double precision FUNCTION DILOGMY(X)
      double precision X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      double precision C(0:18),H,ALFA,B0,B1,B2
      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.289868133696453D0/, PI6 /1.644934066848226D0/
      DATA C( 0) / 0.4299669356081370D0/
      DATA C( 1) / 0.4097598753307711D0/
      DATA C( 2) /-0.0185884366501460D0/
      DATA C( 3) / 0.0014575108406227D0/
      DATA C( 4) /-0.0001430418444234D0/
      DATA C( 5) / 0.0000158841554188D0/
      DATA C( 6) /-0.0000019078495939D0/
      DATA C( 7) / 0.0000002419518085D0/
      DATA C( 8) /-0.0000000319334127D0/
      DATA C( 9) / 0.0000000043454506D0/
      DATA C(10) /-0.0000000006057848D0/
      DATA C(11) / 0.0000000000861210D0/
      DATA C(12) /-0.0000000000124433D0/
      DATA C(13) / 0.0000000000018226D0/
      DATA C(14) /-0.0000000000002701D0/
      DATA C(15) / 0.0000000000000404D0/
      DATA C(16) /-0.0000000000000061D0/
      DATA C(17) / 0.0000000000000009D0/
      DATA C(18) /-0.0000000000000001D0/
      IF(X .EQ. ONE) THEN
         DILOGMY=PI6
         RETURN
      ELSE IF(X .EQ. MONE) THEN
         DILOGMY=MALF*PI6
         RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
         Y=MONE/(ONE+T)
         S=ONE
         A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
         Y=MONE-T
         S=MONE
         A=LOG(-T)
         A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
         Y=(MONE-T)/T
         S=ONE
         A=LOG(-T)
         A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
         Y=-T/(ONE+T)
         S=MONE
         A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
         Y=T
         S=ONE
         A=ZERO
      ELSE
         Y=ONE/T
         S=MONE
         A=PI6+HALF*LOG(T)**2
      END IF
      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO  I = 18,0,-1
         B0=C(I)+ALFA*B1-B2
         B2=B1
 1       B1=B0
      enddo
      DILOGMY=-(S*(B0-H*B2)+A)
c      print*,dilogmy, li2(x)
      RETURN
      END
***********************************************************************************************
      function softcollint_mue(egmax,p1,p2,p3,p4,ionoff)
! this should work also for massless momenta...
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      dimension p(0:3),q(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)
      double precision scales(4,4)
      common/softcollintcmn/ao2pi,masses,eta,scales,ifirst
      data ifirst/0/
      common/parameters/ame,ammu,convfac,alpha,pi

      integer ionoff(4,4)
      
#include "invariants.h"      

! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle

      if (ifirst.eq.0) then

         scales = 0.d0
         
         eta(1) = -1.d0 * Qmu     * QRmu
         eta(2) = -1.d0 * (-1.d0) * QRe
         eta(3) =  1.d0 * Qmu     * QRmu
         eta(4) =  1.d0 * (-1.d0) * QRe

         masses(1) = mm*mm
         masses(2) = me*me
         masses(3) = mm*mm
         masses(4) = me*me
         
         ao2pi = alpha/2.d0/pi
         
         ifirst = 1
      endif

      pmat(1,:) = p1
      pmat(2,:) = p2
      pmat(3,:) = p3
      pmat(4,:) = p4

c      scales(1,1) = p1(0)
c      scales(2,2) = p2(0)
c      scales(3,3) = p3(0)
c      scales(4,4) = p4(0)

c      scales(1,2) = sqrt(abs(dot(p1+p2,p1+p2)))
c      scales(1,3) = sqrt(abs(dot(p1-p3,p1-p3)))
c      scales(1,4) = sqrt(abs(dot(p1-p4,p1-p4)))
c      scales(2,3) = sqrt(abs(dot(p2-p3,p2-p3)))
c      scales(2,4) = sqrt(abs(dot(p2-p4,p2-p4)))
c      scales(3,4) = sqrt(abs(dot(p3+p4,p3+p4)))

c      scales(1,3) = scales(1,2)
c      scales(1,4) = scales(1,2)
c      scales(2,3) = scales(1,2)
c      scales(2,4) = scales(1,2)
c      scales(3,4) = scales(1,2)

cc ! this seeems to mimimize the difference, for mu-      
c      scales(1,1) = sqrt(2.d0*dot(p1,p2))*0.5d0      
      scales(1,1) = sqrt(abs(dot(p1+p2,p1+p2)))*0.5d0 ! slightly better...

      scales(1,2) = scales(1,1)
      scales(1,3) = scales(1,1)
      scales(1,4) = scales(1,1)
      scales(2,3) = scales(1,1)
      scales(2,4) = scales(1,1)
      scales(3,4) = scales(1,1)

      scales(2,2) = scales(1,1)
      scales(3,3) = scales(1,1)
      scales(4,4) = scales(1,1)      
      
      npartl = npart

! off-diagonal contributions
      softint = 0.d0
      do i = 1,npartl-1
         do j = i+1,npartl
            etaij = eta(i) * eta(j) * ionoff(i,j)
            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta_mue(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)
               vl = 0.5d0*(rho*rho*masses(i) - q2)
               v  = vl/(p(0) - q(0)) 

               arglog = 1.d0+2.d0*vl/q2

               terminfra = 0.d0              
               if (arglog.gt.0.d0) then 
                  terminfra = log(arglog)
                  tot = terminfra !+ termfinite_mue(p,q,v,rho,i,j)
               else
                  tot = 0.d0
               endif
               
               tot = -tot * 2.d0*dot(p,q)/vl* etaij
               tot =  tot * 2.d0 ! this is the double product when
                                 ! squaring the eikonal
               tot = tot*(-2.d0*log(egmax/scales(i,j))-1.5d0)
               softint  = softint + tot
            endif
         enddo
      enddo      
      softcollint_mue = softint

! diagonal contributions

      const = - 1.5d0
      
      softint = -4.d0 * (
     .     eta(1)**2*ionoff(1,1)*(-2.d0*log(egmax/scales(1,1))+const)+
     .     eta(2)**2*ionoff(2,2)*(-2.d0*log(egmax/scales(2,2))+const)+
     .     eta(3)**2*ionoff(3,3)*(-2.d0*log(egmax/scales(3,3))+const)+
     .     eta(4)**2*ionoff(4,4)*(-2.d0*log(egmax/scales(4,4))+const)
     .     )
      softcollint_mue = softcollint_mue + softint
      softcollint_mue = softcollint_mue*0.125d0 ! = 1/8

      return
      end
***********************************************************************************************
      subroutine get_eps(k,e1,e2)
      implicit double precision (a-h,o-z)
      double precision k(0:3),e1(0:3),e2(0:3),ku(0:3)
      double precision e0(0:3),e3(0:3)

      ku = k/k(0)

      e1(0) = 0.d0
      if (abs(ku(3)).gt.1d-12) then
         e1(1) = ku(1)
         e1(2) = ku(2)
         e1(3) = (-e1(1)*e1(1) - e1(2)*e1(2))/ku(3)
      elseif (abs(ku(1)).gt.1d-12) then
         e1(3) = ku(3)
         e1(2) = ku(2)
         e1(1) = (-e1(3)*e1(3) - e1(2)*e1(2))/ku(1)
      else
         e1(1) = ku(1)
         e1(3) = ku(3)
         e1(2) = (-e1(3)*e1(3) - e1(1)*e1(1))/ku(2)
      endif
      e1 = e1/(sqrt(tridot(e1,e1)))

      call vecprod(e1,ku,e2)
      e2(0) = 0.d0
      return
      end
************************************
      subroutine vecprod(p1,p2,p3)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3)
      p3(0) =  p1(0)+p2(0)
      p3(1) =  p1(2)*p2(3) - p1(3)*p2(2)
      p3(2) = -p1(1)*p2(3) + p1(3)*p2(1)
      p3(3) =  p1(1)*p2(2) - p1(2)*p2(1)
      return
      end
****************************************************

***** from TEUBNER et al. routine
      double precision function alf_lep_carlo(gen_lep, s)
      implicit NONE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculates Leptonic Contribution \Delta \alpha_lep(s) 
c at LO, NLO and NNLO
c
c   Pi0 := ( Pi(qsq) - Pi(0) )_LO ,                   
c   Pi1 := ( Pi(qsq) - Pi(0) )_NLO ,                  
c
c (Taking real part of Pi is assumed.)
c
c LO and NLO results are checked against the results given in 
c Steinhauser, hep-ph/9803313 
c
c Test at s = Mz^2 using PDG97 values for a comparison with 
c Table 1 of Steinhauser, hep-ph/9803313.  Perfect agreement!
c
c          e             mu            tau
c LO  0.0174346533  0.00917841908  0.00480593422  *1
c LO  0.0174346533  0.00917841908  0.00480592861  *2
c LO  0.0174346531  0.00917841908  0.00480592861  (full formula)
c NLO 0.0000379829  0.00002359991  0.00001603387  (full formula)
c NNLO0.0000003821  0.00000036614  0.00000031449  *3
c
c *1 An expansion up to and including M^2/s terms is used.
c *2 An expansion up to and including M^4/s^2 terms is used.
c *3 The expansion given in hep-ph/9803313 (up to and including 
c    M^2/s terms) is used.
c
c Aug 2009: expression for Pi_1 at q^2/m^2 << 1 improved.
c Mar 2010: CERNLIB no longer necessary to evaluate ReF and ImG
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c input ----------------------------------------------------
      integer gen_lep
      double precision s
c const ----------------------------------------------------
      double precision Me, Mmu, Mtau, PI, ZETA2, ZETA3, ZETA5
c      Parameter(Me=0.510998928d-3, Mmu=0.1056583715d0, Mtau=1.77682d0)! PDG14
c      Parameter( ALPHA = 1.d0/137.035999074d0 ) ! PDG14
      Parameter( PI = 3.141592653589793238d0, ZETA2 = PI ** 2 / 6.d0 )
      Parameter( ZETA3 = 1.2020569032d0,      ZETA5 = 1.0369277551d0 )
c function, tmp var. ---------------------------------------
      double precision Pi0, Pi1, eta, mlep
      double precision PiA, Piell, PiF, Pih, eta1, eta2, stepf, y
      double precision alf_lepLO, alf_lepNLO, alf_lepNNLO

*******CARLO
      double precision ame,ammu,convfac,alpha,piii
      common/parameters/ame,ammu,convfac,alpha,piii
      double precision  AMASSES(9)
      common/massesVP/amasses
      double precision dalo(3),danlo(3),dannlo(3)
      common/dalonlonnlo/dalo,danlo,dannlo
      
c begin ====================================================
c define functions -----------------------------------------
c NNLO -----------------------------------------------------      
      stepf( y ) = ( 1.d0 + abs( y ) / y ) / 2.d0

      PiA( eta ) = 1.d0/(16.d0 * PI **2) * 
     &   ( -121.d0/48.d0 + ( -5.d0 + 8.d0 * Log(2.d0) ) * ZETA2 
     &     - 99.d0/16.d0 * ZETA3 + 10.d0 * ZETA5
     &     +  1.d0/ 8.d0 * Log( abs( 1.d0/eta ) ) )

      Piell( eta1, eta2 ) = 1.d0/(16.d0 * PI**2) * 
     & ( -116.d0/27.d0 + 4.d0/3.d0 * ZETA2 + 38.d0/9.d0 * ZETA3
     &  +  14.d0/9.d0 * Log( abs(1.d0/eta1) ) 
     &  + ( 5.d0/18.d0 - 4.d0/3.d0 * ZETA3 ) * Log( abs(1.d0/eta2) )  
     &  +   1.d0/6.d0 * ( Log( abs(1.d0/eta1) ) **2 
     &                    - PI ** 2 * stepf( - eta1 ) )
     &  -   1.d0/3.d0 * ( Log( abs(1.d0/eta1) ) * Log( abs(1.d0/eta2) )  
     &                    - PI ** 2 * stepf( - eta1 ) ) ) 

      PiF( eta ) = 1.d0/(16.d0 * PI **2) * 
     &  ( -307.d0/216.d0 - 8.d0/3 * ZETA2 + 545.d0/144.d0 * ZETA3
     &  + ( 11.d0/6.d0 - 4.d0/3.d0 * ZETA3 ) * Log( abs(1.d0/eta) ) 
     &  -    1.d0/6.d0 * ( Log( abs(1.d0/eta) ) **2 
     &                     - PI **2 * stepf( - eta ) ) )

      Pih( eta ) = 1.d0/(16.d0 * PI **2) * 
     & ( -37.d0/6.d0 + 38.d0/9.d0 * ZETA3 
     &  +  ( 11.d0/6.d0 - 4.d0/3.d0 * ZETA3 ) * Log( abs(1.d0/eta) )  
     &  -  1.d0/6.d0 * ( Log( abs(1.d0/eta) ) **2 
     &                   - PI **2 * stepf( - eta ) ) ) 

c calc contributions to \Delta \alpha ----------------------
      Me   = amasses(1)
      Mmu  = amasses(2)
      Mtau = amasses(3)

      if ( gen_lep .eq. 1 ) mlep = Me 
      if ( gen_lep .eq. 2 ) mlep = Mmu 
      if ( gen_lep .eq. 3 ) mlep = Mtau 

      alf_lepLO = - 4.d0 * PI * ALPHA * Pi0( - mlep **2 / s ) 

      alf_lepNLO = - 4.d0 * PI * ALPHA * ( ALPHA / PI ) *
     &             Pi1( - mlep **2 / s ) 

      if ( gen_lep .eq. 1 ) then
         alf_lepNNLO = - 4.d0 * PI * ALPHA * ( ALPHA / PI )**2 * ( 
     &         PiA( - me**2/ s ) 
     &      +  PiF( - me**2/ s ) 
     &      +  Pih( - mmu**2 / s ) 
     &      +  Pih( - mtau**2/ s ) 
     &         )
      else if ( gen_lep .eq. 2 ) then
         alf_lepNNLO = - 4.d0 * PI * ALPHA * ( ALPHA / PI )**2 * ( 
     &        PiA( - mmu**2/ s ) 
     &      + Piell( - mmu**2/ s, - me**2/ s ) 
     &      + PiF( - mmu**2/ s ) 
     &      + Pih( - mtau**2/ s ) 
     &        )
      else if ( gen_lep .eq. 3 ) then
         alf_lepNNLO = - 4.d0 * PI * ALPHA * ( ALPHA / PI )**2 * ( 
     &        PiA( - mtau**2/ s ) 
     &      + Piell( - mtau**2/ s, - me**2 / s ) 
     &      + Piell( - mtau**2/ s, - mmu**2/ s ) 
     &      + PiF( - mtau**2/ s ) 
     &        ) 
      end if

      dalo(gen_lep)   = alf_lepLO
      danlo(gen_lep)  = alf_lepNLO
      dannlo(gen_lep) = alf_lepNNLO
c default
      alf_lep_carlo = 1.d0*alf_lepLO+1.d0*alf_lepNLO+1.d0*alf_lepNNLO
c default, end
      return
      end

***** from TEUBNER et al. routine, only nlo blob
      double precision function alf_lep_carlo_nlo(gen_lep, s)
      implicit NONE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculates Leptonic Contribution \Delta \alpha_lep(s) 
c at LO, NLO and NNLO
c
c   Pi0 := ( Pi(qsq) - Pi(0) )_LO ,                   
c   Pi1 := ( Pi(qsq) - Pi(0) )_NLO ,                  
c
c (Taking real part of Pi is assumed.)
c
c LO and NLO results are checked against the results given in 
c Steinhauser, hep-ph/9803313 
c
c Test at s = Mz^2 using PDG97 values for a comparison with 
c Table 1 of Steinhauser, hep-ph/9803313.  Perfect agreement!
c
c          e             mu            tau
c LO  0.0174346533  0.00917841908  0.00480593422  *1
c LO  0.0174346533  0.00917841908  0.00480592861  *2
c LO  0.0174346531  0.00917841908  0.00480592861  (full formula)
c NLO 0.0000379829  0.00002359991  0.00001603387  (full formula)
c NNLO0.0000003821  0.00000036614  0.00000031449  *3
c
c *1 An expansion up to and including M^2/s terms is used.
c *2 An expansion up to and including M^4/s^2 terms is used.
c *3 The expansion given in hep-ph/9803313 (up to and including 
c    M^2/s terms) is used.
c
c Aug 2009: expression for Pi_1 at q^2/m^2 << 1 improved.
c Mar 2010: CERNLIB no longer necessary to evaluate ReF and ImG
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c input ----------------------------------------------------
      integer gen_lep
      double precision s
c const ----------------------------------------------------
      double precision Me, Mmu, Mtau, PI, ZETA2, ZETA3, ZETA5
c      Parameter(Me=0.510998928d-3, Mmu=0.1056583715d0, Mtau=1.77682d0)! PDG14
c      Parameter( ALPHA = 1.d0/137.035999074d0 ) ! PDG14
      Parameter( PI = 3.141592653589793238d0, ZETA2 = PI ** 2 / 6.d0 )
      Parameter( ZETA3 = 1.2020569032d0,      ZETA5 = 1.0369277551d0 )
c function, tmp var. ---------------------------------------
      double precision Pi0, Pi1, eta, mlep
      double precision PiA, Piell, PiF, Pih, eta1, eta2, stepf, y
      double precision alf_lepLO, alf_lepNLO, alf_lepNNLO

*******CARLO
      double precision ame,ammu,convfac,alpha,piii
      common/parameters/ame,ammu,convfac,alpha,piii
      double precision  AMASSES(9)
      common/massesVP/amasses
      double precision dalo(3),danlo(3),dannlo(3)
      common/dalonlonnlo/dalo,danlo,dannlo
      
c begin ====================================================
c define functions -----------------------------------------

c calc contributions to \Delta \alpha ----------------------
      mlep = amasses(gen_lep)

      alf_lep_carlo_nlo = - 4.d0 * PI * ALPHA * ( ALPHA / PI ) *
     &     Pi1( - mlep **2 / s ) 
      return
      end
ctt
c Routine for \Delta\alpha_{top}
ctt
      double precision function alf_top_carlo( s )
      implicit NONE
      double precision  AMASSES(9)
      common/massesVP/amasses
      double precision ame,ammu,convfac,alpha,piii
      common/parameters/ame,ammu,convfac,alpha,piii
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Calculates Top Contribution \Delta \alpha_top(s) at LO and NLO
c
c   Pi0 := ( Pi(qsq) - Pi(0) )_LO ,                   
c   Pi1 := ( Pi(qsq) - Pi(0) )_NLO ,                  
c
c (Taking real part of Pi is assumed.)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c input ----------------------------------------------------
      double precision s
c const ----------------------------------------------------
      double precision Mtop, ALPHAS, PI, Nc, Qt
cc      Parameter( Mtop = 173.21d0 )              ! PDG14
cc      Parameter( ALPHA = 1.d0/137.035999074d0 ) ! PDG14
      Parameter( ALPHAS = 0.1185d0 )            ! PDG14
      Parameter( PI = 3.141592653589793238d0, Nc = 3.d0, Qt = 2.d0/3.d0)
c function, tmp var. ---------------------------------------
      double precision Pi0t, Pi1t, ss, alf_topLO, alf_topNLO

      double precision dalotop,danlotop
      common/dalonlotop/dalotop,danlotop
c begin ====================================================
c LO -------------------------------------------------------
      Pi0t( ss ) = Nc * Qt ** 2 * (
     &    4.d0/  15.d0 *  ss / Mtop **2 
     &  + 1.d0/  35.d0 * ( ss / Mtop **2 ) **2
     &  + 1.d0/ 945.d0 * ( ss / Mtop **2 ) **3 
     &  + 1.d0/5544.d0 * ( ss / Mtop **2 ) **4 
     &  + 1.d0/30030.d0 * ( ss / Mtop **2 ) **5 
     &  + 1.d0/154440.d0 * ( ss/ Mtop **2 ) **6 
     &  + 1.d0/765765.d0 * ( ss/ Mtop **2 ) **7 
     &   ) / ( 16.d0 * PI **2 )

c NLO ------------------------------------------------------

      Pi1t( ss ) =  
     & (        82.d0/81.d0 * ss / Mtop **2  
     &     +   449.d0/2700.d0 * ( ss / Mtop **2 ) **2
     &     + 62479.d0/1984500.d0 * ( ss / Mtop **2 ) **3 
     &  ) / ( 16.d0 * PI **2 )

c calc contributions to \Delta \alpha ----------------------
      Mtop = amasses(4)

      alf_topLO = - 4.d0 * PI * ALPHA * Pi0t( s ) 

      alf_topNLO = - 4.d0 * PI * ALPHA * ( Nc**2 - 1.d0 ) / 2.d0
     &             * Qt**2 * ( ALPHAS / PI ) * Pi1t( s ) 

      dalotop  = alf_topLO
      danlotop = alf_topNLO
      
      alf_top_carlo = 1.d0*alf_topLO + 1.d0 * alf_topNLO

      return
      end
      

**************************************************
      function SigAA(q2) ! from Hto4l
*** B.1
      implicit none
      double complex   SigAA,SigAAtmp
      double precision q2,cq2

      double precision Qd2,Qu2,czero
      integer k
      integer ionlytop,omt,lon
      common/onlytop/ionlytop,omt

      double precision m2gend(3),m2genu(3)
      double precision mmrcl,mercl,murcl,mdrcl,msrcl,mcrcl,mbrcl,mtrcl,
     .     mtaurcl
      common/fermassesrcl/mmrcl,mercl,mtaurcl,murcl,mdrcl,msrcl,mcrcl,
     .    mbrcl,mtrcl

      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif

      cq2   = q2
      SigAA = (0.d0,0.d0)

      czero = 0.d0
      
      omt = 1
** UV test
c      goto 111

** leptonic contribution
      SigAAtmp  = (0.d0,0.d0)
      m2gend(1) = mercl*mercl
      m2gend(2) = mmrcl*mmrcl
      m2gend(3) = mtaurcl*mtaurcl
      do k = 1,3
         SigAAtmp = SigAAtmp + q2/3.d0
     .        -(q2   + omt * m2gend(k)*2.d0) * B0(uquad*cq2,
     .        uquad*m2gend(k),uquad*m2gend(k))
     .        + 2.d0 * omt * m2gend(k)* B0(uquad*czero,uquad*m2gend(k),
     .        uquad*m2gend(k))
      enddo
      SigAAtmp = SigAAtmp*4.d0/3.d0
      SigAA    = SigAA+SigAAtmp
      
*** quark contribution
 111  m2genu(1) = murcl*murcl
      m2genu(2) = mcrcl*mcrcl
      m2genu(3) = mtrcl*mtrcl
      m2gend(1) = mdrcl*mdrcl
      m2gend(2) = msrcl*msrcl
      m2gend(3) = mbrcl*mbrcl
      Qd2       = 1.d0/9.d0 ! down charge squared
      Qu2       = 4.d0/9.d0 ! up charge squared

** down quark contribution
      SigAAtmp   = (0.d0,0.d0)
      do k = 1,2
         SigAAtmp = SigAAtmp + q2/3.d0
     .        -(q2 + omt * m2gend(k)*2.d0)*B0(uquad*cq2,uquad*m2gend(k),
     .        uquad*m2gend(k))
     .        + 2.d0 * omt * m2gend(k)*B0(uquad*czero,uquad*m2gend(k),
     .        uquad*m2gend(k))
      enddo
** bottom
      k = 3
      SigAAtmp = SigAAtmp + q2/3.d0
     .     -  (q2 + m2gend(k)*2.d0)*B0(uquad*cq2,uquad*m2gend(k),
     .     uquad*m2gend(k))
     .     + 2.d0 * m2gend(k)      *B0(uquad*czero,uquad*m2gend(k),
     .     uquad*m2gend(k))
*
      SigAAtmp = Qd2 * SigAAtmp * 4.d0  ! that's 2*2 ! Nc = 3, it cancels 1/3
      SigAA    = SigAA+SigAAtmp !* (1.d0 - ionlytop)
      
** up quark contribution
      SigAAtmp   = (0.d0,0.d0)
      do k = 1,2
         SigAAtmp = SigAAtmp + q2/3.d0
     .        -(q2+ omt * m2genu(k)*2.d0)*B0(uquad*cq2,uquad*m2genu(k),
     .        uquad*m2genu(k))
     .        + 2.d0 * omt * m2genu(k) * B0(uquad*czero,uquad*m2genu(k),
     .        uquad*m2genu(k))
      enddo
* TOP
      k = 3
      SigAAtmp = SigAAtmp + q2/3.d0
     .     -(q2+ m2genu(k)*2.d0) * B0(uquad*cq2,uquad*m2genu(k),
     .     uquad*m2genu(k))
     .     + 2.d0 * m2genu(k)    * B0(uquad*czero,uquad*m2genu(k),
     .     uquad*m2genu(k))
**

      SigAAtmp = Qu2 * SigAAtmp * 4.d0 ! that's 2*2 ! Nc = 3, it cancels 1/3
      SigAA    = SigAA+SigAAtmp

*** rescaling by -alpha/4.d0/pi
      SigAA = -SigAA*alpha*0.25d0/pi

      return
      end

**************************************************
      function dSigAA(q2) ! From Hto4l
*** derivative of B.1
      implicit none
      integer k
      double precision dSigAAtmp
      double complex dSigAA
      double precision mh2,q2,cq2
      double precision Qd2,Qu2
      double precision gfm,gfp,czero

      double precision m2gend(3),m2genu(3)

      integer ionlytop,omt
      common/onlytop/ionlytop,omt

      double precision mmrcl,mercl,murcl,mdrcl,msrcl,mcrcl,mbrcl,mtrcl,
     .     mtaurcl
      common/fermassesrcl/mmrcl,mercl,mtaurcl,murcl,mdrcl,msrcl,mcrcl,
     .    mbrcl,mtrcl

      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi

#include "looptools.h"
      RealType uquad
#ifdef QUAD
      parameter (uquad=1.q0)
#else
      parameter (uquad=1.d0)
#endif
      
      cq2    = q2
      dSigAA = (0.d0,0.d0)
      czero = 0.d0
      omt = 1
** UV test
c      goto 111

** leptonic contribution
      dSigAAtmp = (0.d0,0.d0)
      m2gend(1) = mercl*mercl
      m2gend(2) = mmrcl*mmrcl
      m2gend(3) = mtaurcl*mtaurcl
      do k = 1,3
         dSigAAtmp = dSigAAtmp + 1.d0/3.d0
     .        -(q2+ omt * m2gend(k)*2.d0)*DB0(uquad*cq2,uquad*m2gend(k),
     .        uquad*m2gend(k))
     .     -(1.d0)* B0(uquad*cq2,uquad*m2gend(k),uquad*m2gend(k))
      enddo
      dSigAAtmp = 2.d0 * dSigAAtmp*2.d0/3.d0
      dSigAA    = dSigAA+dSigAAtmp
c      goto 111

*** quark contribution
      m2genu(1) = murcl*murcl
      m2genu(2) = mcrcl*mcrcl
      m2genu(3) = mtrcl*mtrcl
      m2gend(1) = mdrcl*mdrcl
      m2gend(2) = msrcl*msrcl
      m2gend(3) = mbrcl*mbrcl
      Qd2       = 1.d0/9.d0 ! down charge squared
      Qu2       = 4.d0/9.d0 ! up charge squared

** down quark contribution
      dSigAAtmp   = (0.d0,0.d0)
      do k = 1,2
         dSigAAtmp = dSigAAtmp + 1.d0/3.d0
     .        -(q2+ omt * m2gend(k)*2.d0) *
     .        DB0(uquad*cq2,uquad*m2gend(k),uquad*m2gend(k))
     .     -(1.d0)* B0(uquad*cq2,uquad*m2gend(k),uquad*m2gend(k))
      enddo
* bottom 
      k = 3
      dSigAAtmp = dSigAAtmp + 1.d0/3.d0
     .     -(q2+ m2gend(k)*2.d0) * DB0(uquad*cq2,
     .     uquad*m2gend(k),uquad*m2gend(k))
     .     -(1.d0)* B0(uquad*cq2,uquad*m2gend(k),uquad*m2gend(k))
**
      dSigAAtmp = Qd2 * dSigAAtmp * 4.d0  ! that's 2*2 ! Nc = 3, it cancels 1/3
      dSigAA    = dSigAA+dSigAAtmp
** up quark contribution
      dSigAAtmp   = (0.d0,0.d0)
      do k = 1,2
         dSigAAtmp = dSigAAtmp + 1.d0/3.d0
     .        -(q2+ omt * m2genu(k)*2.d0) * DB0(uquad*cq2,
     .        uquad*m2genu(k),uquad*m2genu(k))
     .     -(1.d0) * B0(uquad*cq2,uquad*m2genu(k),uquad*m2genu(k))
      enddo

** TOP 
      k = 3
      dSigAAtmp = dSigAAtmp + 1.d0/3.d0
     .     -(q2 + m2genu(k)*2.d0) * DB0(uquad*cq2,
     .     uquad*m2genu(k),uquad*m2genu(k))
     .     -(1.d0) * B0(uquad*cq2,uquad*m2genu(k),uquad*m2genu(k))
**
      dSigAAtmp = Qu2 * dSigAAtmp * 4.d0 ! that's 2*2 ! Nc = 3, it cancels 1/3
      dSigAA    = dSigAA+dSigAAtmp

*** rescaling by -alpha/4.d0/pi
 111  dSigAA = -dSigAA*alpha*0.25d0/pi

      return
      end
**************************************************
      double precision function cts1g1L(p1,p2,p3,p4,k)
c      implicit none
      implicit double precision (Z) ! posso metterli reali perche' ho raccolto in form tutto cio' che e' imm...
#include "invariants.h"
      double precision t13m1,t24m1,p1km1,p2km1,p3km1,p4km1
c      double precision muemueg1L ! this is the variable in form
      double complex muemueg1L ! this is the variable in form
      double precision k(0:3)
      double precision cover,el8
*      double precision p1k,p2k,p3k,p4k
      common/ifirstcts1g1L/cover,el8,DELTA,ifirst
      double complex vpolc
      external vpolc
      integer ifirst,iii
      data ifirst /0/

      double precision DELTA

      double precision p1_eps, p2_eps, p3_eps, p4_eps, k_eps
      double precision p1_epss,p2_epss,p3_epss,p4_epss,k_epss
      common/gaugeinv/p1_eps, p2_eps, p3_eps, p4_eps, k_eps,
     .     p1_epss,p2_epss,p3_epss,p4_epss,k_epss

      double complex Sigrlii,dSigrlii,dZiie,dZiimu,dSigAA
      double complex SigSii,dSigSii,dZAA
      common/ctscommon/dmasse,dmassmu,dZiie,dZiimu,dZAA
      double precision ctstmp
      integer ifirstauuauuau
      common/qedcts/ctstmp,ifirstauuauuau

      if (ifirst.eq.0) then
         DELTA  = getdelta()
         el8    = ec4*ec4
         cover  = el8*pi216m1
         ifirst = 1

         dmasse  = 0.5d0*me*(2.d0*Sigrlii(me2,me2)+2.d0*SigSii(me2,me2))   
         dmassmu = 0.5d0*mm*(2.d0*Sigrlii(mm2,mm2)+2.d0*SigSii(mm2,mm2))   

         dZiie = -Sigrlii(me2,me2)
     .        - 2.d0*me2*(dSigrlii(me2,me2)+dSigSii(me2,me2))
         dZiimu = -Sigrlii(mm2,mm2)
     .        - 2.d0*mm2*(dSigrlii(mm2,mm2)+dSigSii(mm2,mm2))

** spegnendo dmass[e,mu] ottengo i cts esterni come li avevo fatti naively
c         dmasse = 0.d0
c         dmassmu = 0.d0

! Questi cts corrispondono esattamente a quelli gia' implementati
         dZAA = -dSigAA(0.d0)         
c        dZAA = 0.d0
         
c         print*,'dSigAA_0',dSigAA(0.d0)
c         stop         
c         print*,ctstmp,2.d0*(dZiie+dZiimu)
c         stop
      endif
      
c      call loadinvariants1g(p1,p2,p3,p4,k) ! if first the tree level is called, this can be avoided
c      call loadQEDfuns1g1L ! loads all needed LT functions

      t13m1 = 1.d0/t13
      t24m1 = 1.d0/t24
      
      p1km1 = 1.d0/p1k
      p2km1 = 1.d0/p2k
      p3km1 = 1.d0/p3k
      p4km1 = 1.d0/p4k

c      print*,t13,t24,s12,s34,u14,u23,p1k,p2k,p3k,
c     -     p4k,p1mk2,p2mk2,p3pk2,p4pk2
      

c      print*,'Reenable includes in muemue1g1Lnoud.F!!! stopping'
c      stop
      include 'cts1g1L.f'

c      print*,' '
c      print*,el8,cover
c      print*,muemueg1L
c      print*,t13m1,t24m1
c      print*,p1km1,p2km1,p3km1,p4km1

      
      cts1g1L = cts1g1L*ec4*ec2

      return
      end
****************************************************
