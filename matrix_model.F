      subroutine squared_matrix(model,ng,npairs,ecms,p1,p2,p3,p4,
     >     pbeam1,pbeam2,qph,ie,imtx,summt2)
      implicit double precision (a-h,o-z)
      character*10 model
      character*6 eikng ! multip single
      parameter (imaxph = 40)
      dimension iclose(imaxph)
      double precision csi(1)
      double precision mass,lambda
      dimension p1(0:3),p2(0:3),pin1(0:3),pin2(0:3),qmost(0:3)
      dimension p3(0:3),p4(0:3)
      dimension p1r(0:3),p2r(0:3),p1b(0:3),p2b(0:3),q(0:3)
      dimension q1(0:3),q2(0:3),q3(0:3),q4(0:3),qq(0:3),qr(0:3)
      dimension p1tmp(0:3),p2tmp(0:3),p3tmp(0:3),p4tmp(0:3)
      dimension pin1bl(0:3),pin2bl(0:3),p1bl(0:3),p2bl(0:3)
      dimension ptmp(0:3),ptmp2(0:3)
      dimension pbeam1(0:3),pbeam2(0:3)
      dimension qmostpt(imaxph,0:3)
      dimension p1rr(0:3),p2rr(0:3),pin1rr(0:3),pin2rr(0:3),prif(0:3)
      dimension qph(imaxph,0:3),qphtmp(imaxph,0:3),randvec(3)
      character*2 fs
      common/finalstate/fs
      common/iclosest/iclose
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3),
     >     pin2r(0:3)
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeq2,bornmez,
     >     bornmenovp 
      common/radpattern/nph(4)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/epssoft/eps
      common/fortests/prod
      common/for_debug/ppp(0:3),denphsp,dir3(0:3),dir4(0:3),br,bq,onm
      common/singlefactors/singleeikfact(imaxph)
      double complex vpolc
      external vpolc
      double precision muemuegg
      external muemuegg
      double complex epsfun
      external epsfun
!     imtx is to say that the matrix el. has been calculated..
! summt2 is the squared ME summed over all all spins
      
      common/firstsqmatrix/icount,ifirst
      data icount,ifirst /0,0/

      integer irecola
      common/recolamode/irecola

      common/photoncutoff/egmin,egmin2

      double precision muemue1g1Lnoupdown,muemue1g1Lupdown
      double precision muemue1g1Lupdown_ten
      external muemue1g1Lnoupdown,muemue1g1Lupdown,muemue1g1Lupdown_ten
      double precision udno,udyes,udyesten
      
      character*6 ord
      common/qedORDER/ord
      
***   for cfr recola
      double precision dalp24,dalp13,radmuefermionloop
      common/dalphasfor1g1Lcfrrecola/dalp24,dalp13,radmuefermionloop
* can be commented ************************************      
      double precision sffdebugIR(4,4),lt(13),cll(13)
      common/softFFij/sffdebugIR
      integer Amp(24)
      common/switchampls/Amp! for debugging IR...
      common/invariants/p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,tm1
      common/invariants1g/p1k,p2k,p3k,p4k,s12,s34,t24,t13,u14,u23,
     .     p1mk2,p2mk2,p3pk2,p4pk2
      double precision epsk(0:3),epsk1(0:3),epsk2(0:3)
      double precision q5(0:3),q6(0:3)
      double precision p1_eps,p2_eps,p3_eps,p4_eps,k_eps
      double precision p1_epss,p2_epss,p3_epss,p4_epss,k_epss
      common/gaugeinv/p1_eps, p2_eps, p3_eps, p4_eps, k_eps,
     .     p1_epss,p2_epss,p3_epss,p4_epss,k_epss
* end can be commented ************************************      
      double precision facIS,facFS
      integer iaxial,itnir1,itir1,itnir2,itir2,itnir3,itir3,itnir4,itir4
      integer ib1dir,ib1dnir,ib2dir,ib2dnir
      integer ib3dir,ib3dnir,ib4dir,ib4dnir
      integer ib1cir,ib1cnir,ib2cir,ib2cnir
      integer ib3cir,ib3cnir,ib4cir,ib4cnir
      integer ip1dir,ip1dnir,ip2dir,ip2dnir
      integer ip1cir,ip1cnir,ip2cir,ip2cnir            
      common/IRAXFACTORStest/iaxial,itnir1,itir1,itnir2,
     .     itir2,itnir3,itir3,itnir4,itir4,
     .     ib1dir,ib1dnir,ib2dir,ib2dnir,
     .     ib3dir,ib3dnir,ib4dir,ib4dnir,
     .     ib1cir,ib1cnir,ib2cir,ib2cnir,
     .     ib3cir,ib3cnir,ib4cir,ib4cnir,
     .     ip1dir,ip1dnir,ip2dir,ip2dnir,
     .     ip1cir,ip1cnir,ip2cir,ip2cnir

      double complex rclscale
      common/rescalecouplingrecola/rclscale
      
**** * for crossing symmetry
      double precision Qmu
      common/muoncharge/Qmu
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu      
****
***   for reweighenting
      integer nextraweights
      parameter (nextraweights=4)
      double precision weightdefault,extraweights(0:nextraweights)
      double precision reweightLO,reweightNLO,wnorm
      common/weights/reweightLO,reweightNLO,weightdefault,extraweights
**************************
      double precision p1c(0:3),p2c(0:3),p3c(0:3),p4c(0:3),k1c(0:3),
     .     k2c(0:3)
      common/momentaforloopcalculation/p1c,p2c,p3c,p4c,k1c,k2c
**************************
      double precision dah13,dal13,dah24,dal24
      common/vp1gand1g1L/vp13,vp24,dah13,dal13,dah24,dal24
! coefficients for VP reweighting, from  invariants.h
      double precision cwvp(0:10),cwvpl(0:10)
      common/coeffforvpreweighting/cwvp,cwvpl
****
c      integer xx(4)
c      xx    = 0
c      xx(1) = 1
c      xx(2) = 1
c      do i = 1,4
c         do j=i,4
c            print*,i,j,xx(i)*xx(j)
c         enddo
c      enddo
c      stop
      if (ifirst.eq.0) then
         iaxial = 1 ! actually not used anymore
         vp13   = 1.d0
         vp24   = 1.d0
         dah13  = 0.d0
         dal13  = 0.d0
         dah24  = 0.d0
         dal24  = 0.d0
         ifirst = 1
      endif
      
      cwvp     = 0.d0
      imtx     = 0
      summt2   = 0.d0
      bornme   = 0.d0 !1.d-13
      bornmeq2 = 0.d0 !1.d-13
      bornmez  = 0.d0 !1.d-13
      
      if (ie.gt.0) return
      
**** needed in the case I want to use 1-loop tensors instead of 1-loop tensors coefficients,
**** which now it's the default for up-down boxes and pentagons      
**** the common momentaforloopcalculation is brought around in invariants.h
      p1c = pbeam1
      p2c = pbeam2
      p3c = p1
      p4c = p2
      k1c = qph(1,:)
      k2c = qph(2,:)
***********************************      
      imtx = 1
c     nspinconf = 2**(4+ng)
      
      if (model.eq.'matched') then
****  EXACT REAL PAIRS
         if (npairs.gt.0) then            
cphsp
ctest phspace volume
c            summt2 = 1000.d0
c            return
            q1 = pbeam1
            q2 = pbeam2
            q3 = p1
            q4 = p2
            q5 = p3
            q6 = p4
            
            if (npairs.eq.10) then
***   pi0
c               summt2 = 1.d0
c               return
               call ggpi0ff(dot(q1-q3,q1-q3),dot(q2-q4,q2-q4),ff,ffnorm) ! I use mine...

c               write(*,'(A,F20.18,A)')'  t13 = ',dot(q1-q3,q1-q3),'D0'
c               write(*,'(A,F20.18,A)')'  t24 = ',dot(q2-q4,q2-q4),'D0'
c               write(*,'(A,F20.18,A)')'  FF  = ',ff,'d0'
c               write(*,'(A,F20.18,A)')'  FFN = ',ffnorm,'d0'
               
c               call elmatpi0quad(q1,q2,q3,q4,q5+q6,ff,summt2q)
               call elmatpi0(q1,q2,q3,q4,q5+q6,ff,summt2)

c               r = summt2/summt2q
c               if (abs(r-1.d0).gt.1d-4) then
c                  print*,r
c               endif
c     call elmatczyz(q1,q2,q3,q4,q5,ff,czyz)         
***
               return
            endif
            
            if (Qmu.gt.0.d0) then
               if (npairs.eq.1) then
                  call realpairs(-q3,q2,-q1,q4,q5,q6,ame,ammu,summt2)                                    
               else
!     it works for mu pairs!
                  call realpairs(q2,-q3,q4,-q1,q5,q6,ammu,ame,summt2)
               endif
            else
               if (npairs.eq.1) then
                  call realpairs(q1,q2,q3,q4,q5,q6,ame,ammu,summt2)
               else
!     it works for mu pairs!
                  call realpairs(q2,q1,q4,q3,q5,q6,ammu,ame,summt2)
               endif
            endif
#ifdef RECOLA
            IF (npairs.eq.1) THEN
cc               if (summt2.lt.0.d0) then
!     without identical particles
c         if (Qmu.gt.0.d0) then
c            call recola_ampl(7,pbeam1,pbeam2,p1,p2,p3,p4,
c     .            rclborn,rclloop)            
c         else
c            call recola_ampl(6,pbeam1,pbeam2,p1,p2,p3,p4,
c     .           rclborn,rclloop)
c         endif
c         rclborn = rclborn / dble(rclscale*rclscale)
c         print*,rclborn/summt2*4.d0-1.d0
c         summt2 = 4.d0*rclborn
! with identical particles         
c             if (Qmu.gt.0.d0) then
c                call recola_ampl(9,pbeam1,pbeam2,p1,p2,p3,p4,
c     .               rclborn,rclloop)
c             else
c                call recola_ampl(8,pbeam1,pbeam2,p1,p2,p3,p4,
c     .               rclborn,rclloop)
c             endif
c             print*,4.d0 * 2.d0*rclborn,summt2,
c     .            (rclborn*8.d0-summt2)/summt2
c             endif
          ELSE
!     without identical particles
c         if (Qmu.gt.0.d0) then
c            call recola_ampl(7,pbeam1,pbeam2,p1,p2,p3,p4,
c     .           rclborn,rclloop)
c         else
c            call recola_ampl(6,pbeam1,pbeam2,p1,p2,p3,p4,
c     .           rclborn,rclloop)
c         endif
c         print*,rclborn/summt2*4.d0!/dble(rclscale*rclscale)
c         summt2 = 4.d0*rclborn
! with identical particles         
c             if (Qmu.gt.0.d0) then
c                call recola_ampl(11,pbeam1,pbeam2,p1,p2,p3,p4,
c     .               rclborn,rclloop)
c             else
c                call recola_ampl(10,pbeam1,pbeam2,p1,p2,p3,p4,
c     .               rclborn,rclloop)
c           endif
c             print*,rclborn/summt2*4.d0 * 2.d0,
c     .            (rclborn*8.d0-summt2)/summt2             
          ENDIF
          return
#endif
          return
       endif
****  END EXACT REAL PAIRS
       
       if (ng.eq.0) then
          call twototwo(pbeam1,pbeam2,p1,p2,tttmeq2)
          summt2 = tttmeq2
          imtx = 1
! to get the new weight, cwwp(x)*f(vpt(x)) * default_weight          
          cwvp = cwvpl/tttmeq2 ! normalized to the default weight, usually with full running
c          print*,sum(cwvp)/extraweights(1)
          return
       endif

       reweightLO  = 0.d0
       reweightNLO = 1.d0
       if (ng.eq.1) then
          q = qph(1,:)
          summt2 = radmue(pbeam1,pbeam2,p1,p2,q)          
          imtx   = 1
******************************************
          if (ord.eq.'born'.or.ord.eq.'alpha') then
             cwvp = cwvpl/summt2 ! normalized to the default (usually the full running)
c             print*,sum(cwvp)/extraweights(1)
             return
          endif
****************************
***** 1G1L
          sff = soft_integral_mue(egmin,pbeam1,pbeam2,p1,p2)
          tree = summt2
            
          imtx = 1
#ifdef RECOLA
cccccccccccccccccccccccccccccc
c            if (Qmu.lt.0.d0) then
c               call recola_ampl(3,pbeam1,pbeam2,p1,p2,q,p1b,
c     .              rclborn,rclloop)
c            else
c               call recola_ampl(4,pbeam1,pbeam2,p1,p2,q,p1b,
c     .              rclborn,rclloop)
c            endif
c            rclborn = 4.d0*rclborn
c            rclloop = 4.d0*rclloop
c            rclloop = rclloop + (1.d0+sff)*rclborn
c            summt2  = rclloop
c            return
ccccccccccccccccccccccccccccc
#endif
          call getQEDcts(cts)
          
          if (Qmu.gt.0.d0) then
*** crossing symmetry
! for Qmu = -1 it's called inside radmuesum with correct momenta, no need to recall it here
             call loadinvariants1g(-p1,pbeam2,-pbeam1,p2,q)
             p1c = -p1
             p3c = -pbeam1
          endif
***            
#ifndef COLLIER
          call clearcache
          call loadQEDfuns1g1L
#else
** for collier the cache is cleared inside the functions
          call loadQEDfuns1g1LCOLLIER
#endif
! using axial gauge for external photon, see e.g. https://arxiv.org/pdf/1403.6557.pdf
! no, defaulting to no-axial, the file from form has less operations
          udno  = 0.d0
          udyes = 0.d0
          cwvp = cwvpl * (1.d0+cts+sff)
          if (Qmu.lt.0.d0) then
             udno  = muemue1g1Lnoupdown(pbeam1,pbeam2,p1,p2,q)
             cwvp  = cwvp + cwvpl             
             udyes = muemue1g1Lupdown(pbeam1,pbeam2,p1,p2,q)
             cwvp  = cwvp + cwvpl
          else
             udno  = muemue1g1Lnoupdown(-p1,pbeam2,-pbeam1,p2,q)
             cwvp  = cwvp + cwvpl             
             udyes = muemue1g1Lupdown(-p1,pbeam2,-pbeam1,p2,q)
             cwvp  = cwvp + cwvpl             
          endif
          xxx    = 0.d0 ! -1.5d0*alpha/pi ! XXX IS NOT NEEDED ANYMORE
          summt2 = (1.d0+cts+sff+xxx)*tree + udno + udyes
************
          cwvp = cwvp/summt2
************
          extraweights(0) = 1.d0
          extraweights(1) = cwvp(1)+cwvp(2)+cwvp(5)
          extraweights(2) =
     .         cwvp(1)*1.d0/(1.d0-dal24)/(1.d0-dal24) +
     .         cwvp(2)*1.d0/(1.d0-dal24)/(1.d0-dal13) +
     .         cwvp(5)*1.d0/(1.d0-dal13)/(1.d0-dal13)
***********          
c          print*,sum(cwvp)/extraweights(1) ! OK! check it again, for QRe,QRmu=1. 1 0 and 0 1 are ok
          reweightNLO = tree/summt2
          imtx = 1
          return
       endif
**********************************************************************************************
       if (ng.eq.2) then
          reweightNLO = 0.d0
          icount = icount + 1
          q1 = qph(1,:)
          q2 = qph(2,:)
          if (Qmu.gt.0.d0) then            
             summt2 = muemuegg(pbeam1,pbeam2,p1,p2,q1,q2)
          else
             summt2 = muemuegg(-p1,pbeam2,-pbeam1,p2,q1,q2)
          endif
          cwvp = cwvpl / summt2
!         print*,sum(cwvp)/extraweights(1) ! OK
**
c            if (Qmu.lt.0.d0) then
c               call recola_ampl(1,pbeam1,pbeam2,p1,p2,q1,q2,rcl,rcl1l) ! mu-
c            else
c               call recola_ampl(2,pbeam1,pbeam2,p1,p2,q1,q2,rcl,rcl1l) ! mu+
c            endif
c            rcl = 8.d0*rcl
c            summt2 = rcl
c            print*,summt2,rcl,(summt2-rcl)/summt2
c            stop
c            print*,summt2,rcl,(summt2-rcl)/rcl
c            if (icount.eq.40) stop
**
          imtx = 1
          return
       endif         
      endif
*******************
**      PS       **
*******************
      if (model.eq.'ps') then
         print*,'wrong model. Stopping.'
         stop
      endif
******
      if (model.eq.'eikonal') then
         print*,'wrong model. Stopping.'
         stop
      endif
******
      return
      end
****************************************************************************
      subroutine get_us(p,m,u)
      implicit double precision (a-h,m,o-z)
      double precision p(0:3)
      double complex u(2,4)!,utmp(4)
      double complex im
      common/getus/im,ifirst
      data ifirst/0/
      if (ifirst.eq.0) then
         im = (0.d0,1.d0)
         ifirst = 1
      endif
      den     =  sqrt(p(0)+m)
! u1 for e- ! spin +1 along z in the rest frame
      u(1,1) =  den
      u(1,2) =  0.d0
      u(1,3) =  p(3)/den
      u(1,4) =  (p(1)+im*p(2))/den

!     u2 for e- ! spin -1 along z in the rest frame
      u(2,1) =  0.d0
      u(2,2) =  den
      u(2,3) =  conjg(u(1,4)) ! (p1(1)-im*p1(2))/(p1(0)+me)
      u(2,4) = -u(1,3)        ! (-p1(3)/(p1(0)+me)
      return
      end
****************************************
      subroutine get_vs(p,m,v)
      implicit double precision (a-h,m,o-z)
      double precision p(0:3)
      double complex v(2,4)
      double complex im
      common/getvs/im,ifirst
      data ifirst/0/
      if (ifirst.eq.0) then
         im = (0.d0,1.d0)
         ifirst = 1
      endif
      den     =  sqrt(p(0)+m)
! v1 for e-  ! spin +1 along z in the rest frame
      v(1,1) =  (p(1)-im*p(2))/den
      v(1,2) = -p(3)/den
      v(1,3) =  0.d0
      v(1,4) =  den
! v2 for e-  ! spin -1 along z in the rest frame
      v(2,1) = -v(1,2)        ! p2(3)/(p2(0)+me)
      v(2,2) =  conjg(v(1,1)) !(p2(1)+im*p2(2))/(p2(0)+me)
      v(2,3) =  den
      v(2,4) =  0.d0
      return
      end
***************************************
      subroutine get_spinors(p1,p2,p3,p4,s1,s2,s3,s4)
c p1 = e
c p2 = e
c p3 = mu
c p4 = mu
      implicit double precision (a-h,o-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double complex s1(2,4),s2(2,4),s3(2,4),s4(2,4)
      double complex im
      common/parameters/ame,ammu,convfac,alpha,pi

      double precision me,mm
      common/getspinors/im,me,mm,ifirst
      data ifirst/0/
      
      if (ifirst.eq.0) then
         im = (0.d0,1.d0)
         me = ame
         mm = ammu
         ifirst = 1
      endif

      den     =  sqrt(p1(0)+me)
! u1 for e-
      s1(1,1) =  den
      s1(1,2) =  0.d0
      s1(1,3) =  p1(3)/den
      s1(1,4) =  (p1(1)+im*p1(2))/den
! u2 for e-
      s1(2,1) =  0.d0
      s1(2,2) =  den
      s1(2,3) =  conjg(s1(1,4)) ! (p1(1)-im*p1(2))/(p1(0)+me)
      s1(2,4) = -s1(1,3)        ! (-p1(3)/(p1(0)+me)

      den     =  sqrt(p3(0)+mm)
! u1 for mu-
      s3(1,1) =  den
      s3(1,2) =  0.d0
      s3(1,3) =  p3(3)/den
      s3(1,4) =  (p3(1)+im*p3(2))/den
! u2 for mu-
      s3(2,1) =  0.d0
      s3(2,2) =  den
      s3(2,3) =  conjg(s3(1,4)) !(p3(1)-im*p3(2))/(p3(0)+mm)
      s3(2,4) = -s3(1,3)        !-p3(3)/(p3(0)+mm)

      den     =  sqrt(p2(0)+me)
! v1 for e-
      s2(1,1) =  (p2(1)-im*p2(2))/den
      s2(1,2) = -p2(3)/den
      s2(1,3) =  0.d0
      s2(1,4) =  den
! v2 for e-
      s2(2,1) = -s2(1,2)        ! p2(3)/(p2(0)+me)
      s2(2,2) =  conjg(s2(1,1)) !(p2(1)+im*p2(2))/(p2(0)+me)
      s2(2,3) =  den
      s2(2,4) =  0.d0

      den     =  sqrt(p4(0)+mm)
! v1 for mu-
      s4(1,1) =  (p4(1)-im*p4(2))/den
      s4(1,2) = -p4(3)/den
      s4(1,3) =  0.d0
      s4(1,4) =  den
! v2 for mu-
      s4(2,1) = -s4(1,2)        ! p4(3)/(p4(0)+mm)
      s4(2,2) =  conjg(s4(1,1)) !(p4(1)+im*p4(2))/(p4(0)+mm)
      s4(2,3) =  den
      s4(2,4) =  0.d0
      return
      end
*************************************************************
      subroutine twototwo(p1,p2,p3,p4,smtq2)
#ifdef COLLIER
      use collier
#endif      
      implicit double precision (a-h,l,m,o-z)
      dimension ptmp(0:3),p1b(0:3),p2b(0:3),p3b(0:3),p4b(0:3)
      double complex chis,chit,xi(1)
      double complex vpolc,ctmp,cmuceZ
      double precision leptonicbubble
      double precision Pi02L_Boncianietal
      external Pi02L_Boncianietal
      external vpolc,leptonicbubble
      integer hh
      character*2 fs
      common/finalstate/fs
      common/parameters/ame,ammu,convfac,alpha,pi
      common/tautopmasses/taumass,topmass
      
      double complex sp1(2,4),sp2(2,4),sp3(2,4),sp4(2,4)
      double complex s1l(4),s2l(4),s3l(4),s4l(4),sa(4),sb(4),sa5(4)
      double complex amplborn,ampv(0:15),amplbsm,atmp,amplowz
      double complex c1,c2,boxamp,vertamp,tmp2l,ctmp1,ctmp2,ctmp3
      double complex ctest!, c5mu(0:3),c5e(0:3)

      double complex cmut(0:3),cet(0:3),c5mut(0:3),c5et(0:3),cImut,cIet

      double precision ch(4),amassess(4)
      double complex yfsij(4,4),cyfs,cyfse,cyfsmu,cyfsbox
      
      common/photoncutoff/egmin,egmin2

      character*6 ord
      common/qedORDER/ord
      
      integer ifirst
      common/ifirstborn/ifirst,icount
      data ifirst,icount /0,0/
      
      common/bsmparams/bsmmass,yae,yam,bsmlambda,gagg,ibsm      

      double precision mz,mz2,gz,gz2,mw,mw2,stw,ctw
      common/zparameters/mz,mz2,gz,gz2,mw,mw2,stw,ctw,gv,ga

      double precision phmass
      common/photonmasslambda/phmass

      common/sffudandnoud/sffud,sffnoud

      double complex virtvupvdwn      
      common/testvirtvevmu/virtvupvdwn      
*****************
      double precision sffij(4,4)
      common/softFFij/sffij
*****************
      logical is2L
      double precision p3lab(0:3),p4lab(0:3)
****
      double precision dalo(3),danlo(3),dannlo(3)
      common/dalonlonnlo/dalo,danlo,dannlo

      common/FsPi02Ldebug/F1,F2,F3,F4,F5,F6,F7,F8
      
      common/forbubbles/ctandvmph0,boxmph0
      integer bubblemask(5)
      double precision colch2(5)
      common/colorandch2factors/colch2,bubblemask      
*********
*** for 2 loops, calculated in getQEDvertexes
      double complex amplve,amplvmu
      common/vertexeslegs/amplve,amplvmu
! calculated in getQEDcts
      double precision ctsmu,ctse
      common/qedcts_e_mu/ctsmu,ctse
*****
#include "invariants.h"

      integer ionoff(4,4)
      
      if (ifirst.eq.0) then

         extraweight = 1.d0

         mz2 = mz*mz
         mw2 = mw*mw
         gz2 = gz*gz

         stw = sqrt(1.d0-mw2/mz2)
         ctw = sqrt(1.d0-stw*stw)

         gv = -0.5d0+2.d0*stw*stw
         gv = gv*0.5d0/stw/ctw
         ga = -0.25d0/stw/ctw

         gv = ec*gv
         ga = ec*ga

         ifirst = 1
      endif

      smtq2 = 0.d0
      
      ptmp = p1+p2
      s    = dot(ptmp,ptmp)
      ptmp = p1-p3
      t    = dot(ptmp,ptmp)

      tm1   = 1.d0/t
      vpt   = 1.d0
      sff   = 0.d0
      cwvpl = 0.d0
***   analytic check with form
c      ptmp = p1-p4
c      u    = dot(ptmp,ptmp)
c      s2   = s*s
c      t2   = t*t
c      u2   = u*u
c         xform = !! OK !! already averaged over IS
c     & 2.D0*t2 + 4.D0*s*t + 4.D0*s2 + 4.D0*mm4 - 8.D0*mm2*s + 4.D0*
c     & me4 - 8.D0*me2*s + 8.D0*me2*mm2
c         xform = 4.d0*ec2*ec2/t2*xform                  
c         smtq2 = xform
c         return
******************************

      if (ioneloop.eq.1) then
         call loadinvariants(p1,p2,p3,p4)
         call getQEDcts(cts)
#ifdef COLLIER
         call InitEvent_cll(1)
#else
         call clearcache        ! this is after getQEDcts
! because there only at first event Bs are called...
#endif
         call getVP1loop(t,dal,dah) ! it's in common at the beginning!
         if (ifixedorder.eq.0) then
            vpt  = 1.d0/(1.d0-dal-dah)
            dalw = dal          ! for weights
            dahw = dah          ! for weights
            dal  = 0.d0
            dah  = 0.d0
         endif
         call loadQEDvertexfuns
         call loadQEDboxesfuns
         sff = soft_integral_mue(egmin,p1,p2,p3,p4)
      endif

      if (ioneloop.eq.0.and.ifixedorder.eq.0) then
         call getVP1loop(t,dal,dah)
         dalw = dal             ! for weights
         dahw = dah             ! for weights
         vpt  = 1.d0/(1.d0-dal-dah)
      endif
         
      ampm2     = 0.d0
      born      = 0.d0
      oneloop   = 0.d0
      oneloopm2 = 0.d0
      reweightLO  = 1.d0
      reweightNLO = 1.d0
      boxmph0  = 0.d0
      ctandvmph0  = 0.d0
      
      is2L = .false.
      if (ord.eq.'alpha2') is2L = .true.
     
      yfs2Lee   = 0.d0
      yfs2Lmumu = 0.d0
      yfs2Lemu  = 0.d0
      yfs2Lbb   = 0.d0
      yfs2Lbe   = 0.d0
      yfs2Lbmu  = 0.d0
      yfs2L     = 0.d0
      if (is2L) then
         call YFSvirt(yfsij)         
c         call YFSvirt_onlyIR(yfsij)
         cyfse   = yfsij(2,2)+yfsij(2,4)+yfsij(4,4)
         cyfsmu  = yfsij(1,1)+yfsij(1,3)+yfsij(3,3)
         cyfsbox = yfsij(1,2)+yfsij(3,4)+yfsij(1,4)+yfsij(2,3)
         cyfse   = 0.5d0*cyfse  ! 1/2 perche' nella routine era per 2Re
         cyfsmu  = 0.5d0*cyfsmu
         cyfsbox = 0.5d0*cyfsbox
      else
         cyfse   = (0.d0,0.d0)
         cyfsmu  = (0.d0,0.d0)
         cyfsbox = (0.d0,0.d0)
      endif

ccc TEST YFS oal
c         call YFSvirt(yfsij)
c         cyfse   = yfsij(2,2)+yfsij(2,4)+yfsij(4,4)
c         cyfsmu  = yfsij(1,1)+yfsij(1,3)+yfsij(3,3)
c         cyfsbox = yfsij(1,2)+yfsij(3,4)+yfsij(1,4)+yfsij(2,3)
c         cyfse   = 0.5d0*cyfse  ! 1/2 perche' nella routine era per 2Re
c         cyfsmu  = 0.5d0*cyfsmu
c         cyfsbox = 0.5d0*cyfsbox
ccc END TEST YFS oal
      
      cyfs    = cyfse+cyfsmu+cyfsbox
      do hh = 0,7              ! only half, below the factor of 2!!
c         do hh = 0, 15           ! only half, below the factor of 2!!            
c         do hh = 8, 15 !15 ! only half, below the factor of 2!!
c         do hh = 15, 15 !15 ! only half, below the factor of 2!!
         h1  = mod(hh,16)/8 + 1
         h2  = mod(hh, 8)/4 + 1
         h3  = mod(hh, 4)/2 + 1
         h4  = mod(hh, 2)   + 1
c            s1l = sp1(h1,:)
c            s2l = sp2(h2,:)
c            s3l = sp3(h3,:)
c            s4l = sp4(h4,:)
c            call getcurrent(s4l,s2l,ce,cIe,c5e)    ! c5x put in common!
c            call getcurrent(s3l,s1l,cmu,cImu,c5mu)

         if (Qmu.lt.0.d0) then
            call getcurrentnew(Qmu,h3,h1,p3,p1,mm,mm,cmu,cImu,c5mu)
         else
            call getcurrentnew(Qmu,h3,h1,p1,p3,mm,mm,cmu,cImu,c5mu)
         endif            
         call getcurrentnew(-1.d0,h4,h2,p4,p2,me,me,ce,cIe,c5e)
         
         cmuce   = dotcmom(cmu,ce)
            
!! for the Z I need to use dotcmom(v*cmu+a*c5mu,v*ce+a*c5e)
c         cmuceZ = dotcmom(gv*cmu+ga*c5mu,gv*ce+ga*c5e)
         amplborn = cmuce*im*ec2*tm1 
!! for the Z
c           amplborn = amplborn + im*cmuceZ/(t-mz2+im*mz*gz)
         bbar    = amplborn*conjg(amplborn)
         born    = born + bbar
         if (ioneloop.eq.1) then
            call loadspinvariants(p1,p2,p3,p4,cmu,ce)
            call getQEDvertexes(vertamp)
            call getQEDboxes(boxamp)
            
            oneloop = bbar * (1.d0*cts + 1.d0*sff + 2.d0*(dal+dah))
     .           + 2.d0*dreal(amplborn*conjg(vertamp))
     .           + 2.d0*dreal(amplborn*conjg(boxamp))
            
            boxmph0 = boxmph0 + 2.d0*dreal(amplborn*conjg(boxamp))
            ctandvmph0 = ctandvmph0 + bbar * cts
     .           + 2.d0*dreal(amplborn*conjg(vertamp))
ccc TEST YFS oal
c            oneloop = bbar * (1.d0*cts + 1.d0*sff + 2.d0*(dal+dah))
c     .           + 2.d0*dreal(amplborn*conjg(vertamp))
c     .           + 2.d0*dreal(amplborn*conjg(amplborn*cyfsbox))
ccc END TEST YFS oal            
*** for 2L **********************************************************
            if (is2L) then
               ctmp1 = amplve  + 0.5d0*ctse*amplborn
               ctmp2 = amplvmu + 0.5d0*ctsmu*amplborn
               
               tmp2l     = vertamp + 0.5d0*cts*amplborn + boxamp
               oneloopm2 = oneloopm2 + tmp2l*conjg(tmp2l)
! approximates the 2 loop vertex e                  
               ctmp3 = cyfse*cyfse*0.5d0 * amplborn
     .              + (ctmp1-cyfse*amplborn)*cyfse ! funziona!!                             
               yfs2Lee = yfs2Lee + 2.d0*dreal(amplborn*conjg(ctmp3))                  
! approximates the 2 loop vertex mu
               ctmp3 = cyfsmu*cyfsmu*0.5d0 * amplborn
     .              + (ctmp2-cyfsmu*amplborn)*cyfsmu ! funziona!!                             
               yfs2Lmumu = yfs2Lmumu
     .              + 2.d0*dreal(amplborn*conjg(ctmp3))
! approximates the 2 loop with 1 vertex on e and one on mu                  
               ctmp3 = cyfse*cyfsmu * amplborn
     .              + (ctmp1-cyfse *amplborn)*cyfsmu
     .              + (ctmp2-cyfsmu*amplborn)*cyfse
               yfs2Lemu = yfs2Lemu
     .              + 2.d0*dreal(amplborn*conjg(ctmp3))
! approximates the full two loop, including boxes          
               ctmp3 = cyfs*cyfs*0.5d0*amplborn
     .              + (ctmp1+ctmp2+boxamp - cyfs*amplborn)*cyfs ! full 2L approximation
               yfs2L = yfs2L + 2.d0*dreal(amplborn*conjg(ctmp3))
            endif
** end for 2L ******************************************************
         endif
         ampm2 = ampm2 + bbar + oneloop
      enddo

****  Z LO always added!
      bornwz = 0.d0
      DZ     = 0.d0
cc      goto 42
      do hh = 0, 15 ! notice, I must sum over all helicities because of CP violation
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
         cmuce   = dotcmom(cmu,ce)
         amplowz = cmuce*im*ec2*tm1 ! *vpt ! no vpt here, added below
!! for Z I need to use dotcmom(v*cmu+a*c5mu,v*ce+a*c5e)
         cmuceZ  = dotcmom(gv*cmu+ga*c5mu,gv*ce+ga*c5e)
         amplowz = amplowz + im*cmuceZ/(t-mz2+im*mz*gz)
         bornwz  = bornwz + amplowz*conjg(amplowz)
      enddo
      DZ = bornwz - 2.d0*born ! notice: DZ should have a vpt factorized in front of 2*M^*_z*M_qed
**** end Z LO
 42   continue

      
      vpt2 = vpt*vpt
      
      ampm2      = 2.d0 * ampm2 * vpt2 ! because summing only on half of helicities
      born       = 2.d0 * born  * vpt2
!      DZ         = DZ * vpt ! I assume there is no vpt dependency, but there is in front of 2*M^*_z*M_qed
      
      ampm2      = ampm2 + DZ ! adding Z at LO

      boxmph0    = 2.d0 * boxmph0 ! no vp here
      ctandvmph0 = 2.d0 * ctandvmph0 ! no vp here
      ampm2NLO   = ampm2
      blobsexact = 0.d0

      extrareweightLO = 1.d0
! [[ --- THE FOLLOWING WAS USED FOR NNLO PAIR PAPER
c      IF (ionlybubbles.gt.0) THEN
c         vp    = dot_product(vpvect,bubblemask)         
c         extrareweightLO = 1.d0/vpt/vpt        
c         xx1      = alf_lep_carlo(1,t) ! to fill the common
c         xx1      = alf_lep_carlo(2,t) ! to fill the common
c         xx1      = alf_lep_carlo(3,t) ! to fill the common
c         danlo(1) = danlo(1)*bubblemask(1)
c         danlo(2) = danlo(2)*bubblemask(2)
c         danlo(3) = danlo(3)*bubblemask(3)
c         vpnlo    = danlo(1) + danlo(2) + danlo(3)
c         
c         bornclean = born/vpt/vpt
c         boxbubble  = 0.d0
c         vertbubble = 0.d0
c         ibubblev  = ibubbleonv
c         ibubbleb  = ibubbleonboxsub
c         if (QRe.eq.0.or.QRmu.eq.0) ibubbleb = 0
c         call leptonicbubbles(p1,p2,p3,p4,ibubblev,ibubbleb,
c     .        vertbubble,boxbubble)
c!! that's good!            
c         ampm2 = (ctandvmph0 + boxmph0 + bornclean*sff)*vp*2.d0
c     .        + boxbubble + vertbubble
c     .        + 1.d0 * bornclean * (3.d0*vp*vp + 2.d0*vpnlo)
c      ENDIF
! --- ]]

      if (iexnnlobubbles.gt.1) then
! this is the "non-factorizable" NNLO virtual pair contribution
         boxbubble  = 0.d0
         vertbubble = 0.d0
         ibubblev  = ibubbleonv
         ibubbleb  = ibubbleonboxsub
         if (QRe.eq.0.or.QRmu.eq.0) ibubbleb = 0
         call leptonicbubbles(p1,p2,p3,p4,ibubblev,ibubbleb,
     .        vertbubble,boxbubble)
!! that's good!            
         blobsexact = (boxbubble + vertbubble)*vpt2 ! with vpt^2 I add VP on non-loop photon...
      endif
      
      ampm2 = ampm2 + blobsexact ! added to the 1L contribution...
*******************************************************
      if (is2L) then
*******************************************************        
         oneloopm2 = 2.d0 * oneloopm2 * vpt2 ! |1L|^2 

         yfs2Lee   = 2.d0*yfs2Lee   * vpt2
         yfs2Lmumu = 2.d0*yfs2Lmumu * vpt2
         yfs2Lemu  = 2.d0*yfs2Lemu  * vpt2
         yfs2L     = 2.d0*yfs2L     * vpt2
         
         true2L = twoloop_twobody(s,t,u) * vpt2
         sv1L   = (ampm2-born*(1.d0+sff))*sff
         
         DeltaYFS = yfs2L - yfs2Lee - yfs2Lmumu - yfs2Lemu ! eureka!
         
         dueloop = true2L+DeltaYFS+oneloopm2+sv1L+0.5d0*sff*sff*born
         
         ampm2     = ampm2 + dueloop ! ottimo!
         
         ampm2NNLO = ampm2
      endif     
         
      extraweights(0) = 1.d0
      extraweights(1) = ((ampm2 - DZ )*1.d0/vpt2 + DZ) / ampm2
c      extraweights(2) = 1.d0/vpt2/(1.d0-dalw)/(1.d0-dalw)
      extraweights(2) = ((ampm2 - DZ )*1.d0/vpt2/(1.d0-dalw)/(1.d0-dalw)
     .     + DZ) / ampm2

      reweightLO  = (born + DZ)/ampm2  * extrareweightLO
      
      reweightNLO = ampm2NLO/ampm2

      if (ord.eq.'born') reweightNLO = 0.d0
         
      smtq2 = ampm2
      
****** VP reweighting stuff
      cwvpl(0) = DZ
      cwvpl(1) = (ampm2 - DZ) / vpt2
******
      return
      end
**************************************************************************
      subroutine eikonal_factor(eikng,ng,p1,p2,p3,p4,q,qph,eikonal)
      implicit double precision (a-h,o-z)
      character*6 eikng
      parameter (imaxph = 40)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension qph(imaxph,0:3),q(0:3),peik(0:3)
      character*2 fs
      common/finalstate/fs
      common/singlefactors/singleeikfact(imaxph)
      common/parameters/ame,ammu,convfac,alpha,pi

      eikonal = 1.d0
      test = 1.d0
      if (ng.eq.0) return

      nglocal = ng
      if (eikng.eq.'single') nglocal = 1

      ec2 = 4.d0*pi*alpha ! electron charge^2
      do k = 1,nglocal
         if (eikng.eq.'multip') then
            q(0) = qph(k,0)
            q(1) = qph(k,1)
            q(2) = qph(k,2)
            q(3) = qph(k,3)
         endif
         oop1dk = 1.d0/dot(p1,q)
         oop2dk = 1.d0/dot(p2,q)
         if (fs.eq.'ee'.or.fs.eq.'mm') then
            oop3dk = 1.d0/dot(p3,q)
            oop4dk = 1.d0/dot(p4,q)
         elseif (fs.eq.'gg') then
            oop3dk = 0.d0
            oop4dk = 0.d0
         endif
c         do i = 0,3
c            peik(i)=p1(i)*oop1dk-p2(i)*oop2dk-p3(i)*oop3dk+p4(i)*oop4dk
c         enddo
c         singlephoton = -ec2*dot(peik,peik)
         singlephoton = -(
     .          dot(p1,p1)*oop1dk**2
     .        + dot(p2,p2)*oop2dk**2
     .        + dot(p3,p3)*oop3dk**2
     .        + dot(p4,p4)*oop4dk**2
     .        - 2.d0*dot(p1,p2)*oop1dk*oop2dk
     .        - 2.d0*dot(p3,p4)*oop3dk*oop4dk
     .        - 2.d0*dot(p1,p3)*oop1dk*oop3dk
     .        + 2.d0*dot(p1,p4)*oop1dk*oop4dk
     .        + 2.d0*dot(p2,p3)*oop2dk*oop3dk
     .        - 2.d0*dot(p2,p4)*oop2dk*oop4dk)*ec2               
         
         eikonal = eikonal * singlephoton
         if (eikng.eq.'multip') singleeikfact(k) = singlephoton
      enddo

      eikonal = eikonal/factorial(nglocal)      
      return
      end
**************************************************************************
      subroutine approxNg(ng,p1,p2,p3,p4,qph,me)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),qph(40,0:3),qt(40,0:3)
      dimension tmp(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi
      e2 = alpha*4.d0*pi
      nph = ng + 2
      me = 0.d0
      do i = 0,3
         tmp(i)  = p1(i) + p2(i)
         qt(1,i) = p3(i)
         qt(2,i) = p4(i)
         if (ng.gt.0) then
            do k = 1,ng
               qt(k+2,i) = qph(k,i)
            enddo
         endif
      enddo
      s = dot(tmp,tmp)
      sum  = 0.d0
      prod = 1.d0
      do k = 1,nph
         do i = 0,3
            tmp(i) = qt(k,i)
         enddo
         xi = dot(tmp,p1)
         yi = dot(tmp,p2)
         sum = sum + xi*yi * (xi**2 + yi**2)
         prod = prod * xi*yi
      enddo
      me = sum/prod * s**(nph-2) * e2**nph
      return
      end
*************************************************************************
      subroutine approxNg_SINGLE(ng,p1,p2,p3,p4,qph,me)
! only one term
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),qph(40,0:3),qt(40,0:3)
      dimension tmp(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi
      e2 = alpha*4.d0*pi
      nph = ng + 2
      me = 0.d0
      do i = 0,3
         tmp(i)  = p1(i) + p2(i)
         qt(1,i) = p3(i)
         qt(2,i) = p4(i)
         if (ng.gt.0) then
            do k = 1,ng
               qt(k+2,i) = qph(k,i)
            enddo
         endif
      enddo
      s = dot(tmp,tmp)
      sum  = 0.d0
      prod = 1.d0
      do k = 2,2
         do i = 0,3
            tmp(i) = qt(k,i)
         enddo
         xi = dot(tmp,p1)
         yi = dot(tmp,p2)
         sum = sum + xi*yi * (xi**2 + yi**2)
      enddo
      do k = 1,nph
         do i = 0,3
            tmp(i) = qt(k,i)
         enddo
         xi = dot(tmp,p1)
         yi = dot(tmp,p2)
         prod = prod * xi*yi
      enddo
      me = sum/prod * s**(nph-2) * e2**nph  * 8.d0
      return
      end
