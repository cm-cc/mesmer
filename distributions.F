      subroutine setdistributionfilesname(nw)
! nw is the number of weights
#include "distributions_inc.F"
      character*300 programma
      character*300 distrstring
      character*100   ec,ac
      character*7   bw
      character*6   ord
      common/qedORDER/ord
      character*200  path
      common/cmnpath/path
      
      do k = 1,7
         bw(k:k) = ' '
      enddo
      do k = 1,ndistr
         do i = 1,100
            tmp(k)(i:i) = ' '
            tmp2(k)(i:i) = ' '
         enddo
      enddo

      ebeam = ecmsnom/2.d0

      call itoa(nbin,bw)
      bw = trim(bw)//'.txt'
      
******************
      tmp(1)  = 'ethlab_'
      tmp(2)  = 'eenlab_'
      tmp(3)  = 'mthlab_'
      tmp(4)  = 'menlab_'
      tmp(5)  = 't24_'
      tmp(6)  = 't25_'
      tmp(7)  = 't13_'
      tmp(8)  = 'aco_'
      tmp(9)  = 't26_'
      tmp(10) = 'ethcm_'
      tmp(11) = 'eencm_'
      tmp(12) = 'mthcm_'
      tmp(13) = 'mencm_'
      tmp(14) = 'm34_'

      tmp(15) = 'gthlab_'
      tmp(16) = 'genlab_'
      tmp(17) = 'gthcm_'
      tmp(18) = 'gencm_'
      tmp(19) = 'gecthcm_'
      tmp(20) = 'm45_'
      tmp(21) = 'acollcm_'
      tmp(22) = 'm456_'
      tmp(23) = 'eenlab_p_'
      tmp(24) = 'ethlab_p_'
      tmp(25) = 't2in_'
      tmp(26) = 'g2thcm_'
      tmp(27) = 'g2encm_'

      tmp(28) = 'th5lab_'
      tmp(29) = 'th6lab_'
      tmp(30) = 'en5lab_'
      tmp(31) = 'en6lab_'
      tmp(32) = 'm56_'

      tmp(33) = 'mubeam_th_'
      tmp(34) = 'mubeam_E_'
c      tmp(35) = 'xi3_'
c      tmp(36) = 'xi4_'
      
c      tmp(32) = 'x24_'
      
c      tmp(33) = 'emuin'
c      tmp(34) = 'thmuin'

c      tmp(35) = 'tripleprod_'
c      tmp(36) = 'acoplGA_'
c      tmp(37) = 'acoplCM_'
      
******** 2d distributions
      if (i2d.gt.0) then
         tmp2(1) = '2d_thm_thin_'
         tmp2(1) = 'xi1xi4_'
         tmp2(2) = 'xi3xi4_'
c         tmp2(2) = '2d_the_thm_nohad_'
c         tmp2(3) = '2d_t13z_s4_'
c         tmp2(4) = '2d_thm_aco_'
      endif
******************
      do k=1,ndistr
         lung(k) = len_trim(tmp(k))
      enddo
      do k=1,ndistr2
         if (i2d.gt.0) lung2(k) = len_trim(tmp2(k))
      enddo
      
      if (ord.eq.'alpha') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'oal_'
            tmp(k)(lung(k)+5:) = bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'oal_'
            if (i2d.gt.0) tmp2(k)(lung2(k)+4:) = '.txt'
         enddo
      endif

      if (ord.eq.'alpha2') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'oal2_'
            tmp(k)(lung(k)+6:) = bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'oal2_'
            if (i2d.gt.0) tmp2(k)(lung2(k)+5:) = '.txt'
         enddo
      endif
      
      if (ord.eq.'exp') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'exp_'
            tmp(k)(lung(k)+5:) = bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'exp_'
            if (i2d.gt.0) tmp2(k)(lung2(k)+4:) = '.txt'
         enddo
      endif
      
      if (ord.eq.'born') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'born_'
            tmp(k)(lung(k)+6:) = bw
         enddo
         do k=1,ndistr2
            if (i2d.gt.0) tmp2(k)(lung2(k)+1:) = 'born_'
            if (i2d.gt.0) tmp2(k)(lung2(k)+5:) = '.txt'
         enddo
      endif

      programma = path
      do k = 1,len(programma)
         if (programma(k:k).ne.' ') lpath=k
      enddo
      
      do k = 1,ndistr
         distname(k,1)(1:lpath)  = programma
         distname(k,1)(lpath+1:) = trim(tmp(k))
         if (nw.ge.2) then
            distname(k,2)(1:lpath)  = programma
            distname(k,2)(lpath+1:) = trim(tmp(k))//"_LO"
         endif
         if (nw.eq.3) then
            distname(k,3)(1:lpath)  = programma
            distname(k,3)(lpath+1:) = trim(tmp(k))//"_NLO"
         endif
      enddo
      do k = 1,ndistr2
         if (i2d.gt.0) distname2(k)(1:lpath)  = programma
         if (i2d.gt.0) distname2(k)(lpath+1:) = tmp2(k)
      enddo

***** printing
      print*,'******************************'
      write(ec,*)ndistr
      write(ac,*)nbin
      print*,'Writing also ',trim(adjustl(ec)),' distributions with ',
     .trim(adjustl(ac)),' bins within appropriate ranges.'
      print*,'The distributions names are:'
      i1 = 1
      do k = 1,ndistr
         j = len(tmp(k)(1:index(tmp(k),'_')-1))
         i2 = i1 + j - 1
         distrstring(i1:i2+2) = tmp(k)(1:j)//', '
c         print*,trim(distrstring(1:i2+2))
         i1 = i2 + 3
      enddo
      print*,trim(distrstring(1:i2))//'.'
      if (i2d.gt.0) then
         print*,'Also some 2D distributions are written.'
      endif
      print*,'******************************'
***** 
      return
      end
*******************************************************************
      subroutine distributions(sd,p3,p4,p5,p6,qph,nw)
! nw is the number of weights, sd has dimension nw
#include "distributions_inc.F"
      dimension pin1(0:3),pin2(0:3),csi(1),p5(0:3),p6(0:3)
      double precision kone(0:3)

      double precision qphlab(40,0:3)
      double precision p3lab(0:3),p4lab(0:3),p5lab(0:3),p6lab(0:3)
      common/fslabmomenta/p3lab,p4lab,p5lab,p6lab,qphlab

      double precision sd(nw)
      dimension p4test1(0:3),p3test(0:3),p4test2(0:3)
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b(0:3),p2b(0:3),iref
      character*2 fs
      common/finalstate/fs
      common/ionlyfirsttime/ifirst,icount
      data ifirst,icount /1,0/ 

      common/tmintmaxfordistributions/tmin,tmax
      common/mueexpsetup/emulab,eemin,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela

      common/beamspreadsigmas/esig1,esig2
      
      common/par/pi,radtodeg,ebeam
      
      common/photoncutoff/egminmue,egmin2

      common/parameters/ame,ammu,convfac,alpha,picommon

      integer*8 nc

      double precision lambda
      external lambda
      
      common/rrrr/r,Ei,pmui
      
      double precision qa(0:3),qb(0:3),qa2(0:3),qb2(0:3)
      double precision pluto(0:3)

      double precision p3lo(0:3),p4lo(0:3)
      double precision ptmp(0:3), ptmp2(0:3)

      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      double precision truep1lab(0:3)
      common/labmomenta/p1lab,p2lab,truep1lab,p12lab

      double precision m1,m2,m12,m22
      common/processmasses/m1,m2,m12,m22

      common/na7variable/dthna7 ! filled in cuts!
      
      common/minmaxeng/genmin,genmax
      common/udistr/u

      common/el_limits_here/theminL,themaxL,eemax
      

      common/minmaxfun/funmin,funmax
      data funmin,funmax/1d200,-1d200/

      common/xcrf/xcrf,ycrf,amodk,ckcrf
      
      common/sumdebug/summm,ptmax
      data summm,ptmax/0.d0,-1.d0/
*
***   for reweighenting
      integer nextraweights
      parameter (nextraweights=4)
      double precision weightdefault,extraweights(0:nextraweights)
      double precision reweightLO,reweightNLO,wnorm
      common/weights/reweightLO,reweightNLO,weightdefault,extraweights
***
**************************      
      double precision invmass,invmasswb
      common/invariantmasses/invmass,invmasswb ! calculated in cuts

*** for pairs, momenta pinside_* already in the lab (calculated in cuts routine)
      integer inside_e(3),inside_mu(3)
      double precision pinside_e(0:3),pinside_mu(0:3),pout_1(0:3)
      double precision pout_2(0:3)
      common/insidecuts/inside_e,inside_mu,pinside_e,pinside_mu,
     .     pout_1,pout_2
****      

      common/elasticitiesdistr/dist,distf
      
*** from phasepace for pairs      
      double precision csid(8)
      common/csifordistr/csid
***
      common/nphot_mode/nphotmode

****************      
      if (nw.eq.0) return
****************      
      
      if (ifirst.eq.1) then

         call setdistributionfilesname(nw)
         
         ebeam = ecmsnom/2.d0

         Ei = emulab + m2

         s  = ecmsnom*ecmsnom
         
         pmui = sqrt(emulab*emulab-m1*m1)
         
         pcm = 0.5d0*
     .        sqrt(lambda(ecmsnom*ecmsnom,m1*m1,m2*m2)/ecmsnom/ecmsnom)
         emumaxcm = sqrt(pcm*pcm+m1*m1)
         eemaxcm  = sqrt(pcm*pcm+m2*m2)
                  
         genmin = 1d10
         genmax = 0.d0
                  
         pi = 4.d0 * atan(1.d0)
         radtodeg = 1d3 ! 180.d0/pi

         be = sqrt(tridot(pin2,pin2))/pin2(0)
         ga = 1.d0/sqrt(1.d0-be*be)
         a  = ga*(1.d0-be)*egminmue
         b  = ga*(1.d0+be)*egminmue   

         atleasteg = b
         
         bmi(4) = 0.d0
         bma(4) = emulab
         
         bmi(5) =  -0.143d0 !tmin
         bma(5) =  0.d0 !tmax
         bmi(5) =  tmin
         
         bma(5) =  tmax

ccc         bma(5) = 2.d0*ame*ame - 2.d0*1.d0*ame
         
         bmi(25) =  bmi(5)
         bma(25) =  bma(5)
         
calculating xmin, xmax, theminL, themaxL in the lab (L for local)

         xmin = xvar(tmax)
         xmax = xvar(tmin)
         
         r2 = emulab**2-m12
         r2 = r2 / (emulab+m2)/ (emulab+m2)

         r = sqrt(r2)
         
         eemax = m2 - 0.5d0*tmin/m2         

         theminL = acos(min(1.d0,1.d0/r*sqrt((eemax-m2)/(eemax+m2))))
         themaxL = acos(1.d0/r*sqrt((eemin-m2)/(eemin+m2)))
         
         bmi(1) = nint(theminL*1000d0)
         bma(1) = nint(themaxL*1000d0)
         
         bmi(2) = eemin
         bma(2) = 140.d0        !eemax
         
         bmi(3) = 0.d0 !bmi(1)
         bma(3) = bma(1)/7.d0
         bma(3) = 5.d0

         bmi(7) = -lambda(s,ame*ame,ammu*ammu)/s !bmi(5)
         bma(7) = bma(5)
         bma(7) = 0.d0

         bmi(6) = bmi(5)
         bma(6) = bma(5)

         bmi(9) = bmi(6)
         bma(9) = bma(6)

         
         bmi(8) = -1d-4
         bma(8) = pi - bmi(8)

         bmi(8) = -pi
         bma(8) =  3.d0*pi

         bmi(8) =  -pi*1d3
         bma(8) =  +pi*1d3

         bmi(8) = 0.d0!-pi*1d3
         bma(8) = pi*1d3

c         bmi(8) = -1.d0
c         bma(8) =  1.d0
         
         
c         bmi(8) = -5d0 !bmi(8)*1d-2
c         bma(8) = 5.d0 !bma(8)*1d-2
         
         
         bmi(10) = -1d0 !  0.d0
         bma(10) = 1.d0!-0.95d0!1.d0 !180.d0
         bmi(12) = -1.d0 ! 0.d0
         bma(12) = 1.d0 !180.d0

         bmi(10) = 0.d0
         bma(10) = 180.d0
         bmi(12) = bmi(10)
         bma(12) = bma(10)
         
         bmi(11) = m2
         bma(11) = eemaxcm

c         bma(11) = sqrt(s)*0.5d0
         
         bmi(13) = m1
         bma(13) = emumaxcm

         bmi(14) = bmi(11)+ m1 !bmi(13)
         bma(14) = bma(11)+bma(13)
         bma(14) = bma(14)!*1.15d0
         
         bmi(15) = 0.d0
         bma(15) = 40.d0 !1.25d0*bma(1) !pi * 1d3
         bmi(16) = max(atleasteg,0.005d0*emulab)
         bma(16) = emulab
         bmi(17) = 0.d0
         bma(17) = 180.d0

         bmi(18) = 0.d0 !0.005d0*ecmsnom
         bma(18) = ecmsnom*0.48d0

         bmi(19) = 0.d0
         bma(19) = 2.d0*bma(15) !180.d0
         bmi(19) = -1.d0!0.d0
         bma(19) = 1.d0 !180.d0

         bmi(26) = bmi(17)
         bma(26) = bma(17)
         bmi(27) = bmi(18)
         bma(27) = bma(18)
         
         cmin = max(-1.d0,-1.d0-tmax*0.5d0/pcm/pcm)
         
         bmi(21) = 0.d0
         bma(21) = 180.d0
         bmi(22) = 0.d0 !bmi(21)
         bma(22) = sqrt(s) - ame - ammu!bma(21)
         

         bmi(24) =  bmi(1)
         bma(24) =  bma(1)

         bmi(28) =  0.d0
         bma(28) =  pi * 1d3*0.5d0
         bmi(29) =  0.d0
         bma(29) =  bma(28)
         if (nphotmode.eq.1010) then
            bma(28) = 18.d0
            bma(29) = 18.d0     
         endif
            
c         bmi(1) = bmi(28)
c         bma(1) = bma(28)
         
         bmi(30) =  0.d0
         bma(30) = emulab
         bmi(31) =  0.d0
         bma(31) = emulab
         bmi(32) = 0.d0
         bma(32) = sqrt(s) - ame - ammu
         
         bmi(20) = bmi(32)
         bma(20) = bma(32)

c         bmi(32) = 0.2d0
c         bma(32) = 0.93d0
         
         
         bmi(23) = bmi(2)
         bma(23) = bma(2)

         bmi(33) = 0.d0
         bma(33) = 10.d0

         bmi(34) = 80.d0
         bma(34) = 190.d0
         
c         bmi(33:36) = 0.d0
c         bma(33:36) = 1.d0
         
c         bmi(33) = 50.d0
c         bma(33) = 180.d0
c         bmi(34) = 0d0
c         bma(34) = 10.d0
         
         
c         do k = 1,100
c         print*,'RESET limits for: 2,5,11 '
c         enddo

c         bmi(35) = -3d-6
c         bma(35) = -bmi(35)         
c         bmi(36) = -pi*1d3
c         bma(36) = -bmi(36)
c         bmi(37) = -pi*1d3 /1000.d0
c         bma(37) = -bmi(37)
         
!     filling distr. components when entering the first time
         distr   = 0.d0
         edistr = 0.d0
         do k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               xaxis(k,i)   = y - d
            enddo
         enddo
*****************************
****  2d distributions
         if (i2d.gt.0) then
            bmi2(1,1) = 0.d0
            bmi2(1,2) = 0.d0
            bma2(1,1) = 1.d0
            bma2(1,2) = 1.d0
            bmi2 = 0.d0
            bma2 = 1.d0

            
            dis2d  = 0.d0         
            edis2d = 0.d0
            do k = 1,ndistr2
               dx = (bma2(k,1) - bmi2(k,1))/nbin2
               dy = (bma2(k,2) - bmi2(k,2))/nbin2
               x = bmi2(k,1)
               y = bmi2(k,2)
               do i = 1,nbin2
                  x = x+dx
                  y = y+dy
                  xaxis2(k,i) = x - dx
                  yaxis2(k,i) = y - dy
               enddo  
            enddo
         endif
**********************************         
         ifirst = 0
      endif
      
cc this is more general
c      call new_boost(p12lab,p3,p3lab,-1)
c      call new_boost(p12lab,p4,p4lab,-1)
** the two coincide anyway if p2lab is at rest
      
      p3mcm = sqrt(tridot(p3,p3))
      p4mcm = sqrt(tridot(p4,p4))
      
      p3mlab = sqrt(tridot(p3lab,p3lab))      
      p4mlab = sqrt(tridot(p4lab,p4lab))
      
      fun27 = p4lab(3)/p4mlab
      cccc  = p3lab(3)/p3mlab

      fun(1) = acos(fun27)*radtodeg
      fun(3) = acos(cccc)*radtodeg
      
      fun(2) = p4lab(0)
      fun(4) = p3lab(0)
      
c      t24 = dot(p2lab-p4lab,p2lab-p4lab)
      t24 = dot(p2lab-p4lab,p2lab-p4lab)
      fun(5) = t24
c      t13 = dot(p1lab-p3lab,p1lab-p3lab)
      t13 = dot(p1lab-p3lab,p1lab-p3lab)
      fun(7) = t13

****  extra      
      
ccccccc      fun(8) = testplanarity(p3lab,p4lab)
      fun(8) = acoplanarityPSIdef(p3lab,p4lab)*1d3

c      call acoplanaritiesGA(p3lab,p4lab,fun(35),fun(36),fun(37))
      
      pm4 = sqrt(tridot(p4,p4))
      pm3 = sqrt(tridot(p3,p3))
      fun(10) = bmi(10)-1.d0
      fun(12) = bmi(12)-1.d0
      if (pm4.gt.0.d0) fun(10) = acos(p4(3)/pm4)*180.d0/pi
      if (pm3.gt.0.d0) fun(12) = acos(p3(3)/pm3)*180.d0/pi
      
      fun(11) = p4(0)
      fun(13) = p3(0)
      
      fun(14) = sqrt(dot(p3+p4,p3+p4))
      
      fun(21) = acollinearity(p3,p4)
      fun(22) = -1d3

      fun(33) = acos(truep1lab(3)/sqrt(tridot(truep1lab,truep1lab)))
      fun(33) = fun(33)*1d3
      fun(34) = truep1lab(0)
      
**** photon and pairs variables
      fun(15:19) = bmi(15:19) - 10d0
      fun(24) = fun(1)
      fun(25:27) = bmi(25:27) - 10d0

      qa  = 0.d0                 ! in the lab
      qa2 = 0.d0
      if (qph(1,0).ge.qph(2,0)) then
! ordered in energy!         
         qb  = qph(1,:)
         qb2 = qph(2,:)
      else
         qb2 = qph(1,:)
         qb  = qph(2,:)
      endif      
      if (qb(0).gt.0.d0)  call new_boost(pin2,qb,qa,1)
      if (qb2(0).gt.0.d0) call new_boost(pin2,qb2,qa2,1)
      
      if (qb(0).gt.0.d0) then
         if (qa(0).gt.bmi(16)) then
            fun(15) = radtodeg*acos(qa(3)/qa(0))
            fun(16) = qa(0)            
            fun(19) = tridot(qb,p4)/qb(0)/sqrt(tridot(p4,p4))
         endif
         if (qb(0).gt.bmi(18)) then
            fun(17) = acos(qb(3)/qb(0))*180.d0/pi
            fun(18) = qb(0)
         endif
         
         endif

         if (qb2(0).gt.bmi(27)) then
            fun(26) = acos(qb2(3)/qb2(0))*180.d0/pi
            fun(27) = qb2(0)
         endif
                  
         fun(23) = pinside_e(0)
***********pairs
         fun(28:32) = bmi(28:32) -1d3
         fun(20) = bmi(20) - 1d3
         fun(6) = bmi(6) - 1d3
         fun(9) = bmi(9) - 1d3         
         if (p5(0).gt.0.5d0*ame) then
            fun(24) = acos(pinside_e(3)/
     .           sqrt(tridot(pinside_e,pinside_e)))*1d3
            fun(25) = dot(p2lab - pinside_e,p2lab - pinside_e)
            
c            call new_boost(p12lab,p5,p5lab,-1)
c            call new_boost(p12lab,p6,p6lab,-1)
            fun(28) = acos(p5lab(3)/sqrt(tridot(p5lab,p5lab)))*1d3
            fun(29) = acos(p6lab(3)/sqrt(tridot(p6lab,p6lab)))*1d3
            fun(30) = p5lab(0)
            fun(31) = p6lab(0)

            if (dot(p5,p5).lt.0.5d0*ame*ame) then ! 2 photons from pion decay
               xxx = acos(p5lab(3)/sqrt(tridot(p5lab,p5lab)))*1d3
               yyy = acos(p6lab(3)/sqrt(tridot(p6lab,p6lab)))*1d3
               fun(28) = max(xxx,yyy)
               fun(29) = min(xxx,yyy)
               fun(30) = max(p5lab(0),p6lab(0))
               fun(31) = min(p5lab(0),p6lab(0))
            endif
            
            fun(32) = sqrt(dot(p5+p6,p5+p6))

cc            fun(32) = (p5(0)+p6(0))/sqrt(dot(p12lab,p12lab))

            fun(20) = sqrt(dot(p5+p4,p5+p4))
            fun(22) = sqrt(dot(p4+p5+p6,p4+p5+p6))
            fun(6)  = dot(p2lab-p5lab,p2lab-p5lab)
            fun(9)  = dot(p2lab-p6lab,p2lab-p6lab)
            fun(8)  = acoplanarityPSIdef(p3lab,pinside_e)*1d3

c            call acoplanaritiesGA(p3lab,pinside_e,
c     .           fun(35),fun(36),fun(37))
            
         endif

ccc         fun(32) = xvar(t24)

c         fun(33:36) = csid(1:4)
         
***********************************************************************
         DO k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
** NEW ***
            if (fun(k).le.bma(k).and.fun(k).ge.bmi(k)) then
               i = int((fun(k) - bmi(k))/d)
               i = min(i+1,nbin) ! to avoid roundoffs (?)
               do ii = 1,nw
                  distr(k,i,ii)  = distr(k,i,ii)  + sd(ii)
                  edistr(k,i,ii) = edistr(k,i,ii) + sd(ii)*sd(ii)
               enddo
            endif
** NEW ***
         ENDDO

         IF (i2d.gt.0) then
*** 2D distributions      
c            fun2(1,1) = fun(3)
c            fun2(1,2) = acos(pinside_e(3)/
c     .           sqrt(tridot(pinside_e,pinside_e)))*1d3
            fun2(1,1) = csid(1)
            fun2(1,2) = csid(4)
            fun2(2,1) = csid(3)
            fun2(2,2) = csid(4)
*******************      
            DO k = 1,ndistr2
               dx = (bma2(k,1) - bmi2(k,1))/nbin2
               dy = (bma2(k,2) - bmi2(k,2))/nbin2
** NEW ***
            if ((fun2(k,1).le.bma2(k,1).and.fun2(k,1).ge.bmi2(k,1)).and.
     .         (fun2(k,2).le.bma2(k,2).and.fun2(k,2).ge.bmi2(k,2))) then
               i = int((fun2(k,1) - bmi2(k,1))/dx)
               j = int((fun2(k,2) - bmi2(k,2))/dy)
               i = min(i+1,nbin2) ! to avoid roundoffs (?)
               j = min(j+1,nbin2) ! to avoid roundoffs (?)
            
*** in the 2nd 2d distribution I use the weights without \Delta\alpha_had
** no... 
c            if (k.eq.1) then
               dis2d(k,i,j)  = dis2d(k,i,j)  + sd(1)
               edis2d(k,i,j) = edis2d(k,i,j) + sd(1)*sd(1)
c            else
c               sdw = sd*extraweights(2)
c               dis2d(k,i,j)  = dis2d(k,i,j)  + sdw
c               edis2d(k,i,j) = edis2d(k,i,j) + sdw*sdw
c            endif            
            endif
** NEW ***
         ENDDO
      ENDIF
      return
      end
*********************************************************************
      subroutine writedistributions(realnc,nw)
! nw is the number of weights
#include "distributions_inc.F"
      common/ionlyfirsttime/ifirst,icount
      data ifirst /1/
      common/mueexpsetup/emulab,eemin,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela

      integer*8 nc

      character*200  path
      common/cmnpath/path
            
      character*6 ord
      common/qedORDER/ord

****************
      if (nw.eq.0) return
****************

      if (ifirst.eq.1) return
      
      DO j = 1,nw
      
      do k = 1,ndistr    
         open(12,file=distname(k,j),status='unknown')
         do i = 1, nbin

            av   = 0.d0
            stdd = 0.d0
            stddroot = 0.d0

            if (abs(distr(k,i,j)).gt.0.d0) then
               av  = distr(k,i,j)/realnc
               av2 = edistr(k,i,j)/realnc
               
               stdd     = sqrt(abs(av2 -av*av)/realnc) ! this is the MC error
               stddroot = sqrt(av2/realnc) !! this is equal to root. This is the statistical error with weitghts.
            endif
            
            sezd = av/(bma(k)-bmi(k))*nbin
            esezd = stdd/(bma(k)-bmi(k))*nbin                        
            esezdroot = stddroot/(bma(k)-bmi(k))*nbin                        

            write(12,*)xaxis(k,i),sezd,esezd,esezdroot
         enddo
         close(12)
      enddo

      ENDDO
      
*** 2D distributions      
      IF (i2d.gt.0) THEN
      do k = 1,ndistr2
         open(12,file=distname2(k),status='unknown')
         do i = 1, nbin2
         do j = 1, nbin2

            av   = 0.d0
            stdd = 0.d0
            stddroot = 0.d0

            if (abs(dis2d(k,i,j)).gt.0.d0) then
               av  = dis2d(k,i,j)/realnc
               av2 = edis2d(k,i,j)/realnc
               
               stdd     = sqrt(abs(av2 -av*av)/realnc)
               stddroot = sqrt(av2/realnc)
            endif
            
            sezd = av/(bma2(k,1)-bmi2(k,1))*nbin2
            sezd = sezd/(bma2(k,2)-bmi2(k,2))*nbin2

            esezd = stdd/(bma2(k,1)-bmi2(k,1))*nbin2                        
            esezd = esezd/(bma2(k,2)-bmi2(k,2))*nbin2
            
            esezdroot = stddroot/(bma2(k,1)-bmi2(k,1))*nbin2
            esezdroot = esezdroot/(bma2(k,2)-bmi2(k,2))*nbin2
            
            write(12,*)xaxis2(k,i),yaxis2(k,j),sezd,esezd,esezdroot
         enddo
         write(12,*)' '
         enddo
         close(12)
      enddo
      ENDIF      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function testplanarity(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),ptmp(0:3),pz(0:3),ptmp2(0:3)
      common/zaxis/pz
      data pz /1.d0,0.d0,0.d0,1.d0/
      parameter (pi = 4.d0*atan(1.d0))
      
      call vecprod(p1/sqrt(tridot(p1,p1)),p2/sqrt(tridot(p2,p2)),ptmp)
      tp = tridot(pz,ptmp)
      
cc my definition, or A_2 by Giovanni      
c      call vecprod(p1,p2,ptmp)
c      ptmp = ptmp/sqrt(tridot(ptmp,ptmp))      
c      testplanarity = tridot(pz,ptmp)

cc A_1 by Giovanni
cc the following is identical to acoplanarityPSIdef      

      call vecprod(pz,p1,ptmp)
      ptmp = ptmp/sqrt(tridot(ptmp,ptmp))
      call vecprod(pz,p2,ptmp2)
      ptmp2 = ptmp2/sqrt(tridot(ptmp2,ptmp2))
      
      testplanarity = tridot(ptmp2,ptmp)
      testplanarity = max(-1.d0,testplanarity)      
      testplanarity = (pi - acos(testplanarity))*1d3 ! in mrad
      testplanarity = testplanarity*sign(1.d0,tp)
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine acoplanaritiesGA(p1,p2,trp,gio,cla)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),ptmp(0:3),pz(0:3),ptmp2(0:3)
      common/zaxis/pz
      data pz /1.d0,0.d0,0.d0,1.d0/
      parameter (pi = 4.d0*atan(1.d0))
      
      call vecprod(p1/sqrt(tridot(p1,p1)),p2/sqrt(tridot(p2,p2)),ptmp)
      trp = tridot(pz,ptmp)

cc A_\phi by Giovanni, which is identical to acoplanarityPSIdef
      call vecprod(pz,p1,ptmp)
      ptmp = ptmp/sqrt(tridot(ptmp,ptmp))
      call vecprod(pz,p2,ptmp2)
      ptmp2 = ptmp2/sqrt(tridot(ptmp2,ptmp2))
      
      testplanarity = tridot(ptmp,ptmp2)
      testplanarity = max(-1.d0,testplanarity)      
      testplanarity = (pi - acos(testplanarity))*1d3 ! in mrad
      testplanarity = testplanarity*sign(1.d0,trp)
      gio = testplanarity
cc
      call vecprod(p1/sqrt(tridot(p1,p1)),p2/sqrt(tridot(p2,p2)),ptmp)      
      ptmp = ptmp/sqrt(tridot(ptmp,ptmp))      
      testplanarity = tridot(pz,ptmp)
      cla = pi*0.5d0 - acos(testplanarity)
      cla = cla * 1d3
      return
      end      
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function acoplanarity(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),ptmp(0:3),pz(0:3)
      common/zaxis/pz
      data pz /1.d0,0.d0,0.d0,1.d0/
      parameter (pi = 4.d0*atan(1.d0))

      phi1 = datan2(p1(2),p1(1))
      phi2 = datan2(p2(2),p2(1))

      acoplanarity = pi - phi2 + phi1
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function acoplanarityPSIdef(p1,p2)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      phi1 = datan2(p1(2),p1(1))
      phi2 = datan2(p2(2),p2(1))
      acoplanarityPSIdef = abs(pi - abs(phi2 - phi1))
      return
      end
*****************************************
      function get_dphi(p1,p2)
      implicit double precision (a-h,o-z)
      double precision p1(0:3),p2(0:3)
! from ALPGEN
      dphi = (p1(1)*p2(1)+p1(2)*p2(2))/
     .     sqrt(p1(1)**2+p1(2)**2)/sqrt(p2(1)**2+p2(2)**2)
      if (abs(dphi).gt.1.d0) then
         dphi = dphi/abs(dphi)
      endif
      get_dphi = acos(dphi)
      return
      end
***********************
      function getphi_here(ppp)
      implicit double precision (a-h,o-z)
      dimension ppp(0:3)
      parameter (pi = 3.141592653589793238462643383279502884197169399d0)
      pm = sqrt(ppp(1)**2+ppp(2)**2+ppp(3)**2)
      c  = ppp(3)/pm
      s  = sqrt(1.d0-c**2)
      if (s.eq.0.d0) then
         getphi_here = 0.d0
         return
      else   
         arg = ppp(1)/pm/s
*  avoiding numerical problems......
         if (abs(arg).ge.1.d0) then
            iarg = arg
            arg  = iarg
         endif   
         if (ppp(2).ge.0.d0) getphi_here = acos(arg)
         if (ppp(2).lt.0.d0) getphi_here = 2.d0*pi-acos(arg)
      endif             
      return
      end
****
      subroutine get_Emu(eb,me,mm,thm,the,E1,E2,E)
      implicit double precision (a-h,m,o-z)
c      implicit real*16 (a-h,m,o-z)
      
      c   = cos(thm)      
      me2 = me*me
      mm2 = mm*mm

!!!!!
      if (c*c.lt.(1.d0-me2/mm2)) then
         print*,'NO SOLUTIONS! setting them to 0.'
         E1 = 0.d0
         E2 = 0.d0
         E  = 0.d0
         return
      endif
!!!!       
      
      b = eb+me
      g = sqrt(eb*eb-mm2)*c
      a = mm2+eb*me
      
ccc      print*,'>>>>',(g*g-b*b)*mm2+a*a,abs(b*b-g*g),g,a*b
      if (abs(b*b-g*g).lt.1d-16) then
         E1 = (a*a+g*g*mm2)*0.5d0/a/b
         E2 = E1
      else
         E1 = a*b+abs(g)*sqrt((g*g-b*b)*mm2+a*a)
         E2 = a*b-abs(g)*sqrt((g*g-b*b)*mm2+a*a)
         E1 = E1/(b*b-g*g)
         E2 = E2/(b*b-g*g)
c         E2p = (g*g*mm2 + a*a)/(b*b-g*g)/E1
c         E1p = (g*g*mm2 + a*a)/(b*b-g*g)/E2
      endif

c      print*,'tentativo a capocchia (pmu ~ Emu)'
c      print*,a/(b-g),E1+E2 ! SONO MOLTO SIMILI, interessante...
c      print*,2.d0*b/(b+g) ! ~ 1. , ie b~g....
      
      if (E1.gt.eb) then
         print*,'one solution larger than ',eb
         print*,'setting it to ',eb
         print*,'delta ',abs(eb-E1)
         E1 = eb
      endif
      if (E2.gt.eb) then
         print*,'one solution larger than ',eb
         print*,'delta ',abs(eb-E2)
         print*,'setting it to ',eb
         E2 = eb
      endif
      if (E1.lt.mm) then
         print*,'one solution smaller muon mass '
         print*,'delta ',abs(mm-E1)
         print*,'setting it to muon mass'
         E1 = mm
      endif
      if (E2.lt.mm) then
         print*,'one solution smaller muon mass '
         print*,'delta ',abs(mm-E2)
         print*,'setting it to muon mass'
         E2 = mm
      endif

*** selecting the solution according to th_e      
      r = g/c/b
      ee = b/a*eb*me
      theb = acos(1.d0/r*sqrt((ee-me)/(ee+me)))
      
      if (the.lt.theb) then
         E = E2
      else
         E = E1
      endif

      return
      
*** https://en.wikipedia.org/wiki/Bisection_method      
** first root      
c$$$      xi = mm
c$$$      xs = eb*0.5d0
c$$$
c$$$      xi = 30.d0
c$$$      xs = 36d0
c$$$      uno = 1.d0
c$$$      do k = 1,1000
c$$$         xc = (xs+xi)*0.5d0
c$$$!ERRATA!         fi = xi*xi*b+xi*g+a + sqrt(xi*xi-mm2)*d
c$$$!ERRATA!         fc = xc*xc*b+xc*g+a + sqrt(xc*xc-mm2)*d
c$$$c     if (abs(fc).lt.1d-5) then
c$$$c     E1 = xc
c$$$c         else
c$$$c     endif
c$$$         if (nint(sign(uno,fc)).eq.nint(sign(uno,fi))) then
c$$$            xi = xc
c$$$         else
c$$$            xs = xc
c$$$         endif
c$$$         print*,xc,fc
c$$$      enddo            
c$$$      return
      end
