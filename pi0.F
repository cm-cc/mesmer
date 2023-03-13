*********** PI_0 stuff
      subroutine phasespacemuepi0(s,m1,m2,m3,m4,p1,p2,p3,p4,p5,p6,
     >     np,phsp,w,ier)
! derived from phasespacemuepairsv4
      implicit double precision (a-h,l,m,o-z)
      parameter (nphmax = 2)
      integer np
      parameter (pi = 4.d0*atan(1.d0))
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p5(0:3),p6(0:3)
      dimension q13(0:3),q24(0:3),q13a(0:3),q24a(0:3)
      dimension csi1(1),csi2(2),csi3(3),csi4(4),csi5(5),dir4(0:3)
      dimension pe(0:3),pmu(0:4)
      dimension p1a(0:3),p2a(0:3),p3a(0:3),p4a(0:3)
      dimension p4ka(0:3), p4k(0:3),ps(0:3),ps2(0:3)
      dimension Q(0:3),Qstar(0:3),Qp(0:3),Qpstar(0:3)
      double precision l1,l2
      dimension ptmp(0:3),pbar(0:3),dir3(0:3),ptmp2(0:3)
      double precision kone(0:3),lastk(0:3),nu,kmod,l(0:3),k2
      common/photoncutoff/egmin,egmin2
      common/tmintmax/tmincut,tmaxcut

      common/mueexpsetup/emulab,eemin,eemax,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela

      double precision Qmu
      common/muoncharge/Qmu
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      double precision truep1lab(0:3)
      common/labmomenta/p1lab,p2lab,truep1lab,p12lab
      logical nge2,sol1,sol2,cgt0,clt0
      common/tcutcommon/icut
      common/nfotonicommon/nfotoni
      common/secondphotonmaxenergy/enmaxg2
      common/icount/icount
      common/identicalandperipheral/iper,noper,iid
      dimension twopis(0:nphmax),factorial(0:nphmax)
      common/pispairs_pi0/twopis,factorial,rmax,nmax,nmaxb,ifirst
      data ifirst /0/

      common/xcrf/x,y,modk,ck

      common/pionparameters/pi0m,fpi0,wpi0gg,gpi0
      
*** for distributions      
      double precision csid(8)
      common/csifordistr/csid
      
      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         factorial(0) = 1.d0
         do k = 1,nphmax
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
            factorial(k) = factorial(k-1)*k
         enddo
         rmax = -100d0
         nmax = nphmax            
         ifirst = 1
         icount = 0
      endif
      if (np.ne.10) then
         print*,'Should not reach here! Stopping!'
         stop
      endif      
      m5 = pi0m
      m6 = 0.d0
***********************
      p3 = 0.d0
      p4 = 0.d0
      p5 = 0.d0
      p6 = 0.d0
      ier = 0
      sqs = sqrt(s)      
      if (sqs.lt.(m3+m4+m5+m6)) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      phsp = 1.d0
      w    = 1.d0
      m12  = m1*m1
      m22  = m2*m2
      m32  = m3*m3
      m42  = m4*m4
      m52  = m5*m5
      m62  = m6*m6

      alam = lambda(s,m12,m22)
      
      pm = 0.5d0*sqrt(alam/s)
      
      p1(0) =  sqrt(pm*pm + m12)      
      p1(1) =  0.d0
      p1(2) =  0.d0
      p1(3) =  pm
      
      p2(0) =  sqrt(pm*pm + m22)
      p2(1) =  0.d0
      p2(2) =  0.d0
      p2(3) = -pm

***   depending on the radiation channel, if radiation from mu (ich=2) or e line (ich=1)
      pch1 = 0.d0
      pch2 = 1.d0 - pch1
      ich  = 2                  ! default, t2[46] is anyway sampled...
! verificare se serve      
c      if (np.eq.1) then
c         call getrnd(csi1,1)
c         ich = 1
c         pch1 = 0.005d0
c         if (iid.eq.0) pch1 = 1.d0
c         pch2 = 1.d0 - pch1
c         if (csi1(1).gt.pch1) ich = 2
c      endif
         
      if (ich.eq.1) then
         m3f = m3
         m4f = m4
      else
         m3f = m4
         m4f = m3
      endif
      m3f2 = m3f*m3f
      m4f2 = m4f*m4f
      
****  min and max for (p4+p5+p6)**2
      m4562min = (m4f+m5+m6)*(m4f+m5+m6)
      m4562max = (sqs-m3f)*(sqs-m3f)
****
c flat
      call getrnd(csi1,1)
      an    = m4562max - m4562min
      m4562 = an * csi1(1) + m4562min
      w     = w * an
c  end flat
      
****
** its the energy in the cm of the t13 or t24      
      Es  = 0.5d0*(s + m3f2 - m4562 )/sqs
      psm = sqrt(Es*Es - m3f2)

** t_sampling limits
      if (ich.eq.1) then
         tsmin = m12 + m3f2 - 2.d0*(p1(0)*Es + psm*pm)
         tsmax = m12 + m3f2 - 2.d0*(p1(0)*Es - psm*pm)
      else
         tsmin = m22 + m3f2 - 2.d0*(p2(0)*Es + psm*pm)
         tsmax = m22 + m3f2 - 2.d0*(p2(0)*Es - psm*pm)
      endif

c as 1/t13 THIS SEEMS TO BE BETTER, DEFAULT
      an = log(tsmax/tsmin) 
      call getrnd(csi1,1)
      ts = tsmin*exp(csi1(1)*an)
      w = w * an * ts      
c as 1/t13^2
c      call getrnd(csi1,1)
c      an  = 1.d0/tsmin - 1.d0/tsmax
c      ts = tsmin / (1.d0-tsmin*an*csi1(1))
c      w   = w * an * ts*ts
      
      if (ich.eq.1) then
         cs = ts - m12 - m3f2 + 2.d0*p1(0)*Es
         cs = cs*0.5d0/psm/pm
      else
         cs =  ts - m22 - m3f2 + 2.d0*p2(0)*Es
         cs =  cs*0.5d0/psm/pm
         cs = -cs ! because p2 is along -z
      endif
      if (abs(cs).gt.1.d0) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      ss = sqrt(abs(1.d0-cs*cs))
      
      call getrnd(csi1,1)
      phis = 2.d0*pi*csi1(1)
      w    = w * 2.d0*pi

! this is p3 (ich=1) or p4 (ich=2)
      ps(0) = Es
      ps(1) = psm * ss * sin(phis)
      ps(2) = psm * ss * cos(phis)
      ps(3) = psm * cs
*****
      phsp = phsp * 0.25d0/sqrt(alam)
*****       

! q13 is actually q13 for ich = 1 or q24 for ich = 2      
      if (ich.eq.1) then
         q13 = p1 - ps
         sp  = dot(q13+p2,q13+p2)
      else
         q13 = p2 - ps
         sp  = dot(q13+p1,q13+p1)
      endif
         
***   now I'm in the q13 + p2 RF
      sqsp = sqrt(sp)

c      mk2min = (m5+m6)*(m5+m6)
c      mk2max = (sqsp - m4f)*(sqsp - m4f)
c      if (mk2min.ge.mk2max) then
c         phsp = 0.d0
c         w    = 0.d0
c         ier  = 1
c         return
c      endif
cc flat      
c      call getrnd(csi1,1)      
c      an  = mk2max - mk2min
c      mk2 = an*csi1(1) + mk2min
c      w   = w * an
cc end flat

      mk2 = pi0m*pi0m
      
! this is the other p with respect to ps, ie p4 for ich = 1 or p3 for ich 2, in the q13+p[1 or 2] rest frame
      Es2  = (sp + m4f2 - mk2)*0.5d0/sqsp
      ps2m = sqrt(Es2*Es2-m4f2)

      ek   = (sp - m4f2 + mk2)*0.5d0/sqsp   
      modk = sqrt(ek*ek-mk2)      
***
      phsp = phsp * modk/sqsp * 0.5d0
***      
** now I bring q13 and p2 where q13+p2 is at rest. Notice that here q13a and p2a are always back-to-back.
      if (ich.eq.1) then
         call new_boost(q13+p2,q13,q13a,1)
         call new_boost(q13+p2,p2,p2a,1)
      else
         call new_boost(q13+p1,q13,q13a,1)
         call new_boost(q13+p1,p1,p2a,1) ! notice that if ich = 2 in p2a there is actually p1 boosted
      endif
** DEFAULT ck to be sampled later as 1/((p2a-k)^2-m42) and/or 1/((q13a-k)^2)
      A = mk2-2.d0*p2a(0)*ek
      A = A / 2.d0 /sqrt(tridot(p2a,p2a))/modk
      B = ts + mk2-2.d0*q13a(0)*ek
      B = B / 2.d0 /sqrt(tridot(q13a,q13a))/modk
      call anothercollinear(abs(A),abs(B),1.d0,q13a,p2a,lastk,wcoll)
      w         = w * wcoll
      kone(0)   = ek
      kone(1:3) = lastk(1:3)*modk

      p4a(0)   =  Es2
      p4a(1:3) = -kone(1:3)

******
      q24a = p2a - p4a

*** now I bring q13a and q24a in the kone RF. Notice that here q13a and q24a are always back-to-back.
c      call new_boost(kone,q13a,q13a,1)
c      call new_boost(kone,q24a,q24a,1)

c this is when p5 was the pi0 in the q13+p2(p1) rest frame !      p5 = kone

** I decay the pi0 in the kone rest frame
      call getrnd(csi1,1)
      cphot = 2.d0 * csi1(1) - 1.d0
      sphot = sqrt(abs(1.d0-cphot*cphot))
      call getrnd(csi1,1)
      phph = 2.d0*pi*csi1(1)
! NO!! because it cancels with the rest of the d\Phi_2...      w = w * 4.d0*pi
      p5(0) = sqrt(mk2)*0.5d0
      p5(1) = p5(0)*sin(phph)*sphot
      p5(2) = p5(0)*cos(phph)*sphot
      p5(3) = p5(0)*cphot
      p6(0) = p5(0)
      p6(1:3) = -p5(1:3)
      
      if (ich.eq.1) then
         call new_boost(q13+p2,p4a,p4,-1)
         call new_boost(q13+p2,kone,kone,-1)
c         call new_boost(q13+p2,p5,p5,-1)
c         call new_boost(q13+p2,p6,p6,-1)
         p3 = ps
      else
         call new_boost(q13+p1,p4a,p3,-1)
         call new_boost(q13+p1,kone,kone,-1)
c         call new_boost(q13+p1,p5,p5,-1)
c         call new_boost(q13+p1,p6,p6,-1)
         p4 = ps
      endif
* and bring the photons back in the p1+p2 rest frame
      call new_boost(kone,p5,p5,-1)
      call new_boost(kone,p6,p6,-1)
***
      phsp = phsp * twopis(1)*0.5d0

      return
      end
*********************************************
*********** form ggpi0 factor
      subroutine ggpi0ff(t1,t2,ff,ffnorm)
      implicit double precision (a-h,o-z)
      integer nc
      parameter (nc = 3)
      double precision hv(3),fv(3),fw(3)
      double precision hw(3),api0(3),Fphi(3)
      double precision sv(3),aphiwpi0(3)
      double precision mrho2(3),mphi2(3),mw2(3)
      common/ffggpi0cmn/pi,hv,fv,fw,
     .     hw,api0,Fphi,
     .     sv,aphiwpi0,fpi,sq2,pi2,
     .     mrho2,mphi2,mw2,
     .     ifirst
      data ifirst /0/
      double precision ce(10),my(10)
      common/cecoeff/ce
      double precision pi0m,fpi0,wpi0gg,gpi0
      common/pionparameters/pi0m,fpi0,wpi0gg,gpi0

      D(t,am2) = 1.d0/(t-am2)

      if (ifirst.eq.0) then
         pi = 4.d0*atan(1.d0)
         sq2 = sqrt(2.d0)
         pi2 = pi*pi

         fpi = 0.092388d0
         fpi = fpi0 ! from userinterface

! init         
         hv = 0.d0
         fv = 0.d0
         hw = 1.d0
         sv = 0.d0
         fw       = 1.d0
         fphi     = 1.d0
         api0     = 0.d0
         aphiwpi0 =0.d0
!!
         
         hv(1) = 0.0335d0
         hv(1) = 3.3497294616090301d-2 ! ekhara (9)
         
         fv(1) = 0.2022d0
         fv(1) = 0.20218353880783102d0 ! ekhara (9)

         hv(2) =  0.00184d0
         hv(3) = -0.485d0
         hv(2) =  1.d0 ! any value is ok, because only hv * fv enter the ff
         hv(3) =  1.d0 ! any value is ok, because only hv * fv enter the ff

         fv(2) = -0.0013/hv(2)
         fv(2) = -1.3130539042110201d-3/hv(2) ! ekhara (9)

         hw(1) = 1.160d0
         hw(1) = 1.1596831889368309d0 ! ekhara (9)
         fw(1) = 0.881d0
         fw(1) = 0.88067796038122648d0 ! ekhara (9)
         
         fphi(1) =  0.783d0
         fphi(1) =  0.78285303775172110d0 ! ekhara (9)
         api0(1) = -0.094d0
         api0(1) = -9.3536295969627758d-2 ! ekhara (9)
         api0(2) = -12.04d0
         api0(2) = -12.035890860025670d0 ! ekhara (9)
         api0(3) =  0.08d0
         api0(3) =  7.6455441519080125d-2 ! ekhara (9)
         
         aphiwpi0(1) = -0.23d0
         aphiwpi0(1) = -0.23212353281028397d0 ! ekhara (9)
         
         mrho2(1) = 0.77526d0**2
         mrho2(2) = 1.465d0**2
         mrho2(3) = 1.720d0**2  !  ekhara (9)
         
         mphi2(1) = 1.01946d0**2 ! ekhara (9)
         mphi2(2) = 1.680d0**2
         mphi2(3) = 2.188d0**2   ! ekhara (9)
         
         mw2(1) = 0.78265d0**2  ! ekhara (9)                  
         mw2(2) = 1.425d0**2    ! ekhara (9)
         mw2(3) = 1.670d0**2
         
! la (19)         
         fv3hv3=nc/16.d0/sqrt(2.d0)/pi/pi
         fv3hv3=fv3hv3-(hv(1)*fv(1)*(1.d0+fw(1)*hw(1)+api0(1)*fphi(1)))
         fv3hv3=fv3hv3-(hv(2)*fv(2)*(1.d0+fw(2)*hw(2)+api0(2)*fphi(2)))

         fv3hv3=fv3hv3/(1+fw(3)*hw(3)+api0(3)*fphi(3))
         
         fv(3)=fv3hv3/hv(3)

! la (20)         
         sv1fv12=sqrt(2.d0)*fv(1)*hv(1)
         sv2fv22=sqrt(2.d0)*fv(2)*hv(2)
         sv3fv32=sqrt(2.d0)*fv(3)*hv(3)

         sv(1) = sv1fv12/fv(1)**2 
         sv(2) = sv2fv22/fv(2)**2 
         sv(3) = sv3fv32/fv(3)**2 
         
         ifirst = 1
      endif

      ff = -nc / 12.d0/fpi/pi2
      do i = 1,3
         ft1 = D(t1,mrho2(i)) + fw(i)*hw(i)*D(t1,mw2(i))
     .        +api0(i)*fphi(i)*D(t1,mphi2(i))
         ft1 = t1*ft1

         ft2 = D(t2,mrho2(i)) + fw(i)*hw(i)*D(t2,mw2(i))
     .        +api0(i)*fphi(i)*D(t2,mphi2(i))
         ft2 = t2*ft2

         ff = ff + 4.d0*sq2*hv(i)*fv(i)/3.d0/fpi*(ft1+ft2)

         ft12 = D(t2,mrho2(i))*D(t1,mw2(i)) +
     .          D(t1,mrho2(i))*D(t2,mw2(i)) +
     .  (api0(i)*fphi(i)-aphiwpi0(i))  * D(t1,mphi2(i))*D(t2,mphi2(i)) +
     .  (fw(i)*hw(i)-1.d0-aphiwpi0(i)) * D(t1,mw2(i))*D(t2,mw2(i)) +
     .        aphiwpi0(i)*
     .      (D(t1,mphi2(i))*D(t2,mw2(i)) + D(t2,mphi2(i))*D(t1,mw2(i))) 
         ft12 = t1*t2*ft12

         ff = ff - 4.d0*sv(i)*fv(i)**2/3.d0/fpi* ft12
      enddo

      ffnorm = ff/ (-nc) * 12.d0*fpi*pi2
      
      return
      end
********************************************************************
* matrix element for the reaction mu- e- --> to mu- e- pi0
*
      subroutine elmatpi0(p1,p2,p3,p4,p5,formfactor,elmatsq)
*
*     p1 is the incoming muon
*     p2 is the incoming electron
*     p3 is the muon momentum
*     p4 is the electron momentum
*     p5 is the pion momentum
*
*      implicit none
c      implicit real*16 (Z)
      implicit double precision (a-z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      COMMON/CONST/EMUBEAM,EEBEAM,ALPHAa,CONVFACa,PIa,AMEa,AMUa
      common/pionparameters/pi0m,fpi0,wpi0gg,gpi0
      common/parameters/ame,amu,convfac,alpha,pi
      common/pi0cmn/ampi0,ampi02
      integer k
      double precision pezzo(16)      
      double precision p4m2(0:3),p3m1(0:3)      
      real*16 pezzoq(16),elmqcmn,constcmn,diff47q,diff1013q
      common/piecesquad/pezzoq,elmqcmn,constcmn,diff47q,diff1013q
      real*16 qdot
      external qdot
*
      ampi0 = pi0m
      ampi02 = pi0m*pi0m

      me2= ame*ame
      me4= me2*me2
      mm2= amu*amu
      mm4= mm2*mm2
      
      p1_p2 = dot(p1,p2)
      p1_p3 = dot(p1,p3)
      p1_p4 = dot(p1,p4)
      p2_p3 = dot(p2,p3)
      p2_p4 = dot(p2,p4)
      p3_p4 = dot(p3,p4)

      t13 = dot(p3-p1,p3-p1)
      t24 = dot(p4-p2,p4-p2)
      
! questo per fattore di forma costante (tipo Brodsky)
c      zff = -1.d0/4.d0/pi/pi/fpi0 ! Czyz ff at t1=t2=0
c      constant = (4*pi*alpha)**2/t13**2/t24**2  !mu-e diagram part
c     +        *(4*pi*alpha)**2*zff**2    !feynman rule from czyz
      
      constant = (4*pi*alpha)**2/t13**2/t24**2  !mu-e diagram part
     +        *(4*pi*alpha)**2*formfactor**2    !feynman rule from czyz

***** BEST ONE - pi0/matrix-element-mine.frm
      a = p1_p2
      b = p1_p3
      c = p1_p4
      d = p2_p3
      e = p2_p4
      f = p3_p4
      x = (a+f)**2 - b*e
      y = (c+d)**2 - b*e
      z = -2.d0*((a*f-c*d)**2)
      xpy = x+y

      Z1_= - 2*e**2
      Z2_=2*f
      Z3_=c - d
      Z4_=Z2_ - Z3_
      Z4_=a*Z4_
      Z5_=c*d
      Z6_= - xpy + 2*Z5_
      Z7_=Z3_*f
      Z4_=Z4_ + Z7_ + Z6_
      Z4_=2*e*Z4_
      Z8_= - 2*b**2
      Z9_=12*b
      Z3_=Z2_ + Z3_
      Z3_=a*Z3_
      Z3_=Z3_ - Z7_ + Z6_
      Z3_=2*b*Z3_
      Z6_=12*e
      Z7_=c + d
      Z2_= - Z2_ - Z7_
      Z2_=a*Z2_
      Z7_=f*Z7_
      Z2_=Z7_ - Z2_
      Z7_=b*e
      Z2_= - 6*Z7_ + 5*xpy - 8*Z5_ - 4*Z2_
      Z5_=xpy*Z7_
      Z5_=z + Z5_
      elmatsq=mm4*Z1_ + mm2*Z4_ + me4*Z8_ - 10*me4*mm4 + me4*mm2*Z9_ + 
     & me2*Z3_ + me2*mm4*Z6_ + me2*mm2*Z2_ + Z5_
      
      elmatsq = 32*elmatsq*constant

c      pezzo(1) = mm4*Z1_
c      pezzo(2) =  + mm2*Z4_
c      pezzo(3) =  + me4*Z8_
c      pezzo(4) =  - 10*me4*mm4
c      pezzo(5) =  + me4*mm2*Z9_
c      pezzo(6) =  + me2*Z3_
c      pezzo(7) =  + me2*mm4*Z6_
c      pezzo(8) =  + me2*mm2*Z2_
c      pezzo(9) =  Z5_     
c      call elmatpi0quad(p1,p2,p3,p4,p5,formfactor,A8)
c      r = elmatsq / A8
c      if (abs(r-1.d0).gt.5d-8) then
c         print*,' '
c         print*,r-1.d0,elmatsq/32/constant-A8/32/constant
c         print*,(1q0*p1_p2-qdot(p1,p2))/qdot(p1,p2)
c         print*,(1q0*p1_p3-qdot(p1,p3))/qdot(p1,p3)
c         print*,(1q0*p1_p4-qdot(p1,p4))/qdot(p1,p4)
c         print*,(1q0*p2_p3-qdot(p2,p3))/qdot(p2,p3)
c         print*,(1q0*p2_p4-qdot(p2,p4))/qdot(p2,p4)
c         print*,(1q0*p3_p4-qdot(p3,p4))/qdot(p3,p4)
c         print*,sum(pezzo(1:9)),sum(pezzoq(1:9))
c         do k = 1,9
c            print*,k,pezzo(k),dble(pezzoq(k)),pezzo(k)-pezzoq(k)
c         enddo
c      endif     
      return

*** with p5 and mpi0: this seems to be much better
c      mp  = pi0m
c      mp2 = mp*mp
c      mp4 = mp2*mp2
c      p1_p5 = qdot(p1,p5)
c      p2_p5 = qdot(p2,p5)
c      p3_p5 = qdot(p3,p5)
c      p4_p5 = qdot(p4,p5)
c$$$      Z1_=p1_p2
c$$$      Z2_=p1_p3
c$$$      Z3_=p1_p5
c$$$      Z4_=p2_p5
c$$$      Z5_=p2_p3
c$$$      Z6_=p3_p5
c$$$      Z7_=Z1_*Z4_
c$$$      Z8_=Z5_*Z4_
c$$$      Z9_= - Z7_ + 3*Z8_
c$$$      Z9_=mm2*Z9_
c$$$      Z10_=2*Z1_
c$$$      Z11_= - Z5_*Z10_
c$$$      Z12_=Z4_ + Z5_ - Z1_
c$$$      Z13_=me2 + Z12_
c$$$      Z14_=mm2*Z13_
c$$$      Z11_=Z11_ + Z14_
c$$$      Z11_=Z3_*Z11_
c$$$      Z14_=Z7_ + Z8_
c$$$      Z15_=Z3_*me2
c$$$      Z15_=Z15_ - Z14_
c$$$      Z15_=Z2_*Z15_
c$$$      Z9_=Z15_ + Z9_ + Z11_
c$$$      Z11_=3*me2 + 2*Z12_
c$$$      Z11_=Z11_*mm2
c$$$      Z15_=Z1_**2
c$$$      Z16_=Z2_*Z13_
c$$$      Z16_=Z16_ + 2*Z15_ - Z11_
c$$$      Z16_=Z6_*Z16_
c$$$      Z9_=2*Z9_ + Z16_
c$$$      Z9_=Z6_*Z9_
c$$$      Z16_=4*mm2
c$$$      Z17_= - Z13_*Z16_
c$$$      Z18_=Z10_ + Z5_
c$$$      Z18_=Z5_*Z18_
c$$$      Z19_=Z2_*Z12_
c$$$      Z17_=Z19_ + Z17_ + Z15_ + Z18_
c$$$      Z17_=Z2_*Z17_
c$$$      Z10_=Z10_ - 3*Z5_
c$$$      Z10_=Z5_*Z10_
c$$$      Z10_= - 3*Z15_ + Z10_
c$$$      Z10_=mm2*Z10_
c$$$      Z12_=3*Z12_ + 4*me2
c$$$      Z12_=mm4*Z12_
c$$$      Z10_=Z12_ + Z10_ + Z17_
c$$$      Z10_=mp2*Z10_
c$$$      Z7_=3*Z7_ - Z8_
c$$$      Z7_=mm2*Z7_
c$$$      Z8_=Z5_**2
c$$$      Z8_=2*Z8_ - Z11_
c$$$      Z8_=Z3_*Z8_
c$$$      Z7_=2*Z7_ + Z8_
c$$$      Z7_=Z3_*Z7_
c$$$      Z8_=Z4_**2
c$$$      Z11_=Z8_*Z16_
c$$$      Z12_=Z3_*Z13_
c$$$      Z12_= - 2*Z14_ + Z12_
c$$$      Z12_=Z3_*Z12_
c$$$      Z11_=Z11_ + Z12_
c$$$      Z11_=Z2_*Z11_
c$$$      Z8_=mm4*Z8_
c$$$      Z7_=Z10_ + Z9_ - 4*Z8_ + Z7_ + Z11_
c$$$      B=32*Z7_
c$$$      B = -B*constant
c$$$      elmatsq = B
c$$$      return
      end
**********************************************************************
      subroutine elmatpi0quad(p1,p2,p3,p4,p5,formfactor,elmatquad)
*     implicit none
      implicit real*16 (Z)
      real*16 B,C
      real*16 a,d,e,f,x,y,z,xpy
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      real*16 elmatsq
      double precision formfactor
      real*16 elmatq
      real*16 ampi0,ampi02
      double precision pi0m,fpi0,wpi0gg,gpi0
      common/pionparameters/pi0m,fpi0,wpi0gg,gpi0
      double precision ame,amu,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/pi0cmnquad/ampi0,ampi02
      real*16 me2,me4,mm2,mm4
      real*16 p1_p2,p1_p3,p1_p4,p2_p3,p2_p4,p3_p4
      real*16 t13,t24
      real*8 dot,elmatquad
      external dot
      real*16 fpi,constant
      real*16 piwidth
      real*16 qdot
      external qdot

      real*16 pezzoq(16),elmqcmn,constcmn,diff47q,diff1013q
      common/piecesquad/pezzoq,elmqcmn,constcmn,diff47q,diff1013q
*
      ampi0 = pi0m
      ampi02 = pi0m*pi0m

      me2= ame*ame
      me4= me2*me2
      mm2= amu*amu
      mm4= mm2*mm2

      p1_p2= qdot(p1,p2)
      p1_p3= qdot(p1,p3)
      p1_p4= qdot(p1,p4)
      p2_p3= qdot(p2,p3)
      p2_p4= qdot(p2,p4)
      p3_p4= qdot(p3,p4)

c      t24= 2.q0*me2-2.q0*p2_p4
c      t13= 2.q0*mm2-2.q0*p1_p3
      t13 = qdot(p3-p1,p3-p1)
      t24 = qdot(p4-p2,p4-p2)

! questo per fattore di forma costante (tipo Brodsky)
c      zff = -1.d0/4.d0/pi/pi/fpi0 ! Czyz ff at t1=t2=0
c      constant = (4*pi*alpha)**2/t13**2/t24**2  !mu-e diagram part
c     +        *(4*pi*alpha)**2*zff**2    !feynman rule from czyz
      
      constant = (4*pi*alpha)**2/t13**2/t24**2  !mu-e diagram part
     +        *(4*pi*alpha)**2*formfactor**2    !feynman rule from czyz

***** BEST ONE - pi0/matrix-element-mine.frm
      a = p1_p2
      b = p1_p3
      c = p1_p4
      d = p2_p3
      e = p2_p4
      f = p3_p4
      x = (a+f)**2 - b*e
      y = (c+d)**2 - b*e
      z = -2.d0*((a*f-c*d)**2)
      xpy = x+y

      Z1_= - 2*e**2
      Z2_=2*f
      Z3_=c - d
      Z4_=Z2_ - Z3_
      Z4_=a*Z4_
      Z5_=c*d
      Z6_= - xpy + 2*Z5_
      Z7_=Z3_*f
      Z4_=Z4_ + Z7_ + Z6_
      Z4_=2*e*Z4_
      Z8_= - 2*b**2
      Z9_=12*b
      Z3_=Z2_ + Z3_
      Z3_=a*Z3_
      Z3_=Z3_ - Z7_ + Z6_
      Z3_=2*b*Z3_
      Z6_=12*e
      Z7_=c + d
      Z2_= - Z2_ - Z7_
      Z2_=a*Z2_
      Z7_=f*Z7_
      Z2_=Z7_ - Z2_
      Z7_=b*e
      Z2_= - 6*Z7_ + 5*xpy - 8*Z5_ - 4*Z2_
      Z5_=xpy*Z7_
      Z5_=z + Z5_
      elmatsq=mm4*Z1_ + mm2*Z4_ + me4*Z8_ - 10*me4*mm4 + me4*mm2*Z9_ + 
     & me2*Z3_ + me2*mm4*Z6_ + me2*mm2*Z2_ + Z5_

      elmatquad = 32*elmatsq*constant

c      pezzoq(1) = mm4*Z1_
c      pezzoq(2) =  + mm2*Z4_
c      pezzoq(3) =  + me4*Z8_
c      pezzoq(4) =  - 10*me4*mm4
c      pezzoq(5) =  + me4*mm2*Z9_
c      pezzoq(6) =  + me2*Z3_
c      pezzoq(7) =  + me2*mm4*Z6_
c      pezzoq(8) =  + me2*mm2*Z2_
c      pezzoq(9) =  Z5_
      return
      end

********************
      subroutine elmatczyz(p1,p2,p3,p4,p5,formfactor,elmatsq)
      implicit none
      real*8 B,C
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      real*8 elmatsq,formfactor
      REAL*8 EMUBEAM,EEBEAM,ALPHAa,CONVFACa,PIa,AMEa,AMUa
      COMMON/CONST/EMUBEAM,EEBEAM,ALPHAa,CONVFACa,PIa,AMEa,AMUa
      real*8 ampi0,ampi02
      double precision pi0m,fpi0,wpi0gg,gpi0
      common/pionparameters/pi0m,fpi0,wpi0gg,gpi0
      double precision ame,amu,convfac,alpha,pi
      common/parameters/ame,amu,convfac,alpha,pi
      common/pi0cmn/ampi0,ampi02
      real*8 me2,me4,mm2,mm4
      real*8 p1_p2,p1_p3,p1_p4,p2_p3,p2_p4,p3_p4
      real*8 p1_p5,p2_p5,p3_p5,p4_p5,mp,mp2,mp4
      real*8 t13,t24,s12,s34,u14,u23
      real*8 dot
      external dot
      real*8 fpi,constant,btmp
      real*8 piwidth,eltmp

      double precision  k(0:3),g1(0:3),l1(0:3),l2(0:3)
      double precision q1(0:3),q2(0:3)
      
      double precision it1,it2
      
      double precision p4m2(0:3),p3m1(0:3),zff
      double precision p1_p3m1,p1_p4m2,p2_p3m1,p2_p4m2,p3m1_p4m2

      real*16 pezzoq(16),elmqcmn,constcmn,diff47q,diff1013q
      common/piecesquad/pezzoq,elmqcmn,constcmn,diff47q,diff1013q

      double precision pezzo(16)

      double precision  DT_kl1, DT_kl2, 
     &         DT_k1k1, DT_k1l1, DT_k1l2, 
     &         DT_q1k, DT_q1k1, DT_q1q2,  DT_q1l1, DT_q1l2,
     &         DT_q2k, DT_q2k1,  DT_q2l1, DT_q2l2,
     &         DT_p1k, DT_p1k1, DT_p1q1, DT_p1q2, DT_p1p2,
     &         DT_p1l1, DT_p1l2,
     &         DT_p2k,  DT_p2k1, DT_p2q1, DT_p2q2,
     &         DT_p2l1, DT_p2l2,
     &         DT_l1l1, DT_l1l2,  DT_l2l2, MsqT

      me2= ame*ame
      me4= me2*me2
      mm2= amu*amu
      mm4= mm2*mm2

      mp  = pi0m
      mp2 = mp*mp
      mp4 = mp2*mp2

      q1 = p3
      q2 = p4

      l1 = p1-q1
      l2 = p2-q2

      it1 = dot(l1,l1)
      it2 = dot(l2,l2)
            
      DT_q1q2 = dot(q1,q2)
      DT_q1l1 = dot(q1,l1)
      DT_q1l2 = dot(q1,l2)
      
      DT_q2l1 = dot(q2,l1)
      DT_q2l2 = dot(q2,l2)
      
      DT_p1q1 = dot(p1,q1)
      DT_p1q2 = dot(p1,q2)
      DT_p1p2 = dot(p1,p2)
      DT_p1l1 = dot(p1,l1)
      DT_p1l2 = dot(p1,l2)
      
      DT_p2q1 = dot(p2,q1)
      DT_p2q2 = dot(p2,q2)
      DT_p2l1 = dot(p2,l1)
      DT_p2l2 = dot(p2,l2)
      
      DT_l1l1 = it1
      DT_l1l2 = dot(l1,l2)
      DT_l2l2 = it2

      MsqT = 0.d0
      MsqT = -( 
     &     - 32.Q0*DT_p1p2*DT_q1q2*DT_l1l1*DT_l2l2 + 32.Q0*
     &     DT_p1p2*DT_q1q2*DT_l1l2*DT_l1l2 + 32.Q0*DT_p1p2*
     &     DT_q1l1*DT_q2l1*DT_l2l2 - 32.Q0*DT_p1p2*DT_q1l1*
     & DT_q2l2*DT_l1l2 - 32.Q0*DT_p1p2*DT_q1l2*DT_q2l1*
     & DT_l1l2 + 32.Q0*DT_p1p2*DT_q1l2*DT_q2l2*DT_l1l1
     &  + 32.Q0*DT_p1q1*DT_p2q2*DT_l1l1*DT_l2l2 - 32.Q0*
     & DT_p1q1*DT_p2q2*DT_l1l2*DT_l1l2 - 32.Q0*DT_p1q1*
     & DT_p2l1*DT_q2l1*DT_l2l2 + 32.Q0*DT_p1q1*DT_p2l1*
     & DT_q2l2*DT_l1l2 + 32.Q0*DT_p1q1*DT_p2l2*DT_q2l1*
     & DT_l1l2 - 32.Q0*DT_p1q1*DT_p2l2*DT_q2l2*DT_l1l1
     &  - 32.Q0*DT_p1q2*DT_p2q1*DT_l1l1*DT_l2l2 + 32.Q0*
     & DT_p1q2*DT_p2q1*DT_l1l2*DT_l1l2 + 32.Q0*DT_p1q2*
     & DT_p2l1*DT_q1l1*DT_l2l2 - 32.Q0*DT_p1q2*DT_p2l1*
     & DT_q1l2*DT_l1l2 - 32.Q0*DT_p1q2*DT_p2l2*DT_q1l1*
     & DT_l1l2 + 32.Q0*DT_p1q2*DT_p2l2*DT_q1l2*DT_l1l1
     &  + 32.Q0*DT_p1l1*DT_p2q1*DT_q2l1*DT_l2l2 - 32.Q0*
     & DT_p1l1*DT_p2q1*DT_q2l2*DT_l1l2 - 32.Q0*DT_p1l1*
     & DT_p2q2*DT_q1l1*DT_l2l2 + 32.Q0*DT_p1l1*DT_p2q2*
     & DT_q1l2*DT_l1l2 + 32.Q0*DT_p1l1*DT_p2l1*DT_q1q2*
     & DT_l2l2 - 32.Q0*DT_p1l1*DT_p2l1*DT_q1l2*DT_q2l2
     &  - 32.Q0*DT_p1l1*DT_p2l2*DT_q1q2*DT_l1l2 + 64.Q0*
     & DT_p1l1*DT_p2l2*DT_q1l1*DT_q2l2 - 32.Q0*DT_p1l1*
     & DT_p2l2*DT_q1l2*DT_q2l1 - 32.Q0*DT_p1l1*DT_q1l1*
     & DT_l2l2*me2 + 32.Q0*DT_p1l1*DT_q1l2*DT_l1l2*me2
     &  - 32.Q0*DT_p1l2*DT_p2q1*DT_q2l1*DT_l1l2 + 32.Q0*
     & DT_p1l2*DT_p2q1*DT_q2l2*DT_l1l1 + 32.Q0*DT_p1l2*
     & DT_p2q2*DT_q1l1*DT_l1l2 - 32.Q0*DT_p1l2*DT_p2q2*
     & DT_q1l2*DT_l1l1 - 32.Q0*DT_p1l2*DT_p2l1*DT_q1q2*
     & DT_l1l2 - 32.Q0*DT_p1l2*DT_p2l1*DT_q1l1*DT_q2l2
     &  + 64.Q0*DT_p1l2*DT_p2l1*DT_q1l2*DT_q2l1 + 32.Q0*
     & DT_p1l2*DT_p2l2*DT_q1q2*DT_l1l1 - 32.Q0*DT_p1l2*
     & DT_p2l2*DT_q1l1*DT_q2l1 + 32.Q0*DT_p1l2*DT_q1l1*
     & DT_l1l2*me2 - 32.Q0*DT_p1l2*DT_q1l2*DT_l1l1*me2
     &  - 32.Q0*DT_p2l1*DT_q2l1*DT_l2l2*me2 + 32.Q0*
     & DT_p2l1*DT_q2l2*DT_l1l2*me2 + 32.Q0*DT_p2l2*
     & DT_q2l1*DT_l1l2*me2 - 32.Q0*DT_p2l2*DT_q2l2*
     & DT_l1l1*me2 + 32.Q0*DT_l1l1*DT_l2l2*me2**2
     &  - 32.Q0*DT_l1l2*DT_l1l2*me2**2
     &           )
      MsqT =  MsqT / (it1 * it2)**2
! now account for the transition form factor
      MsqT =  MsqT * abs(formfactor)**2

      MsqT = MsqT * (4.d0*pi*alpha)**4
      
      elmatsq = MsqT
      return
      end
