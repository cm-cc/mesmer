      subroutine getomax(Einfo,Esupo,m,ua,ub,uc,ud,ue,Emax,omax)
*** https://en.wikipedia.org/wiki/Golden-section_search
      implicit double precision (a-h,m,o-z)
      parameter (grm1=1.d0/((sqrt(5.d0)+1.d0)*0.5d0))
      parameter (eps = 1d-5)
      a = Einfo
      b = Esupo
      c = b - (b-a)*grm1
      d = a + (b-a)*grm1
c      ic = 0
      fstart = ofun(a,m,ua,ub,uc,ud,ue)
      fc = ofun(c,m,ua,ub,uc,ud,ue)
      fd = ofun(d,m,ua,ub,uc,ud,ue)
      do while (abs(fc-fd).gt.eps*fstart)
c         ic = ic + 1
         if (fc.gt.fd) then
            b = d
         else
            a = c
         endif
         c = b - (b-a)*grm1
         d = a + (b-a)*grm1
         fc = ofun(c,m,ua,ub,uc,ud,ue)
         fd = ofun(d,m,ua,ub,uc,ud,ue)
      enddo
      Emax = (a+b)*0.5d0
      omax = ofun(Emax,m,ua,ub,uc,ud,ue) * (1.d0+1.50*eps)
      return
      end
      function ofun(E,m,ua,ub,uc,ud,ue)
      implicit double precision (a-h,m,o-z)
      p    = sqrt(abs(E*E-m*m))
      ofun = (ua+ub*E+uc*p)/(ud-E+ue*p)*0.5d0
      return
      end
***************************************************************************      
      subroutine phasespacemue_v2(s,m1,m2,m3,m4,nph,p1,p2,p3,p4,qphext,
     >     phsp,w,ier)
** in this one I sample c4 and E4. It should be generalized easily to n-gammas      
      implicit double precision (a-h,l,m,o-z)
      parameter (nphmax = 40)
      parameter (pi = 4.d0*atan(1.d0))
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),v3(0:3),q3(0:3)
      dimension csi1(1),csi2(2),csi3(3),csi4(4),csi5(5),dir4(0:3)
      dimension qph(0:3,max(1,nph)),qphext(nphmax,0:3),pe(0:3),pmu(0:4)
      dimension ptmp(0:3),pbar(0:3),dir3(0:3)
      double precision kone(0:3),lastk(0:3),nu,kmod,l(0:3),k2
      common/photoncutoff/egmin,egmin2
      common/s3eps/epss3_ch1,epss3_ch2,sbefore,egminbefore
      data sbefore,egminbefore /1d12,1d12/
      common/tmintmax/tmincut,tmaxcut
      integer*8 count
      common/kappa/count
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab

      logical nge2,sol1,sol2
      
      common/tcutcommon/icut

      common/nfotonicommon/nfotoni

      common/secondphotonmaxenergy/enmaxg2

      common/icount/icount
      
      double precision Erange(2,2)

      common/mueexpsetup/emulab,eemin,semu,thmumin,themax,thmumax,
     .     dthna7max,cutela,ina7,iela
      
      dimension twopis(0:nphmax),factorial(0:nphmax)
      common/pis_new/twopis,factorial,rmax,nmax,nmaxb,ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
cc         print*,'**factorial for photons in the phsp!**'
         twopis(0) = (2.d0*pi)**(4-6)
         factorial(0) = 1.d0
         do k = 1,nphmax
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
            factorial(k) = factorial(k-1)*k
         enddo
         rmax = -100d0
         nmax = nphmax

         enmaxg2 = -1.d0
         
         ifirst = 1
         icount = 0
      endif

      nfotoni = nph
      nmaxb   = nmax
      
      p3 = 0.d0
      p4 = 0.d0
      qph(:,1:max(1,nph))      = 0.d0
      qphext(1:max(1,nmaxb),:) = 0.d0
      nmax = nph
      
      ier = 0

      sqs = sqrt(s)      

      if (sqs.lt.(m3+m4+egmin*(nph-1)+egmin2)) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif

      egminl = egmin

      phsp = 1.d0
      w    = 1.d0

**      using  two channels, for mu- or e-line radiation (for nph > 0)
** no, mu radiation is really negligible

      m12 = m1*m1
      m22 = m2*m2
      
      pm = 0.5d0*sqrt(lambda(s,m12,m22)/s)
     
      p1(0) =  sqrt(pm*pm + m12)
      p1(1) =  0.d0
      p1(2) =  0.d0
      p1(3) =  pm
      
      p2(0) =  sqrt(pm*pm + m22)
      p2(1) =  0.d0
      p2(2) =  0.d0
      p2(3) = -pm

c      sfreeze = m12 + m22 + 2.d0*m2*150d0
c      pmfreeze = 0.5d0*sqrt(lambda(sfreeze,m12,m22)/sfreeze)
      
      if (nph.eq.0) then
** 2 body phasepace
*** non la faccio tanto lunga assumendo 4 masse diverse...         
         e3 = sqrt(pm*pm + m12)
         e4 = sqrt(pm*pm + m22)

         tmax = 2.d0*m22 - 2.d0*e4*e4 + 2.d0*pm*pm
         tmin = 2.d0*m22 - 2.d0*e4*e4 - 2.d0*pm*pm
         tmax = min(tmax,tmaxcut)
         tmin = max(tmin,tmincut)

***** let's see if themax allows for a stricter range in the 2 -> 2
         r  = sqrt(emulab**2-m3**2)/(emulab+m2)
         r2 = r*r
         ce = cos(themax*1d-3)
         ce2 = ce*ce
         Eeminfromthemax = m4*(1.d0+r2*ce2)/(1.d0-r2*ce2)
         tmaxfromthemax  = 2.d0*m4*m4-2.d0*m4*Eeminfromthemax
         tmax = min(tmax,tmaxfromthemax)
         if (tmin.gt.tmax) then
            phsp = 0.d0
            w    = 0.d0
            ier  = 1
            return
         endif
************         
         
         if (thmumin.gt.1d-12) then         
**** [[ implementing the thmumin cut on t13
** it's from the condition |pmu_cm|*sin(thmucm) = |pmu_lab|*sin(thmulab)
** and |pmu_lab|*cos(thmulab) = gb * (|pmu_cm|*cos(thmucm) - bb*emucm)      
            tgl = tan(thmumin*1d-3)
            bb  = -p12lab(3)/p12lab(0) ! boost beta
            gb  = 1.d0/sqrt(1.d0-bb*bb) ! boost gamma
            bcm = pm/E3         ! muon beta in the cm
            
            A =  gb*gb*tgl*tgl + 1.d0
            B = -gb*gb*tgl*tgl * bb/bcm
            C =  gb*gb*tgl*tgl*bb/bcm*bb/bcm - 1.d0
            
            uu1 = -B + sqrt(B*B-A*C)
            uu1 = uu1 /A
            uu2 = -B - sqrt(B*B-A*C)
            uu2 = uu2 /A
! here in the 2 --> t24 = t13
            t13max = 2.d0*m12 - 2.d0*(p1(0)*E3 - uu1*pm*pm)

c            print*,t13max,tmax
c            stop
            
            tmax = min(tmax,t13max)
            if (tmin.gt.tmax) then
               phsp = 0.d0
               w    = 0.d0
               ier  = 1
               return
            endif
***  ]]
         endif
      
         c4max = min( 1.d0,-1.d0-tmin*0.5d0/pm/pm)
         c4min = max(-1.d0,-1.d0-tmax*0.5d0/pm/pm)

         call get_cth(c4,c4min,c4max,s,pm,m12,m22,wc)         
         w  = w * wc
         s4 = sqrt(1.d0-c4*c4)

         call getrnd(csi1,1)
         phi = 2.d0*pi*csi1(1)
         w   = w * 2.d0*pi

         p4(0) = e4
         p4(1) = pm * sin(phi)*s4
         p4(2) = pm * cos(phi)*s4
         p4(3) = pm * c4

         p3   = p1+p2-p4
         
         phsp = 0.25d0*phsp*twopis(nph)*pm/sqs
         return
      endif

      Erange = -1.d0
      
      nge2=.FALSE.
      if (nph.ge.2) nge2=.TRUE.

      qph(:,1:nph)    = 0.d0
      qphext(1:nph,:) = 0.d0
      kone            = 0.d0
      
      e3 = sqrt(pm*pm + m12)
      e4 = sqrt(pm*pm + m22)
      
      tmax =  0.d0 !2.d0*m22 - 2.d0*e4*e4 + 2.d0*pm*pm
      tmin = -4.d0*pm*pm
      tmax = min(tmax,tmaxcut)
      tmin = max(tmin,tmincut)
         
      c4max = min( 1.d0,-1.d0-tmin*0.5d0/pm/pm)
      c4min = max(-1.d0,-1.d0-tmax*0.5d0/pm/pm)
               
      al    = m22 - tmax*0.5d0
      be    = p2(0)

      ilarge = 0
      if (al.le.be*m2) then
         ilarge = 1  ! it must go down to -1... there is below the equation ga*p > al - be*E,
                     ! which can be true for any c4 if p = 0
      endif
      call get_cth_larger(c4,c4min,c4max,s,pm,m12,m22,ilarge,wc)!attenzione, vedi considerazioni sotto

      w = w * wc
      s4 = sqrt(1.d0-c4*c4)

      call getrnd(csi1,1)
      phi = 2.d0*pi*csi1(1)
      w   = w * 2.d0*pi

      dir4(0) = p2(0)
      dir4(1) = pm * sin(phi)*s4
      dir4(2) = pm * cos(phi)*s4
      dir4(3) = pm * c4

      dir3(0) = sqrt(pm*pm + m12)
      dir3(1:3) = -dir4(1:3)
** I calculate E4min, which comes by requiring t24 < t_max(LAB)
*** al and be calculated above      
      ga   = pm*c4

      sol1 = .FALSE.
      sol2 = .FALSE.
      gagabebe = ga*ga-be*be    ! this is always < 0
                                ! al and be always > 0, ga depends on c4
      disc = ga*ga * (al*al+gagabebe*m22)
      if (disc.ge.0.d0) then
c         E_1 = (-al*be+sqrt(disc))/gagabebe
         E_2 = (-al*be-sqrt(disc))/gagabebe         
         E_1 = (al*al*be*be-disc)/gagabebe/gagabebe/E_2
         
         p_1 = sqrt(abs(E_1*E_1 - m22))
         p_2 = sqrt(abs(E_2*E_2 - m22))
         
         if (sign(1.d0,ga*p_1).eq.sign(1.d0,al-be*E_1)
     .        .and.E_1.ge.m2.and.E_1.le.be) then
            sol1 = .TRUE.
         endif
         if (sign(1.d0,ga*p_2).eq.sign(1.d0,al-be*E_2)
     .        .and.E_2.ge.m2.and.E_2.le.be) then
            sol2 = .TRUE.
         endif
      endif

      if (al.le.be*m2) then
         Erange(1,1) = m2
         Erange(1,2) = be
         if (sol1.and..not.sol2) then
            Erange(1,1) = m2
            Erange(1,2) = E_1
         elseif (sol2.and..not.sol1) then
            Erange(1,1) = m2
            Erange(1,2) = E_2
         elseif (sol2.and.sol1) then
            Erange(1,1) = m2
            Erange(1,2) = min(E_1,E_2)
            Erange(2,1) = max(E_1,E_2)
            Erange(2,2) = be
         endif      
      else ! (if al > be*m2), here there should be always only one solution, or none
         Erange(1,1) = m2
         Erange(1,2) = be
         if (sol1.and..not.sol2) then
            Erange(1,1) = E_1
            Erange(1,2) = be
         elseif (sol2.and..not.sol1) then
            Erange(1,1) = E_2
            Erange(1,2) = be
         elseif (sol1.and.sol2) then
            print*,'should not happen! sol[12]',sol1,sol2
         endif
      endif
      
      E4min = Erange(1,1)
      p4min = sqrt(abs(E4min*E4min - m22))
      
      if (nge2) then
         osup = (s - m12 - m22)/sqs*0.5d0
         omink = egminl         
         do k = 2,nph
            call collinearmue_v2(p1,p2,m12,m22,dir3,dir4,ptmp,wcoll,ie)
            w = w * wcoll
            
            k2    = dot(kone,kone)
            qc4q  = tridot(kone,dir4)/pm
            nu    = tridot(kone,ptmp)
            cg4   = tridot(ptmp,dir4)/pm
            omaxk = (s+m22-m12+k2+2.d0*E4min*kone(0)-2.d0*p4min*qc4q
     .              -2.d0*E4min*sqs-2.d0*kone(0)*sqs)
            omaxk = 0.5d0*omaxk/(sqs-E4min-kone(0)+p4min*cg4+nu)

            omaxo = omaxk

            if (cg4.lt.0.d0) then
               ua = s+m22-m12+k2-2.d0*sqs*kone(0)
               ub = -2.d0*(sqs-kone(0))
               uc = -2.d0*pm4*qc4q
               ud = sqs-kone(0)+nu
               ue = cg4

c ie omax =      ooo = (ua+ub*E4min+uc*p4min)/(ud-E4min+ue*p4min)*0.5d0
               oatme = (ua+ub*m2)/(ud-m2)*0.5d0
               if (omaxk.gt.oatme) then
                  orun = 1.001*oatme ! I checked in gnuplot this is safe, it might be not for different me
c                  goto 135
***************************************************************************************
                  E4inf = m2
** I calculate other E4 where omega(E4) = oatme
                  al   = 2.d0*oatme*cg4  + 2.d0*qc4q
                  be   = s+m22-m12+k2-2.d0*kone(0)*sqs
     .                 -2.d0*oatme*(sqs-kone(0)+nu)
                  ga   = 2.d0*(oatme-sqs + kone(0))
                  alalgaga = al*al-ga*ga
                  disc = al*al*(be*be+m22*alalgaga)
                  if (disc.lt.0.d0) then
!     it can happen only for nph > 1
!     Shoud not happen, but if it does orun is the defaut above
                  else
                     E4sup = (be*ga - sqrt(disc))/alalgaga
                     call getomax(E4inf,E4sup,m2,ua,ub,uc,ud,ue,
     .                    E4atmax,orun)
                  endif
***********************************************************************************               
 135              omaxk = orun
               endif
            endif           
c           if (omaxk.lt.egminl) then
c               phsp = 0.d0
c               w    = 0.d0
c               ier  = 2
c               return
c            endif

c here is sampled as 1/om, I need also some 1/(osup-om). I cheched omaxk < p2(0) always
c            call getrnd(csi1,1)
c            anomk = log(omaxk/omink)
c            omk = omink*exp(csi1(1)*anomk)
c            w  = w * anomk*omk
**************************************************************************
c            goto 777
            an1 = log(omaxo/omink)
            an2 = log((osup-omink)/(osup-omaxo))
            pr1 = an1/(an1+an2) ! 0.5d0 !<- less efficient
            pr2 = 1.d0 - pr1

            eps   = (omaxk - omaxo)/omaxo
            fomax = 1d0*(pr1/an1/omaxo + pr2/an2/(osup-omaxo))
            ant   = 1.d0 + fomax*omaxo*eps

            csi1(1) = 0.d0            
            if (omaxk.gt.omaxo) call getrnd(csi1,1)
            if (csi1(1)*ant.le.1.d0) then
               call getrnd(csi1,1)
               if (csi1(1).lt.pr1) then
                  call getrnd(csi1,1)
                  omk = omink*exp(csi1(1)*an1)
               else
                  call getrnd(csi1,1)
                  omk = osup - (osup-omink)*exp(-an2*csi1(1))
               endif
               w = w  / (pr1/an1/omk + pr2/an2/(osup-omk)) * ant
            else
               omk = omaxo + (ant*csi1(1)-1.d0)/fomax
               w = w  / fomax * ant                
            endif

            qph(:,k) = ptmp*omk
            kone = kone + qph(:,k)
         enddo
         
      endif
***   per tagliare uno dei due fotoni sopra soglia e' meglio tagliare solo k=1, ovvero qui,
***  anche se e' tutto simmetrico. Ripensare cmq bene la cosa.
c     c      egminl = 0.5d-1
*     * per l'altro fotone, lo genero > egmin2
      egminl = egmin2
      
      call collinearmue_v2(p1,p2,m12,m22,dir3,dir4,lastk,wcoll,ie)      
      w = w * wcoll

      cg4  = tridot(lastk,dir4)/pm
      qc4q = tridot(kone,dir4)/pm
      nu   = tridot(kone,lastk)
      k2   = dot(kone,kone)

** I calculate the maximum E4, which happens for omega = egminl
      al   = 2.d0*egminl*cg4  + 2.d0*qc4q
      be   = s+m22-m12+k2-2.d0*kone(0)*sqs-2.d0*egminl*(sqs-kone(0)+nu)
      ga   = 2.d0*(egminl-sqs + kone(0))
      
      alalgaga = al*al-ga*ga
c      discold = be*be*ga*ga + (al*al-ga*ga)*(m22*al*al+be*be)
      disc = al*al*(be*be+m22*alalgaga)

      if (disc.lt.0.d0) then
! it can happen only for nph > 1         
         phsp = 0.d0
         w    = 0.d0
         ier  = 3
         return
      else
c         E_1 = (be*ga + sqrt(disc))/alalgaga
         E_2 = (be*ga - sqrt(disc))/alalgaga
         E_1 = (be*be*ga*ga-disc)/alalgaga/alalgaga/E_2
         
         p_1 = sqrt(abs(E_1*E_1 - m22))
         p_2 = sqrt(abs(E_2*E_2 - m22))

         sol1 = .FALSE.
         sol2 = .FALSE.
         if (sign(1.d0,al*p_1).eq.sign(1.d0,be+ga*E_1).and.E_1.ge.m2
     .        .and.E_1.le.p2(0)) sol1 = .TRUE.
         if (sign(1.d0,al*p_2).eq.sign(1.d0,be+ga*E_2).and.E_2.ge.m2
     .        .and.E_2.le.p2(0)) sol2 = .TRUE.

         if (sol1.and..not.sol2) then
c            E4max = E_1
            E4inf = m2
            E4sup = E_1
         elseif(.not.sol1.and.sol2) then
c            E4max = E_2
            E4inf = m2
            E4sup = E_2
         elseif(sol1.and.sol2) then
! it can happen only for nph > 1
c            print*,'two solutions',E_1,E_2,Erange(1,:),Erange(2,:)
c            E4max = max(E_1,E_2) ! this and the commented below seem to give same integrals...

            E4inf = min(E_1,E_2)
            E4sup = max(E_1,E_2)

c            E4max = E_1
c            call getrnd(csi1,1)
c            if (csi1(1).lt.0.5d0) E4max = E_2
c            w = w * 2.d0
         else                   ! sol[12] = .FALSE.
! it can happen only for nph > 1
c     print*,'no solutions'
            phsp = 0.d0
            w    = 0.d0
            ier  = 4
            return    
         endif
      endif
*********
      if (Erange(2,1).gt.-0.5d0) then
         if (     E4sup.lt.Erange(1,1)
     .        .or.E4inf.gt.Erange(2,2)
     .        .or.(E4inf.gt.Erange(1,2).and.E4sup.lt.Erange(2,1))) then
!     no intersection
            phsp = 0.d0
            w    = 0.d0
            ier  = 10
            return 
         elseif (E4inf.ge.Erange(1,1).and.E4sup.le.Erange(1,2)) then
! [E4inf:E4sup] contained in Erange(1,:)
            E4min = max(Erange(1,1),E4inf)
            E4max = min(Erange(1,2),E4sup)
         elseif (E4inf.ge.Erange(2,1).and.E4sup.le.Erange(2,2)) then
! [E4inf:E4sup] contained in Erange(2,:)
            E4min = max(Erange(2,1),E4inf)
            E4max = min(Erange(2,2),E4sup)
         elseif (E4inf.ge.Erange(1,1).and.E4sup.le.Erange(2,1)) then
! E4inf > Erange(1,1) and E4sup between Erange(1,2) and Erange(2,1)
            E4min = max(Erange(1,1),E4inf)
            E4max = min(Erange(2,1),E4sup)
         elseif (E4inf.ge.Erange(1,2).and.E4sup.le.Erange(2,2)) then
! E4inf > Erange(1,2) and E4sup < Erange(2,2)
            E4min = max(Erange(1,2),E4inf)
            E4max = min(Erange(2,2),E4sup)
         else
!     there are two ranges
            Einf1 = max(E4inf,Erange(1,1))
            Esup1 = Erange(1,2)
            Einf2 = Erange(2,1)
            Esup2 = min(E4sup,Erange(2,2))

            r1 = Esup1 - Einf1
            r2 = Esup2 - Einf2

            pr1 = r1/(r1+r2)
            pr2 = 1.d0 - pr1

            call getrnd(csi1,1)
            if (csi1(1).lt.pr1) then
               E4min = Einf1
               E4max = Esup1
               w     = w / pr1
            else
               E4min = Einf2
               E4max = Esup2
               w     = w / pr2
            endif            
         endif                    
c         print*,'NOT READY YET! line ~463 phasespacemue.F'
c         stop
      else  ! when only Erange(1,:) is active
         E4min = max(Erange(1,1),E4inf)
         E4max = min(Erange(1,2),E4sup)
      endif
**** 
      if (E4min.ge.E4max) then  ! ge is essential, because it can be they are the same (eq me)
         phsp = 0.d0
         w    = 0.d0
         ier  = 5
         return
      endif

*** I then generate E4 ~ 1/(s+m22-m12-2*E*sqs) + sth else, see the subroutine
c     al = (s + m22 - m12)*0.5d0/sqs ! for 1 photon
      bmax = sqrt(1.d0-m22/E4max/E4max)
      al = s + m22 - m12 + k2-2.d0*kone(0)*sqs
      al = al*0.5d0/(sqs - kone(0) + bmax*qc4q) ! bmax here is essential to grant (almost) always al > E4max
      
      call get_Ee_new(0,anE,E4,E4min,E4max,al,s,m1,m2,cg4,wE4)
      w = w * wE4

 123  pm4   = sqrt(E4*E4-m22)
      p4(0) = E4
      p4(1:3) = dir4(1:3)/pm*pm4

      den  = sqs-kone(0)-E4+pm4*cg4 + nu      
      
      qph(0,1)   = 0.5d0*
     .     (s+m22-m12+k2-2.d0*sqs*kone(0)
     .     -2.d0*E4*(sqs-kone(0))-2.d0*pm4*qc4q)
     .     /den

c      if (den.lt.0.d0.or.qph(0,1).lt.egminl) then
c seems that den is never < 0. I put abs when used below in phsp
      if (qph(0,1).lt.egminl) then
         ier = 6
         phsp = 0.d0
         w = 0.d0         
         return
      endif

***************************************
      if (nge2) then
c the following if qph(0,1) is generated from egmin2 and qph(0,2) from egmin        
         if (qph(0,2).lt.egmin2) w = w *2.d0 ! statistical weight, if qph(0,1) generated > egmin2!!
cc the commented one only if both photons are generated from egmin         
c      if (qph(0,1).gt.egmin2.or.qph(0,2).gt.egmin2) then
c            continue
c         else
c            ier = 7
c            phsp = 0.d0
c            w = 0.d0         
c            return
c         endif
      endif
***************************************
      
      qph(1:3,1) = lastk(1:3)*qph(0,1)

      kone = kone + qph(:,1)
            
      p3 = p1 + p2 - kone - p4

      if (p3(0).lt.m1) then
         ier = 8
         phsp = 0.d0
         w = 0.d0         
         return
      endif
      
      call tcut(p1,p2,p3,p4,tmin,tmax,icut)
      if (icut.gt.0) then
         ier = 9
         phsp = 0.d0
         w = 0.d0         
         return
      endif
      
      prodomega = 1.d0
      do k = 1,nph
         prodomega = prodomega * 0.5d0*qph(0,k)
         do i = 0,3
            qphext(k,i) = qph(i,k)
         enddo
      enddo

      phsp = phsp * 0.25d0 * twopis(nph) * prodomega * pm4/abs(den)
      phsp = phsp/factorial(nph)

      if (qph(0,2).gt.enmaxg2) then
         enmaxg2 = qph(0,2)
c         print*,enmaxg2
      endif
      
      if (phsp.lt.0.d0) then
         print*,'Should not happen!'
      endif
***   regulator
!     not needed for nph = 1, think about nph > 1
! yes, it's needed for energy sampling, I think it's better if it's symmetric
c$$$      if (nph.ge.2) then
c$$$         prodwe = 1.d0
c$$$         do k = 1,nph
c$$$            prodwe = prodwe/qph(0,k)
c$$$         enddo
c$$$
c$$$         den = 0.d0
c$$$         do k = 1,nph
c$$$            den = den + prodwe * qph(0,k)/(p2(0)-qph(0,k))
c$$$         enddo
c$$$      endif
***   end regulator      

      return
      end      
********************************************************************
      subroutine phasespacemue_v1(s,m1,m2,m3,m4,nph,p1,p2,p3,p4,qphext,
     >     phsp,w,ier)
* this is the one used up to July 2019. It works very well, except that many points are outside bounds
      implicit double precision (a-h,l,m,o-z)
      parameter (pi = 4.d0*atan(1.d0))
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),v3(0:3),q3(0:3)
      dimension csi1(1),csi2(2),csi3(3),csi4(4),csi5(5)
      dimension qph(0:3,max(1,nph)),ptmp(0:3),qphext(40,0:3)
      dimension ptmp2(0:3),ptmp3(0:3)
      double precision kone(0:3)
      common/photoncutoff/egmin,egmin2
      common/s3eps/epss3_ch1,epss3_ch2,sbefore,egminbefore
      data sbefore,egminbefore /1d12,1d12/
      common/tmintmax/tmincut,tmaxcut

      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab
      
      common/tcutcommon/icut

      common/nfotonicommon/nfotoni
      
      dimension twopis(0:40)
      common/pis/twopis,rmax,ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         do k = 1,40
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
         enddo
         rmax = -100d0
         ifirst = 1
      endif

      nfotoni = nph
      
      p3 = 0.d0
      p4 = 0.d0
      qph(:,1) = 0.d0
      qphext(1,:) = 0.d0

      ier = 0

      if (nph.gt.1) then
         print*,'for now ONLY for up-to-1 photon!!'
         print*,'for now ONLY for up-to-1 photon!!'
         print*,'for now ONLY for up-to-1 photon!!'
         stop
      endif
         
      if (s.lt.(m3+m4)*(m3+m4)) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      
      phsp = 1.d0
      w    = 1.d0
      
**      using  two channels, for mu- or e-line radiation (for nph > 0)
      if (nph.gt.0) then
         ich = 2
c         pchmu = 0.001d0 !0.25d0 !0.00d0 !0.25d0
         pchmu = 0.d0 !!1.d-4
         pche  = 1.d0 - pchmu
         call getrnd(csi1,1)
         if (csi1(1).le.pchmu) ich = 1
      else
         ich = 1
      endif
      if (ich.eq.1) then
         m12 = m1*m1
         m22 = m2*m2
         m32 = m3*m3
         m42 = m4*m4
      else
         m12 = m2*m2
         m22 = m1*m1
         m32 = m4*m4
         m42 = m3*m3
      endif
      
      pm = 0.5d0*sqrt(lambda(s,m12,m22)/s)
     
cc      if (ich.eq.2) then
      p1(0) =  sqrt(pm*pm + m12)
      p1(1) =  0.d0
      p1(2) =  0.d0
      p1(3) =  pm
      
      p2(0) =  sqrt(pm*pm + m22)
      p2(1) =  0.d0
      p2(2) =  0.d0
      p2(3) = -pm
cc      endif
      
      if (nph.eq.0) then
** 2 body phasepace
*     per mettere i limiti a t assumo assenza di fotoni
         e3  = 0.5d0*(s + m32 - m42)/sqrt(s)
         
         pm3 = sqrt(abs(e3*e3 - m32))
         pm1 = pm
         
         e1  = p1(0)
         tmax = m12 + m32 - 2.d0*e1*e3 + 2.d0*pm1*pm3
         tmin = m12 + m32 - 2.d0*e1*e3 - 2.d0*pm1*pm3
         
         tmax = min(tmax,tmaxcut)
         tmin = max(tmin,tmincut)
      
c         call get_t(t,tmin,tmax,s,m12,m32,wt)
** warning! must be used with muon and electron mass^2 !!!
         call get_t(t,tmin,tmax,s,m32,m42,wt)
         w = w * wt

         call getrnd(csi1,1)
         phi = 2.d0*pi*csi1(1)
         w   = w * 2.d0*pi
         
         c3 = t - m12 - m32 + 2.d0*e1*e3
         c3 = c3*0.5d0/pm1/pm3

         arg = 1.d0-c3*c3
         if (arg.lt.0.d0) then
            ier  = 1
            phsp = 0.d0
            w    = 0.d0
            return
         endif      
         s3 = sqrt(arg)

         p3(0) = e3
         p3(1) = pm3 * sin(phi)*s3
         p3(2) = pm3 * cos(phi)*s3
         p3(3) = pm3 * c3

         p4 = p1+p2-p3

         call tcut(p1,p2,p3,p4,tmincut,tmaxcut,icut)
         if (icut.gt.0) then
            phsp = 0.d0
            w    = 0.d0
            ier  = 1
            return
         endif


! I checked that sqrt(lambda) below is correct also if m22 is not me^2         
         phsp = 0.25d0*phsp/sqrt(lambda(s,m12,m22))*twopis(nph)
         
         return
      endif

      sqs = sqrt(s)      
      
** prima genero s3
      if (ich.eq.1) then
         s3max = (sqs-m4)*(sqs-m4)
      else
         s3max = (sqs-m3)*(sqs-m3)
      endif         
c      call getlimitss3fromt(s3max,s3min,s,t,m12,m42,e1,pm1,
c     .     s3maxnew,s3minnew)
c      s3max = min(s3max,s3maxnew)
      
      s3min = m32 ! poi questo verra' modificato quando il fotone > omin!!

cc      probabilmente non funziona bene!
      if (abs(s-sbefore).gt.1d-8.or.abs(egmin-egminbefore).gt.1d-8) then
ccc         call getdeltas3(m32,m42,s,egmin,1.d0,epss3)
         call getdeltas3(m3*m3,m4*m4,s,egmin,1.d0,epss3_ch1)
         call getdeltas3(m4*m4,m3*m3,s,egmin,1.d0,epss3_ch2)
c         print*,s/sbefore
         sbefore = s
         egminbefore = egmin         
      endif
      if (ich.eq.1) then
         s3min = s3min + epss3_ch1
      else
         s3min = s3min + epss3_ch2
      endif
         
c      call getrnd(csi1,1)
cc      s3 = (s3max - s3min) * csi1(1) + s3min
cc      w  = (s3max - s3min) * w
c      anom = log((s3max-m32)/(s3min-m32))
c      s3   = m32 + (s3min-m32)*exp(anom*csi1(1))
c      w    = w * anom * (s3-m32)

      call gets3(s3max,s3min,s,m32,m42,s3,ws3)
      w = w*ws3
      
*     per mettere i limiti a t assumo assenza di fotoni
      e3  = 0.5d0*(s + s3 - m42)/sqs
      pm3 = sqrt(abs(e3*e3 - s3))
      pm1 = pm
      e1  = sqrt(pm1*pm1+m12)
      tmax = m12 + s3 - 2.d0*e1*e3 + 2.d0*pm1*pm3
      tmin = m12 + s3 - 2.d0*e1*e3 - 2.d0*pm1*pm3
      
      tmax = min(tmax,tmaxcut)
      tmin = max(tmin,tmincut)

c      tmax = 0.5d0*tmax
      
c      call get_t(t,tmin,tmax,s,m12,s3,wt)
** warning! must be used with muon and electron mass^2 !!!      
** maybe for radiative events is better m22 or m12 = s3, it must be checked which is the most
** efficient      
      call get_t(t,tmin,tmax,s,m12,m22,wt)
      w = w * wt
      
** ora posso calcolare q3. q3m e' |\vec{q3}| e q30 e' l'energia
      q3m = pm3  !0.5d0/sqs*sqrt(lambda(s,s3,m42))      
      q30 = e3                  !sqrt(s3+q3m*q3m)

      cq3 = t - m12 - s3 + 2.d0*e1*q30
      cq3 = cq3*0.5d0/pm1/q3m

      arg = 1.d0-cq3*cq3
      if (abs(cq3).gt.1.d0) then      
         ier  = 1
         phsp = 0.d0
         w    = 0.d0
         goto 666
      endif      
      sq3 = sqrt(arg)
      
      call getrnd(csi1,1)      
      phi = 2.d0 * pi * csi1(1)
      w   = w * 2.d0*pi
      
      q3(0) = q30
      q3(1) = q3m*sin(phi)*sq3
      q3(2) = q3m*cos(phi)*sq3
      q3(3) = q3m*cq3

      
*** ora calcolo il fotone
      cmax =  1.d0
      cmin = -1.d0

c      om = 0.5*(s3-m32)/(q30-q3m*c) deve essere sempre om > egmin --> limite su cmin
c      q30 - q3m*c = 0.5*(s3-m32)/egmin
      cc = q30 - 0.5d0*(s3-m32)/egmin
      cc = cc/q3m !! this and cmin are wrt to q3!!!!
      if (abs(cc).lt.1.d0) then
         cmin = cc
      endif
**      
c      cth = (cmax - cmin) * csi2(2) + cmin
c      w   = w * (cmax - cmin) 
c      sth = sqrt(1.d0 - cth*cth)
** per ora e' kone/om      
c      kone = 0.d0
c      call getrnd(csi2,2)
c      phi = 2.d0 * pi * csi2(1)
c      w   = w * 2.d0*pi
c      kone(0) = 1.d0
c      kone(1) = sth*sin(phi)
c      kone(2) = sth*cos(phi)
c      kone(3) = cth
c      call rot(-1,q3,kone,kone)
****

c      ptmp(0) = p1(0)
c      ptmp(1:3) = sqrt(p1(0)**2-m12)*q3(1:3)/sqrt(tridot(q3,q3))
      
      call collinearmue(cmin,p1,m12,q3,s3,kone,wcoll,ie)
c      call collinearmue(cmin,p1,m12,ptmp,m12,kone,wcoll,ie)
      w = w*wcoll

      if (ie.gt.0) then
         ier  = 1
         phsp = 0.d0
         w    = 0.d0
         goto 666
      endif      
      
*     * occhio qui cth e' rispettto a q3... quando ci sara' da campionare
*     * anche lungo p1 o p2 bisogna pensare come fare      
      
      dq3km1 = 1.d0/dot(q3,kone)
      om = 0.5d0*(s3-m32)*dq3km1
      
      if (om.lt.egmin) then
         ier  = 1
         phsp = 0.d0
         w    = 0.d0
         goto 666
      endif      
      
! .03125 = 1/32
      phsp = 0.03125d0/sqrt(lambda(s,m12,m22))*(s3-m32)*dq3km1*dq3km1
      phsp = phsp * twopis(nph)
      
      kone = om * kone
      
      p3 = q3 - kone
      p4 = p1 + p2 - q3
      
      qph(:,1) = kone
      
 666  if (ich.eq.1) then
         if (ier.gt.0) return
      else !if (ich.eq.2) then
! mirroring for ich = 2 and exchanging
         ptmp = p1
         p1 = p2
         p2 = ptmp
         p1(1:3)    = -p1(1:3)
         p2(1:3)    = -p2(1:3)
         if (ier.gt.0) return
         
         ptmp = p3
         p3 = p4
         p4 = ptmp
         
         p3(1:3)    = -p3(1:3)
         p4(1:3)    = -p4(1:3)
         qph(1:3,1) = -qph(1:3,1)
      endif

*** t cut      
      call tcut(p1,p2,p3,p4,tmincut,tmaxcut,icut)
      if (icut.gt.0) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      
cc  regulator
c$$$      if (ich.eq.1) then
c$$$         anum = (1.d0/dot(p1,qph(:,1))*p1(0)*qph(0,1) + 
c$$$     .                   1.d0/dot(p3,qph(:,1))*p3(0)*qph(0,1))
c$$$         den = pchmu*anum + pche*(1.d0/dot(p2,qph(:,1))*p2(0)*qph(0,1)
c$$$     .                          + 1.d0/dot(p4,qph(:,1))*p4(0)*qph(0,1))
c$$$      else
c$$$         anum = (1.d0/dot(p2,qph(:,1))*p2(0)*qph(0,1) + 
c$$$     .                       1.d0/dot(p4,qph(:,1))*p4(0)*qph(0,1))
c$$$         den = pche*anum + pchmu*(1.d0/dot(p1,qph(:,1))*p1(0)*qph(0,1) 
c$$$     .                          + 1.d0/dot(p3,qph(:,1))*p3(0)*qph(0,1))
c$$$      endif      
c$$$      w = w * anum/den

      eg2 = qph(0,1)*qph(0,1)
      if (ich.eq.1) then
         anum = 1.d0/dot(p1,qph(:,1))/dot(p3,qph(:,1))*eg2
         den  = pchmu*anum + pche/dot(p2,qph(:,1))/dot(p4,qph(:,1))*eg2
      else
         anum = 1.d0/dot(p2,qph(:,1))/dot(p4,qph(:,1))*eg2
         den  = pche*anum + pchmu/dot(p1,qph(:,1))/dot(p3,qph(:,1))*eg2
      endif      
      w = w * anum/den
** end regulator
      
      do i = 0,3
         qphext(1,i) = qph(i,1)
      enddo

c     1/(ab) = 1/(a+b) * (1/a + 1/b) 
c$$$      if (ich.eq.2) then
c$$$         q3 = p4 + kone
c$$$         a = dot(q3,kone)/q3(0)/kone(0)
c$$$         b = dot(p2,kone)/p2(0)/kone(0)
c$$$
c$$$         x = 1.d0/a + 1.d0/b
c$$$         y = 1.d0/a/b
c$$$         z = 1.d0/(a+b)
c$$$         r = y/x
c$$$         if (r.gt.rmax) then
c$$$            rmax = r
c$$$            print*,y/x,x/y,rmax,p4(0),kone(0)
c$$$
c$$$         endif                  
c$$$      endif
      return
      end      
********************************************************************
      subroutine phasespacemue_new(s,m1,m2,m3,m4,nph,p1,p2,p3,p4,qphext,
     >     phsp,w,ier)
*     * in this I follow messpv, the independent variable is not omega but E4. I do not attempt
*     * any sample along the muon line because it's really negligible. This works only for up to nph = 1      
      implicit double precision (a-h,l,m,o-z)
      parameter (pi = 4.d0*atan(1.d0))
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision kvers(0:3)
      dimension csi1(1),csi2(2),csi3(3),csi4(4),csi5(5)
      dimension qph(0:3,max(1,nph)),ptmp(0:3),qphext(40,0:3)
      dimension ptmp2(0:3),ptmp3(0:3)
      double precision kone(0:3),dir4(0:3)
      common/photoncutoff/egmin,egmin2
      common/tmintmax/tmincut,tmaxcut

      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab
      
      common/tcutcommon/icut

      common/nfotonicommon/nfotoni
      
      dimension twopis(0:40)
      common/pis/twopis,rmax,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         do k = 1,40
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
         enddo
         rmax = -100d0
         ifirst = 1
      endif

      nfotoni = nph

      p3 = 0.d0
      p4 = 0.d0
      qph(:,1) = 0.d0
      qphext(1,:) = 0.d0

      ier = 0

      if (nph.gt.1) then
         print*,'for now ONLY for up-to-1 photon!!'
         print*,'for now ONLY for up-to-1 photon!!'
         print*,'for now ONLY for up-to-1 photon!!'
         stop
      endif
         
      if (s.lt.(m3+m4)*(m3+m4)) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      
      phsp = 1.d0
      w    = 1.d0
      
      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      
      pm = 0.5d0*sqrt(lambda(s,m12,m22)/s)
     
      p1(0) =  sqrt(pm*pm + m12)
      p1(1) =  0.d0
      p1(2) =  0.d0
      p1(3) =  pm
      
      p2(0) =  sqrt(pm*pm + m22)
      p2(1) =  0.d0
      p2(2) =  0.d0
      p2(3) = -pm
      
      if (nph.eq.0) then
** 2 body phasepace
*     per mettere i limiti a t assumo assenza di fotoni
         e3  = 0.5d0*(s + m32 - m42)/sqrt(s)
         
         pm3 = sqrt(abs(e3*e3 - m32))
         pm1 = pm
         
         e1  = p1(0)
         tmax = m12 + m32 - 2.d0*e1*e3 + 2.d0*pm1*pm3
         tmin = m12 + m32 - 2.d0*e1*e3 - 2.d0*pm1*pm3
         
         tmax = min(tmax,tmaxcut)
         tmin = max(tmin,tmincut)
      
c         call get_t(t,tmin,tmax,s,m12,m32,wt)
** warning! must be used with muon and electron mass^2 !!!
         call get_t(t,tmin,tmax,s,m32,m42,wt)
         w = w * wt

         call getrnd(csi1,1)
         phi = 2.d0*pi*csi1(1)
         w   = w * 2.d0*pi
         
         c3 = t - m12 - m32 + 2.d0*e1*e3
         c3 = c3*0.5d0/pm1/pm3

         arg = 1.d0-c3*c3
         if (arg.lt.0.d0) then
            ier  = 1
            phsp = 0.d0
            w    = 0.d0
            return
         endif      
         s3 = sqrt(arg)

         p3(0) = e3
         p3(1) = pm3 * sin(phi)*s3
         p3(2) = pm3 * cos(phi)*s3
         p3(3) = pm3 * c3

         p4 = p1+p2-p3

         call tcut(p1,p2,p3,p4,tmincut,tmaxcut,icut)
         if (icut.gt.0) then
            phsp = 0.d0
            w    = 0.d0
            ier  = 1
            return
         endif

! I checked that sqrt(lambda) below is correct also if m22 is not me^2         
         phsp = 0.25d0*phsp/sqrt(lambda(s,m12,m22))*twopis(nph)
         
         return
      endif

      sqs = sqrt(s)

!    these are e4 and pm4 in the case there are no photons... I use these also
!    to find the limits on c4 (in the cm)      
      e4  = 0.5d0*(s + m42 - m32)/sqs
      pm4 = sqrt(e4*e4-m42)

ccc   cmin = (tmincut - m22 - m42 + 2.d0*p2(0)*e4)/pm/pm4*0.5d0
      cmin = -1.d0
      cmax = (tmaxcut - m22 - m42 + 2.d0*p2(0)*e4)/pm/pm4*0.5d0

      cveromax = -cmin
      cveromin = -cmax !cos(pi-acos(cmax))
      
      cmax = cveromax
      cmin = cveromin
      
      call getrnd(csi1,1)
      phi4 = 2.d0*pi*csi1(1)
c      call getrnd(csi1,1)
c      c4   = 2.d0 * csi1(1) - 1.d0

      call getc4(cmin,cmax,c4,wc4)

      s4   = sqrt(1.d0 - c4*c4)
      sphi4 = sin(phi4)
      cphi4 = cos(phi4)

      w = w * 2.d0*pi * wc4

      
      dir4(0) = e4
      dir4(1) = pm4 * s4 * sphi4
      dir4(2) = pm4 * s4 * cphi4
      dir4(3) = pm4 * c4

*** photon angles
      call collinearmue(-1.d0,p2,m22,dir4,m42,kone,wang,ierr)
      if (ierr.gt.0) then
         ier = 1
         phsp = 0.d0
         w = 0.d0
         return
      endif
      w = w * wang
      
*** limits on E4      
      zed = tridot(dir4,kone)/pm4
      E4min = m4

      dm2 = m42 - m32
      a   = 2.d0*(sqs-egmin)
      b   = 2.d0*egmin * sqs - s - dm2
      c   = 2.d0*egmin*zed

      disc = b*b - m42 * (a*a - c*c)
      if (disc.lt.0.d0) then
         ier = 1
         phsp = 0.d0
         w = 0.d0
         return
      endif
      
      E4max_1 = -a*b + abs(c)*sqrt(disc)
      E4max_2 = -a*b - abs(c)*sqrt(disc)
      E4max_1 = E4max_1/(a*a-c*c)
c      E4max_2 = E4max_2/(a*a-c*c)      
      E4max_2 = (b*b + m42*c*c)/(a*a-c*c)/E4max_1 ! more precise (usually)

      x1 = a*E4max_1 + b ! this is the left side of the equation i solved to get E4max_[12]
      y1 = -c*sqrt(E4max_1**2 - m42) ! this is the right side
      x2 = a*E4max_2 + b
      y2 = -c*sqrt(E4max_2**2 - m42)

!     the solution must have the same sign of left and right side
      if (x1/y1.lt.0.d0.and.x2/y2.lt.0.d0) then
!     no solutions
         ier = 1
         phsp = 0.d0
         w = 0.d0
         return         
      elseif (x1/y1.gt.0.d0.and.x2/y2.lt.0.d0) then
         E4max = E4max_1
      elseif (x2/y2.gt.0.d0.and.x1/y1.lt.0.d0) then
         E4max = E4max_2
      elseif (x1/y1.gt.0.d0.and.x2/y2.gt.0.d0) then
!     two solutions
         print*,'two solutions, i think thats strange...'
         call getrnd(csi1,1)
         if (csi1(1).lt.0.5d0) then
            E4max = E4max_1
         else            
            E4max = E4max_2
         endif
         w = w * 2.d0
      endif

      if (E4max.lt.E4min) then
         ier = 1
         phsp = 0.d0
         w = 0.d0
         return         
      endif
      
c      call getrnd(csi1,1)
c      E4 = (E4max - E4min)*csi1(1) + E4min
c      w = w * (E4max - E4min)
      call getE4(E4max,E4min,s,m32,m42,E4,we4)
      w = w * we4
      
      pm4   = sqrt(E4*E4-m42)
      p4(0) = E4
      p4(1) = pm4 * s4 * sphi4
      p4(2) = pm4 * s4 * cphi4
      p4(3) = pm4 * c4

      denom = sqs - E4 + pm4*zed
      omega = 0.5d0*(s-2.d0*sqs*E4 + dm2)/denom

      kone = omega * kone

      p3 = p1 + p2 - p4 - kone 

      phsp = 0.125d0*omega*pm4/denom*twopis(nph)
            
      do i = 0,3
         qphext(1,i) = kone(i)
      enddo

      return
      end
******************************************************************************
      subroutine collinearmue_new_ok(cmin,p,mp2,q,mq2,vers,w,ie) ! new_ok
      implicit double precision (a-h,m,o-z)
      dimension q(0:3),vers(0:3),csi2(2),csi1(1),p(0:3)
      dimension pr(0:3),qr(0:3)
      parameter (pi = 4.d0*atan(1.d0))

      common/forcsampling/ak,almax,almin

      double precision idefint_dc,idefint_dphidc
      external idefint_dc,idefint_dphidc
      
      ie = 0
      
      qm = sqrt(tridot(q,q))
      bq = qm/q(0)      
      pm = sqrt(tridot(p,p))
      bp = pm/p(0)

      wphi = 2.d0*pi

      anq = 1.d0/bq*log((1.d0-bq*cmin)/(1.d0-bq))
      if (cmin.eq.-1.d0) then
         anp = 1.d0/bp*log((1.d0-bp*cmin)/(1.d0-bp))
      else
         call rot(1,q,p,pr)
         vx   = pr(1)/pm
         vy   = pr(2)/pm
         vz   = pr(3)/pm
         anp = idefint_dphidc(bp,vz,cmin,1.d0)/wphi
      endif
      
      rq = 0.5d0         
      rp = 1.d0 - rq
      
      ant = rq * anq + rp * anp
      
      prq = rq * anq/ant
      prp = 1.d0 - prq

      w = 1.d0
      
      call getrnd(csi1,1)      
      if (csi1(1).le.prq) then

         call getrnd(csi2,2)
         
         cth = (1.d0 - (1.d0-bq*cmin)*exp(-bq*anq*csi2(1)))/bq
         sth = sqrt(abs(1.d0-cth*cth))

         phi = 2.d0*pi * csi2(2)! - pi
         
         vers(0) = 1.d0
         vers(1) = sth*sin(phi)
         vers(2) = sth*cos(phi)
         vers(3) = cth         
         call rot(-1,q,vers,vers)
         anum = 1.d0            !prq/anq/(1.d0 - tridot(q,vers)/q(0))         
      else

         call getrnd(csi2,2)
         
         phi = 2.d0*pi * csi2(1)! - pi

         if (cmin.eq.-1.d0) then
            cth = (1.d0 - (1.d0-bp*cmin)*exp(-bp*anp*csi2(2)))/bp
            sth = sqrt(abs(1.d0-cth*cth))
            vers(0) = 1.d0
            vers(1) = sth*sin(phi)
            vers(2) = sth*cos(phi)
            vers(3) = cth            
         else
            Esp = almin*exp(bp*anp*csi2(2))

            cth = (vz+0.5d0*(Esp-ak/Esp))/bp
            sth = sqrt(1.d0-cth*cth)

            vers(0) = 1.d0
            vers(1) = sth*sin(phi)
            vers(2) = sth*cos(phi)
            vers(3) = cth
            
            call rot(-1,q,vers,vers)
         endif
         anum = 1.d0 !prp/anp/(1.d0 - tridot(p,vers)/p(0))
      endif         

      cgq = tridot(q,vers)/qm
      if (cgq.lt.cmin) then
         w  = 0.d0
         ie = 1
         return
      endif

      if (cmin.eq.-1d0) then
         cgp = tridot(p,vers)/pm
         funp = 1.d0/(1.d0 - bp*cgp)
      else
         sgq = sqrt(1.d0-cgq*cgq)
         argsq = bp*bp*(sgq*sgq*(vx*vx+vy*vy)-cgq*cgq*vz*vz)
         argsq = 1.d0 - argsq - 2.d0*bp*cgq*vz         
         funp = 1.d0/sqrt(argsq)
      endif
      den = rq * 1.d0/(1.d0 - bq*cgq)
      den = den + rp * funp
      den = den/ant

      w   = w * wphi * anum/den

      return
      end
**********************************************      
      subroutine collinearmue_WORK(cmin,p,mp2,q,mq2,vers,w,ie) ! ok but slow
! this is to be completed. look for "stop"
      implicit double precision (a-h,m,o-z)
      dimension q(0:3),vers(0:3),csi2(2),csi1(1),p(0:3),s(0:3)
      dimension pr(0:3),qr(0:3)
      double precision kmin(0:3),kmax(0:3),kminr(0:3),kmaxr(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      double precision u
      common/udistrcmncmn/u
      double precision ntrue,napp
      common/normalization/ntrue,napp

      double precision idefint_dc,idefint_dphidc
      external idefint_dc,idefint_dphidc
      
      ie = 0
      
      qm = sqrt(tridot(q,q))
      bq = qm/q(0)      
      pm = sqrt(tridot(p,p))
      bp = pm/p(0)
      
      phimax =  pi
      phimin = -pi
      
      anq = 1.d0/bq*log((1.d0-bq*cmin)/(1.d0-bq))


      call rot(1,q,p,pr)
      pr(1:3) = pr(1:3)/pm
      anp = idefint_dphidc(bp,pr(3),cmin,1.d0)*0.5d0/pi

c      if (isnan(anp)) print*,'NAN!!',anp,bp,cmin,pr(3)/pm
      
cc      print*,anp/anp_false
      
c###      rq = anq/(anq+anp)
      rq = 0.5d0
      rp = 1.d0 - rq
      
      ant = rq * anq + rp * anp
      
      prq = rq * anq/ant
      prp = 1.d0 - prq

      call getrnd(csi1,1)
      
      if (csi1(1).le.prq) then
         call getrnd(csi2,2)
         
         cth = (1.d0 - (1.d0-bq*cmin)*exp(-bq*anq*csi2(1)))/bq
         sth = sqrt(abs(1.d0-cth*cth))

         phi = 2.d0*pi * csi2(2) - pi
         wphi = 2.d0*pi
         
         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth         
         call rot(-1,q,vers,vers)
         anum = 1.d0            !prq/anq/(1.d0 - tridot(q,vers)/q(0))         
      else

         if (cmin.gt.0.1d0.and.cmin.lt.0.5d0) then

            call getrnd(csi1,1)
            phi = 2.d0*pi * csi1(1) - pi
            sph = sin(phi)
            cph = cos(phi)
            
            smin = sqrt(1.d0 - cmin*cmin)
            
            kmin(0) = 1.d0
            kmin(1) = smin * sph
            kmin(2) = smin * cph
            kmin(3) = cmin

            kmax(0) = 1.d0
            kmax(1) = 0.d0
            kmax(2) = 0.d0
            kmax(3) = 1.d0

            
            pt = sqrt(pr(1)*pr(1)+pr(2)*pr(2)) ! pr already normalized

            kminr(1) =  pr(3)*pr(1)/pt*kmin(1)
     .                + pr(3)*pr(2)/pt*kmin(2)
     .                - pt*kmin(3)  
            kminr(2) = - pr(2)/pt*kmin(1)
     .                 + pr(1)/pt*kmin(2)
            kminr(3) = pr(1)*kmin(1)+pr(2)*kmin(2)+pr(3)*kmin(3)

            kmaxr(1) =  pr(3)*pr(1)/pt*kmax(1)
     .                + pr(3)*pr(2)/pt*kmax(2)
     .                - pt*kmax(3)  
            kmaxr(2) = - pr(2)/pt*kmax(1)
     .                 + pr(1)/pt*kmax(2)
            kmaxr(3) = pr(1)*kmax(1)+pr(2)*kmax(2)+pr(3)*kmax(3)

            print*,kmin(1:3)
            print*,kmax(1:3)
            
            call rot(1,pr,kmax,kmax)
            call rot(1,pr,kmin,kmin)
            
            
            print*,kmin(1:3)
            print*,kminr(1:3)
            print*,kmax(1:3)
            print*,kmaxr(1:3)
            print*,pr(1:3)

            print*,'stop at line 869 phasespacemue.F'
            
            
            stop
            
         endif
         
         c0 = -1.d0
         c1 =  1.d0      
         anp_false = 1.d0/bp*log((1.d0-bp*c0)/(1.d0-bp*c1))

         wphi = phimax - phimin
         ccontrol = cmin-1.d0
         jjj = 0
         do while(ccontrol.lt.cmin)
            call getrnd(csi2,2)
            jjj = jjj + 1
            cth = (1.d0 - (1.d0-bp*c0)*exp(-bp*anp_false*csi2(1)))/bp
            sth = sqrt(abs(1.d0-cth*cth))
            
            phi  = (phimax - phimin) * csi2(2) + phimin
            
            vers(0) = 1.d0
            vers(1) = sth*cos(phi)
            vers(2) = sth*sin(phi)
            vers(3) = cth

            ccontrol = tridot(vers,q)/qm

         enddo
         
         if (p(3).lt.0.d0) call rot(-1,p,vers,vers)

         anum = 1.d0 !prp/anp/(1.d0 - tridot(p,vers)/p(0))
      endif         

      cgq = tridot(q,vers)/qm
      if (cgq.lt.cmin) then
         w  = 0.d0
         ie = 1
         return
      endif
      cgp = tridot(p,vers)/pm
      
      den = rq * 1.d0/(1.d0 - bq*cgq)
      den = den + rp * 1.d0/(1.d0 - bp*cgp)
      den = den/ant
      
      w   = wphi * anum/den
      return
      end
**********************************************      
******************************************************************************
      subroutine collinearmue(cmin,p,mp2,q,mq2,vers,w,ie) ! default
      implicit double precision (a-h,m,o-z)
      dimension q(0:3),vers(0:3),csi2(2),csi1(1),p(0:3),s(0:3)
      dimension pr(0:3),qr(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      double precision u
      common/udistrcmncmn/u
      double precision ntrue,napp
      common/normalization/ntrue,napp

      double precision idefint_dc,idefint_dphidc
      external idefint_dc,idefint_dphidc
      
      ie = 0
      
      qm = sqrt(tridot(q,q))
      bq = qm/q(0)      
      pm = sqrt(tridot(p,p))
      bp = pm/p(0)
      
** as 1/(1-b*c)

      call getrnd(csi1,1)      
      call getrnd(csi2,2)
***************
c$$$      qr = q/qm
c$$$      phi0 = atan2(qr(2),qr(1))
c$$$      if (phi0.le.0.d0) phi0 = phi0 + 2.d0*pi
c$$$      Aphi0 = qr(1)*cos(phi0)+qr(2)*sin(phi0)
c$$$c      test = atan2(sin(phi),cos(phi))
c$$$c      if (test.lt.0.d0) test = test + 2.d0*pi
c$$$c      print*,phi/test
c$$$      Aphi = qr(1)*cos(phi)+qr(2)*sin(phi)
c$$$      c0 = cos(atan(Aphi/qr(3)))
***************
      c0 = -1.d0
      c1 =  1.d0
      phimax =  pi
      phimin = -pi
**************

c      print*,atan2(q(2),q(1))*180.d0/pi
      
c### IT SEEMS TO WORK!!      
c      thcmin = acos(cmin)
c      thq    = acos(q(3)/qm)
c      c0 = cos(min(pi,thcmin+thq))
c      c1 = cos(max(0.d0,thq-thcmin))

** for phi doesn't work...      
c      phiq   = atan2(q(2),q(1))
c      phimin = max(-pi,phiq - thcmin)
c      phimax = min(pi,phiq + thcmin)
c      print*,phimin,phimax
      
c#######
c      print*,''
c      print*,thq*180.d0/pi,thcmin*180.d0/pi
c      print*,(thcmin+thq)*180.d0/pi,(-thcmin+thq)*180.d0/pi
c      print*,c0,c1
      
      anq = 1.d0/bq*log((1.d0-bq*cmin)/(1.d0-bq))
      anp = 1.d0/bp*log((1.d0-bp*c0)/(1.d0-bp*c1))

c###      rq = anq/(anq+anp)
      rq = 0.5d0
      rp = 1.d0 - rq
      
      ant = rq * anq + rp * anp
      
      prq = rq * anq/ant
      prp = 1.d0 - prq

      
      if (csi1(1).le.prq) then
         cth = (1.d0 - (1.d0-bq*cmin)*exp(-bq*anq*csi2(1)))/bq
         sth = sqrt(abs(1.d0-cth*cth))

         phi = 2.d0*pi * csi2(2) - pi
         wphi = 2.d0*pi
         
         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth         
         call rot(-1,q,vers,vers)
         anum = 1.d0            !prq/anq/(1.d0 - tridot(q,vers)/q(0))         
      else
         cth = (1.d0 - (1.d0-bp*c0)*exp(-bp*anp*csi2(1)))/bp
         sth = sqrt(abs(1.d0-cth*cth))

         phi  = (phimax - phimin) * csi2(2) + phimin
         wphi = phimax - phimin
                  
         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth

         if (p(3).lt.0.d0) call rot(-1,p,vers,vers)

         anum = 1.d0 !prp/anp/(1.d0 - tridot(p,vers)/p(0))
      endif         

      cgq = tridot(q,vers)/qm
      if (cgq.lt.cmin) then
         w  = 0.d0
         ie = 1
         return
      endif
      cgp = tridot(p,vers)/pm
      
      den = rq * 1.d0/(1.d0 - bq*cgq)
      den = den + rp * 1.d0/(1.d0 - bp*cgp)
      den = den/ant
      
      w   = wphi * anum/den
      return
      end
**********************************************
******************************************************************************
      subroutine collinearmue_v2(p1,p2,m12,m22,p3,p4,vers,w,ie)
      implicit double precision (a-h,m,o-z)
      dimension q(0:3),vers(0:3),csi2(1),csi1(1),p3(0:3),p4(0:3)
      dimension pr(0:3),qr(0:3),pa(0:3),pb(0:3)
      dimension p1(0:3),p2(0:3),q1(0:3),q2(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      double precision u
      common/udistrcmncmn/u
      double precision ntrue,napp
      common/normalization/ntrue,napp

      double precision idefint_dc,idefint_dphidc
      external idefint_dc,idefint_dphidc

**** charges to switch on radiation on legs (0 or 1), from invariants.h
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
********************************************************************      
      ie = 0
** as 1/(1-b*c)

      call getrnd(csi2,1)
      phi  = 2.d0*pi * csi2(1) - pi
      wphi = 2.d0*pi

      pa = p2
      pb = p4

      if (QRmu.eq.1.and.QRe.eq.0) then
! flat ! ok
c         call getrnd(csi1,1)
c         cth = 2.d0*csi1(1)-1.d0
c         sth = sqrt(1.d0-cth*cth)         
c         vers(0) = 1.d0
c         vers(1) = sth*cos(phi)
c         vers(2) = sth*sin(phi)
c         vers(3) = cth
c         w = wphi*2.d0
c         return         
         pa = p1
         pb = p3
      endif
      
      pem = sqrt(tridot(pa,pa))
      be  = pem/pa(0)
      
      ane = 1.d0/be*log((1.d0+be)/(1.d0-be))
      ant = 2.d0*ane
      
      pe4 = ane / ant
      pe2 = ane / ant
      
      call getrnd(csi1,1)     
      if (csi1(1).le.pe4) then
         call getrnd(csi2,1)
         cth = (1.d0 - (1.d0+be)*exp(-be*ane*csi2(1)))/be
         sth = sqrt(abs(1.d0-cth*cth))
         
         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth         
         call rot(-1,pb,vers,vers)
      else
         call getrnd(csi2,1)
         cth = (1.d0 - (1.d0+be)*exp(-be*ane*csi2(1)))/be
         sth = sqrt(abs(1.d0-cth*cth))

         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth
         call rot(-1,pa,vers,vers)         
      endif         

      cg4 = tridot(pb,vers)/pem
      cg2 = tridot(pa,vers)/pem
      
      w   = wphi * 1.d0/(
     .       pe4/ane/(1.d0-be*cg4)
     .     + pe2/ane/(1.d0-be*cg2)
     .     )
      return
      end
**********************************************
******************************************************************************
      subroutine collinearmueCRF(p2,p4,m12,m22,cmin,cmax,vers,w,ie)
! I exploit here that p2 is along -z and p4 along z
      implicit double precision (a-h,m,o-z)
      dimension q(0:3),vers(0:3),csi2(1),csi1(1),p3(0:3),p4(0:3)
      dimension pr(0:3),qr(0:3),pa(0:3),pb(0:3)
      dimension p1(0:3),p2(0:3),q1(0:3),q2(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      double precision u
      common/udistrcmncmn/u
      double precision ntrue,napp
      common/normalization/ntrue,napp

      double precision idefint_dc,idefint_dphidc
      external idefint_dc,idefint_dphidc

**** charges to switch on radiation on legs (0 or 1), from invariants.h
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
********************************************************************      
      ie = 0
** as 1/(1-b*c)

      call getrnd(csi2,1)
      phi  = 2.d0*pi * csi2(1) - pi
      wphi = 2.d0*pi

      pa = p2
      pb = p4

      pem = sqrt(tridot(pa,pa))
      be  = pem/pa(0)
      
      an4 = 1.d0/be*log((1.d0-be*cmin)/(1.d0-be*cmax))

      cmin2 = -cmax
      cmax2 = -cmin
           
      an2 = 1.d0/be*log((1.d0-be*cmin2)/(1.d0-be*cmax2))
      
      ant = an2 + an4
      
      pe4 = an4 / ant
      pe2 = an2 / ant
      
      call getrnd(csi1,1)     
      if (csi1(1).le.pe4) then
         call getrnd(csi2,1)
         cth = (1.d0 - (1.d0-be*cmin)*exp(-be*an4*csi2(1)))/be
         sth = sqrt(abs(1.d0-cth*cth))
         
         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth         
         call rot(-1,pb,vers,vers)
      else
         call getrnd(csi2,1)
         cth = (1.d0 - (1.d0-be*cmin2)*exp(-be*an2*csi2(1)))/be
         sth = sqrt(abs(1.d0-cth*cth))

         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth
         call rot(-1,pa,vers,vers)         
      endif         

      cg4 = tridot(pb,vers)/pem
      cg2 = tridot(pa,vers)/pem
      
      w   = wphi * 1.d0/(
     .       pe4/an4/(1.d0-be*cg4)
     .     + pe2/an2/(1.d0-be*cg2)
     .     )
      return
      end
**********************************************
      
******************************************************************************
      subroutine collinearpairs(p1,p2,m12,m22,p3,p4,ek,mk2,bk,vers,w,ie)
      implicit double precision (a-h,m,o-z)
      dimension q(0:3),vers(0:3),csi2(1),csi1(1),p3(0:3),p4(0:3)
      dimension pr(0:3),qr(0:3),pa(0:3),pb(0:3)
      dimension p1(0:3),p2(0:3),q1(0:3),q2(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      double precision u
      common/udistrcmncmnpairs/u
      double precision ntrue,napp
      common/normalizationpairs/ntrue,napp

      double precision idefint_dc,idefint_dphidc
      external idefint_dc,idefint_dphidc

**** charges to switch on radiation on legs (0 or 1), from invariants.h
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
********************************************************************      
      ie = 0
** as 1/(1-b*c)

      call getrnd(csi2,1)
      phi  = 2.d0*pi * csi2(1) - pi
      wphi = 2.d0*pi

      if (bk.eq.0.d0) then
         call getrnd(csi2,1)
         c = 2.d0*csi2(1)-1.d0
         s = sqrt(abs(1.d0 -c*c))
         vers(0) = 1.d0
         vers(1) = sin(phi)*s
         vers(2) = cos(phi)*s
         vers(3) = c

         w = 2.d0*wphi
         
         return         
      endif         
      
      pa = p2
      pb = p4

      if (QRmu.eq.1.and.QRe.eq.0) then
! flat ! ok
c         call getrnd(csi1,1)
c         cth = 2.d0*csi1(1)-1.d0
c         sth = sqrt(1.d0-cth*cth)         
c         vers(0) = 1.d0
c         vers(1) = sth*cos(phi)
c         vers(2) = sth*sin(phi)
c         vers(3) = cth
c         w = wphi*2.d0
c         return         
         pa = p1
         pb = p3
      endif
      
      pem = sqrt(tridot(pa,pa))
      be  = pem/pa(0)

      A  = 0.5d0*mk2/pa(0)/ek
      rp = (1.d0 + A)/be/bk
      rm = (1.d0 - A)/be/bk
      
      anp = log((rp+1.d0)/(rp-1.d0))
      anm = log((rm+1.d0)/(rm-1.d0))
      ant = anp + anm
      
      pe4 = anp / ant
      pe2 = anm / ant
      
      call getrnd(csi1,1)     
      if (csi1(1).le.pe4) then
         call getrnd(csi2,1)

         cth = rp - (rp+1.d0)*exp(-anp*csi2(1))         
         sth = sqrt(abs(1.d0-cth*cth))
         
         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth         
         call rot(-1,pb,vers,vers)
      else
         call getrnd(csi2,1)
         cth = rm - (rm+1.d0)*exp(-anm*csi2(1))         
         sth = sqrt(abs(1.d0-cth*cth))

         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth
         call rot(-1,pa,vers,vers)         
      endif         

      cg4 = tridot(pb,vers)/pem
      cg2 = tridot(pa,vers)/pem
      
      w   = wphi * 1.d0/(
     .       pe4/anp/(rp-cg4)
     .     + pe2/anm/(rm-cg2)
     .     )
      return
      end
**********************************************
******************************************************************************
      subroutine collinearmue_wint_v1(p1,p2,m12,m22,p3,p4,vers,w,ie)
! it is not efficient, despite sampling exactly as 1/p2.k/p4.k...
      implicit double precision (a-h,m,o-z)
      dimension vers(0:3),csi2(1),csi1(1)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p1h(0:3),p2h(0:3),p3h(0:3),p4h(0:3)      
      dimension p24h(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      double precision u
      
      ie = 0

      me2h = m22/p2(0)/p2(0)
      
      p2h = p2/p2(0)
      p4h = p4/p4(0)
      
      be = sqrt(1.d0-me2h)

      p2p4h = dot(p2h,p4h)

      p24h  = p2h + p4h

      betone24 = sqrt(tridot(p24h,p24h))/p24h(0)

      sqp = p2p4h + me2h
      sqp = sqrt(sqp)
      sqm = p2p4h - me2h
      sqm = sqrt(sqm)      
      antrue = 2.d0*2.d0*p2p4h/sqp/sqm*log((sqp+sqm)/(sqp-sqm)) ! il 2 e' per 24 e 42

!      antrue = antrue - 2.d0 - 2.d0
      
** as 1/(1-b*c)
            
      ane = 1.d0/be*log((1.d0+be)/(1.d0-be))

      ant = 2.d0*ane
      
      pe4 = ane / ant
      pe2 = ane / ant

      istop = 0
      ic = 0
      rmax = 1.d0/(1.d0-betone24)/p24h(0) ! p24h(0) is always 2.d0
      wphi = 2.d0*pi
      do while(istop.eq.0)

         ic = ic + 1
         
         call getrnd(csi2,1)
         cth = (1.d0 - (1.d0+be)*exp(-be*ane*csi2(1)))/be
         sth = sqrt(abs(1.d0-cth*cth))
         call getrnd(csi2,1)
         phi  = 2.d0*pi * csi2(1) - pi
         
         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth
         
         call getrnd(csi1,1)
         if (csi1(1).le.pe4) then
            call rot(-1,p4h,vers,vers)
         else
            call rot(-1,p2h,vers,vers)
         endif

         p2hv = dot(p2h,vers)
         p4hv = dot(p4h,vers)
         
         ftrue = 2.d0*p2p4h/p2hv/p4hv
         
!         ftrue = ftrue - me2h/p2hv**2 - me2h/p4hv**2

         ffalse = 2.d0*p2p4h*(1.d0/p2hv + 1.d0/p4hv)
         
         call getrnd(csi1,1)
         if (csi1(1)*rmax.lt.ftrue/ffalse) then
            istop = 1
         endif
      enddo

c      w = ftrue ! i.e. 2d0*p2p4h/p2hv/p4hv
      
      w = wphi*antrue/ftrue
      
      return
      
      end
**********************************************
******************************************************************************
      subroutine collinearmue_wint(p1,p2,m12,m22,p3,p4,vers,w,ie)
* adding also mu-like radiation, it's more efficient but quite much slower due to internal unweightening
      implicit double precision (a-h,m,o-z)
      dimension vers(0:3),csi2(1),csi1(1)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      dimension p1h(0:3),p2h(0:3),p3h(0:3),p4h(0:3)      
      dimension p24h(0:3),p13h(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      double precision u
      common/rmaxcmn/rmaxmax
      data rmaxmax /-1.d0/
      
      ie = 0

      me2h = m22/p2(0)/p2(0)
      mm2h = m12/p1(0)/p1(0)
      
      p2h = p2/p2(0)
      p4h = p4/p4(0)
      p1h = p1/p1(0)
      p3h = p3/p3(0)
      
      be = sqrt(1.d0-me2h)
      bm = sqrt(1.d0-mm2h)

      p2p4h = dot(p2h,p4h)
      p1p3h = dot(p1h,p3h)

      p24h  = p2h + p4h
      p13h  = p1h + p3h

      betone24 = sqrt(tridot(p24h,p24h))/p24h(0)
      betone13 = sqrt(tridot(p13h,p13h))/p13h(0)

      sqp = p2p4h + me2h
      sqp = sqrt(sqp)
      sqm = p2p4h - me2h
      sqm = sqrt(sqm)      
      antrue = 2.d0*2.d0*p2p4h/sqp/sqm*log((sqp+sqm)/(sqp-sqm)) ! il 2 e' per 24 e 42

      sqp = p1p3h + mm2h
      sqp = sqrt(sqp)
      sqm = p1p3h - mm2h
      sqm = sqrt(sqm)      
      antrue = antrue + 
     .     2.d0*2.d0*p1p3h/sqp/sqm*log((sqp+sqm)/(sqp-sqm))
      
      antrue = antrue - 2.d0 - 2.d0 - 2.d0 - 2.d0

c      antrue = colllogwint_mue(p1,p2,p3,p4)
      
** as 1/(1-b*c)
            
      ane = 1.d0/be*log((1.d0+be)/(1.d0-be))
      anm = 1.d0/bm*log((1.d0+bm)/(1.d0-bm))

      ant = 2.d0*ane*p2p4h
      ant = ant + 2.d0*anm*p1p3h
      
      pe4 = ane*p2p4h / ant
      pe2 = ane*p2p4h / ant
      pm3 = anm*p1p3h / ant
      pm1 = anm*p1p3h / ant

      istop = 0
      ic = 0
      rmax = 1.d0/(1.d0-betone24)/p24h(0) ! p24h(0) is always 2.d0
      rmax = rmax + 1.d0/(1.d0-betone13)/p13h(0)
      rmax = rmax * 0.95d0
      wphi = 2.d0*pi
      do while(istop.eq.0)
***
cccccc         istop = 1
***
         ic = ic + 1

         call getrnd(csi1,1)
         iel = 0
         if (csi1(1).lt.(pe4+pe2)) iel = 1
         
         call getrnd(csi2,1)
         if (iel.eq.1)  then
            cth = (1.d0 - (1.d0+be)*exp(-be*ane*csi2(1)))/be
         else
            cth = (1.d0 - (1.d0+bm)*exp(-bm*anm*csi2(1)))/bm
         endif
         sth = sqrt(abs(1.d0-cth*cth))
         call getrnd(csi2,1)
         phi  = 2.d0*pi * csi2(1) - pi
         
         vers(0) = 1.d0
         vers(1) = sth*cos(phi)
         vers(2) = sth*sin(phi)
         vers(3) = cth
         
         call getrnd(csi1,1)
         if (csi1(1).le.0.5d0) then
            if (iel.eq.1) then
               call rot(-1,p4h,vers,vers)
            else
               call rot(-1,p3h,vers,vers)
            endif
         else
            if (iel.eq.1) then
               call rot(-1,p2h,vers,vers)
c            else
               ! not needded, p1h is along z call rot(-1,p1h,vers,vers)
            endif
         endif

         p2hv = dot(p2h,vers)
         p4hv = dot(p4h,vers)
         p1hv = dot(p1h,vers)
         p3hv = dot(p3h,vers)
         
         ftrue = 2.d0*p2p4h/p2hv/p4hv
         ftrue = ftrue + 2.d0*p1p3h/p1hv/p3hv
c         ftrue = ftrue + 2.d0*dot(p1h,p2h)/p1hv/p2hv
c         ftrue = ftrue - 2.d0*dot(p1h,p4h)/p1hv/p4hv
c         ftrue = ftrue - 2.d0*dot(p2h,p3h)/p2hv/p3hv
c         ftrue = ftrue + 2.d0*dot(p3h,p4h)/p3hv/p4hv

         
         ftrue = ftrue - me2h/p2hv**2 - me2h/p4hv**2
         ftrue = ftrue - mm2h/p1hv**2 - mm2h/p3hv**2
         
         ffalse = 2.d0*p2p4h*(1.d0/p2hv + 1.d0/p4hv)
         ffalse = ffalse + 2.d0*p1p3h*(1.d0/p1hv + 1.d0/p3hv)

         r = ftrue/ffalse / rmax

         if (r.gt.rmaxmax) then
            rmaxmax = r
            print*,rmaxmax
         endif
         
         if (r.gt.1.d0) print*,'OOOOOPS',r,rmax
         
         call getrnd(csi1,1)
         if (csi1(1).lt.r) then
            istop = 1
         endif
      enddo

c      w = ftrue ! i.e. 2d0*p2p4h/p2hv/p4hv
      
      w = wphi*antrue/ftrue

      
      return
      end
**********************************************
***********************************************************************************************
      function colllogwint_mue(p1,p2,p3,p4)
! this should work also for massless momenta...
      implicit double precision (a-h,l,m,o-z)
      parameter (npart = 4)
      dimension p(0:3),q(0:3),pmat(npart,0:3)
      double precision mass(5)
      dimension masses(npart),eta(npart)
      common/collogmue/ao2pi,masses,eta,phmass,alog2,ifirst
      data ifirst/0/
      common/parameters/ame,ammu,convfac,alpha,pi

#include "invariants.h"      
      
! ch are the charges of the field (not anti-field!!)
! Use this convention: the integer factor in front of the charge
! must be
! -1 --> for incoming particle
! -1 --> for outgoing anti-particle
! +1 --> for outgoing particle
! +1 --> for incoming anti-particle

      if (ifirst.eq.0) then

         eta(1) = -1.d0 * Qmu
         eta(2) = -1.d0 * (-1.d0)
         eta(3) =  1.d0 * Qmu
         eta(4) =  1.d0 * (-1.d0)

         masses(1) = mm*mm
         masses(2) = me*me
         masses(3) = mm*mm
         masses(4) = me*me
         
         ao2pi = 1.d0
         
         ifirst = 1
      endif
      
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

c to check only vertexes ! ok!            
c            if (i.eq.1.and.j.eq.2) etaij = 0.d0
c            if (i.eq.1.and.j.eq.4) etaij = 0.d0
c            if (i.eq.2.and.j.eq.3) etaij = 0.d0
c            if (i.eq.3.and.j.eq.4) etaij = 0.d0
c to check only boxes
c            if (i.eq.1.and.j.eq.3) etaij = 0.d0
c            if (i.eq.2.and.j.eq.4) etaij = 0.d0
c only direct box !  ok!
c            if (i.eq.1.and.j.eq.4) etaij = 0.d0
c            if (i.eq.2.and.j.eq.3) etaij = 0.d0
c only crossed box ! ok !
c            if (i.eq.1.and.j.eq.2) etaij = 0.d0
c            if (i.eq.3.and.j.eq.4) etaij = 0.d0
            
            if (abs(etaij).gt.1.d-3) then
               call rescale_momenta_mue(npart,i,j,pmat,masses,p,q,rho)
               q2 = masses(j)
               vl = 0.5d0*(rho*rho*masses(i) - q2)
               v  = vl/(p(0) - q(0)) 
               arglog = 1.d0+2.d0*vl/q2

               terminfra = 0.d0              
               if (arglog.gt.0.d0) then 
                  tot = log(arglog)
               else
                  tot = 0.d0
               endif
               
               tot = -tot * 2.d0*dot(p,q)/vl* etaij

               softint  = softint + tot
            endif
         enddo
      enddo      
      soft_integral_mue = softint

! diagonal contributions
cc      softint = 0.d0

c      only boxes
c      goto 666
 666  soft_integral_mue = (soft_integral_mue - 8.d0)
      colllogwint_mue = soft_integral_mue
      
      return
      end
***********************************************************************************************
*******************************************************************************
      function idefint_dphidc(b,vz,cmin,cmax)
*** IT WORKS!!!      
*this should be the integral between  -pi and pi, and between cmin and cman in dc
      implicit double precision (a-h,m,o-z)
      double precision idefint_dphidc
      parameter (pi = 4.d0*atan(1.d0))
      parameter (eps2lim = 1d-8)

      common/forcsampling/ak,almax,almin
      
      if (vz.eq.1.d0) then
         ak = 0.d0
         almax = 1.d0-b*cmin
         almin = 1.d0-b*cmax
         idefint_dphidc = log(almax/almin)
         idefint_dphidc = idefint_dphidc/b*2.d0*pi
         return
      endif

      discrmax = b*cmax-vz
      discrmin = b*cmin-vz
      ak       = (1.d0-vz*vz)*(1.d0-b*b)
      epsmax   = ak/discrmax/discrmax
      epsmin   = ak/discrmin/discrmin
      
      sqmax = abs(discrmax)*sqrt(1.d0+epsmax)
      sqmin = abs(discrmin)*sqrt(1.d0+epsmin)

      almax = sqmax+discrmax
      almin = sqmin+discrmin
      
!! I need an expansion when eps's are small and discr's are negative.... 
      eps2max = epsmax*epsmax
      if (discrmax.lt.0.d0.and.eps2max.lt.eps2lim) then
! almax = abs(discrmax)(sqrt(1.d0+epsmax) - 1.d0)

         almax = 0.5d0     * epsmax
     .        -0.125d0     * eps2max
     .        +0.0625d0    * epsmax*eps2max
     .        -0.0390625d0 * eps2max*eps2max
     .        +7d0/256d0   * eps2max*eps2max*epsmax

         almax = almax * abs(discrmax)
      endif
      eps2min = epsmin*epsmin
      if (discrmin.lt.0.d0.and.eps2min.lt.eps2lim) then
! almin = abs(discrmin)(sqrt(1.d0+epsmin) - 1.d0)

         almin = 0.5d0     * epsmin
     .        -0.125d0     * eps2min
     .        +0.0625d0    * epsmin*eps2min
     .        -0.0390625d0 * eps2min*eps2min
     .        +7d0/256d0   * eps2min*eps2min*epsmin

         almin = almin * abs(discrmin)
      endif
      
      idefint_dphidc = log(almax/almin)
      idefint_dphidc = idefint_dphidc/b*2.d0*pi

c      if (isnan(idefint_dphidc)) then
c         print*,''
c         print*,b,vz,cmin,cmax
c         print*,sqmax,sqmin
c         print*,sqrt((b*cmax-vz)**2+(1.d0-vz*vz)*(1.d0-b*b))/sqmax
c         print*,sqrt((b*cmin-vz)**2+(1.d0-vz*vz)*(1.d0-b*b))/sqmin
c         print*,almax
c         print*,sqrt((b*cmax-vz)**2+(1.d0-vz*vz)*(1.d0-b*b))+b*cmax-vz
c         print*,almin
c         print*,sqrt((b*cmin-vz)**2+(1.d0-vz*vz)*(1.d0-b*b))+b*cmin-vz
c      endif
      return
      end
*******************************************************************************
      function thfun(x)
      double precision thfun,x
      thfun = 0.d0
      if (x.gt.0.d0) thfun = 1.d0
      return
      end
******************************************************************************
      subroutine collinearmue_BCK(cmin,p,mp2,q,mq2,vers,w,ie)
      implicit double precision (a-h,m,o-z)
      dimension q(0:3),vers(0:3),csi2(2),csi1(1),p(0:3),s(0:3)
      dimension pr(0:3),qr(0:3)
      parameter (pi = 4.d0*atan(1.d0))
      double precision u
      common/udistrcmncmn/u
      double precision ntrue,napp
      common/normalization/ntrue,napp
      
***   for the moment I assume cmin = -1.d0, I think this reduces a lot usable points
***   but I don't know how to set limits when vers is sampled along p (q?)

      ie = 0
      
      qm = sqrt(tridot(q,q))
      bq = qm/q(0)      
      pm = sqrt(tridot(p,p))
      bp = pm/p(0)

      
c      goto 666
cc      goto 555      
      
** as 1/(1-b*c)
      anq = 1.d0/bq*log((1.d0+bq)/(1.d0-bq))
      anp = 1.d0/bp*log((1.d0+bp)/(1.d0-bp))
      
      prq = anq/(anq+anp)
      prp = 1.d0 - prq
      
      call getrnd(csi1,1)      
      call getrnd(csi2,2)
      phi = 2.d0*pi * csi2(2)
      if (csi1(1).le.prq) then
         cth = (1.d0 - (1.d0+bq)*exp(-bq*anq*csi2(1)))/bq
         sth = sqrt(abs(1.d0-cth*cth))
         vers(0) = 1.d0
         vers(1) = sth*sin(phi)
         vers(2) = sth*cos(phi)
         vers(3) = cth         
         call rot(-1,q,vers,vers)
         anum = 1.d0!prq/anq/(1.d0 - tridot(q,vers)/q(0))
      else
         cth = (1.d0 - (1.d0+bp)*exp(-bp*anp*csi2(1)))/bp
         sth = sqrt(abs(1.d0-cth*cth))
         vers(0) = 1.d0
         vers(1) = sth*sin(phi)
         vers(2) = sth*cos(phi)
         vers(3) = cth         
         anum = 1.d0 !prp/anp/(1.d0 - tridot(p,vers)/p(0))
      endif         

      den = prq/anq/(1.d0 - tridot(q,vers)/q(0))
      den = den + prp/anp/(1.d0 - tridot(p,vers)/p(0))

      w   = 2.d0*pi * anum/den
      
      return

 666  continue

      s = p + q
      sm = sqrt(tridot(s,s))
      bs  = sm/s(0)

      anq = 1.d0/bq*log((1.d0+bq)/(1.d0-bq))/q(0)
      anp = 1.d0/bp*log((1.d0+bp)/(1.d0-bp))/p(0) 
      ans = 1.d0/bs*log((1.d0+bs)/(1.d0-bs))/s(0)

      ansumm1 = 1.d0/(anq+anp+2.d0*ans)
      prq = anq*ansumm1
      prp = anp*ansumm1
      prs = 1.d0 - prq - prp
      
      call getrnd(csi1,1)
      call getrnd(csi2,2)
      phi = 2.d0*pi * csi2(2)
      if (csi1(1).le.prq) then
         cth = (1.d0 - (1.d0+bq)*exp(-bq*anq*q(0)*csi2(1)))/bq
         sth = sqrt(abs(1.d0-cth*cth))
         vers(0) = 1.d0
         vers(1) = sth*sin(phi)
         vers(2) = sth*cos(phi)
         vers(3) = cth         
         call rot(-1,q,vers,vers)
         anum = 1.d0
      elseif(csi1(1).gt.prq.and.csi1(1).le.(prq+prp)) then
         cth = (1.d0 - (1.d0+bp)*exp(-bp*anp*p(0)*csi2(1)))/bp
         sth = sqrt(abs(1.d0-cth*cth))
         vers(0) = 1.d0
         vers(1) = sth*sin(phi)
         vers(2) = sth*cos(phi)
         vers(3) = cth         
         anum = 1.d0
      else
         cth = (1.d0 - (1.d0+bs)*exp(-bs*ans*s(0)*csi2(1)))/bs
         sth = sqrt(abs(1.d0-cth*cth))
         vers(0) = 1.d0
         vers(1) = sth*sin(phi)
         vers(2) = sth*cos(phi)
         vers(3) = cth         
         call rot(-1,s,vers,vers)
         anum = 1.d0
      endif

c      den = prq/anq/(1.d0 - tridot(q,vers)/q(0))
c      den = den + prp/anp/(1.d0 - tridot(p,vers)/p(0))

      den = prq/anq/dot(q,vers)
      den = den + prp/anp/dot(p,vers)
      den = den + prs/ans/dot(s,vers)
      
      w   = 2.d0*pi * anum/den
      return

 555  continue  ! uses feynman parameterization!!
                ! it's very cool
                ! it's a sort of multichannel in integral form!!
                ! it's less efficient because there is an extra integration, i guess
      w = 1.d0

c      qr = q/q(0)
c      pr = p/p(0)      
            
      
c      call getrnd(csi1,1)
c      u = 2.d0 * csi1(1) - 1.d0
c      w = 1.d0 ! viene da 1/2 * 2

      call getufeyn(p+q,p-q,u,wu)
      w = w * wu

      s = p+q + u*(p-q)
            
      sm = sqrt(tridot(s,s))
      bs  = sm/s(0)

c      print*,p(0),q(0),sqrt(dot(q,q))
c      stop
      
!      ans = 1.d0/bs*(1.d0/(1.d0 - bs) - 1.d0/(1.d0+bs)) ! i.e.
      ans = 2.d0/(1.d0 - bs*bs)

      call getrnd(csi2,2)
      phi = 2.d0*pi * csi2(2)

      cth = 1.d0 - 1.d0/(bs*ans*csi2(1) + 1.d0/(1.d0+bs))
      cth = cth/bs
      
      sth = sqrt(abs(1.d0-cth*cth))
      vers(0) = 1.d0
      vers(1) = sth*sin(phi)
      vers(2) = sth*cos(phi)
      vers(3) = cth         
      call rot(-1,s,vers,vers)
      
      w = w * 2.d0*pi * ans * dot(s,vers)**2/s(0)/s(0)*0.5d0
c      w = w * 2.d0*pi * ans * dot(p,vers)*dot(q,vers)/p(0)/q(0)*0.5d0

      
      end
*******************************************************************************
      subroutine getdeltas3(m32,m42,s,egmin,c,eps)      
      implicit double precision (a-h,m,o-z)
      eps  = 0.d0
      if (egmin.eq.0.d0) return
      epsb = 1000.d0
c      ic = 0
c      print*,''
      do while(abs(epsb-eps).gt.1d-10)
c         ic = ic + 1
         epsb = eps
         ee  = 0.5d0/sqrt(s)*(s+m32-m42 + eps)
         b   = sqrt(1.d0-m32/ee/ee)
         eps = 2.d0*egmin*ee*(1.d0-b*c)
c         print*,eps
      enddo
      
      return
      end
*******************************************************************************
      subroutine gets3max(m32,m42,s,egmin,c,s3max,s3maxnew)
      implicit double precision (a-h,m,o-z)
      s3maxnew  = s3max
      if (egmin.eq.0.d0) return
      s3maxb = 1000.d10
      print*,'no!'
      do while(abs(s3maxb-s3maxnew).gt.1d-10)
         s3maxb = s3maxnew
         ee  = 0.5d0/sqrt(s)*(s+s3maxnew-m42)
         b   = sqrt(1.d0-s3maxnew/ee/ee)
         s3maxnew = m32 + 2.d0 * egmin * (ee+b*c)
      enddo
      
      return
      end
*******************************************************************************
      subroutine tcut(p1,p2,p3,p4,tmin,tmax,icut)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      icut = 0
!!!!!!!!!!!
      return
!!!!!!!!!!!
      t13 = dot(p1-p3,p1-p3)
      t24 = dot(p2-p4,p2-p4)
      
      if (t13.lt.tmin.or.t13.gt.tmax) icut = 1
      if (t24.lt.tmin.or.t24.gt.tmax) icut = 1

      return
      end

******
*******************************************************************************
      subroutine getufeyn(s,d,u,w)
! samples u so that 1/(1-b^2) is flat
! integrals from https://en.wikipedia.org/wiki/List_of_integrals_of_rational_functions
      implicit double precision (a-i,m,n,o-z)
      dimension s(0:3),d(0:3)
      dimension csi1(1),csi2(2)
      common/normalization/ntrue,napp

      w = 1.d0
      
      a = dot(d,d)
      b = 2.d0*dot(s,d)
      c = dot(s,s)
      
      det = 4.d0*a*c - b*b
      
      sdet = sqrt(-det)

      imin = intind(-1.d0,sdet,a,b)
      norm = intind( 1.d0,sdet,a,b) - imin
      

      call getrnd(csi1,1)
      u = (tanh(-(norm * csi1(1) + imin)*0.5d0*sdet)*sdet-b)*0.5d0/a
      
c         ee = exp(sdet*(napp*csi2(1)+imin))
c         y = sdet*(1.d0-ee)/(1.d0+ee)
c         u = (y - b)*0.5d0/a
c         print*,u/ub
         
         
cc      w = napp*(a*u*u + b*u + c) * 0.5d0

      w = norm*(a*u*u + b*u + c)

      
      return
      end

*********************************************
      function intind(x,sdet,a,b)
      implicit double precision (a-i,m,n,o-z)

c      arg = (2.d0*a*x + b - sdet)/(2.d0*a*x + b + sdet)      
c      intindappa = 1.d0/sdet*log(abs(arg))

      arg = (2.d0*a*x+b)/sdet
      intind = -2.d0/sdet*atanh(arg)
      
      return
      end
*******************************************************************
      function intindtrue(x,sdet,a,b,c,al,be,ga)
      implicit double precision (a-i,m,n,o-z)

c      if (al.ne.0.d0) then
      m = be - b*al/a
      n = ga - c*al/a
      arg = (2.d0*a*x+b)/sdet
      intindtrue =
     .     x*al
     .     + m*0.5d0*log(abs(a*x*x + b*x + c))
     .     - (2.d0*a*n-b*m)/sdet*atanh(arg)
      intindtrue = intindtrue/a

      return
      end
***************************************************************************      
      subroutine phasespacemuepairsCRF(s,m1,m2,m3,m4,p1,p2,p3,p4,p5,p6,
     >     np,phsp,w,ier)      
! COMOVING REFERENCE FRAME
      implicit double precision (a-h,l,m,o-z)
      parameter (nphmax = 2)
      integer np
      parameter (pi = 4.d0*atan(1.d0))
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),v3(0:3),q3(0:3)
      dimension p5(0:3),p6(0:3)
      dimension csi1(1),csi2(2),csi3(3),csi4(4),csi5(5),dir4(0:3)
      dimension pe(0:3),pmu(0:4)
      dimension p1a(0:3),p2a(0:3),p3a(0:3),p4a(0:3)
      dimension p4ka(0:3), p4k(0:3)
      dimension Q(0:3),Qstar(0:3),Qp(0:3),Qpstar(0:3)
      double precision l1,l2
      dimension ptmp(0:3),pbar(0:3),dir3(0:3),ptmp2(0:3)
      double precision kone(0:3),lastk(0:3),nu,kmod,l(0:3),k2
      common/photoncutoff/egmin,egmin2
      data sbefore,egminbefore /1d12,1d12/
      common/tmintmax/tmincut,tmaxcut
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab

      logical nge2,sol1,sol2,cgt0,clt0
      
      common/tcutcommon/icut

      common/nfotonicommon/nfotoni

      common/secondphotonmaxenergy/enmaxg2

      common/icount/icount

      common/identicalandperipheral/iper,noper,iid

      common/t13rangeCRF/t13max
      data t13max/-1d20/
      
      dimension twopis(0:nphmax),factorial(0:nphmax)
      common/pispairsCRF/twopis,factorial,rmax,nmax,nmaxb,ifirst
      data ifirst /0/

      common/xcrf/x,y,modk,ck
      
      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         factorial(0) = 1.d0
         do k = 1,nphmax
            print*,'iid and iper driven in userinterface!'
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
            factorial(k) = factorial(k-1)*k
         enddo
         rmax = -100d0
         nmax = nphmax

c         iper = 1               ! peripheral diagrams in pairs
c         noper = 1
c         iid  = 0 ! identical particles in pairs

         ifirst = 1
         icount = 0
      endif
      icount = icount + 1
******mass of the pair
      if (np.eq.1) then
         m5 = m2  ! electron pair
      else
         m5 = m1  ! muon pair
      endif
      m6 = m5
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

      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      m52 = m5*m5
      m62 = m6*m6

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

      phsp = phsp * 0.25d0/sqrt(alam)

*** min and max for kone2
      mk2min = (m5+m6)*(m5+m6)
      mk2max = (sqs - m3 - m4)*(sqs - m3 - m4)

c flat      call getrnd(csi1,1)
c      an  = mk2max - mk2min
c      mk2 = an * csi1(1) + mk2min
c      w = w * an

      an  = log(mk2max/mk2min)
      call getrnd(csi1,1)
      mk2 = mk2min*exp(csi1(1)*an)
      w   = w * an*mk2
            
*****************************************
      sqmk2 = sqrt(mk2)
c      sqmk2 = m5+m6
*****************************************      
*     limiti su t24 e t13
! wrong
c      sp     = s + (m5+m6-2.d0*sqs)*(m5+m6)
c      p4m    = 0.5d0*sqrt(lambda(sp,m42,m32)/sp)
c      E4max  = sqrt(p4m**2+m42)
c      tmine  = m22+m42 - 2.d0*(p2(0)*E4max+p4m*pm)
c      tmaxe  = m22+m42 - 2.d0*(p2(0)*E4max-p4m*pm)
********
! ok      
      E3max  = 0.5d0*(s+m32-(m4+sqmk2)**2 )/sqs
      p3max  = sqrt(E3max**2-m32)
      tminmu = m12 + m32 - 2.d0*(p1(0)*E3max+p3max*pm)
      tmaxmu = m12 + m32 - 2.d0*(p1(0)*E3max-p3max*pm)
      
      E4max = 0.5d0*(s+m42-(m3+sqmk2)**2 )/sqs
      p4max = sqrt(E4max**2-m42)
      tmine = m22 + m42 - 2.d0*(p2(0)*E4max+p4max*pm)
      tmaxe = m22 + m42 - 2.d0*(p2(0)*E4max-p4max*pm)
      
***   now I generate t24 as 1/t24^2, maybe I can take into account mk2
      t24min = tmine
      t24max = min(tmaxcut,tmaxe)
      t24max = tmaxe
      if (t24min.gt.t24max) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      
cphsp 
      call getrnd(csi1,1)
      an  = 1.d0/t24min - 1.d0/t24max
      t24 = t24min / (1.d0-t24min*an*csi1(1))
      w   = w * an * t24*t24

c flat      an = t24max - t24min
c      call getrnd(csi1,1)
c      t24 = an * csi1(1) + t24min
c      w   = w * an

      sqt24 = sqrt(-t24)      
***   
      kone = 0.d0

** in the Qpstar rest frame, Breit frame [p2,p4]
      alam24 =  lambdanoabs(t24,m22,m42)
      if (alam24.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      p2a(0) =  0.5d0*sqrt(alam24)/sqt24
      p2a(1) =  0.d0
      p2a(2) =  0.d0
      p2a(3) =  0.5d0*(t24+m22-m42)/sqt24
      p4a(0) =  p2a(0)
      p4a(1) =  0.d0
      p4a(2) =  0.d0
      p4a(3) = -0.5d0*(t24+m42-m22)/sqt24
!! notice: I'm in the p4-p2 BF. here p2 goes along -z and Qp = p4-p2 along z

      Qp3 = p4a(3)-p2a(3) ! i.e. Qp3 = sqt24
** now tminmu - t24 < mk2 - 2.d0*modk*lastk(3)*Qp3 < tmaxmu - t24
** I also have (p4+kone)^2 < (sqs-m3)^2
***   [[[ and (p4+kone)^2 > s + m32 - 2*sqs*E3max ]]]

      alim = (t24 - tmaxmu + mk2) * 0.5d0/Qp3
      blim = (t24 - tminmu + mk2) * 0.5d0/Qp3
      if (alim.ge.0.d0.and.blim.ge.0.d0) then
         cgt0 = .true.
         clt0 = .false.
      elseif (alim.le.0.d0.and.blim.ge.0.d0) then
         cgt0 = .true.
         clt0 = .true.
      elseif (alim.le.0.d0.and.blim.le.0.d0) then
         cgt0 = .false.
         clt0 = .true.
      endif
      
      R     = ((sqs - m3)**2 - m42 - mk2)/p4a(0)*0.5d0      
      b4    = p4a(3)/p4a(0)
      c2lim = 1.d0/b4/b4*(1.d0-R*R/mk2)
     
      if (c2lim.lt.0.d0) then
         cmax =  1.d0
         cmin = -1.d0
         if (cgt0.and..not.clt0) cmin = 0.d0
         if (.not.cgt0.and.clt0) cmax = 0.d0         
      elseif (c2lim.gt.1.d0) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,icount,'1'
         return
      else
         if (cgt0.and.clt0) then
            cmax = 1.d0
            cmin = sqrt(c2lim)
            call getrnd(csi1,1)
            if (csi1(1).lt.0.5d0) then
               cmax = -sqrt(c2lim)
               cmin = -1.d0
            endif
            w = 2.d0*w
         elseif(cgt0) then
            cmax = 1.d0
            cmin = sqrt(c2lim)
         elseif(clt0) then
            cmax = -sqrt(c2lim)
            cmin = -1.d0
         endif
      endif
      
cphsp
      call collinearmueCRF(p2a,p4a,m22,m42,cmin,cmax,lastk,wcoll,ie)
      w = w * wcoll
c      call getrnd(csi1,1)
c      ck = (cmax-cmin)*csi1(1) + cmin
c      sk = sqrt(1.d0-ck*ck)
c      call getrnd(csi1,1)
c      phik = 2.d0*pi*csi1(1)
c      w    = w *2.d0*pi * (cmax-cmin)
c      lastk(0) = 1.d0
c      lastk(1) = sin(phik)*sk
c      lastk(2) = cos(phik)*sk
c      lastk(3) = ck

      ck = lastk(3) 
      
***   from t13 limits
      mkinf = 0.d0
      if (alim.gt.0.d0.and.blim.gt.0.d0) then
         mkinf = alim/ck
         mksup = blim/ck
      elseif (alim.le.0.d0.and.blim.gt.0.d0) then
         if (ck.gt.0.d0) then
            mksup = blim/ck
         else
            mksup = alim/ck
         endif
      elseif (alim.lt.0.d0.and.blim.lt.0.d0) then
         mkinf = blim/ck
         mksup = alim/ck
      endif
************************
      AA =  1.d0-b4*b4*ck*ck
      BB = -2.d0*R*b4*ck
      CC =  mk2-R*R
      
      disc = 4.d0*(R*R - mk2*AA) ! this should be guaranteeded > 0 with cmin cmax above
      
      mkm1 = (-BB + sqrt(disc))*0.5d0/AA
      mkm2 = (-BB - sqrt(disc))*0.5d0/AA
      sol1 = .false.
      sol2 = .false.
      if (mkm1.ge.0.d0.and.
     .    sign(1.d0,sqrt(mkm1**2+mk2)).eq.sign(1.d0,R+b4*mkm1*lastk(3)))
     .     sol1 = .true.
      if (mkm2.ge.0.d0.and.
     .    sign(1.d0,sqrt(mkm2**2+mk2)).eq.sign(1.d0,R+b4*mkm2*lastk(3)))
     .     sol2 = .true.

      if (.not.sol1.and..not.sol2) then                 
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,icount,'2'
         return
      elseif (.not.sol1.and.sol2) then
         modkmin = max(0.d0,mkinf)
         modkmax = min(mkm2,mksup)         
      elseif (sol1.and..not.sol2) then
         modkmin = max(0.d0,mkinf)
         modkmax = min(mkm1,mksup)         
      else
         modkmin = min(mkm1,mkm2)
         modkmax = max(mkm1,mkm2)
         modkmin = max(modkmin,mkinf)
         modkmax = min(modkmax,mksup)
      endif

c      modkmax = min(modkmax,sqrt((R+b4*blim)**2-mk2))
               
      if (modkmin.gt.modkmax) then         
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,icount,'3'
         return
      endif

! integrating in ek      
      ekmax = sqrt(modkmax**2 + mk2)
      ekmin = sqrt(modkmin**2 + mk2)
c flat      an    = ekmax - ekmin
c      call getrnd(csi1,1)
c      ek = an*csi1(1) + ekmin
c     w  = w * an
      an  = log(ekmax/ekmin)
      call getrnd(csi1,1)
      ek = ekmin*exp(csi1(1)*an)
      w  = w * an*ek

      modk = sqrt(abs(ek*ek-mk2))

! integrating in modk      
c      an = modkmax - modkmin
c      call getrnd(csi1,1)
c      modk = an*csi1(1) + modkmin
c      w   = w * an
c      ek  = sqrt(mk2+modk*modk)      
****************

      kone(0)   = ek
      kone(1:3) = modk*lastk(1:3)

      p4ka  = p4a + kone
      mp4k2 = dot(p4ka,p4ka)

      t13a = t24 + mk2 - 2.d0*Qp3*modk*lastk(3)
      if (t13a.lt.tminmu.or.t13a.gt.tmaxmu) then
! included above
         
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,icount,'4'
         return
      endif
      
      A = m12
      B = (t13a-m12-m32)
      C = 0.25d0*B*B + pm*pm*m32
      B = p1(0)*B
      disc = B*B-4.d0*A*C
      if (disc.ge.0.d0) then
         x1 = -B + sqrt(disc)
         x1 = x1*0.5d0/A
c         x2 = -B - sqrt(B*B-4.d0*A*C)
c         x2 = x2*0.5d0/A
         x2 = C/A/x1
      else
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      
! deve essere E3= (s+m32-mp4k2)*0.5/sqs > x2
c i.e       -mp4k2 > x2*2.d0*sqs -s -m32
c i.e       mp4k2 < s+m32-x2*2.d0*sqs
      if (mp4k2.gt.(s+m32-x2*2.d0*sqs)) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1         
c         print*,icount,'5'
         return
      endif
            
      if (mp4k2.gt.(s+m32-2.d0*sqs*m3)) then
! included above
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,icount,'6'
         return
      endif      
      
** I can now build p3 (in the CM of p1+p2)
      alam = lambdanoabs(s,m32,mp4k2)
      if (alam.lt.0.d0) then ! if the previous condition is met also this one is met
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,icount,'7'
         return
      endif      
      modp3 = 0.5d0*sqrt(alam/s)
      p3(0) = sqrt(modp3*modp3+m32)
! or  p3(0) = 0.5d0(s+m32-mp4k2)/sqs


! not using the x - y trick
!      phsp = phsp * 0.5d0*(sqt24*0.5d0)**3/ek*y*y/(1.d0-y)**4
!     notice y**2/(1.d0-y)**4 = (1+x)**2*x**2
c      phsp = phsp * 0.5d0*(sqt24*0.5d0)**3/ek*(1d0+x)**2*x*x

      phsp = phsp * modk * 0.5d0 ! *modk/ek <- this if integrating in modk
      
      t13 = dot(p2a-p4ka,p2a-p4ka)
      
c      A = m12
c      B = p1(0)*(t13-m12-m32)
c      C = 0.25d0*(t13-m12-m32)**2 + pm*pm*m32
c      disc = sqrt(B*B-4.d0*A*C)
c      if (disc.gt.0.d0) then
c         x1 = -B + sqrt(B*B-4.d0*A*C)
c         x1 = x1*0.5d0/A
c         x2 = -B - sqrt(B*B-4.d0*A*C)
c         x2 = x2*0.5d0/A
c         print*,'HERE',x1.lt.E3max,x2.lt.E3max
c      endif
                  
      c3 = t13-m12-m32 + 2.d0*p1(0)*p3(0)
      c3 = c3/pm/modp3*0.5d0
      if (abs(c3).gt.1.d0) then
!     this cuts a lot
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,icount,'8'
         return
      endif
c      print*,'after cut',p3(0).gt.x2 ! quindi E3 = 0.5*(s+m32-mp4k2)/sqs deve essere > x2, che significa
                                     ! un limite inferiore su mp4k2, che posso mettere in qualsiasi frame.
                                     ! Nota pero' che deve gia' esserci un t13 nei limiti [tminmu,tmaxmu]

c      tminmu < t13 = (p4+k-p2)^2 = t24 + k^2 - 2*(p4a(3)-p2a(3))*modk*ck < tmaxmu
c      mp4k2 = m42 + k^2 + 2*(p4a(0)*sqrt(modk**2+k^2) - p4a(3)*modk*ck) < XXX
      
      s3   = sqrt(abs(1.d0-c3*c3))
      call getrnd(csi1,1)
      phi3 = 2.d0*pi*csi1(1)
      w    = w * 2.d0*pi

      p3(1) = modp3*s3*sin(phi3)
      p3(2) = modp3*s3*cos(phi3)
      p3(3) = modp3*c3
*
      p4k   = p1+p2 - p3
*
      al_a = lambdanoabs(t24,m22,m42) ! ie al_a = alam24
      al_b = lambdanoabs(t13,m22,mp4k2)

      if (al_b.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,icount,'9'
         return
      endif

      
*      phsp = phsp * sqrt(-t24*al_a)/sqrt(-t13*al_b)
      phsp = phsp * sqrt(al_a)/sqrt(al_b) ! AAAHHRRGGG... questo funziona, c'e' un errorino da qualche parte
! the jacobian for \delta((p2-p4-k)^2-t13) is 1, because I used it to eliminate d_t13

c      l1 = -t24 + m22 - m42
c      l2 = -t24 - m22 + m42
c      l1 = l1/sqrt(lambda(t24,m22,m42))
c      l2 = l2/sqrt(lambda(t24,m22,m42))
c      Qpstar = l1 * p2a + l2 * (p4a)
c      print*,Qpstar ! ok this is the Qpstar rest frame
      
** Now I must do the various rotations and boosts and then the decay k --> p5 p6
c      p2a,p4a and kone are at this point in the Qpstar rest frame
      l1 = -t13 + m22 - mp4k2
      l2 = -t13 - m22 + mp4k2
      l1 = l1/sqrt(al_b) !sqrt(lambda(t13,m22,mp4k2))
      l2 = l2/sqrt(al_b) !sqrt(lambda(t13,m22,mp4k2))

      
! as seen in the Qpstar rest frame
      Qstar = l1 * p2a + l2 * p4ka

! I go now in the Qstar rest frame
      call new_boost(Qstar,p2a,p2a,1)
      call new_boost(Qstar,p4a,p4a,1)
      call new_boost(Qstar,kone,kone,1)
      
** ruoto dove Q = (0,0,0,sqt24) !! CRUCIALE !!
      call loadrot(1,p4a + kone - p2a)
      call rotafterload(p2a,p2a)
      call rotafterload(p4a,p4a)
      call rotafterload(kone,kone)
***
      Qstar = l1 * p2 + l2 * p4k !! p2 and p4k are in the p1+p2 CM frame!
      
      call new_boost(Qstar,p2,ptmp,1)
** ptmp e' p2 come visto in Bf-Qs. ora devo ruotare p2 in modo che coincidano e poi boostare back
      
      call getrotmatrix(p2a,ptmp)
! not needed      call rotatevec(p2a)
      call rotatevec(p4a)
      call rotatevec(kone)

c back in the p1+p2 rest frame!

!not needed      call new_boost(Qstar,p2a,p2,-1)
      call new_boost(Qstar,p4a,p4,-1)
      call new_boost(Qstar,kone,kone,-1)
      
***   decaying kone where it's at rest ***
      p5m = 0.5d0*sqrt(lambda(mk2,m52,m62)/mk2)
      call getrnd(csi1,1)
      c5 = 2.d0*csi1(1) - 1.d0
      s5 = sqrt(abs(1.d0-c5*c5))
      call getrnd(csi1,1)
      phi5 = 2.d0*pi*csi1(1)
      w    = 4.d0*pi*w
      p5(0)   =  sqrt(p5m*p5m + m52)
      p5(1)   =  p5m * sin(phi5)*s5
      p5(2)   =  p5m * cos(phi5)*s5
      p5(3)   =  p5m * c5
      p6(0)   =  sqrt(p5m*p5m + m62)
      p6(1:3) = -p5(1:3)

      phsp2b = p5m/sqrt(mk2)

!     boosting in the cm
      call new_boost(kone,p5,p5,-1)
      call new_boost(kone,p6,p6,-1)

      phsp = phsp * phsp2b * twopis(2) * 0.25d0
      
      return
      end
***************************************************************************      
***************************************************************************      
      subroutine phasespacemuepairs(s,m1,m2,m3,m4,p1,p2,p3,p4,p5,p6,np,
     >     phsp,w,ier)
** in this one I sample c4 and E4. It should be generalized easily to n-gammas      
      implicit double precision (a-h,l,m,o-z)
      parameter (nphmax = 2)
      integer np
      parameter (pi = 4.d0*atan(1.d0))
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),v3(0:3),q3(0:3)
      dimension p5(0:3),p6(0:3)
      dimension csi1(1),csi2(2),csi3(3),csi4(4),csi5(5),dir4(0:3)
      dimension pe(0:3),pmu(0:4)
      dimension ptmp(0:3),pbar(0:3),dir3(0:3),ptmp2(0:3)
      double precision kone(0:3),lastk(0:3),nu,kmod,l(0:3),k2
      common/photoncutoff/egmin,egmin2
      common/s3eps/epss3_ch1,epss3_ch2,sbefore,egminbefore
      data sbefore,egminbefore /1d12,1d12/
      common/tmintmax/tmincut,tmaxcut
      integer*8 count
      common/kappa/count
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab

      logical nge2,sol1,sol2
      
      common/tcutcommon/icut

      common/nfotonicommon/nfotoni

      common/secondphotonmaxenergy/enmaxg2

      common/icount/icount

      common/identicalandperipheral/iper,noper,iid

      common/t13range/t13max
      data t13max/-1d20/
      
      dimension twopis(0:nphmax),factorial(0:nphmax)
      common/pispairs/twopis,factorial,rmax,nmax,nmaxb,ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         factorial(0) = 1.d0
         do k = 1,nphmax
            print*,'iid and iper driven in userinterface!'
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
            factorial(k) = factorial(k-1)*k
         enddo
         rmax = -100d0
         nmax = nphmax

c         do k = 1,1000
c            print*,'CHECKING PHASE SPACE VOLUME'
c            print*,'look for ''cphsp'''
c            print*,'in cuts, main, phasespacemue, matrix_model'
c         enddo

c         iper = 1 ! peripheral diagrams in pairs
c         noper = 1
c         iid  = 0 ! identical particles in pairs
         
         ifirst = 1
         icount = 0
      endif

      nfotoni = 1 ! for cos_theta sampling
      
******mass of the pair
      if (np.eq.1) then
         m5 = m2  ! electron pair
      else
         m5 = m1  ! muon pair
      endif
      m6 = m5      
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

      m12 = m1*m1
      m22 = m2*m2
      m52 = m5*m5
      m62 = m6*m6
      
      pm = 0.5d0*sqrt(lambda(s,m12,m22)/s)
      
      p1(0) =  sqrt(pm*pm + m12)
      p1(1) =  0.d0
      p1(2) =  0.d0
      p1(3) =  pm
      
      p2(0) =  sqrt(pm*pm + m22)
      p2(1) =  0.d0
      p2(2) =  0.d0
      p2(3) = -pm

*** min and max for kone
      mk2min = (m5+m6)*(m5+m6)
      mk2max = (sqs - m3 - m4)*(sqs - m3 - m4)
***
***   I generate kone^2
      an  = log(mk2max/mk2min)
      call getrnd(csi1,1)
      mk2 = mk2min*exp(csi1(1)*an)
      w   = w * an*mk2
      
c      print*,'frozen!'
c      mk2 = 0.99d0*mk2max
      
c      an = 1.d0/mk2min - 1.d0/mk2max
c      call getrnd(csi1,1)
c      mk2 = mk2min / (1.d0 - an*csi1(1)*mk2min)
c      w = w * an*mk2*mk2
      
*** min and max for |vec{kone}| = modkmin,modkmax
      modkmin = 0.d0
      modkmax = 0.5d0*sqrt(lambda(s,(m3+m4)*(m3+m4),mk2)/s)
******
      kone = 0.d0

      e3 = sqrt(pm*pm + m12)
      e4 = sqrt(pm*pm + m22)
      
      tmax =  0.d0              !2.d0*m22 - 2.d0*e4*e4 + 2.d0*pm*pm
      tmin = -4.d0*pm*pm
      tmax = min(tmax,tmaxcut)
      tmin = max(tmin,tmincut)
      
      c4max = min( 1.d0,-1.d0-tmin*0.5d0/pm/pm)
      c4min = max(-1.d0,-1.d0-tmax*0.5d0/pm/pm)

      al    = m22 - tmax*0.5d0
      be    = p2(0)

      ilarge = 1
      call get_cth_larger(c4,c4min,c4max,s,pm,m12,m22,ilarge,wc) !attenzione, vedi considerazioni sotto
cphsp      
c phspace volume test
c      call getrnd(csi1,1)
c      c4 = 2.d0*csi1(1) - 1.d0
c      wc = 2.d0
      
      
c      anc = 1.d0/(1.d0+c4min) - 1.d0/(1.d0+c4max)
c      call getrnd(csi1,1)
c      c4 = (1.d0+c4min)/(1.d0-(1.d0+c4min)*anc*csi1(1)) - 1.d0
c      wc = anc*(1.d0+c4)**2
c mh to be checked      call get_cth_pairs(c4,c4min,c4max,s,pm,m12,m22,wc)!attenzione, vedi considerazioni sotto
      w = w * wc
      s4 = sqrt(1.d0-c4*c4)
      
      call getrnd(csi1,1)
      phi = 2.d0*pi*csi1(1)
      w   = w * 2.d0*pi

      dir4(0) = p2(0)
      dir4(1) = pm * sin(phi)*s4
      dir4(2) = pm * cos(phi)*s4
      dir4(3) = pm * c4

      dir3(0)   =  sqrt(pm*pm + m12)
      dir3(1:3) = -dir4(1:3)

** I generate the energy of kone and its direction
      ekmin = sqrt(modkmin*modkmin + mk2)
      ekmax = sqrt(modkmax*modkmax + mk2)
      
cphsp      
c phspace volume test      
c      an = ekmax - ekmin
c      call getrnd(csi1,1)
c      ek = (ekmax - ekmin)*csi1(1) + ekmin
c      w = w * an      
      an    = log(ekmax/ekmin)
      call getrnd(csi1,1)
      ek    = ekmin*exp(csi1(1)*an)
      w     = w * an*ek
      pk    = sqrt(ek*ek-mk2)

***   instead of ek, I generate pk
c      an = modkmax - modkmin
c      call getrnd(csi1,1)
c      pk = (modkmax - modkmin)*csi1(1) + modkmin
c      w = w * an

c      a   = 0.5d0 ! must be < 1.
c      uma = 1.d0 - a 
c      an = (modkmax -modkmin)**uma/uma
c      call getrnd(csi1,1)
c      pk = modkmax - (an*uma*(1.d0-csi1(1)))**(1.d0/uma)
c      w = w * an * (modkmax-pk)**a

      
c      pkma = modkmax*1.01d0
c      an = log(pkma - modkmi) - log(pkma - modkmax)
c      call getrnd(csi1,1)
c      pk = pkma - (pkma - modkmin)*exp(-an*csi1(1))
c      w = w * an * (pkma - pk)
      
      
      ek = sqrt(pk*pk+mk2)
      
      
c      an = 1.d0/ekmin - 1.d0/ekmax
c      call getrnd(csi1,1)
c      ek = ekmin / (1.d0 - an*csi1(1)*ekmin)
c      w = w * an*ek*ek

c      call getrnd(csi1,1)
c      ek = (ekmax - ekmin)*csi1(1) + ekmin
c      w  = w * (ekmax - ekmin)

c      call collinearmue_v2(p1,p2,m12,m22,dir3,dir4,lastk,wcoll,ie)
      bek = pk/ek
      call collinearpairs(p1,p2,m12,m22,dir3,dir4,
     .     ek,mk2,bek,lastk,wcoll,ie)
cphsp
c phspace volume test
c      call getrnd(csi1,1)
c      cp = 2.d0*csi1(1)-1.d0
c      sp = sqrt(1.d0-cp*cp)
c      call getrnd(csi1,1)
c      phip = 2.d0*pi*csi1(1)
c      wcoll = 4.d0*pi
c      lastk(0) = 1.d0
c      lastk(1) = sp * cos(phip)
c      lastk(2) = sp * sin(phip)
c      lastk(3) = cp

      w       = w * wcoll
      kone(0) = ek
      kone(1) = lastk(1)*pk
      kone(2) = lastk(2)*pk
      kone(3) = lastk(3)*pk
***      
      
**  calculate p4mod
      B = -2.d0*(sqs-kone(0))
      A =  (B*0.5d0)**2 - pk*pk + m22 - m12
      C = -2.d0*tridot(dir4,kone)/pm

      arg = B*B*(A*A + m22 *(C*C-B*B))
      if (arg.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      p4mod_1 = (A*C + sqrt(arg))/(B*B-C*C)
      p4mod_2 = (A*C - sqrt(arg))/(B*B-C*C)
      E4_1    = sqrt(p4mod_1*p4mod_1 + m22)
      E4_2    = sqrt(p4mod_2*p4mod_2 + m22)
      
      sol1 = .FALSE.
      sol2 = .FALSE.
      if (sign(1.d0,E4_1*B).eq.sign(1.d0,-A-C*p4mod_1)
     .     .and.p4mod_1.ge.0.d0) sol1 = .TRUE.
      if (sign(1.d0,E4_2*B).eq.sign(1.d0,-A-C*p4mod_2)
     .     .and.p4mod_2.ge.0.d0) sol2 = .TRUE.
      
      if (sol1.and..not.sol2) then
         E4 = E4_1         
      elseif (sol2.and..not.sol1) then
         E4 = E4_2
      elseif (sol1.and.sol2) then
         call getrnd(csi1,1)
         aaa = csi1(1)
         if (csi1(1).lt.0.5d0) then
            E4 = E4_1
         else
            E4 = E4_2
         endif
         w = 2.d0 * w
      else !if (.not.sol1.and..not.sol2) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif

      pm4     = sqrt(E4*E4 - m22)
      p4(0)   = E4
      p4(1:3) = dir4(1:3)/pm*pm4
      
      p3 = p1 + p2 - kone - p4
      
***   decaying kone where it's at rest ***
      p5m = 0.5d0*sqrt(lambda(mk2,m52,m62)/mk2)
      call getrnd(csi1,1)
      c5 = 2.d0*csi1(1) - 1.d0
      s5 = sqrt(1.d0-c5*c5)
      call getrnd(csi1,1)
      phi5 = 2.d0*pi*csi1(1)
      w    = 4.d0*pi*w
      p5(0)   =  sqrt(p5m*p5m + m52)
      p5(1)   =  p5m * sin(phi5)*s5
      p5(2)   =  p5m * cos(phi5)*s5
      p5(3)   =  p5m * c5
      p6(0)   =  sqrt(p5m*p5m + m62)
      p6(1:3) = -p5(1:3)

      phsp2b = p5m/sqrt(mk2)

!     boosting in the cm
      call new_boost(kone,p5,p5,-1)
      call new_boost(kone,p6,p6,-1)

c      call new_boost(kone,p5,ptmp,-1)
c      call boostmatrix(kone,p5,ptmp2,-1) ! this is new, and works
c      print*,'a',ptmp
c      print*,'b',ptmp2
c      print*,'', ptmp-ptmp2
      
!
      phsp = twopis(2)*phsp2b*0.03125d0 ! 0.03125 = 1/2^5
      E3   = p3(0)
      phsp = phsp * pm4*pm4 * pk / abs( E3*pm4 + E4*(pm4  - C*0.5d0) )
c      phsp = phsp * pk/ek ! this jacobian because I'm generating pk instead of ek

!     for identical particles
      if (iid.eq.1)  phsp = 0.5d0*phsp

! this only if iid means electrons      
      if (iid.eq.1.and.np.eq.1) then
c      if (iid.eq.0) then

         t24  = abs(dot(p2-p4,p2-p4))
         t26  = abs(dot(p2-p6,p2-p6))
         mk45 = dot(p4+p5,p4+p5)
         mk56 = dot(p6+p5,p6+p5)
                  
         e45 = p4(0)+p5(0)
         e56 = p6(0)+p5(0)
         
         anum1 = 1.d0/t24/t24/mk56/e56
     .        *(1.d0/dot(p2,p5+p6) + 1.d0/dot(p4,p5+p6))
         
         anum2 = 1.d0/t26/t26/mk45/e45
     .        *(1.d0/dot(p2,p5+p4) + 1.d0/dot(p6,p5+p4))
         
         den = 0.5d0*(anum1 + anum2)
         
         call getrnd(csi1,1)
         if (csi1(1).lt.0.5d0) then
            reg = anum1 / den
         else
            ptmp = p4
            p4   = p6
            p6   = ptmp
         
            reg = anum2 / den
         endif
         
         w = w * reg

      endif
      
      return
      end
***************************************************************************      
      subroutine phasespacemuep_phspcheck(s,m1,m2,m3,m4,p1,p2,p3,p4,
     >     p5,p6,phsp,np,w,ier)
** in this one I sample c4 and E4. It should be generalized easily to n-gammas      
      implicit double precision (a-h,l,m,o-z)
      parameter (nphmax = 2)
      parameter (pi = 4.d0*atan(1.d0))
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),v3(0:3),q3(0:3)
      dimension p5(0:3),p6(0:3)
      dimension csi1(1),csi2(2),csi3(3),csi4(4),csi5(5),dir4(0:3)
      dimension pe(0:3),pmu(0:4)
      dimension ptmp(0:3),pbar(0:3),dir3(0:3)
      double precision kone1(0:3),nu,kmod,l(0:3),k2,kone2(0:3)
      common/photoncutoff/egmin,egmin2
      common/s3eps/epss3_ch1,epss3_ch2,sbefore,egminbefore
      data sbefore,egminbefore /1d12,1d12/
      common/tmintmax/tmincut,tmaxcut
      integer*8 count
      common/kappa/count
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab

      logical nge2,sol1,sol2
      
      common/tcutcommon/icut

      common/nfotonicommon/nfotoni

      common/secondphotonmaxenergy/enmaxg2

      common/icount/icount

      common/identicalandperipheral/iper,noper,iid
      
      dimension twopis(0:nphmax),factorial(0:nphmax)
      common/pispairs/twopis,factorial,rmax,nmax,nmaxb,ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         factorial(0) = 1.d0
         do k = 1,nphmax
            print*,'iid and iper driven in userinterface'
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
            factorial(k) = factorial(k-1)*k
         enddo
         rmax = -100d0
         nmax = nphmax

         do k = 1,1000
            print*,'CHECKING PHASE SPACE VOLUME'
            print*,'look for ''cphsp'''
            print*,'in cuts, main, phasespacemue, matrix_model'
            print*,'*** this phase space is only for volume check!! ***'
         enddo

c         iper = 1 ! peripheral diagrams in pairs
c         noper = 1
c         iid  = 0 ! identical particles in pairs
         
         ifirst = 1
         icount = 0
      endif

      nfotoni = 1 ! for cos_theta sampling

******mass of the pair
      if (np.eq.1) then
         m5 = m2  ! electron pair
      else
         m5 = m1  ! muon pair
      endif
      m6 = m5      
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

      m12 = m1*m1
      m22 = m2*m2
      m52 = m5*m5
      m62 = m6*m6
      
      pm = 0.5d0*sqrt(abs(lambda(s,m12,m22)/s))
      if ((lambda(s,m12,m22)/s).le.0.d0) print*,'1-',lambda(s,m12,m22),s
      
      p1(0) =  sqrt(pm*pm + m12)
      p1(1) =  0.d0
      p1(2) =  0.d0
      p1(3) =  pm
      
      p2(0) =  sqrt(pm*pm + m22)
      p2(1) =  0.d0
      p2(2) =  0.d0
      p2(3) = -pm

*** min and max for kone1
      mk2min = (m5+m6)*(m5+m6)
      mk2max = (sqs - m3 - m4)*(sqs - m3 - m4)
***   I generate kone1^2
c      an  = log(mk2max/mk2min)
c      call getrnd(csi1,1)
c      mk2 = mk2min*exp(csi1(1)*an)
c      w   = w * an*mk2
      an    = mk2max - mk2min
      call getrnd(csi1,1)
      mk2   = an * csi1(1) + mk2min
      w = w * an
      mk21 = mk2

***   I generate kone2^2
*** min and max for kone2
      mk2min = (m3+m4)*(m3+m4)
      mk2max = (sqs - sqrt(abs(mk21)))*(sqs - sqrt(abs(mk21)))
      if (mk21.le.0.d0) print*,'2-',mk21
c      an  = log(mk2max/mk2min)
c      call getrnd(csi1,1)
c      mk2 = mk2min*exp(csi1(1)*an)
c     w   = w * an*mk2
      an = mk2max - mk2min
      call getrnd(csi1,1)
      mk2 = an * csi1(1) + mk2min
      w = w * an
      mk22 = mk2

      pk1 = 0.5d0*sqrt(abs(lambda(s,mk21,mk22)/s))
c      if (lambdanoabs(s,mk21,mk22).le.0.d0) print*,'**->',
c     .     lambdanoabs(s,mk21,mk22),mk21,mk22
      pk2 = pk1
      ek1 = sqrt(pk1*pk1 + mk21)
      ek2 = sqrt(pk2*pk2 + mk22)
      
      call getrnd(csi1,1)
      cp = 2.d0*csi1(1) - 1.d0
      sp = sqrt(abs(1.d0-cp*cp))
      call getrnd(csi1,1)
      phip = 2.d0*pi*csi1(1)
      w = w * 4.d0*pi
      kone1(0) = ek1
      kone1(1) = pk1 * sp * cos(phip)
      kone1(2) = pk1 * sp * sin(phip)
      kone1(3) = pk1 * cp

      kone2(0) = ek2
      kone2(1:3) = -kone1(1:3)

      phsp2b = pk1/sqs * 0.5d0
            
***   decaying kone1 where it's at rest ***
      p5m = 0.5d0*sqrt(abs(lambda(mk21,m52,m62)/mk21))
      call getrnd(csi1,1)
      c5 = 2.d0*csi1(1) - 1.d0
      s5 = sqrt(abs(1.d0-c5*c5))
      call getrnd(csi1,1)
      phi5 = 2.d0*pi*csi1(1)
      w    = 4.d0*pi*w
      p5(0)   =  sqrt(p5m*p5m + m52)
      p5(1)   =  p5m * sin(phi5)*s5
      p5(2)   =  p5m * cos(phi5)*s5
      p5(3)   =  p5m * c5
      p6(0)   =  sqrt(p5m*p5m + m62)
      p6(1:3) = -p5(1:3)

      phsp2b = phsp2b*p5m/sqrt(abs(mk21))
      
! boosting in the cm
      call new_boost(kone1,p5,p5,-1)
      call new_boost(kone1,p6,p6,-1)

***   decaying kone2 where it's at rest ***
      p5m = 0.5d0*sqrt(abs(lambda(mk22,m3*m3,m4*m4)/mk22))
      call getrnd(csi1,1)
      c5 = 2.d0*csi1(1) - 1.d0
      s5 = sqrt(abs(1.d0-c5*c5))
      call getrnd(csi1,1)
      phi5 = 2.d0*pi*csi1(1)
      w    = 4.d0*pi*w
      p4(0)   =  sqrt(p5m*p5m + m4*m4)
      p4(1)   =  p5m * sin(phi5)*s5
      p4(2)   =  p5m * cos(phi5)*s5
      p4(3)   =  p5m * c5
      p3(0)   =  sqrt(p5m*p5m + m3*m3)
      p3(1:3) = -p4(1:3)

      phsp2b = phsp2b*p5m/sqrt(abs(mk22))
      
! boosting in the cm
      call new_boost(kone2,p4,p4,-1)
      call new_boost(kone2,p3,p3,-1)      
!
      phsp = twopis(2)*phsp2b*0.03125d0 ! 0.03125 = 1/2^5
      
      return
      end
*************************************************************************
***************************************************************************      
      subroutine phasespacemuepairsCRF_firstworkingversion
     >     (s,m1,m2,m3,m4,p1,p2,p3,p4,p5,p6,np,phsp,w,ier)
! COMOVING REFERENCE FRAME
      implicit double precision (a-h,l,m,o-z)
      parameter (nphmax = 2)
      integer np
      parameter (pi = 4.d0*atan(1.d0))
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),v3(0:3),q3(0:3)
      dimension p5(0:3),p6(0:3)
      dimension csi1(1),csi2(2),csi3(3),csi4(4),csi5(5),dir4(0:3)
      dimension pe(0:3),pmu(0:4)
      dimension p1a(0:3),p2a(0:3),p3a(0:3),p4a(0:3)
      dimension p4ka(0:3), p4k(0:3)
      dimension Q(0:3),Qstar(0:3),Qp(0:3),Qpstar(0:3)
      double precision l1,l2
      dimension ptmp(0:3),pbar(0:3),dir3(0:3),ptmp2(0:3)
      double precision kone(0:3),lastk(0:3),nu,kmod,l(0:3),k2
      common/photoncutoff/egmin,egmin2
      data sbefore,egminbefore /1d12,1d12/
      common/tmintmax/tmincut,tmaxcut
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab

      logical nge2,sol1,sol2
      
      common/tcutcommon/icut

      common/nfotonicommon/nfotoni

      common/secondphotonmaxenergy/enmaxg2

      common/icount/icount

      common/identicalandperipheral/iper,noper,iid

      common/t13rangeCRF/t13max
      data t13max/-1d20/
      
      dimension twopis(0:nphmax),factorial(0:nphmax)
      common/pispairsCRF/twopis,factorial,rmax,nmax,nmaxb,ifirst
      data ifirst /0/

      common/xcrf/x,y,modk,ck
      
      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         factorial(0) = 1.d0
         do k = 1,nphmax
            print*,'iid and iper driven in userinterface'
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
            factorial(k) = factorial(k-1)*k
         enddo
         rmax = -100d0
         nmax = nphmax

c         iper = 1 ! peripheral diagrams in pairs
c         noper = 1
c         iid  = 0 ! identical particles in pairs
         
         ifirst = 1
         icount = 0
      endif

******mass of the pair
      if (np.eq.1) then
         m5 = m2  ! electron pair
      else
         m5 = m1  ! muon pair
      endif
      m6 = m5      
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

      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      m52 = m5*m5
      m62 = m6*m6

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
c      print*,p1(0),(s+m12-m22)/sqs*0.5d0
c      print*,p2(0),(s-m12+m22)/sqs*0.5d0

      phsp = phsp * 0.25d0/sqrt(alam)

*     limiti su t24 e t13
! wrong
c      sp     = s + (m5+m6-2.d0*sqs)*(m5+m6)
c      p4m    = 0.5d0*sqrt(lambda(sp,m42,m32)/sp)
c      E4max  = sqrt(p4m**2+m42)
c      tmine  = m22+m42 - 2.d0*(p2(0)*E4max+p4m*pm)
c      tmaxe  = m22+m42 - 2.d0*(p2(0)*E4max-p4m*pm)
********
! ok      
      E3max  = 0.5d0*(s+m32-(m4+m5+m6)**2 )/sqs
      p3max  = sqrt(E3max**2-m32)
      tminmu = m12 + m32 - 2.d0*(p1(0)*E3max+p3max*pm)
      tmaxmu = m12 + m32 - 2.d0*(p1(0)*E3max-p3max*pm)

      E4max = 0.5d0*(s+m42-(m3+m5+m6)**2 )/sqs
      p4max = sqrt(E4max**2-m42)
      tmine = m22 + m42 - 2.d0*(p2(0)*E4max+p4max*pm)
      tmaxe = m22 + m42 - 2.d0*(p2(0)*E4max-p4max*pm)
            
***   now I generate t24 as 1/t24^2, maybe I can take into account mk2
      t24min = tmine
      t24max = min(tmaxcut,tmaxe)
cphsp 
c      call getrnd(csi1,1)
c      an  = 1.d0/t24min - 1.d0/t24max
c      t24 = t24min / (1.d0-t24min*an*csi1(1))
c      w   = w * an * t24*t24
      an = t24max - t24min
      call getrnd(csi1,1)
      t24 = an * csi1(1) + t24min
      w   = w * an

      
      sqt24 = sqrt(-t24)      
***   
      kone = 0.d0

** in the Qpstar rest frame, Breit frame [p2,p4]
      alam24 =  lambdanoabs(t24,m22,m42)
      p2a(0) =  0.5d0*sqrt(alam24)/sqt24
      p2a(1) =  0.d0
      p2a(2) =  0.d0
      p2a(3) =  0.5d0*(t24+m22-m42)/sqt24
      p4a(0) =  p2a(0)
      p4a(1) =  0.d0
      p4a(2) =  0.d0
      p4a(3) = -0.5d0*(t24+m42-m22)/sqt24
!! notice: I'm in the p4-p2 BF. here p2 goes along -z and Qp = p4-p2 along z
      
**   I now generate kone


*** min and max for kone2
      mk2min = (m5+m6)*(m5+m6)
      mk2max = (sqs - m3 - m4)*(sqs - m3 - m4)

****  I try the barbatrucco, mh...
*     * I write |k| = x * sqt24*0.5
*     * x = y/(1 - y)  
      ymax = 1.d0

*     * altro test
      
      b4 = p4a(3)/p4a(0)
      AA = ((sqs-m3)**2 - m42 - mk2min)/p4a(0)*0.5d0
      c2lim = (1.d0 - AA*AA/mk2min)/b4/b4
      if (c2lim.lt.0.d0) then
         cmax =  1.d0
         cmin = -1.d0
      else
         cmax = 1.d0
         cmin = sqrt(c2lim)
         call getrnd(csi1,1)
         if (csi1(1).lt.0.5d0) then
            cmax = -sqrt(c2lim)
            cmin = -1.d0
         endif
         w = 2.d0*w
      endif

cphsp
c      call collinearmueCRF(p2a,p4a,m22,m42,lastk,wcoll,ie)
c      w = w * wcoll
      call getrnd(csi1,1)
      ck = (cmax-cmin)*csi1(1) + cmin
      sk = sqrt(1.d0-ck*ck)
      call getrnd(csi1,1)
      phik = 2.d0*pi*csi1(1)
      w    = w *2.d0*pi * (cmax-cmin)
      lastk(0) = 1.d0
      lastk(1) = sin(phik)*sk
      lastk(2) = cos(phik)*sk
      lastk(3) = ck

      BB = b4*lastk(3)
      DD = AA*AA + mk2min*(BB*BB-1.d0)
      if (DD.ge.0.d0) then
         XX1 = -AA*BB + sqrt(DD)
         XX1 = XX1/(BB*BB-1.d0)
         XX2 = -AA*BB - sqrt(DD)
         XX2 = XX2/(BB*BB-1.d0)

         if (xx2.lt.0.d0) then
            phsp = 0.d0
            w    = 0.d0
            ier  = 1
            return
         endif

         if (xx1.lt.0.d0) xx1 = 0.d0
         
! it looks like that only XX2 > 0         
         xmax = XX2*2.d0/sqt24
         ymax = xmax /(1.d0 + xmax)
      else
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
*     * fine altro test
      
      call getrnd(csi1,1)
! not using the y-x trick      
c      ymax = min(1.d0,ymax)
c      ymin = 0.d0
c      an = ymax - ymin
c      y = csi1(1)*an + ymin
c      w = w *an      
c      x = y/(1.d0-y) ! can go up to \infty  ** y = x/(1+x)
c      modk = x * sqt24 * 0.5d0
      
      modkmin = xx1
      modkmax = xx2
      
      an = modkmax - modkmin
      modk = an*csi1(1) + modkmin
      w   = w * an
      
! from the condition that at least p3 at rest must exist in the lab      
! ie  mp4k2 must be < s+m32-2.d0*sqs*m3
      u2  = modk*modk
      xi  = (sqs - m3)**2 ! s + m32 -2.d0*sqs*m3
      xip = xi - m42 + 2.d0*tridot(p4a,lastk)*modk
      ee4 = p4a(0)
      
      
!      disc = (xip + 2.d0*ee4**2)**2 - xip**2 + 4.d0*ee4**2*u2
      disc = 4.d0*ee4**2*(ee4**2 + xip + u2)
      
      fxmin = mk2min + 2.d0*ee4*sqrt(mk2min + u2)

      sol1 = .false.
      sol2 = .false.
      if (fxmin.le.xip.and.disc.ge.0.d0) then
         x1 = xip + 2.d0*ee4**2 + sqrt(disc)
c         x2 = xip + 2.d0*ee4**2 - sqrt(disc)
! or better numerical
c     x2 = (xip**2-4.d0*ee4**2*u2)/x1 
         x2 = (xip + 2.d0*ee4*modk)*(xip - 2.d0*ee4*modk)/x1
         
         if (sign(1.d0,xip-x1).eq.sign(1.d0,2.d0*ee4*sqrt(x1+u2))
     .        .and.x1.gt.mk2min) sol1 = .true.
         if (sign(1.d0,xip-x2).eq.sign(1.d0,2.d0*ee4*sqrt(x2+u2))
     .        .and.x2.gt.mk2min) sol2 = .true.

         if (sol1) mk2max = min(x1,mk2max)
         if (sol2) mk2max = min(x2,mk2max)

         uu = tminmu -t24+2.d0*(p4a(3)-p2a(3))*modk*lastk(3)
         zz = tmaxmu -t24+2.d0*(p4a(3)-p2a(3))*modk*lastk(3)
         mk2min = max(mk2min,uu)
         mk2max = min(mk2max,zz)
         if (mk2max.lt.mk2min) then
            phsp = 0.d0
            w    = 0.d0
            ier  = 1
            return
         endif         
         
***
***   I generate kone^2
cphsp 
c      an  = log(mk2max/mk2min)
c      call getrnd(csi1,1)
c      mk2 = mk2min*exp(csi1(1)*an)
c      w   = w * an*mk2
         an  = mk2max - mk2min
         call getrnd(csi1,1)
         mk2 = an*csi1(1) + mk2min
         w   = w *an
      else
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,'1'
         return
      endif         
****************
      ek   = sqrt(mk2+modk*modk)      
      
      kone(0)   = ek
      kone(1:3) = modk*lastk(1:3)
            
      p4ka  = p4a + kone
      mp4k2 = dot(p4ka,p4ka)

      t13a = t24 + mk2 - 2.d0*(p4a(3)-p2a(3))*modk*lastk(3)
      if (t13a.lt.tminmu.or.t13a.gt.tmaxmu) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif

c      print*,mp4k2,
c     . m42 + mk2 + 2.d0*(p4a(0)*sqrt(modk**2+mk2)-p4a(3)*modk*lastk(3))
      
      
      A = m12
      B = (t13a-m12-m32)
      C = 0.25d0*B*B + pm*pm*m32
      B = p1(0)*B
      disc = B*B-4.d0*A*C
      if (disc.ge.0.d0) then
         x1 = -B + sqrt(disc)
         x1 = x1*0.5d0/A
c         x2 = -B - sqrt(B*B-4.d0*A*C)
c         x2 = x2*0.5d0/A
         x2 = C/A/x1         
c         print*,'HERE',x1.lt.E3max,x2.lt.E3max
      endif
! deve essere E3= (s+m32-mp4k2)*0.5/sqs > x2
c i.e       -mp4k2 > x2*2.d0*sqs -s -m32
c i.e       mp4k2 < s+m32-x2*2.d0*sqs
      if (mp4k2.gt.(s+m32-x2*2.d0*sqs)) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
      
      
      if (mp4k2.gt.(s+m32-2.d0*sqs*m3)) then
! included above         
! should not happen, because is checked before, but sometimes it happens...         
c         print*,mp4k2/(s+m32-2.d0*sqs*m3)
c         print*,(xip + 2.d0*ee4**2)/sqrt(disc)

! from the condition that at least p3 at rest must exist in the lab

c         print*,sol1,sol2,fxmin.lt.xip,disc
         
c     print*,x,y,lastk(3),sqrt(mk2)/sqs
         
c         u2 = modk*modk
c         al = s+m32-2.d0*sqs*m3 - m42 + 2.d0*tridot(p4a,lastk)*modk
c         print*,''
c         print*,2.d0*p4a(0)*sqrt(mk2+u2).lt.(al-mk2)
c         print*,2.d0*p4a(0)*sqrt(mk2min+u2).lt.(al-mk2min)
c         print*,2.d0*p4a(0)*sqrt(mk2max+u2).lt.(al-mk2max)

         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,'2'
         return
      endif      
      
** I can now build p3 (in the CM of p1+p2)
      alam = lambdanoabs(s,m32,mp4k2)
      if (alam.lt.0.d0) then ! if the previous condition is met also this one is met
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,'3'
         return
      endif      
      modp3 = 0.5d0*sqrt(alam/s)
      p3(0) = sqrt(modp3*modp3+m32)
! or  p3(0) = 0.5d0(s+m32-mp4k2)/sqs


! not using the x - y trick
!      phsp = phsp * 0.5d0*(sqt24*0.5d0)**3/ek*y*y/(1.d0-y)**4
!     notice y**2/(1.d0-y)**4 = (1+x)**2*x**2
c      phsp = phsp * 0.5d0*(sqt24*0.5d0)**3/ek*(1d0+x)**2*x*x

      phsp = phsp * modk*modk/ek * 0.5d0
      
      t13 = dot(p2a-p4ka,p2a-p4ka)
      
c      A = m12
c      B = p1(0)*(t13-m12-m32)
c      C = 0.25d0*(t13-m12-m32)**2 + pm*pm*m32
c      disc = sqrt(B*B-4.d0*A*C)
c      if (disc.gt.0.d0) then
c         x1 = -B + sqrt(B*B-4.d0*A*C)
c         x1 = x1*0.5d0/A
c         x2 = -B - sqrt(B*B-4.d0*A*C)
c         x2 = x2*0.5d0/A
c         print*,'HERE',x1.lt.E3max,x2.lt.E3max
c      endif
                  
      c3 = t13-m12-m32 + 2.d0*p1(0)*p3(0)
      c3 = c3/pm/modp3*0.5d0
      if (abs(c3).gt.1.d0) then
!     this cuts a lot
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif
c      print*,'after cut',p3(0).gt.x2 ! quindi E3 = 0.5*(s+m32-mp4k2)/sqs deve essere > x2, che significa
                                     ! un limite inferiore su mp4k2, che posso mettere in qualsiasi frame.
                                     ! Nota pero' che deve gia' esserci un t13 nei limiti [tminmu,tmaxmu]


c      tminmu < t13 = (p4+k-p2)^2 = t24 + k^2 - 2*(p4a(3)-p2a(3))*modk*ck < tmaxmu
c      mp4k2 = m42 + k^2 + 2*(p4a(0)*sqrt(modk**2+k^2) - p4a(3)*modk*ck) < XXX
      
      s3 = sqrt(1.d0-c3*c3)
      call getrnd(csi1,1)
      phi3 = 2.d0*pi*csi1(1)
      w    = w * 2.d0*pi

      p3(1) = modp3*s3*sin(phi3)
      p3(2) = modp3*s3*cos(phi3)
      p3(3) = modp3*c3
*
      p4k   = p1+p2 - p3
*
      al_a = lambdanoabs(t24,m22,m42)
      al_b = lambdanoabs(t13,m22,mp4k2)

      if (al_b.lt.0.d0) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
c         print*,'5'
         return
      endif

c      disc = (m12+m32-t13)**2 - 4.d0*m32*(p1(0)**2-p1(3)**2*c3**2)
c      disc = disc * p1(3)**2*c3**2
c      uu1 = p1(0)*(m12+m32-t13)+sqrt(disc)
c      uu1 = uu1/(p1(0)**2-p1(3)**2*c3**2)*0.5d0
c      uu2 = p1(0)*(m12+m32-t13)-sqrt(disc)
c      uu2 = uu2/(p1(0)**2-p1(3)**2*c3**2)*0.5d0      
!disc > 0 --> 0.25d0*(m12+m32-t13)**2 > m32*(p1(0)**2-p1(3)**2*c3**2)
!     p1(0)**2 - 0.25d0*(m12+m32-t13)**2/m32 < p1(3)**2*c3**2
!     c3**2 > (p1(0)**2 - 0.25d0*(m12+m32-t13)**2/m32)/p1(3)**2
c      print*,''
c      print*,E3max,p3(0)
c      print*,uu1,uu2
c      print*,t13-m12-m32+2.d0*p1(0)*uu1,2.d0*p1(3)*sqrt(uu1**2-m32)*c3
c      print*,t13-m12-m32+2.d0*p1(0)*uu2,2.d0*p1(3)*sqrt(uu2**2-m32)*c3
c      print*,c3**2,(p1(0)**2 - 0.25d0*(m12+m32-t13)**2/m32)/p1(3)**2
      
*      phsp = phsp * sqrt(-t24*al_a)/sqrt(-t13*al_b)
      phsp = phsp * sqrt(al_a)/sqrt(al_b) ! AAAHHRRGGG... questo funziona, c'e' un errorino da qualche parte
! the jacobian for \delta((p2-p4-k)^2-t13) is 1, because I used it to eliminate d_t13

c      l1 = -t24 + m22 - m42
c      l2 = -t24 - m22 + m42
c      l1 = l1/sqrt(lambda(t24,m22,m42))
c      l2 = l2/sqrt(lambda(t24,m22,m42))
c      Qpstar = l1 * p2a + l2 * (p4a)
c      print*,Qpstar ! ok this is the Qpstar rest frame
      
** Now I must do the various rotations and boosts and then the decay k --> p5 p6
c      p2a,p4a and kone are at this point in the Qpstar rest frame
      l1 = -t13 + m22 - mp4k2
      l2 = -t13 - m22 + mp4k2
      l1 = l1/sqrt(al_b) !sqrt(lambda(t13,m22,mp4k2))
      l2 = l2/sqrt(al_b) !sqrt(lambda(t13,m22,mp4k2))

      
! as seen in the Qpstar rest frame
      Qstar = l1 * p2a + l2 * p4ka

! I go now in the Qstar rest frame
      call new_boost(Qstar,p2a,p2a,1)
      call new_boost(Qstar,p4a,p4a,1)
      call new_boost(Qstar,kone,kone,1)
      
** ruoto dove Q = (0,0,0,sqt24) !! CRUCIALE !!
      call loadrot(1,p4a + kone - p2a)
      call rotafterload(p2a,p2a)
      call rotafterload(p4a,p4a)
      call rotafterload(kone,kone)
***
      Qstar = l1 * p2 + l2 * p4k !! p2 and p4k are in the p1+p2 CM frame!
      
      call new_boost(Qstar,p2,ptmp,1)
** ptmp e' p2 come visto in Bf-Qs. ora devo ruotare p2 in modo che coincidano e poi boostare back
      
      call getrotmatrix(p2a,ptmp)
! not needed      call rotatevec(p2a)
      call rotatevec(p4a)
      call rotatevec(kone)

c back in the p1+p2 rest frame!

!not needed      call new_boost(Qstar,p2a,p2,-1)
      call new_boost(Qstar,p4a,p4,-1)
      call new_boost(Qstar,kone,kone,-1)
      
***   decaying kone where it's at rest ***
      p5m = 0.5d0*sqrt(lambda(mk2,m52,m62)/mk2)
      call getrnd(csi1,1)
      c5 = 2.d0*csi1(1) - 1.d0
      s5 = sqrt(1.d0-c5*c5)
      call getrnd(csi1,1)
      phi5 = 2.d0*pi*csi1(1)
      w    = 4.d0*pi*w
      p5(0)   =  sqrt(p5m*p5m + m52)
      p5(1)   =  p5m * sin(phi5)*s5
      p5(2)   =  p5m * cos(phi5)*s5
      p5(3)   =  p5m * c5
      p6(0)   =  sqrt(p5m*p5m + m62)
      p6(1:3) = -p5(1:3)

      phsp2b = p5m/sqrt(mk2)

!     boosting in the cm
      call new_boost(kone,p5,p5,-1)
      call new_boost(kone,p6,p6,-1)

      phsp = phsp * phsp2b * twopis(2) * 0.25d0

      return
      end
***************************************************************************      
** another phase space for pairs, it's now the default
***************************************************************************      
      subroutine phasespacemuepairsv3(s,m1,m2,m3,m4,p1,p2,p3,p4,p5,p6,
     >     np,phsp,w,ier)      
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
      dimension p4ka(0:3), p4k(0:3)
      dimension Q(0:3),Qstar(0:3),Qp(0:3),Qpstar(0:3)
      double precision l1,l2
      dimension ptmp(0:3),pbar(0:3),dir3(0:3),ptmp2(0:3)
      double precision kone(0:3),lastk(0:3),nu,kmod,l(0:3),k2
      common/photoncutoff/egmin,egmin2
      common/tmintmax/tmincut,tmaxcut

      common/mueexpsetup/emulab,eemin,semu,thmumin,themax,thmumax,
     .     dthna7max,cutela,ina7,iela
      
      double precision Qmu
      common/muoncharge/Qmu
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab

      logical nge2,sol1,sol2,cgt0,clt0
      
      common/tcutcommon/icut

      common/nfotonicommon/nfotoni

      common/secondphotonmaxenergy/enmaxg2

      common/icount/icount

      common/identicalandperipheral/iper,noper,iid

      dimension twopis(0:nphmax),factorial(0:nphmax)
      common/pispairs_v3/twopis,factorial,rmax,nmax,nmaxb,ifirst
      data ifirst /0/

      common/xcrf/x,y,modk,ck
      
      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         factorial(0) = 1.d0
         do k = 1,nphmax
c            print*,'iid and iper driven in userinterface'
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
            factorial(k) = factorial(k-1)*k
         enddo
         rmax = -100d0
         nmax = nphmax

c         iper = 1 ! peripheral diagrams in pairs
c         noper = 1
c         iid  = 0 ! identical particles in pairs
         do k = 1,100
cc            print*,'iid = 0 and sampling flat and summt2 = 1000!!'
         enddo
            
         ifirst = 1
         icount = 0
      endif
      icount = icount + 1
******mass of the pair
      if (np.eq.1) then
         m5 = m2  ! electron pair
      else
         m5 = m1  ! muon pair
      endif
      m6 = m5
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

      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      m52 = m5*m5
      m62 = m6*m6

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
      
****  min and max for (p4+p5+p6)**2
      m4562min = (m4+m5+m6)*(m4+m5+m6)
      m4562max = (sqs-m3)*(sqs-m3)
****
*** this must be sampled as 1/(m4562 - m22)
c flat
c      call getrnd(csi1,1)
c      an    = m4562max - m4562min
c      m4562 = an * csi1(1) + m4562min
c      w     = w * an
c  end flat
c as 1/m4562  DEFAULT
      an    = log( (m4562max-m42) / (m4562min-m42) )
      call getrnd(csi1,1)
      m4562 = (m4562min-m42)*exp(csi1(1)*an) + m42
      w     = w * an * (m4562 - m42)
c as 1/m4562^2 mh, better no
c      call getrnd(csi1,1)
c      an  = 1.d0/m4562min - 1.d0/m4562max
c      m4562 = m4562min / (1.d0-m4562min*an*csi1(1))
c      w   = w * an * m4562*m4562

      
****
      E3  = 0.5d0*(s + m32 - m4562 )/sqs
      p3m = sqrt(E3*E3 - m32)
      
*t13 limits
      t13min = m12 + m32 - 2.d0*(p1(0)*E3 + p3m*pm)
      t13max = m12 + m32 - 2.d0*(p1(0)*E3 - p3m*pm)
     
      if (thmumin.gt.1d-12.and.np.eq.1.and.iid.eq.0) then
c      if (thmumin.gt.1d12) then
**** [[ implementing the thmumin cut on t13
** it's from the condition |pmu_cm|*sin(thmucm) = |pmu_lab|*sin(thmulab)
** and |pmu_lab|*cos(thmulab) = gb * (|pmu_cm|*cos(thmucm) - bb*emucm)      
         tgl = tan(thmumin*1d-3)
         bb  = -p12lab(3)/p12lab(0) ! boost beta
         gb  = 1.d0/sqrt(1.d0-bb*bb) ! boost gamma
         bcm = p3m/E3           ! muon beta in the cm
         
         A =  gb*gb*tgl*tgl + 1.d0
         B = -gb*gb*tgl*tgl * bb/bcm
         C =  gb*gb*tgl*tgl*bb/bcm*bb/bcm - 1.d0
         
         disc = B*B-A*C
         if (disc.lt.0.d0) then
            phsp = 0.d0
            w    = 0.d0
            ier  = 1
            return
         endif
         
         uu1 = -B + sqrt(disc)
         uu1 = uu1 /A
         uu2 = -B - sqrt(disc)
         uu2 = uu2 /A
         
         t13max = m12 + m32 - 2.d0*(p1(0)*E3 - uu1*p3m*pm)         
         if (t13min.gt.t13max) then
            phsp = 0.d0
            w    = 0.d0
            ier  = 1
            return
         endif
***   ]]
      endif

** t13 must be sampled as 1/t13^2
c flat      
c      call getrnd(csi1,1)
c      an  = t13max - t13min
c      t13 = an * csi1(1) + t13min
c      w   = w * an
c end flat
c as 1/t13^2
c      call getrnd(csi1,1)
c      an  = 1.d0/t13min - 1.d0/t13max
c      t13 = t13min / (1.d0-t13min*an*csi1(1))
c      w   = w * an * t13*t13
c as 1/t13 THIS SEEMS TO BE BETTER, DEFAULT
      an = log(t13max/t13min) 
      call getrnd(csi1,1)
      t13 = t13min*exp(csi1(1)*an)
      w = w * an *t13
     
      c3 = t13 - m12 - m32 + 2.d0*p1(0)*E3
      c3 = c3*0.5d0/p3m/pm
      s3 = sqrt(abs(1.d0-c3*c3))

      call getrnd(csi1,1)
      phi3 = 2.d0*pi*csi1(1)
      w    = w * 2.d0*pi

      p3(0) = E3
      p3(1) = p3m * s3 * sin(phi3)
      p3(2) = p3m * s3 * cos(phi3)
      p3(3) = p3m * c3

*****
      phsp = phsp * 0.25d0/sqrt(alam)
***** 

      q13 = p1 - p3
      
*** now I'm in the q13 + p2 RF
      sp   = dot(q13+p2,q13+p2)
      sqsp = sqrt(sp)
      
      mk2min = (m5+m6)*(m5+m6)
      mk2max = (sqsp - m4)*(sqsp - m4)
      if (mk2min.ge.mk2max) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif

! as 1/mk2      
c      call getrnd(csi1,1)
c      an  = log(mk2max/mk2min)
c      mk2 = mk2min*exp(csi1(1)*an)
c      w   = w * an * mk2
! DEFAULT as 1/mk2**2, seems to be better (pay attention to the regulator below)
      an  = 1.d0/mk2min - 1.d0/mk2max
      mk2 = mk2min / (1.d0-mk2min*an*csi1(1))
      w   = w * an * mk2*mk2     
c flat      
c      call getrnd(csi1,1)      
c      an  = mk2max - mk2min
c      mk2 = an*csi1(1) + mk2min
c      w   = w * an
c end flat
      
      E4  = (sp + m42 - mk2)*0.5d0/sqsp
      p4m = sqrt(E4*E4-m42)

      ek   = (sp - m42 + mk2)*0.5d0/sqsp   
      modk = sqrt(ek*ek-mk2)

***
      phsp = phsp * modk/sqsp * 0.5d0
***
      
** now I bring q13 and p2 where q13+p2 is at rest. Notice that here q13a and p2a are always back-to-back.
      call new_boost(q13+p2,q13,q13a,1)
      call new_boost(q13+p2,p2,p2a,1)

** DEFAULT ck to be sampled later as 1/((p2a-k)^2-m42) and/or 1/(q13a-k)^2
      A = mk2-2.d0*p2a(0)*ek
      A = A / 2.d0 /sqrt(tridot(p2a,p2a))/modk
      B = t13 + mk2-2.d0*q13a(0)*ek
      B = B / 2.d0 /sqrt(tridot(q13a,q13a))/modk
! originally      call anothercollinear(abs(A),abs(B),-0.5d0,q13a,p2a,lastk,wcoll)
      call anothercollinear(abs(A),abs(B),0.5d0,q13a,p2a,lastk,wcoll)
      w         = w * wcoll
      kone(0)   = ek
      kone(1:3) = lastk(1:3)*modk
c flat
c      call getrnd(csi1,1)      
c      an  = 2.d0
c      ck = an*csi1(1) + (-1.d0)
c      sk = sqrt(abs(1.d0-ck*ck))
c      w   = w * an
c      call getrnd(csi1,1)      
c      phik = 2.d0*pi*csi1(1)
c      w    = w * 2.d0*pi
c      sphik = sin(phik)
c      cphik = cos(phik)
c      kone(0) = ek
c      kone(1) = modk * sk * sphik
c      kone(2) = modk * sk * cphik
c      kone(3) = modk * ck
c end flat
           
      p4a(0)   =  E4
      p4a(1:3) = -kone(1:3)

******
      q24a = p2a - p4a

*** now I bring q13a and q24a in the kone RF. Notice that here q13a and q24a are always back-to-back.
      call new_boost(kone,q13a,q13a,1)
      call new_boost(kone,q24a,q24a,1)
      
      sqmk2 = sqrt(mk2)
      p5(0) = sqmk2*0.5d0
      p6(0) = p5(0)
      p5m   = sqrt(mk2*0.25d0-m52)

      IF (iper.gt.-1) THEN ! i.e always
cc      IF (iper.ge.1) THEN
cccc      IF (iper.gt.10) THEN ! i.e never
cc here A = B !      
         A = dot(q13a,q13a) - 2.d0*q13a(0)*p5(0)
         A = A * 0.5d0 /sqrt(tridot(q13a,q13a))/p5m
         B = dot(q24a,q24a) - 2.d0*q24a(0)*p5(0)
         B = B * 0.5d0 /sqrt(tridot(q24a,q24a))/p5m
         call anothercollinear(abs(A),abs(B),0.5d0,q13a,q24a,lastk,
     .        wcoll)
         w = w * wcoll

*** symmetrizing p5 and p6      
         call getrnd(csi1,1)      
         if (csi1(1).le.0.5d0) then      
            p5(1:3) =  p5m * lastk(1:3)
            p6(1:3) = -p5(1:3)
            anum = 1.d0/abs(dot(q13a-p5,q13a-p5))
     .           + 1.d0/abs(dot(q24a-p5,q24a-p5))
            den  = anum +
     .             1.d0/abs(dot(q13a-p6,q13a-p6))
     .           + 1.d0/abs(dot(q24a-p6,q24a-p6))
         else
            p6(1:3) =  p5m * lastk(1:3)
            p5(1:3) = -p6(1:3)
            anum = 1.d0/abs(dot(q13a-p6,q13a-p6))
     .           + 1.d0/abs(dot(q24a-p6,q24a-p6))
            den = anum +
     .             1.d0/abs(dot(q13a-p5,q13a-p5))
     .           + 1.d0/abs(dot(q24a-p5,q24a-p5))
         endif
         w = w * 2.d0*anum/den
      ELSE
         call getrnd(csi1,1)
         an  = 2.d0
         c5 = an*csi1(1) + (-1.d0)
         s5 = sqrt(abs(1.d0-c5*c5))
         w   = w * an
         call getrnd(csi1,1)      
         phi5 = 2.d0*pi*csi1(1)
         w    = w * 2.d0*pi
         sphi5 = sin(phi5)
         cphi5 = cos(phi5)
         p5(1) = p5m * sphi5 * s5
         p5(2) = p5m * cphi5 * s5
         p5(3) = p5m * c5
         p6(1:3) = -p5(1:3)
      ENDIF
*** c5 or c6 later to be sampled according to 1/((q24a-p5[6])^2-m5[6]2) and/or 1/((q13a-p5[6])^2-m5[6]2)
c flat
c      call getrnd(csi1,1)      
c      an  = 2.d0
c      c5 = an*csi1(1) + (-1.d0)
c      s5 = sqrt(abs(1.d0-c5*c5))
c      w   = w * an
c      call getrnd(csi1,1)      
c      phi5 = 2.d0*pi*csi1(1)
c      w    = w * 2.d0*pi
c      sphi5 = sin(phi5)
c      cphi5 = cos(phi5)
c      p5(1) = p5m * sphi5 * s5
c      p5(2) = p5m * cphi5 * s5
c      p5(3) = p5m * c5
c      p5(1:3) = p5m * lastk(1:3)
c      p6(1:3) = -p5(1:3)
c end flat

****
      phsp = phsp * p5m/sqmk2 * 0.5d0
****      

*** now p5 and p6 back in the q13+p2 RF, ie where kone is kone      
      call new_boost(kone,p5,p5,-1)
      call new_boost(kone,p6,p6,-1)

*** and now back in the p1+p2 RF
      call new_boost(q13+p2,p4a,p4,-1)
      call new_boost(q13+p2,p5,p5,-1)
      call new_boost(q13+p2,p6,p6,-1)
***
      phsp = phsp * twopis(2)*0.25d0
*************************
      
      if (iid.gt.0) then
****
         phsp = 0.5d0*phsp
****         
         if (np.eq.1) then
!     for electron pairs

            mk2_56 = dot(p5+p6,p5+p6)
            mk2_45 = dot(p4+p5,p4+p5)
            f1_56 = 1.d0/mk2_56 **2 ! **2 if sampling as 1/mk2^2
            f1_45 = 1.d0/mk2_45 **2 
            
            a_56 = dot(p2-p5-p6,p2-p5-p6) - m42
            a_45 = dot(p2-p5-p4,p2-p5-p4) - m42
c            b_56 = dot(p1-p3-p5-p6,p1-p3-p5-p6)
c            b_45 = dot(p1-p3-p5-p4,p1-p3-p5-p4) ! i.e.
            b_56 = dot(p4-p2,p4-p2)
            b_45 = dot(p6-p2,p6-p2)
            f2_56 = 1.d0/a_56 + 1.d0/b_56
            f2_45 = 1.d0/a_45 + 1.d0/b_45

            f3_56 = 1.d0/(dot(p2-p4-p5,p2-p4-p5)-m52)
     .            + 1.d0/(dot(p2-p4-p6,p2-p4-p6)-m62)
     .            + 1.d0/(dot(p1-p3-p5,p1-p3-p5)-m52)
     .            + 1.d0/(dot(p1-p3-p6,p1-p3-p6)-m62)      

            f3_45 = 1.d0/(dot(p2-p6-p5,p2-p6-p5)-m52)
     .            + 1.d0/(dot(p2-p6-p4,p2-p6-p4)-m62)
     .            + 1.d0/(dot(p1-p3-p5,p1-p3-p5)-m52)
     .            + 1.d0/(dot(p1-p3-p4,p1-p3-p4)-m62)      
            
            anum = f1_56 * f2_56 * f3_56
            den  = anum + f1_45 * f2_45 * f3_45
            call getrnd(csi1,1)
            if (csi1(1).le.0.5d0) then
               ptmp = p4
               p4   = p6
               p6   = ptmp
            endif
            w = w * anum/den * 2.d0
            
         elseif(np.eq.2) then
!     for muon pairs
            if (Qmu.lt.0.d0) then
               t13 = dot(p1-p3,p1-p3)
               t16 = dot(p1-p6,p1-p6)

               anum = 1.d0/t13!**2
               den  = anum + 1.d0/t16!**2
               call getrnd(csi1,1)
               if (csi1(1).le.0.5d0) then
                  ptmp = p3
                  p3   = p6
                  p6   = ptmp
               endif
               w = w * anum/den * 2.d0
            else
               t13 = dot(p1-p3,p1-p3)
               t15 = dot(p1-p5,p1-p5)
               anum = 1.d0/t13!**2
               den  = anum + 1.d0/t15!**2
               call getrnd(csi1,1)
               if (csi1(1).le.0.5d0) then
                  ptmp = p3
                  p3   = p5
                  p5   = ptmp
               endif
               w = w * anum/den * 2.d0
            endif
         endif
      endif
*************************      
      return
      end
***************************************************************************      
      subroutine phasespacemuepairsv4(s,m1,m2,m3,m4,p1,p2,p3,p4,p5,p6,
     >     np,phsp,w,ier)
! derived from phasespacemuepairsv3
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

      common/mueexpsetup/emulab,eemin,semu,thmumin,themax,thmumax,
     .     dthna7max,cutela,ina7,iela
      
      double precision Qmu
      common/muoncharge/Qmu
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab
      logical nge2,sol1,sol2,cgt0,clt0
      common/tcutcommon/icut
      common/nfotonicommon/nfotoni
      common/secondphotonmaxenergy/enmaxg2
      common/icount/icount
      common/identicalandperipheral/iper,noper,iid
      dimension twopis(0:nphmax),factorial(0:nphmax)
      common/pispairs_v3/twopis,factorial,rmax,nmax,nmaxb,ifirst
      data ifirst /0/

      common/xcrf/x,y,modk,ck
      
      if (ifirst.eq.0) then
         twopis(0) = (2.d0*pi)**(4-6)
         factorial(0) = 1.d0
         do k = 1,nphmax
c            print*,'iid and iper driven in userinterface'
            twopis(k) = twopis(k-1) * (2.d0*pi)**(-3)
            factorial(k) = factorial(k-1)*k
         enddo
         rmax = -100d0
         nmax = nphmax
            
         ifirst = 1
         icount = 0
      endif
      icount = icount + 1
******mass of the pair
      if (np.eq.1) then
         m5 = m2  ! electron pair
      else
         m5 = m1  ! muon pair
      endif
      m6 = m5
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

      m12 = m1*m1
      m22 = m2*m2
      m32 = m3*m3
      m42 = m4*m4
      m52 = m5*m5
      m62 = m6*m6

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
      pch1 = 1.d0
      pch2 = 1.d0 - pch1
      ich  = 1                   ! default, t2[46] is anyway sampled...
      if (np.eq.1) then
         call getrnd(csi1,1)
         ich = 1
         pch1 = 0.005d0
         if (iid.eq.0) pch1 = 1.d0
         pch2 = 1.d0 - pch1
         if (csi1(1).gt.pch1) ich = 2
      endif
         
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
*** this must be sampled as 1/(m4562 - m22)
c flat
c      call getrnd(csi1,1)
c      an    = m4562max - m4562min
c      m4562 = an * csi1(1) + m4562min
c      w     = w * an
c  end flat
c as 1/m4562  DEFAULT
      an    = log( (m4562max-m4f2) / (m4562min-m4f2) )
      call getrnd(csi1,1)
      m4562 = (m4562min-m4f2)*exp(csi1(1)*an) + m4f2
      w     = w * an * (m4562 - m4f2)
      
c as 1/m4562^2 mh, better no
c      call getrnd(csi1,1)
c      an  = 1.d0/m4562min - 1.d0/m4562max
c      m4562 = m4562min / (1.d0-m4562min*an*csi1(1))
c      w   = w * an * m4562*m4562

****
** its the energy in the cm of the t13 or t24      
      Es  = 0.5d0*(s + m3f2 - m4562 )/sqs
      psm = sqrt(Es*Es - m3f2)
      
** t_sampling limits
      if (ich.eq.1) then
         tsmin = m12 + m3f2 - 2.d0*(p1(0)*Es + psm*pm)
         tsmax = m12 + m3f2 - 2.d0*(p1(0)*Es - psm*pm)

         if (thmumin.gt.1d-12.and.np.eq.1) then
**** [[ implementing the thmumin cut on t13
** it's from the condition |pmu_cm|*sin(thmucm) = |pmu_lab|*sin(thmulab)
** and |pmu_lab|*cos(thmulab) = gb * (|pmu_cm|*cos(thmucm) - bb*emucm)      
            tgl = tan(thmumin*1d-3)
            bb  = -p12lab(3)/p12lab(0) ! boost beta
            gb  = 1.d0/sqrt(1.d0-bb*bb) ! boost gamma
            bcm = psm/Es        ! muon beta in the cm
            
            A =  gb*gb*tgl*tgl + 1.d0
            B = -gb*gb*tgl*tgl * bb/bcm
            C =  gb*gb*tgl*tgl*bb/bcm*bb/bcm - 1.d0
            
            disc = B*B-A*C
            if (disc.lt.0.d0) then
               phsp = 0.d0
               w    = 0.d0
               ier  = 1
               return
            endif
            
            uu1 = -B + sqrt(disc)
            uu1 = uu1 /A
            uu2 = -B - sqrt(disc)
            uu2 = uu2 /A
            
            tsmax = m12 + m3f2 - 2.d0*(p1(0)*Es - uu1*psm*pm)         
            if (tsmin.gt.tsmax) then
               phsp = 0.d0
               w    = 0.d0
               ier  = 1
               return
            endif
***   ]]
         endif
      else
         tsmin = m22 + m3f2 - 2.d0*(p2(0)*Es + psm*pm)
         tsmax = m22 + m3f2 - 2.d0*(p2(0)*Es - psm*pm)
      endif

** t13 must be sampled as 1/t13^2
c flat      
c      call getrnd(csi1,1)
c      an  = tsmax - tsmin
c      ts = an * csi1(1) + tsmin
c      w   = w * an
c end flat
c as 1/t13^2
c      call getrnd(csi1,1)
c      an  = 1.d0/tsmin - 1.d0/tsmax
c      ts = tsmin / (1.d0-tsmin*an*csi1(1))
c      w   = w * an * ts*ts

c as 1/t13 THIS SEEMS TO BE BETTER, DEFAULT
      an = log(tsmax/tsmin) 
      call getrnd(csi1,1)
      ts = tsmin*exp(csi1(1)*an)
      w = w * an * ts

c      an1 = 1.d0/tsmin - 1.d0/tsmax
c      an2 = log(tsmax/tsmin) 
c      ptta = 0.25d0
c      pttb = 1.d0 - ptta
c      call getrnd(csi1,1)
c      if (csi1(1).lt.ptta) then
c         call getrnd(csi1,1)
c         ts = tsmin / (1.d0-tsmin*an1*csi1(1))
c      else        
c         call getrnd(csi1,1)
c         ts = tsmin*exp(csi1(1)*an2)
c      endif
c      w = w * 1.d0/(ptta/an1/ts/ts + pttb/an2/ts)

c      an1 = tsmax - tsmin
c      an2 = log(tsmax/tsmin) 
c      pta = 0.0d0
c      ptb = 1.d0 - pta
c      call getrnd(csi1,1)
c      if (csi1(1).lt.pta) then
c         call getrnd(csi1,1)
c         ts = an1 * csi1(1) + tsmin
c      else        
c         call getrnd(csi1,1)
c         ts = tsmin*exp(csi1(1)*an2)
c      endif
c      w = w * 1.d0/(pta/an1 + ptb/an2/ts)
      
      
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
      
      mk2min = (m5+m6)*(m5+m6)
      mk2max = (sqsp - m4f)*(sqsp - m4f)

      if (mk2min.ge.mk2max) then
         phsp = 0.d0
         w    = 0.d0
         ier  = 1
         return
      endif

! as 1/mk2 DEFAULT
      call getrnd(csi1,1)
      an  = log(mk2max/mk2min)
      mk2 = mk2min*exp(csi1(1)*an)
      w   = w * an * mk2
! as 1/mk2**2, seems to be better (pay attention to the regulator below)
c      call getrnd(csi1,1)
c      an  = 1.d0/mk2min - 1.d0/mk2max
c      mk2 = mk2min / (1.d0-mk2min*an*csi1(1))
c      w   = w * an * mk2*mk2     
c flat      
c      call getrnd(csi1,1)      
c      an  = mk2max - mk2min
c      mk2 = an*csi1(1) + mk2min
c      w   = w * an
c end flat

c      prm1 = 0.2d0
c      prm2 = 1.d0 - prm1
c      an1 = mk2max - mk2min
c      an2 = log(mk2max/mk2min)
c      call getrnd(csi1,1)
c      if (csi1(1).lt.prm1) then
c         call getrnd(csi1,1)
c         mk2 = an1*csi1(1) + mk2min
c      else
c         call getrnd(csi1,1)
c         mk2 = mk2min*exp(csi1(1)*an2)
c      endif
c      w = w / (prm1/an1 + prm2/an2/mk2)

      
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
      call anothercollinear(abs(A),abs(B),0.75d0,q13a,p2a,lastk,wcoll)

      w         = w * wcoll
      kone(0)   = ek
      kone(1:3) = lastk(1:3)*modk
c flat
c      call getrnd(csi1,1)      
c      an  = 2.d0
c      ck = an*csi1(1) + (-1.d0)
c      sk = sqrt(abs(1.d0-ck*ck))
c      w   = w * an
c      call getrnd(csi1,1)      
c      phik = 2.d0*pi*csi1(1)
c      w    = w * 2.d0*pi
c      sphik = sin(phik)
c      cphik = cos(phik)
c      kone(0) = ek
c      kone(1) = modk * sk * sphik
c      kone(2) = modk * sk * cphik
c      kone(3) = modk * ck
c end flat
      p4a(0)   =  Es2
      p4a(1:3) = -kone(1:3)

******
      q24a = p2a - p4a

*** now I bring q13a and q24a in the kone RF. Notice that here q13a and q24a are always back-to-back.
      call new_boost(kone,q13a,q13a,1)
      call new_boost(kone,q24a,q24a,1)
      
      sqmk2 = sqrt(mk2)
      p5(0) = sqmk2*0.5d0
      p6(0) = p5(0)
      p5m   = sqrt(mk2*0.25d0-m52)

      IF (iper.gt.-1) THEN ! i.e always
cc      IF (iper.ge.1) THEN
cccc      IF (iper.gt.10) THEN ! i.e never
cc here A = B !      
         A = dot(q13a,q13a) - 2.d0*q13a(0)*p5(0)
         A = A * 0.5d0 /sqrt(tridot(q13a,q13a))/p5m
         B = dot(q24a,q24a) - 2.d0*q24a(0)*p5(0)
         B = B * 0.5d0 /sqrt(tridot(q24a,q24a))/p5m
         call anothercollinear(abs(A),abs(B),0.5d0,q13a,q24a,lastk,
     .        wcoll)
         w = w * wcoll

*** symmetrizing p5 and p6      
         call getrnd(csi1,1)      
         if (csi1(1).le.0.5d0) then      
            p5(1:3) =  p5m * lastk(1:3)
            p6(1:3) = -p5(1:3)
            anum = 1.d0/abs(dot(q13a-p5,q13a-p5))
     .           + 1.d0/abs(dot(q24a-p5,q24a-p5))
            den  = anum +
     .             1.d0/abs(dot(q13a-p6,q13a-p6))
     .           + 1.d0/abs(dot(q24a-p6,q24a-p6))
         else
            p6(1:3) =  p5m * lastk(1:3)
            p5(1:3) = -p6(1:3)
            anum = 1.d0/abs(dot(q13a-p6,q13a-p6))
     .           + 1.d0/abs(dot(q24a-p6,q24a-p6))
            den = anum +
     .             1.d0/abs(dot(q13a-p5,q13a-p5))
     .           + 1.d0/abs(dot(q24a-p5,q24a-p5))
         endif
         w = w * 2.d0*anum/den
      ELSE
         call getrnd(csi1,1)
         an  = 2.d0
         c5 = an*csi1(1) + (-1.d0)
         s5 = sqrt(abs(1.d0-c5*c5))
         w   = w * an
         call getrnd(csi1,1)      
         phi5 = 2.d0*pi*csi1(1)
         w    = w * 2.d0*pi
         sphi5 = sin(phi5)
         cphi5 = cos(phi5)
         p5(1) = p5m * sphi5 * s5
         p5(2) = p5m * cphi5 * s5
         p5(3) = p5m * c5
         p6(1:3) = -p5(1:3)
      ENDIF
*** c5 or c6 later to be sampled according to 1/((q24a-p5[6])^2-m5[6]2) and/or 1/((q13a-p5[6])^2-m5[6]2)
c flat
c      call getrnd(csi1,1)
c      an  = 2.d0
c      c5 = an*csi1(1) + (-1.d0)
c      s5 = sqrt(abs(1.d0-c5*c5))
c      w   = w * an
c      call getrnd(csi1,1)      
c      phi5 = 2.d0*pi*csi1(1)
c      w    = w * 2.d0*pi
c      sphi5 = sin(phi5)
c      cphi5 = cos(phi5)
c      p5(1) = p5m * sphi5 * s5
c      p5(2) = p5m * cphi5 * s5
c      p5(3) = p5m * c5
c      p5(1:3) = p5m * lastk(1:3)
c      p6(1:3) = -p5(1:3)
c end flat

****
      phsp = phsp * p5m/sqmk2 * 0.5d0
****      

*** now p5 and p6 back in the q13+p2 RF, ie where kone is kone      
      call new_boost(kone,p5,p5,-1)
      call new_boost(kone,p6,p6,-1)

***   and now back in the p1+p2 RF
      if (ich.eq.1) then
         call new_boost(q13+p2,p4a,p4,-1)
         call new_boost(q13+p2,p5,p5,-1)
         call new_boost(q13+p2,p6,p6,-1)
         p3 = ps
      else
         call new_boost(q13+p1,p4a,p3,-1)
         call new_boost(q13+p1,p5,p5,-1)
         call new_boost(q13+p1,p6,p6,-1)
         p4 = ps
      endif
***
      phsp = phsp * twopis(2)*0.25d0
c      if (iid.gt.0) phsp = 0.5d0*phsp
c     i.e.
      phsp = phsp * (iid * 0.5d0 + (1-iid))
      
*************************
*** since the are two channel also for iid=0, there is a regulator
      
      f0_1 = 1.d0/dot(p1-p3,p1-p3)
      f0_2 = 1.d0/dot(p2-p4,p2-p4)

      f1_1 = 1.d0/(dot(p4+p5+p6,p4+p5+p6)-m42)
      f1_2 = 1.d0/(dot(p3+p5+p6,p3+p5+p6)-m32)

      f2_1 = 1.d0/dot(p5+p6,p5+p6)!**2
      f2_2 = f2_1
      
      a = dot(p2-p5-p6,p2-p5-p6) - m42
c      b = dot(p1-p3-p5-p6,p1-p3-p5-p6) !i.e.
      b = dot(p4-p2,p4-p2)
      f3_1 = 1.d0/a + 1.d0/b

      a = dot(p1-p5-p6,p1-p5-p6) - m32
      b = dot(p1-p3,p1-p3)
      f3_2 = 1.d0/a + 1.d0/b

      if (np.eq.1) then
         if (iid.eq.0) then
      
            den = pch1*f0_1*f1_1*f2_1*f3_1
            den = den + pch2*f0_2*f1_2*f2_1*f3_2
            
            if (ich.eq.1) then
               anum = f0_1*f1_1*f2_1*f3_1
            else
               anum = f0_2*f1_2*f2_2*f3_2
            endif
            w = w * anum/den
         else            
            
            f0_1f = 1.d0/dot(p1-p3,p1-p3)
            f0_2f = 1.d0/dot(p2-p6,p2-p6)
            
            f1_1f = 1.d0/(dot(p4+p5+p6,p4+p5+p6)-m42)
            f1_2f = 1.d0/(dot(p3+p5+p4,p3+p5+p4)-m32)
            
            f2_1f = 1.d0/dot(p5+p4,p5+p4)!**2
            f2_2f = f2_1f
            
            a = dot(p2-p5-p4,p2-p5-p4) - m42
c      b = dot(p1-p3-p5-p6,p1-p3-p5-p6) !i.e.
            b = dot(p6-p2,p6-p2)
            f3_1f = 1.d0/a + 1.d0/b
            
            a = dot(p1-p5-p4,p1-p5-p4) - m32
            b = dot(p1-p3,p1-p3)
            f3_2f = 1.d0/a + 1.d0/b

            den = pch1*f0_1*f1_1*f2_1*f3_1
            den = den + pch2*f0_2*f1_2*f2_1*f3_2

            den = den + pch1*f0_1f*f1_1f*f2_1f*f3_1f
            den = den + pch2*f0_2f*f1_2f*f2_1f*f3_2f

            den = den * 0.5d0

            iswap = 0
            call getrnd(csi1,1)
            if (csi1(1).lt.0.5d0) iswap = 1

            if (iswap.eq.1) then
               ptmp = p4
               p4   = p6
               p6   = ptmp
            endif
            
            if (ich.eq.1) then
!no! because the weight is calculated before swapping 4 6!  anum = (1-iswap) * f0_1*f1_1*f2_1*f3_1
!no!               anum = anum + iswap * f0_1f*f1_1f*f2_1f*f3_1f
               anum = f0_1*f1_1*f2_1*f3_1
            else
!no!               anum = (1-iswap) * f0_2*f1_2*f2_2*f3_2
!no!               anum = anum + iswap * f0_2f*f1_2f*f2_2f*f3_2f
               anum = f0_2*f1_2*f2_2*f3_2
            endif
            w = w * anum/den
            
         endif
      endif
******************************************************************
      if (iid.gt.0.and.np.eq.2) then
!     for muon pairs
!     here ich = 1 always         
         if (Qmu.lt.0.d0) then
            t13 = dot(p1-p3,p1-p3)
            t16 = dot(p1-p6,p1-p6)
            
            anum = 1.d0/t13     !**2               
            den  = anum + 1.d0/t16 !**2
            call getrnd(csi1,1)
            if (csi1(1).le.0.5d0) then
               ptmp = p3
               p3   = p6
               p6   = ptmp
            endif
            w = w * anum/den * 2.d0
         else
            t13 = dot(p1-p3,p1-p3)
            t15 = dot(p1-p5,p1-p5)
            anum = 1.d0/t13     !**2            
            den  = anum + 1.d0/t15 !**2
            call getrnd(csi1,1)
            if (csi1(1).le.0.5d0) then
               ptmp = p3
               p3   = p5
               p5   = ptmp
            endif
            w = w * anum/den * 2.d0
         endif
      endif
*************************         
      return
      end
*********************************************      
      subroutine anothercollinear(A,B,proba,pa,pb,vers,w)
      implicit double precision (a-h,o-z)
      dimension vers(0:3),pa(0:3),pb(0:3),xi(1)
      parameter (pi = 4.d0*atan(1.d0))

ccc      if (A.lt.1.d0.or.B.lt.1.d0) print*,A,B
      
      ana = log( (A + 1.d0) / (A - 1.d0) )
      anb = log( (B + 1.d0) / (B - 1.d0) )
      ant = ana + anb
      
      pra = ana / ant
      prb = anb / ant

      if (proba.gt.0.d0) then
! governed by proba if it's > 0
         pra = proba
         prb = 1.d0 - pra
      endif

      call getrnd(xi,1)
      if (xi(1).le.pra) then
         call getrnd(xi,1)
         c =  A - (1.d0+A) * exp(-ana*xi(1))
      else
         call getrnd(xi,1)
         c =  B - (1.d0+B) * exp(-anb*xi(1))
         c = -c ! because anyway pa and pb are always back-to-back
      endif

      call getrnd(xi,1)
      phi = 2.d0*pi * xi(1)

      s       = sqrt(abs(1.d0-c*c))
      vers(0) = 1.d0
      vers(1) = s*sin(phi)
      vers(2) = s*cos(phi)
      vers(3) = c

      call rot(-1,pa,vers,vers)
      
      ca = tridot(vers,pa)/sqrt(tridot(pa,pa))
      cb = tridot(vers,pb)/sqrt(tridot(pb,pb))

      w = 2.d0*pi / (pra/ana/(A-ca) + prb/anb/(B-cb))
      return
      end
