      subroutine cuts(p1,p2,p3,p4,qph,nph,npairs,icut)
! p1(0...3): final state electron four-momentum
! p2(0...3): final state positron four-momentum
! qph(0...40,0...3): four-momenta of 40 emitted photons
! icut: icut = 1 event rejected, icut = 0 event accepted
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),ptmp(0:3)
      dimension p3(0:3),p4(0:3)
      dimension p3lab(0:3),p4lab(0:3),p5lab(0:3),p6lab(0:3)
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi

      common/mueexpsetup/emulab,eemin,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela

      common/ifirstcuts/ifirst
      data ifirst /0/
      dimension ptvector(0:3)
      dimension pelabin(0:3),pmu(0:3),pe(0:3),pmulabin(0:3),pmulo(0:3)
      dimension peid(0:3),pmuid(0:3),ppid(0:3),pmatlab(0:3,4),thslab(4)
      dimension pelo(0:3)
      common/cutcommon/pelabin,pmulabin,r,s,cthminmu,thmaxmukin,
     .     thmaxekin,emulabmin,iacademic

      common/na7variable/dthna7 ! filled here to be passed in distributions!

      double precision lambda
      external lambda

      integer djb_hash
      external djb_hash
      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      double precision truep1lab(0:3)
      common/labmomenta/p1lab,p2lab,truep1lab,p12lab

      common/identicalandperipheral/iper,noper,iid

      logical lcond1,lcond2
      
      integer inside_e(3),inside_mu(3)
      double precision pinside_e(0:3),pinside_mu(0:3),pout_1(0:3)
      double precision pout_2(0:3)
      common/insidecuts/inside_e,inside_mu,pinside_e,pinside_mu,
     .     pout_1,pout_2
      integer itriggers(4)
      common/tmintmax/tmin,tmax

      common/elasticityparam/emubeam0,ame0,amu0,ame20,amu20
      common/elasticitiesdistr/dist,distf      
      
      if (ifirst.eq.0) then

         ame0 = ame
         amu0 = ammu
         ame02 = ame0*ame0
         amu02 = amu0*amu0
         
         inside_e  = 0
         inside_mu = 0

         pelabin(0) = ame
         pelabin(1:3) = 0.d0

         pmulabin(0)   = emulab
         pmulabin(1:3) = 0.d0
         pmulabin(3)   = sqrt(emulab*emulab-ammu*ammu)
         
         r = sqrt(emulab**2-ammu**2)/(emulab+ame)

         cthminmu = sqrt((1.d0-ame/ammu)*(1.d0+ame/ammu))
         thmaxmukin = acos(cthminmu)*1d3

         cthmine   = 1.d0/r*sqrt((Eemin-ame)/(Eemin+ame))
         thmaxekin = acos(cthmine)*1d3
         
         pmatlab = 0.d0
         thslab  = 0.d0

         s  = dot(pelabin+pmulabin,pelabin+pmulabin)
         px = 0.5d0*sqrt(lambda(s,ame**2,ammu**2)/s)
         ex = (s+ammu*ammu-ame*ame)/sqrt(s)*0.5d0
         pinside_mu = (/ex,0.d0,0.d0,-px/)
         call new_boost(p12lab,pinside_mu,pinside_mu,-1)
         emulabmin = pinside_mu(0)

         iacademic = 1
         if (iid.eq.1) iacademic = 0
         
         ifirst = 1
      endif
** for elasticity routine      
      emubeam0 = p1lab(0)
**
      
cphsp
c phspace volume test
c      icut = 0
c      return

      icut = 1 ! rejected by default
      
***
* I put in p[3456]lab the final state lepton momenta...      
      call new_boost(p12lab,p1,p3lab,-1) 
      call new_boost(p12lab,p2,p4lab,-1)
      if (p3(0).gt.0.d0) then
         call new_boost(p12lab,p3,p5lab,-1)
         call new_boost(p12lab,p4,p6lab,-1)
      else
         p5lab = 0.d0
         p6lab = 0.d0
      endif

      pmatlab(:,1) = p3lab
      pmatlab(:,2) = p4lab
      pmatlab(:,3) = p5lab
      pmatlab(:,4) = p6lab

      nlept = 2
      if (npairs.gt.0.and.iacademic.eq.0) nlept = nlept + 2
      call triggerdetector(nlept,pmatlab,Ethr,ththr,thslab,itriggers)

c     nocut  itriggers = 1
      if (sum(itriggers(1:nlept)).ne.2) then
         return
      endif

      
***  now checking how many pairs mimic mu/e signal
      ncoppie = 0 
      do i = 1,nlept-1
         do j = i+1,nlept
c            if (j.ne.i.and.itriggers(i)*itriggers(j).gt.0) then
            if (itriggers(i)*itriggers(j).gt.0) then

               lcond1 = thslab(i).ge.thmumin.and.thslab(i).le.thmaxmukin
     .              .and.pmatlab(0,i).ge.emulabmin
               lcond1 = lcond1.and.pmatlab(0,j).gt.eemin
     .              .and.thslab(j).ge.themin 
               
               lcond2 = thslab(j).ge.thmumin.and.thslab(j).le.thmaxmukin
     .              .and.pmatlab(0,j).ge.emulabmin
               lcond2 = lcond2.and.pmatlab(0,i).gt.eemin
     .              .and.thslab(i).ge.themin

               if (lcond1.or.lcond2) then
                  ncoppie = ncoppie + 1
               endif
            endif
         enddo
      enddo
      
      if (ncoppie.eq.1) then
         if (itriggers(2).eq.1) pinside_e = pmatlab(:,2)
         if (itriggers(3).eq.1) pinside_e = pmatlab(:,3)
         if (itriggers(4).eq.1) pinside_e = pmatlab(:,4)
         icut = 0
         if (ina7.eq.1) then 
            aco = acoplanarityPSIdef(pmatlab(:,1),pinside_e)*1d3
            if (aco.gt.dthna7max) icut = 1
         endif         
      endif

c no cut 
c      pinside_e = pmatlab(:,2)
c      icut = 0

      dist  = 0.d0
      if (icut.eq.0.and.(nph+npairs).gt.0.and.iela.gt.0) then
         thein = pinside_e(3)/sqrt(tridot(pinside_e,pinside_e))
         thein = acos(thein)*1d3
         thmin = p3lab(3)/sqrt(tridot(p3lab,p3lab))
         thmin = acos(thmin)*1d3

c         call elasticity_distance(thein,thmin,theoutf,thmoutf)
c         dist = 1d3*sqrt((thein-theoutf)**2+(thmin-thmoutf)**2)
c         distf = dist
         
         call myelasticityopt(thein,thmin,nph+npairs,ththr,
     .        cutela,theout,thmout,dist)
         if (dist.gt.cutela) icut = 1

      endif

      return
      end
*******************************************************************
      subroutine triggerdetector(n,p,eth,thth,th,itriggers)
      implicit double precision (a-h,o-z)
      double precision p(0:3,4),th(4),eth,thth
      integer itriggers(4)
      integer k
      itriggers = 0
      do k = 1,n
         th(k) = acos(p(3,k)/sqrt(tridot(p(:,k),p(:,k))))*1d3
         if (p(0,k).ge.eth.and.th(k).le.thth) itriggers(k) = 1
      enddo
      end
******************************************************************
      subroutine beamprofile(e,s,es,truep1lab)
      implicit double precision (a-h,o-z)
      double precision truep1lab(0:3)
      character*150 line
      character*200 path
      common/cmnpath/path      
      common/parameters/ame,ammu,convfac,alpha,pi
      integer ismualongz,ifirst
      common/ISmuondirection/ismualongz
      common/bmprfl/ic,ifirst
      data ic,ifirst /0,0/

      es = e
      if (s.gt.1d-8) then
         es  = gau_spread(e,s)      
      endif

      if (ismualongz.eq.1) then
** Preparing for generic p1lab, with xyz components
         truep1lab(0) = es
         truep1lab(1) = 0.d0
         truep1lab(2) = 0.d0
         truep1lab(3) = sqrt(es*es-ammu*ammu)
      else
         if (ifirst.eq.0) then
            print*,'Creating fifo for feeding in muon beam momenta...'
            call EXECUTE_COMMAND_LINE('mkfifo '
     >           //trim(path)//'beamprofile.fifo 2>/dev/null')

            print*,'Waiting for muon beam momenta from '//
     >           trim(path)//'beamprofile.fifo'
            open(82,file=trim(path)//'beamprofile.fifo',status='old')
c            read(82,'(A)')line ! not needed, dealt with below with the err directive in read statement
c            print*,'Starting!'
            ifirst = 1
         endif
         
 999     continue
c ! this was for the file with E,posx,posy,px,py,pz         
c         read(82,*,err=999)energyinfile,posx,posy,
c     .        truep1lab(1),truep1lab(2),truep1lab(3)
c ! this is for the file with px,py,pz         
         read(82,*,err=999)truep1lab(1),truep1lab(2),truep1lab(3)
c         ic = ic + 1
c         print*,ic,truep1lab(3)
         
! forcing the muon on-shell again
         truep1lab(0)=sqrt(tridot(truep1lab,truep1lab) + ammu*ammu)

         es = truep1lab(0)
c         p = sqrt(es*es-ammu*ammu)
c         truep1lab(0) = es
c         truep1lab(1) = p * sin(45.d0*pi/180.d0)
c         truep1lab(2) = 0.d0
c         truep1lab(3) = p * cos(45.d0*pi/180.d0)
      endif
      return
      end
*****************************************************************      
      subroutine cuts_old(p1,p2,p3,p4,qph,nph,npairs,icut)
! p1(0...3): final state electron four-momentum
! p2(0...3): final state positron four-momentum
! qph(0...40,0...3): four-momenta of 40 emitted photons
! icut: icut = 1 event rejected, icut = 0 event accepted
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),ptmp(0:3)
      dimension p3(0:3),p4(0:3)
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi
      common/mueexpsetup/emulab,eemin,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela
      common/ifirstcuts/ifirst
      data ifirst /0/
      dimension ptvector(0:3)
      dimension pelabin(0:3),pmu(0:3),pe(0:3),pmulabin(0:3),pmulo(0:3)
      dimension peid(0:3),pmuid(0:3),ppid(0:3)
      dimension pelo(0:3)
      common/cutcommon/pelabin,pmulabin,r,s,cthminmu,thmaxmukin,
     .     thmaxekin,emulabmin,iacademic

      common/na7variable/dthna7 ! filled here to be passed in distributions!

      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      double precision truep1lab(0:3)
      common/labmomenta/p1lab,p2lab,truep1lab,p12lab

      common/identicalandperipheral/iper,noper,iid

      integer inside_e(3),inside_mu(3)
      double precision pinside_e(0:3),pinside_mu(0:3),pout_1(0:3)
      double precision pout_2(0:3)
      common/insidecuts/inside_e,inside_mu,pinside_e,pinside_mu,
     .     pout_1,pout_2
            
      if (ifirst.eq.0) then

         inside_e  = 0
         inside_mu = 0

!         ina7 and dthna7max read in userinterface
         
         pelabin(0) = ame
         pelabin(1:3) = 0.d0

         pmulabin(0)   = emulab
         pmulabin(1:3) = 0.d0
         pmulabin(3)   = sqrt(emulab*emulab-ammu*ammu)
         
         r = sqrt(emulab**2-ammu**2)/(emulab+ame)

         s = dot(pmulabin+pelabin,pmulabin+pelabin)


         pmui = sqrt(emulab**2-ammu*ammu)
         emu  = emulab + ame - eemin

         csup = emulab*emu-ammu**2 + ame*(ame-eemin)
         csup = csup /pmui/sqrt(emu**2-ammu**2)

         cthminmu   = sqrt((1.d0-ame/ammu)*(1.d0+ame/ammu))
         thmaxmukin = acos(cthminmu)*1d3
         cthmine    = 1.d0/r*sqrt((Eemin-ame)/(Eemin+ame))
         thmaxekin  = acos(cthmine)*1d3


******test
** from doc/Em_vs_thm.pdf         
c         thtilm = 0.5d0 * 0.001d0
c         c = cos(thtilm)
c         al = ammu*ammu + Emulab*ame
c         be = Emulab + ame
c         ga = sqrt(Emulab*Emulab - ammu*ammu)*c
c         disc = (ga*ga - be*be)*ammu*ammu + al*al
c         Ep = al*be + abs(ga)*sqrt(disc)
c         Ep = - Ep / (ga*ga - be*be)
c         Em = al*be - abs(ga)*sqrt(disc)
c         Em = - Em / (ga*ga - be*be)
c         ptmp(0) = Em
c         ptmp(1) = sqrt(Em*Em - ammu*ammu)*sqrt(1.d0-c*c)
c         ptmp(2) = 0.d0
c         ptmp(3) = sqrt(Em*Em - ammu*ammu)*c
c         print*,dot(pmulabin-ptmp,pmulabin-ptmp)        
c         stop
         
c         cth = 1.d0 + tmax*0.5d0/px/px
c         sth = sqrt(1.d0-cth*cth)         
c         ptmp(0) = ex
c         ptmp(1) = px * sth
c         ptmp(2) = 0.d0
c         ptmp(3) = px * cth
c         call new_boost(p12lab,ptmp,ptmp,-1)
c         th = ptmp(3)/sqrt(tridot(ptmp,ptmp))
c         th = acos(th)*1d3 ! this is the muon lab theta below which LO is tiny. it corresponds to set symmetric
c                           ! theta cut in the cm, ie the angle is between x and pi - x
c         print*,th,tmin,tmax         
c         stop
c$$$         thmin = 0.d0
c$$$         thmax = pi
c$$$         dth = (thmax-thmin)/100000d0
c$$$         open(888)
c$$$         do k = 0,100000
c$$$            th = thmin + k * dth
c$$$            cth = cos(th)
c$$$            sth = sqrt(1.d0-cth*cth)
c$$$            
c$$$            ptmp(0) = ex
c$$$            ptmp(1) = px * sth
c$$$            ptmp(2) = 0.d0
c$$$            ptmp(3) = px * cth
c$$$            call new_boost(p12lab,ptmp,pinside_mu,-1)
c$$$
c$$$            cth = pinside_mu(3)/sqrt(tridot(pinside_mu,pinside_mu))
c$$$            thl = acos(cth)
c$$$
c$$$            write(888,*)th,thl*1d3,pinside_mu(0),pinside_mu(3)
c$$$
c$$$         enddo
c$$$         close(888)
c$$$         stop
***** test

         
         ifirst = 1
      endif

cphsp      
c phspace volume test
c      icut = 0
c      return

      inside_e  = 0
      inside_mu = 0
      
      icut = 1 ! event rejected      

c      icut = 0
c      return
      
      dthna7 = 0.d0
****************************************************+
      
cc this doesnt work if not boosting where pin2 is at rest (p2lab is not at rest in general)
c      call new_boost(pin2,p1,pmu,1) 
c      call new_boost(pin2,p2,pe, 1)
cc this is more general
      call new_boost(p12lab,p1,pmu,-1) 
      call new_boost(p12lab,p2,pe, -1)
      thmu = acos(pmu(3)/sqrt(tridot(pmu,pmu)))*1d3
      the  = acos(pe(3) /sqrt(tridot(pe,pe)))  *1d3      

c      if (iid.eq.1.and.npairs.gt.0) then
      if (npairs.gt.0) then         
         if (npairs.eq.1) then
            call new_boost(p12lab,p4,peid, -1)
            call new_boost(p12lab,p3,ppid, -1)
            theid  = acos(peid(3) /sqrt(tridot(peid,peid)))*1d3
            thpid  = acos(ppid(3) /sqrt(tridot(ppid,ppid)))*1d3
            pmuid  = 0.d0
            thmuid = pi*0.5d0*1d3 * (1.1d0)            
         elseif (npairs.eq.2) then
            call new_boost(p12lab,p3,pmuid, -1)
            thmuid = acos(pmuid(3) /sqrt(tridot(pmuid,pmuid)))*1d3
            pmuid  = 0.d0
            thmuid = pi*0.5d0*1d3 * (1.1d0)
         endif
      else
         peid   = 0.d0
         ppid   = 0.d0
         pmuid  = 0.d0
         theid  = pi*0.5d0*1d3 * (1.1d0)
         thpid  = pi*0.5d0*1d3 * (1.1d0)         
         thmuid = pi*0.5d0*1d3 * (1.1d0)
      endif

! muon cut      
      if (thmu.ge.thmumax.and.thmuid.ge.thmumax) then
         return 
      endif      
      if (thmu.lt.thmumin) then
         return 
      endif
      
! electron cut
c      if (the .ge.100d0.and.theid.gt.100d0) return
c      if (pe(0).lt.eemin.and.peid(0).lt.eemin) return

      if (npairs.eq.1) then
         if (the  .lt.themax.and.pe(0)  .gt.eemin) inside_e(1) = 1
c         if (theid.lt.themax.and.peid(0).gt.eemin) inside_e(2) = 1
c         if (thpid.lt.themax.and.ppid(0).gt.eemin) inside_e(3) = 1

         if (sum(inside_e).eq.1) then
            continue
         else
            return
         endif
         pinside_e=inside_e(1)*pe+inside_e(2)*peid+inside_e(3)*ppid
         
      else
c         if (the.lt.themax.and.pe(0).gt.eemin) then
         if (the.lt.themax.and.the.lt.31.851481694841610d0) then
c     .        and.pe(0).gt.Eemin) then
            continue
         else
            return
         endif
      endif

      
      if (ina7.eq.1) then
* else if ina7 = 2 acoplanarity cut. the max acoplanarity is still dthna7max      
c         aco  = acoplanarityPSIdef(pmu,pe)*1d3
c     if (aco.gt.dthna7max) return
         if (npairs.eq.0) then
***   acoplanarity cut
            aco  = acoplanarity(pmu,pe)
            aco1 = abs(aco)*1d3
            aco2 = abs(2.d0*pi-aco)*1d3      
            if (aco1.gt.dthna7max.and.aco2.gt.dthna7max) return
         else
            aco  = acoplanarityPSIdef(pmu,pinside_e)*1d3
            if (aco.gt.dthna7max) return
         endif
      endif
      
      icut = 0                  ! event accepted

      return
      end
*******************************************************************
      subroutine getna7angles(emubeam,m1,m2,p1cm,p2cm,p1,p2,th1,th2,dth)
      implicit double precision (a-h,m,o-z)
      double precision p1(0:3),p2(0:3),csi(1)
      double precision p1lo(0:3),p2lo(0:3)
      double precision p1cm(0:3),p2cm(0:3)
      double precision ptmp1(0:3),ptmp2(0:3)
      common/na7vars/r,xi,al,pmui,Ei,E1th1max,c1min,th1max,E2th1max,
     .     c2atmax,th2atmax,pi,ifirst
      data ifirst /0/
      double precision pin1(0:3),pin2(0:3)
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      double precision truep1lab(0:3)
      common/labmomenta/p1lab,p2lab,truep1lab,p12lab
      common/momentainitial/pin1,pin2

      if (ifirst.eq.0) then
         pi = 4.d0*atan(1.d0)
         
         Ei   = emubeam + m2
         pmui = sqrt(emubeam*emubeam - m1*m1)

         r2 = emubeam**2-m1*m1
         r2 = r2 / (emubeam+m2)/ (emubeam+m2)
         r  = sqrt(r2)

         xi = Ei/pmui
         al = (emubeam*m2+m1*m1)/pmui

! this is the muon lab energy where thm has a maximum as a function of the         
         E1th1max = xi*m1*m1/al
! this is the cosine of the thmmax angle
         c1min    = (xi*E1th1max-al)/sqrt(E1th1max*E1th1max-m1*m1)
!     i.e., verificato analiticamente!!!
c         c1min = sqrt(1.d0-m2*m2/m1/m1)

** ecco perche' tenchini dice che la posizione del massimo non si sposta e dipende
** solo dal rapporto delle masse
**** thmmax ~ m2/m1 !!!!         
         
!     this is the thmmax angle
         th1max   = acos(c1min)
! this is the electron energy corresponding to thm = th1max
         E2th1max = Ei - E1th1max
! this is the electron cos corresponding to the the max in the the,thm plane        
         c2atmax = 1.d0/r*sqrt((E2th1max-m2)/(E2th1max+m2))
!     and this is the angle
         th2atmax = acos(c2atmax)

c         print*,th1max*1d3,th2atmax*1.d3,m2/m1*1d3
c         print*,E1th1max,E2th1max
c         stop

***   NOTA: E2th1max puo' essere usata come discriminante: al LO nel piano the,thmu
***   quando the < th2atmax si ha E_e > E2th1max mentre
***      per the > th2atmax si ha E_e < E2th1max

         ifirst = 1
      endif

      dth  = 0.d0
      
      pm2  = sqrt(tridot(p2,p2))
      ce   = p2(3)/pm2
      pm1  = sqrt(tridot(p1,p1))
      cmu  = p1(3)/pm1         
      the  = acos(ce)
      thmu = acos(cmu)

      goto 555
      
      call getrnd(csi,1)
      if (csi(1).lt.0.5d0) then
caso 1: ricostruisco p1lo a partire da p2           
         Eelo = m2 * (1.d0+r*r*ce*ce)/(1.d0-r*r*ce*ce)
         p2lo(0)   = Eelo
         p2lo(1:3) = p2(1:3)/pm2 * sqrt(Eelo*Eelo - m2*m2)
         p1lo = p1lab+p2lab-p2lo         
      else
caso 2: ricostruisco p2lo a partire da p1: occhio che si controlla comunque p2(0) o the
         disc  = xi*xi*al*al - (xi*xi-cmu*cmu)*(al*al+m1*m1*cmu*cmu)
         emuf1 = xi*al+sqrt(disc)
         emuf1 = emuf1/(xi*xi-cmu*cmu)
         emuf2 = xi*al-sqrt(disc)
         emuf2 = emuf2/(xi*xi-cmu*cmu)
         
c        if (p2(0).lt.E2th1max) then ! or
         if (the.gt.th2atmax) then
            emuf = emuf1
         else
            emuf = emuf2
         endif
         
         p1lo(0)   = emuf
         p1lo(1:3) = p1(1:3)/pm1*sqrt(emuf*emuf-m1*m1)
         
         p2lo = p1lab+p2lab-p1lo
      endif      
c      if (p2(0).gt.E2th1max) then
c         print*,' '
c         print*,p2
c         print*,p2lo
c      endif

      cmulo = p1lo(3)/sqrt(tridot(p1lo,p1lo))
      celo  = p2lo(3)/sqrt(tridot(p2lo,p2lo))
      thmulo = acos(cmulo)
      thelo  = acos(celo)

      th1 = thmulo
      th2 = thelo
      
      dth=sqrt((the-thelo)**2+(thmu-thmulo)**2)

      if (thmu.gt.thmulo) dth = -dth

***** I try to do everything in the CM
 555  continue
      
      pm2cm = sqrt(tridot(p2cm,p2cm))
      pm1cm = sqrt(tridot(p1cm,p1cm))
      cthecm = p2cm(3)/pm2cm
      cthmcm = p1cm(3)/pm1cm
      thecm = acos(cthecm)
      thmcm = acos(cthmcm)

      cphie = p2cm(1)/pm2cm/sqrt(1.d0-cthecm*cthecm)
      sphie = p2cm(2)/pm2cm/sqrt(1.d0-cthecm*cthecm)
      
      thelo = 0.5d0*(pi + thecm - thmcm)
      thmlo = pi - thelo

      cthelo = cos(thelo)
      sthelo = sqrt(1.d0-cthelo*cthelo)
      be   = sqrt(tridot(pin2,pin2))/pin2(0)
      p2lo(0) = pin2(0)
      p2lo(1) = pin2(0)*be*cphie*sthelo
      p2lo(2) = pin2(0)*be*sphie*sthelo
      p2lo(3) = pin2(0)*be*cthelo
      
      p1lo = pin1 + pin2 - p2lo

      call new_boost(pin2,p1lo,p1lo,1) 
      call new_boost(pin2,p2lo,p2lo,1)
      
      th1lo = p1lo(3)/sqrt(tridot(p1lo,p1lo))
      th2lo = p2lo(3)/sqrt(tridot(p2lo,p2lo))
      th1lo = acos(th1lo)
      th2lo = acos(th2lo)

      dth = sqrt((thmu-th1lo)**2+(the-th2lo)**2)

**************TEST
c$$$      celo  = cos(th2lo)
c$$$      cphie = p2lo(1)/sqrt(tridot(p2lo,p2lo))/sqrt(1.d0-celo*celo)
c$$$      sphie = p2lo(2)/sqrt(tridot(p2lo,p2lo))/sqrt(1.d0-celo*celo)
c$$$      thalo = th2lo*0.999d0
c$$$      cthelo = cos(thalo)
c$$$      sthelo = sqrt(1.d0-cthelo*cthelo)
c$$$      be   = sqrt(tridot(p2lo,p2lo))/p2lo(0)
c$$$      
c$$$      p2lo(1) = p2lo(0)*be*cphie*sthelo
c$$$      p2lo(2) = p2lo(0)*be*sphie*sthelo
c$$$      p2lo(3) = p2lo(0)*be*cthelo
c$$$
c$$$      
c$$$      p1lo = p1lab + p2lab - p2lo
c$$$
c$$$      tha = acos(p2lo(3)/sqrt(tridot(p2lo,p2lo)))
c$$$      thb = acos(p1lo(3)/sqrt(tridot(p1lo,p1lo)))
c$$$
c$$$      dth2 = sqrt((thmu-thb)**2+(the-tha)**2)
c$$$
c$$$      print*,dth2/dth
c$$$      
c$$$
      
c      print*,dth
c      dth = 0.d0
      
      return
      end
********************************************************************
      subroutine na7cut(p3lab,p4lab,icut)
      implicit double precision (a-h,m,o-z)
      double precision p3(0:3),p4(0:3),csi(1)
      double precision p3lab(0:3),p4lab(0:3)
      double precision p3lo(0:3),p4lo(0:3)

      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      double precision truep1lab(0:3)
      common/labmomenta/p1lab,p2lab,truep1lab,p12lab
      
      common/processmasses/m1,m2,m12,m22
      common/mueexpsetup/emulab,eemin,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela
      common/rrrre/r,radtodeg,ifirst
      dimension pin1(0:3),pin2(0:3)
      common/momentainitial/pin1,pin2
      data ifirst /0/
      if (ifirst.eq.0) then
         
         radtodeg = 1d3
         r2 = emulab**2-m12
         r2 = r2 / (emulab+m2)/ (emulab+m2)
         
         r = sqrt(r2)
         ifirst = 1
      endif
      
      icut = 1

      p3mlab = sqrt(tridot(p3lab,p3lab))      
      p4mlab = sqrt(tridot(p4lab,p4lab))

      fun1 = acos(p4lab(3)/p4mlab)*radtodeg
      fun3 = acos(p3lab(3)/p3mlab)*radtodeg

*** NA7 variable D\theta
      ce = cos(fun1/radtodeg)
      Eelo = m2 * (1.d0+r*r*ce*ce)/(1.d0-r*r*ce*ce)
      
      p4lo(0) = Eelo
      pelo = sqrt(Eelo**2 - m22)
      p4lo(1:3) = p4lab(1:3)/p4mlab * pelo
      
      p3lo = p1lab+p2lab-p4lo


c$$$      thLlo = 0.d0
c$$$      thRlo = 0.d0      
c$$$      phi3lo = atan2(p3lo(2),p3lo(1))
c$$$      phi4lo = atan2(p4lo(2),p4lo(1))      
c$$$      if (phi3lo.ge.0.d0.and.phi4lo.le.0.d0) then
c$$$         thLlo = acos(p4lo(3)/sqrt(tridot(p4lo,p4lo)))*radtodeg
c$$$         thRlo = acos(p3lo(3)/sqrt(tridot(p3lo,p3lo)))*radtodeg
c$$$      elseif (phi3lo.le.0.d0.and.phi4lo.ge.0.d0) then
c$$$         thRlo = acos(p4lo(3)/sqrt(tridot(p4lo,p4lo)))*radtodeg
c$$$         thLlo = acos(p3lo(3)/sqrt(tridot(p3lo,p3lo)))*radtodeg
c$$$      endif
c$$$      phi3 = atan2(p3lab(2),p3lab(1))
c$$$      phi4 = atan2(p4lab(2),p4lab(1))
c$$$      thL = 10000.d0
c$$$      thR = 10000.d0
c$$$      if (phi3.ge.0.d0.and.phi4.le.0.d0) then
c$$$         thL = fun1
c$$$         thR = fun3
c$$$      elseif (phi3.le.0.d0.and.phi4.ge.0.d0) then
c$$$         thR = fun1
c$$$         thL = fun3
c$$$      endif

**** ramdomized      
      call getrnd(csi,1)
      if (csi(1).lt.0.5d0) then
         thLlo = acos(p4lo(3)/sqrt(tridot(p4lo,p4lo)))*radtodeg
         thRlo = acos(p3lo(3)/sqrt(tridot(p3lo,p3lo)))*radtodeg
         thL   = fun1
         thR   = fun3
      else
         thRlo = acos(p4lo(3)/sqrt(tridot(p4lo,p4lo)))*radtodeg
         thLlo = acos(p3lo(3)/sqrt(tridot(p3lo,p3lo)))*radtodeg
         thR   = fun1
         thL   = fun3
      endif

      cut = sqrt((thR-thRlo)**2+(thL-thLlo)**2)

      if (cut.gt.dthna7max) return
      
      icut = 0
      
      return
      end
************************************************************************
      subroutine get_dangleUM(E,da)
      implicit double precision (a-h,m,o-z)
      common/gdaUM/Emu,dT,X0T,X0A,c0,c1,ifirst
      data ifirst /0/
c  3 const double Emu    = 150;
c  4 const double me  = 0.5109989461 *1e-3; // GeV
c  5 const double mmu = 105.6583745 * 1e-3; // GeV
c  6 const double dT  = 2.; // cm
c  7 const double X0T = 35.28; // cm
c  8 const double X0A = 9.37; // cm
c  9 const double c0  = 1.2 * 1e-3;   // rad * GeV : MSC at 1 GeV         
c 10 const double c1  = 0.03 * 1e-3;  // rad 
c 11 const double r = sqrt(Emu*Emu - mmu*mmu)/(Emu + me);
      if (ifirst.eq.0) then
         Emu = 150d0
         dT  = 2.d0
         X0T = 35.28d0
         X0A = 9.37d0
         c0  = 1.2d-3
         c1  = 0.03d-3
         ifirst = 1
      endif
      sT = 13.6d-3/E*sqrt(0.5d0*dT/X0T)*(1.d0+0.038d0*log(0.5d0*dT/X0T))
      sA = sqrt( c1*c1 + c0*c0/E/E)
      da = 5.d0*sqrt(sT * sT + sA * sA)
c      sigma   = sqrt(sT * sT + sA * sA)
c      da      = 5.d0 * sigma
      return
      end

***********************
      function pseudoacoll(emub,mm,me,pml,pel)
      implicit double precision (a-h,m,o-z)
      dimension pml(0:3),pel(0:3)
      dimension pmlo(0:3),pelo(0:3)
      dimension ptmp(0:3)

      a = mm*mm+emub*me
      b = emub+me

      pimu = sqrt(emub*emub - mm*mm)

      ptmp(0)   = emub+me
      ptmp(1:2) = 0.d0
      ptmp(3)   = pimu
      
      r = pimu/b
      
      pelm = sqrt(tridot(pel,pel))
      pmlm = sqrt(tridot(pml,pml))
      
      cthe = pel(3)/pelm
      cthm = pml(3)/pmlm

      g = cthm * pimu
      
      the = acos(cthe)
      thm = acos(cthm)

      thebar = acos(b/sqrt(b*b+mm*mm-me*me)) ! eq. 7 notes
      
      Eelo = me * (1.d0+r*r*cthe*cthe)/(1.d0-r*r*cthe*cthe)

      pelo    = pel/pelm * sqrt(Eelo*Eelo-me*me)
      pelo(0) = Eelo

      b2mg2 = b*b - g*g
      sq = sqrt(a*a - b2mg2*mm*mm)
      if (the.le.thebar) then
         emulo = (a*b - abs(g)*sq)/b2mg2
      else
         emulo = (a*b + abs(g)*sq)/b2mg2
      endif
      pmlo    = pml/pmlm * sqrt(emulo*emulo-mm*mm)
      pmlo(0) = emulo

      call new_boost(ptmp,pelo,pelo,1) 
      call new_boost(ptmp,pmlo,pmlo,1)
            
      pseudoacoll = acollinearity(pelo,pmlo)

      return
      end
***********************
