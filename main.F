      program main
      implicit double precision (a-h,o-z)
      parameter (mph=40,maxnw=3)
      double precision distwgts(maxnw)
      double precision csi(1)
      dimension p1(0:3),p2(0:3),qph(mph,0:3)
      dimension p3(0:3),p4(0:3)
      dimension pin1(0:3),pin2(0:3),ptmp(0:3)
      integer isvec(25),irngstate(104)
      common/rlxstatus/isvec,irngstate
*** filled in the subroutine userinterface
      character*6   ord
      character*10  model
      character*100 outfile,storefile
      character*3   eventlimiter,store
      character*2 fs
      character*20 hstnm
      common/hostandpid/hstnm,idproc
      common/finalstate/fs
      common/ecms/ecms,ecmsnom
      common/nphot_mode/nphotmode
      common/epssoft/eps
      common/parameters/ame,ammu,convfac,alpha,pi
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync,
     .     istorver,icwvp,ivpfl
      common/qedORDER/ord
      common/charinput/model,eventlimiter,store,storefile,outfile
      character*200 path
      common/cmnpath/path
      common/realinput/anpoints,sdifmax
      common/iseedoninput/iseedext
      common/ialpharunning/iarun
      common/teubner/iteubn
      common/nskvp/inskvp
*********************************************
      double precision Qmu
      common/muoncharge/Qmu
      
      common/momentainitial/pin1,pin2

      common/tmintmax/tmin,tmax
      common/tmintmaxfordistributions/tmind,tmaxd
      common/mueexpsetup/emulab,eemin,eemax,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela

      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      double precision truep1lab(0:3)
      common/labmomenta/p1lab,p2lab,truep1lab,p12lab
      
      double precision lambda,m1,m2,m12,m22
      external lambda
      common/processmasses/m1,m2,m12,m22
      
***   for reweighenting
      integer nextraweights
      parameter (nextraweights=4)
      double precision weightdefault,extraweights(0:nextraweights)
      double precision reweightLO,wnorm
      common/weights/reweightLO,reweightNLO,weightdefault,extraweights
      common/weightnormalization/wnorm
***      
      logical writecond
      real starttime,endtime,truestarttime,tmptime
      
      dimension sump(0:mph-1),sum2p(0:mph-1)
      dimension fractions(0:mph-1)
      dimension xsecp(0:mph-1),varp(0:mph-1)
      double precision sumpairs,sumpairs2,sumpairsv(2),sumpairsv2(2)      
      integer*8 npoints,k,naccepted,iwriteout,kt,itime,itime2,nh
      integer*8 kwritenow,kwritelast
      integer ucs4
      parameter (ucs4  = selected_char_KIND ('ISO_10646'))
      character(len=1,kind=ucs4) greekmu
      integer*8 foonhpm,foonh
      double precision ntomicro
      integer rlxstate(104)      
      common/forfooter/xsw,exsw,xsunw,exsunw,
     .     xsbias,exsbias,xsbiasn,exsbiasn,truemax,sum,sum2,sumneg,
     .     sum2neg,foonhpm,foonh,nabove,nlt0
      common/onceinmain/truestarttime,starttime,s,sumw1,sumw2,
     .     sumw12,xsec,var,hmxsect,hmerr,
     .     sumw22,sumLO,sumLO2,sumNLO,sumNLO2,sumpairs,sumpairs2,
     .     sumpairsv,sumpairsv2,sump,sum2p,sdiftruemax,sdiftruemin,hit,
     .     hitpmiss,
     .     istopsearch,nneg,sumover,sum2over,fmax,
     .     ntomicro,k,npoints,naccepted,kt,itime,rlxstate,nphmax,nover, 
     .     in_conf_spin
      integer ifirstmain,ifirstsub
      common/firstmainsub/ifirstmain,ifirstsub
      data ifirstmain,ifirstsub /0,0/
***
      character*500 datacard
      integer isstandalone
      common/standalone/isstandalone,datacard ! set in main.F or init_mesmer() !!
***
      integer ismualongz,itmp
      common/ISmuondirection/ismualongz
      logical needswarmup,imwarmingup
      common/warmuppery/needswarmup,imwarmingup
      logical areweighted
      common/weightedeventsornot/areweighted

**** for the header of the events      
      common/tobewrittenintheheader/wnormalization,sdmax
****
      isstandalone = 1
****************************************
      call userinterface
****************************************
      if (.not.needswarmup.and.store.eq.'yes') then
         wnormalization = wnorm
         sdmax = 1.1d0*sdifmax
         call initstorage(storefile)
      endif
****************************************      
      if (needswarmup) then
         imwarmingup = .true.      
      endif
      itmp = ismualongz

*------------------------------------      
 4242 continue

      if (imwarmingup) then
         ismualongz = 1
      else
         ismualongz = itmp
      endif
      
      if (ifirstmain.eq.0) then
         greekmu = ucs4_'\u03BC'
         greekmu = char(int (z'03BC'), ucs4) 

****  RNG initialization
         call init_rnd_sequence ! uses iseedext and sets iseed
****  end RNG initialization
      
         itime = time8()
         call cpu_time(starttime)
         truestarttime = starttime
      
         ntomicro = 0.001d0
      
** [[begin initialization]]
** [[ this is only to set distribution limits...      
         s = ame*ame + ammu*ammu + 2.d0*ame*emulab
         pcmm  = 0.5d0 * sqrt(lambda(s,ame*ame,ammu*ammu)/s)
         tmind = -4.d0*pcmm*pcmm
         tmaxd =  2.d0*ame*ame - 2.d0*Eemin*ame
** ]]                  
         ecmsnom = ecms
         npoints = anpoints

         in_conf_spin = 4

         sumw1   = 0.d0
         sumw2   = 0.d0
         sumw12  = 0.d0
         sumw22  = 0.d0
         sumLO   = 0.d0
         sumLO2  = 0.d0
         sumNLO  = 0.d0
         sumNLO2 = 0.d0
         
         sum = 0.d0
         sum2 = 0.d0
         sumpairs = 0.d0
         sumpairs2 = 0.d0
         sumpairsv  = 0.d0
         sumpairsv2 = 0.d0
         do kk = 0,mph-1
            sump(kk) = 0.d0
            sum2p(kk) = 0.d0
         enddo

         sdiftruemax = -1d30
         sdiftruemin =  1d30
      
         nphmax       = 0
         xsec         = 0.d0
         xsec2        = 0.d0
         var          = 0.1d0
         naccepted    = 0
         nwhenmax     = 0
         nover        = 0
         hitpmiss     = 0.d0
         hit          = 0.d0
         istopsearch  = 0
         nneg         = 0
         sumover      = 0.d0
         sum2over     = 0.d0
         sumneg       = 0.d0
         sum2neg      = 0.d0
         ie           = 0

         xsec2v2 = 0.d0
         
         if (wnorm.gt.0.d0) istopsearch = 1
         
         k  = 0
         kt = 0
         kwritelast = 0
         qph = 0.d0
         qphlab = 0.d0
         ng = 0
         fmax = 1.1d0*sdifmax

***   gamma & beta boost
         bboost = (emulab**2-ammu**2)/(emulab+ame)**2
         bboost = sqrt(bboost)
         gg  = 1.d0/sqrt(1.d0-bboost*bboost)
         gumb = gg*(1-bboost)
         gupb = gg*(1+bboost)
         
         
         ifirstmain = 1
** [[end initialization]]
      endif
      
* [[[[[ LOOP OVER EVENTS ]]]]]      
      iregenerate = 1
      
      do while(kt.lt.npoints.or.imwarmingup.or.iregenerate.eq.1)
         
** [[begin initialization]]
         if (iverbose.eq.1) call getrndstatus(irngstate)
************************
         call syncrndseq
************************
         emuspr = emulab        ! to be spreaded
         call beamprofile(emulab,semu,emuspr,truep1lab)
         
         m1  = ammu
         m2  = ame
         m12 = m1*m1
         m22 = m2*m2
         
         snom = m12 + m22 + 2.d0*m2*emulab
         ecmsnom = sqrt(snom)
        
         plab = sqrt(emuspr*emuspr - m12)

         p1lab(0)   = emuspr
         p1lab(1:2) = 0.d0
         p1lab(3)   = plab
         
         p2lab(0)   = m2
         p2lab(1:3) = 0.d0
         
         watom = 1.d0
c        call getpeinverysimple(p1lab,p2lab,watom)
c        call getpeinverysimple(ptmp,p2lab,watommmm)
c        call getpein(p1lab,p2lab,watom)
c        call getpein(p1lab,ptmp,watommmm)
        
         p12lab = p1lab+p2lab        
         s      = dot(p12lab,p12lab)
         
         pcmm   = 0.5d0 * sqrt(lambda(s,m12,m22)/s)
         ee     = sqrt(pcmm*pcmm + m22)
         emu    = sqrt(pcmm*pcmm + m12)
         
         ecms   = emu + ee
        
cc notice: pin[12] ALWAYS along z
         pin1(0) = emu
         pin1(1) = 0.d0
         pin1(2) = 0.d0
         pin1(3) = pcmm
         
         pin2(0) =  ee
         pin2(1) =  0.d0
         pin2(2) =  0.d0
         pin2(3) = -pcmm
******        
         tmin = -4.d0*pcmm*pcmm
         tmax =  2.d0*m22 - 2.d0*Eemin*m2
** [[end initialization]]

** analytic integrated LO xsection
c        cumt2 = 2.d0*((ame*ame+ammu*ammu)**2
c     .        + s*(s-2.d0*(ame*ame+ammu*ammu)))
c        print*,cumt2
c        print*,2.d0*(s-ame*ame-ammu*ammu)**2
c        cumt1 = 2.d0*s
c        cumt0 = 1.d0
c        aint = cumt0 * (tmax - tmin) + cumt1 * log(tmax/tmin)
c        aint = aint - cumt2 * (1.d0/tmax - 1.d0/tmin)
c        aint = aint * convfac*4.d0*pi*alpha*alpha*ntomicro
cc        aint = aint * 0.25d0/sqrt(dot(pin1,pin2)**2-ame*ame*ammu*ammu)
c        aint = aint * 0.25d0/ame/sqrt(emuspr*emuspr-ammu*ammu)        
c        aint = aint /sqrt(lambda(s,ame*ame,ammu*ammu))
c        print*,aint
c        stop       
****************************************
         k = k + 1

         call multiplicity(eps,ecms,1.d0,ng,npairs,wnphot)
         sdif = wnphot
****
         if (npairs.eq.10) then
            call phasespacemuepi0(s,m1,m2,ammu,ame,pin1,pin2,
     .           p1,p2,p3,p4,npairs,phsp,w,ie)            
         elseif (npairs.gt.0) then
            qph = 0.d0
c            call phasespacemuepairs(s,m1,m2,ammu,ame,pin1,pin2,
c     .           p1,p2,p3,p4,npairs,phsp,w,ie)
c DEFAULT
            call phasespacemuepairsv4(s,m1,m2,ammu,ame,pin1,pin2,
     .           p1,p2,p3,p4,npairs,phsp,w,ie)            
c            call phasespacemuep_phspcheck(s,m1,m2,ammu,ame,pin1,pin2,
c     .           p1,p2,p3,p4,phsp,npairs,w,ie)
c            call phasespacemuepairsCRF(s,m1,m2,ammu,ame,pin1,pin2,
c     .           p1,p2,p3,p4,npairs,phsp,w,ie)
         else
            p3 = 0.d0
            p4 = 0.d0
            call phasespacemue_v2(s,m1,m2,ammu,ame,ng,pin1,pin2,
     .           p1,p2,qph,phsp,w,ie)
c         call phasespacemue_v1(s,m1,m2,ammu,ame,ng,pin1,pin2,
c     .        p1,p2,qph,phsp,w,ie)
         endif
***
         ieb = 0
         ie  = ie + ieb ! ieb if beam energy spread fails
         if (ieb.gt.0) phsp = 0.d0
***
         sdif = sdif * phsp * w * watom
         
*****************************************
         if (ie.ge.1) ie = 1
         if (ie.lt.1) then
            call momentatolab(p1,p2,p3,p4,ng,qph)
            call cuts(p1,p2,p3,p4,qph,ng,npairs,icut)
         else
            icut = 1
         endif
         ie = ie + icut
         
         if (ng.gt.nphmax.and.ie.eq.0) nphmax = ng
         if (icut.eq.0) naccepted = naccepted + 1

! notice: this is calculated where pin[12] are along z, while cuts are applied in the rotated frame
         call squared_matrix(model,ng,npairs,ecms,p1,p2,p3,p4,pin1,pin2,
     .         qph,ie,icalc,emtx)

cphsp
ctest phspace volume
c         emtx = emtx / ntomicro * in_conf_spin*flux/convfac

         if (ie.eq.0) then
c            call svfactor(model,ng,ecms,p1,p2,eps,sv,deltasv)
            sv   = 1.d0
            sdif = sdif * sv
         else
            sdif = 0.d0
         endif
         sdif = sdif * emtx

         flux = 4.d0*sqrt(dot(pin1,pin2)**2-m12*m22)
!     also this ! flux = 2.d0*sqrt(lambda(s,me*me,mm*mm)) !!
****  Converting to microbarn and standard units and divided by initial spin conf and divided by the flux
         sdif = sdif * ntomicro * convfac/in_conf_spin/flux
****
         if (isnan(sdif)) then
            if (iverbose.gt.0) then
               call printstatus(5,k,p1,p2,qph,ng,xsec,var,varbefore,
     .              sdif,sdifmax,fmax)
            endif
         endif

         if (sdif.gt.sdiftruemax) sdiftruemax = sdif
         if (sdif.lt.sdiftruemin) sdiftruemin = sdif
**************************************************************
         sum  = sum  + sdif
         sum2 = sum2 + sdif**2

         sumw1  = sumw1 + sdif*extraweights(1)
         sumw2  = sumw2 + sdif*extraweights(2)
         sumw12 = sumw12 + (sdif*extraweights(1))**2
         sumw22 = sumw22 + (sdif*extraweights(2))**2

         sumLO   = sumLO + sdif*reweightLO
         sumLO2  = sumLO2 + (sdif*reweightLO)**2
         sumNLO  = sumNLO + sdif*reweightNLO
         sumNLO2 = sumNLO2 + (sdif*reweightNLO)**2

         if (npairs.eq.0) then
            sump(ng)  = sump(ng)  + sdif
            sum2p(ng) = sum2p(ng) + sdif*sdif
         else
            sumpairs  = sumpairs  + sdif
            sumpairs2 = sumpairs2 + sdif*sdif
            if (npairs.le.2) then
               sumpairsv(npairs)  = sumpairsv(npairs)  + sdif
               sumpairsv2(npairs) = sumpairsv2(npairs) + sdif*sdif
            endif
         endif
         varbefore = var

**** https://diego.assencio.com/?index=c34d06f4f4de2375658ed41f70177d59
         xsec  = xsec  + 1.d0/dble(k)*(sdif - xsec)
         xsec2 = xsec2 + 1.d0/dble(k)*(sdif*sdif - xsec2) ! media dei quadrati
cc         print*,(xsec-sum/k)/xsec,(sum2/k-xsec2)/xsec2,sum2/k,xsec2
         
c     original xsec = sum/k
c     original var  = sqrt(abs((sum2/k-xsec*xsec)/k))
         var  = sqrt(abs((xsec2-xsec*xsec)/dble(k)))

c         xsec2v2 = xsec2v2 + (sdif - xsec)**2
c         var     = sqrt(xsec2v2/k/k) ! this seems better when var ~ 0., but underestimates for small k
*************************************************************
         
         abssdif = abs(sdif)
         if (abssdif.gt.sdifmax) then 
            if (nsearch.gt.0) sdifmax = abssdif ! changing only if nsearch > 0
            nwhenmax = ng
         endif

!!     unweightening for unweighted events...
         if (k.ge.nsearch.and.imwarmingup) then
            print*,'Warming-up completed.'
            ifirstmain = 0 ! so next time everything is reset
            imwarmingup = .false.
            istopsearch = 1     ! not used anymore here
            iregenerate = 0
            if (hitpmiss.lt.1.d0) then
               fmax = 1.1d0*sdifmax
               if (wnorm.lt.0.d0) wnorm = xsec ! changing only if not explicitly set in input
               print*,'Starting now event generation.'
               if (store.eq.'yes') then
                  wnormalization = wnorm
                  sdmax = fmax
                  call initstorage(storefile)
               endif
            endif
            goto 4242
         endif

         if (.not.imwarmingup) then
            hitpmiss = hitpmiss + 1.d0
            call getrnd(csi,1)
            if (fmax*csi(1).lt.abssdif) then
               uno = nint(sdif/abssdif) ! added with +-1
               hit = hit +  uno
               nh = hit
               iregenerate = 0
!     !---- storing in case of unweighted generation ----
               if (store.eq.'yes') then
                  if (.not.areweighted) then
                     call eventstorage(1.d0*uno,k,ie,p1,p2,p3,p4,qph)
                  endif
               endif
!     !--------------------------------------------------
            endif
            
            if (sdif.lt.-1.d-20) then 
               nneg  = nneg  + 1
               sumneg  = sumneg + abs(sdif)
               sum2neg = sum2neg + sdif**2
            endif
            
            if (sdif.gt.fmax) then
               nover = nover + 1
               sumover  = sumover  + sdif - fmax
               sum2over = sum2over + (sdif - fmax)**2
               if (iverbose.gt.0) then
                  call printstatus(2,k,p1,p2,qph,ng,xsec,var,varbefore,
     .                 sdif,sdifmax,fmax)
               endif
            endif
         endif

********************
         if (areweighted) iregenerate = 0 ! never regenerate!
********************
****************************************************************
         ie = ie * (1-iregenerate) ! (1-irigenerate) = 1 always for weighted generation
                                   ! = 1 if hit, = 0 if miss for unweighted generation
****************************************************************
         if (ie.eq.0.and..not.imwarmingup) kt = kt + 1

         tollerate = 1.8d0 ! 1.3d0
         if (var.gt.tollerate*varbefore
     .        .and.varbefore.gt.0.d0.and.ie.eq.0) then            
            if (iverbose.gt.0) then
               call printstatus(3,k,p1,p2,qph,ng,xsec,var,varbefore,
     .              sdif,sdifmax,fmax)
            endif
         endif

         if (ie.eq.0.and..not.imwarmingup) then
            distwgts(1) = sdif
            distwgts(2) = sdif*reweightLO
            distwgts(3) = sdif*reweightNLO
            call distributions(distwgts(1:nw),p1,p2,p3,p4,qph,nw)
         endif
*
!!---- storing in case of weighted generation ----
         if (store.eq.'yes'
     .        .and.areweighted
     .        .and..not.imwarmingup) then
            call eventstorage(sdif/wnorm,k,ie,p1,p2,p3,p4,qph)
         endif
!!--------------------------------------------------         
**[[WRITING OUTPUT]]
         writecond = kt.eq.npoints.or.nint(hit).eq.npoints
         if (iwriteout.gt.0) then
            writecond = mod(k,iwriteout).eq.0.or.writecond
         else
            if (mod(k,30000).eq.0) then
               call cpu_time(tmptime)
               if ((tmptime-starttime).ge.-iwriteout) then
                  writecond = .true.
! done inside printoutput!  starttime = tmptime                  
               endif
            endif
         endif
         if (writecond) then
            hmxsect = 0.d0
            hmerr   = 0.d0
            if (hitpmiss.gt.0.d0)  then
               hmeff   = hit/hitpmiss
               hmxsect = fmax*hmeff
               hmerr   = fmax* sqrt(hmeff*(1-hmeff)/hitpmiss)
            endif
            call printoutput
            call writedistributions(1.d0*k,nw)
cc            call writedistributions(hit/xsec) ! for unweighted events
         endif
**[[END WRITING OUTPUT]]
      enddo
* [[[[[ END LOOP OVER EVENTS ]]]]]      
c      close(33)
      foonhpm  = hitpmiss
      foonh    = hit
      xsw      = xsec
      exsw     = var
      xsunw    = hmxsect
      exsunw   = hmerr
      nabove   = nover
      nlt0     = nneg
      xsbias   = sezover
      exsbias  = errsezover
      truemax  = sdiftruemax/wnorm
      xsbiasn  = sezneg
      exsbiasn = errsezneg
      if (store.eq.'yes') then
         if (isstandalone.eq.1) print*,'...finalizing storage'
         call finalizestorage(storefile)
      endif
      call exitlooplibraries
      call exitrecola
      print*,'Generation finished'
      stop
      end
