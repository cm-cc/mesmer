      subroutine generate_event(pmubeam,nfspart,mcids,pmat,wout,itag,
     .     ievtnr,wnovp,wnohad,wLO,wNLO,cwvpout,ierr)
      implicit double precision (a-h,o-z)
*** inputs & outputs *****************
      integer mcids(20),nfspart,itag,ierr
      integer*8 ievtnr
c      double precision pmu3beam(1:3) ! only tri-momentum!
      double precision pmubeam(0:3)
      double precision pmat(0:3,20),cwvpout(0:10)
      double precision wout,wnovp,wnohad,wLO,wNLO
**************************************
      parameter (mph=40,maxnw=3)
      double precision distwgts(maxnw)
      double precision csi(1)
      dimension p1(0:3),p2(0:3),qph(mph,0:3),qphlab(mph,0:3)
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
      integer ismualongz
      common/ISmuondirection/ismualongz
*********************************************
      double precision Qmu
      common/muoncharge/Qmu
      
      common/momentainitial/pin1,pin2

      common/tmintmax/tmin,tmax
      common/tmintmaxfordistributions/tmind,tmaxd
      common/mueexpsetup/emulab,eemin,semu,thmumin,themin,themax,
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

** from cuts      
      dimension p3lab(0:3),p4lab(0:3),p5lab(0:3),p6lab(0:3)
      common/fslabmomenta/p3lab,p4lab,p5lab,p6lab
***
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
     .     hitpmiss,istopsearch,nneg,sumover,sum2over,fmax,
     .     ntomicro,k,npoints,naccepted,kt,itime,rlxstate,nphmax,nover, 
     .     in_conf_spin
      integer ifirstmain,ifirstsub
      common/firstmainsub/ifirstmain,ifirstsub
      data ifirstmain,ifirstsub /0,0/
      common/varversion2/xsec2v2
      
****** coefficients for VP reweighting, from  invariants.h, to be put in output
      double precision cwvp(0:10),cwvpl(0:10)
      common/coeffforvpreweighting/cwvp,cwvpl
******************************************************************
      logical needswarmup,imwarmingup
      common/warmuppery/needswarmup,imwarmingup
      logical areweighted
      common/weightedeventsornot/areweighted

******* for the header      
      common/tobewrittenintheheader/wnormalization,sdmax
****
      
      if (ifirstsub.eq.0) then
****** RESETTING store IF USING MESMER AS SUBROUTINE
         if (store.eq.'yes') then
            store = 'no'
            print*,' :::'
            print*,' ::: Using MESMER as a subroutine, thus'
            print*,' ::: ''store'' variable is reset to ''no'','
            print*,' ::: because presumably events are stored'
            print*,' ::: independently...'
            print*,' :::'
            print*,' '
         endif
******
         greekmu = ucs4_'\u03BC'
         greekmu = char(int (z'03BC'), ucs4) 

****  RNG initialization. I do it here, so I redo it also after warmup
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
         
         sum        = 0.d0
         sum2       = 0.d0
         sumpairs   = 0.d0
         sumpairs2  = 0.d0
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
         xsec2v2      = 0.d0
         var          = 0.1d0
         naccepted    = 0
         nwhenmax     = 0
         nover        = 0
         hitpmiss     = 0.d0
         hit          = 0.d0
         nneg         = 0
         sumover      = 0.d0
         sum2over     = 0.d0
         sumneg       = 0.d0
         sum2neg      = 0.d0
         ie           = 0
         
         istopsearch  = 0 ! not used anymore here
cc         if (wnorm.gt.0.d0) istopsearch = 1
         
         k          = 0
         kt         = 0
         kwritelast = 0
         qph        = 0.d0
         qphlab     = 0.d0
         ng         = 0

! needs to be done outside after calling userinterface if (needswarmup) imwarmingup = .true.
         fmax = 1.1d0*sdifmax
         
         ifirstsub = 1
** [[end initialization]]
      endif
      
* [[[[[ LOOP OVER EVENTS ]]]]]      
cc      do while(kt.lt.npoints)
cc[[[[[ LOOP OVER NWARMUP !! ]]]]]
      iregenerate = 1
      do while(imwarmingup.or.iregenerate.eq.1)
** [[begin initialization]]
         if (iverbose.eq.1) call getrndstatus(irngstate)
************************
         call syncrndseq
************************
         if (imwarmingup) then
            call beamprofile(emulab,semu,emuspr,truep1lab)
         else
            truep1lab(0)   = 0.d0
            truep1lab(1:3) = pmubeam(1:3)
! making sure it's on shell
            truep1lab(0) = sqrt(ammu*ammu + tridot(truep1lab,truep1lab))            
         endif
c         emuspr         = emulab ! to be spreaded
c     MUST BE DONE OUTSIDE         call beamprofile(emulab,semu,emuspr,truep1lab)

         emuspr = truep1lab(0)

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

         k = k + 1

         call multiplicity(eps,ecms,1.d0,ng,npairs,wnphot)
         sdif = wnphot
****
         if (npairs.eq.10) then
            call phasespacemuepi0(s,m1,m2,ammu,ame,pin1,pin2,
     .           p1,p2,p3,p4,npairs,phsp,w,ie)            
         elseif (npairs.gt.0) then
            qph = 0.d0
c DEFAULT
            call phasespacemuepairsv4(s,m1,m2,ammu,ame,pin1,pin2,
     .           p1,p2,p3,p4,npairs,phsp,w,ie)            
         else
            p3 = 0.d0
            p4 = 0.d0
            call phasespacemue_v2(s,m1,m2,ammu,ame,ng,pin1,pin2,
     .           p1,p2,qph,phsp,w,ie)
         endif
***
         ieb = 0
         ie  = ie + ieb          ! ieb if beam energy spread fails
         if (ieb.gt.0) phsp = 0.d0
***
         sdif = sdif * phsp * w * watom         
*****************************************
         if (ie.ge.1) ie = 1
         if (ie.lt.1) then 
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
         
         flux = 4.d0*sqrt(dot(pin1,pin2)**2-m12*m22)
!     also this ! flux = 2.d0*sqrt(lambda(s,me*me,mm*mm)) !!

         emtx = emtx/in_conf_spin/flux ! divided by initial spin conf and divided by the flux

         if (ie.eq.0) then
c            call svfactor(model,ng,ecms,p1,p2,eps,sv,deltasv)
            sv   = 1.d0
            sdif = sdif * sv
         else
            sdif = 0.d0
         endif
         sdif = sdif * emtx

****  Converting to microbarn and standard units
         sdif = sdif * ntomicro * convfac
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
c     original xsec = sum/k
c     original var  = sqrt(abs((sum2/k-xsec*xsec)/k))
         var  = sqrt(abs((xsec2-xsec*xsec)/k))

c         xsec2v2 = xsec2v2 + (sdif - xsec)**2
c         var     = sqrt(xsec2v2/k/k) ! this seems better when var ~ 0.
**************************************************************

         abssdif = abs(sdif)
         if (abssdif.gt.sdifmax) then 
            if (nsearch.gt.0) sdifmax = abssdif ! changing only if nsearch > 0
            nwhenmax = ng
         endif
         
         if (k.ge.nsearch.and.imwarmingup) then
            print*,'Warming-up completed.'
            ifirstsub = 0 ! so next time everything is reset
            imwarmingup = .false.
            istopsearch = 1     ! not used anymore here
            iregenerate = 0
            if (hitpmiss.lt.1.d0) then
               fmax = 1.1d0*sdifmax
               if (wnorm.lt.0.d0) wnorm = xsec ! changing only if not explicitly set in input
               print*,'Starting now event generation.'
               wnormalization = wnorm
               sdmax = 1.1d0*sdifmax
               if (store.eq.'yes') then
                  call initstorage(storefile)
               endif
            endif
         endif
         
         if (.not.imwarmingup) then
            hitpmiss = hitpmiss + 1.d0
            call getrnd(csi,1)
!!     unweightening for unweighted events...
            if (fmax*csi(1).lt.abssdif) then
               uno = nint(sdif/abssdif) ! added with +-1
               hit = hit + uno
               nh = hit
               iregenerate = 0
!!---- storing in case of unweighted generation ----
               if (store.eq.'yes') then
                  if (.not.areweighted) then
                     call eventstorage(1.d0*uno,k,ie,p1,p2,p3,p4,qph)
                  endif
               endif
!!--------------------------------------------------
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
c         writecond = k.eq.npoints.or.hit.eq.npoints
         if (iwriteout.gt.0) then
            writecond = mod(k,iwriteout).eq.0.or.writecond
         else
            if (mod(k,30000).eq.0) then
               call cpu_time(tmptime)
               if ((tmptime-starttime).ge.-iwriteout) then
                  writecond = .true.
! done inside printoutput                  starttime = tmptime
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
            if (.not.imwarmingup) call writedistributions(1.d0*k,nw)
cc            call writedistributions(hit/xsec) ! for unweighted events
         endif
**[[END WRITING OUTPUT]]
         
*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
********  output momenta etc!
         mcids   = 0
         wout    = 0.d0
         pmat    = 0.d0
         wnovp   = 0.d0
         wnohad  = 0.d0
         wLO     = 0.d0
         wNLO    = 0.d0
         cwvpout = 0.d0
         ierr    = 1
         
!         pmat(0:3,1) = p1lab ! this must become truep1lab
         pmat(0:3,1) = truep1lab
         pmat(0:3,2) = p2lab
         itag   = iseedext
         ievtnr = k
         
         mcids(1) = 13
         if (Qmu.gt.0.d0) mcids(1) = -13
         mcids(2) = 11
         mcids(3) = mcids(1)
         mcids(4) = mcids(2)
         
         nfspart = 2 ! default as agreed with Giovanni
         if (ie.eq.0) then
            ierr = ie
            pmat(0:3,3) = p3lab ! from cuts
            pmat(0:3,4) = p4lab
            if (ng.gt.0) then
               do ki = 1,ng
                  call new_boost(p12lab,qph(ki,:),qphlab(ki,:),-1)
                  pmat(0:3,ki+4) = qphlab(ki,0:3)
                  mcids(ki+4)    = 22
               enddo
            endif
            nppp = 0
            if (npairs.gt.0) then
               pmat(0:3,5+ng) = p5lab
               pmat(0:3,6+ng) = p6lab
               if (npairs.eq.1) then
                  mcids(5+ng) = -mcids(2)
                  mcids(6+ng) =  mcids(2)
                  pmat(0:3,6+ng) = p6lab                  
                  nppp = 2
               elseif (npairs.eq.2) then
                  mcids(5+ng) = -mcids(1)
                  mcids(6+ng) =  mcids(1)
                  pmat(0:3,6+ng) = p6lab                  
                  nppp = 2
               elseif (npairs.eq.10) then
                  mcids(5+ng) = 22
                  mcids(6+ng) = 22
                  nppp = 2
               endif
            endif
            nfspart = 2 + ng + nppp

            iclean  = 1
            wout    = sdif/wnorm
            if (.not.areweighted) then
               wout    = sdif/abssdif
               iclean  = 0
            endif
            
            wnovp   = wout*extraweights(1) *iclean
            wnohad  = wout*extraweights(2) *iclean
            wLO     = wout*reweightLO      *iclean
            wNLO    = wout*reweightNLO     *iclean
            cwvpout = wout * cwvp          *iclean
****  rotating for beam profile
            if (ismualongz.eq.0) then
               do ii = 1,(nfspart+2)
                  call rot(-1,truep1lab,pmat(:,ii),pmat(:,ii))
               enddo
            endif
*******************************
         endif
*************************************************************************
*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      enddo
 100  continue
* [[[[[ END LOOP OVER WARMING UP EVENTS ]]]]]      
c     close(33)
      return
      end
