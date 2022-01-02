      program main
      implicit double precision (a-h,o-z)
      parameter (mph=40,maxnw=3)
      double precision distwgts(maxnw)
      double precision csi(1)
      dimension p1(0:3),p2(0:3),qph(mph,0:3)
      dimension p3(0:3),p4(0:3)
      dimension pin1(0:3),pin2(0:3),ptmp(0:3)
      integer isvec(25)
      common/rlxstatus/isvec
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
      common/mueexpsetup/emulab,eemin,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela

      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      double precision truep1lab(0:3)
      common/labmomenta/p1lab,p2lab,truep1lab,p12lab
      
      double precision lambda,m1,m2,m12,m22
      external lambda

c      double precision m1,m2,m12,m22
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
      common/oldmain/greekmu,truestarttime,s,sumw1,sumw2,sumw12,sumw22,
     .     sumLO,sumLO2,sumNLO,sumNLO2,sumpairs,sumpairs2,sumpairsv,
     .     sumpairsv2,sump,sum2p,sdiftruemax,sdiftruemin,hit,hitpmiss,
     .     istopsearch,nneg,sumover,sum2over,fmax,wnpoints,
     .     ntomicro,k,npoints,naccepted,kt,itime,rlxstate,nphmax,nover, 
     .     in_conf_spin
****************************************      
      call userinterface
****************************************      

      greekmu = ucs4_'\u03BC'
      greekmu = char (int (z'03BC'), ucs4) 
      
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
      ieb          = 0
      ie           = 0

      if (wnorm.gt.0.d0) istopsearch = 1
      
      k  = 0
      kt = 0
      kwritelast = 0
      wnpoints = npoints
      do j = 0,3
         do ki = 1,mph
            qph(ki,j) = 0.d0
         enddo
      enddo
      ng = 0
            
      s = ammu**2 + ame**2 + 2.d0*ame*emulab
      call initmuescattering(s)

***   gamma & beta boost
      bboost = (emulab**2-ammu**2)/(emulab+ame)**2
      bboost = sqrt(bboost)
      gg  = 1.d0/sqrt(1.d0-bboost*bboost)
      print*,'*********************************************************'
      print*,'CM-lab boost parameters '
      print*,'gamma, beta of the boost',gg,bboost
      gumb = gg*(1-bboost)
      gupb = gg*(1+bboost)
      print*,'gamma*(1-beta), gamma*(1+beta)',gumb,gupb
      print*,'*********************************************************'
****      
      dm = abs(ammu - ame)
      dm = dm/ammu
      if (dm.lt.1d-3) then
         print*,'mmu = me !'
         print*,'mmu = me !'
      endif
***
      if (eventlimiter.eq.'unw') wnpoints = 1.d13

* [[[[[ LOOP OVER EVENTS ]]]]]      
c      do while(k.lt.wnpoints)
      do while(kt.lt.wnpoints)
************************
        call syncrndseq
************************         
        emuspr = emulab     ! to be spreaded
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

         k = k + 1

         call multiplicity(eps,ecms,1.d0,ng,npairs,wnphot)
         sdif = wnphot
****
         if (npairs.gt.0) then
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
         ie = ie + ieb          ! ieb if beam energy spread fails
         if (ieb.gt.0) phsp = 0.d0
***         
         sdif = sdif * phsp * w * watom
         
         if (isnan(phsp).or.isnan(w)) then
            print*,'NAN found (pos 1)!!!'
            print*,isvec
            print*,'end ----!!!'
            ie = 1
         endif
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

cphsp
ctest phspace volume
c         emtx = emtx * 1d12 / ntomicro * in_conf_spin*flux/convfac

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
         iii = 0
         
         if (sdif.gt.sdiftruemax) sdiftruemax = sdif
         if (sdif.lt.sdiftruemin) sdiftruemin = sdif

         if (sdif.gt.sdifmax) then 
            if (nsearch.gt.0) sdifmax = sdif ! changing only if nsearch > 0
            nwhenmax = ng
            iii = 1
         endif
         if (istopsearch.eq.0) then            
            if (iverbose.gt.0.and.iii.eq.1) then
               call printstatus(1,k,p1,p2,qph,ng,xsec,var,varbefore,
     .              sdif,sdifmax,fmax)
            endif
         endif         
         
         if (isnan(sdif)) then
            if (iverbose.gt.0) then
               call printstatus(5,k,p1,p2,qph,ng,xsec,var,varbefore,
     .              sdif,sdifmax,fmax)
            endif
         endif
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
            sumpairsv(npairs)  = sumpairsv(npairs)  + sdif
            sumpairsv2(npairs) = sumpairsv2(npairs) + sdif*sdif
         endif
         varbefore = var
         xsec = sum/k

**** https://diego.assencio.com/?index=c34d06f4f4de2375658ed41f70177d59
         xsec  = xsec  + 1.d0/dble(k)*(sdif - xsec)
         xsec2 = xsec2 + 1.d0/dble(k)*(sdif*sdif - xsec2) ! media dei quadrati
         
c original var  = sqrt(abs((sum2/k-xsec*xsec)/k))
         var  = sqrt(abs((xsec2-xsec*xsec)/k))         
**************************************************************
         
!!     unweightening for unweighted events...
         if (k.gt.nsearch) then
            istopsearch = 1
            if (hitpmiss.lt.1.d0) then
               fmax = 1.05d0*sdifmax
               if (wnorm.lt.0.d0) wnorm = xsec ! changing only if not explicitly set in input
               print*,'Starting now also unweighted generation!'
               if (store.eq.'yes') then
                  call initstorage(storefile,wnorm,fmax)
               endif
            endif
            hitpmiss = hitpmiss + 1.d0
            call getrnd(csi,1)
            if (fmax*csi(1).lt.sdif) then
               hit = hit + 1.d0
               nh = hit
!!---- storing in case of unweighted generation ----
               if (eventlimiter.eq.'unw') then
                  if (store.eq.'yes') then
                   call eventstorage(1.d0,k,ie,p1,p2,p3,p4,qph)
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

         if (ie.eq.0.and.istopsearch.eq.1) kt = kt + 1
         
         tollerate = 1.8d0            
cc         tollerate = 1.3d0
         if (var.gt.tollerate*varbefore
     .        .and.varbefore.gt.0.d0.and.ie.eq.0) then            
            if (iverbose.gt.0) then
               call printstatus(3,k,p1,p2,qph,ng,xsec,var,varbefore,
     .              sdif,sdifmax,fmax)
            endif
         endif

         if (ie.eq.0) then
c     if (ie.eq.0.and.eventlimiter.eq.'w') then
            distwgts(1) = sdif
            distwgts(2) = sdif*reweightLO
            distwgts(3) = sdif*reweightNLO
          call distributions(distwgts(1:nw),p1,p2,p3,p4,qph,nw)
         endif
*
!!---- storing in case of weighted generation ----
         if (store.eq.'yes'
     .        .and.eventlimiter.eq.'w'
     .        .and.istopsearch.eq.1) then
            call eventstorage(sdif/wnorm,k,ie,p1,p2,p3,p4,qph)
         endif
!!--------------------------------------------------
         
**[[WRITING OUTPUT]]
         writecond = kt.eq.npoints.or.hit.eq.npoints
         if (iwriteout.gt.0) then
            writecond = mod(k,iwriteout).eq.0.or.writecond
         else
            if (mod(k,30000).eq.0) then
               call cpu_time(tmptime)
               if ((tmptime-starttime).ge.-iwriteout) then
                  writecond = .true.
               endif
            endif
         endif
         
         if (writecond) then
ccc   call printoutput
            
            open(42,file=trim(path)//'iwriteout',status='unknown')
            read(42,*)iwriteout
            close(42)            
            
            kwritenow = k
            idk   = kwritenow-kwritelast
            kwritelast = kwritenow
            itime2 = time8()
            call cpu_time(endtime)
            
            call getrndstatus(rlxstate)            
!! hit or miss cross section...
            hmxsect = 0.d0
            hmerr   = 0.d0
            if (hitpmiss.gt.0.d0)  then 
               hmeff   = hit/hitpmiss
               hmxsect = fmax*hmeff
               hmerr   = fmax* sqrt(hmeff*(1-hmeff)/hitpmiss)
            endif
!!            
            xsec = sum/k
            var  = sqrt((abs(sum2/k-xsec**2))/k)            
            open(10,file=outfile,status='unknown',encoding='UTF-8')

            write(10,*)'RNG status at this point:'
            write(10,*)rlxstate
            write(10,*)' '
            
            write(10,*)' Process PID:',idproc,' running on ',hstnm

            dtime     = endtime-starttime
            dtimetrue = endtime-truestarttime
            write(10,*)' file written last time ~',itime2-itime,'s ago'

            write(10,*)' CPU time to generate last         ',idk,
     .           'events',dtime,'s (',dtime/idk*1d6,
     .           greekmu,'s/ev)' 
            write(10,*)' CPU time to generate     ',k,'events',
     .           dtimetrue,'s (',dtimetrue/k*1d6,greekmu,'s/ev)' 
            starttime = endtime
            itime = itime2
            write(10,*)' '
            write(10,*)'Inputs & parameters for this run: '
            iqmu = qmu
            write(10,'(1x,A,I4)')'Incoming muon charge',iqmu
            write(10,'(1x,A,f12.4,A)')
     .           'nominal muon beam energy =',emulab,' GeV'
            write(10,'(1x,A,f12.4,A)')
     .           'muon beam energy  spread =',semu/emulab*100d0,' %'
            write(10,'(1x,A,f12.4,A)')
     .           'sqrt(s)   =',ecmsnom,' GeV'
            write(10,'(1x,A,f12.4,A)')
     .           'minimum electron LAB energy =',eemin,' GeV'
            write(10,'(1x,A,f12.4,A)')
     .           'maximum electron LAB angle  =',themax,' mrad'
            write(10,'(1x,A,f12.4,A)')
     .           'minimum muon LAB angle      =',thmumin,' mrad'
            write(10,'(1x,A,f12.4,A)')
     .           'maximum muon LAB angle      =',thmumax,' mrad'
            write(10,'(1x,A,i2,1x,f12.4,A)')
     .           'applying acopl. cut and value =',ina7,dthna7max,
     .           ' mrad'
            write(10,'(1x,A)')
     .           'ord    = '//ord
            write(10,'(1x,A)')
     .           'model  = '//model
            write(10,'(1x,A,i5)')
     .           'nphot mode =',nphotmode
            write(10,'(1x,A,i9,i9)')
     .           'seeds   =',iseedext,iseed
            write(10,'(1x,A,i5)')
     .           'alpha running  =',iarun
            if (iarun.eq.1) then
               if (iteubn.eq.1) then
                  write(10,*)'(Teubner & al. alpha running)'
               elseif (inskvp.eq.1) then
                  write(10,*)'(Novosibirsk alpha running)'
               else
                  write(10,*)'(Jegerlehner alpha running)'
               endif
            endif
            write(10,'(1x,A,f10.9)')
     .           'eps    = ',eps
            write(10,*)' '
            write(10,*)'Weights normalization ',wnorm,' (',greekmu,'b)'
            write(10,*)'fmax for unweightening ',fmax
            write(10,*)'true sdifmax ',sdiftruemax
            
            write(10,*)' '
            if (eventlimiter.eq.'w') then
               write(10,'(A,f12.0,A)')'~ Generating ',wnpoints,
     .              ' weighted events ~'
            else
               write(10,'(A,i12,A)')'~ Generating ',npoints,
     .              ' unweighted events ~'
            endif
            write(10,*)' '
            write(10,*)'::::::>>>>>> weighted events <<<<<<::::::'
            do i = 0,nphmax
               xsecp(i) = sump(i)/k
               varp(i)  = sqrt((abs(sum2p(i)/k-xsecp(i)**2))/k)
               fractions(i) = xsecp(i)/xsec * 100.d0
               iii = i
               if (fs.eq.'gg') iii = i + 2
               write(10,'(i2,A,f22.11,A,f22.11,A,f8.4,A)')
     :              iii,'  photons: ',xsecp(i),' +-',varp(i),
     :              '     (',fractions(i),' %)'
            enddo
            if (sumpairs.gt.0.d0) then
               xsecpairs = sumpairs/k
               xsecpairs1 = sumpairsv(1)/k
               xsecpairs2 = sumpairsv(2)/k
               varpairs  = sqrt((abs(sumpairs2/k-xsecpairs**2))/k)
               varpairs1  = sqrt((abs(sumpairsv2(1)/k-xsecpairs1**2))/k)
               varpairs2  = sqrt((abs(sumpairsv2(2)/k-xsecpairs2**2))/k)
               fracpairs = xsecpairs/xsec * 100.d0
               write(10,'(A,f22.11,A,f22.11,A,f8.4,A)')
     :              ' real pairs: ',xsecpairs,' +-',varpairs,
     :              '     (',fracpairs,' %)'
               write(10,'(A,f17.11,A,f22.11)')
     :              '   [e+e-   pairs: ',xsecpairs1,' +-',varpairs1
               write(10,'(A,f17.11,A,f22.11,A)')
     :              '    mu+mu- pairs: ',xsecpairs2,' +-',varpairs2,']'
            endif
            write(10,'(1x,A,f22.11,A,f22.11,A,A,A)')
     :           'total:      ',xsec,' +-',var,' ',greekmu,'b'

            write(10,*)' '
            write(10,'(1x,A,i12,A,f20.0)')
     :           'n. ',k,' of ',wnpoints
            eff = (1.d0*naccepted)/k
            write(10,'(1x,A,f12.8,A)')
     :           'cut points ',abs(1.d0 - eff)*100d0,' %'
            write(10,*)'::::::>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<::::::'
            write(10,*)' '
            write(10,*)'::::::>>>>>> unweighted events <<<<<<::::::'
            write(10,'(1x,A,f10.6,A)')
     :           'hit or miss efficiency ',hmeff*100,' %'
            write(10,'(1x,A,f12.0)')
     :           'unweighted events generated ',hit
            biashit        = 0.d0
            biashitpmiss   = 0.d0
            biasneghit     = 0.d0
            biasneghitmiss = 0.d0
            sezover        = 0.d0
            errsezover     = 0.d0
            if (hit.gt.0.d0) then 
               biashit      = 1.d0*nover/hit
               biashitpmiss = 1.d0*nover/hitpmiss
               biasneghit   = 1.d0*nneg/hit
               biasneghitmiss = 1.d0*nneg/hitpmiss
               sezover    = sumover/hitpmiss
               errsezover = (sum2over/hitpmiss - sezover**2)/hitpmiss
               errsezover = sqrt(abs(errsezover))
               sezneg    = sumneg/hitpmiss
               errsezneg = (sum2neg/hitpmiss - sezneg**2)/hitpmiss
               errsezneg = sqrt(abs(errsezneg))
            endif
            write(10,*)' '
            write(10,'(1x,A,A,A,f22.11,A,f22.11,A)')
     .           'total (',greekmu,'b):   ',hmxsect,' +-',hmerr,' +'
            write(10,'(1x,A,f22.11,A,f22.11,A)')
     .           '       !!!!   ',sezover,' +-',
     .           errsezover,' (bias over fmax) +'
            write(10,'(1x,A,f22.11,A,f22.11,A)')
     .           '       !!!!   ',-sezneg,' +-',
     .           errsezneg,' (bias negative)'
            write(10,'(1x,A,f22.11,A,f22.11)')
     .           'total + biases: ',hmxsect
     .           +sezover-sezneg,' +-',hmerr+errsezover+errsezneg
            write(10,*)' '
            write(10,'(1x,A,i12)')
     .           'N. points with w > fmax (bias): ',nover
            write(10,'(1x,A,f10.7,A,f10.7,A)')
     .           'bias/hit and bias/(hit+missed):',
     .           biashit*100,' % and ',biashitpmiss*100,' %'
            write(10,'(1x,A,i12)')
     .           'N. points with w < 0:',nneg
            write(10,'(1x,A,f10.5,A,f10.5,A)')
     .           'biases for w < 0:   ',
     .           biasneghit*100,' % and ',biasneghitmiss*100,' %'
            write(10,*)' '
            write(10,*)
     .     'True upper/lower limits:'
            write(10,*)sdiftruemax,sdiftruemin

            nwhenprint = nwhenmax
c            write(10,'(1x,A,i3,A)')
c     :           'when there were ',nwhenprint,' photons'
            write(10,*)'::::::>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<::::::'
            write(10,*)' '
            write(10,*)'reweighted cross sections'
c            xsec = sum/k
c            var  = sqrt((abs(sum2/k-xsec**2))/k)
            xsw1 = sumw1/k
            write(10,*)'W1  ',xsw1,'+-',sqrt((abs(sumw12/k-xsw1**2))/k)
            xsw2 = sumw2/k
            write(10,*)'W2  ',xsw2,'+-',sqrt((abs(sumw22/k-xsw2**2))/k)
            xsLO = sumLO/k
            write(10,*)'LO  ',xsLO,'+-',sqrt((abs(sumLO2/k-xsLO**2))/k)
            xsNLO = sumNLO/k
            write(10,*)'NLO ',xsNLO,'+-',
     .           sqrt((abs(sumNLO2/k-xsNLO**2))/k)
            close(10)

            call writedistributions(1.d0*k,nw)
cc            call writedistributions(hit/xsec) ! for unweighted events
            
            
         endif
**[[END WRITING OUTPUT]]
         if (eventlimiter.eq.'unw'.and.hit.ge.npoints) goto 100

****  fifo to syncronize with other process
c         write(33,*)1d3
c         read(33,*)xxxx
****
      enddo
 100  continue
* [[[[[ END LOOP OVER EVENTS ]]]]]      
c      close(33)
      if (store.eq.'yes') then
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
         print*,'...finalizing storage'
         call finalizestorage(storefile)
      endif

      call exitlooplibraries
      call exitrecola
      print*,'Generation finished'
      stop
      end
