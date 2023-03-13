      function init_mesmer(indatacard) bind(C,name="init_mesmer")
      use iso_c_binding
      character(kind=c_char) :: indatacard(500)
      character*500 datacard
      integer init_mesmer
      integer isstandalone      
      common/standalone/isstandalone,datacard ! set in main.F or init_mesmer() !!
      integer k
      logical needswarmup,imwarmingup
      common/warmuppery/needswarmup,imwarmingup
      integer mcids(20),nfspart,itag,ierr
      integer*8 ievtnr
      double precision pmubeam(0:3)
      double precision pmat(0:3,20),cwvpout(0:10)
      double precision wout,wnovp,wnohad,wLO,wNLO
      double precision emulab,eemin,eemax,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela
      integer ina7,iela      
      common/mueexpsetup/emulab,eemin,eemax,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela
      integer ismualongz,itmp
      common/ISmuondirection/ismualongz
      integer initerr
      common/init_error_code/initerr ! filled in user interface

      isstandalone = 0
      do k = 1,500
         datacard(k:k) = ' '
      enddo
      k = 1
      do while (indatacard(k).ne.C_NULL_CHAR)
         datacard(k:k) = indatacard(k)
         k = k + 1
      enddo
      datacard = trim(datacard)
      call userinterface

      if (initerr.eq.0) then      
!     I do here the warmup
         if (needswarmup) then
            itmp = ismualongz
            ismualongz  = 1 ! for the warming up, the mu is always alongz, and gauss-spread
            imwarmingup = .true.      
            call generate_event(pmubeam,nfspart,mcids,pmat,wout,itag,
     .           ievtnr,wnovp,wnohad,wLO,wNLO,cwvpout,ierr)
            ismualongz = itmp
      endif
      imwarmingup = .false.      
!!!!!!!!!!
      else
         print*,'Something wrong in the initialization.'
         print*,'Exiting with an error.'
      endif
      init_mesmer = initerr
      return
      end
*************************************************
      subroutine IncomingMuonMomentum(pmu)
     .        bind(C,name = "IncomingMuonMomentum_mesmer")
***   inputs & outputs *****************
      implicit double precision (a-h,o-z)
      double precision pmu(0:3)
***
      double precision truep1lab(0:3)
      common/mueexpsetup/emulab,eemin,eemax,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela      
      call beamprofile(emulab,semu,emuspr,truep1lab)
      pmu = truep1lab
      return
      end
******************************************
      subroutine generate_event_mesmer(pmubeam,nfspart,mcids,
     .     pmat,w,itag,ievtnr,wnovp,wnohad,wLO,wNLO,cwvp,ierr)
     .     bind(C,name = "generate_event_mesmer")
      use iso_c_binding
*** inputs & outputs *****************
      integer itag,ierr
      integer*8 ievtnr
      double precision pxmu,pymu,pzmu
      double precision pmubeam(0:3)
      double precision pmat(0:3,20),cwvp(0:10)
      double precision w,wnovp,wnohad,wLO,wNLO
      integer(kind=c_int) :: nfspart
      integer(kind=c_int) :: mcids(20)
      call generate_event(pmubeam,nfspart,mcids,pmat,w,itag,
     .     ievtnr,wnovp,wnohad,wLO,wNLO,cwvp,ierr)
      return
      end
************************************************
      subroutine finalize_mesmer(xsw,exsw,foonhpm,foonh,truemax,nabove,
     .     nlt0,xsbias,exsbias,xsbiasn,
     .     exsbiasn,sumow,sum2ow2,sumnow,sum2now2)
     .     bind(C,name="finalize_mesmer")
      use iso_c_binding
      implicit double precision (a-h,o-z)
      integer*8 foonhpm,foonh,nabove,nlt0
      call get_footer(xsw,exsw,foonhpm,foonh,truemax,nabove,nlt0,xsbias,
     .     exsbias,xsbiasn,
     .     exsbiasn,sumow,sum2ow2,sumnow,sum2now2)
      call exitlooplibraries
      call exitrecola
      return
      end
*************************************************
      subroutine get_footer(xsw,exsw,foonhpm,foonh,truemax,nabove,nlt0,
     .     xsbias,exsbias,xsbiasn,
     .     exsbiasn,sumow,sum2ow2,sumnow,sum2now2)
      implicit double precision (a-h,o-z)
      integer*8 nabove,nlt0
      real starttime,endtime,truestarttime,tmptime
      parameter (mph=40,maxnw=3)
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
      integer*8 foonhpma,foonha
      double precision ntomicro
      integer rlxstate(104)
! i use the forfooter common only for sum and sum2      
      common/forfooter/xswa,exswa,xsunwa,exsunwa ,
     .     xsbiasa,exsbiasa,xsbiasna,exsbiasna,truemaxa,sum,sum2,sumneg,
     .     sum2neg,foonhpma,foonha,nabovea,nlt0a
      common/onceinmain/truestarttime,starttime,s,sumw1,sumw2,
     .     sumw12,xsec,var,hmxsect,hmerr,
     .     sumw22,sumLO,sumLO2,sumNLO,sumNLO2,sumpairs,sumpairs2,
     .     sumpairsv,sumpairsv2,sump,sum2p,sdiftruemax,sdiftruemin,hit,
     .     hitpmiss,istopsearch,nneg,sumover,sum2over,fmax,
     .     ntomicro,k,npoints,naccepted,kt,itime,rlxstate,nphmax,nover, 
     .     in_conf_spin
      common/weightnormalization/wnorm
*************
      integer iseed,nsearch,iverbose,nw,isync,
     .     istorver,icwvp,ivpfl      
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync,
     .     istorver,icwvp,ivpfl
      logical areweighted
      common/weightedeventsornot/areweighted
      
c      if (k.le.0) then
c! in common         xsec = 0.d0
c! in common         var  = 0.d0
c      else
c! in common         xsec = sum/k
c! in common         var  = sqrt((abs(sum2/k-xsec**2))/k)
c      endif
*** writing anyway outputs and distributions      
      call printoutput
      call writedistributions(1.d0*k,nw)
***
      sezover    = 0.d0
      errsezover = 0.d0
      sezneg = 0.d0
      errsezneg = 0.d0
      if (hitpmiss.gt.0) then
         sezover    = sumover/hitpmiss
         errsezover = (sum2over/hitpmiss - sezover**2)/hitpmiss
         errsezover = sqrt(abs(errsezover))
         sezneg    = sumneg/hitpmiss
         errsezneg = (sum2neg/hitpmiss - sezneg**2)/hitpmiss
         errsezneg = sqrt(abs(errsezneg))
      endif
      
      foonhpm = int(hitpmiss,kind=8)      

      foonh   = int(hit,kind=8)
      if (areweighted) foonh = foonhpm

      xsw     = xsec
      exsw    = var
      xsunw   = hmxsect
      exsunw  = hmerr

      nabove   = int(nover,kind=8)
      nlt0     = int(nneg,kind=8)
      xsbias   = sezover
      exsbias  = errsezover
      truemax  = sdiftruemax/wnorm
      xsbiasn  = sezneg
      exsbiasn = errsezneg
      sumow    = sum/wnorm
      sum2ow2  = sum2/wnorm/wnorm
      
      sumnow   = sumneg/wnorm
      sum2now2 = sum2neg/wnorm/wnorm
      return
      end      
*************************************************
      subroutine get_masses_mesmer(mumass,emass)
     .     bind(C,name="get_masses_mesmer")
      implicit double precision (a-h,o-z)
      double precision mumass,emass,taumass
      common/parameters/ame,ammu,convfac,alpha,pi
      common/tautopmasses/taumasscmn,topmasscmn
      mumass  = ammu
      emass   = ame
      taumass = taumasscmn
      return
      end
*************************************************
      subroutine mesmer_setup
     .     (ctagout,
     .     mesmerversion,
     .     hostnameout,
     .     datetimeout,
     .     idprocout,
     .     nevout,
     .     areweightedout,
     .     RCorderout,
     .     includedfsout,
     .     radmuout,
     .     radeout,
     .     iseed1,
     .     iseed2,
     .     emulabout,
     .     spreadout,
     .     extmubeamout,
     .     Qmuout,
     .     mumassout,
     .     elmassout,
     .     invalphaout,
     .     wnormout,
     .     wmaxout,
     .     eeminout,
     .     theminout,
     .     eemaxout,
     .     thmuminout,
     .     thmumaxout,
     .     ethrout,
     .     ththrout,
     .     iacoplout,
     .     acoplout,
     .     ielaout,
     .     elaout,
     .     ivpwgtsout,
     .     ihadonout,
     .     ivpflout,
     .     nwarmupout,
     .     ndistrwout,
     .     isyncout,
     .     epsout,
     .     phmassout)
     .     bind(C,name="mesmer_setup")
      use iso_c_binding
      implicit double precision (a-h,o-z)
*** output vars
      integer ctagout
      character(kind=c_char), dimension(21) :: mesmerversion
c      character :: mesmerversion(40)
      character(kind=c_char),dimension(21) :: hostnameout
      character(kind=c_char),dimension(21) :: datetimeout
      character(kind=c_char),dimension(5)  :: RCorderout

      integer idprocout
      integer*8 nevout
      integer includedfsout
      integer radmuout,radeout
      integer iseed1,iseed2
      integer ivpwgtsout
      double precision emulabout
      double precision spreadout
      double precision Qmuout
      double precision mumassout,elmassout
      double precision invalphaout
      double precision wnormout
      double precision wmaxout
      double precision eeminout
      double precision theminout
      double precision eemaxout
      double precision thmuminout
      double precision thmumaxout
      double precision ethrout
      double precision ththrout
      double precision acoplout
      double precision elaout
      double precision epsout
      double precision phmassout
      integer iacoplout
      integer ihadonout
      integer ielaout
      integer ivpflout
      integer nwarmupout
      integer ndistrwout
      integer isyncout      
      integer areweightedout
      integer extmubeamout
******************************************************      
      integer*8 nev,nstored,kbefore,nevts
      dimension qph(40,0:3),qph1(0:3),qph2(0:3)
      dimension pin1(0:3),pin2(0:3),p5(0:3),p6(0:3),ptmp(0:3)
      common/momentainitial/pin1,pin2
      common/storageunit/istorunit,ilens
      common/storedevents/nstored,nevts,ifirst

      character*30 cnstored,cnp,ctmp,ctag,ctmp2
      character(len=:),allocatable,target :: cstring
      
      common/parameters/ame,ammu,convfac,alpha,pi
      common/mueexpsetup/emulab,eemin,eemax,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela
      integer*8 iwriteout
      integer istorver,icwvp
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync,
     .     istorver,icwvp,ivpfl
      integer ihadoff
      common/hadoff/ihadoff
      
      parameter (ikind = 8)
      double precision p3lab(0:3),p4lab(0:3)
      double precision p5lab(0:3),p6lab(0:3)
      double precision qphlab(40,0:3)

      integer g,mum,mup,em,ep
      common/LHAparticlesids/g,mum,mup,em,ep
      character*4 cg,cmum,cmup,cem,cep
      common/LHAparticlesidschar/cg,cmum,cmup,cem,cep

      double precision mz,mz2,gz,gz2,mw,mw2,stw,ctw
      common/zparameters/mz,mz2,gz,gz2,mw,mw2,stw,ctw,gv,ga
      
      integer id1,id2,id3,id4,id5,id6
      character*3 cid1,cid2,cid3,cid4,cid5,cid6
      
      common/eventcount/kbefore
      
      character*20 programversion
      common/version/programversion
      character*20 hstnm
      common/hostandpid/hstnm,idproc

      common/iseedoninput/iseedext
      character*8 date
      character*10 time
      character*5 zone
      integer values(8)
      
      character*6   ord
      character*10  model
      character*100 outfile,storefile
      character*3   eventlimiter,store
      common/qedORDER/ord
      common/charinput/model,eventlimiter,store,storefile,outfile
      common/photoncutoff/egmin,egmin2
      common/realinput/anpoints,sdifmax

      double precision eps
      common/epssoft/eps
      
      common/nphot_mode/nphotmode
      
      common/weiunwei/iweight

      common/tobewrittenintheheader/wnormalization,sdmax
      
      double precision p0(0:3)
      data p0 /0.d0,0.d0,0.d0,0.d0/

      integer ismualongz
      common/ISmuondirection/ismualongz

      double precision phmass
      common/photonmasslambda/phmass
      
      character*2 c1,c2,c4,c5,c6,c7
      character*4 c3
*     included only to access the common 'weights', QRmu QRe and coeffforreweighting
#include "invariants.h"
      logical areweighted
      common/weightedeventsornot/areweighted

      ctagout = iseedext
      
      call date_and_time(date,time,zone,values)
      write(c1,'(i2)')values(3)
      if (values(3).lt.10) write(c1,'(i1,i1)')0,values(3)
      write(c2,'(i2)')values(2)
      if (values(2).lt.10) write(c2,'(i1,i1)')0,values(2)
      write(c3,'(i4)')values(1)
      write(c4,'(i2)')values(5)
      if (values(5).lt.10) write(c4,'(i1,i1)')0,values(5)
      write(c5,'(i2)')values(6)
      if (values(6).lt.10) write(c5,'(i1,i1)')0,values(6)
      write(c6,'(i2)')values(7)
      if (values(7).lt.10) write(c6,'(i1,i1)')0,values(7)
      write(ctmp,*)c1//'/'//c2//'/'//c3//' '//c4//':'//c5//':'//c6
     
      ctmp = adjustl(trim(ctmp))
      do k = 1,min(20,len_trim(ctmp))
         datetimeout(k:k) = ctmp(k:k)
      enddo
      datetimeout(k:k) = C_NULL_CHAR
      
      ctmp = adjustl(trim(programversion))
      do k = 1,min(20,len_trim(ctmp))
         mesmerversion(k:k) = ctmp(k:k)
      enddo
      mesmerversion(k:k) = C_NULL_CHAR
      
      ctmp = adjustl(trim(hstnm))
      do k = 1,min(20,len_trim(ctmp))
         hostnameout(k:k)   = ctmp(k:k)
      enddo
      hostnameout(k:k) = C_NULL_CHAR
      
      if (ord.eq.'born') then
         RCorderout(1:1) = 'L'
         RCorderout(2:2) = '0'
      elseif (ord.eq.'alpha') then
         RCorderout(1:1) = 'N'
         RCorderout(2:2) = 'L'
         RCorderout(3:3) = 'O'
      elseif (ord.eq.'alpha2') then
         RCorderout(1:1) = 'N'
         RCorderout(2:2) = 'N'
         RCorderout(3:3) = 'L'
         RCorderout(4:4) = 'O'
      endif
      RCorderout(5:5) = C_NULL_CHAR
      
      idprocout     = idproc
      nevout        = int(anpoints,kind=8)
      
      areweightedout = 0
      if (areweighted) areweightedout = 1
            
      includedfsout = nphotmode
      radmuout = QRmu
      radeout  = QRe

      iseed1 = iseedext
      iseed2 = iseed

      emulabout = emulab
      spreadout = semu

      extmubeamout = 1 - ismualongz

      Qmuout = Qmu

      mumassout = ammu
      elmassout = ame
      invalphaout = 1.d0/alpha

c      k0out = egmin
      epsout = eps
      phmassout = phmass
      
      wnormout = wnorm
      wmaxout  = sdmax/wnormout

      eeminout   = eemin
      theminout  = themin
      eemaxout   = eemax
      thmuminout = thmumin
      thmumaxout = thmumax
      ethrout    = ethr
      ththrout   = ththr
      iacoplout  = ina7
      acoplout   = dthna7max
      ielaout    = iela
      elaout     = cutela
***
      ihadonout  = 1 - ihadoff
      ivpflout   = ivpfl
      nwarmupout = nsearch
      ndistrwout = nw
      isyncout   = isync
***
      ivpwgtsout = icwvp      
      return
      end
*************************************************
      subroutine simplechar(stringa) bind(C,name="simplechar")
! test routine
      use iso_c_binding
      character(kind=c_char) stringa
c      type(c_ptr) stringa
      character*200 fortstring
      fortstring = "334abcdefg 1 1"//C_NULL_CHAR
      print*,len(fortstring)      
      do k = 1,len(fortstring)
         stringa(k:k) = fortstring(k:k)
      enddo
      return
      end
