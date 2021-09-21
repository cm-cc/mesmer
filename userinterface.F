      subroutine userinterface
      use, intrinsic :: iso_fortran_env, only : stdin=>input_unit,
     .                                        stdout=>output_unit,
     .                                         stderr=>error_unit
      implicit double precision (a-b,d-h,o-z)
      implicit character*10 (c)
      logical istop
      character*50 r8tostr,r4tostr,itostr,i8tostr
      double precision csi(2)
      character*20  c,ctmp
      character*100 cinput
      character*6   ord,arun,chadoff
      character*10  model
      character*100 outfile,storefile
      character*200 path
      character*100 cx
      character*100 stringa
      character*3   eventlimiter,store
      character*200 programma,corrections
      double precision convfac

      character*2 fs
      common/finalstate/fs
      common/ecms/ecms,ecmsnom
      common/nphot_mode/nphotmode
      
      common/epssoft/eps
      common/parameters/ame,ammu,convfac,alpha,pi
      common/tautopmasses/taumass,topmass
      double precision mz,mz2,gz,gz2,mw,mw2,stw,ctw
      common/zparameters/mz,mz2,gz,gz2,mw,mw2,stw,ctw,gv,ga
      common/qedORDER/ord
      integer*8 iwriteout
      integer isync,seed,nsearch,iverbose,nw
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync
      common/charinput/model,eventlimiter,store,storefile,outfile
      common/cmnpath/path
      common/realinput/anpoints,sdifmax

      common/hadoff/ihadoff
      integer ififo
      common/ififocmn/ififo

******************************
      integer ifredv
      common/fredversion/ifredv
      data ifredv/19/ ! or 12 or 17
*****************************      
*** from invariants.h      
      double precision Qmu
      common/muoncharge/Qmu
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
***
      double precision cutela
      common/mueexpsetup/emulab,eemin,semu,thmumin,themax,thmumax,
     .     dthna7max,cutela,ina7,iela

      double precision phmass
      common/photonmasslambda/phmass
      
      common/ialpharunning/iarun

!!!     should be initialized here! used for real pairs
      integer iper,noper,iid
      common/identicalandperipheral/iper,noper,iid
      data iper,iid /0,0/
!!!      
      common/weightnormalization/wnorm
      character*9 setkey(37)
c$$$      
c$$$ for teubner routine, 2010/05/19
c$$$
      integer iteubn
      common/teubner/iteubn
c$$$
c$$$
c$$$ for Novosibirsk VP routine, 2015/01/22
c$$$
      integer inskvp
      common/nskvp/inskvp
c$$$

      data setkey/
     +     'Ebeam  ', !  1
     +     'bspr   ', !  2
     +     'Eemin  ', !  3
     +     'Qmu    ', !  4
     +     'radchs ', !  5
     +     'nev    ', !  6
     +     'path   ', !  7
     +     'store  ', !  8
     +     'menu2  ', !  9
     +     'seed   ', ! 10
     +     'ord    ', ! 11
     +     'model  ', ! 12
     +     'nphot  ', ! 13
     +     'mode   ', ! 14
     +     'eps    ', ! 15
     +     'nwrite ', ! 16
     +     'nwarmup', ! 17
     +     'verbose', ! 18
     +     'sdmax  ', ! 19
     +     'help   ', ! 20
     +     'arun   ', ! 21
     +     'thmumin', ! 22
     +     'thmumax', ! 23
     +     'themax ', ! 24
     +     'hadoff ', ! 25
     +     'phmass ', ! 26
     +     'storemode', ! 27
     +     'ndistr', ! 28      
     +     'acoplcut ', ! 29
     +     'sync', ! 30 
     +     'elastcut', ! 31
     +     'massmin', ! 32
     +     'massmax', ! 33
     +     'scEdhad', ! 34
     +     'sprb1  ', ! 35
     +     'sprb2  ', ! 36
     +     'wnorm '   ! 37
     +     /
* BSM parameters
      common/bsmparams/bsmmass,yae,yam,bsmlambda,gagg,ibsm
      data ibsm /0/
***
      character*20 programversion
      common/version/programversion
      common/iseedoninput/iseedext
      character*200 welcome(23)

      character*4 bbf,ebf
      character*4 bit,eit
      integer ionlybubbles,ibubbleonv,ibubbleonboxsub,iexnnlobubbles
      common/ionlybubblescmn/iexnnlobubbles,ionlybubbles,ibubbleonv,
     .     ibubbleonboxsub
      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
            
*** opening stdout as utf-8!      
      open(stdout,encoding='UTF-8')
***

*** for real pairs
      iid   = 1
      iper  = 1
      noper = 1
*** for virtual pairs
      ionlybubbles    = 0
      ibubbleonv      = 0
      ibubbleonboxsub = 0
      iexnnlobubbles  = 0
c      bubblemask = (/1,1,1,1,1/)
      bubblemask = (/0,0,0,0,0/)

#ifndef RELEASE
      print*,'iid iper and noper'
      read*,iid,iper,noper
      
**********************      
*******************************

c      print*,'Input ionlybubbles'
c      read*,ionlybubbles
      print*,'Input iexnnlobubbles (0/1/2'
      read*,iexnnlobubbles
      
      print*,'only bubbles on verter and/or boxes'
      read*,ibubbleonv,ibubbleonboxsub
      
      print*,'bubble mask?'
      read*,bubblemask
      print*,bubblemask
c      stop
#endif
      
      colch2(1:3) = 1.d0      * bubblemask(1:3)
      colch2(4)   = 4.d0/3.d0 * bubblemask(4)    ! i.e. for top 3 colors * (2/3 charge)^2 = 12/9 = 4/3
      colch2(5)   = 1.d0      * bubblemask(5)    ! for Rhad!
*******************************      
      
      bbf=achar(27)//'[1m'
      ebf=achar(27)//'[0m'
      bit=achar(27)//'[3m'
      eit=achar(27)//'[0m'
      
      programversion='MESMER 1.0'
      welcome = ' '      
      welcome(1) =
     . '12345678901234567890123456789012345678901234567890123456789012' !1
      welcome(2) =
     .' *************************************************************' !2
      welcome(3) =
     .' ********                                             ********'!4
      welcome(4) =
     . ' ******               Welcome to '//bbf//'MESMER'//ebf//'
     .          ******'!5
      welcome(5) =
     .' ****              ~~~~~~~~~~~~~~~~~~~~~~                 ****'!6
      welcome(6) =
     .' **       A fully exclusive Monte Carlo event generator     **'!7
      welcome(7) =
     .' ****          for '//bbf//'M'//ebf//'uon '//bbf//'E'//ebf//
     .'lectron '//bbf//'S'//ebf//'cattering  with         ****'           !8
      welcome(8) =
     .' ******        '//bbf//'M'//ebf//'ultiple '
     .//bbf//'E'//ebf//'lectromagnetic '
     .//bbf//'R'//ebf//'adiation       ******'!9
      welcome(9) =
     .' ********                                             ********'!10
      welcome(10) =
     .' *************************************************************'!12
      welcome(11) ='   Questo è quel pezzo / di calamita:'
      welcome(12) ='   pietra '//bbf//'mesmerica'//ebf//
     .     ', / ch''ebbe l''origine'
      welcome(13) ='   nell''Alemagna, / Che poi sì celebre'//
     .     ' / là in Francia fu.'
      welcome(14) ='          Mozart/Da Ponte, '//bit//
     .'Così Fan Tutte'//eit//', Act II, Scene XVI'

      do k = 2,14
         print*,trim(welcome(k))
      enddo
**
      convfac = 0.3893793656D6
**
      pi = 3.1415926535897932384626433832795029D0
      do i = 1,10
         model(i:i) = ' '
      enddo
      do k = 1,200
         path(k:k) = ' '
      enddo
      do i = 1,200
         programma(i:i)   = ' '
      enddo
      do i = 1,50
         corrections(i:i) = ' '
      enddo

*** BSM parameter arXiv:1607.01022, Passera et al.
      bsmlambda = 1000d0
      yae       = 0.d0
      yam       = 0.d0
      gagg      = 0.d0
      bsmmass   = 0.5d0
*** end BSM parameters
      
      ame  =   0.51099906d-3
      ammu = 105.65836900d-3

**
      ame  = 0.510998928d-3 ! electron mass
      ammu = 105.6583715d-3 ! muon mass
      taumass = 1.77682D0
      topmass = 175.6D0

** for Z exchange (at LO)      
      mz = 91.1876d0
      gz = 2.4952d0
      mw = 80.398d0
      
c      ame = ammu
c      ammu = ame
      
c      write(*,'(A,A)',advance='no')'test',achar(13)
c      stop
*
      alpha = 1.d0/137.03599907430637d0

***** default input values
      fs   = 'ee' ! 'ee' or 'gg' or 'mm'
      path = 'test-run/'
      ord  = 'alpha'
      nphotmode = -1
      model  = 'matched'
      anpoints = 10000000.d0

      emulab = 150.d0
      spremu = 0.d0
      eemin  = 1.d0
      thmumin = 0.d0
      thmumax = 100d0
      themax  = 100d0

      ecms = ame*ame + ammu*ammu + 2.d0*emulab*ame
      ecms = sqrt(ecms)      
      
      ihadoff = 0
      chadoff = 'no'

      QRmu = 1
      QRe  = 1
      
      iqmu     = 1
      Qmu      = iqmu*1.d0

      phmass = 1d-10
      eps    = 2.d-5

      ina7 = 0
      iela = 0
      dthna7max = 3.5d0
      cutela    = 0.2d0
      
      nw = 1 ! number of distributions

      iseed        = 42
      iwriteout    = 500000
      isync        = 0
      
      nsearch      =  0 ! 1000000
      wnorm        =  1250.d0 ! -1d0
      wnorm        = -1.d0
      
      iverbose     = 0
      eventlimiter = 'w'  ! 'unw' or 'w'
      store        = 'no'
c      ififo = 2 ! for gzipped file on the fly
c      ififo = 1 ! for writing to a fifo, used with root interface      
c      ififo = 0 ! ascii file
      ififo = 0
      
      sdifmax      = 1.d-17
      arun         = 'on'

      iarun = 0
      if (arun.ne.'off') iarun = 1
***********
      
      istop = .false.
      iii   = 0
      do while(.not.istop)
         print*,' '
         if (iii.gt.0) then
c            print*,'>> Overwriting run parameters <<'
         endif
         iii = 1
         print*,bbf//'Principal parameters:'//ebf
         print*,'  [ type "run" to start generation, '//
     .        '"',trim(setkey(20)),'" for help or "quit" to abort ]'
         write(stringa,*)' [ '//setkey(4)//' ] Incoming muon charge = ',
     .        trim(itostr(iqmu)),' (in e+ charge units)'
         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(1)//' ] Muon beam energy = ',
     .        trim(r8tostr(emulab)),' GeV'
         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(2)//' ] Beam energy spread = ',
     .        trim(r8tostr(spremu)),' %'
         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(3)//' ] Minimum electron'//
     .       ' LAB energy = ',trim(r8tostr(eemin)),' GeV'
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(24)//' ] Maximum electron'//
     .       ' LAB angle = ',trim(r8tostr(themax)),' mrad'
         print*,trim(stringa)
         
         write(stringa,*)' [ '//setkey(22)//' ] Minimum muon'//
     .       ' LAB angle = ',trim(r8tostr(thmumin)),' mrad'
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(23)//' ] Maximum muon'//
     .       ' LAB angle = ',trim(r8tostr(thmumax)),' mrad'
         print*,trim(stringa)

         ctmp = 'no'
         if (ina7.eq.1) ctmp = 'yes'
         write(stringa,*)' [ '//setkey(29)//' ] Apply acoplanarity cut'
     .        //' (and value) = '//trim(ctmp),' (',
     .        trim(r8tostr(dthna7max)),' mrad)'
         print*,trim(stringa)

         ctmp = 'no'
         if (iela.eq.1) ctmp = 'yes'
         write(stringa,*)' [ '//setkey(31)//' ] Apply "elasticity" cut'
     .        //' (and value) = '//trim(ctmp),' (',
     .        trim(r8tostr(cutela)),' mrad)'
         print*,trim(stringa)
         
         write(stringa,*)' [ '//setkey(11)
     .//' ] Simulation at "'//trim(ord)//'" order'
         print*,trim(stringa)
         
         write(stringa,*)' [ '//setkey(21)
     .//' ] alpha running is '//arun
         print*,trim(stringa)         
         
         write(stringa,*)' [ '//setkey(25)//' ] Hadronic VP is off: '//
     .       trim(chadoff)
         print*,trim(stringa)
         
         write(stringa,*)' [ '//setkey(6)
     .//' ] n. of events to generate: ',trim(r8tostr(anpoints))
         print*,trim(stringa)
         
         write(stringa,*)' [ '//setkey(8)
     .        //' ] events storage: '//store
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(27)
     .//' ] mode to store: ',trim(itostr(ififo))         
         print*,trim(stringa)

         do k = 1,200
            if (path(k:k).ne.' ') lpath=k
         enddo
         if (path(lpath:lpath).ne.'/') then
            lpath = lpath+1
            path(lpath:lpath)='/'
         endif
         write(stringa,*)' [ '//setkey(7)
     .//' ] files saved in '//path(1:lpath)
         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(10)
     .//' ] seed for pseudo-RNG = ',trim(itostr(iseed))
         print*,trim(stringa)         
*********
         print*,' '
         print*,bbf//'Internal tweaks:'//ebf
         ctmp = 'unweighted'
         if (eventlimiter.eq.'w') ctmp = 'weighted'
         write(stringa,*)' [ '//setkey(14)
     .//' ] requested evts. are '//trim(ctmp)
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(5)
     .        //' ] "radiative" charges: muon ~> ',trim(itostr(QRmu)),
     .   ' electron ~> ',trim(itostr(QRe))
         print*,trim(stringa)
         
         write(stringa,*)' [ '//setkey(15)
     .        //' ] soft photon cutoff in sqrt(s)/2 units = ',
     .        trim(r8tostr(eps))
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(26)
     .        //' ] photon mass = ',trim(r8tostr(phmass)),' GeV'
         print*,trim(stringa)
         
c         write(stringa,*)' [ '//setkey(12)
c     .//' ] model for corrections is '//trim(model)
c         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(13)
     .//' ] max. number of photons mode is ',trim(itostr(nphotmode))
         print*,trim(stringa)
         if (iwriteout.gt.0) then
            write(stringa,*)' [ '//setkey(16)
     .           //' ] file(s) dumped every ',trim(i8tostr(iwriteout)),
     .           ' events'
         else
            write(stringa,*)' [ '//setkey(16)
     .           //' ] file(s) dumped every ',
     .           trim(r8tostr(-iwriteout/60.d0)),
     .           ' minutes'
         endif
         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(17)
     .//' ] events for maximum searching = ',trim(itostr(nsearch))
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(28)
     .//' ] number of distr. at different orders: ',trim(itostr(nw))
         print*,trim(stringa)
         
         write(stringa,*)' [ '//setkey(19)
     .        //' ] starting "sdifmax" ',trim(r8tostr(sdifmax))
         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(37)
     .        //' ] normalization cross section = ',trim(r8tostr(wnorm))
         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(30)
     .        //' ] random numbers sequence syncronization: ',
     .        trim(itostr(isync))//' [0/1]'
         print*,trim(stringa)
*********
         print*,' '
cc 111     write(*,'(''Insert "variable value" or "run" or "quit": '',$)')
 111     write(*,42)trim(bbf//'Insert "parameter value" or'
     .        //' "run" or "quit": '//ebf)
 42      format(A54,$)
         read(*,'(A)')cinput
         k = index(cinput,' ')
         c = cinput(1:k-1)
         if (c.eq.'help') then 
            print*,' '
            call explainvalues(setkey)
            print*,' '
            goto 111
         endif
         if ( c.eq.'ord'.or.
     .        c.eq.'model'.or.
     .        c.eq.'path'.or.
     .        c.eq.'store'.or.
     .        c.eq.'arun'.or.
     .        c.eq.'mode'.or.
     .        c.eq.'hadoff'.or.
     .        c.eq.'acoplcut'.or.
     .        c.eq.'elastcut'.or.
     .        c.eq.'radchs') then
            read(cinput(k+1:),'(bn,A)') cx
            if (c.eq.'ord')     ord     = cx
            if (c.eq.'model')   model   = cx
            if (c.eq.'path')    path    = cx
            if (c.eq.'store')   store   = cx
            if (c.eq.'arun')    arun    = cx
            if (c.eq.'hadoff')  chadoff = cx    
            if (c.eq.'radchs')  read(cx,*)QRmu,QRe
            if (c.eq.'mode') then 
               if (cx.eq.'weighted')   eventlimiter = 'w'
               if (cx.eq.'unweighted') eventlimiter = 'unw'
            endif
            if (c.eq.'acoplcut') then
               read(cx,*)ctmp,dthna7max
               if (ctmp.eq.'yes') ina7 = 1
            endif
            if (c.eq.'elastcut') then
               read(cx,*)ctmp,cutela
               if (ctmp.eq.'yes') iela = 1
            endif
         elseif (c.eq.'seed'.or.
     .           c.eq.'nwarmup'.or.
     .           c.eq.'nwrite'.or.
     .           c.eq.'verbose'.or.
     .           c.eq.'nphot'.or.
     .           c.eq.'storemode'.or.
     .           c.eq.'ndistr'.or.
     .           c.eq.'sync'.or.
     .           c.eq.'Qmu') then
            ctmp = cinput(k+1:)
            kk = index(ctmp,'.')
            if (kk.gt.0) ctmp = ctmp(:kk-1)
            read(ctmp,'(bn,i9)')ix
            if (c.eq.'seed')      iseed = ix
            if (c.eq.'nwarmup')   nsearch = ix
            if (c.eq.'nwrite')    iwriteout = ix
            if (c.eq.'verbose')   iverbose = ix
            if (c.eq.'nphot')     nphotmode = ix
            if (c.eq.'sync')      isync = ix
            if (c.eq.'Qmu')       iqmu = ix
            if (c.eq.'storemode') ififo = ix
            if (c.eq.'ndistr')    nw = ix 
         else
            ctmp = cinput(k+1:)
            kk = index(ctmp,'.')
            read(cinput(k+1:),'(bn,f40.20)')x
            if (kk.eq.0) read(cinput(k+1:),'(bn,f40.0)')x
            if (c.eq.'Ebeam')   emulab = x
            if (c.eq.'bspr')    spremu = x
            if (c.eq.'Eemin')   eemin = x
            if (c.eq.'nev')     anpoints = x
            if (c.eq.'eps')     eps = x
            if (c.eq.'phmass')  phmass = x
            if (c.eq.'sdmax')   sdifmax = x
            if (c.eq.'wnorm')   wnorm = x
            if (c.eq.'thmumin') thmumin = x
            if (c.eq.'thmumax') thmumax = x
            if (c.eq.'themax')  themax = x            
         endif
         if (c(1:3).eq.'run')  then
            istop = .true.
            print*,'                  '
         endif
         if (c(1:4).eq.'quit') stop
      enddo

**********************************************************************
      iarun  = 0
      iteubn = 0
      inskvp = 0
      if (arun.eq.'on'.or.arun.eq.'hadr5') then
         iarun  = 1
         iteubn = 0
      elseif (arun.eq.'knt') then
         iarun  = 1
         iteubn = 1
      elseif (arun.eq.'nsk') then
         iarun  = 1
         inskvp = 1
      else
         iarun = iarun
      endif
      if (iarun.eq.1) then
         if (inskvp.eq.1) then
c            print*,'                  '
            print*, 'Using Fedor Ignatov''s routine for '//
     .           'hadronic \Delta\alpha.'
         else
            if (iteubn.lt.1) then
c               print*,'                  '
               if (ifredv.eq.19)
     .              print*, 'Using Fred Jegerlehner''s HADR5x19'//
     .              ' routine for hadronic \Delta\alpha.'
               if (ifredv.eq.12)
     .              print*, 'Using Fred Jegerlehner''s HADR5n12'//
     .              ' routine for hadronic \Delta\alpha.'
               if (ifredv.eq.17)
     .              print*, 'Using Fred Jegerlehner''s HADR5n17'//
     .              ' routine for hadronic \Delta\alpha.'
            else
c               print*,'                  '
               print*,'Using Keshavarzi-Nomura-Teubner''s routine for'//
     .              ' hadronic \Delta\alpha.'//
     .              ' [VP_KNT (version 3.0.1, Feb 16 2021)]'
            endif
         endif
      endif
***********************************************************************            
      do k = 1,200
         if (path(k:k).ne.' ') lpath=k
      enddo
      if (path(lpath:lpath).ne.'/') then
         lpath = lpath+1
         path(lpath:lpath)='/'
      endif

      call execute_command_line('[ -d '//trim(path)//' ]',exitstat=ierr)      
      if (ierr.eq.0) then ! i.e. if the directory exists
         print*,' '
         print*,'WARNING:'
         print*,'The directory ',path(1:lpath),' already exists.'
         print*,'The files inside it are going to be overwritten...'
      endif
      close(11)
      call execute_command_line('mkdir -p '//trim(path)//' 2>/dev/null')

      corrections=ord
      if (ord.eq.'alpha')  corrections='oal'
      if (ord.eq.'alpha2') corrections='oal2'

      if (corrections.eq.'born') then
         storefile = path(:lpath)//'events-LO.dat'
      elseif (corrections(1:4).eq.'oal2') then
         storefile = path(:lpath)//'events-NNLO.dat'
      else
         storefile = path(:lpath)//'events-NLO.dat'
      endif
      
      outfile  = path(:lpath)//
     .     'stat_'//corrections(1:len_trim(corrections))//'.txt'
c
      if (store.eq.'yes') then
*** checking if nsearch, sdifmax and wnorm are consistent      
         if (nsearch.gt.0) then
            print*,'**************************************************'
            print*,'Using nwarmup =',nsearch
            print*,'and resetting wnorm and sdmax'
            print*,'**************************************************'
            wnorm = -1.d0
            sdifmax = 1d-17
         endif
         if (wnorm.gt.0.d0) then
            print*,'**************************************************'
            print*,'Using wnorm=',wnorm
            print*,'and resetting nwarmup to 0'
            print*,'**************************************************'
            nsearch = 0
         endif
         if (eventlimiter.eq.'unw') then
            if (nsearch.lt.1.and.sdifmax.lt.0.001d0) then
               print*,'You should set sdmax > 0 or nwarmup > 0.'
               print*,'Stopping.'
               stop
            endif            
            if (nsearch.lt.1.and.wnorm.lt.0.d0) then
               print*,'You should set wnorm > 0 or nwarmup > 0.'
               print*,'Stopping.'
               stop
            endif
         endif
      endif

      spremu = spremu * 0.01d0
      semu   = spremu * emulab 
      Eemin = max(Eemin,ame)
      if (abs(iqmu).ne.1) then
         print*,'Wrong charge selected! stopping...'
         stop
      endif
      Qmu = iqmu * 1.d0

      if (QRe.lt.0.or.QRe.gt.1) then
         print*,'Wrong electron "radiative" charge.'
         print*,'Must be 0 or 1. Stopping'
         stop
      endif
      if (QRmu.lt.0.or.QRmu.gt.1) then
         print*,'Wrong muon "radiative" charge.'
         print*,'Must be 0 or 1. Stopping'
         stop
      endif

      if (chadoff.eq.'yes') ihadoff = 1
      if (chadoff.eq.'no')  ihadoff = 0

**** number of distributions in one shot
      if (ord.eq.'born')  nw = 1
      if (ord.eq.'alpha') nw = min(2,nw)
      nw = min(nw,3)
*****

*** switch off any NNLO VP RC if iarun = 0  or ord = born
      if (iarun.eq.0) iexnnlobubbles = 0
      if (ord.eq.'born') iexnnlobubbles = 0
***      
*******************************************************
#ifdef RELEASE
c      if (ord.eq.'alpha2'
c     .     .or.(nphotmode.gt.999.or.nphotmode.lt.-999)) then
      if (nphotmode.gt.999.or.nphotmode.lt.-999) then
         print*,' '
         print*,'You want to run with not-yet-released features.'
         print*,'Please, contact Carlo... stopping execution!'
         print*,' '
         stop
      endif
#endif
**************************************************
      print*,'                  '
      print*,'Starting generation...'
      print*,'                  '
****  RNG initialization
      call initrng(4323432)
      iseedext = iseed
      do k = 1,iseed
         call getrndroutine(csi,2)
      enddo
      iii   = int(csi(1)*1d7)
      iii   = iii + int(csi(2)*1d7)
      iseed = iii
      call initrng(iseed)
****  end RNG initialization

c      a = ammu*ammu + emulab*ame
c      b = emulab + ame
c      ac = sqrt(1.d0-ame*ame/ammu/ammu)
c      g = sqrt(emulab*emulab-ammu*ammu)*(ac)
c      em1 = a*b + abs(g)*sqrt((g*g-b*b)*ammu*ammu+a*a)
c      em1 = em1 / (b*b-g*g)
c      em2 = a*b - abs(g)*sqrt((g*g-b*b)*ammu*ammu+a*a)
c      em2 = em2 / (b*b-g*g)      
c      print*,em1,em2,em1-em2,sqrt((g*g-b*b)*ammu*ammu+a*a)
c      stop
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine explainvalues(setkey)
      character*(*) setkey(*)
      integer ucs4
      parameter (ucs4  = selected_char_KIND ('ISO_10646'))
      character(len=1,kind=ucs4) greekmu
! https://en.wikipedia.org/wiki/Box-drawing_character
      greekmu = char (int (z'2514'), ucs4) 

      print*,'Principal parameters:'
      print*,setkey(4) ,' ---> incoming muon charge in e+ charge'//
     .     ' units [-1/1]'
      print*,setkey(1)  ,' ---> nominal muon beam energy (GeV)'
      print*,setkey(2)  ,' ---> beam energy spread percentage (%)'
      print*,setkey(3)  ,' ---> minimum electron energy in the LAB'//
     .     ' (GeV)'
      print*,setkey(24) ,' ---> maximum electron angle in the LAB'//
     .     ' (mrad)'
      print*,setkey(22) ,' ---> minimum muon angle in the LAB (mrad)'
      print*,setkey(23) ,' ---> maximum muon angle in the LAB (mrad)'
      print*,setkey(29) ,' ---> if applying acoplanarity cut and'//
     .     ' to which value [yes/no mrad]'
      print*,setkey(31) ,' ---> if applying'//
     .     ' distance-from-elasticity-curve cut and'//
     .     ' to which value [yes/no mrad]'
      print*,setkey(11) ,' ---> to which order simulate events '//
     .     '[born/alpha/alpha2] for [LO/NLO/NNLO]'
      print*,setkey(21) ,' ---> if running of alpha must be used'//
     .     ' [off/on/hadr5/nsk/knt]'
      print*,setkey(25) ,' ---> if switching off hadronic VP [yes/no]'
      print*,setkey(6)  ,' ---> number of events to generate'
      print*,setkey(8)  ,' ---> if events have to be stored [yes/no]'
      print*,setkey(27) ,' ---> which mode to use to store events'//
     .     ' [0/1/2]'
      print*,'            ',greekmu,'> [0] plain ascii file '
      print*,'            ',greekmu,'> [1] root file (see README)'
      print*,'            ',greekmu,'> [2] on-the-fly xz'//
     .     ' compressed ascii file'
      print*,setkey(7)  ,' ---> path where to store outputs'
      print*,setkey(10) ,' ---> pseudo-RNG seed ("small" int)'
      print*,' '
      print*,'Internal tweaks:'
      print*,setkey(14)  ,' ---> if stored events are weighted'//
     .     ' or unweighted [weighted/unweighted]'
      print*,setkey(5)  ,' ---> "radiative" charges of'//
     .     ' the muon and electron leg, to switch on/off'//
     .     ' gauge invariants subsets [0/1 0/1]'
      print*,setkey(15)  ,' ---> minimum photon CM energy of the'//
     .     ' emitted photons in sqrt(s)/2 units. A.k.a. soft/hard'//
     .     ' separator'
      print*,setkey(26)  ,' ---> photon mass (IR regulator) (GeV)'
      print*,setkey(13)  ,' ---> drives maximum number of'//
     .     ' hard photons (> eps) and real pair emission [int]'
      print*,'            ',greekmu,'> [ < 0 ] maximum possible'//
     .     ' (up to 2)'
      print*,'            ',greekmu,'> [ 0 | 1 | 2 ] exactly this'//
     .     ' number of photons'
      print*,'            ',greekmu,'> [ 1000 ] both real electron and'
     .     //' muon pairs'
      print*,'            ',greekmu,'> [ 1001 ] only real electron'//
     .     ' pairs'
      print*,'            ',greekmu,'> [ 1002 ] only real muon pairs'
      print*,'            ',greekmu,'> [-1000 ] all possible final'//
     .     ' states'
      print*,setkey(16) ,' ---> files are dumped every nwrite'//
     .     ' generated events [int] (negative integer means files'//
     .' written every -'//trim(setkey(16))//' seconds)'
      print*,setkey(17) ,' ---> after nwarmup events, also unweighted'//
     .     ' generation is started [int]. Plays with ',
     .     trim(setkey(37)),'.'
      print*,setkey(28)  ,' ---> if writing distributions also at'//
     .     ' different orders [1/2/3]'
      print*,setkey(19)  ,' ---> maximum integrand for unweightening'
      print*,setkey(37)  ,' ---> typical integrated cross section'//
     .     ' within applied cuts, used for storage in ROOT format'
      print*,setkey(30)  ,' ---> syncronization mode for'//
     .     ' random numbers (see README) [0/1]'
      return
      end
**************************
      function r8tostr(a)
      character(*) r8tostr
      double precision a
      integer i
      logical is0      
      character*200 stringa

      i = precision(a)
      
      write(stringa,'(f30.18)')a
      stringa = adjustl(trim(stringa))
      is0 = .true.
      k = index(stringa,' ')-1
      do while(is0)
         if (.not.(stringa(k:k).eq.'0')) is0 = .false.
         k = k -1
      enddo
      r8tostr = stringa(1:k+1)
      return
      end
**************************
      function r4tostr(a)
      character(*) r4tostr
      real*4 a
      integer i
      logical is0      
      character*200 stringa
      write(stringa,'(f30.18)')a
      stringa = adjustl(trim(stringa))
      is0 = .true.
      k = index(stringa,' ')-1
      do while(is0)
         if (.not.(stringa(k:k).eq.'0')) is0 = .false.
         k = k -1
      enddo
      r4tostr = stringa(1:k+1)
      return
      end
**************************
      function itostr(a)
      character(*) itostr
      integer a
      integer i
      logical is0      
      character*200 stringa
      write(stringa,'(i30)')a
      stringa = adjustl(trim(stringa))
      itostr  = stringa
      return
      end
**************************
      function i8tostr(a)
      character(*) i8tostr
      integer*8 a
      integer i
      logical is0      
      character*200 stringa
      write(stringa,'(i30)')a
      stringa = adjustl(trim(stringa))
      i8tostr  = stringa
      return
      end

      
