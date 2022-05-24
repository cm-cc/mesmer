      subroutine userinterface
      use, intrinsic :: iso_fortran_env, only : stdin=>input_unit,
     .                                        stdout=>output_unit,
     .                                         stderr=>error_unit
      implicit double precision (a-b,d-h,o-z)
      implicit character*10 (c)
      logical istop
      character*50 r8tostr,r4tostr,itostr,i8tostr
      double precision csi(2)
      character*20  c,ctmp,ctmp2,ctmp3
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
      common/pionparameters/pi0m,fpi0,wpi0gg,gpi0
      double precision mz,mz2,gz,gz2,mw,mw2,stw,ctw
      common/zparameters/mz,mz2,gz,gz2,mw,mw2,stw,ctw,gv,ga
      common/qedORDER/ord
      integer*8 iwriteout
      integer isync,iseed,nsearch,iverbose,nw,istorver,icwvp
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync,
     .     istorver,icwvp,ivpfl
      common/iseedoninput/iseedext
      common/charinput/model,eventlimiter,store,storefile,outfile
      common/cmnpath/path
      common/realinput/anpoints,sdifmax

      logical areweighted
      common/weightedeventsornot/areweighted
      
      common/hadoff/ihadoff
      integer ififo
      common/ififocmn/ififo

      character*20 hstnm
      common/hostandpid/hstnm,idproc      
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
      common/mueexpsetup/emulab,eemin,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela

      double precision phmass
      common/photonmasslambda/phmass
      
      common/ialpharunning/iarun

!!!     should be initialized here! used for real pairs
      integer iper,noper,iid
      common/identicalandperipheral/iper,noper,iid
      data iper,iid /0,0/
!!!      
      common/weightnormalization/wnorm
      character*9 setkey(41)
c for teubner et al, KNT
      integer iteubn
      common/teubner/iteubn
c for Novosibirsk VP routine
      integer inskvp
      common/nskvp/inskvp

cc
      integer initerr
      common/init_error_code/initerr
      data initerr /0/
      
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
     +     'wnorm ',  ! 37
     +     'themin', ! 38
     +     'Ethr', ! 39
     +     'ththr',  ! 40
     +     'extmubeam'  ! 41
     +     /
* BSM parameters
      common/bsmparams/bsmmass,yae,yam,bsmlambda,gagg,ibsm
      data ibsm /0/
***
      character*20 programversion
      common/version/programversion
      character*200 welcome(23)

      character*4 bbf,ebf
      character*4 bit,eit
      integer ionlybubbles,ibubbleonv,ibubbleonboxsub,iexnnlobubbles
      common/ionlybubblescmn/iexnnlobubbles,ionlybubbles,ibubbleonv,
     .     ibubbleonboxsub
      double precision colch2(5)
      integer bubblemask(5)
      common/colorandch2factors/colch2,bubblemask
*****
******* for the header      
      common/tobewrittenintheheader/wnormalization,sdmax
****
      integer ismualongz
      common/ISmuondirection/ismualongz
***** 
      integer isstandalone
      character*500 datacard
      common/standalone/isstandalone,datacard ! set in main.F or init_mesmer() !!

      logical needswarmup,imwarmingup
      common/warmuppery/needswarmup,imwarmingup
      data needswarmup,imwarmingup /.false.,.false./
*******************************************
***   it decides from where to read the input, according to isstandalone.
***   isstandalone is set to 1 in main.F and to 0 in init_mesmer()
      inu = stdin
      if (isstandalone.eq.0) then
         inu = 33
ccc set in init_mesmer in wrapper.F
         open(inu,file=datacard,status='old',iostat=ioperr)
         if (ioperr.gt.0) then
            print*,' The file ',trim(datacard),' does not exists.'
            initerr = 1
            return
c            print*,' Stopping MESMER execution.'
c            stop
         endif
      endif
*******************************************
*** opening stdout as utf-8!      
      open(stdout,encoding='UTF-8')
***
      istorver = 2
      icwvp    = 1
***
      ismualongz = 1
c      ismualongz = 0

*** for real pairs
      iid   = 1
      iper  = 1
      noper = 1
*** for virtual pairs
      ionlybubbles    = 0
      ibubbleonv      = 0
      ibubbleonboxsub = 0
      iexnnlobubbles  = 0
      bubblemask = (/0,0,0,0,0/)
      
#ifndef RELEASE
      print*,'iid iper and noper'
      read(inu,'(I2,I2,I2)')iid,iper,noper
**********************      
*******************************
c      print*,'Input ionlybubbles'
c      read*,ionlybubbles
      print*,'Input iexnnlobubbles (0/1/2)'
!     0 -- no nnlo exact bubbles
!     1 -- only VP NLO blob
!     2 -- add also non factorizable corrections via DR relations
      read(inu,'(I2)')iexnnlobubbles
      print*,'only bubbles on verter and/or boxes'
      read(inu,'(I2,I2)')ibubbleonv,ibubbleonboxsub
      print*,'bubble mask?'
      read(inu,'(I2,I2,I2,I2,I2)')bubblemask
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
      programversion='MESMER 1.1'
#include "mesmerversion.f" 
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
      pi = 3.1415926535897932384626433832795029D0
**
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
***   end BSM parameters

      alpha = 1.d0/137.03599907430637d0

      ame  = 0.510998928d-3     ! electron mass
      ammu = 105.6583715d-3 ! muon mass

c      ammu  = ame
      
      taumass = 1.77682D0
      topmass = 175.6D0
      pi0m    = 0.1349766d0
      fpi0    = 0.092388d0 !!!!!!!!!!!!!!!
      wpi0gg  = alpha*alpha*pi0m**3/64.d0/pi**3/fpi0/fpi0
      gpi0    = 4.d0*pi*wpi0gg/pi0m**3
      gpi0    = sqrt(gpi0)
      
** for Z exchange (at LO)      
      mz = 91.1876d0
      gz = 2.4952d0
      mw = 80.398d0

c      ame = ammu
c      ammu = ame      
c      write(*,'(A,A)',advance='no')'test',achar(13)
c      stop
*

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
      themin  = 0.d0
      thmumax = 100d0
      themax  = 100d0
      Ethr    = 0.2d0
      ththr   = 100.d0
      
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
      
      nsearch      =  0       ! 1000000
      wnorm        = -1.d0
      sdifmax      = 1.d-10
      
      iverbose     = 0
      eventlimiter = 'w'  ! 'unw' or 'w'
      store        = 'no'
c      ififo = 2 ! for gzipped file on the fly
c      ififo = 1 ! for writing to a fifo, used with root interface      
c      ififo = 0 ! ascii file
c      ififo = 3 ! just write on the fifo
      ififo = 0
      
      arun         = 'on'

      iarun = 0
      if (arun.ne.'off') iarun = 1
***********
      istop = .false.
      iii   = 0
      ishow = isstandalone
      do while(.not.istop)
         IF (ishow.eq.1) THEN
         print*,' '
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
         ctmp2 = 'yes'
         if (ismualongz.eq.1) ctmp2 = 'no'
         write(stringa,*)' [ '//setkey(41)//' ] Feed generic'//
     .    ' muon beam externally (beam profile): ',trim(ctmp2)
         print*,trim(stringa)
         write(stringa,*)' [ '//setkey(3)//' ] Minimum electron'//
     .       ' LAB energy = ',trim(r8tostr(eemin)),' GeV'
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(38)//' ] Minimum electron'//
     .       ' LAB angle = ',trim(r8tostr(themin)),' mrad'
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

         write(stringa,*)' [ '//setkey(39)//' ] ''Detectability'''//
     .        ' energy threshold for extra leptons = ',
     .        trim(r8tostr(Ethr)),' GeV'
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(40)//' ] Maximum'//
     .        ' ''detectability'' LAB angle'//
     .        ' for extra leptons = ',
     .        trim(r8tostr(ththr)),' mrad'
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

         if (icwvp.eq.0) ctmp2 = 'no'
         if (icwvp.eq.1) ctmp2 = 'yes'
         if (istorver.eq.1) ctmp3 = 'old'
         if (istorver.eq.2) ctmp3 = 'new'
         write(stringa,*)' [ '//setkey(8)
     .        //' ] events storage: '//
     .        trim(store)//' '//trim(ctmp2)
c         //
c     .        ' '//trim(ctmp3)
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(27)
     .//' ] mode to store: ',trim(itostr(ififo))         
         print*,trim(stringa)

         write(stringa,*)' [ '//setkey(7)
     .//' ] files saved in '//trim(path)
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
      ENDIF
      do k = 1,200
         if (path(k:k).ne.' ') lpath=k
      enddo         
      if (path(lpath:lpath).ne.'/') then
         lpath = lpath+1
         path(lpath:lpath)='/'
      endif
      print*,' '
c 111     write(*,42)trim(bbf//'Insert "parameter value" or'
c     .        //' "run" or "quit": '//ebf)
c 42      format(A54,$)
 111     write(*,42,advance="no")
     .        trim(bbf//'Insert "parameter value" or'
     .        //' "run" or "quit": '//ebf)
 42      format(A54)
         read(inu,'(A)')cinput
c         print*,'Setting ',trim(cinput)
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
     .        c.eq.'extmubeam'.or.
     .        c.eq.'elastcut'.or.
     .        c.eq.'radchs') then
            read(cinput(k+1:),'(bn,A)') cx
            if (c.eq.'ord')     ord   = cx
            if (c.eq.'model')   model = cx
            if (c.eq.'path')    path  = cx
            if (c.eq.'store') then
c               read(cx,*)ctmp,ctmp3,ctmp2
               ctmp3 = 'new'
               read(cx,*)ctmp,ctmp2
               if (ctmp2.eq.'yes') icwvp = 1
               if (ctmp2.eq.'no')  icwvp = 0
c               if (ctmp3.eq.'new') istorver = 2
c               if (ctmp3.eq.'old') istorver = 1
               store   = ctmp
            endif
            if (c.eq.'extmubeam') then
               read(cx,*)ctmp
               if (ctmp.eq.'yes') ismualongz = 0
               if (ctmp.eq.'no')  ismualongz = 1
            endif
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
            if (c.eq.'themin')  themin = x
            if (c.eq.'themax')  themax = x
            if (c.eq.'Ethr')    Ethr = x
            if (c.eq.'ththr')   ththr = x
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
c     below the double \\ in front of \a is because otherwise the -fbackslash compile option
c     interprets \a as an alert
         if (inskvp.eq.1) then
c            print*,'                  '
            print*, 'Using Fedor Ignatov''s routine for '//
     .           'hadronic \Delta\\alpha.'
         else
            if (iteubn.lt.1) then
c               print*,'                  '
               if (ifredv.eq.19)
     .              print*, 'Using Fred Jegerlehner''s HADR5x19'//
     .              ' routine for hadronic \Delta\\alpha.'
               if (ifredv.eq.12)
     .              print*, 'Using Fred Jegerlehner''s HADR5n12'//
     .              ' routine for hadronic \Delta\\alpha.'
               if (ifredv.eq.17)
     .              print*, 'Using Fred Jegerlehner''s HADR5n17'//
     .              ' routine for hadronic \Delta\\alpha.'
            else
c               print*,'                  '
               print*,'Using Keshavarzi-Nomura-Teubner''s routine for'//
     .              ' hadronic \Delta\\alpha.'//
     .              ' [VP_KNT (version 3.0.1, Feb 16 2021)]'
            endif
         endif
      endif

      if (iarun.eq.0) then
         ivpfl = 0
      else
         if (iteubn.eq.0.and.inskvp.eq.0) then            
            if (ifredv.eq.19) ivpfl = 19
            if (ifredv.eq.12) ivpfl = 12
            if (ifredv.eq.17) ivpfl = 17            
         endif
         if (iteubn.eq.1) ivpfl = 20
         if (inskvp.eq.1) ivpfl = 30
      endif
***********************************************************************
      call execute_command_line('[ -d '//trim(path)//' ]',exitstat=ierr)      
      if (ierr.eq.0) then ! i.e. if the directory exists
         print*,' '
         print*,'WARNING:'
         print*,'The directory ',path(1:lpath),' already exists.'
         print*,'The files inside it are going to be overwritten...'
      endif
      close(11)
      call execute_command_line('mkdir -p '//trim(path)//' 2>/dev/null')

      corrections = ord
      if (ord.eq.'alpha')  corrections='oal'
      if (ord.eq.'alpha2') corrections='oal2'

cc
c      if (corrections.eq.'born') then
c         storefile = path(:lpath)//'events-LO.dat'
c      elseif (corrections(1:4).eq.'oal2') then
c         storefile = path(:lpath)//'events-NNLO.dat'
c      else
c         storefile = path(:lpath)//'events-NLO.dat'
c      endif
c useing always events.dat (widthout "order" info...)
      storefile = path(:lpath)//'events.dat'

      outfile  = path(:lpath)//
     .     'stat_'//corrections(1:len_trim(corrections))//'.txt'

****************************
      open(42,file=trim(path)//'iwriteout',status='unknown')
      write(42,*)iwriteout
      close(42)
****************************
***   checking if nsearch, sdifmax and wnorm are consistent
***   setting here needswarmup, imwarmingup
      print*,' '
      if (store.eq.'yes'.or.isstandalone.eq.0
     .     .or.eventlimiter.eq.'unw') then
         if (eventlimiter.eq.'w') then
            print*,'You are requesting non-constant weight events.'
            if (wnorm.le.0.d0) then
               if (nsearch.lt.1) then
                  print*,'You should set wnorm > 0 or nwarmup > 0.'
                  if (isstandalone.eq.1) then
                     print*,'Stopping.'
                     stop
                  else
                     initerr = 1
                     return
                  endif
               else
                  print*,'Number of events for warming-up:',nsearch
                  print*,'Will be used to calculate wnorm'
                  needswarmup = .true.
               endif
            else !if (wnorm.gt.0.d0) then
               print*,'Using wnorm:',wnorm
               print*,'and not warming-up'
            endif
         else !if (eventlimiter.eq.'unw') then
            print*,'You are requesting constant weight events.'
            if (nsearch.lt.1.and.sdifmax.lt.0.0001d0) then
               print*,'You should set sdmax > 0 or nwarmup > 0.'
               if (isstandalone.eq.1) then
                  print*,'Stopping.'
                  stop
               else
                  initerr = 1
                  return
               endif
            else
               if (sdifmax.lt.0.0001d0) then
                  print*,'Number of events for warming-up:',nsearch
                  needswarmup = .true.
               else
                  print*,'Using sdmax:',sdifmax
                  print*,'and not warming-up'
               endif
            endif
         endif
      else
         print*,'No need of warming-up'
         if (isstandalone.eq.0) then
            print*,'Running as a subroutine'
         else
            print*,'Running in standalone mode'
         endif
      endif
***********************************************************************
      
      if (eventlimiter.eq.'unw'.and.icwvp.eq.1.and.store.eq.'yes') then
         print*,' Storing unweighted events selected, resetting writing'
         print*,' coefficients for VP reweight to ''no'''
         icwvp = 0
      endif
      
      spremu = spremu * 0.01d0
      semu   = spremu * emulab 
      Eemin  = max(Eemin,ame)
      if (abs(iqmu).ne.1) then
         print*,'Wrong charge selected! stopping...'
         if (isstandalone.eq.1) then
            stop
         else
            initerr = 1
            return
         endif
      endif
      Qmu = iqmu * 1.d0

      if (QRe.lt.0.or.QRe.gt.1) then
         print*,'Wrong electron "radiative" charge.'
         print*,'Must be 0 or 1. Stopping'
         if (isstandalone.eq.1) then
            stop
         else
            initerr = 1
            return
         endif
      endif
      if (QRmu.lt.0.or.QRmu.gt.1) then
         print*,'Wrong muon "radiative" charge.'
         print*,'Must be 0 or 1. Stopping'
         if (isstandalone.eq.1) then
            stop
         else
            initerr = 1
            return
         endif
      endif

      if (chadoff.eq.'yes') ihadoff = 1
      if (chadoff.eq.'no')  ihadoff = 0

**** number of distributions in one shot
      if (nw.gt.0) then
         if (ord.eq.'born')  nw = 1
         if (ord.eq.'alpha') nw = min(2,nw)
         nw = min(nw,3)
      else
         print*,' Not writing distribution files!'
      endif
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
         if (isstandalone.eq.1) then
            stop
         else
            initerr = 1
            return
         endif
      endif
#endif
**************************************************

*********      
      if (isstandalone.eq.0) close(inu)
*********
      print*,'                  '
      print*,'Initializing '//trim(programversion)//'...'
      print*,'                  '

      idproc = getpid()
      call hostnm(hstnm)

***********************************************      
      areweighted = .true.
      if (eventlimiter.eq.'unw') areweighted = .false.
***********************************************      

      
**** further initializations
      s = ammu**2 + ame**2 + 2.d0*ame*emulab
      call initmuescattering(s)      
****
      call init_rnd_sequence    ! uses iseedext and sets iseed

      wnormalization = wnorm
      sdmax          = 1.1d0*sdifmax
      
      if (abs(ammu-ame).lt.1d-3) then
         do jj = 1,1000
            print*,'ammu ~ ame!!!!'
         enddo
      endif
      
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
      print*,setkey(1) ,' ---> nominal muon beam energy (GeV)'
      print*,setkey(2) ,' ---> beam energy spread percentage (%)'
      print*,setkey(41),' ---> if muon beam 3-momentum is externally'//
     .     ' fed (beam profile) [yes/no].'
      print*,'               If ''yes'', the beam 3-momentum must'//
     .     ' be fed into the pipe ''path''/beamprofile.fifo'
      print*,setkey(3) ,' ---> minimum electron energy in the LAB'//
     .     ' (GeV)'
      print*,setkey(38) ,' ---> minimum electron angle in the LAB'//
     .     ' (mrad)'      
      print*,setkey(24) ,' ---> maximum electron angle in the LAB'//
     .     ' (mrad)'
      print*,setkey(22) ,' ---> minimum muon angle in the LAB (mrad)'
      print*,setkey(23) ,' ---> maximum muon angle in the LAB (mrad)'
      print*,setkey(29) ,' ---> if applying acoplanarity cut and'//
     .     ' to which value [yes/no mrad]'
      print*,setkey(31) ,' ---> if applying'//
     .     ' distance-from-elasticity-curve cut and'//
     .     ' to which value [yes/no mrad]'
      print*,setkey(39) ,' ---> ''detectability'' threshold energy'//
     .     ' for extra leptons (GeV)'
      print*,setkey(40) ,' ---> maximum ''detectability'' angle for'//
     .     ' extra leptons  (mrad)'      
      print*,setkey(11) ,' ---> to which order simulate events '//
     .     '[born/alpha/alpha2] for [LO/NLO/NNLO]'
      print*,setkey(21) ,' ---> if running of alpha must be used'//
     .     ' [off/on/hadr5/nsk/knt]'
      print*,setkey(25) ,' ---> if switching off hadronic VP [yes/no]'
      print*,setkey(6)  ,' ---> number of events to generate'
      print*,setkey(8)  ,' ---> (2 parameters) if events have to be'//
     .     ' stored [yes/no],'
c      print*,'               which format [new/old] (for v2/v1),'
      print*,'               if writing also'//
     .     ' coefficients for VP reweighting [yes/no]'
c      //
c     .     ' (active only for v2)'
      print*,setkey(27) ,' ---> which mode to use to store events'//
     .     ' [0/1/2/3]'
      print*,'            ',greekmu,'> [0] plain ascii file '
      print*,'            ',greekmu,'> [1] root file (see README)'
      print*,'            ',greekmu,'> [2] on-the-fly xz'//
     .     ' compressed ascii file'
      print*,'            ',greekmu,'> [3] just write to the fifo'//
     .     ' and wait for an external process to read events from it'
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
     .     ' states (except 1010)'
      print*,'            ',greekmu,'> [ 1010 ] pi0 production'//
     .     ' (and decay in to two photons)'
      print*,setkey(16) ,' ---> files are dumped every nwrite'//
     .     ' generated events [int] (negative integer means files'//
     .' written every -'//trim(setkey(16))//' seconds)'
      print*,setkey(17) ,' ---> after nwarmup events, also unweighted'//
     .     ' generation is started [int]. Plays with ',
     .     trim(setkey(37)),'.'
      print*,setkey(28)  ,' ---> if writing distributions also at'//
     .     ' different orders [0/1/2/3]'//
     .     ' (0 disables writing distributions)'
      print*,setkey(19)  ,' ---> maximum integrand for unweightening'
      print*,setkey(37)  ,' ---> typical integrated cross section'//
     .     ' within applied cuts, used for storage in ROOT format'
      print*,setkey(30)  ,' ---> syncronization mode for'//
     .     ' random numbers (see README) [0/1]'
      return
      end
**************************
      subroutine printoutput
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
      common/iseedoninput/iseedext
      common/qedORDER/ord
      common/charinput/model,eventlimiter,store,storefile,outfile
      character*200 path
      common/cmnpath/path
      common/realinput/anpoints,sdifmax
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
     .     hitpmiss,istopsearch,nneg,sumover,sum2over,fmax,
     .     ntomicro,k,npoints,naccepted,kt,itime,rlxstate,nphmax,nover, 
     .     in_conf_spin

      common/printoutcommonblock/kwritelast
      data kwritelast /0/
      
      open(42,file=trim(path)//'iwriteout',status='unknown')
      read(42,*)iwriteout
      close(42)

      greekmu = ucs4_'\u03BC'
      greekmu = char(int (z'03BC'), ucs4)
      
      kwritenow = k
      idk   = kwritenow-kwritelast
      kwritelast = kwritenow
      itime2 = time8()
      call cpu_time(endtime)
            
      call getrndstatus(rlxstate)
      
cc in common      xsec = sum/k
cc in common      var  = sqrt((abs(sum2/k-xsec*xsec))/k)

      open(10,file=outfile,status='unknown',encoding='UTF-8')
      
      write(10,*)'RNG status at this point:'
      write(10,*)rlxstate
      write(10,*)' '
      
      write(10,*)' Process PID:',idproc,' running on ',hstnm
      
      dtime     = endtime-starttime
      dtimetrue = endtime-truestarttime
      starttime = endtime
      write(10,*)' file written last time ~',itime2-itime,'s ago'
      
      write(10,*)' CPU time to generate last         ',idk,
     .     'events',dtime,'s (',dtime/idk*1d6,
     .     greekmu,'s/ev)' 
      write(10,*)' CPU time to generate     ',k,'events',
     .     dtimetrue,'s (',dtimetrue/k*1d6,greekmu,'s/ev)' 
      itime = itime2
      write(10,*)' '
      write(10,*)'Inputs & parameters for this run: '
      iqmu = qmu
      write(10,'(1x,A,I4)')'Incoming muon charge',iqmu
      write(10,'(1x,A,f12.4,A)')
     .     'nominal muon beam energy =',emulab,' GeV'
      write(10,'(1x,A,f12.4,A)')
     .     'muon beam energy  spread =',semu/emulab*100d0,' %'
      write(10,'(1x,A,f12.4,A)')
     .     'sqrt(s)   =',ecmsnom,' GeV'
      write(10,'(1x,A,f12.4,A)')
     .     'minimum electron LAB energy =',eemin,' GeV'
      write(10,'(1x,A,f12.4,A)')
     .     'maximum electron LAB angle  =',themax,' mrad'
      write(10,'(1x,A,f12.4,A)')
     .     'minimum muon LAB angle      =',thmumin,' mrad'
      write(10,'(1x,A,f12.4,A)')
     .     'maximum muon LAB angle      =',thmumax,' mrad'
      write(10,'(1x,A,i2,1x,f12.4,A)')
     .     'applying acopl. cut and value =',ina7,dthna7max,
     .     ' mrad'
      write(10,'(1x,A)')
     .     'ord    = '//ord
      write(10,'(1x,A)')
     .     'model  = '//model
      write(10,'(1x,A,i5)')
     .     'nphot mode =',nphotmode
      write(10,'(1x,A,i9,i9)')
     .     'seeds   =',iseedext,iseed
      write(10,'(1x,A,i5)')
     .     'alpha running  =',iarun
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
     .     'eps    = ',eps
      write(10,*)' '
      write(10,*)'Weights normalization ',wnorm,' (',greekmu,'b)'
      write(10,*)'fmax for unweightening ',fmax
      write(10,*)'true sdifmax ',sdiftruemax
      
      write(10,*)' '
      if (eventlimiter.eq.'w') then
         write(10,'(A,i12,A)')'~ Generating ',npoints,
     .        ' weighted events ~'
      else
         write(10,'(A,i12,A)')'~ Generating ',npoints,
     .        ' unweighted events ~'
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
     :        iii,'  photons: ',xsecp(i),' +-',varp(i),
     :        '     (',fractions(i),' %)'
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
     :        ' real pairs: ',xsecpairs,' +-',varpairs,
     :        '     (',fracpairs,' %)'
         write(10,'(A,f17.11,A,f22.11)')
     :        '   [e+e-   pairs: ',xsecpairs1,' +-',varpairs1
         write(10,'(A,f17.11,A,f22.11,A)')
     :        '    mu+mu- pairs: ',xsecpairs2,' +-',varpairs2,']'
      endif
      write(10,'(1x,A,f22.11,A,f22.11,A,A,A)')
     :     'total:      ',xsec,' +-',var,' ',greekmu,'b'
      
      write(10,*)' '
      write(10,'(1x,A,i12,A,i20)')
     :     'n. ',k,' of ',npoints
      eff = (1.d0*naccepted)/k
      write(10,'(1x,A,f12.8,A)')
     :     'cut points ',abs(1.d0 - eff)*100d0,' %'
      write(10,*)'::::::>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<::::::'
      write(10,*)' '
      write(10,*)'::::::>>>>>> unweighted events <<<<<<::::::'
      hmeff = 0.d0
      if (hitpmiss.gt.0.d0) hmeff = hit/hitpmiss
      write(10,'(1x,A,f10.6,A)')
     :     'hit or miss efficiency ',hmeff*100,' %'
      write(10,'(1x,A,f12.0)')
     :     'unweighted events generated ',hit
      biashit        = 0.d0
      biashitpmiss   = 0.d0
      biasneghit     = 0.d0
      biasneghitmiss = 0.d0
      sezover        = 0.d0
      errsezover     = 0.d0
      sezneg         = 0.d0
      errsezneg      = 0.d0
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
     .     'total (',greekmu,'b):   ',hmxsect,' +-',hmerr,' +'
      write(10,'(1x,A,f22.11,A,f22.11,A)')
     .     '       !!!!   ',sezover,' +-',
     .     errsezover,' (bias over fmax) +'
      write(10,'(1x,A,f22.11,A,f22.11,A)')
     .     '       !!!!   ',-sezneg,' +-',
     .     errsezneg,' (bias negative)'
      write(10,'(1x,A,f22.11,A,f22.11)')
     .     'total + biases: ',hmxsect
     .     +sezover-sezneg,' +-',hmerr+errsezover+errsezneg
      write(10,*)' '
      write(10,'(1x,A,i12)')
     .     'N. points with w > fmax (bias): ',nover
      write(10,'(1x,A,f11.7,A,f11.7,A)')
     .     'bias/hit and bias/(hit+missed): ',
     .     biashit*100,' % and ',biashitpmiss*100,' %'
      write(10,'(1x,A,i12)')
     .     'N. points with w < 0:',nneg
      write(10,'(1x,A,f10.5,A,f10.5,A)')
     .     'biases for w < 0:   ',
     .     biasneghit*100,' % and ',biasneghitmiss*100,' %'
      write(10,*)' '
      write(10,'(1x,A,f12.4,A,f7.4,A)')
     .     'Used fmax: ',fmax,' = ',fmax/sdifmax,' * sdmax'
      write(10,*)
     .     'True upper/lower limits (sdmax sdmin):'
      write(10,*)sdiftruemax,sdiftruemin
      
      nwhenprint = nwhenmax
c     write(10,'(1x,A,i3,A)')
c     :           'when there were ',nwhenprint,' photons'
      write(10,*)'::::::>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<::::::'
      write(10,*)' '
      write(10,*)'reweighted cross sections'
      xsw1 = sumw1/k
      write(10,*)'W1  ',xsw1,'+-',sqrt((abs(sumw12/k-xsw1**2))/k)
      xsw2 = sumw2/k
      write(10,*)'W2  ',xsw2,'+-',sqrt((abs(sumw22/k-xsw2**2))/k)
      xsLO = sumLO/k
      write(10,*)'LO  ',xsLO,'+-',sqrt((abs(sumLO2/k-xsLO**2))/k)
      xsNLO = sumNLO/k
      write(10,*)'NLO ',xsNLO,'+-',
     .     sqrt((abs(sumNLO2/k-xsNLO**2))/k)
      close(10)

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
      write(stringa,'(f30.16)')a
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
      write(stringa,'(f30.16)')a
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
