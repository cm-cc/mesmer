*** storage routines, storing to a simple ASCII file
      subroutine initstorage(storfile,wnorm,sdifmax)
      implicit double precision (a-h,o-z)
      character*(*) storfile
      common/storageunit/istorunit,ilens
      common/ififocmn/ififo
      common/tobewritten/wnormalization,sdmax

      integer g,mum,mup,em,ep
      common/LHAparticlesids/g,mum,mup,em,ep
      data g,mum,mup,em,ep/22,13,-13,11,-11/

      print*,'Initializing storage...'
      
      wnormalization = wnorm
      sdmax          = sdifmax
      
      istorunit = 61
c      ififo = 2 ! for gzipped file on the fly
c      ififo = 1 ! for writing to a fifo, used with root interface      
c     ififo = 0 ! ascii file

c      print*,ififo
c      stop

      ilens = lnblnk(storfile)-4
         
      if (ififo.eq.2) then  
         call
     >        EXECUTE_COMMAND_LINE("mkfifo "
     >        //storfile(1:ilens)//".fifo 2>/dev/null")

c XZ         
         call EXECUTE_COMMAND_LINE("xz -9 < "//storfile(1:ilens)//
     >        ".fifo > "//storfile(1:lnblnk(storfile))//".xz",
     >        wait=.false.)
c BZIP2         
c         call EXECUTE_COMMAND_LINE("bzip2 -9 < "//storfile(1:ilens)//
c     >        ".fifo > "//storfile(1:lnblnk(storfile))//".bz2",
c     >        wait=.false.)
c ZSTD
c         call EXECUTE_COMMAND_LINE("zstd -9 -T0 < "//storfile(1:ilens)//
c     >        ".fifo > "//storfile(1:lnblnk(storfile))//".zstd",
c     >        wait=.false.)

         open
     >     (istorunit,file=storfile(1:ilens)//".fifo",status='unknown')

      elseif (ififo.eq.1) then
         call
     >    EXECUTE_COMMAND_LINE("mkfifo "//storfile(1:ilens)//".fifo")         
c         call EXECUTE_COMMAND_LINE("./write_MuE_MCevents.exe "
         call EXECUTE_COMMAND_LINE("./write-root-events "
     >        //storfile(1:ilens)//".fifo "
c     >        //storfile(1:ilens)//".footer "
     >        //storfile(1:ilens)//".root",
     >        wait=.false.)
         open
     >     (istorunit,file=storfile(1:ilens)//".fifo",status='unknown')
      else
         open(istorunit,file=storfile,status='unknown')
      endif

      print*,'...and now starting to store...'
      
      
      return
      end
*********************************************************
      subroutine eventstorage_v1(w,nev,ie,p3,p4,p5,p6,qph)
      implicit double precision (a-h,o-z)
! already in invariats.h      dimension p3(0:3),p4(0:3)
      integer*8 nev,nstored,kbefore,nevts
      dimension qph(40,0:3),qph1(0:3),qph2(0:3)
      dimension pin1(0:3),pin2(0:3),p5(0:3),p6(0:3)
      common/momentainitial/pin1,pin2
      common/storageunit/istorunit,ilens
      common/storedevents/nstored,nevts,ifirst
      data nstored,nevts,ifirst /0,0,0/

      common/parameters/ame,ammu,convfac,alpha,pi     
      common/mueexpsetup/emulab,eemin,semu,thmumin,themax,thmumax,
     .     dthna7max,cutela,ina7,iela
      integer*8 iwriteout
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync

      parameter (ikind = 8)
      double precision p3lab(0:3),p4lab(0:3)
      double precision p5lab(0:3),p6lab(0:3)      
      double precision qphlab(40,0:3)

      common/eventcount/kbefore
      data kbefore /0/
      
      character*20 programversion
      common/version/programversion
      character*20 hstnm
      common/hostandpid/hstnm,idproc

      common/iseedoninput/iseedext
      character*8 date
      character*10time
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

      common/weiunwei/iweight

      common/tobewritten/wnormalization,sdmax
      
      
      double precision p0(0:3)
      data p0 /0.d0,0.d0,0.d0,0.d0/
      
*     included only to access the common 'weights'      
#include "invariants.h"

      if (ifirst.eq.0) then
         write(istorunit,*)'<header>'
         write(istorunit,*)'  Program version: ',programversion
         write(istorunit,*)'  Process ID',idproc,'running on ',hstnm
         call date_and_time(date,time,zone,values)
         write(istorunit,'(A,I2,A,I2,A,I4,A,I2,A,I2,A,I2)')
     >    '   Started on ',values(3),'/',values(2),'/',
     >        values(1),' at ',values(5),':',values(6),':',values(7)
         if (eventlimiter.eq.'w') then
            write(istorunit,*)
     >  '  Storing ',int(anpoints),' weighted (non-zero weight) events'
            iweight = 1
         else
            write(istorunit,*)
     >           '  Storing ',int(anpoints),' unweighted events'
            iweight = 0
         endif
         write(istorunit,*)
     >           '  Unweighted events (1 = TRUE, 0 = FALSE)'
         if (iweight.eq.0) then
            write(istorunit,*)1
         else
            write(istorunit,*)0
         endif
         if (ord.eq.'born') then
            write(istorunit,*)'  Running mode: LO'
         elseif (ord.eq.'alpha') then
            write(istorunit,*)'  Running mode: NLO'
         elseif (ord.eq.'alpha2') then
            write(istorunit,*)'  Running mode: NNLO'
         else
            print*,'  Running mode not valid. Stopping in eventstorage'
            stop
         endif         
         write(istorunit,*)'  External and internal seed for RNG'
         write(istorunit,*)iseedext,iseed
         write(istorunit,*)'  ------------------------------'
         write(istorunit,*)'  Nominal muon beam energy (GeV)'
         write(istorunit,*)emulab
         write(istorunit,*)'  RMS for gaussian energy spread (GeV)'
         write(istorunit,*)semu
         write(istorunit,*)'  Muon charge'
         write(istorunit,*)Qmu
         write(istorunit,*)'  Muon mass (GeV)'
         write(istorunit,*)' ',ammu
         write(istorunit,*)'  Electron mass (GeV)'
         write(istorunit,*)' ',ame
         write(istorunit,*)'  1/alpha(0)'
         write(istorunit,*)' ',1.d0/alpha
         write(istorunit,*)'  CoM photon k_0 (GeV)',egmin
         write(istorunit,*)'  MC generation cut: min E_e (GeV)',eemin
         write(istorunit,*)'  ------------------------------'
         write(istorunit,*)'  Normalization cross section (mub)'
         write(istorunit,*)' ',wnormalization
         write(istorunit,*)'  Maximum weight [guessed]'
         write(istorunit,*)' ',sdmax/wnorm
         write(istorunit,*)'</header>'
         ifirst = 1
      endif
      
      nstored = nstored + 1
      nevts   = nevts + nev - kbefore
      ng = 0                    ! number of photons in the record

      write(istorunit,*)'<event>'
      if (ie.ne.0) then
         write(istorunit,*)iseedext ! RunNr
         write(istorunit,*)nstored  ! EventNr
         write(istorunit,*)0 ! number of final state particles
! p1lab already in common  call new_boost(p12lab,pin1,p1lab,1) 
         write(istorunit,*)0.d0,0.d0,0.d0  ! weights
         write(istorunit,*)0.d0            ! weight for LO
         
         write(istorunit,*)real(p1lab(0),ikind)
         write(istorunit,*)real(p0,ikind)
         write(istorunit,*)real(p0,ikind)    
      else
ccc   counting photons 
         do k = 1,40
            if (qph(k,0).gt.0.d0) then
               ng = ng +1
               call new_boost(p12lab,qph(k,0:3),qphlab(k,0:3),-1)
            endif
         enddo
         
c p1lab already in common         call new_boost(pin2,pin1,p1lab,1)
         call new_boost(p12lab,p3,p3lab,-1)
         call new_boost(p12lab,p4,p4lab,-1)

         write(istorunit,*)iseedext ! RunNr
         write(istorunit,*)nstored  ! EventNr

         npart = 2+ng
         ipair = 0
         if (p5(0).gt.0.d0.and.p6(0).gt.0.d0) then
            ipair = 1
            npart = npart + 2
            call new_boost(p12lab,p5,p5lab,-1)
            call new_boost(p12lab,p6,p6lab,-1)
         endif
         
         write(istorunit,*)npart ! number of final state particles

         kbefore = nev
         write(istorunit,*)real(w,ikind),
     >        real(w*extraweights(1),ikind),
     >        real(w*extraweights(2),ikind)
         write(istorunit,*)w*reweightLO  ! weight for LO         
         write(istorunit,*)real(p1lab(0),ikind)

         write(istorunit,*)real(p3lab,ikind)
         write(istorunit,*)real(p4lab,ikind)
         if (ng.gt.0) then
            do k = 1, ng
               write(istorunit,*)real(qphlab(k,:),ikind)
            enddo
         endif
         if (ipair.eq.1) then
            write(istorunit,*)real(p5lab,ikind)
            write(istorunit,*)real(p6lab,ikind)
         endif
      endif
      write(istorunit,*)'</event>'
ccc      print*,nevts
      return
      end
*********************************************************
      subroutine eventstorage(w,nev,ie,p3,p4,p5,p6,qph)
      implicit double precision (a-h,o-z)
! already in invariats.h      dimension p3(0:3),p4(0:3)
      integer*8 nev,nstored,nevts
      dimension qph(40,0:3),qph1(0:3),qph2(0:3)
      dimension pin1(0:3),pin2(0:3),p5(0:3),p6(0:3)
      common/momentainitial/pin1,pin2
      common/storageunit/istorunit,ilens
      common/storedevents/nstored,nevts,ifirst
      data nstored,nevts,ifirst /0,0,0/

      parameter (ikind = 8)
      double precision p3lab(0:3),p4lab(0:3)
      double precision qphlab(40,0:3)
* included only to access the common 'weights'      
#include "invariants.h"      
c      call eventstorage_newversion(w,nev,ie,p3,p4,p5,p6,qph)
      call eventstorage_v1(w,nev,ie,p3,p4,p5,p6,qph)
      return
      end
*******************************************
      subroutine finalizestorage(storfile)
      implicit double precision (a-h,o-z)
      character*(*) storfile
      integer istorunit,ilens
      logical itisthere,iopen
      common/storageunit/istorunit,ilens
      integer*8 foonhpm,foonh
      common/forfooter/xsw,exsw,xsunw,exsunw,
     .     xsbias,exsbias,xsbiasn,exsbiasn,truemax,sum,sum2,sumneg,
     .     sum2neg,foonhpm,foonh,nabove,nlt0
      common/weightnormalization/wnorm      
      ifoo = istorunit
      write(ifoo,*)'<footer>'
      write(ifoo,*)'(sub-)sample weighted XS (mub)'
      write(ifoo,*)xsw,' +- ',exsw
c      write(ifoo,*)'(sub-)sample unweighted XS (mub)'
c      write(ifoo,*)xsunw,' +- ',exsunw
      write(ifoo,*)'total shut points (weighted events)'
      write(ifoo,*)foonhpm
c      write(ifoo,*)'total hit points (unweighted events)'
c      write(ifoo,*)foonh
      write(ifoo,*)'true wmax'
      write(ifoo,*)truemax
      write(ifoo,*)'N. points above guessed wmax'
      write(ifoo,*)nabove
      write(ifoo,*)'N. points with negative weight'
      write(ifoo,*)nlt0
      write(ifoo,*)'Estimated XS bias for w > wmax (mub)'
      write(ifoo,*)xsbias,' +- ',exsbias
      write(ifoo,*)'Estimated XS bias for w < 0 (mub)'
      write(ifoo,*)xsbiasn,' +- ',exsbiasn
      write(ifoo,*)'Sum of weights and weights^2'
      write(ifoo,*)sum/wnorm,sum2/wnorm/wnorm
      write(ifoo,*)'Sum of weights < 0 and (weights < 0)^2'
      write(ifoo,*)sumneg/wnorm,sum2neg/wnorm/wnorm
      write(ifoo,*)'</footer>'      
***************************
      close(istorunit)
      end
*************************************************
*********************************************************
      subroutine eventstorage_newversion(w,nev,ie,p3,p4,p5,p6,qph)
! after agreeing new format with Giovanni (5/2/21)
      implicit double precision (a-h,o-z)
! already in invariats.h      dimension p3(0:3),p4(0:3)
      integer*8 nev,nstored,kbefore,nevts
      dimension qph(40,0:3),qph1(0:3),qph2(0:3)
      dimension pin1(0:3),pin2(0:3),p5(0:3),p6(0:3)
      common/momentainitial/pin1,pin2
      common/storageunit/istorunit,ilens
      common/storedevents/nstored,nevts,ifirst
      data nstored,nevts,ifirst /0,0,0/

      common/parameters/ame,ammu,convfac,alpha,pi     
      common/mueexpsetup/emulab,eemin,semu,thmumin,themax,thmumax,
     .     dthna7max,cutela,ina7,iela
      integer*8 iwriteout
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync

      parameter (ikind = 8)
      double precision p3lab(0:3),p4lab(0:3)
      double precision p5lab(0:3),p6lab(0:3)
      double precision qphlab(40,0:3)

      integer g,mum,mup,em,ep
      common/LHAparticlesids/g,mum,mup,em,ep

      integer id1,id2,id3,id4,id5,id6
      
      common/eventcount/kbefore
      data kbefore /0/
      
      character*20 programversion
      common/version/programversion
      character*20 hstnm
      common/hostandpid/hstnm,idproc

      common/iseedoninput/iseedext
      character*8 date
      character*10time
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

      common/weiunwei/iweight

      common/tobewritten/wnormalization,sdmax
      
      double precision p0(0:3)
      data p0 /0.d0,0.d0,0.d0,0.d0/
      
*     included only to access the common 'weights'      
#include "invariants.h"

      if (ifirst.eq.0) then
         write(istorunit,*)'<header>'
         write(istorunit,*)'  Program version: ',programversion
         write(istorunit,*)'  Process ID',idproc,'running on ',hstnm
         call date_and_time(date,time,zone,values)
         write(istorunit,'(A,I2,A,I2,A,I4,A,I2,A,I2,A,I2)')
     >    '   Started on ',values(3),'/',values(2),'/',
     >        values(1),' at ',values(5),':',values(6),':',values(7)
         if (eventlimiter.eq.'w') then
            write(istorunit,*)
     >  '  Storing ',int(anpoints),' weighted (non-zero weight) events'
            iweight = 1
         else
            write(istorunit,*)
     >           '  Storing ',int(anpoints),' unweighted events'
            iweight = 0
         endif
         write(istorunit,*)
     >           '  Unweighted events (1 = TRUE, 0 = FALSE)'
         if (iweight.eq.0) then
            write(istorunit,*)1
         else
            write(istorunit,*)0
         endif
         if (ord.eq.'born') then
            write(istorunit,*)'  Running mode: LO'
         elseif (ord.eq.'alpha') then
            write(istorunit,*)'  Running mode: NLO'
         elseif (ord.eq.'alpha2') then
            write(istorunit,*)'  Running mode: NNLO'
         else
            print*,'  Running mode not valid. Stopping in eventstorage'
            stop
         endif         
         write(istorunit,*)'  External and internal seed for RNG'
         write(istorunit,*)iseedext,iseed
         write(istorunit,*)'  ------------------------------'
         write(istorunit,*)'  Nominal muon beam energy (GeV)'
         write(istorunit,*)emulab
         write(istorunit,*)'  RMS for gaussian energy spread (GeV)'
         write(istorunit,*)semu
         write(istorunit,*)'  Muon charge'
         write(istorunit,*)Qmu
         write(istorunit,*)'  Muon mass (GeV)'
         write(istorunit,*)' ',ammu
         write(istorunit,*)'  Electron mass (GeV)'
         write(istorunit,*)' ',ame
         write(istorunit,*)'  1/alpha(0)'
         write(istorunit,*)' ',1.d0/alpha
         write(istorunit,*)'  CoM photon k_0 (GeV)',egmin
         write(istorunit,*)'  MC generation cut: min E_e (GeV)',eemin
         write(istorunit,*)'  ------------------------------'
         write(istorunit,*)'  Normalization cross section (mub)'
         write(istorunit,*)' ',wnormalization
         write(istorunit,*)'  Maximum weight [guessed]'
         write(istorunit,*)' ',sdmax/wnorm
         write(istorunit,*)'</header>'
         ifirst = 1
      endif
      
      nstored = nstored + 1
      nevts   = nevts + nev - kbefore
      ng = 0      ! number of photons in the record
      npairs = 0  ! number of pairs
      
      id2 = em
      id4 = id2
      id1 = mup
      if (Qmu.lt.0.d0) id1 = mum
      id3 = id1

      id5 = ep
      id6 = em
      
      write(istorunit,*)'<event>'
      if (ie.ne.0) then
c     write(istorunit,*)0
         write(istorunit,*)iseedext ! RunNr
         write(istorunit,*)nstored      ! EventNr
         write(istorunit,*)0 ! number of final state particles
! p1lab already in common         call new_boost(pin2,pin1,p1lab,1)
         write(istorunit,*)0.d0,0.d0,0.d0  ! weights
         write(istorunit,*)0.d0,0.d0       ! weight for LO and NLO
         
         write(istorunit,*)real(p1lab(0),ikind)
         write(istorunit,*)id3,real(p0(1:3),ikind)
         write(istorunit,*)id4,real(p0(1:3),ikind)    
      else
ccc   counting photons 
         do k = 1,40
            if (qph(k,0).gt.0.d0) then
               ng = ng +1
               call new_boost(p12lab,qph(k,0:3),qphlab(k,0:3),-1)
            endif
         enddo
ccc   counting pairs
         if (p5(0).gt.0.d0) then            
            npairs = npairs + 1
            call new_boost(p12lab,p5,p5lab,-1)
            call new_boost(p12lab,p6,p6lab,-1)

            if ((abs(dot(p5,p5)-ammu**2)/ammu/ammu).lt.1d-4) then
               id5 = mup
               id6 = mum
            endif
         endif
         
c p1lab already in common         call new_boost(pin2,pin1,p1lab,1)
         call new_boost(p12lab,p3,p3lab,-1)
         call new_boost(p12lab,p4,p4lab,-1)

         write(istorunit,*)iseedext ! RunNr
         write(istorunit,*)nstored  ! EventNr
            
         write(istorunit,*)2+ng  + 2*npairs! number of final state particles

         kbefore = nev
         write(istorunit,*)real(w,ikind),
     >        real(w*extraweights(1),ikind),
     >        real(w*extraweights(2),ikind)
         write(istorunit,*)w*reweightLO,w*reweightNLO  ! weight for LO and NLO
         write(istorunit,*)real(p1lab(0),ikind)

         write(istorunit,*)id3,real(p3lab(1:3),ikind)
         write(istorunit,*)id4,real(p4lab(1:3),ikind)
         
         if (npairs.gt.0) then
            write(istorunit,*)id5,real(p5lab(1:3),ikind)
            write(istorunit,*)id6,real(p6lab(1:3),ikind)
         endif

         if (ng.gt.0) then
            do k = 1, ng
               write(istorunit,*)g,real(qphlab(k,1:3),ikind)
            enddo
         endif
      endif
      write(istorunit,*)'</event>'
ccc      print*,nevts
      return
      end
*********************************************************
*** end of storage routines
*********************************************************
      subroutine twomostenergetic(ng,q,q1,q2)
      implicit double precision (a-h,o-z)
      integer iloc(1)
      dimension q(1:ng,0:3),q1(0:3),q2(0:3),ql(1:ng,0:3)
      q1 = 0.d0
      q2 = 0.d0
      ql = q      
      iloc = maxloc(ql(:,0))
      q1 = ql(iloc(1),:)
      ql(iloc(1),:) = 0.d0      
      iloc = maxloc(ql(:,0))
      q2 = ql(iloc(1),:)
      return
      end
*********************************************************      
      subroutine sortmomenta(ng,qm)
      implicit double precision (a-h,o-z)
      dimension pa(0:3),pb(0:3),qm(40,0:3),q1(0:3),q2(0:3),qmord(40,0:3)
      if (ng.eq.0) return
      ngl = ng
      call orderinverseqm(ngl,qm,qmord)
! in qmord photons are from the less energetic (k=1) to the most energetic 
! (k=ng). In qm I want to sort from the most to the less energetic...
      do k = 1,ng
         qm(ng+1-k,0) = qmord(k,0)
         qm(ng+1-k,1) = qmord(k,1)
         qm(ng+1-k,2) = qmord(k,2)
         qm(ng+1-k,3) = qmord(k,3)
      enddo
      return
      end
***
****************************************************
      subroutine mixandsortmomenta(ng,pa,pb,qm)
      implicit double precision (a-h,o-z)
      dimension pa(0:3),pb(0:3),qm(40,0:3),q1(0:3),q2(0:3),qmord(40,0:3)
      common/beforesort/p1o(0:3),p2o(0:3),qpho(40,0:3)
      if (ng.eq.0) return
      do k = 0,3
         p1o(k) = pa(k)
         p2o(k) = pb(k)
         do i = 1,ng
            qpho(i,k) = qm(i,k)
         enddo
      enddo
      ngl = ng + 2
      do k = 0,3
         qm(ng+1,k) = pa(k)
         qm(ng+2,k) = pb(k)
      enddo
      call orderinverseqm(ngl,qm,qmord)
      do k = 0,3
         pa(k)= qmord(ng+2,k)
         pb(k)= qmord(ng+1,k)
         qm(ng+1,k) = 0.d0
         qm(ng+2,k) = 0.d0
      enddo
! in qmord photons are from the less energetic (k=1) to the most energetic 
! (k=ng). In qm I want to sort from the most to the less energetic...
      do k = 1,ng
         qm(ng+1-k,0) = qmord(k,0)
         qm(ng+1-k,1) = qmord(k,1)
         qm(ng+1-k,2) = qmord(k,2)
         qm(ng+1-k,3) = qmord(k,3)
      enddo
      return
      end
c
      subroutine orderinverseqm(nph,qm,qmord)
      implicit double precision (a-h,o-z)
      parameter (nmax = 40)
      dimension qm(nmax,0:3),qmord(nmax,0:3),en(nph),sen(nph)
      dimension indexes(nph)
      do k = 1,nph
         en(k) = qm(k,0)
      enddo
      call sort_vector_inverse(en,sen,indexes,nph)
      do k = 1,nph
         qmord(k,0) = qm(indexes(k),0)
         qmord(k,1) = qm(indexes(k),1)
         qmord(k,2) = qm(indexes(k),2)
         qmord(k,3) = qm(indexes(k),3)
      enddo
      do k = nph+1,nmax
         qmord(k,0) = 0.d0
         qmord(k,1) = 0.d0
         qmord(k,2) = 0.d0
         qmord(k,3) = 0.d0
      enddo
      return
      end
******************************************
      subroutine sort_vector(v,sv,iold,n)
! pensare quick sort!!
      implicit double precision (a-h,o-z)
      dimension iold(n)
      dimension v(n),sv(n)
      do k = 1,n
         iold(k) = k
         sv(k) = v(k)
      enddo
      do i = 1,n
         vmax = sv(i)
         do j=i+1,n
            if (v(j).gt.vmax) then
               vmax = sv(j)
               vtmp = sv(i)
               sv(i) = sv(j)
               sv(j) = vtmp

               ktmp = iold(i)
               iold(i) = iold(j)
               iold(j) = ktmp
            endif
         enddo
      enddo
      return
      end

      subroutine sort_vector_inverse(v,sv,iold,n)
! pensare quick sort!!
      implicit double precision (a-h,o-z)
      dimension iold(n)
      dimension v(n),sv(n)
C this uses the quick sort algorithm reported at the end of the file
c      call SORTRX(n,v,iold)
c      do k = 1,n
c         sv(k) = v(iold(k))
c      enddo
c      return
      do k = 1,n
         iold(k) = k
         sv(k) = v(k)
      enddo
      do i = 1,n
         vmin = sv(i)
         do j=i+1,n
            if (v(j).lt.vmin) then
               vmin = sv(j)
               vtmp = sv(i)
               sv(i) = sv(j)
               sv(j) = vtmp
               ktmp = iold(i)
               iold(i) = iold(j)
               iold(j) = ktmp
            endif
         enddo
      enddo
      return
      end
