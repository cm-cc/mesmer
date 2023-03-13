*******************************************************************************
      subroutine initrng(iseed)
      integer iseed,lux,k1,k2,ilux,n
      parameter (lux = 4,k1 = 0,k2 = 0) ! for RANLUX
      integer rnlxstatesize_c
      external rnlxstatesize_c
c      do n = 1,100
c         print*,'using old ranlux!!'
c      enddo
c      call rluxgo(lux,iseed,k1,k2) ! commenta
c      return  ! commenta
      ilux = 2
      call drnlx_c_init(ilux,iseed) ! for C version of ranlux, double precision
c      call srnlx_c_init(ilux,iseed) ! for C version of ranlux, single precision
c      print*,rnlxstatesize_c()
c     stop
      return
      end
************************************************
      subroutine getrndroutine(csi8,n)
! it can be modified by the user
      integer n
      double precision csi8(n)
      real*8 csitruedouble(n)
      real*4 csi4(n)
c      call ranlux(csi4,n)  ! commenta
c      csi8 = 1.d0 *csi4  ! commenta
c      return  ! commenta
c      call srnlx_c(csi4,n) ! for C version of ranlux, single precision
c      csi8 = 1.d0 *csi4
c      return
#ifndef QUAD
      call drnlx_c(csi8,n)      ! for C version of ranlux, double precision
#else
      call drnlx_c(csitruedouble,n) ! for C version of ranlux, double precision
      csi8 = csitruedouble
#endif
      return
      end
*******************************************************************************
      subroutine syncrndseq
      implicit double precision (a-h,o-z)
      integer*8 iwriteout
      integer isync,iseed,nsearch,iverbose,nw,istorver,icwvp
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync,
     .     istorver,icwvp,ivpfl
      parameter (nvec=1000)
      double precision r(nvec)
      common/getrndwvcmn/r
      common/lastusedrnd/last,ifirst
*****
      if (isync.eq.0) return ! does nothing in this case
*****
      last = nvec      
      return
      end
*******************************************************************************      
      subroutine getrnd(xi,n)
      implicit double precision (a-h,o-z)
      integer*8 iwriteout
      integer isync,iseed,nsearch,iverbose,nw,istorver,icwvp
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync,
     .     istorver,icwvp,ivpfl
      dimension xi(n)
      parameter (nvec=1000)
      parameter (nvecsync=30) ! must be <= nvec
      double precision r(nvec)
      common/getrndwvcmn/r
      common/lastusedrnd/last,ifirst
      data ifirst,last /0,nvec/
c      call getrndroutine(xi,n)
c      return      
      if (last.eq.nvec) then
         if (isync.eq.0) then
            call getrndroutine(r,nvec)
         else ! I do not need nvec random per event...
            call getrndroutine(r(1:nvecsync),nvecsync)
         endif
         last  = 0
      endif
      if ((last+n).le.nvec) then
         xi(1:n)=r((last+1):(last+n))
         last = last + n
      else         
         xi(1:(nvec-last))=r(last+1:nvec)
         call getrndroutine(r,nvec)
         xi((nvec-last+1):n) = r(1:(n-(nvec-last)))
         last = n - (nvec - last)
      endif
      return
      end
************************************************      
      subroutine getrndstatus(n)
      integer n(104)
      call rnlxdstate_c(n)
      return
      end
*******************************************************************************      
      subroutine getrnd_firstone(csi8,n)
! it can be modified by the user
      integer n
      double precision csi8(n)
      real*8 csitruedouble(n)
      real*4 csi4(n)
c      call ranlux(csi4,n)  ! commenta
c      csi8 = 1.d0 *csi4  ! commenta
c      return  ! commenta
c      call srnlx_c(csi4,n) ! for C version of ranlux, single precision
c      csi8 = 1.d0 *csi4
c      return
#ifndef QUAD
      call drnlx_c(csi8,n)      ! for C version of ranlux, double precision
#else
      call drnlx_c(csitruedouble,n) ! for C version of ranlux, double precision
      csi8 = csitruedouble
#endif
      return
      end
***************************************************************************
      subroutine anotherrng(csi,n)
      integer n,ifirst,k
      double precision csi(n)
      integer seed(33)      
      real*4 csi4(n)
      common/yetanotherrng/ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         seed = 0
         seed(1) =  1956357262
         seed(2) =  1751479972
         seed(3) = -485056068
         seed(4) = -1880652247
         seed(5) =  147017111
         seed(6) =  1123015416
         seed(7) = -506977708
         seed(8) = -1870866328
c         call random_seed(put = seed)
         ifirst = 1
         call srnlx_c_init(2,424242) ! for C version of ranlux, double precision

c         call gslrndinit(seed(1))
         
      endif
c      call gslrnd(csi,n)
c      return
      igo = 0
      do while (igo.eq.0)
c         call random_number(csi4)
         call srnlx_c(csi4,n)
         if (csi4(1).ne.0..and.csi4(1).ne.1.) igo = 1
      enddo
      csi = csi4
      do k = 1,n
         csi(k) = min(1.d0,csi(k))
         csi(k) = max(0.d0,csi(k))
      enddo
      return
      end
**************************************************************************      
      subroutine init_rnd_sequence
! I enclose it here so that can be used after warmup      
      implicit none
      character*20 a
      integer iii,k,icantor,djb_hash,n,idim
      parameter(idim = 1000)
      integer vec(idim),fl(1)
      double precision xis(idim)
      integer*8 iwriteout,ii,i8
      integer isync,iseed,nsearch,iverbose,nw,istorver,icwvp
      integer iseedext,ivpfl
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync,
     .     istorver,icwvp,ivpfl
      common/iseedoninput/iseedext
      double precision csi(2),ri
      real*4 csis(2)
      integer last,nvec, ifirst
      parameter (nvec=1000)
      common/lastusedrnd/last,ifirst
      integer ifirstinit
      common/initrndcmn/ifirstinit
      data ifirstinit /0/
! I call this routine in userinterface and at least once more (twice if warming up) in generete_event or main.
      if (ifirstinit.eq.0) then
c this was the first
         iseedext = iseed
         call getanotherseed(iseedext,iseed)
c         print*,iseed
c         stop
c         call initrng(4323432)
c         do k = 1,iseedext
c            call getrndroutine(csi,2)
c         enddo
c         iii   = int(csi(1)*1d7)
c         iii   = iii + int(csi(2)*1d7)
c         iseed = iii
         ifirstinit = 1
      endif
****  RNG initialization      
      call initrng(iseed)
****  end RNG initialization
      last = nvec ! to be sure that the random array is reloaded when init_rnd_sequence is called
      return
      end
****************************************************
      function djb_hash(str) result(hash)
      implicit none
      character(len=*),intent(in) :: str
      integer :: hash
      integer :: i
      hash = 5381
      do i=1,len(str)
         hash = (ishft(hash,5) + hash) + ichar(str(i:i))
      end do
      end function djb_hash
****************************************************
      subroutine converttobinary(i)
      implicit double precision (a-h,o-z)
      integer bin(31),shuffled(31)
      double precision xi(1)
      
      ic = 0
      j  = i
      bin = 0
      shuffled = 0
      if (j.eq.0) then
         ic = 1
         bin(1) = 0
      else
         do while(j.ne.0)
            ic = ic + 1
            if (mod(j,2).eq.0) then
               bin(ic) = 0
            else
               bin(ic) = 1
            endif
            j = j/2
         enddo
      endif

      print*,i,ic
      print*,'direct  ',bin(1:ic)

      j = 0
      do k = 1,ic
         j = j + bin(k)*2**(k-1)
      enddo
      print*,j
* shuffling bits
*      https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
      shuffled = bin
      do k = 1,ic-1
         call getrndroutine(xi,1)
         j = k + (ic-k)*xi(1) + 1
         itmp = shuffled(k)
         shuffled(k) = shuffled(j)
         shuffled(j) = itmp         
      enddo

      print*,'shuffled',shuffled(1:ic)
      j = 0
      do k = 1,ic
         j = j + shuffled(k)*2**(k-1)
      enddo
      print*,j

      call initrng(i)

      do k = 1,31
         call getrndroutine(xi,1)
         if (xi(1).lt.0.5d0) then
            shuffled(k) = 1
         else
            shuffled(k) = 0
         endif         
      enddo
      
      j = 0
      do k = 1,31
         j = j + shuffled(k)*2**(k-1)
      enddo
      print*,j      
      return
      end
****************************************************
      subroutine getanotherseed(i,j)
      implicit double precision (a-h,o-z)
      integer rndbit(31)
      integer bin(31)
      double precision xi(1)
c      goto 123
****
c      call initrng(i)
****
c      do k = 1,31
c         call getrndroutine(xi,1)
c         if (xi(1).ge.0.5d0) then
c            rndbit(k) = 1
c         else
c            rndbit(k) = 0
c         endif         
c      enddo
c      j = 0
c      do k = 1,31
c         j = j + rndbit(k)*2**(k-1)
c      enddo
c 123  continue ! just flipping the bits
      ic = 0
      j  = i
      bin = 0
      if (j.eq.0) then
         ic = 1
         bin(1) = 0
      else
         do while(j.ne.0)
            ic = ic + 1
            if (mod(j,2).eq.0) then
               bin(ic) = 0
            else
               bin(ic) = 1
            endif
            j = j/2
         enddo
      endif
      j = 0
      do k = 1,31
         bin(k) = 1-bin(k)
         j = j + bin(k)*2**(k-1)
      enddo
      return
      end
