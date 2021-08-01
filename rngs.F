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
c      stop
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
      integer isync,seed,nsearch,iverbose,nw
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync
      parameter (nvec=1000)
      double precision r(nvec)
      common/getrndwvcmn/r,last,ifirst
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
      integer isync,seed,nsearch,iverbose,nw
      common/intinput/iwriteout,iseed,nsearch,iverbose,nw,isync
      dimension xi(n)
      parameter (nvec=1000)
      parameter (nvecsync=30) ! must be <= nvec
      double precision r(nvec)
      common/getrndwvcmn/r,last,ifirst
      data ifirst,last /0,nvec/
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

      
