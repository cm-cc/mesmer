      double precision function muemue1g1Lnoupdown(p1,p2,p3,p4,k)
c      implicit none
      implicit double precision (Z) ! posso metterli reali perche' ho raccolto in form tutto cio' che e' imm...
#include "invariants.h"
#include "funsdeccmn1g1L.h"      
      double precision t13m1,t24m1,p1km1,p2km1,p3km1,p4km1
      double precision muemueg1L ! this is the variable in form
      double precision k(0:3),kf(0:3),p1_kf,p2_kf,p3_kf,p4_kf,k_kf
      double precision cover,el8
*      double precision p1k,p2k,p3k,p4k
      common/ifirstmuemue1g1L/cover,el8,DELTA,ifirst
      double complex vpolc
      external vpolc
      integer ifirst,iii
      data ifirst /0/
      integer iarun
      common/ialpharunning/iarun
      
      double precision DELTA
#ifdef QUAD
      real*16 getdelta
#else
      double precision getdelta
#endif
      external getdelta
      
      if (ifirst.eq.0) then
         do iii = 1,2
c            print*,'Did you call getQEDcts,loadQEDfuns1g1L'//
c     .           ' and loadivariants1g first?'
         enddo
cc         DELTA  = getdelta()
         el8    = ec4*ec4
         cover  = el8*pi216m1
         ifirst = 1

         vp13 = 1.d0
         vp24 = 1.d0
      endif

      muemue1g1Lnoupdown = 0.d0

***   VP stuff
      if (ifixedorder.eq.0.and.iarun.eq.1) then
         call getVP1loop(t13,dal,dah)
         vp13 = 1.d0/(1.d0-dal-dah)         
         call getVP1loop(t24,dal,dah)
         vp24 = 1.d0/(1.d0-dal-dah)
      endif      
***   end VP stuff      

      t13m1 = 1.d0/t13
      t24m1 = 1.d0/t24
      
      p1km1 = 1.d0/p1k
      p2km1 = 1.d0/p2k
      p3km1 = 1.d0/p3k
      p4km1 = 1.d0/p4k

      kf(0)   = k(0)
      kf(1:3) = -k(1:3)

      p1_kf = dot(p1,kf)
      p2_kf = dot(p2,kf)
      p3_kf = dot(p3,kf)
      p4_kf = dot(p4,kf)
      k_kf  = dot(k,kf)
      
c      print*,'Reenable includes in muemue1g1Lnoud.F!!! stopping'
c      stop
! the default now includes the axial gauge!
      include 'muemue1g1Lnoupdown.f'
      
      muemue1g1Lnoupdown = muemueg1L*cover

      return
      end
****************************************************
