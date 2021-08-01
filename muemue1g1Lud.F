      double precision function muemue1g1Lupdown(p1,p2,p3,p4,k)
c      implicit none
      implicit double precision (Z) ! posso metterli reali perche' ho raccolto in form tutto cio' che e' imm...
#include "invariants.h"
#include "funsdeccmn1g1L.h"      
      double precision t13m1,t24m1,p1km1,p2km1,p3km1,p4km1
      double precision muemueg1L ! this is the variable in form
      double precision k(0:3)
      double precision cover,el8
*      double precision p1k,p2k,p3k,p4k
      common/ifirstmuemue1g1Lupdown/cover,el8,ifirst
      double complex vpolc
      external vpolc
      integer ifirst,iii
      data ifirst /0/
      integer iarun
      common/ialpharunning/iarun

      integer Amp(24)
      common/switchampls/Amp! for debugging IR...

      double precision p1_eps, p2_eps, p3_eps, p4_eps, k_eps
      double precision p1_epss,p2_epss,p3_epss,p4_epss,k_epss
      common/gaugeinv/p1_eps, p2_eps, p3_eps, p4_eps, k_eps,
     .     p1_epss,p2_epss,p3_epss,p4_epss,k_epss

      integer mulined,mulinec,elined,elinec
      common/subsetsBP/mulined,mulinec,elined,elinec

      double precision kf(0:3),p1_kf,p2_kf,p3_kf,p4_kf,k_kf

      double precision pik
      external pik

      integer iaxial,itnir1,itir1,itnir2,itir2,itnir3,itir3,itnir4,itir4
      integer ib1dir,ib1dnir,ib2dir,ib2dnir
      integer ib3dir,ib3dnir,ib4dir,ib4dnir
      integer ib1cir,ib1cnir,ib2cir,ib2cnir
      integer ib3cir,ib3cnir,ib4cir,ib4cnir
      integer ip1dir,ip1dnir,ip2dir,ip2dnir
      integer ip1cir,ip1cnir,ip2cir,ip2cnir            
      common/IRAXFACTORStest/iaxial,itnir1,itir1,itnir2,
     .     itir2,itnir3,itir3,itnir4,itir4,
     .     ib1dir,ib1dnir,ib2dir,ib2dnir,
     .     ib3dir,ib3dnir,ib4dir,ib4dnir,
     .     ib1cir,ib1cnir,ib2cir,ib2cnir,
     .     ib3cir,ib3cnir,ib4cir,ib4cnir,
     .     ip1dir,ip1dnir,ip2dir,ip2dnir,
     .     ip1cir,ip1cnir,ip2cir,ip2cnir            

      
      
      if (ifirst.eq.0) then

         mulined = 1
         mulinec = 1
         elined  = 1
         elinec  = 1

         vp13 = 1.d0
         vp24 = 1.d0
         
         el8 = ec4*ec4
         cover  = el8*pi216m1
         ifirst = 1
         Amp = 1
c         do iii = 1,100
c            print*,'Remember to set Amp(i)=1 here and in form later on!'
c         enddo
      endif

      muemue1g1Lupdown = 0.d0
      if (QRmu.eq.0.or.QRe.eq.0) return
      
      kf(0)   = k(0)
      kf(1:3) = -k(1:3)
      p1_kf   = dot(p1,kf)
      p2_kf   = dot(p2,kf)
      p3_kf   = dot(p3,kf)
      p4_kf   = dot(p4,kf)
      k_kf    = dot(k,kf)

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

c      print*,'Reenable includes in muemue1g1Lud.F!!! stopping'
      include 'muemue1g1Lupdown.f'
      
      muemue1g1Lupdown = muemueg1L*cover

      return
      end
****************************************************
      double precision function muemue1g1Lupdown_ten(p1,p2,p3,p4,k)
*      implicit none
      implicit double precision (Z) ! posso metterli reali perche' ho raccolto in form tutto cio' che e' imm...
#include "invariants.h"
#include "funsdeccmn1g1L.h"      
      double precision t13m1,t24m1,p1km1,p2km1,p3km1,p4km1
      double precision muemueg1L ! this is the variable in form
      double precision k(0:3)
      double precision cover,el8
*      double precision p1k,p2k,p3k,p4k
      common/ifirstmuemue1g1Lupdownten/cover,el8,zerop,ifirst
      double complex vpolc
      external vpolc
      integer ifirst,iii
      data ifirst /0/

      double precision zerop(0:3)
      double complex dotten
      external dotten
      double complex test1,test2,test3,test4,test5

      integer Amp(24)
      common/switchampls/Amp! for debugging IR...
      
      if (ifirst.eq.0) then
         Amp = 1
         el8 = ec4*ec4
         cover  = el8*pi216m1
         zerop = 0.d0
         ifirst = 1
      endif

      muemue1g1Lupdown_ten = 0.d0
      if (QRmu.eq.0.or.QRe.eq.0) return
      
      t13m1 = 1.d0/t13
      t24m1 = 1.d0/t24
      
      p1km1 = 1.d0/p1k
      p2km1 = 1.d0/p2k
      p3km1 = 1.d0/p3k
      p4km1 = 1.d0/p4k

      print*,'reenable print of muemueg1L in the form!!!'
c      include 'muemue1g1Lupdownwten.f'
      
      muemue1g1Lupdown_ten = muemueg1L*cover
      
      return
      end
************************************************************************
      double precision function pik(p,k)
      double precision p(0:3), k(0:3)
      double precision dot
      external dot
      pik = dot(p,k)
      return
      end
