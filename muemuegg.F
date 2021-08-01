****************************************************
      double precision function muemuegg(p1,p2,p3,p4,k1,k2)
c      implicit none !double precision (a-h,o-z)
      implicit double precision (Z)
#include "invariants.h"
      double precision k1(0:3),k2(0:3),k1k2,p2mp4(0:3)
      double precision p1k1,p2k1,p3k1,p4k1
      double precision p1k2,p2k2,p3k2,p4k2
      double precision p1k1m1,p2k1m1,p3k1m1,p4k1m1
      double precision p1k2m1,p2k2m1,p3k2m1,p4k2m1
      integer iarun
      common/ialpharunning/iarun

      double precision tm100,tm110,tm102,tm112
      
      double precision el8
      integer ifirst
      common/ifirstmuemuegg/el8,ifirst
      data ifirst /0/

!      double precision vp13,vp24 already in invariants
      double precision vp13k1,vp24k1
      double precision vp13k2,vp24k2
      double precision vp13k1k2,vp24k1k2
      common/vps2gtree/!already in invariants and a common there vp13,vp24,
     .     vp13k1,vp24k1,
     .     vp13k2,vp24k2,
     .     vp13k1k2,vp24k1k2

      double precision vp13k1_1,vp24k1_1
      double precision vp13k2_1,vp24k2_1
      double precision vp13k1k2_1,vp24k1k2_1
      double precision vp13k1_2,vp24k1_2
      double precision vp13k2_2,vp24k2_2
      double precision vp13k1k2_2,vp24k1k2_2
      double precision vp13_1,vp24_1
      double precision vp13_2,vp24_2
      common/vpsmuemueggforreweight/vp13k1_1,vp24k1_1,vp13k2_1,vp24k2_1,
     .     vp13k1k2_1,vp24k1k2_1,vp13k1_2,vp24k1_2,
     .     vp13k2_2,vp24k2_2,vp13k1k2_2,vp24k1k2_2,
     .     vp13_1,vp24_1,vp13_2,vp24_2



      
      if (ifirst.eq.0) then
         vp13     = 1.d0
         vp24     = 1.d0
         vp13k1   = 1.d0
         vp24k1   = 1.d0
         vp13k2   = 1.d0
         vp24k2   = 1.d0
         vp13k1k2 = 1.d0
         vp24k1k2 = 1.d0


         vp13k1_1 = 1.d0
         vp24k1_1 = 1.d0
         vp13k2_1 = 1.d0
         vp24k2_1 = 1.d0
         vp13k1k2_1 = 1.d0
         vp24k1k2_1 = 1.d0
         vp13k1_2 = 1.d0
         vp24k1_2 = 1.d0
         vp13k2_2 = 1.d0
         vp24k2_2 = 1.d0
         vp13k1k2_2 = 1.d0
         vp24k1k2_2 = 1.d0
         vp13_1 = 1.d0
         vp24_1 = 1.d0
         vp13_2 = 1.d0
         vp24_2 = 1.d0
                  
         el8 = ec4*ec4
         ifirst = 1
      endif

      reweightLO  = 0.d0
      reweightNLO = 0.d0
      
      call loadinvariants(p1,p2,p3,p4)
      p1k1 = dot(p1,k1)
      p2k1 = dot(p2,k1)
      p3k1 = dot(p3,k1)
      p4k1 = dot(p4,k1)
      p1k2 = dot(p1,k2)
      p2k2 = dot(p2,k2)
      p3k2 = dot(p3,k2)
      p4k2 = dot(p4,k2)
      k1k2 = dot(k1,k2)
      
      p1k1m1 = 1.d0/p1k1
      p2k1m1 = 1.d0/p2k1
      p3k1m1 = 1.d0/p3k1
      p4k1m1 = 1.d0/p4k1
      p1k2m1 = 1.d0/p1k2
      p2k2m1 = 1.d0/p2k2
      p3k2m1 = 1.d0/p3k2
      p4k2m1 = 1.d0/p4k2

      p2mp4 = p2 - p4
      
      tm100 = 1.d0/dot(p2mp4,p2mp4)
      tm110 = 1.d0/dot(p2mp4-k1,p2mp4-k1)
      tm102 = 1.d0/dot(p2mp4-k2,p2mp4-k2)
      tm112 = 1.d0/dot(p2mp4-k1-k2,p2mp4-k1-k2)

***   VP STUFF
      if (ifixedorder.eq.0.and.iarun.eq.1) then
         call getVP1loop(1.d0/tm100,dal,dah)
         vp24   = 1.d0/(1.d0-dal-dah)
         vp24_1 = 1.d0
         vp24_2 = 1.d0/(1.d0-dal)
         call getVP1loop(1.d0/tm110,dal,dah)
         vp24k1   = 1.d0/(1.d0-dal-dah)
         vp24k1_1 = 1.d0
         vp24k1_2 = 1.d0/(1.d0-dal)
         call getVP1loop(1.d0/tm102,dal,dah)
         vp24k2   = 1.d0/(1.d0-dal-dah)
         vp24k2_1 = 1.d0
         vp24k2_2 = 1.d0/(1.d0-dal)
         call getVP1loop(1.d0/tm112,dal,dah)      
         vp24k1k2   = 1.d0/(1.d0-dal-dah)
         vp24k1k2_1 = 1.d0
         vp24k1k2_2 = 1.d0/(1.d0-dal)
      endif
********
c      print*,'Reenable includes in muemuegg.F!!! stopping'
c      stop
      muemuegg = 0.d0

cc always the matrix element for mu+ used, by crossing symmetry in matrix_model.F      
      include 'muemuegg-plus.f'

!! this part is commented in the file above      
      muemuegg=tm112**2*vp24k1k2**2*Z7_ + tm102*tm112*vp24k2*vp24k1k2*
     & Z12_ + tm102**2*vp24k2**2*Z5_ + tm110*tm112*vp24k1*vp24k1k2*Z25_
     &  + tm110*tm102*vp24k1*vp24k2*Z13_ + tm110**2*vp24k1**2*Z8_ + 
     & tm100*tm112*vp24*vp24k1k2*Z11_ + tm100*tm102*vp24*vp24k2*Z9_ + 
     & tm100*tm110*vp24*vp24k1*Z6_ + tm100**2*vp24**2*Z10_

      extraweights(0) = 1.d0

      extraweights(1) =tm112**2*vp24k1k2_1**2*Z7_
     &     + tm102*tm112*vp24k2_1*vp24k1k2_1*Z12_ 
     &     + tm102**2*vp24k2_1**2*Z5_
     &     + tm110*tm112*vp24k1_1*vp24k1k2_1*Z25_
     &     + tm110*tm102*vp24k1_1*vp24k2_1*Z13_
     &     + tm110**2*vp24k1_1**2*Z8_ + 
     &     tm100*tm112*vp24_1*vp24k1k2_1*Z11_
     &     + tm100*tm102*vp24_1*vp24k2_1*Z9_ + 
     &     tm100*tm110*vp24_1*vp24k1_1*Z6_ + tm100**2*vp24_1**2*Z10_

      extraweights(2) =tm112**2*vp24k1k2_2**2*Z7_
     &     + tm102*tm112*vp24k2_2*vp24k1k2_2*Z12_ 
     &     + tm102**2*vp24k2_2**2*Z5_
     &     + tm110*tm112*vp24k1_2*vp24k1k2_2*Z25_
     &     + tm110*tm102*vp24k1_2*vp24k2_2*Z13_
     &     + tm110**2*vp24k1_2**2*Z8_ + 
     &     tm100*tm112*vp24_2*vp24k1k2_2*Z11_
     &     + tm100*tm102*vp24_2*vp24k2_2*Z9_ + 
     &     tm100*tm110*vp24_2*vp24k1_2*Z6_ + tm100**2*vp24_2**2*Z10_
      

      extraweights(1:2) = extraweights(1:2)/muemuegg
      
ccccccccccccc
ccc not used to speedup compilation time, in matrix_model.F crossing symmetry is used
c      if (Qmu.lt.0.d0) then
c         include 'muemuegg-minus.f'
c      else
c         include 'muemuegg-plus.f'
c     endif
cccccccccccccccccccccccccccc
      
      muemuegg = el8 * muemuegg
      return
      end
