      subroutine realpairs(p1,p2,p3,p4,p5,p6,m1,m2,ampl2)
c     implicit none
      implicit double precision (Z)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      double precision m1,m2,m12,m22,ampl2,el8
      double precision ame,ammu,convfac,alpha,pi
      common/parameters/ame,ammu,convfac,alpha,pi
      integer ifirst,iper,noper,iid

      common/identicalandperipheral/iper,noper,iid
**** charges to switch on radiation on legs (0 or 1), from invariants.h
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
**********************************************************************      
      common/realpairscmn/el8,m12,m22,ifirst
      data ifirst/0/
      double precision p1p2,p1p3,p1p4,p1p5,p1p6
      double precision      p2p3,p2p4,p2p5,p2p6
      double precision           p3p4,p3p5,p3p6
      double precision                p4p5,p4p6
      double precision                     p5p6,dot
      external dot

      if (ifirst.eq.0) then

         el8 = (4.d0*pi*alpha)**4


** SET IN phasespacepairs
c     iper = 1
c     iid  = 1
         
         ifirst = 1
      endif
      
c         m1 = ame
c         m2 = ammu
      m12 = m1*m1
      m22 = m2*m2

      p1p2 = dot(p1,p2)
      p1p3 = dot(p1,p3)
      p1p4 = dot(p1,p4)
      p1p5 = dot(p1,p5)
      p1p6 = dot(p1,p6)
      p2p3 = dot(p2,p3)
      p2p4 = dot(p2,p4)
      p2p5 = dot(p2,p5)
      p2p6 = dot(p2,p6)
      p3p4 = dot(p3,p4)
      p3p5 = dot(p3,p5)
      p3p6 = dot(p3,p6)
      p4p5 = dot(p4,p5)
      p4p6 = dot(p4,p6)
      p5p6 = dot(p5,p6)
      
      include 'realpairs_ampl2.f'
cccccccc      include 'realpairs_ampl2iid1iper1.f'
      
      ampl2 = el8 * ampl2

      return
      end
