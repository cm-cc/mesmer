      subroutine myelasticityopt(the0,thm0,nphpnc,ththr,dthr,
     .     theout,thmout,dist)
! all angles in mrad and distances!
      implicit double precision (a-h,m,o-z)
      logical lcond
      integer ifirst
      common/elasticityparam/emubeam0,ame,amu,ame2,amu2
      common/myelasticityoptcmn/thmumax,gr,grm1,pi,eps,ifirst
      data ifirst /0/
      
      if (ifirst.eq.0) then
         gr   = (sqrt(5.d0) + 1.d0)*0.5d0 ! golden ratio
         grm1 = 1.d0/gr
         pi   = 4.d0*atan(1.d0)
         ame2 = ame*ame
         amu2 = amu*amu
         thmumax = 1.d3*acos(sqrt(1.d0-ame2/amu2))
         eps     = 1.d-7 * 1.d3 ! because here is mrad
         ifirst  = 1
      endif      

      dist = 0.d0
      if (nphpnc.eq.0) return

      it1 = 0
      it2 = 0

      dthr2 = dthr*dthr

      Emu  = emubeam0
      beta = Emu + ame
** th_e where the correlation th_e/th_mu has the maximum, th_mu = thmumax
      themax = 1d3*acos(beta/sqrt(beta*beta + amu2-ame2))

cc https://en.wikipedia.org/wiki/Golden-section_search
***   first branch
      a = 0.d0   !themax
      b = themax !ththr ! pi*0.5d0

      c = b - (b-a)*grm1
      d = a + (b-a)*grm1
c      thmc = funthmmrad(c)
c      thmd = funthmmrad(d)
c      dc   = (thmc-thm0)**2 + (c-the0)**2
c      dd   = (thmd-thm0)**2 + (d-the0)**2
c      lcond = dc.gt.dthr2.or.dd.gt.dthr2

c      do while (abs(a-b).gt.eps.and.lcond)
      do while (abs(a-b).gt.eps)
         it1  = it1+1
         thmc = funthmmrad(c)
         thmd = funthmmrad(d)
         dc = (thmc-thm0)**2 + (c-the0)**2
         dd = (thmd-thm0)**2 + (d-the0)**2
c         lcond = dc.gt.dthr2.or.dd.gt.dthr2
         
         if (dc.lt.dd) then
            b = d
         else
            a = c
         endif
         c = b - (b-a)*grm1
         d = a + (b-a)*grm1
      enddo
      themindist1 = (a+b)*0.5d0
      thmmindist1 = funthmmrad(themindist1)
      dist1 = sqrt((thmmindist1-thm0)**2 + (themindist1-the0)**2)

      dist2 = 1d3 * dist1

      if (dist1.gt.dthr) then
***   second branch
         a = themax ! 0.d0
         b = ththr  ! pi*0.5d0

         c = b - (b-a)*grm1
         d = a + (b-a)*grm1
c         thmc = funthmmrad(c)
c         thmd = funthmmrad(d)
c         dc   = (thmc-thm0)**2 + (c-the0)**2
c         dd   = (thmd-thm0)**2 + (d-the0)**2
c         lcond = dc.gt.dthr2.or.dd.gt.dthr2
      
c         do while (abs(a-b).gt.eps.and.lcond)
         do while (abs(a-b).gt.eps)
            it2 = it2 + 1
            thmc = funthmmrad(c)
            thmd = funthmmrad(d)
            dc = (thmc-thm0)**2 + (c-the0)**2
            dd = (thmd-thm0)**2 + (d-the0)**2
c            lcond = dc.gt.dthr2.or.dd.gt.dthr2
            
            if (dc.lt.dd) then
               b = d
            else
               a = c
            endif
            c = b - (b-a)*grm1
            d = a + (b-a)*grm1
         enddo
         themindist2 = (a+b)*0.5d0
         thmmindist2 = funthmmrad(themindist2)
         dist2 = sqrt((thmmindist2-thm0)**2 + (themindist2-the0)**2)
      endif

      if (dist1.lt.dist2) then
         theout = themindist1
         thmout = thmmindist1
         dist   = dist1
      else
         theout = themindist2
         thmout = thmmindist2
         dist   = dist2
      endif      
      return
      end
***************************************************************
***************************************************************
      subroutine myelasticity2(the0,thm0,theout,thmout,dist1,dist2,dist)
      implicit double precision (a-h,m,o-z)
      double precision invphi,invphi2
      common/elasticityparam/emubeam0,ame,amu,ame2,amu2

      common/myelasticitycmnv2/thmumax,gr,grm1,pi,eps,rd,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         gr = (sqrt(5.d0) + 1.d0)*0.5d0 ! golden ratio
         grm1 = 1.d0/gr
         rd = 1.d0/3.d0
         
         invphi  = (sqrt(5.d0) - 1.d0)*0.5d0
         invphi2 = (3.d0 - sqrt(5.d0))*0.5d0
         
         pi = 4.d0*atan(1.d0)
         
         ame2 = ame*ame
         amu2 = amu*amu

         thmumax = acos(sqrt(1.d0-ame2/amu2))
         eps = 1d-7         
         ifirst = 1
      endif

      Emu  = emubeam0
      beta = Emu + ame
** th_e where the correlation th_e/th_mu has the maximum, th_mu = thmumax
      themax = acos(beta/sqrt(beta*beta + amu2-ame2))
      
cc https://en.wikipedia.org/wiki/Ternary_search
***   first branch
      a = 0.d0
      b = themax
      it1 = 0
      do while (abs(a-b).gt.eps)
         it1 = it1 + 1
         a3rd = a + (b-a)*rd
         b3rd = b - (b-a)*rd

         thma = funthm(a3rd)
         thmb = funthm(b3rd)
         da = (thma-thm0)**2 + (a3rd-the0)**2
         db = (thmb-thm0)**2 + (b3rd-the0)**2
         
         if (da.lt.db) then
            b = b3rd
         else
            a = a3rd
         endif
      enddo
         
      themindist1 = (a+b)*0.5d0
      thmmindist1 = funthm(themindist1)
      dist1 = sqrt((thmmindist1-thm0)**2 + (themindist1-the0)**2)
      
***   second branch
      a = themax
      b = pi*0.5d0
      it2 = 0
      do while (abs(a-b).gt.eps)
         it2 = it2 + 1
         a3rd = a + (b-a)*rd
         b3rd = b - (b-a)*rd

         thma = funthm(a3rd)
         thmb = funthm(b3rd)
         da = (thma-thm0)**2 + (a3rd-the0)**2
         db = (thmb-thm0)**2 + (b3rd-the0)**2
         
         if (da.lt.db) then
            b = b3rd
         else
            a = a3rd
         endif
      enddo
      
      themindist2 = (a+b)*0.5d0
      thmmindist2 = funthm(themindist2)
      dist2 = sqrt((thmmindist2-thm0)**2 + (themindist2-the0)**2)

 333  if (dist1.lt.dist2) then
         theout = themindist1
         thmout = thmmindist1
         dist   = dist1
      else
         theout = themindist2
         thmout = thmmindist2
         dist   = dist2
      endif
      return
      end
******************
      subroutine myelasticity(the0,thm0,theout,thmout,dist1,dist2,dist)
      implicit double precision (a-h,m,o-z)
      common/elasticityparam/emubeam0,ame,amu,ame2,amu2

      common/myelasticitycmn/thmumax,gr,grm1,pi,eps,ifirst
      data ifirst /0/

      if (ifirst.eq.0) then
         gr = (sqrt(5.d0) + 1.d0)*0.5d0 ! golden ratio
         grm1 = 1.d0/gr
         pi = 4.d0*atan(1.d0)
         ame2 = ame*ame
         amu2 = amu*amu
         thmumax = acos(sqrt(1.d0-ame2/amu2))
         eps = 1d-7
         ifirst = 1
      endif

      Emu  = emubeam0
      beta = Emu + ame
** th_e where the correlation th_e/th_mu has the maximum, th_mu = thmumax
      themax = acos(beta/sqrt(beta*beta + amu2-ame2))

cc https://en.wikipedia.org/wiki/Golden-section_search
***   first branch
      a = 0.d0
      b = themax

      c = b - (b-a)*grm1
      d = a + (b-a)*grm1
      it1 = 0

      dc = 100.d0*eps
      dd = 0.d0      
      do while (abs(b-a).gt.eps)
c      do while (abs(dc-dd).gt.eps)
         it1  = it1+1
         thmc = funthm(c)
         thmd = funthm(d)
         dc = (thmc-thm0)**2 + (c-the0)**2
         dd = (thmd-thm0)**2 + (d-the0)**2
         if (dc.lt.dd) then
            b = d
         else
            a = c
         endif
         c = b - (b-a)*grm1
         d = a + (b-a)*grm1
      enddo
      themindist1 = (a+b)*0.5d0
      thmmindist1 = funthm(themindist1)
      dist1 = sqrt((thmmindist1-thm0)**2 + (themindist1-the0)**2)

***   second branch
      a = themax
      b = pi*0.5d0

      c = b - (b-a)*grm1
      d = a + (b-a)*grm1
      it2 = 0
      dc = 100.d0*eps
      dd = 0.d0
      
      do while (abs(b-a).gt.eps)
c      do while (abs(dc-dd).gt.eps)
         it2 = it2 + 1
         thmc = funthm(c)
         thmd = funthm(d)
         dc = (thmc-thm0)**2 + (c-the0)**2
         dd = (thmd-thm0)**2 + (d-the0)**2
         if (dc.lt.dd) then
            b = d
         else
            a = c
         endif
         c = b - (b-a)*grm1
         d = a + (b-a)*grm1
      enddo
      themindist2 = (a+b)*0.5d0
      thmmindist2 = funthm(themindist2)
      dist2 = sqrt((thmmindist2-thm0)**2 + (themindist2-the0)**2)

 333  if (dist1.lt.dist2) then
         theout = themindist1
         thmout = thmmindist1
         dist   = dist1
      else
         theout = themindist2
         thmout = thmmindist2
         dist   = dist2
      endif      
c      print*,it1,it2
      
      return
      end
******************
      function funthm(the)
      implicit double precision (a-h,m,o-z)
      common/elasticityparam/emubeam0,ame,amu,ame2,amu2
      Emu = emubeam0
      pmu = sqrt(Emu*Emu-amu2)
      r   = sqrt(Emu*Emu-amu2)/(Emu+ame)
      ce  = cos(the)
      Eef  = ame*(1. + r*r*ce*ce)/(1. - r*r*ce*ce)
      pzmu = pmu - sqrt(Eef*Eef-ame2) * ce
      Emuf = Emu + ame - Eef
      cmu  = pzmu/sqrt(Emuf*Emuf-amu2)
      funthm = acos(cmu)
      return
      end
******************
******************
      function funthmmrad(the)
      implicit double precision (a-h,m,o-z)
      common/elasticityparam/emubeam0,ame,amu,ame2,amu2
      Emu = emubeam0
      pmu = sqrt(Emu*Emu-amu2)
      r   = sqrt(Emu*Emu-amu2)/(Emu+ame)
      ce  = cos(the*1d-3)
      Eef  = ame*(1. + r*r*ce*ce)/(1. - r*r*ce*ce)
      pzmu = pmu - sqrt(Eef*Eef-ame2) * ce
      Emuf = Emu + ame - Eef
      cmu  = pzmu/sqrt(Emuf*Emuf-amu2)
      funthmmrad = 1d3*acos(cmu)
      return
      end
******************      
      subroutine elasticity_distance(thetaein,thetamuin,
     +                               thetaeout,thetamuout)
      implicit none
      real*8 thetaein,thetamuin
      real*8 thetaeout,thetamuout
      real*8 thetae,thetamu
      real*8 thetaebar,thetamubar
      common/currentpoint/thetaebar,thetamubar
      real*8 thetamup
      external thetamup
      real*8 Dp
      external Dp
      real*8 D
      external D
      real*8 rtbis
      real*8 rtsafe
      external rtbis,rtsafe
      integer itest
      real*8 xlow,xup,xup0,precision
      real*8 solution
      integer iterations
      integer iter,i
      common/calls/iter
      real*8 thetae_elasticity_maximum
      real*8 thetae_elasticity_maximum_1,thetae_elasticity_maximum_2
      real*8 thetae_1,thetae_2
      real*8 delta,h,dp2
      real*8 emubeam0,ame,amu,ame2,amu2
      common/elasticityparam/emubeam0,ame,amu,ame2,amu2
      integer initial
c      parameter (initial=0)
      save initial
      save thetae_elasticity_maximum
*
      thetaebar= thetaein
      thetamubar= thetamuin
      precision= 1.d-11
*
      thetaeout = 0.d0
      thetamuout = 0.d0      
*
c      emubeam0= 150.d0
c      ame= 5.d-4
c      amu= 0.1d0
      ame2= ame*ame
      amu2= amu*amu


      
c      emubeam0= 150.d0
c      ame= 5.d-4
c      amu= 0.1d0
c      ame2= ame*ame
c      amu2= amu*amu
*
      if(initial.ne.1) then
*     
*     First step: 
*     Looking for the theta_e corresponding to maximum theta_mu
*     in the elasticity curve
*
         xlow= 0.d0
         xup= 0.1d0
         delta= xup-xlow
*
c$$$         iter= 0
c$$$         solution= rtbis(thetamup,xlow,xup,precision)
c$$$         print*,'looking for the root with bisection' 
c$$$         print*,'x_0= ',solution,' with ',iter,' bisections'
c$$$         print*,'f(x_0)= ',thetamup(solution)
*     
 10      iter= 0
         itest= 0
         solution= rtsafe(thetamup,xlow,xup,precision,itest)
         if(itest.gt.0) then
            xup= xup-delta/10.d0
            goto 10
         endif
c$$$         print*,'looking for the root with Newton-Raphson'
c$$$         print*,'x_0= ',solution,'  with ',iter,' iterations'
c$$$         print*,'f(x_0)= ',thetamup(solution)
c$$$         print*,'thetamu(solution)= ',thetamu(solution)
*
         thetae_elasticity_maximum= solution
         initial= 1
      endif
*     
*     Second step:
*     Looking for the point theta_e on the curve which minimizes the distance from
*     the point at hand (tehtaebar,thetamubar)
*
*     particular case: (thetaebar,thetamubar) on the elasticity curve
*
      if(abs(thetamu(thetaebar)-thetamubar).lt.1.d-10) then
         thetaeout= thetaebar
         thetamuout= thetamu(thetaebar)
         return
      endif
*
*     particular case: thetamubar > thetamu(thetae_elasticity_maximum)
*
      if(thetamubar.gt.thetamu(thetae_elasticity_maximum)) then
         thetaeout= thetae_elasticity_maximum
         thetamuout= thetamu(thetaeout)
c$$$         print*,' '
c$$$         print*,'thetamubar > maximum thetamu'
c$$$         print*,'thetae = ',thetae,'  thetamu= ',thetamu(thetae)
         return
      endif
*
* general case: first look at one solution over the whole range [0,0.1] rad
*
      xlow= 0.d0
      xup= 0.1d0
*
      xup0= xup      
      delta= xup-xlow

c$$$      print*,' '
c$$$      print*,'solution 1 over the whole range '
c$$$      print*,' '
c$$$      print*,'xlow    = ',xlow
c$$$      print*,'xup     = ',xup
c$$$      print*,'Dp(xlow)= ',Dp(xlow)
c$$$      print*,'Dp(xup) = ',Dp(xup)
*      
c      iter= 0
c      solution= rtbis(Dp,xlow,xup,precision)
c      print*,'looking for the root with bisection' 
c      print*,'x_0= ',solution,' with ',iter,' bisections'
c      print*,'f(x_0)= ',Dp(solution)
*     
 20   iter= 0
      itest= 0
      solution= rtsafe(Dp,xlow,xup,precision,itest)
      if(itest.gt.0) then
         xup= xup-delta/10.d0
         if(xup.gt.xlow) then
            goto 20
         else
            xup= xup0
 25         iter= 0
            itest= 0
            solution= rtsafe(Dp,xlow,xup,precision,itest)
            if(itest.gt.0) then
               xup= xup - delta/100.d0
               if(xup.gt.xlow) then
                  goto 25
               else
                  if(D(xlow).lt.D(xup0)) solution= xlow
                  if(D(xlow).ge.D(xup0)) solution= xup0
                  goto 28
               endif
            endif
         endif
      endif
      
c      print*,'looking for the root with Newton-Raphson'
c      print*,'x_0= ',solution,'  with ',iter,' iterations'
c      print*,'f(x_0)= ',Dp(solution)
*
 28   continue
      thetae_1= solution
*     
c      print*,' '
c      print*,'point     = ',thetaein,thetamuin
c      print*,'solution 1= ',thetae_1,thetamu(thetae_1),D(thetae_1)
c      print*,' '
c      print*,'thetae_1= ',thetae_1,'  D(thetae_1)= ',D(thetae_1)
c      print*,'thetae_1= ',thetae_1,'  thetamu= ',thetamu(thetae_1)
*
* then look at the second solution in one of the two sides 
* around thetae_elasticity_maximum
*
      if(thetaebar.gt.thetae_elasticity_maximum.and.
     +     thetamubar.gt.thetamu(thetaebar).and.
     +   thetae_1.gt.thetae_elasticity_maximum) then
         thetaeout= thetae_1
         thetamuout= thetamu(thetaeout)
*
* check the second derivative at the solution
*
         h= 1.d-6
         dp2= (Dp(thetaeout+h)-Dp(thetaeout-h))/2.d0/h
         if(dp2.lt.0) then
            print*,'second derivative at the solution < 0!!! ',dp2
            stop
         endif
         return
      endif
      if(thetaebar.lt.thetae_elasticity_maximum.and.
     +   thetae_1.lt.thetae_elasticity_maximum) then
         thetaeout= thetae_1
         thetamuout= thetamu(thetaeout)
*
* check the second derivative at the solution
*
         h= 1.d-6
         dp2= (Dp(thetaeout+h)-Dp(thetaeout-h))/2.d0/h
         if(dp2.lt.0) then
            print*,'second derivative at the solution < 0!!! ',dp2
            stop
         endif
         return
      endif
*
      if(thetaebar.gt.0.004d0.and.
     +   thetae_1.gt.thetae_elasticity_maximum) then
         thetaeout= thetae_1
         thetamuout= thetamu(thetaeout)
*
* check the second derivative at the solution
*
         h= 1.d-6
         dp2= (Dp(thetaeout+h)-Dp(thetaeout-h))/2.d0/h
         if(dp2.lt.0) then
            print*,'second derivative at the solution < 0!!! ',dp2
            stop
         endif

         return
      endif
      if(thetae_1.gt.thetae_elasticity_maximum) then
c     
         xlow= 0.d0
         xup= thetae_elasticity_maximum !- 1.d-5 !this 1.d-5 is to avoid
*                                        D(xlow) * D(xup) > 0
         xup0= xup
         delta= xup-xlow
c$$$         print*,' '
c$$$         print*,'second solution '
c$$$         print*,' '
c$$$         print*,'xlow    = ',xlow
c$$$         print*,'xup     = ',xup
c$$$         print*,'Dp(xlow)= ',Dp(xlow)
c$$$         print*,'Dp(xup) = ',Dp(xup)
*      
c         iter= 0
c         solution= rtbis(Dp,xlow,xup,precision)
c         print*,'looking for the root with bisection' 
c         print*,'x_0= ',solution,' with ',iter,' bisections'
c         print*,'f(x_0)= ',Dp(solution)
*
 30      iter= 0
         itest= 0
         solution= rtsafe(Dp,xlow,xup,precision,itest)
         if(itest.gt.0) then
            xup= xup-delta/10.d0
            if(xup.gt.xlow) then
               goto 30
            else
               xup= xup0
 35            iter= 0
               itest= 0
               solution= rtsafe(Dp,xlow,xup,precision,itest)
               if(itest.gt.0) then
                  xup= xup - delta/100.d0
                  if(xup.gt.xlow) then
                     goto 35
                  else
                     if(D(xlow).lt.D(xup0)) solution= xlow
                     if(D(xlow).ge.D(xup0)) solution= xup0
                     goto 50
                  endif
               endif
            endif
         endif
c         print*,'looking for the root with Newton-Raphson'
c         print*,'x_0= ',solution,'  with ',iter,' iterations'
c         print*,'f(x_0)= ',Dp(solution)
*
 50      continue
         thetae_2= solution
*     
c         print*,' '
c         print*,'thetae_2= ',thetae_2,'  D(thetae_2)= ',D(thetae_2)
c         print*,'thetae_2= ',thetae_2,'  thetamu= ',thetamu(thetae_2)
      else
         xlow= thetae_elasticity_maximum
         xup= 0.1d0
*
         xup0= xup
         delta= xup-xlow

c$$$         print*,' '
c$$$         print*,'second solution '
c$$$         print*,' '
c$$$         print*,'xlow    = ',xlow
c$$$         print*,'xup     = ',xup
c$$$         print*,'Dp(xlow)= ',Dp(xlow)
c$$$         print*,'Dp(xup) = ',Dp(xup)
*      
c         iter= 0
c         solution= rtbis(Dp,xlow,xup,precision)
c         print*,'looking for the root with bisection' 
c         print*,'x_0= ',solution,' with ',iter,' bisections'
c         print*,'f(x_0)= ',Dp(solution)
*
 40      iter= 0
         itest= 0
         solution= rtsafe(Dp,xlow,xup,precision,itest)
         if(itest.gt.0) then
            xup= xup-delta/10.d0
            if(xup.gt.xlow) then
               goto 40
            else
               xup= xup0
 45            iter= 0
               itest= 0
               solution= rtsafe(Dp,xlow,xup,precision,itest)
               if(itest.gt.0) then
                  xup= xup - delta/100.d0
                  if(xup.gt.xlow) then
                     goto 45
                  else
                     if(D(xlow).lt.D(xup0)) solution= xlow
                     if(D(xlow).ge.D(xup0)) solution= xup0
                     goto 60
                  endif
               endif
            endif
         endif
c         print*,'looking for the root with Newton-Raphson'
c         print*,'x_0= ',solution,'  with ',iter,' iterations'
c         print*,'f(x_0)= ',Dp(solution)
*
 60      continue
         thetae_2= solution
*     
c         print*,' '
c         print*,'thetae_2= ',thetae_2,'  D(thetae_2)= ',D(thetae_2)
c         print*,'thetae_2= ',thetae_2,'  thetamu= ',thetamu(thetae_2)
      endif
*
      if(thetae_2.ge.0.d0) then
         if(D(thetae_1).lt.D(thetae_2)) then
            thetaeout= thetae_1
            thetamuout= thetamu(thetaeout)
*
* check the second derivative at the solution
*
            h= 1.d-6
            dp2= (Dp(thetaeout+h)-Dp(thetaeout-h))/2.d0/h
            if(dp2.lt.0) then
               print*,'second derivative at the solution < 0!!! ',dp2
               stop
            endif
         else
            thetaeout= thetae_2
            thetamuout= thetamu(thetaeout)
*
* check the second derivative at the solution
*
            h= 1.d-6
            dp2= (Dp(thetaeout+h)-Dp(thetaeout-h))/2.d0/h
            if(dp2.lt.0) then
               print*,'second derivative at the solution < 0!!! ',dp2
               stop
            endif
         endif
c$$$         print*,' '
c$$$         print*,'point     = ',thetaein,thetamuin
c$$$         print*,'solution 1= ',thetae_1,thetamu(thetae_1),D(thetae_1)
c$$$         print*,'solution 2= ',thetae_2,thetamu(thetae_2),D(thetae_2)
      else
         thetaeout= thetae_1
         thetamuout= thetamu(thetaeout)
*
* check the second derivative at the solution
*
         h= 1.d-6
         dp2= (Dp(thetaeout+h)-Dp(thetaeout-h))/2.d0/h
         if(dp2.lt.0) then
            print*,'second derivative at the solution < 0!!! ',dp2
            stop
         endif
c$$$         print*,' '
c$$$         print*,'point     = ',thetaein,thetamuin
c$$$         print*,'solution 1= ',thetae_1,thetamu(thetae_1),D(thetae_1)
c$$$         print*,'solution 2 not found'
      endif
*
      return
*
      end
*
*     distance of the point (thetaebar,thetamubar)
*     with respect to the curve thetamu(thetae)
*     D= sqrt((thetaebar-thetae)^2+(thetamubar-thetamu(thetae))^2)
*
      real*8 function D(thetae)
      implicit none
      real*8 thetae
      real*8 thetamu
      external thetamu
      real*8 thetaebar,thetamubar
      common/currentpoint/thetaebar,thetamubar
      real*8 num,den
*
      D= sqrt((thetaebar-thetae)**2 + (thetamubar-thetamu(thetae))**2)
*
      return
      end
*
*
*     derivative of the distance of the point (thetaebar,thetamubar)
*     with respect to the curve thetamu(thetae)
*     D= sqrt((thetaebar-thetae)^2+(thetamubar-thetamu(thetae))^2)
*     We look for the zero of D(thetae)/dthetae, neglecting
*     the factor 1/D which is always > 0 by construction
*
      real*8 function Dp(thetae)
      implicit none
      real*8 thetae
      real*8 thetamup,thetamu
      external thetamu,thetamup
      real*8 thetaebar,thetamubar
      common/currentpoint/thetaebar,thetamubar
*
      Dp= thetae - thetaebar + thetamu(thetae)*thetamup(thetae)
     +   -thetamubar*thetamup(thetae)
      Dp= Dp 
     +  / sqrt((thetaebar-thetae)**2 + (thetamubar-thetamu(thetae))**2)
*
      return
      end
*
      real*8 function thetamup(thetae)
      implicit none
      real*8 thetae,thetamu
      real*8 r,pmu0
      real*8 num,den
      real*8 emubeam0,ame,amu,ame2,amu2
      common/elasticityparam/emubeam0,ame,amu,ame2,amu2
*
      pmu0= sqrt(emubeam0*emubeam0-amu2)
      r= pmu0/(emubeam0+ame)
      num= 2.d0*ame*r*(pmu0-r*(pmu0*r+2.d0*ame)*(cos(thetae))**2)
     +     *(cos(2.d0*thetae))
     +     -ame*r**2*(pmu0*r+2.d0*ame)*(sin(2.d0*thetae))**2
      den= (pmu0-r*(pmu0*r+2.d0*ame)*(cos(thetae))**2)**2
     +    +r*r*ame2*(sin(2.d0*thetae))**2
      thetamup= num/den
*
      return
      end
*
      real*8 function thetamu(thetae)
      implicit none
      real*8 thetae
      real*8 r,pmu0
      real*8 num,den,arg
      real*8 emubeam0,ame,amu,ame2,amu2
      common/elasticityparam/emubeam0,ame,amu,ame2,amu2
*
      pmu0= sqrt(emubeam0*emubeam0-amu2)
      r= pmu0/(emubeam0+ame)
      num= ame*r*(sin(2.d0*thetae))
      den= pmu0-r*(pmu0*r+2.d0*ame)*(cos(thetae))**2
      arg= num/den
      thetamu= atan(num/den)
*
      return
      end
*
c
c Using a combination of Newton-Raphson and bisection, 
c find the root of a function bracketed between x1 and x2. 
c The root, returned as the function value RTSAFE, will be refined 
c until its accuracy is known within +-xacc. funcd is a user-supplied 
c subroutine which returns both the function value and the first derivative 
c of the function.
c
      REAL*8 function rtsafe(funcd,x1,x2,xacc,icheck)
      implicit none
      integer maxit
      integer icheck
      real*8 x1,x2,xacc
      real*8 funcd
      external funcd
      parameter (maxit=100)      !maximum allowed number of iterations
      integer j
      real*8 df,dx,dxold,f,fh,fl,temp,xh,xl
      real*8 h
      integer iter
      common/calls/iter
*
      h= (x2-x1)/1.d6          !step for the derivative calculation
      call funczf(funcd,x1,h,fl,df)
      call funczb(funcd,x2,h,fh,df)
      if((fl.gt.0.d0.and.fh.gt.0.d0).or.(fl.lt.0.d0.and.fh.lt.0.d0)) 
     + then
c$$$         print*,'root must be bracketed in rtsafe'
         rtsafe= -x1  !just to remind that it has no meaning
         icheck= 1
         return
      endif
*
      if(abs(fl).eq.0.d0) then
         rtsafe=x1
         return
      elseif(abs(fh).eq.0.d0) then
         rtsafe=x2
         return
      elseif(fl.lt.0.d0) then     !orient the search so that f(x1)<0
         xl= x1
         xh= x2
      else
         xh= x1
         xl= x2
      endif
      rtsafe= 0.5d0*(x1+x2)       !initialize the guess for root
      dxold= abs(x2-x1)           !the "stepsize before last"
      dx= dxold                   !and the last step
      h= dxold/1.d6
      call funcz(funcd,rtsafe,h,f,df)
      iter= 0
      do j=1,maxit                !loop over allowed iterations
         iter= iter+1
         if(((rtsafe-xh)*df-f)
     +     *((rtsafe-xl)*df-f).gt.0.d0 .or. !bisect if Newton out of range
     +     abs(2.d0*f).gt.abs(dxold*df) ) then !or not decreasing fast enough  
           dxold= dx
           dx= 0.5d0*(xh-xl)
           rtsafe= xl+dx
           if(abs(xl-rtsafe).lt.1.d-14) return  !change in root is negligible
        else
           dxold= dx
           dx= f/df
           temp= rtsafe
           rtsafe= rtsafe - dx
           if(abs(temp-rtsafe).lt.1.d-14) return !Newto step acceptable: take it
        endif
        if(abs(dx).lt.xacc) return        !Convergence criterion
        h= dx/100.d0
        call funcz(funcd,rtsafe,h,f,df)   !new function evaluation per iteration
        if(f.lt.0.d0) then
           xl= rtsafe
        else
           xh= rtsafe
        endif
      enddo
      print*, 'rtsafe exceeding maximum iterations'
      return
      end
*
      subroutine funcz(func,x,h,f,df)
      implicit none
      real*8 func
      external func
      real*8 x,f,df
      real*8 h,xup,xlow
      f= func(x)
      xup= x+h
      xlow= x-h
      df= (func(xup)-func(xlow))/2.d0/h      !derivative with o(h^2)
      return
      end
*
      subroutine funczf(func,x,h,f,df)
      implicit none
      real*8 func
      external func
      real*8 x,f,df
      real*8 h,xup,xlow
      f= func(x)
      xup= x+h
      xlow= x
      df= (func(xup)-func(xlow))/h      !derivative with o(h^2)
      return
      end
*
      subroutine funczb(func,x,h,f,df)
      implicit none
      real*8 func
      external func
      real*8 x,f,df
      real*8 h,xup,xlow
      f= func(x)
      xup= x-h
      xlow= x
      df= (func(xup)-func(xlow))/(-h)      !derivative with o(h^2)
      return
      end
*
c 
c Using bisection, finds the root of a function func known 
c to lie between x1 and x2. The root, returned as rtbis, 
c will be refined until its accuracy is +- xacc
c 
      REAL*8 function rtbis(func,x1,x2,xacc)
      implicit none
      integer jmax   !maximum number of iterations for bisection
      real*8 x1,x2,xacc
      real*8 func
      external func
      parameter (jmax=40)      !2^-40 ~ 10^-12
      integer j
      real*8 dx,f,fmid,xmid
      integer iter
      common/calls/iter
*
      fmid= func(x2)
      f= func(x1)

      if(f*fmid.ge.0.d0) then
         print*,'root must be bracketed in rtbis'
         stop
      endif
      if(f.lt.0.d0) then !orient the search so that f>0 is at x+dx
         rtbis= x1
         dx= x2-x1
      else
         rtbis= x2
         dx= x1-x2
      endif
      do j=1,jmax
         iter= j
         dx= dx/2.d0
         xmid= rtbis+dx
         fmid= func(xmid)
         if(fmid.le.0.d0) rtbis=xmid
         if(abs(dx).lt.xacc.or.fmid.eq.0.d0) return  !Convergence
      enddo
      print*,'too many bisections in rtbis'

      return
      end
