      subroutine svfactor(model,ng,ecms,p1,p2,eps,sv,dsv)
      implicit double precision (a-h,l,o-z)
      dimension p1(0:3),p2(0:3),pin1(0:3),pin2(0:3),p1b(0:3),p2b(0:3)
      character*10 model
      character*6 ord
      double precision masses(4),charges(4)
      character*2 fs
      common/finalstate/fs
      common/qedORDER/ord
      common/parameters/ame,ammu,convfac,alpha,pi
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3)
     >     ,pin2r(0:3)
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeq2,bornmez,
     >     bornmenovp
      common/ialpharunning/iarun
      common/icountsv/icsvf,ifirst
      data ifirst,icsvf /0,0/
! not needed (yet) for mue scattering!      
      sv = 1.d0
      return
      
      icsvf = icsvf + 1

      dsv = 0.d0

      if (ifirst.eq.0) then
         do k = 1,10
         enddo
         ifirst = 1
      endif

      if (ord.eq.'born') then
         sv = 1.d0
         return
      endif

      if (ord.eq.'alpha') then
         if (ng.ge.1) then
            sv = 1.d0
            return
         endif         
         if (ng.eq.0) then
            s = ecms**2
            c = p1(3)/sqrt(tridot(p1,p1))
            t = -s/2.d0*(1.d0 - c)            

!            sigma0 = sdif_BCK(s,c)      
            sigma0  = bornme*phsp2b/flux2b
            sigmaq2 = bornmeq2*phsp2b/flux2b

            if (fs.eq.'gg') then
               sv = aintsvgg(s,c,eps)/sigma0
            else
c               sv = aintsv(s,t,eps)/sigma0
               sv = aintsv_new(pin1,pin2,p1,p2,eps)/sigma0
            endif

            if (model.eq.'ps') then
               aieps = -2.d0*log(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
! see below
               dp1p2 = log(2.d0*dot(pin1,pin2)/ame**2)
               dp1p3 = log(2.d0*dot(pin1,p1)/ame**2)
               dp1p4 = log(2.d0*dot(pin1,p2)/ame**2)
               dp2p3 = log(2.d0*dot(pin2,p1)/ame**2)
               dp2p4 = log(2.d0*dot(pin2,p2)/ame**2)
               dp3p4 = log(2.d0*dot(p1,p2)/ame**2)

               if (fs.eq.'ee'.or.fs.eq.'mm') then
                  lcollnew = dp1p2+dp1p3-dp1p4-dp2p3+dp2p4+dp3p4
                  lcollnew = 2.d0 * lcollnew - 4.d0
                  lcollnew = lcollnew/4.d0
                  arg = 2.d0*alpha/pi * lcollnew * aieps
**** NEW ***************
                  charges(1) = -1.d0
                  charges(2) = -1.d0
                  charges(3) = -1.d0
                  charges(4) = -1.d0
                  masses(1)  = ame
                  masses(2)  = ame
                  masses(3)  = ame
                  masses(4)  = ame
                  if (fs.eq.'mm') then
                     masses(3)  = ammu
                     masses(4)  = ammu
                  endif
                  lcoll1=eikonalintegral(pin1,pin2,p1,p2,charges,masses)
                  arg2  = 2.d0*alpha/pi * lcoll1 * aieps 
                  arg = arg2
**************************

               elseif (fs.eq.'gg') then
ccc                  return
                  lcollnew = dp1p2
                  lcollnew = lcollnew - 1.d0
                  arg = alpha/pi * lcollnew * aieps
               endif
               sv = 1.d0 - arg
            endif
            return
         endif
      endif

      if (ord.eq.'exp') then
         aieps = -2.d0*log(eps) -1.5d0+2.d0*eps-0.5d0*eps**2
! with reduced to born momenta...
         aopi = alpha/pi
         s = 4.d0*pin1b(0)*pin1b(0)

! for gg not-reduced center of mass..
ccccc         if (fs.eq.'gg') s = ecms**2

         p1mod = sqrt(tridot(p1b,p1b))
         c = p1b(3)/p1mod
         t = -s/2.d0*(1.d0 - c)
         u = -s/2.d0*(1.d0 + c)

         ame2mu = 1.d0/ame/ame
         sdif = bornme*phsp2b/flux2b

** VP STUFF
*** changed for no VP in 1-loop
         sdif = bornmenovp*phsp2b/flux2b
***
         
         if (fs.eq.'ee'.or.fs.eq.'mm') then
            scale = s*t/u
            lcoll = log(scale*ame2mu) - 1.d0
            arg = 2.d0 * aopi * lcoll * aieps

**** NEW ***************
            charges(1) = -1.d0
            charges(2) = -1.d0
            charges(3) = -1.d0
            charges(4) = -1.d0
            masses(1)  = ame
            masses(2)  = ame
            masses(3)  = ame
            masses(4)  = ame
            if (fs.eq.'mm') then
               masses(3)  = ammu
               masses(4)  = ammu
            endif
            lcoll1 = eikonalintegral(pin1b,pin2b,p1b,p2b,charges,masses)
            arg2   = 2.d0*aopi * lcoll1 * aieps 
            arg = arg2
**************************
** VP STUFF
ccc   changed for no VP in virtual
            iarunbck = iarun ! added
            iarun = 0        ! added ! ON OFF
            svpureoal = aintsv_new(pin1b,pin2b,p1b,p2b,eps)/sdif-1.d0
            iarun = iarunbck ! added
         elseif (fs.eq.'gg') then
            scale = s
            lcoll = log(scale*ame2mu) - 1.d0
            arg = 1.d0 * aopi * lcoll * aieps
            svpureoal = aintsvgg(s,c,eps)/sdif - 1.d0
         endif
         deltasv = svpureoal + arg      !!! NB

*** NEW SFF, formula 2.23 PhD. Thesis. in the limit dot(pi,pj) >> m**2
         dp1p2 = log(2.d0*dot(pin1,pin2)*ame2mu)
         dp1p3 = log(2.d0*dot(pin1,p1)*ame2mu)
         dp1p4 = log(2.d0*dot(pin1,p2)*ame2mu)
         dp2p3 = log(2.d0*dot(pin2,p1)*ame2mu)
         dp2p4 = log(2.d0*dot(pin2,p2)*ame2mu)
         dp3p4 = log(2.d0*dot(p1,p2)*ame2mu)

cc         dp1p2 = log(2.d0*dot(pin1b,pin2b)*ame2mu)
cc         dp1p3 = log(2.d0*dot(pin1b,p1b)*ame2mu)
cc         dp1p4 = log(2.d0*dot(pin1b,p2b)*ame2mu)
cc         dp2p3 = log(2.d0*dot(pin2b,p1b)*ame2mu)
cc         dp2p4 = log(2.d0*dot(pin2b,p2b)*ame2mu)
cc         dp3p4 = log(2.d0*dot(p1b,p2b)*ame2mu)

         if (fs.eq.'ee'.or.fs.eq.'mm') then
            lcollnew = dp1p2+dp1p3-dp1p4-dp2p3+dp2p4+dp3p4
            lcollnew = 2.d0 * lcollnew - 4.d0
            lcollnew = lcollnew/4.d0
            arg = 2.d0 * aopi * lcollnew * aieps

**** NEW ***************
            charges(1) = -1.d0
            charges(2) = -1.d0
            charges(3) = -1.d0
            charges(4) = -1.d0
            masses(1)  = ame
            masses(2)  = ame
            masses(3)  = ame
            masses(4)  = ame
            if (fs.eq.'mm') then
               masses(3)  = ammu
               masses(4)  = ammu
            endif
            lcoll1 = eikonalintegral(pin1,pin2,p1,p2,charges,masses)
            arg2   = 2.d0*aopi * lcoll1 * aieps 
            arg = arg2
**************************

         elseif (fs.eq.'gg') then
            lcollnew = dp1p2
            lcollnew = lcollnew - 1.d0
            arg = 1.d0 * aopi * lcollnew * aieps
         endif

!         deltasv = svpureoal + arg      !!! NB
         ffs = dexp(-arg)
         if (model.eq.'ps'.or.model.eq.'eikonal') then 
            deltasv = 0.d0
         endif
         sv = ffs * (1.d0 + deltasv)
         dsv = deltasv
         return
      endif

      if (ord.eq.'alpha2') then
         aieps = -2.d0*log(eps) -1.5d0+2.d0*eps-0.5d0*eps**2
         s = (pin1b(0) + pin2b(0))**2
         p1mod = sqrt(tridot(p1b,p1b))
         c = p1b(3)/p1mod

         t = -s/2.d0*(1.d0 - c)
         u = -s/2.d0*(1.d0 + c)

         scale = s*t/u
         lcoll = log(scale/ame**2) - 1.d0

         arg = 2.d0*alpha/pi * lcoll * aieps

         sdif = bornme*phsp2b/flux2b

         svpureoal = aintsv(s,t,eps)/sdif - 1.d0
         deltasv = svpureoal + arg      !!! NB

*** NEW SFF, formula 2.23 PhD. Thesis. in the limit dot(pi,pj) >> m**2
         dp1p2 = log(2.d0*dot(pin1,pin2)/ame**2)
         dp1p3 = log(2.d0*dot(pin1,p1)/ame**2)
         dp1p4 = log(2.d0*dot(pin1,p2)/ame**2)
         dp2p3 = log(2.d0*dot(pin2,p1)/ame**2)
         dp2p4 = log(2.d0*dot(pin2,p2)/ame**2)
         dp3p4 = log(2.d0*dot(p1,p2)/ame**2)
         lcollnew = dp1p2 + dp1p3 - dp1p4 - dp2p3 + dp2p4 + dp3p4
         lcollnew = 2.d0 * lcollnew - 4.d0
         lcollnew = lcollnew/4.d0
         arg = 2.d0*alpha/pi * lcollnew * aieps

!         deltasv = svpureoal + arg      !!! NB

         ffs = dexp(-arg)
         if (model.eq.'ps') deltasv = 0.d0

         if (ng.le.2) then
            sv = 1.d0
         endif
         if (ng.le.1) then
            sv = sv - arg + deltasv
         endif
         if (ng.eq.0) then
            sv = sv + arg**2/2.d0 - arg*deltasv
         endif
         return
      endif
      return
      end
