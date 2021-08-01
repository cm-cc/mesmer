* COMPLEX VACUUM POLARIZATION: ALPHA(0)->ALPHA(0)*VPOL(Q2)
      FUNCTION VPOLC(Q2) 
      IMPLICIT double precision (A-H,O-Z)
      double precision QIN4,ST24,DER4,ERRDER4,DEG4,ERRDEG4
      DIMENSION AMASSES(9)
      common/hadoff/ihadoff
      common/parameters/ame,ammu,convfac,alpha,pi     
      common/tautopmasses/taumass,topmass
      common/ialpharunning/iarun
      common/teubner/iteubn
      integer inskvp
      common/nskvp/inskvp
      common/resonances/ires,ionDHAD
      common/vpolcomc/ifirst
      common/jpsiparameters/amjpsi,gtjpsi,gejpsi,effcjpsi
      common/jpsiparameters2/amjpsi2s,gtjpsi2s,gejpsi2s,effcjpsi2s
      common/jpsiparameters3/amjpsi3s,gtjpsi3s,gejpsi3s,effcjpsi3s
      common/jpsiparameters4/amjpsi4s,gtjpsi4s,gejpsi4s,effcjpsi4s
      common/jpsiparameters5/amjpsi5s,gtjpsi5s,gejpsi5s,effcjpsi5s

      common/tmintmax/tmin,tmax

*** from invariants.h      
      double complex vpdah
      double precision vpdal
      double precision vpvect(5)
      common/vp_dalphas/vpdal,vpdah,vpvect
      integer ionlybubbles,ibubbleonv,ibubbleonboxsub
      common/ionlybubblescmn/ionlybubbles,ibubbleonv,ibubbleonboxsub
***      
      data ifirst /0/

      double complex vpolcnsk,vpolchlmnt
      external vpolcnsk,vpolchlmnt

      double complex vpolc,deltahad,im,cder
      parameter (nres=15)
      double precision AM(nres),GT(nres),GE(nres)
      common/resonancesparams/am,gt,ge
*     AM=HADRON MASSES, GT=TOTAL WIDTH, GE=PARTIAL E+E- OR Mu+Mu- WIDTH.
* 1=RHO, 2=OMEGA, 3=PHI, 4=J/PSI, 5=PSI(2S), 6=PSI(3770), 7=PSI(4040),
* 8=PSI(4160), 9=PSI(4415), 10=Y, 11=Y(2S), 12=Y(3S), 13=Y(4S), 14=Y(10860),
* 15=Y(11020).
      double precision EAM(nres),EGT(nres),EGE(nres)
      common/resonancesparamserrors/eam,egt,ege
** errors on the previous parameters, from PDG
      
****  TO BE UPDATED (e.g. from PDG)!!
      DATA AM/ ! GeV
     >     0.769D0,    ! 1.  rho
     >     0.78265D0,  ! 2.  omega
     >     1.019455D0, ! 3.  phi
     >     3.096916D0, ! 4.  J/Psi
     >     3.686109D0, ! 5.  Psi(2S)
     >     3.77315D0,  ! 6.  Psi(3770)
     >     4.039D0,    ! 7.  Psi(4040)
     >     4.191D0,    ! 8.  Psi(4160)
     >     4.421D0,    ! 9.  Psi(4415)
     >     9.4603D0,   ! 10. Y
     >     10.02326D0, ! 11. Y(2S)
     >     10.35520,   ! 12. Y(3S)
     >     10.5794D0,  ! 13. Y(4S)
     >     10.876D0,   ! 14. Y(10860)
     >     11.019D0/   ! 15. Y(11020)
      DATA GT/ ! GeV
     >     0.149D0,   ! 1.  rho
     >     8.49D-3,   ! 2.  omega
     >     4.266D-3,  ! 3.  phi
     >     92.9D-6,   ! 4.  J/Psi
     >     299.D-6,   ! 5.  Psi(2S)
     >     27.2D-3,   ! 6.  Psi(3770)
     >     80.D-3,    ! 7.  Psi(4040)
     >     70.D-3,    ! 8.  Psi(4160)
     >     62.D-3,    ! 9.  Psi(4415)
     >     54.02D-6,  ! 10. Y
     >     31.98D-6,  ! 11. Y(2S)
     >     20.32D-6,  ! 12. Y(3S)
     >     20.5D-3,   ! 13. Y(4S)
     >     55.D-3,    ! 14. Y(10860)
     >     79.D-3/    ! 15. Y(11020)
      DATA GE/ ! GeV
     >     7.04D-6,   ! 1.  rho
     >     0.60D-6,   ! 2.  omega
     >     1.27D-6,   ! 3.  phi
     >     5.55D-6,   ! 4.  J/Psi
     >     2.36D-6,   ! 5.  Psi(2S)
     >     0.262D-6,  ! 6.  Psi(3770)
     >     0.86D-6,   ! 7.  Psi(4040)
     >     0.48D-6,   ! 8.  Psi(4160)
     >     0.58D-6,   ! 9.  Psi(4415)
     >     1.34D-6,   ! 10. Y
     >     0.612D-6,  ! 11. Y(2S)
     >     0.443D-6,  ! 12. Y(3S)
     >     0.272D-6,  ! 13. Y(4S)
     >     0.31D-6,   ! 14. Y(10860)
     >     0.130D-6/  ! 15. Y(11020)
*** ERRORS ON THE PREVIOUS PARAMETERS, from PDG
      DATA EAM/                  ! GeV
     >     0.0009D0,   ! 1.  rho
     >     0.00012D0,  ! 2.  omega
     >     0.000019D0, ! 3.  phi
     >     0.000011D0, ! 4.  J/Psi
     >     0.000013D0, ! 5.  Psi(2S)
     >     0.00033D0,  ! 6.  Psi(3770)
     >     0.001D0,    ! 7.  Psi(4040)
     >     0.005D0,    ! 8.  Psi(4160)
     >     0.004D0,    ! 9.  Psi(4415)
     >     0.00026D0,  ! 10. Y
     >     0.00031D0,  ! 11. Y(2S)
     >     0.0005D0,   ! 12. Y(3S)
     >     0.0012D0,   ! 13. Y(4S)
     >     0.011D0,    ! 14. Y(10860)
     >     0.008D0/    ! 15. Y(11020)
      DATA EGT/ ! GeV
     >     0.001D0,  ! 1.  rho
     >     0.08D-3,  ! 2.  omega
     >     0.031D-3, ! 3.  phi
     >     2.8D-6,   ! 4.  J/Psi
     >     8.0D-6,   ! 5.  Psi(2S)
     >     1.0D-3,   ! 6.  Psi(3770)
     >    10.D-3,    ! 7.  Psi(4040)
     >    10.D-3,    ! 8.  Psi(4160)
     >    20.D-3,    ! 9.  Psi(4415)
     >     1.25D-6,  ! 10. Y
     >     2.63D-6,  ! 11. Y(2S)
     >     1.85D-6,  ! 12. Y(3S)
     >     2.5D-3,   ! 13. Y(4S)
     >    28.D-3,    ! 14. Y(10860)
     >    16.D-3/    ! 15. Y(11020)
      DATA EGE/ ! GeV
     >     0.06D-6,  ! 1.  rho
     >     0.02D-6,  ! 2.  omega
     >     0.04D-6,  ! 3.  phi
     >     0.16D-6,  ! 4.  J/Psi
     >     0.04D-6,  ! 5.  Psi(2S)
     >     0.018D-6, ! 6.  Psi(3770)
     >     0.07D-6,  ! 7.  Psi(4040)
     >     0.22D-6,  ! 8.  Psi(4160)
     >     0.07D-6,  ! 9.  Psi(4415)
     >     0.018D-6, ! 10. Y
     >     0.011D-6, ! 11. Y(2S)
     >     0.008D-6, ! 12. Y(3S)
     >     0.029D-6, ! 13. Y(4S)
     >     0.07D-6,  ! 14. Y(10860)
     >     0.03D-6/  ! 15. Y(11020)

*
* E,MU,TAU,TOP MASSES: 
*
      common/massesVP/amasses
c     DATA AMASSES /0.510999D-3,.1056584D0,1.777D0,175.6D0,      
!     Fael tau
      DATA AMASSES /0.510999D-3,.1056584D0,1.77682D0,175.6D0,
     . 0.066d0,! u
     . 0.066d0,! d
     . 1.2d0,  ! c
     . 0.150d0,! s 
     . 4.6d0  /! b
*     
      double precision dalo(3),danlo(3),dannlo(3)
      common/dalonlonnlo/dalo,danlo,dannlo
      double precision dalotop,danlotop
      common/dalonlotop/dalotop,danlotop

      if (ifirst.eq.0) then
! filling the common for importance sampling...

*** taken from the common ***     
         amasses(1) = ame
         amasses(2) = ammu
         amasses(3) = taumass
         amasses(4) = topmass
*****************************         
         ifirst = 1

         amjpsi   = am(4)
         gtjpsi   = gt(4)
         gejpsi   = ge(4)
         effcjpsi = 3.d0*gejpsi/amjpsi/alpha

         amjpsi2s   = am(5)
         gtjpsi2s   = gt(5)
         gejpsi2s   = ge(5)
         effcjpsi2s = 3.d0*gejpsi2s/amjpsi2s/alpha

         amjpsi3s   = am(10)
         gtjpsi3s   = gt(10)
         gejpsi3s   = ge(10)
         effcjpsi3s = 3.d0*gejpsi3s/amjpsi3s/alpha
         amjpsi4s   = am(11)
         gtjpsi4s   = gt(11)
         gejpsi4s   = ge(11)
         effcjpsi4s = 3.d0*gejpsi4s/amjpsi4s/alpha
         amjpsi5s   = am(12)
         gtjpsi5s   = gt(12)
         gejpsi5s   = ge(12)
         effcjpsi5s = 3.d0*gejpsi5s/amjpsi5s/alpha

         ires = 0

         ionDHAD = 1
         if (ires.eq.1) ionDHAD = 0

      endif

      vpolc    = (1.d0,0.d0)
      vpdal    = 0.d0
      vpdah    = 0.d0
      vpvect = 0.d0

      if (iarun.eq.0) then
         vpolc = 1.d0
         return
      endif

      if (inskvp.eq.1) then
         vpolc = VPOLCnsk(Q2)

         datotnsk = 1.d0 - 1.d0/vpolc

         SOMMA=0.D0
         DO I=1,4
            vpvect(i) = SUMMA(AMASSES(I),Q2,I)
            SOMMA = SOMMA + vpvect(i)
         ENDDO
         dalphaQED = ALPHA/PI*SOMMA !+(alpha/pi)**3*(log(-q2/ame/ame)-1.d0)**2
         vpvect  = vpvect*alpha/pi
         
         vpdal = dalphaQED
         if (ihadoff.eq.0) vpdah = datotnsk - vpdal

         vpvect(5) = vpdah
         
         vpolc = 1.d0/(1.d0-vpdal-vpdah)
         
         return
      endif
         
      SOMMA=0.D0
      DO I=1,4 !3 only leptons
         vpvect(i) = SUMMA(AMASSES(I),Q2,I)
         SOMMA = SOMMA + vpvect(i)
      ENDDO
      dalphaQED = ALPHA/PI*SOMMA  !+(alpha/pi)**3*(log(-q2/ame/ame)-1.d0)**2
      vpvect(1:4)  = vpvect(1:4) * alpha/pi
     
c***********************************************************************************
c!! adding leptonic nlo and nnlo contributions to  \Delta\alpha. Top only NLO (tiny)
c! this just to load the commons dalonlonnlo and dalonlotop      
c      ateubnlonnlo = alf_lep_carlo(1,q2) + alf_lep_carlo(2,q2)
c     .     + alf_lep_carlo(3,q2)+alf_top_carlo(q2)
c*
c      dalphaQED = dalphaQED + sum(danlo) + sum(dannlo) + danlotop
c***********************************************************************************
      vpdal     = dalphaQED
      
      if (ihadoff.eq.1) then
         vpolc = 1.d0/(1.d0 - dalphaQED)
         return
      endif
      
      qin = q2/abs(q2) * sqrt(abs(q2))
      der = 0.d0
      if (iteubn.eq.0) then
         call wrap_Fred_hadVP(qin,der,errder)                        
      elseif (iteubn.eq.1) then
c     old            call d_alpha(dabs(qin),das,ddas,dat,ddat)
         call d_alpha_reso_carlo(dabs(qin),das,ddas,dat,ddat)            
         if (qin.lt.0.d0) then
            der = alpha* das
         else
            der = alpha* dat
         endif
      endif
      cder  =  der
      vpdah = cder
      vpolc = 1.d0/(1.d0 - dalphaQED - cder)
      vpvect(5) = vpdah
      RETURN
      END
********************************************************
      subroutine wrap_Fred_hadVP(qin,der,errder)
      implicit double precision (a-h,o-z)
      parameter (ST2 = 0.2322d0)
      integer ifredv
      common/fredversion/ifredv
      if (ifredv.eq.12) then
         CALL DHADR5n12(QIN,ST2,DER,ERRDER,DEG,ERRDEG)
      elseif (ifredv.eq.17) then
         call dhadr5n17cmcc(qin,st2,der,errdersta,errdersys,deg,
     .        errdegsta,errdegsys) ! it's the '17 version
         errder = sqrt(errdersta**2+errdersys**2)
      elseif (ifredv.eq.19) then
         call dhadr5x19cmcc(qin,st2,der,errdersta,errdersys,deg,
     .        errdegsta,errdegsys) ! it's the '19 version
         errder = sqrt(errdersta**2+errdersys**2)
      else
         print*,'Wrong version for Jegerlehner''s VP. Stopping.'
         stop
      endif
      return
      end
******************************************************
* COMPLEX VACUUM POLARIZATION from NSBSK DATA: ALPHA(0)->ALPHA(0)*VPOL(Q2)
* from    http://cmd.inp.nsk.su/~ignatov/vpl/
      FUNCTION VPOLCnsk(Q2) 
      implicit double precision (A-H,L,O-Z)
c      parameter (np=7186)
      parameter (np=6883)
      dimension s(0:np-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)
      double complex vpolcnsk,deltahad,im
      common/vpcnskv2/s,pres,pre,pim,epres,epre,epim,cv,im,ifirst
      common/nskimid/imidguess
      common/tmintmax/tmin,tmax
      data ifirst/0/
***      
      if (ifirst.eq.0) then
         im = (0.d0,1.d0)
         open(45,file='vpol_novosibirsk_v2.dat',status='old')
         do k = 0,np-1
            read(45,*)s(k),pres(k),pre(k),pim(k),epres(k),epre(k),
     .           epim(k),cv(k)
         enddo
         close(45)
         ifirst=1
*** guessing where to start the binary search....         
         k = 0
         aq2 = -tmax*2d0
         do while(aq2.gt.s(k))
            k = k + 1         
         enddo
         k = k - 1
         imidguess = k ! this is correct for born s-channel and no energy spread
***************
      endif

c      k   = 0
c      sk  = 0.d0
c      aq2 = abs(q2)
c      do while(aq2.gt.s(k))
c         k = k + 1         
c      enddo
c      k = k - 1
******** binary search ****
      aq2     = abs(q2)
      istart  = 0
      ifinish = np-1
c      imid = (istart + ifinish)/2
      imid = imidguess
      do while(.not.(s(imid+1).gt.aq2.and.s(imid).le.aq2))
         if (aq2.gt.s(imid)) then
            istart  = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
******************
      if (q2.gt.0.d0) then
         vpre = pre(k) + (pre(k+1)-pre(k))/(s(k+1)-s(k)) * (q2-s(k))
         vpim = pim(k) + (pim(k+1)-pim(k))/(s(k+1)-s(k)) * (q2-s(k))

         evpre = epre(k) + (epre(k+1)-epre(k))/(s(k+1)-s(k)) * (q2-s(k))
         evpim = epim(k) + (epim(k+1)-epim(k))/(s(k+1)-s(k)) * (q2-s(k))

c         rho = cv(k) + (cv(k+1)-cv(k))/(s(k+1)-s(k)) * (q2-s(k))
c         cov = rho * evpre*evpim
      else
         vpre = pres(k) + (pres(k+1)-pres(k))/(s(k+1)-s(k))*(aq2-s(k))
         vpim = 0.d0

         evpre = epres(k)+(epres(k+1)-epres(k))/(s(k+1)-s(k))*(aq2-s(k))
         evpim = 0.d0

c         rho = 0.d0
c         cov = 0.d0
      endif
*********************************************************************************************
      VPOLCnsk = 1.d0/(1.d0  - vpre - im*vpim)
      return
      end
*==============================================================================================
*     LEPTONIC AND TOP CONTRIBUTION TO VACUUM POLARIZATION
** with correct expansion for |q2| <<< m2!!
      FUNCTION SUMMA(AM,Q2,I)
      IMPLICIT double precision (A-H,O-Z)
      double precision NC(9),QF2(9)
      common/ncqf2/nc,qf2
      data nc  /1.d0,1.d0,1.d0,
     . 3.d0, ! top, which is perturbative
     . 3.d0, !u
     . 3.d0, !d
     . 3.d0, !c
     . 3.d0, !s
     . 3.d0 /!b
      data qf2 /1.d0,1.d0,1.d0,0.44444444444444444444d0,
     . 0.44444444444444444444d0, !u
     . 0.11111111111111111111d0, !d
     . 0.44444444444444444444d0, !c
     . 0.11111111111111111111d0, !s
     . 0.11111111111111111111d0 /!b
* NC AND QF ARE COLOR FACTOR (1 FOR LEPTONS, 3 FOR TOP) AND CHARGE

      if (Q2.eq.0.d0) then
         summa = 0.d0
         return
      endif      
      
      AM2=AM**2

      z6lim = 1d-8

      IF (Q2.GT.0.D0.AND.Q2.LT.(4.D0*AM2)) THEN
         x = 4.D0*AM2/Q2
         z = sqrt(1.d0/x)
         z2 = z*z
         z4 = z2*z2
         z6 = z2*z4

         if (z6.gt.z6lim) then
            SQ = SQRT(1d0 - z2)
            s1 = -5.D0/3.D0 - 1.d0/z2
            s2 = (1.d0 + z2- 2.D0*z4)/z/z2*1.d0/SQ*ATAN(z/SQ)
            s1ps2 = s1 + s2
         else
*     with wxmaxima
* s2 = taylor((1+z^2-2*z^4)/sqrt(1-z^2)/z^3*atan(z/sqrt(1-z^2)), z , 0 , 20);
            s1ps2 = -z2*
     .       (4.d0/5.d0+12.d0/35.d0*z2+64.d0/315.d0*z4+32.d0/231.d0*z6)
         endif
         SUMMA=1.d0/3.d0*NC(I)*QF2(I)*s1ps2
      ELSE
           x = -4.D0*AM2/Q2
           if (x.gt.0.d0) then ! space-like
              z  = sqrt(1.d0/x)
              z2 = z*z
              z4 = z2*z2
              z6 = z2*z4
              
              if (z6.gt.z6lim) then
                SQ     = SQRT(1.D0 + z2)
                ARGLOG = ABS( (z - sq)/(z+sq) )              
                s1 = -5.D0/3.D0 + 1.d0/z2
                s2 = 0.5d0*(1.d0 - z2 - 2.D0*z4)/SQ*LOG(ARGLOG)/z2/z
                
                s1ps2 = s1 + s2
              else
* with wxmaxima
* s2 = taylor(log( abs ((z - sqrt(1+z^2))/(z + sqrt(1+z^2)) ))/sqrt(1+z^2)*(1-z^2-2*z^4)/2/z^3,z,0,20); 
                 s1ps2 = z2*
     .        (4.d0/5.d0-12.d0/35.d0*z2+64.d0/315.d0*z4-32.d0/231.d0*z6)
              endif
           else ! for x > 0 ie q2 > 4m^2
              x = -x
              z  = sqrt(1.d0/x)
              z2 = z*z
              z4 = z2*z2

              SQ     = SQRT(z2-1.d0)
              ARGLOG = ABS((z - sq)/(z+sq) )            
              s1 = -5.D0/3.D0 - 1.d0/z2
              s2 =  0.5d0*(1.d0 + z2 - 2.D0*z4)/SQ*LOG(ARGLOG)/z2/z
              s1ps2 = s1 + s2
           endif
           SUMMA = 1.d0/3.d0 * NC(I)*QF2(I)*s1ps2
       ENDIF
       RETURN
       END

*****************
      function wrapRhad(s)
      implicit double precision (a-h,o-z)

      parameter (np = 3885)     ! number of points in Rhad-scan.dat
      double precision rhens(np),rhad(np,2)
      common/rhadscans/rhens,rhad

      double precision wrapRhad,s,eR,sin
      double precision en,rout
      integer ifirst
      common/wraprhadcmn/ifirst,i
      data ifirst /0/

      integer iteubn
      common/teubner/iteubn
      
      if (ifirst.eq.0) then
         call r_intpl_reso(1.d0,rrrr,errrr) ! just to fill the common wKNT...

         i = 1    ! which rhad to choose [1 or 2 for fred or knt]
         if (iteubn.eq.1) i = 2

         print*,'Loading Rhad data, obtained from'//
     .        ' alphaQEDc19 and KNT 3.0.1. Stored in Rhad-scan.dat.'
         if (i.eq.1) print*,'Using alphaQEDc19...'
         if (i.eq.2) print*,'Using KNT...'
         open(300,file='Rhad-scan.dat',status='old')
         do k = 1,np
            read(300,*)ii,rhens(k),rhad(k,1),rhad(k,2) ! ie rhad(:,1) = rfred, rhad(:,2) = rknt
         enddo
         close(300)
         ifirst = 1
         print*,'... done!'
      endif
      
      en = sqrt(s)

      if (en.ge.rhens(np)) then
         wrapRhad = rhad(np,i)
         return
      endif
      if (en.le.rhens(1)) then
         wrapRhad = rhad(1,i)
         return
      endif

***   binary search
      istart  = 1
      ifinish = np
      imid = (ifinish+istart)/2
      do while(.not.(rhens(imid+1).gt.en
     .     .and.rhens(imid).le.en))
         if (en.gt.rhens(imid)) then
            istart = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid

** not binary search...      
c      k = 1
c      do while(en.gt.rhens(k))
c         k = k + 1
c      enddo
c      k = k - 1      
****
      wrapRhad = (rhad(k+1,i)- rhad(k,i))
     .     /(rhens(k+1) - rhens(k)) * (en - rhens(k)) + rhad(k,i)
      return
      end
*****************
      function wrapRhadv1(s)
      double precision wrapRhadv1,s,eR,sin
      double precision en,rout
      en = sqrt(s)
      call r_intpl_reso(en,wrapRhadv1,eR)
c      call wrapRfred(s,wrapRhad)
      return
      end
***********************************************************
c$$$      subroutine wrapRfred(s,Rout)
c$$$      implicit double precision (a-h,o-z)      
c$$$      double precision Rout
c$$$      double complex cone,calp,cvpt_new,cvpt_hig,cvpt_low,alphac,alphah
c$$$     &     ,alphal,alphact,alphaht,alphalt,calept1,cahadr1,caltop1
c$$$
c$$$      
c$$$      double precision myalpha,mypi
c$$$      common/carlocommons/myalpha,mypi
c$$$      mypi = 4.d0*atan(1.d0)
c$$$      myalpha = 1.d0/137.03599907430637d0
c$$$
c$$$      call alphaQEDcxcmcc(s,cvpt_new,cvpt_hig,cvpt_low,alphac,
c$$$     &    alphah,alphal,alphact,alphaht,alphalt,calept1,cahadr1,caltop1)
c$$$
c$$$      Rout = (0.d0,1.d0)*cahadr1*mypi/myalpha
c$$$
c$$$      return
c$$$      end
c$$$      
c$$$*;;;;;;;;;;;;;;;;;;;;;;;;;;;;; -*- Mode: Fortran -*- ;;;;;;;;;;;;;;;;;;;;;;;;;;;
c$$$*;; alphaQEDc_sub.f --- 
c$$$*;; Author          : Fred Jegerlehner
c$$$*;; Created On      : Wed Jan 11 03:26:27 2006
c$$$*;; Last Modified By: Friedrich Jegerlehner
c$$$*;; Last Modified On: Tue Jul 18 12:57:38 2017
c$$$*;; RCS: $Id$
c$$$*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
c$$$*;; Copyright (c) 2010 Fred Jegerlehner
c$$$*;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
c$$$*;; 
c$$$      SUBROUTINE alphaQEDcxcmcc(s,cvpt_new,cvpt_hig,cvpt_low,alphac,
c$$$     &    alphah,alphal,alphact,alphaht,alphalt,calept1,cahadr1,caltop1)
c$$$C Link files dggvapx,dhadr5n17,cggvapx,chadr5n17
c$$$C Calculating the running QED coupling alpha(E) as a complex function
c$$$C Leptons: 1--loop, 2--loop exact, 3--loop in high energy approximation
c$$$C Quarks form routine dhadr5n including effective s--channel shifts in low energy range
c$$$C r1....,r2....,r3.... : one-loop, two-loop, three-loop result
c$$$C  ..real: real-part, ..imag: imaginary part;
c$$$C ......l: light fermions only (high energy approximation for 3--loops) 
c$$$C dallepQED1n=r1real,dallepQED2n=r2real,dallepQED3l=r3reall
c$$$C dallepQED3l also returns one-loop and two-loop high energy (light fermion) 
c$$$C approximations 
c$$$c 21/09/2009 common pQCD parameter setting via constants_qcd.f      
c$$$c FJ 28/01/2012 recalculated Nall=2055 --> 2177
c$$$c 29/01/2012 common resonance parameters setting via resonances_dat.f
c$$$c 21/09/2012 new option for W contribution: iOSW=0 strict MSbar, iOSW=1 including 
c$$$c            high energy asymptotic constant of the OS W contribution at s>>MW2
c$$$c 10/07/2017 providing statistical and systematic hadronic errors separately now
c$$$c            allows to set up the covariance matrix
c$$$      implicit none
c$$$      double precision s,null
c$$$c      double precision dalept,dahadr,daltop,Dalphaweak1MSb
c$$$      double complex cggvapx,cerror,cerrorsta,cerrorsys,
c$$$     &     calept,cahadr,caltop,cDalphaweak1MSb
c$$$      double complex cone,calp,cvpt_new,cvpt_hig,cvpt_low,alphac,alphah
c$$$     &     ,alphal,alphact,alphaht,alphalt,calept1,cahadr1,caltop1
c$$$      external cggvapx
c$$$      integer ifirst
c$$$      common/complexdalphacarlo/ifirst
c$$$      data ifirst /0/
c$$$********************* part to be included and edited in main program **************
c$$$c check cuts CHPTCUT,EC,ECUT   ! specify where to use CHPT and pQCD in constants_qcd
c$$$c decommment next line in main
c$$$      double precision pi1,ALINP,ERRAL,EINP,MTOP1 ! get QCD parameters as set in constants_qcd
c$$$c QCD parameters required to calculare imaginary part (R) and top contribution 
c$$$      include 'alphaQEDc19/common.h'
c$$$c     decomment next 3 lines in main
c$$$      double precision st2,als,mtop,sin2ell
c$$$      integer Ismooth
c$$$      common/var1/pi1,ALINP,ERRAL,EINP,MTOP1
c$$$      common /parm/st2,als,mtop
c$$$      COMMON/SMOO/Ismooth
c$$$c decomment next line if needed (provided also by cggvap)
c$$$c      common /resu/dalept,dahadr,daltop,Dalphaweak1MSb   ! real results
c$$$      common /cres/calept,cahadr,caltop,cDalphaweak1MSb  ! complex results
c$$$c decomment next 2 lines in main
c$$$      if (ifirst.eq.0) then
c$$$         LEPTONflag='had'       ! choose had,lep,top,ele,muo,tau
c$$$         iLEP=LFLAG(LEPTONflag)
c$$$*         
c$$$         call constants()
c$$$         call constants_qcd()
c$$$         sin2ell=0.23153        ! pm 0.00016 LEPEEWG Phys Rep 427 (2006) 257
c$$$c decommment next line and change default as appropriate in main      
c$$$         st2=0.23153d0          ! Reference value for weak mixing parameter
c$$$c decomment next 2 lines in main
c$$$         als=ALINP              ! alpha strong
c$$$         mtop=MTOP1             ! top quark mass
c$$$         Ismooth=0              ! flag for imaginary part: 0= R(s) data; 1= R(s) fits
c$$$c following two entried set as default in chad5n09.f
c$$$         IRESON=1               ! include narrow resonances in fits
c$$$         iresonances=1          ! include narrow resonances in data sets
c$$$c decomment next line in main
c$$$         call Rdata()
c$$$c decomment next line in main
c$$$         call resonances_renocomplex()
c$$$c      call resonances_renoreal()
c$$$C LEPTONflag=all,had,lep,ele,muo,tau -> iLEP=-3,-2,-1,1,2,3
c$$$C Default: 
c$$$c       LEPTONflag='all'
c$$$c       iLEP  = -3  ! for sum of leptons + quarks              
c$$$C decomment in main the following 2 lines for changing the default
c$$$c ***********************************************************************
c$$$c
c$$$      endif
c$$$      null=0.d0
c$$$      cone=DCMPLX(1.d0,null)
c$$$      calp=DCMPLX(alp,null)
c$$$      cvpt_new= cggvapx(s,cerror,cerrorsta,cerrorsys)           ! central value
c$$$      cvpt_hig= cvpt_new+cerror            !    ...        + total error 
c$$$      cvpt_low= cvpt_new-cerror            !    ...        - total error
c$$$      alphac  =calp/(cone-cvpt_new)        ! alpha no top
c$$$      alphah  =calp/(cone-cvpt_hig)
c$$$      alphal  =calp/(cone-cvpt_low)
c$$$      alphact =calp/(cone-cvpt_new-caltop-cDalphaweak1MSb) ! alpha incl top + W
c$$$      alphaht =calp/(cone-cvpt_hig-caltop-cDalphaweak1MSb)
c$$$      alphalt =calp/(cone-cvpt_low-caltop-cDalphaweak1MSb)
c$$$      calept1 =calept                      ! lepton VP
c$$$      cahadr1 =cahadr                      ! hadronic VP 5 flavors
c$$$      caltop1 =caltop                      ! top contribution                 
c$$$      return
c$$$      END
