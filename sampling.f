      subroutine getE4(E4max,E4min,s,m32,m42,E4,w)
      implicit double precision (a-h,m,o-z)
      double precision csi1(1)
      common/E4get/sbefore,dm2,sqs,anmin,ifirst
      data sbefore /-1d20/

      if (abs(s-sbefore).gt.1d-10) then
         sbefore = s
         dm2 = m42 - m32
         sqs = sqrt(s)
         anmin = -0.5d0/sqs*log(s-2.d0*sqs*E4min+dm2)
      endif
      anorm = -0.5d0/sqs*log(s-2.d0*sqs*E4max+dm2)-anmin

      call getrnd(csi1,1)
      E4 = s+dm2 - exp(-2.d0*sqs*(anorm*csi1(1)+anmin))
      E4 = E4*0.5d0/sqs

      w = (s - 2.d0*sqs*E4 + dm2)*anorm
      
      return
      end
********************************************************
      subroutine gets3(s3max,s3min,s,m32,m42,s3,w)
      implicit double precision (a-h,m,o-z)
      double precision csi1(1),csi2(2)
c      goto 114
***   as 1/(s3-m32)      
      call getrnd(csi1,1)
c      s3 = (s3max - s3min) * csi1(1) + s3min
c      w  = (s3max - s3min) * w
      anom = log((s3max-m32)/(s3min-m32))
      s3   = m32 + (s3min-m32)*exp(anom*csi1(1))
      w    = anom * (s3-m32)
      return

******   as 1/(s3)**4 + 1/(s-s3)**4
 110  continue ! no decisamente no
c      s3 = (s3max - s3min) * csi1(1) + s3min
c      w  = (s3max - s3min) * w
      anom1 = 1.d0/(s3min)**2 - 1.d0/(s3max)**2
      anom1 = anom1/2.d0

      anom2 = 1.d0/(s-s3max)**3 - 1.d0/(s-s3min)**3
      anom2 = anom2/3.d0

      call getrnd(csi1,1)
      if (csi1(1).lt.anom1/(anom1+anom2)) then
         call getrnd(csi1,1)

         aa = 1.d0/(s3min)**2 - 2.d0*anom1*csi1(1)
         s3 = 1.d0/aa**(1.d0/2.d0)
         
      else
         call getrnd(csi1,1)
         aa = anom2*3.d0*csi1(1)+1.d0/(s-s3min)**3
         s3 = s - 1.d0/aa**(1.d0/3.d0)
      endif
      
      w = (1.d0/(s3)**3 + 1.d0/(s-s3)**4)/(anom1+anom2)
      w = 1.d0/w
      return
      
***   as p1/(s3-m32)  + p2/(s3 + s - m42) = p1/(s3-m32)  + p2/(s3 -(m42 - s))
 111  anom1 = log((s3max-m32)/(s3min-m32))
      anom2 = log((s3max-(m42-s))/(s3min-(m42-s)))

      p1 = 0.95d0
      p2 = 1.d0 - p1

      call getrnd(csi2,2)
      if (csi2(2).lt.p1) then
         s3   = m32 + (s3min-m32)*exp(anom1*csi2(1))
      else
         s3   = (m42-s) + (s3min-(m42-s))*exp(anom2*csi2(1))
      endif

      w = p1/anom1/(s3-m32)
      w = w + p2/anom2/(s3-(m42-s))

      w = 1.d0/w
      return
***   as p1/(s3-m32)  + p2/sqrt(s3max-s3)
 112  anom1 = log((s3max-m32)/(s3min-m32))
      anom2 = 2.d0 * sqrt(s3max-s3min)

      p1 = 0.5d0
      p2 = 1.d0 - p1

      call getrnd(csi2,2)
      if (csi2(2).lt.p1) then
         s3   = m32 + (s3min-m32)*exp(anom1*csi2(1))
      else
c         s3   = (m42-s) + (s3min-(m42-s))*exp(anom2*csi2(1))

         s3 = s3max - (anom2*0.5d0*csi2(1))**2
         
      endif

      w = p1/anom1/(s3-m32)
      w = w + p2/anom2/sqrt(s3max -s3)

      w = 1.d0/w
      return
***   as p1/(s3-m32)  + p2 * sqrt(s3max-s3)
 113  anom1 = log((s3max-m32)/(s3min-m32))
      anom2 = 1.d0/1.5d0 * (s3max-s3min) ** 1.5d0

      p1 = 0.5d0
      p2 = 1.d0 - p1

      call getrnd(csi2,2)
      if (csi2(2).lt.p1) then
         s3   = m32 + (s3min-m32)*exp(anom1*csi2(1))
      else

       s3 = s3max - (1.5d0 * anom2*csi2(1))**(2.d0/3.d0)
         
      endif

      w = p1/anom1/(s3-m32)
      w = w + p2/anom2*sqrt(s3max -s3)

      w = 1.d0/w
      
      return

***   1/(s3-m32)^2 no, proprio no
 114  anom = 1.d0/(s3min-m32) - 1.d0/(s3max-m32)

      call getrnd(csi1,1)
      
      x = 1.d0/(s3min-m32) - anom*csi1(1)

      s3 = m32 + 1.d0/x

      w = anom * (s3-m32)*(s3-m32)

      
      return

***   1/s3 no, proprio no
 115  anom = log(s3max/s3min)

      call getrnd(csi1,1)

      s3 = s3min * exp(anom*csi1(1))

      w = anom * s3
      
      return

***   as 1/(s3-m32)**2
 116  call getrnd(csi1,1)
      anom = 1.d0/(s3min-m32) - 1.d0/(s3max-m32)

c      anom * csi1(1)= 1.d0/(s3min-m32) - 1.d0/(s3-m32)

c      1.d0/(s3min-m32) - anom * csi1(1)=  1.d0/(s3-m32)      


      ww = 1.d0/(s3min-m32) - anom * csi1(1)

c      ww =  1.d0/(s3-m32)      
      
      s3 = 1.d0/ww + m32

c      print*,s3
      w    = anom * (s3-m32)**2
      return

      
      end
************************************************************
      subroutine get_t(t,tmin,tmax,s,m12,m22,w)
      implicit double precision (a-h,m,o-z)
      double precision csi(1),csi2(1)
      common/nfotonicommon/nfotoni

      ant  = 1.d0/tmin - 1.d0/tmax
c      call getrnd(csi,1)
c      t    = tmin / (1.d0 - ant * csi(1) * tmin)
c      w    = ant*t*t
c      return

      A = 2.d0*s
      B = s*s + m12*m12-2.d0*m12*s+m22*m22-2.d0*m22*s+2.d0*m12*m22
      B = 2.d0 * B
      an1 = tmax - tmin
      an2 = log(tmax/tmin)*A
      an3 = (1.d0/tmin - 1.d0/tmax)*B
      an23 = an2 + an3

      fmax = 1.d0+A/tmax+B/tmax/tmax
      fmax = fmax*0.25d0

c      an4  = - fmax * tmax      ! flat
c      an4 = -fmax * tmax * 0.5d0 ! as fmax*t/tmax

      
      sqtmax = sqrt(-tmax)
      fsqtmax = fmax*sqtmax
      
      an4 = 2.d0*fsqtmax * sqtmax ! as fmax*sqrt(-tmax)/sqrt(-t)

      if (nfotoni.eq.0) an4 = 0.d0
      
      an1234 = an1 + an23 + an4
      p1 = an1/an1234
      p2 = p1 + an23/an1234
      p3 = p2 + an4/an1234
      
      call getrnd(csi,1)
      if (csi(1).lt.p1) then
         call getrnd(csi,1)
         t = (tmax - tmin)*csi(1)+tmin
      elseif(csi(1).ge.p1.and.csi(1).lt.p2) then
         istop = 0
         rmin = A/B*tmax + 1.d0
         do while(istop.eq.0)
            call getrnd(csi,1)
            t    = tmin / (1.d0 - ant * csi(1) * tmin)
            call getrnd(csi2,1)            
            if (csi2(1).lt.(A/B*t  + 1.d0)/rmin) istop = 1
         enddo
      else
         call getrnd(csi,1)
c         t = tmax * csi(1) ! flat

c         t = an4*csi(1)*2.d0*tmax/fmax+tmax*tmax
c         t = -sqrt(t) ! as fmax*t/tmax

         t = sqtmax-an4*csi(1)*0.5d0/fsqtmax
         t = -t*t

c         print*,t
         
      endif

      if (t.le.tmax) then
         w = (an1+an23+an4)/(1.d0+A/t+B/t/t)
      else
c         w = (an1+an23+an4)/fmax ! flat
c         w = (an1+an23+an4)/fmax/t*tmax ! as fmax*t/tmax
         w = (an1+an23+an4)/fsqtmax*sqrt(-t)
      endif
       
      return
      end

************************************************************
      subroutine get_Ee_new(ionlyan,an,E,Emin,Emax,A,s,mmu,m,c,w)
** v1
** better than v2
! that's ok, i try to improve it with the other one
      implicit double precision (a-h,m,o-z)
      double precision csi(1),csi2(1)
      common/nfotonicommon/nfotoni

      if (A.le.Emax.or.Emax.lt.Emin*(1.d0+1d-8)) then
!     flat if A < Emax...
         an = Emax -Emin
         if (ionlyan.eq.1) return
         call getrnd(csi,1)
         w = an
         E = an * csi(1) + Emin
         return
      endif
      
      ml = m
      if (Emin.le.1.001d0*m) then
         ml = m*0.1d0
      endif
c      ml = 2.d0*m !!!!!!!!
      
      an1 = log((A-Emin)/(A-Emax))

      x = 1.5d0

c      x = 1.6d0
      
      umx = 1.d0 - x
      
      an2 = 1.d0/umx*((Emax-ml)**umx - (Emin-ml)**umx)

      an = an1 + an2

      if (ionlyan.eq.1) return
      
c     p1 = 0.4d0
c      if (nfotoni.ge.1) then
c         p1 = 1.d0
c      else
         p1 = an1/an
c      endif
      p2 = 1.d0 - p1

      call getrnd(csi,1)
      if (csi(1).le.p1) then      
         call getrnd(csi,1)
         E = A - (A-Emin)*exp(-an1*csi(1))
      else
         call getrnd(csi,1)
         E = ml + (umx*an2*csi(1)+(Emin-ml)**umx)**(1.d0/umx)         
      endif

      if (E.eq.A) then
         w = 0.d0
         return
      endif

      w = 1.d0 /(p1/(A-E)/an1 + p2/(E-ml)**x/an2)

      return
      end      
************************************************************
      subroutine get_Ee_v2(E,Emin,Emax,A,s,mmu,m,c,w)
** v2
** v1 is better      
      implicit double precision (a-h,m,o-z)
      double precision csi(1),csi2(1)
c         call get_Ee_v1(E,Emin,Emax,A,s,mmu,m,c,w)
c         return

      ml = m
      
      an1 = log((A-Emin)/(A-Emax))

      dm2 = ml*ml-mmu*mmu
      sqs = sqrt(s)
      omax = s+dm2-2.d0*Emin*sqs
      omax = omax/(sqs - Emin + sqrt(Emin*Emin-ml*ml)*c)
      omax = 0.5d0*omax
      al = 0.5d0*(sqs+dm2/sqs)-omax
      rp = (1.d0-c)/sqs
      be = 1.d0-omax*rp
      aob = al/be

      if (Emin.lt.aob) aob = 0.9999d0*Emin

      if (Emin.lt.aob) then
         ! should not happen
         call get_Ee(E,Emin,Emax,A,s,mmu,m,c,w)
         return
      endif
      
      an2=(1.d0 - aob*rp)*log((Emax-aob)/(Emin-aob)) - rp*(Emax - Emin)

      
c      print*,an2,aob,rp,Emin,Emin-aob,c

c      p1 = 0.4d0
      p1 = an1/(an1+an2)
c      p1 = 0.2d0
      p2 = 1.d0 - p1

      call getrnd(csi,1)
      if (csi(1).le.p1) then      
         call getrnd(csi,1)
         E = A - (A-Emin)*exp(-an1*csi(1))
      else

         antmp = log((Emax-aob)/(Emin-aob))

         istop = 0
         do while(istop.eq.0)
            call getrnd(csi,1)
            E = Emin + (Emin-aob)*exp(antmp*csi(1))

            f     = (1.d0-aob*rp)/(E-aob)
            ftrue = (1.d0-aob*rp)/(E - aob) - rp

            call getrnd(csi2,1)
            if (csi2(1).lt.ftrue/f) istop = 1
            
         enddo
      endif      

      w = 1.d0 /(p1/(A-E)/an1 + p2/an2 * ((1.d0-aob*rp)/(E-aob) - rp))      

      return
      end      
************************************************************
      subroutine get_Ee(E,Emin,Emax,A,s,mmu,m,c,w)
** v1
** better than v2
! that's ok, i try to improve it with the other one
      implicit double precision (a-h,m,o-z)
      double precision csi(1),csi2(1)

      ml = m
      
      an1 = log((A-Emin)/(A-Emax))
c      an2 = log((Emax-ml)/(Emin-ml))

      x   = 1.5d0
      umx = 1.d0 - x 
      an2 = 1.d0/umx*((Emax-ml)**umx - (Emin-ml)**umx)

c      p1 = 0.4d0
      p1 = an1/(an1+an2)
      p2 = 1.d0 - p1

      call getrnd(csi,1)
      if (csi(1).le.p1) then      
         call getrnd(csi,1)
         E = A - (A-Emin)*exp(-an1*csi(1))
      else
         call getrnd(csi,1)
c         E = ml + (Emin-ml) * exp(an2*csi(1))      

         E = ml + (umx*an2*csi(1)+(Emin-ml)**umx)**(1.d0/umx)
         
      endif      

c      w = 1.d0 /(p1/(A-E)/an1 + p2/(E-ml)/an2)

      w = 1.d0 /(p1/(A-E)/an1 + p2/(E-ml)**x/an2)
      

      return
      end      
************************************************************
      subroutine get_Ee_v3(E,Emin,Emax,A,s,mmu,m,c,w)
** v3
      implicit double precision (a-h,m,o-z)
      double precision csi(1),csi2(1)

      ml = m*m/sqrt(s)
      
      an1 = log((A-Emin)/(A-Emax))
c      an2 = log((Emax-ml)/(Emin-ml))

      x   = 1.5d0
      umx = 1.d0 - x 
      an2 = 1.d0/umx*((Emax-ml)**umx - (Emin-ml)**umx)

c      p1 = 0.4d0
      p1 = an1/(an1+an2)
c      p1 = 0.5d0
      p2 = 1.d0 - p1

      call getrnd(csi,1)
      if (csi(1).le.p1) then      
         call getrnd(csi,1)
         E = A - (A-Emin)*exp(-an1*csi(1))
      else
         call getrnd(csi,1)
c         E = ml + (Emin-ml) * exp(an2*csi(1))      

         E = ml + (umx*an2*csi(1)+(Emin-ml)**umx)**(1.d0/umx)
         
      endif      

c      w = 1.d0 /(p1/(A-E)/an1 + p2/(E-ml)/an2)

      w = 1.d0 /(p1/(A-E)/an1 + p2/(E-ml)**x/an2)
      

      return
      end      
************************************************************
      subroutine get_cth(c,cmin,cmax,s,pm,m12,m22,w)
      implicit double precision (a-h,m,o-z)
      double precision csi(1),csi2(1)
      double precision lambda
      external lambda
      common/nfotonicommon/nfotoni

      an = 1.d0/(1.d0+cmin) - 1.d0/(1.d0+cmax) ! needed below!

c      if (nfotoni.ge.1) then
c      call getrnd(csi,1)
c      upcmin = 1.d0+cmin
c      c = upcmin/(1.d0-upcmin*an*csi(1)) - 1.d0      
c      w = an * (1.d0+c)**2
c      return
c      endif

      A  = 2.d0*s / (-2.d0*pm*pm)
      CF = 1.d0 ! coefficient for FLAT, 0 or 1

      if (nfotoni.ge.1) then
         A  =  A*0.d0 ! A = 0, CF > 0 seems to be the best choice for integrated xs
         CF = 50.d0 ! it seems that raising this is better
      endif
      
!      B = s*s + m12*m12 +m22*m22 -2.d0*m12*s-2.d0*m22*s + 2.d0*m12*m22
!or      
      B = (s-m12-m22)**2

!      B   = 2.d0 * B / 4.d0 /pm**4
      B   = 0.5d0 * B / pm**4
      an1 = CF * (cmax - cmin)
      an2 = A * log((1.d0+cmax)/(1.d0+cmin))
      an3 = B * ( 1.d0/(1.d0+cmin) - 1.d0/(1.d0+cmax) )

      an23 = an2 + an3

      an4 = 0
      CL  = 0.d0
      if (cmin.lt.-0.9999d0) then
         CL = B/(1.d0+cmin)/(1.d0+cmin)*6d-6
         an4 = CL*2.d0
      endif
      
      an1234 = an1 + an23 + an4
      p1  = an1/an1234
      p23 = p1 + an23/an1234
      p4 = p1 + p23 + an4/an1234
      
      call getrnd(csi,1)
      if (csi(1).lt.p1) then
         call getrnd(csi,1)
         c = (cmax - cmin)*csi(1)+cmin
      elseif (csi(1).ge.p1.and.(csi(1)).lt.p23) then
         istop = 0
         rmin = A/B*(1.d0+cmin) + 1.d0
         do while(istop.eq.0)
            call getrnd(csi,1)
            upcmin = 1.d0+cmin
            c = upcmin/(1.d0-upcmin*an*csi(1)) - 1.d0
            if (abs(A).gt.0) then
               call getrnd(csi2,1)
               if (csi2(1).lt.(A/B*(1.d0+c)  + 1.d0)/rmin) istop = 1
            else
               istop = 1
            endif
         enddo
      else ! if p23 < csi < p4
         call getrnd(csi,1)
         c = 2.d0*csi(1) - 1.d0
      endif

      if (c.ge.cmin) then
         w = an1234/(CF+A/(1.d0+c)+B/(1.d0+c)/(1.d0+c) + CL)
      else
         w = an1234/CL
      endif
         
      return
      end
************************************************************
      subroutine get_cth_larger(c,cmin,cmax,s,pm,m12,m22,il,w)
      implicit double precision (a-h,m,o-z)
      double precision csi(1),csi2(1)
      common/nfotonicommon/nfotoni

      an = 1.d0/(1.d0+cmin) - 1.d0/(1.d0+cmax) ! needed below!

c      if (nfotoni.ge.1) then
c      call getrnd(csi,1)
c      upcmin = 1.d0+cmin
c      c = upcmin/(1.d0-upcmin*an*csi(1)) - 1.d0      
c      w = an * (1.d0+c)**2
c      return
c      endif
      
      A  = 2.d0*s / (-2.d0*pm*pm)
      CF = 1.d0 ! coefficient for FLAT, 0 or 1
      
      if (nfotoni.ge.1) then
         A  =  A*0.075d0 ! A = 0, CF > 0 seems to be the best choice for integrated xs
         CF =  250.d0    ! it seems that raising this is better
c         A = 0.d0
c         CF = 50.d0
      endif
      
      B   = s*s + m12*m12-2.d0*m12*s+m22*m22-2.d0*m22*s+2.d0*m12*m22
      B   = 2.d0 * B / 4.d0 /pm**4
      an1 = CF * (cmax - cmin)
      an2 = A * log((1.d0+cmax)/(1.d0+cmin))
      an3 = B * an
      
      an23 = an2 + an3

      an4 = 0.d0
      CL  = 0.d0
      if (il.gt.0) then
         CL = B/(1.d0+cmin)/(1.d0+cmin) *6d-6         
         an4 = CL*2.d0
      endif
      
      an1234 = an1 + an23 + an4
      p1     = an1/an1234
      p23    = p1 + an23/an1234
!      p4     = p1 + p23 + an4/an1234 !!! WRONG !!!
      p4     = p23 + an4/an1234

      call getrnd(csi,1)
      if (csi(1).lt.p1) then
         call getrnd(csi,1)
         c = (cmax - cmin)*csi(1)+cmin
      elseif (csi(1).gt.p1.and.csi(1).lt.p23) then
         istop = 0
         rmin = A/B*(1.d0+cmin) + 1.d0
         do while(istop.eq.0)
            call getrnd(csi,1)
            upcmin = 1.d0+cmin
            c = upcmin/(1.d0-upcmin*an*csi(1)) - 1.d0
            if (abs(A).gt.0) then
               call getrnd(csi2,1)
               if (csi2(1).lt.(A/B*(1.d0+c)  + 1.d0)/rmin) istop = 1
            else
               istop = 1
            endif
         enddo
      else ! if p23 < csi < p4
         call getrnd(csi,1)
         c = 2.d0*csi(1) - 1.d0
      endif

      if (c.ge.cmin) then
         w = an1234/(CF+A/(1.d0+c)+B/(1.d0+c)/(1.d0+c) + CL)
      else
         w = an1234/CL
      endif
      return
      end
************************************************************
      subroutine get_cth_pairs(c,cmin,cmax,s,pm,m12,m22,w)
      implicit double precision (a-h,m,o-z)
      double precision csi(1),csi2(1)
      common/nfotonicommon/nfotoni

      cminthr = -0.999d0
      if (cmin.lt.cminthr) cmin = cminthr
      
      an = 1.d0/(1.d0+cmin) - 1.d0/(1.d0+cmax) ! needed below!

c      if (nfotoni.ge.1) then
c      call getrnd(csi,1)
c      upcmin = 1.d0+cmin
c      c = upcmin/(1.d0-upcmin*an*csi(1)) - 1.d0      
c      w = an * (1.d0+c)**2
c      return
c      endif

      A  = 2.d0*s / (-2.d0*pm*pm)
      B  = s*s + m12*m12-2.d0*m12*s+m22*m22-2.d0*m22*s+2.d0*m12*m22
      B  = 2.d0 * B / 4.d0 /pm**4
      CF = 1.d0 ! coefficient for FLAT, 0 or 1
      
      if (nfotoni.ge.1) then
         A  =  A*0.075d0 ! A = 0, CF > 0 seems to be the best choice for integrated xs
         CF =  250.d0    ! it seems that raising this is better
c         A = 0.d0
c     CF = 50.d0

         A  = 0.d0
         B  = 1.d0
         CF = 0.d0
         
      endif
      
      an1 = CF * (cmax - cmin)
      an2 = A * log((1.d0+cmax)/(1.d0+cmin))
      an3 = B * an
      
      an23 = an2 + an3

      CL = B/(1.d0+cmin)/(1.d0+cmin) * 1.d0
      an4 = CL*2.d0
      
      an1234 = an1 + an23 + an4
      p1     = an1/an1234
      p23    = p1 + an23/an1234
!      p4     = p1 + p23 + an4/an1234 !!! WRONG !!!
      p4     = p23 + an4/an1234

c      print*,p1,p23,p4
      
      call getrnd(csi,1)
      if (csi(1).lt.p1) then
         call getrnd(csi,1)
         c = (cmax - cmin)*csi(1)+cmin
      elseif (csi(1).gt.p1.and.csi(1).lt.p23) then
         istop = 0
         rmin = A/B*(1.d0+cmin) + 1.d0
         do while(istop.eq.0)
            call getrnd(csi,1)
            upcmin = 1.d0+cmin
            c = upcmin/(1.d0-upcmin*an*csi(1)) - 1.d0
            if (abs(A).gt.0) then
               call getrnd(csi2,1)
               if (csi2(1).lt.(A/B*(1.d0+c)  + 1.d0)/rmin) istop = 1
            else
               istop = 1
            endif
         enddo
      else ! if p23 < csi < p4
         call getrnd(csi,1)
         c = 2.d0*csi(1) - 1.d0
      endif

      if (c.ge.cmin) then
         w = an1234/(CF+A/(1.d0+c)+B/(1.d0+c)/(1.d0+c) + CL)
      else
         w = an1234/CL
      endif
      return
      end
***************************************************************************            
      function gau_spread(e,s)
      implicit double precision (a-h,o-z)
      double precision csi(2)
      common/gauspreadfun/tpi,x1,x2,ifirst,icalc
      data ifirst /0/
      if (ifirst.eq.0) then
         tpi    = 8.d0 * atan(1.d0)
         icalc  = 1
         ifirst = 1
      endif
***   using Box-Mueller algorithm to generate 2 normal numbers
      x = 0.d0
      if (icalc.eq.1) then
        call getrnd(csi,2)
c         call anotherrng(csi,2)         
         r  = sqrt(-2.d0*log(csi(1))) 
         x1 = r*sin(tpi*csi(2)) 
         x2 = r*cos(tpi*csi(2)) 
         x  = x1
c     * since it's used also for "detector smearing", where s is not always
c     * constant, I always recalculate it at price of throwing one away...
cc NOT REALLY: here x1 and x2 are distributed as Gaussian with mean 0 and sigma = 1.
cc the translation to true sigma and mean is done in the last line...
         icalc = 0
cc         ical = 1
      else
         x     = x2
         icalc = 1
      endif
      gau_spread = e + s*x
      return
      end
*************************************************************************
      subroutine get_pattern(ncharged,n,ep)
      integer n,ep(n),ncharged
      double precision csi(1)
      do k = 1,n
         call getrnd(csi,1)
         ep(k) = 1.d0*ncharged*csi(1) + 1
         if (ep(k).gt.ncharged) ep(k) = ncharged ! to avoid round-offs...
      enddo
      return
      end 
*************************************************************************
      subroutine get_patternNEW(ncharged,n,ep,w)
      implicit double precision (a-h,o-z)
      integer n,ep(n),ncharged
      double precision csi(1)
      double precision pc(0:4)
      common/emispattern/pis,pfs,pc,ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         pis = 0.5d0 ! 0.96d0
         pfs = 1.d0 - pis
c         do k = 1,5
c            print*,pis*100,' % photons from IS!!'
c         enddo
         pc(0) = 0.d0
         pc(1) = pc(0) + pis*0.5d0
         pc(2) = pc(1) + pis*0.5d0
         pc(3) = pc(2) + pfs*0.5d0
         pc(4) = pc(3) + pfs*0.5d0
         ifirst = 1
      endif
      w    = 1.d0
      do k = 1,n
         call getrnd(csi,1)
         i = 1
         do while(csi(1).gt.pc(i))
            i = i + 1
         enddo
         ep(k) = i
         w = w * 1.d0/(pc(i)-pc(i-1))
      enddo
      return
      end 
*-----------------------------------------------------
      subroutine scanvpols
      implicit double precision (a-h,o-z)
      parameter (n=100000)
      double complex vpolc,result
      external vpolc
      common/resonances/ires,ionHAD
      common/jpsiparameters/amjpsi,gtjpsi,gejpsi,effcjpsi
      common/jpsiparameters2/am2s,gt2s,ge2s,coup2s ! from vpolc
      
      Emin = 0.1d0
      Emax = 5.5d0

      Emin = (amjpsi - 50.d0*gtjpsi)*.5d0
      Emax = (amjpsi + 50.d0*gtjpsi)*.5d0

      emin = 1.5d0
      emax = 1.55d0
      
      d = Emax - Emin
      d = d / n
      E = Emin
c      if (ires.eq.1) open(33,file='vpolnsk',status='unknown')
c      if (ires.eq.0) open(33,file='vpolres0',status='unknown')
      open(33,file='arunres2',status='unknown')
      do k = 0,n
         s = 4.d0*E*E
         if (E.lt.0.d0) s = -s
         result = vpolc(s)
         write(33,*)2.d0*E,dReal(result),Imag(result)
         E = E + d
      enddo
      close(33)
      print*,'STOP IN SCANVPOLS'
      stop
      
      return
      end
***************************************************
      subroutine ap_vertex(x,omx)
      implicit double precision (a-h,o-z)
*  x generation according to ap splitting function
*  (1+x^2)/(1-x), 0 <= x <= 1-eps
      double precision eps
      common/epssoft/eps
      double precision r(2)

      integer ifirst
      data ifirst /0/
      save ifirst

      double precision alne
      save alne

      if (ifirst.eq.0) then
        ifirst = 1
        alne = dlog(eps)
      endif

      irigenera = 1
      do while(irigenera.eq.1)
         call getrnd(r,2)
         cx  = r(1)
         omx = dexp(cx*alne)
         x   = 1.d0 - omx
         rx  = r(2)*2.d0/omx
         px  = (1.d0 + x**2)/omx
         if (rx.lt.px) irigenera = 0
      enddo
      return
      end
***
      subroutine ap_vertexNEW(eps,x,omx)
      implicit double precision (a-h,o-z)
*  x generation according to ap splitting function
*  (1+x^2)/(1-x), 0 <= x <= 1-eps
      double precision eps
      double precision r(2)
      alne = log(eps)
      irigenera = 1
      do while(irigenera.eq.1)
         call getrnd(r,2)
         cx  = r(1)
         omx = dexp(cx*alne)
         x   = 1.d0 - omx
         rx  = r(2)*2.d0/omx
         px  = (1.d0 + x**2)/omx
         if (rx.lt.px) irigenera = 0
      enddo
      return
      end
************************************************************
      subroutine multiplicity(eps,ecms,cth,n,np,w)
      implicit double precision (a-h,l,o-z)
      double precision csi(1)
      dimension vect(0:8)
      character*6 ord
      common/qedORDER/ord
      common/nphot_mode/nphotmode
      common/parameters/ame,ammu,convfac,alpha,pi
      common/photoncutoff/egmin,egmin2
      double precision phmass
      common/photonmasslambda/phmass

*** cuts and inputs            
      common/mueexpsetup/emulab,eemin,eemax,semu,thmumin,themin,themax,
     .     thmumax,ththr,Ethr,dthna7max,cutela,ina7,iela
********
      
**** charges to switch on radiation on legs (0 or 1), from invariants.h
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
********************************************************************
      common/molteplicita/arg,a,b,c,ppairs,ialsopairs,ifirst
      data ifirst /0/
      if (ifirst.eq.0) then
         arg = alpha/pi*(-log(2.d0*eps))*log(4.d0*ecms*ecms/ame/ame)

         if (QRmu.eq.1.and.QRe.eq.0) then
            arg=alpha/pi*(-log(2.d0*eps))*log(4.d0*ecms*ecms/ammu/ammu)
         endif
         
         a = abs(1.d0 - arg + 0.5d0*arg*arg)
         b = abs(arg - 0.5d0*arg*arg)
         c = 0.5d0*arg*arg

         if (ord.eq.'alpha') a = abs(1.d0 - arg)
         
         c = c * 2.d0 ! this seems to reduce the final error
!     I manipulate the probabilities
c         a = 1.5d0 * a
c         c = 2.d0  * c

*********************
         ialsopairs = 0
         ppairs = 0.d0
         if (nphotmode.gt. 999) ialsopairs = 1
         if (nphotmode.lt.-999) ialsopairs = 1
         if (ialsopairs.eq.1) then
            ppairs = 0.05d0
            ppairs = 0.1d0
            if (thmumin.gt.1d-5) then
               ppairs = 0.001d0
            endif
         endif
*******************
         
         ifirst =1
      endif
***
      np = 0  ! regulates pair production: np=0 no pairs, =1 electron pairs, =2 muon pairs, =10 pi0
*** dealing with nphotmode >= 0 cases
      if (nphotmode.gt.999) then
         w  = 1.d0
         n  = 0
         if (nphotmode.eq.1001) then
            np = 1              ! electron pairs
            return
         endif
         if (nphotmode.eq.1002) then
            np = 2              ! muon pairs
            return
         endif

!! pion_0         
         if (nphotmode.eq.1010) then
            np = 10 !pion_0
            return
         endif
****
         
! here if nphomode = 1000         
         pelp = 0.999d0
         pmup = 1.d0 - pelp
         call getrnd(csi,1)
         if (csi(1).lt.pelp) then
            np = 1
            w  = 1.d0/pelp
         else
            np = 2
            w  = 1.d0/pmup
         endif
         return
      endif
**      
      if (nphotmode.ge.0) then
         w = 1.d0
         n = nphotmode
         return
      endif

****** up-to-here all nphotmode >= 0 cases are dealt with

***** and now dealing with nphomode < 0
*****************************************************************
      if (ord.eq.'born') then
         if (ialsopairs.eq.0) then
            n = 0
            w = 1.d0
         else
            n      = 0
            p0     = 1.d0 - ppairs
            call getrnd(csi,1)
            if (csi(1).lt.p0) then
               np = 0
               w = 1.d0/p0
            else
               w = 1.d0/ppairs
               pelp = 0.999d0
               pmup = 1.d0 - pelp               
               call getrnd(csi,1)
               if (csi(1).lt.pelp) then
                  np = 1
                  w = w / pelp
               else
                  np = 2
                  w = w / pmup
               endif
            endif
         endif
         return
      endif
*****************************************************************
      if (ord.eq.'alpha') then
         
         p0 = a * 0.5d0 ! better efficiency         
         p1 = 1.d0 - p0

         p0 = p0 * (1.d0-ppairs)
         p1 = p1 * (1.d0-ppairs)

         np = 0
            
         call getrnd(csi,1)
         if (csi(1).lt.p0) then
            n = 0
            w = 1.d0 / p0
         elseif (csi(1).lt.(p0+p1)) then
            n = 1
            w = 1.d0 / p1
         else
            n = 0
            w = 1.d0/ppairs
            
            pelp = 0.999d0
            pmup = 1.d0 - pelp
            call getrnd(csi,1)
            if (csi(1).lt.pelp) then
               np = 1
               w  = w /pelp
            else
               np = 2
               w  = w /pmup
            endif
         endif
         return
      endif
*****************************************************************
      if (ord.eq.'alpha2') then
         p0 = a/(a+b+c)
         p1 = b/(a+b+c)
         p2 = c/(a+b+c)


         p0 = p0 * (1.d0-ppairs)
         p1 = p1 * (1.d0-ppairs)
         p2 = p2 * (1.d0-ppairs)
         
         np = 0

         call getrnd(csi,1)
         if (csi(1).lt.p0) then 
            n = 0
            w = 1.d0 / p0
         elseif (csi(1).gt.p0.and.csi(1).lt.(p0+p1)) then 
            n = 1
            w = 1.d0 / p1
         elseif(csi(1).gt.(p0+p1).and.csi(1).lt.(p0+p1+p2)) then 
            n = 2
            w = 1.d0 / p2
         else
            n = 0
            w = 1.d0/ppairs

            pelp = 0.999d0
            pmup = 1.d0 - pelp
            call getrnd(csi,1)
            if (csi(1).lt.pelp) then
               np = 1
               w  = w /pelp
            else
               np = 2
               w  = w /pmup
            endif
         endif
         return
      endif
      end
**********************************************************************************************
      subroutine getpeinverysimple(pp,p,w)
      implicit double precision (a-h,m,o-z)
      double precision p(0:3),csi(1),xi(1),pp(0:3)

      double precision Eioniz(4)
      
      common/parameters/ame,ammu,convfac,alpha,pi
      common/ifgetpeinverysimple/a0,gamma,zm,zmax,anorm1,anorm2,
     .     anorm,anormtrue,Eioniz,ifirst

      data ifirst /0/

      if (ifirst.eq.0) then

         Z = 4.
         n = 1         
         a0 = 1.d0/ame/alpha ! Bohr radius
         gamma = Z/n/a0

         
         Eioniz(1) = 1.0364d-2 * 899.5d0
         Eioniz(2) = 1.0364d-2 * 1757.1d0
         Eioniz(3) = 1.0364d-2 * 14848.7d0
         Eioniz(4) = 1.0364d-2 * 21006.6d0
         Eioniz = Eioniz * 1d-9 ! from eV to GeV         

c         print*,gamma/Z * 4
c         print*,gamma/Z * 3
c         print*,gamma/Z * 2* 0.5d0
c         print*,gamma/Z * 1* 0.5d0
c         print*,sqrt(2.d0*ame*Eioniz)
c         stop

         
         do k = 1,10
            print*,'ATOMIC EFFECTS INCLUDED!'
         enddo
         
         ifirst = 1
      endif

      w = 1.d0
      
      call getrnd(csi,1)

      work = Eioniz(int(4.d0*csi(1)+1))
      
      pp(0) = pp(0) - work
      pp(3) = sqrt(pp(0)**2-ammu*ammu)
      
      return
      end
********************************************************
      subroutine getpein(pp,p,w)
      implicit double precision (a-h,m,o-z)
      double precision p(0:3),csi(1),xi(1),pp(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi

      double precision Eioniz(4),gammas(4)
      common/ifgetpein/a0,gamma,zm,zmax,anorm1,anorm2,anorm,anormtrue,
     .     Eioniz,gammas,ifirst

      data ifirst /0/
************************************************************
***   Podolsky Pauling 1929
***   https://journals.aps.org/pr/pdf/10.1103/PhysRev.34.109      
************************************************************
      truefun(xx) = xx*xx/(1.d0+xx*xx)**4

      if (ifirst.eq.0) then

         Z = 4.d0
         n = 1
         
         a0 = 1.d0/ame/alpha ! Bohr radius
         gamma = Z/n/a0

         zm = 2.83996566d0
         zmax = 1d20

         Eioniz(1) = 1.0364d-2 * 899.5d0
         Eioniz(2) = 1.0364d-2 * 1757.1d0
         Eioniz(3) = 1.0364d-2 * 14848.7d0
         Eioniz(4) = 1.0364d-2 * 21006.6d0
         Eioniz = Eioniz * 1d-9 ! from eV to GeV         

         gammas(4) = gamma/Z * Z * 1.0000489713477514d0
         gammas(3) = gamma/Z * (Z-1.d0) * 1.1210534966798504d0
         gammas(2) = gamma/Z * (Z-2.d0) * 0.5d0 * 1.1569154349729127d0
         gammas(1) = gamma/Z * (Z-3.d0) * 0.5d0 * 1.6555181841076323d0
c     print*,Eioniz/(gammas**2*0.5d0/ame)

c         print*,Eioniz*1000d0
c     print*,Eioniz*1d6
c         print*,'stop line 95 sampling.f'
c         stop
         

!     ! I do not have relativistic effects. I tolerate a maximum lorentz gamma
!     ! for the IS electron as given by gmax (which corresponds to a given zmax)
         gmax = 1.05d0
c         zmax  = sqrt(gmax*gmax - 1.d0)*ame/gamma
         
         anorm1 = 0.125d0*(1.d0 - 1.d0/(zm**2+1.d0)**2)
         anorm2 = 0.2d0*(1.d0/zm**5-1.d0/zmax**5)
         anorm = anorm1 + anorm2
         anormtrue =
     .        zmax*(3.d0*zmax**4+8.d0*zmax**2-3.d0)/(zmax**2+1.d0)**3
     .        +3.d0*atan(zmax)
         anormtrue = anormtrue - 3.d0*atan(0.d0)
         anormtrue = anormtrue/48.d0

c         print*,anormtrue
c         stop

         do k = 1,10
            print*,'ATOMIC EFFECTS INCLUDED!'
         enddo
         
         ifirst = 1
      endif

c      ic = 0
 111  continue
c      ic = ic + 1
      call getrnd(csi,1)
      if (csi(1)*anorm.lt.anorm1) then
         z = sqrt(1.d0/sqrt(1.d0-8.d0*anorm*csi(1)) -1.d0)
      else
         z = 5.d0*(anorm*csi(1)-anorm1) + 1.d0/zmax**5
         z = 1.d0/z**0.2d0
      endif
      
      call getrnd(xi,1)
      if (xi(1)*falsefun(z,zm).gt.truefun(z)) goto 111
      
c      w = anorm/anormtrue/falsefun(z,zm)*truefun(z)
      w = 1.d0

      call getrnd(csi,1)
      phi = 2.d0*pi*csi(1)
      call getrnd(csi,1)
      c = 2.d0*csi(1)-1.d0
      s = sqrt(1.d0-c*c)

      gammarun = gamma
      call getrnd(csi,1)
      gammarun = gammas(int(4.d0*csi(1)+1.d0))
      
      p(0) = sqrt(z*z*gammarun*gammarun + ame*ame)
! I leave it on-shell for the moment
c      p(0) = ame - gammarun*gammarun*0.5d0/ame
      
      p(1) = gammarun * z * s * sin(phi)
      p(2) = gammarun * z * s * cos(phi)
      p(3) = gammarun * z * c

      return
      end
********************************************************
      function falsefun(x,xm)
      implicit double precision (a-h,o-z)
      if (x.lt.xm) then      
         falsefun = 0.5d0*x/(1.d0+x*x)**3
      else
         falsefun = 1.d0/x**6
      endif
      return
      end
*********************************************************
      subroutine getpeindummy(pp,p,w)
      implicit double precision (a-h,m,o-z)
      double precision p(0:3),csi(1),xi(1),pp(0:3)
      common/parameters/ame,ammu,convfac,alpha,pi

      double precision Eioniz(4)
      common/ifgetpeindummy/a0,gamma,zm,zmax,anorm1,anorm2,anorm,
     .     anormtrue,Eioniz,ifirst

      data ifirst /0/
************************************************************
***   Podolsky Pauling 1929
***   https://journals.aps.org/pr/pdf/10.1103/PhysRev.34.109      
************************************************************
ccc      truefun(xx) = xx*xx/(1.d0+xx*xx)**4

      if (ifirst.eq.0) then

         Z = 4.
         n = 1
         
         a0 = 1.d0/ame/alpha ! Bohr radius
         gamma = Z/n/a0

         zm = 2.83996566d0
         zmax = 1d20

         Eioniz(1) = 1.0364d-2 * 899.5d0
         Eioniz(2) = 1.0364d-2 * 1757.1d0
         Eioniz(3) = 1.0364d-2 * 14848.7d0
         Eioniz(4) = 1.0364d-2 * 21006.6d0
         Eioniz = Eioniz * 1d-9 ! from eV to GeV         
         

         anorm1 = 0.125d0*(1.d0 - 1.d0/(zm**2+1.d0)**2)
         anorm2 = 0.2d0*(1.d0/zm**5-1.d0/zmax**5)
         anorm = anorm1 + anorm2
         anormtrue =
     .        zmax*(3.d0*zmax**4+8.d0*zmax**2-3.d0)/(zmax**2+1.d0)**3
     .        +3.d0*atan(zmax)
         anormtrue = anormtrue - 3.d0*atan(0.d0)
         anormtrue = anormtrue/48.d0

c         print*,anormtrue
c         stop

         do k = 1,10
            print*,'ATOMIC EFFECTS INCLUDED!'
         enddo
         
         ifirst = 1
      endif

      w = 1.d0

****  I get z at max only
** this is the maximum of the function      
      z = 1.d0/sqrt(3.d0)
*     * this should be the average
c      z = 1.d0
      

      call getrnd(csi,1)
      work = Eioniz(int(4.d0*csi(1)+1))
      z = sqrt(2.d0*ame*Eioniz(int(4.d0*csi(1)+1)))
      z = z/gamma
      
      
      call getrnd(csi,1)
      phi = 2.d0*pi*csi(1)
      call getrnd(csi,1)
      c = 2.d0*csi(1)-1.d0
      s = sqrt(1.d0-c*c)

      p(0) = sqrt(z*z*gamma*gamma + ame*ame)
      p(0) = ame - work
      p(1) = gamma * z * s * sin(phi)
      p(2) = gamma * z * s * cos(phi)
      p(3) = gamma * z * c

      p(0) = ame - work      
      p(1) = sqrt(2.d0*ame*work) * s * sin(phi)
      p(2) = sqrt(2.d0*ame*work) * s * cos(phi)
      p(3) = sqrt(2.d0*ame*work) * c
      
      return
      end
****
      subroutine getc4(cmin,cmax,c,w)
      implicit double precision (a-h,m,o-z)
      double precision csi1(1)

      anmin = -1.d0/(1.d0+cmin)
      anorm = -1.d0/(1.d0+cmax) - anmin

      call getrnd(csi1,1)
      
      c = -1.d0 - 1.d0/(anmin+anorm*csi1(1))

      w = (1.d0 + c)*(1.d0 + c)*anorm
      
      return
      end
********************************************************
      
