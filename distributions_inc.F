      implicit double precision (a-h,o-z)
      parameter (ndistr  = 34)
#ifdef RELEASE
      parameter (nbin   = 100)	
* 2D
      parameter (ndistr2 = 1, i2d = 0) ! i2d = 1 calculate 2d distr, = 0 no
      parameter (nbin2 = 100)
#else
      parameter (nbin   = 100)
* 2D
      parameter (ndistr2 = 2, i2d = 0) ! i2d = 1 calculate 2d distr, = 0 no
      parameter (nbin2 = 20)
#endif	
      parameter(maxwgt=3)
      
      common/ecms/ecms,ecmsnom
      common/epssoft/eps
      dimension p3(0:3),p4(0:3),qph(40,0:3)
      integer*8 nevinbin(ndistr,nbin)
      common/nevxbin/nevinbin

      character*100 distname(ndistr,maxwgt),tmp(ndistr),
     .     distname2(ndistr2)
      character*100 tmp2(ndistr)
      integer lung(ndistr),lung2(ndistr2)
      dimension dis2d(ndistr2,nbin2,nbin2),edis2d(ndistr2,nbin2,nbin2)
      dimension bmi(ndistr),bma(ndistr)
      dimension bmi2(ndistr2,2),bma2(ndistr2,2),xaxis(ndistr,nbin)
      dimension distr(ndistr,nbin,maxwgt),fun(ndistr),fun2(ndistr2,2)
      dimension edistr(ndistr,nbin,maxwgt)
      dimension xaxis2(ndistr2,nbin2),yaxis2(ndistr2,nbin2)
      
      common/filenames/distname,distname2
      common/distr/distr,edistr,dis2d,edis2d,bmi,bma,bmi2,bma2,
     .     xaxis,xaxis2,yaxis2
