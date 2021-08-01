      double complex cmu(0:3),ce(0:3),c5mu(0:3),c5e(0:3)
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,tm1
      double precision p1k,p2k,p3k,p4k,s12,s34,t24,t13,u14,u23
      double precision p1mk2,p2mk2,p3pk2,p4pk2
      common/invariants/p1p2,p1p3,p1p4,p2p3,p2p4,p3p4,tm1
      common/invariants1g/p1k,p2k,p3k,p4k,s12,s34,t24,t13,u14,u23,
     .     p1mk2,p2mk2,p3pk2,p4pk2 

      double complex p1cmu,p2cmu,p3cmu,p4cmu,p1ce,p2ce,p3ce,p4ce
      double complex cImu,cIe,cmuce,c5muc5e
      double complex e13cmuc5e, e14cmuc5e, e14cec5mu, e34cec5mu
      double complex e12cmuc5e, e12cec5mu, e23cec5mu
      double complex emy12cmuc5e, emy12cec5mu, emy13cmuc5e, emy14cmuc5e,
     .     emy14cec5mu, emy23cec5mu, emy34cec5mu
      double complex p1c5mu,p2c5mu,p3c5mu,p4c5mu, 
     .     p1c5e, p2c5e, p3c5e, p4c5e

      common/variants/p1cmu,p2cmu,p3cmu,p4cmu,p1ce,p2ce,p3ce,p4ce,
     .     cImu,cIe,cmuce,c5muc5e,e13cmuc5e, e14cmuc5e, e14cec5mu,
     .     e34cec5mu,e12cmuc5e, e12cec5mu, e23cec5mu,
     .     p1c5mu,p2c5mu,p3c5mu,p4c5mu, 
     .     p1c5e, p2c5e, p3c5e, p4c5e,
     .     c5mu,c5e
      double complex cuno,im,czero
      parameter (cuno  = (1.d0,0.d0))
      parameter (im    = (0.d0,1.d0))
      parameter (czero = (0.d0,0.d0))
      double complex dotcmom,epsfuncmom,dot1cmom,epsfun2cmom
      double precision dot
      external dotcmom,dot,epscfuncmom,dot1cmom,epsfun2cmom

      double precision ec,ec2,ec4,mm,mm2,me,me2,mm4,me4,mm6,me6
      common/elchmasses/ec,ec2,ec4,mm,mm2,me,me2,mm4,me4,mm6,me6
      double precision pi216m1,pi216
      common/pi216m1/pi216m1,pi216


******* calculated in main      
      double precision p1lab(0:3),p2lab(0:3),p12lab(0:3)
      common/labmomenta/p1lab,p2lab,p12lab
******* calculated in main      
      
      integer h1,h2,h3,h4
      common/spinindex/h1,h2,h3,h4
      double precision Qmu
      integer ioneloop
      common/muoncharge/Qmu
      common/oneloopcalc/ioneloop

**** charges to switch on radiation on legs (0 or 1)      
      integer QRe,QRmu
      common/radiationcharges/QRe,QRmu
      
*********** for "cross" (Denner A.5)
      double complex Crl,CSpm,Crlmu,CSpmmu,Crle,CSpme
      double complex Crlmuoe2,CSpmmuoe2,Crleoe2,CSpmeoe2
      common/CrlCSpmoe2_mu_e/Crlmuoe2,CSpmmuoe2,Crleoe2,CSpmeoe2 ! these enter muemue1g1L
      double complex dZrlii,dmass

      double complex dZrliimu,dmassmu
      double complex dZrliie ,dmasse
c      common/dZsdm/dZrliimu,dmassmu,dZrliie,dmasse      
****************************      
      
***   for \Delta\alpha
      double precision dal,dah,vpt
      integer ifixedorder
      common/deltaalphas/dal,dah,vpt,ifixedorder
      double precision vpdal
      double complex vpdah
      double precision vpvect(5)
      common/vp_dalphas/vpdal,vpdah,vpvect
      double precision vp13,vp24
      common/vp1g1L/vp13,vp24

      integer ionlybubbles,ibubbleonv,ibubbleonboxsub
      common/ionlybubblescmn/ionlybubbles,ibubbleonv,ibubbleonboxsub
***
***   for reweighenting
      integer nextraweights
      parameter (nextraweights=4)
      double precision weightdefault,extraweights(0:nextraweights)
      double precision reweightLO,reweightNLO,wnorm
      common/weights/reweightLO,reweightNLO,weightdefault,extraweights
      common/weightnormalization/wnorm
***
