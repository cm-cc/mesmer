      Z1_=p1
      Z2_=p3
      Z3_=p2
      Z4_=p4
      Z5_=k
      Z6_=cmu
      Z7_=ce
      Z8_=evkcec5mu
      Z9_=evkcmuc5e
      Z10_=dp3km1 - dp1km1
      Z11_=dp1km1 + dp3km1
      Z12_=dp4km1 - dp2km1
      Z13_=dp2km1 + dp4km1
      Z14_=dp1km1*Z1_
      Z15_= - Z2_*dp3km1
      Z14_=Z14_ + Z15_
      Z14_=2*Z14_
      Z15_=dp2km1*Z3_
      Z16_= - Z4_*dp4km1
      Z15_=Z15_ + Z16_
      Z15_=2*Z15_
      rad=t24m1*Z8_*Z10_ + t24m1*kcmu*Z7_*Z11_ - t24m1*kce*Z6_*Z11_ + 
     & t13m1*Z9_*Z12_ - t13m1*kcmu*Z7_*Z13_ + t13m1*kce*Z6_*Z13_ + 
     & cmuce*t24m1*Z14_ - cmuce*t24m1*Z5_*Z11_ + cmuce*t13m1*Z15_ - 
     & cmuce*t13m1*Z5_*Z13_


