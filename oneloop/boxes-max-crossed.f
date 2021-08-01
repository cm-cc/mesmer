
      boxes =
     &  + D33ltc * ( 160*me2*p1p4**2 - 32*me2*p1p3*p2p4 + 64*me2*p1p2*
     &    p1p4 + 32*me2*p1p2**2 - 96*me2**2*p1p3 - 32*mm2*me2*p2p4 + 
     &    160*mm2*me2**2 )
      boxes = boxes + D23ltc * ( 256*p1p4**3 - 320*p1p2*p1p4**2 + 64*
     &    p1p2*p1p3*p2p4 + 128*p1p2**2*p1p4 - 64*p1p2**3 - 128*me2*p1p3
     &    *p1p4 + 64*me2*p1p2*p1p3 - 128*mm2*p1p4*p2p4 + 64*mm2*p1p2*
     &    p2p4 + 256*mm2*me2*p1p4 - 192*mm2*me2*p1p2 )
      boxes = boxes + D13ltc * ( 256*p1p4**3 + 128*p1p2**2*p1p4 - 128*
     &    me2*p1p3*p1p4 - 64*me2*p1p2*p1p3 - 128*mm2*p1p4*p2p4 - 64*mm2
     &    *p1p2*p2p4 + 384*mm2*me2*p1p4 + 128*mm2*me2*p1p2 )
      boxes = boxes + D22ltc * (  - 256*p1p3*p1p4**2 + 128*me2*p1p3**2
     &     + 256*mm2*p1p4**2 + 128*mm2*p1p3*p2p4 - 384*mm2*me2*p1p3 - 
     &    128*mm2**2*p2p4 + 256*mm2**2*me2 )
      boxes = boxes + D12ltc * (  - 256*p1p3*p1p4**2 + 128*me2*p1p3**2
     &     + 256*mm2*p1p4**2 + 128*mm2*p1p3*p2p4 - 384*mm2*me2*p1p3 - 
     &    128*mm2**2*p2p4 + 256*mm2**2*me2 )
      boxes = boxes + D11ltc * ( 160*mm2*p1p4**2 - 32*mm2*p1p3*p2p4 + 
     &    64*mm2*p1p2*p1p4 + 32*mm2*p1p2**2 - 32*mm2*me2*p1p3 - 96*
     &    mm2**2*p2p4 + 160*mm2**2*me2 )
      boxes = boxes + D00ltc * ( 512*p1p4**2 + 128*p1p2**2 - 320*me2*
     &    p1p3 - 320*mm2*p2p4 + 640*mm2*me2 )
      boxes = boxes + D3ltc * ( 128*p1p4**3 + 128*p1p2**2*p1p4 + 64*me2
     &    *p1p4**2 - 64*me2*p1p3*p1p4 - 64*me2*p1p2*p1p3 + 64*me2*
     &    p1p2**2 - 64*me2**2*p1p3 - 128*mm2*p1p4*p2p4 - 64*mm2*me2*
     &    p2p4 + 64*mm2*me2*p1p4 + 64*mm2*me2*p1p2 + 128*mm2*me2**2 )
      boxes = boxes + D2ltc * ( 128*p1p4**3 - 128*p1p3*p1p4**2 - 192*
     &    p1p2*p1p4**2 + 64*p1p2*p1p3*p2p4 + 128*p1p2**2*p1p4 - 64*
     &    p1p2**3 - 64*me2*p1p3*p1p4 + 64*me2*p1p3**2 - 64*mm2*p1p4*
     &    p2p4 + 128*mm2*p1p4**2 + 64*mm2*p1p3*p2p4 + 128*mm2*me2*p1p4
     &     - 192*mm2*me2*p1p3 - 64*mm2*me2*p1p2 - 64*mm2**2*p2p4 + 128*
     &    mm2**2*me2 )
      boxes = boxes + D1ltc * ( 128*p1p4**3 + 128*p1p2**2*p1p4 - 128*
     &    me2*p1p3*p1p4 - 64*mm2*p1p4*p2p4 + 64*mm2*p1p4**2 - 64*mm2*
     &    p1p2*p2p4 + 64*mm2*p1p2**2 + 64*mm2*me2*p1p4 - 64*mm2*me2*
     &    p1p3 + 64*mm2*me2*p1p2 - 64*mm2**2*p2p4 + 128*mm2**2*me2 )
      boxes = boxes + D0ltc * ( 128*p1p4**3 + 128*p1p2**2*p1p4 - 128*
     &    me2*p1p3*p1p4 - 128*mm2*p1p4*p2p4 + 256*mm2*me2*p1p4 )

