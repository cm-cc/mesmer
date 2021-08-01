RELEASE=yes
RELEASE=

EXE = mesmer

F77 = gfortran
CC  = gcc
# tune for your processor, if you want
ifdef RELEASE
  FFLAGS = -O3 -ffast-math # safer and more portable
else
  FFLAGS = -O3 -march=native -mtune=native -ffast-math
endif
#-freal-8-real-16   ! this promotes all real*8 to real*16 (quad precision)
### for debugging  FFLAGS = -g -fbounds-check -fbacktrace

SMH=
ifdef RELEASE
  DEFINERELEASE=-DRELEASE
else
  DEFINERELEASE=
  SMH=-DSMH
endif
SMH=

RECOLA=-DRECOLA
RECOLA=

COLLIER=
COLLIER=-DCOLLIER

## C ranlux optimizations, AVX2 faster but only on recent CPUS, SSE2 slower but present also on
## old hardware. If both are specified, -DSSE2 is ignored (and won't run on old hardware).
ifdef RELEASE
  RLXOPT=-DSSE2  # safer and more portable
else
  RLXOPT=-DSSE2 -DAVX2
endif

QUAD=-DRealType=real*16 -DComplexType=complex*32 -DQUAD
QUAD=

ROOTINTERFACE=
ROOTINTERFACE=yes

GSLLIBS= -lgsl -lgslcblas

F77 += $(FFLAGS)

default: $(EXE)

VPKNT=vp_knt_v3_0_1
##old one vp_hlmnt_v2_2
HPLOG=hplog

LTVER = 2.16
LTDIR = LoopTools-$(LTVER)
LTSTRING =-looptools
ifdef QUAD
  LTSTRING =-looptools-quad
endif

HANDYGDIR=handyG/
CHAPLINDIR=chaplin/
HANDYG=
CHAPLIN=
ifeq ($(SMH),-DSMH)
  HANDYG=-I$(HANDYGDIR)include -L$(HANDYGDIR)lib -lhandyg
##HANDYG=-IhandyGquad/include -LhandyGquad/lib -lhandyg
  CHAPLIN=-L$(CHAPLINDIR)lib -lchaplin
endif

ifeq ($(COLLIER),-DCOLLIER)
  CLLDIR = collier/COLLIER-1.2.5/
  CLLMOD = -I$(CLLDIR)include/
  CLLLIB = -L$(CLLDIR)lib -lcollier
endif

ifeq ($(RECOLA),-DRECOLA)
  RCLVERSION=1.4.0
  RCLDIR = recolas/recola-collier-$(RCLVERSION)
  RCLMOD = -I$(RCLDIR)/include/ -I$(RCLDIR)/recola-$(RCLVERSION)/include/
  RCLLIB = -L$(RCLDIR)/lib/ -lrecola -lcollier
endif

LIBFILES  = $(LTDIR)/lib64/libooptools.a
ifeq ($(ROOTINTERFACE),yes)
  LIBFILES += root-interface/write_MuE_MCevents.exe
endif
ifeq ($(COLLIER),-DCOLLIER)
  LIBFILES += $(CLLDIR)/lib/libcollier.a
endif
ifeq ($(SMH),-DSMH)
  LIBFILES += $(HANDYGDIR)/lib/libhandyg.a
  LIBFILES += $(CHAPLINDIR)/lib/libchaplin.a
endif

OBJECTS = main.o cuts.o sv.o matrix_model.o loops.o muemuegg.o realpairs.o vacuumpolarization.o\
          rngs.o routines.o sampling.o phasespacemue.o distributions.o $(VPKNT).o\
          hadr5n12.o hadr5n17.o hadr5x19.o hadr5n.o ranlux.o userinterface.o muemue1g1Lnoud.o\
          muemue1g1Lud.o storage.o c_rnlx_interface.o ranlux_common.o ranlxd.o ranlxs.o recola_int.o\
          twoloop_virtual.o $(HPLOG).o quadpack.o gsl_random.o elasticity.o
ifeq ($(SMH),-DSMH)
  OBJECTS += light-heavy-2LFF.o
endif

#### packaging
SAVEDIR    = release
RELEASEDIR = MESMER

STRINGSMH=
ifeq ($(SMH),-DSMH)
  STRINGSMH  = mkdir -p $(RELEASEDIR)/$(CHAPLINDIR)/ && cp -ra handyG-clean/ $(RELEASEDIR)/handyG &&
  STRINGSMH += cp -ra $(CHAPLINDIR)/chaplin-1.2-clean/ $(RELEASEDIR)/$(CHAPLINDIR)/chaplin-1.2 &&
endif

pack: # use only to release MESMER
	mkdir -p $(RELEASEDIR)/oneloop &&\
	mkdir -p $(RELEASEDIR)/c_ranlux &&\
	mkdir -p $(RELEASEDIR)/collier/ &&\
	cp -ra LoopTools-$(LTVER)-clean/ $(RELEASEDIR)/LoopTools-$(LTVER) &&\
	cp -ra collier/COLLIER-1.2.5-clean/ $(RELEASEDIR)/collier/COLLIER-1.2.5 &&\
        cp -ra root-interface/ write-root-events MuEtreeDict_rdict.pcm Makefile README.md input-example\
        distributions.F distributions_inc.F invariants.h muemue1g1Lnoud.F funsdeccmn1g1L.h\
        main.F matrix_model.F vpol_novosibirsk.dat vpol_novosibirsk_v2.dat muemue1g1Lud.F\
        vacuumpolarization.F cuts.F sv.F routines.f sampling.f phasespacemue.F $(VPKNT).f recola_int.F\
        rngs.F hadr5n12.f hadr5n17.f hadr5x19.f hadr5n.f ranlux.f userinterface.F\
        loops.F storage.F muemuegg-minus.f muemuegg-plus.f muemuegg.F muemue1g1Lnoupdown.f realpairs.F\
        realpairs_ampl2.f muemue1g1Lupdown.f funssettozero.f printltfun cts1g1L.f\
        twoloop_virtual.F dalhadslow17.f dalhadshigh17.f dalhadt17.f quadpack.F $(HPLOG).f\
        f1_light_heavy.f f2_light_heavy.f f1_light_heavy_quad.f f2_light_heavy_quad.f\
        light-heavy-2LFF.F constgpl_defs.f vargpl_defs.f gsl_random.c Rhad-scan.dat\
        elasticity.F\
        $(RELEASEDIR) &&\
	cp oneloop/*.f $(RELEASEDIR)/oneloop/ &&\
	cp c_ranlux/*.* $(RELEASEDIR)/c_ranlux/ && $(STRINGSMH) \
	rm -f $(RELEASEDIR)/root-interface/MuEtreeDict* $(RELEASEDIR)/root-interface/write_MuE_MCevents.exe
	mkdir $(SAVEDIR) ;\
	tar cjvf $(SAVEDIR)/$(RELEASEDIR).tar.bz2 $(RELEASEDIR)/ &&\
	rm -rf $(RELEASEDIR)
##### end packaging

clean:
	rm -f $(OBJECTS) *.a $(EXE) *~

# source files
main.o: main.F
	$(F77) -c main.F
cuts.o: cuts.F
	$(F77) -c cuts.F
$(HPLOG).o: $(HPLOG).f
	$(F77) -c $(HPLOG).f
matrix_model.o: matrix_model.F invariants.h funsdeccmn1g1L.h Makefile $(LIBFILES)
	$(F77) -c $(RECOLA) $(COLLIER) $(CLLMOD) matrix_model.F
sv.o: sv.F
	$(F77) -c sv.F
loops.o: loops.F Makefile oneloop/*.f invariants.h funsdeccmn1g1L.h $(LIBFILES)
	$(F77) -c $(QUAD) $(COLLIER) $(RECOLA) -I$(LTDIR)/include loops.F $(CLLMOD)
muemuegg.o: muemuegg.F invariants.h muemuegg-minus.f muemuegg-plus.f
	$(F77) -c muemuegg.F
muemue1g1Lud.o: muemue1g1Lud.F invariants.h funsdeccmn1g1L.h muemue1g1Lupdown.f
	$(F77) $(QUAD) -c muemue1g1Lud.F
muemue1g1Lnoud.o: muemue1g1Lnoud.F  invariants.h funsdeccmn1g1L.h muemue1g1Lnoupdown.f
	$(F77) $(QUAD) -c muemue1g1Lnoud.F
userinterface.o: userinterface.F Makefile
	$(F77) $(QUAD) $(DEFINERELEASE) -c userinterface.F
realpairs.o: realpairs.F realpairs_ampl2.f 
	$(F77) -c realpairs.F
phasespacemue.o: phasespacemue.F 
	$(F77) -c phasespacemue.F
vacuumpolarization.o: vacuumpolarization.F
	$(F77) -c vacuumpolarization.F
hadr5n12.o: hadr5n12.f 
	$(F77) -c hadr5n12.f
hadr5n17.o: hadr5n17.f 
	$(F77) -c hadr5n17.f
hadr5x19.o: hadr5x19.f 
	$(F77) -c hadr5x19.f
hadr5n.o: hadr5n.f 
	$(F77) -c hadr5n.f
$(VPKNT).o: $(VPKNT).f 
	$(F77) -c $(VPKNT).f
quadpack.o: quadpack.F 
	$(F77) -c quadpack.F
sampling.o: sampling.f
	$(F77) -c sampling.f
routines.o: routines.f 
	$(F77) $(QUAD) -c routines.f
distributions.o: distributions.F distributions_inc.F Makefile
	$(F77) -c distributions.F
ranlux.o: ranlux.f 
	$(F77) -c ranlux.f
rngs.o: rngs.F Makefile
	$(F77) $(QUAD) -c rngs.F
elasticity.o: elasticity.F
	$(F77) $(QUAD) -c elasticity.F
storage.o: storage.F  invariants.h
	$(F77) -c storage.F
recola_int.o: recola_int.F Makefile $(LIBFILES)
	$(F77) $(RECOLA) $(QUAD) -c $(RCLMOD) recola_int.F
twoloop_virtual.o: twoloop_virtual.F Makefile $(LIBFILES)
	$(F77) $(QUAD) $(COLLIER) -I$(LTDIR)/include $(SMH) -c twoloop_virtual.F $(CLLMOD)
light-heavy-2LFF.o: light-heavy-2LFF.F f1_light_heavy.f f2_light_heavy.f constgpl_defs.f vargpl_defs.f $(LIBFILES)
	$(F77) $(QUAD) -I$(LTDIR)/include $(HANDYG) $(SMH) -c light-heavy-2LFF.F -ffixed-line-length-85

# C version of ranlux by Martin Luscher, http://luscher.web.cern.ch/luscher/ranlux/index.html
c_rnlx_interface.o: c_ranlux/c_rnlx_interface.c  Makefile
	$(CC) $(FFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/c_rnlx_interface.c
ranlxd.o: c_ranlux/ranlxd.c  Makefile 
	$(CC) $(FFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlxd.c
ranlxs.o: c_ranlux/ranlxs.c  Makefile
	$(CC) $(FFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlxs.c
ranlux_common.o: c_ranlux/ranlux_common.c  Makefile
	$(CC) $(FFLAGS) $(RLXOPT) -std=c99 -Ic_ranlux/ -c c_ranlux/ranlux_common.c
##########
gsl_random.o: gsl_random.c
	$(CC) $(FFLAGS) -c gsl_random.c

### external libraries (LoopTools, Collier, HandyG, Chaplin, ROOT interface by Giovanni Abbiendi ####
extlibs: looptools collier handyg chaplin rootinterface

looptools: $(LTDIR)/lib64/libooptools.a
$(LTDIR)/lib64/libooptools.a:
	@echo " "
	@echo "Building LoopTools"
	@echo " "
	cd $(LTDIR) && ./configure --prefix=. && make && make install

collier: $(CLLDIR)/lib/libcollier.a
$(CLLDIR)/lib/libcollier.a:
	@echo " "
	@echo "Building Collier"
	@echo " "
	cd $(CLLDIR)/build/ && cmake .. -DCMAKE_INSTALL_PREFIX=.. -Dstatic=ON && make && make install

handyg: $(HANDYGDIR)/lib/libhandyg.a
$(HANDYGDIR)/lib/libhandyg.a:
	@echo " "
	@echo "Building handyG"
	@echo " "
	cd $(HANDYGDIR)/ && ./configure --prefix=. && make && make install

chaplin: $(CHAPLINDIR)/lib/libchaplin.a
$(CHAPLINDIR)/lib/libchaplin.a:
	@echo " "
	@echo "Building Chaplin"
	@echo " "
	cd $(CHAPLINDIR)/chaplin-1.2/ && ./configure --prefix=`pwd`/../ --disable-shared && make && make install

rootinterface: root-interface/write_MuE_MCevents.exe
root-interface/write_MuE_MCevents.exe:
	@echo " "
	@echo "Building ROOT interface (by G. Abbiendi)"
	@echo " "
	cd root-interface && ./compile_writer.sh
##########

# MESMER library
libmesmer.a: $(OBJECTS)
	@echo " "
	@echo "Creating MESMER library"
	@echo " "
	ar cr libmesmer.a $(OBJECTS)

# MESMER executable
$(EXE): $(LIBFILES) libmesmer.a
	@echo " "
	@echo "Creating MESMER program"
	@echo " "
	$(F77) main.o -L. -lmesmer -L$(LTDIR)/lib64 $(LTSTRING) $(HANDYG) $(RCLLIB) $(CLLLIB) $(CHAPLIN) $(RCLEXTRA) $(GSLLIBS) -o $(EXE)
