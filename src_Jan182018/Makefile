CC=gcc
CXX=g++
lCOPT=-DOFFLINE -DDEBUG -DNT -DRAW
CXXOPT = $(COPT) -Wall

CONFIGS	= -DUSE_ROOT #-DDEBUG

LD	= $(CXX)

#CFLAGS	= -O2 -g -pg -fPIC $(MACROS) $(CONFIGS)
CFLAGS	= -O4 -g -fPIC $(MACROS) $(CONFIGS)
ifeq ($(shell uname),Darwin)
UNDEFOPT = dynamic_lookup
SOFLAGS = -dynamiclib -single_module -undefined $(UNDEFOPT) -install_name $(CURDIR)
else
SOFLAGS = -shared
endif

#ROOT	= /cern/ROOT/pro
#ROOT	= /w/e15/common/root
ROOT	= $(ROOTSYS)

ROOTLIB	= $(shell $(ROOT)/bin/root-config --libs) -lMinuit
ROOTINC	= $(shell $(ROOT)/bin/root-config --cflags)

CERN	= /cern/pro

INC	= $(ROOTINC)
MATH	= -lm
#LIBF77	= -lg2c
LIBLINK	= $(MATH) $(ROOTLIB)
#LIBLINK	= $(MATH) $(LIBF77) $(ROOTLIB)

BINDIR	= ../bin
VPATH	= $(BINDIR)

RTLIBDIR = ./lib

.SUFFIXES: .cc .o

.c.o:
	$(CC) -MMD $(CFLAGS) $(INC) -c $<

.cpp.o:
	$(CXX) -MMD $(CFLAGS) $(INC) $(CXXOPT) -c $<

OBJS	= 	TKO.o TKO_Dict.o\
		EventHeader.o EventHeader_Dict.o \
		HodoscopeLikeHit.o HodoscopeLikeHit_Dict.o \
		CherenkovLikeHit.o CherenkovLikeHit_Dict.o \
		ChamberLikeHit.o ChamberLikeHit_Dict.o	\
		BeamLineHitMan.o BeamLineHitMan_Dict.o	\
		CDCHit.o CDCHit_Dict.o	\
		CDSHitMan.o CDSHitMan_Dict.o	\
		ScalerMan.o ScalerMan_Dict.o	\
		SDDHit.o SDDHit_Dict.o  \
		SDDHitCom.o SDDHitCom_Dict.o    \
		SDDHitMan.o SDDHitMan_Dict.o    \
		GlobalVariables.o GlobalVariables_Dict.o \
		TransferMatrixMan.o TransferMatrixMan_Dict.o	\
		ConfMan.o ConfMan_Dict.o	\
		CounterMapMan.o CounterMapMan_Dict.o	\
		GainMapMan.o GainMapMan_Dict.o	\
		GateMapMan.o GateMapMan_Dict.o	\
		SlewingMapMan.o SlewingMapMan_Dict.o \
		GeomMapMan.o GeomMapMan_Dict.o	\
		XTMapMan.o XTMapMan_Dict.o	\
		CDCWireMapMan.o CDCWireMapMan_Dict.o	\
		BLDCWireMapMan.o BLDCWireMapMan_Dict.o	\
		ScalerMapMan.o ScalerMapMan_Dict.o	\
		ReslMapMan.o ReslMapMan_Dict.o	\
		CDSFittingParamMan.o CDSFittingParamMan_Dict.o	\
		BLDCFittingParamMan.o BLDCFittingParamMan_Dict.o	\
		FADCParamMan.o FADCParamMan_Dict.o	\
		SimDataMan.o SimDataMan_Dict.o \
		Display.o Display_Dict.o \
		DisplayFADC.o DisplayFADC_Dict.o \
		HelixFit.o HelixFit_Dict.o \
		CircleFit.o CircleFit_Dict.o \
		CDSTrackingMan.o CDSTrackingMan_Dict.o	\
		CDSTrack.o CDSTrack_Dict.o	\
		HitCluster.o HitCluster_Dict.o	\
		MathTools.o MathTools_Dict.o	\
		BeamLineTrackMan.o BeamLineTrackMan_Dict.o	\
		LinearTrack.o LinearTrack_Dict.o	\
		BeamSpectrometer.o BeamSpectrometer_Dict.o	\
		EventAlloc.o EventTemp.o \
		mtdc_decoder.o mtdc_off_ana.o
#		Particle.o Particle_Dict.o	\
#		AnalysisMan.o AnalysisMan_Dict.o	\
#		BeamLineChamber.o BeamLineChamber_Dict.o	\
#		Display3D.o Display3D_Dict.o \

OBJSMN  =       mngld.o
OBJSMNA =       libMnGld.a
OBJSA	=	libObjs.a
OBJSLIB	=	-L. -lObjs -lMnGld

RAWDATA	= 	TKO.o TKO_Dict.o\
		EventHeader.o EventHeader_Dict.o
RAWDATASO = 	libData.so
RAWDATALIB= 	-L$(RTLIBDIR) -lData

RAWHIT	=	HodoscopeLikeHit.o HodoscopeLikeHit_Dict.o \
		CherenkovLikeHit.o CherenkovLikeHit_Dict.o \
		ChamberLikeHit.o ChamberLikeHit_Dict.o	\
		BeamLineHitMan.o BeamLineHitMan_Dict.o	\
		CDCHit.o CDCHit_Dict.o	\
		CDSHitMan.o CDSHitMan_Dict.o	\
		ScalerMan.o ScalerMan_Dict.o	\
		SDDHit.o SDDHit_Dict.o  \
		SDDHitCom.o SDDHitCom_Dict.o    \
		SDDHitMan.o SDDHitMan_Dict.o    \
		GlobalVariables.o GlobalVariables_Dict.o \
		ConfMan.o ConfMan_Dict.o	\
		CounterMapMan.o CounterMapMan_Dict.o	\
		GainMapMan.o GainMapMan_Dict.o	\
		GateMapMan.o GateMapMan_Dict.o	\
		SlewingMapMan.o SlewingMapMan_Dict.o	\
		GeomMapMan.o GeomMapMan_Dict.o	\
		XTMapMan.o XTMapMan_Dict.o	\
		CDCWireMapMan.o CDCWireMapMan_Dict.o	\
		BLDCWireMapMan.o BLDCWireMapMan_Dict.o	\
		ScalerMapMan.o ScalerMapMan_Dict.o	\
		ReslMapMan.o ReslMapMan_Dict.o	\
		CDSFittingParamMan.o CDSFittingParamMan_Dict.o	\
		BLDCFittingParamMan.o BLDCFittingParamMan_Dict.o	\
		FADCParamMan.o FADCParamMan_Dict.o	\
		SimDataMan.o SimDataMan_Dict.o \
		Display.o Display_Dict.o \
		DisplayFADC.o DisplayFADC_Dict.o \
		HelixFit.o HelixFit_Dict.o \
		CircleFit.o CircleFit_Dict.o \
		CDSTrackingMan.o CDSTrackingMan_Dict.o \
		CDSTrack.o CDSTrack_Dict.o	\
		HitCluster.o HitCluster_Dict.o	\
		MathTools.o MathTools_Dict.o	\
		BeamLineTrackMan.o BeamLineTrackMan_Dict.o	\
		LinearTrack.o LinearTrack_Dict.o	\
		TransferMatrixMan.o TransferMatrixMan_Dict.o	\
		BeamSpectrometer.o BeamSpectrometer_Dict.o	
#		Particle.o Particle_Dict.o	\
#		AnalysisMan.o AnalysisMan_Dict.o	\
#		BeamLineChamber.o BeamLineChamber_Dict.o	\
#		Display3D.o Display3D_Dict.o

RAWHITSO=	libHit.so
RAWHITLIB=	-L$(RTLIBDIR) -lHit

ALLSO	= libAll.so
ALLLIB	= -L$(RTLIBDIR) -lAll

all:
	make lib;
#	make lib; make ana;make evanaraw;make evdisp;make evtracking;

lib:	$(OBJSA) $(RAWDATASO) $(ALLSO) $(DISPSO) $(OBJSMNA)

$(OBJSA): $(OBJS)
	ar r $(OBJSA) $(OBJS)

$(OBJSMNA): $(OBJSMN)
	ar r $(OBJSMNA) $(OBJSMN)

$(RAWDATASO): $(RAWDATA)
	$(LD) $(SOFLAGS) -O2 $^ -o $(RTLIBDIR)/$@ $(EXPLLINKLIBS)

$(RAWHITSO): $(RAWHIT)
	$(LD) $(SOFLAGS) -O2 $^ -o $(RTLIBDIR)/$@ $(EXPLLINKLIBS)

$(ALLSO): $(RAWDATA) $(RAWHIT) $(OBJSMNA)
	$(LD) $(SOFLAGS) -O2 $^ -o $(RTLIBDIR)/$@ $(EXPLLINKLIBS)

%_Dict.cpp: %.h
	@echo ------------------------- Making $@
	$(ROOTSYS)/bin/rootcint -f $@ -c $(INCLUDE) $< LinkDef/$(<D)/$*_LinkDef.h
	@echo ------------------------- Done "\(^o^)/" Done

SRCS	= $(wildcard *.cpp)
-include $(SRCS:%.cpp=%.d)

#mktree: rawread.o mktree.o $(OBJS)
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ rawread.o mktree.o $(LIBLINK) $(OBJSLIB)

#mktree: rawread2.o mktree.o $(OBJS)
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ rawread2.o mktree.o $(LIBLINK) $(OBJSLIB)

#cds: UserCDS.o $(ALLSO)
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ UserCDS.o $(LIBLINK) $(RTLIBDIR)/$(ALLSO)

test: UserTest.o $(ALLSO)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ UserTest.o $(LIBLINK) $(RTLIBDIR)/$(ALLSO)

#cdstree: Main.o File.o FileType.o VEvent.o $(OBJS)
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o VEvent.o $(LIBLINK) $(OBJSLIB)

#tkotree: Main.o File.o FileType.o TKOEvent.o $(OBJS)
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o TKOEvent.o $(LIBLINK) $(OBJSLIB)

evanaraw: Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisRawAll.o $(OBJS)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisRawAll.o $(LIBLINK) $(OBJSLIB)
evananc: Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisNCcalib.o $(OBJS)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisNCcalib.o $(LIBLINK) $(OBJSLIB)


evanatko: Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisTKO.o $(OBJS)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisTKO.o $(LIBLINK) $(OBJSLIB)


evanaall: Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisAll.o $(OBJS)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisAll.o $(LIBLINK) $(OBJSLIB)

bhdt0: Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisBHDT0.o $(OBJS)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisBHDT0.o $(LIBLINK) $(OBJSLIB)

t0tof: Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisT0TOF.o $(OBJS)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisT0TOF.o $(LIBLINK) $(OBJSLIB)

t0cdh: Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisT0CDH.o $(OBJS)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o EventAnalysisT0CDH.o $(LIBLINK) $(OBJSLIB)

include Makefile.user

clean:
	rm -f $(TARGET) dbtest libe471db.a *.o *.d *_Dict.h *_Dict.cpp *~ LinkDef/*~

libclean:
	rm -f $(RTLIBDIR)/lib* libObjs.a libMnGld.a

aclean:
	make clean; make libclean
