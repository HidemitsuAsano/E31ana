CC=gcc
CXX=g++
#CXX=clang
lCOPT=-DOFFLINE -DDEBUG -DNT -DRAW
CXXOPT = $(COPT) -Wall #-Werror

CONFIGS	= -DUSE_ROOT #-DDEBUG

LD	= $(CXX)

#CFLAGS	= -O2 -g -pg -fPIC $(MACROS) $(CONFIGS)
CFLAGS	= -O2 -g -fPIC $(MACROS) $(CONFIGS)
ifeq ($(shell uname),Darwin)
UNDEFOPT = dynamic_lookup
#SOFLAGS = -dynamiclib -single_module -undefined $(UNDEFOPT) -install_name $(CURDIR)
SOFLAGS = -dynamiclib -single_module -undefined $(UNDEFOPT)
else
SOFLAGS = -shared
endif

#ROOT	= /cern/ROOT/pro
#ROOT	= /w/e15/common/root
ROOT	= $(ROOTSYS)
TMPDIR  = ./tmp

ROOTLIB	= $(shell $(ROOT)/bin/root-config --libs) -lMinuit -lGeom -lEG
ROOTINC	= $(shell $(ROOT)/bin/root-config --cflags)

EXPLLINKLIBS = $(shell $(ROOT)/bin/root-config --libs) -lMinuit -lGeom -lEG
CERN	= /cern/pro

KFINC	= -I/gpfs/home/had/hiasano/source/

INC	= $(ROOTINC) $(KFINC)
MATH	= -lm
#LIBF77	= -lg2c
LIBLINK	= $(MATH) $(ROOTLIB)
#LIBLINK	= $(MATH) $(LIBF77) $(ROOTLIB)

BINDIR	= ../bin
VPATH	= $(BINDIR):$(TMPDIR)

RTLIBDIR = ./lib

.SUFFIXES: .cc .o

$(TMPDIR)/%_Dict.o: $(TMPDIR)/%_Dict.cpp
	$(CXX) -MMD $(CFLAGS) $(INC) $(CXXOPT) -c -o $@ $<

$(TMPDIR)/%.o: %.c
	$(CC) -MMD $(CFLAGS) $(INC) -c -o $@ $<

$(TMPDIR)/%.o: %.cpp
	$(CXX) -MMD $(CFLAGS) $(INC) $(CXXOPT) -c -o $@ $<



OBJS	= 	TKO.o TKO_Dict.o\
		EventHeader.o EventHeader_Dict.o \
		HodoscopeLikeHit.o HodoscopeLikeHit_Dict.o \
		CherenkovLikeHit.o CherenkovLikeHit_Dict.o \
		ChamberLikeHit.o ChamberLikeHit_Dict.o	\
		MTDCHit.o MTDCHit_Dict.o	\
		TDCHit.o TDCHit_Dict.o	\
		ADCHit.o ADCHit_Dict.o	\
		HitMan.o HitMan_Dict.o	\
		BeamLineHitMan.o BeamLineHitMan_Dict.o	\
		CDCHit.o CDCHit_Dict.o	\
		CDSHitMan.o CDSHitMan_Dict.o	\
		ScalerMan.o ScalerMan_Dict.o	\
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
		Display.o Display_Dict.o \
		HelixFit.o HelixFit_Dict.o \
		LineFit.o LineFit_Dict.o \
		CircleFit.o CircleFit_Dict.o \
		CDSTrackingMan.o CDSTrackingMan_Dict.o	\
		CDSTrack.o CDSTrack_Dict.o	\
		HitCluster.o HitCluster_Dict.o	\
		MathTools.o \
		ELossTools.o \
		TrackTools.o \
		GeomTools.o GeomTools_Dict.o \
		BeamLineTrackMan.o BeamLineTrackMan_Dict.o	\
		LocalTrack.o LocalTrack_Dict.o	\
		BLDCClusterMan.o BLDCClusterMan_Dict.o	\
		HodoClusterMan.o HodoClusterMan_Dict.o	\
		BeamSpectrometer.o BeamSpectrometer_Dict.o	\
		DCTimeCorrMan.o DCTimeCorrMan_Dict.o	\
		UshiwakaFieldMapMan.o UshiwakaFieldMapMan_Dict.o	\
		USWKTrack.o USWKTrack_Dict.o	\
		UshiwakaTableMan.o UshiwakaTableMan_Dict.o	\
		EventAlloc.o EventTemp.o	\
		Tools.o Tools_Dict.o	\
		DB_Manager.o DB_Manager_Dict.o	\
		Particle.o Particle_Dict.o	\
		DetectorList.o	DetectorList_Dict.o 
#		SimDataMan.o SimDataMan_Dict.o
#		mtdc_decoder.o mtdc_off_ana.o
#		AnalysisMan.o AnalysisMan_Dict.o	\
#		ProtonArm.o ProtonArm_Dict.o	\
#		mtdc_decoder.o mtdc_off_ana.o
#		Particle.o Particle_Dict.o	\
#		BeamLineChamber.o BeamLineChamber_Dict.o	\
#		Display3D.o Display3D_Dict.o \
#		SDDHit.o SDDHit_Dict.o  \
#		SDDHitCom.o SDDHitCom_Dict.o    \
#		SDDHitMan.o SDDHitMan_Dict.o    \
#		DisplayFADC.o DisplayFADC_Dict.o \

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
		MTDCHit.o MTDCHit_Dict.o	\
		TDCHit.o TDCHit_Dict.o	\
		ADCHit.o ADCHit_Dict.o	\
		HitMan.o HitMan_Dict.o	\
		BeamLineHitMan.o BeamLineHitMan_Dict.o	\
		CDCHit.o CDCHit_Dict.o	\
		CDSHitMan.o CDSHitMan_Dict.o	\
		ScalerMan.o ScalerMan_Dict.o	\
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
		Display.o Display_Dict.o \
		HelixFit.o HelixFit_Dict.o \
		CircleFit.o CircleFit_Dict.o \
		LineFit.o LineFit_Dict.o \
		CDSTrackingMan.o CDSTrackingMan_Dict.o \
		CDSTrack.o CDSTrack_Dict.o	\
		HitCluster.o HitCluster_Dict.o	\
		MathTools.o \
		ELossTools.o \
		GeomTools.o \
		TrackTools.o \
		Tools.o Tools_Dict.o	\
		BeamLineTrackMan.o BeamLineTrackMan_Dict.o	\
		LocalTrack.o LocalTrack_Dict.o	\
		BLDCClusterMan.o BLDCClusterMan_Dict.o	\
		HodoClusterMan.o HodoClusterMan_Dict.o	\
		TransferMatrixMan.o TransferMatrixMan_Dict.o	\
		DCTimeCorrMan.o DCTimeCorrMan_Dict.o	\
		BeamSpectrometer.o BeamSpectrometer_Dict.o	\
		UshiwakaFieldMapMan.o UshiwakaFieldMapMan_Dict.o	\
		UshiwakaTableMan.o UshiwakaTableMan_Dict.o	\
		DataContainer.o DataContainer_Dict.o	\
		USWKTrack.o USWKTrack_Dict.o DB_Manager.o DB_Manager_Dict.o	\
		Particle.o Particle_Dict.o	\
		DetectorList.o	DetectorList_Dict.o 
#		SimDataMan.o SimDataMan_Dict.o 
#		AnalysisMan.o AnalysisMan_Dict.o	\
#		BeamLineChamber.o BeamLineChamber_Dict.o	\
#		Display3D.o Display3D_Dict.o
#		ProtonArm.o ProtonArm_Dict.o	\
#		AnalysisMan.o AnalysisMan_Dict.o	\
#		SDDHit.o SDDHit_Dict.o  \
#		SDDHitCom.o SDDHitCom_Dict.o    \
#		SDDHitMan.o SDDHitMan_Dict.o    \
#		DisplayFADC.o DisplayFADC_Dict.o \

RAWHITSO=	libHit.so
RAWHITLIB=	-L$(RTLIBDIR) -lHit

ALLSO	= libAll.so
ALLLIB	= -L$(RTLIBDIR) -lAll

EVENT	=	Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o

all:
	make link;
	make lib;
#	make lib; make tra; make evcdc;
#	make lib; make ana;make evanaraw;make evdisp;make evtracking;

LIBLINKDIR := ../lib
link:
	@if [ -d $(LIBLINKDIR) ]; then\
		rm -r $(LIBLINKDIR);\
	fi
	@if [ ! -e $(LIBLINKDIR) ]; then\
		ln -s src/lib ../;\
	fi

lib:	$(OBJSA) $(RAWDATASO) $(ALLSO) $(DISPSO) $(OBJSMNA)

$(OBJSA): $(patsubst %,$(TMPDIR)/%,$(OBJS))
	ar r $(OBJSA) $(patsubst %,$(TMPDIR)/%,$(OBJS))

$(OBJSMNA): $(OBJSMN)
	ar r $(OBJSMNA) $(OBJSMN)

$(RAWDATASO): $(patsubst %,$(TMPDIR)/%,$(RAWDATA))
	$(LD) $(SOFLAGS) -O2 $^ -o $(RTLIBDIR)/$@ $(EXPLLINKLIBS)

$(RAWHITSO): $(patsubst %,$(TMPDIR)/%,$(RAWHIT))
	$(LD) $(SOFLAGS) -O2 $^ -o $(RTLIBDIR)/$@ $(EXPLLINKLIBS)

$(ALLSO): $(patsubst %,$(TMPDIR)/%,$(RAWDATA)) $(patsubst %,$(TMPDIR)/%,$(RAWHIT)) $(OBJSMNA)
	$(LD) $(SOFLAGS) -O2 $^ -o $(RTLIBDIR)/$@ $(EXPLLINKLIBS)

$(TMPDIR)/%_Dict.cpp: %.h	
	@echo ------------------------- Making $@
	$(ROOTSYS)/bin/rootcint -f $@ -c $(INCLUDE) $(PWD)/$< LinkDef/$(<D)/$*_LinkDef.h
	@echo ------------------------- Done "_(:3vU)_" Done

SRCS	= $(wildcard *.cpp)
-include $(SRCS:%.cpp=%.d)

evanaraw:  $(TMPDIR)/EventAnalysisRawAll.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisRawAll.o $(LIBLINK) $(OBJSLIB)

evanaxt:  $(TMPDIR)/EventAnalysisXT.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisXT.o $(LIBLINK) $(OBJSLIB)

evanadctimecorr:  $(TMPDIR)/EventAnalysisDCTimeCorr.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisDCTimeCorr.o $(LIBLINK) $(OBJSLIB)

evanatko:  $(TMPDIR)/EventAnalysisTKO.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(TMPDIR)/EventAnalysisTKO.o $(patsubst %,$(TMPDIR)/%,$(EVENT))  $(LIBLINK) $(OBJSLIB)

evanabpc:  $(TMPDIR)/EventAnalysisBPC.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisBPC.o $(LIBLINK) $(OBJSLIB)

evanabldc:  $(TMPDIR)/EventAnalysisBLDC.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisBLDC.o $(LIBLINK) $(OBJSLIB)

evanadc:  $(TMPDIR)/EventAnalysisDC.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisDC.o $(LIBLINK) $(OBJSLIB)

evanacdc:  $(TMPDIR)/EventAnalysisCDC.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisCDC.o $(LIBLINK) $(OBJSLIB)

# source file does not exist. (asano)
#mktree:  $(TMPDIR)/EventAnalysisTree.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
#	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(TMPDIR)/EventAnalysisTree.o $(patsubst %,$(TMPDIR)/%,$(EVENT))  $(LIBLINK) $(OBJSLIB)


include Makefile.user

clean:
	rm -rf $(TARGET) dbtest libe471db.a *.o *.d *_Dict.h *_Dict.cpp $(TMPDIR)/*.o $(TMPDIR)/*.d $(TMPDIR)/*.h *~ LinkDef/*~ $(BINDIR)/*

libclean:
	rm -f $(RTLIBDIR)/lib* libObjs.a libMnGld.a

aclean:
	make clean; make libclean
