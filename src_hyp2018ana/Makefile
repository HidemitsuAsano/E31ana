CC=gcc
CXX=g++
lCOPT=-DOFFLINE -DDEBUG -DNT -DRAW
CXXOPT = $(COPT) -Wall #-Werror

CONFIGS	= -DUSE_ROOT #-DDEBUG

LD	= $(CXX)

CFLAGS	= -O4 -g -fPIC $(MACROS) $(CONFIGS)
ifeq ($(shell uname),Darwin)
UNDEFOPT = dynamic_lookup
SOFLAGS = -dynamiclib -single_module -undefined $(UNDEFOPT)
else
SOFLAGS = -shared
endif

ROOT	= $(ROOTSYS)
TMPDIR  = ./tmp

ROOTLIB	= $(shell $(ROOT)/bin/root-config --libs) -lMinuit -lGeom -lEG
ROOTINC	= $(shell $(ROOT)/bin/root-config --cflags)

EXPLLINKLIBS = $(shell $(ROOT)/bin/root-config --libs) -lGeom -lEG
CERN	= /cern/pro

KFINC   = -I/w/e15/common/KinFitter/
KFLIB   = -L/w/e15/common/KinFitter/lib/Linux26SL6_x86_64_gcc447/5.34.13/ -lTKinFitter

INC	= $(ROOTINC) $(KFINC)
MATH	= -lm

LIBLINK	= $(MATH) $(ROOTLIB)

BINDIR	= ../bin
VPATH	= $(BINDIR):$(TMPDIR)

RTLIBDIR = ./lib

.SUFFIXES: .cc .o


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
		GeomTools.o GeomTools_Dict.o \
		SimTools.o SimTools_Dict.o \
		TrackTools.o \
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
		DetectorList.o	DetectorList_Dict.o	\
		SlewingData.o SlewingData_Dict.o	\
		NpipiData.o NpipiData_Dict.o	\
		NppimData.o NppimData_Dict.o	\
		SimDataMan.o SimDataMan_Dict.o \
		KinFitMan.o KinFitMan_Dict.o	\
		KnuclRootData.o KnuclRootData_Dict.o \
		ComCrossSectionTable.o ComCrossSectionTable_Dict.o \
                ChargeParticle.o ChargeParticle_Dict.o    \
		ProtonArm.o ProtonArm_Dict.o	\
		BeamInfo.o BeamInfo_Dict.o      \
		CDSInfo.o CDSInfo_Dict.o        \
		CDS2Info.o CDS2Info_Dict.o      \
		ForwardNeutralInfo.o ForwardNeutralInfo_Dict.o  \
		ForwardChargeInfo.o ForwardChargeInfo_Dict.o    \
		HistInfo.o HistInfo_Dict.o        \
		AnaInfo.o AnaInfo_Dict.o	\
		SlewingTools.o SlewingTools_Dict.o \
		TemplateFitter.o TemplateFitter_Dict.o \

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
		SimTools.o \
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
		DetectorList.o	DetectorList_Dict.o	\
		SlewingData.o SlewingData_Dict.o	\
		NpipiData.o NpipiData_Dict.o	\
		NppimData.o NppimData_Dict.o	\
		SimDataMan.o SimDataMan_Dict.o \
		KinFitMan.o KinFitMan_Dict.o	\
		KnuclRootData.o KnuclRootData_Dict.o \
		ComCrossSectionTable.o ComCrossSectionTable_Dict.o \
		ProtonArm.o ProtonArm_Dict.o	\
		ChargeParticle.o ChargeParticle_Dict.o  \
		BeamInfo.o BeamInfo_Dict.o      \
		CDSInfo.o CDSInfo_Dict.o        \
		CDS2Info.o CDS2Info_Dict.o      \
		ForwardNeutralInfo.o ForwardNeutralInfo_Dict.o  \
		ForwardChargeInfo.o ForwardChargeInfo_Dict.o    \
		AnaInfo.o AnaInfo_Dict.o	\
		HistInfo.o HistInfo_Dict.o        \
		SlewingTools.o SlewingTools_Dict.o \
		TemplateFitter.o TemplateFitter_Dict.o \

RAWHITSO=	libHit.so
RAWHITLIB=	-L$(RTLIBDIR) -lHit

ALLSO	= libAll.so
ALLLIB	= -L$(RTLIBDIR) -lAll

EVENT	=	Main.o File.o FileType.o VEvent.o BLEvent.o rawread3.o

all:
	make link;
	make lib;

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

%_Dict.cpp: %.h	
	@echo ------------------------- Making $@
	$(ROOTSYS)/bin/rootcint -f $@ -c $(INCLUDE) $< LinkDef/$(<D)/$*_LinkDef.h
	@echo ------------------------- Done "\(^o^)/" Done

SRCS	= $(wildcard *.cpp)
-include $(SRCS:%.cpp=%.d)

evtracking: $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysis_Tracking.o $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysis_Tracking.o $(LIBLINK) $(OBJSLIB)

evanaraw:  $(TMPDIR)/EventAnalysisRawAll.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisRawAll.o $(LIBLINK) $(OBJSLIB)

evanatko:  $(TMPDIR)/EventAnalysisTKO.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(TMPDIR)/EventAnalysisTKO.o $(patsubst %,$(TMPDIR)/%,$(EVENT))  $(LIBLINK) $(OBJSLIB)

mktree:  $(TMPDIR)/EventAnalysisTree.o $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@  $(TMPDIR)/EventAnalysisTree.o $(patsubst %,$(TMPDIR)/%,$(EVENT))  $(LIBLINK) $(OBJSLIB)

include Makefile.user
include Makefile.calib

clean:
	rm -f $(TARGET) dbtest libe471db.a *.o *.d *_Dict.h *_Dict.cpp $(TMPDIR)/*.o $(TMPDIR)/*.d $(TMPDIR)/*.h *~ LinkDef/*~ $(BINDIR)/*~
libclean:
	rm -f $(RTLIBDIR)/lib* libObjs.a libMnGld.a

aclean:
	make clean; make libclean
