evall:
	make lib; make evtracking;

evdisp: $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventDisplay.o $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventDisplay.o $(LIBLINK) $(OBJSLIB)

ev1: $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysis1.o $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysis1.o $(LIBLINK) $(OBJSLIB)

momentum: $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisMomentum.o $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventAnalysisMomentum.o $(LIBLINK) $(OBJSLIB)

xt: $(TMPDIR)/UserXT.o $(ALLSO)
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(TMPDIR)/UserXT.o $(LIBLINK) $(RTLIBDIR)/$(ALLSO)

mtdc_decoder: $(TMPDIR)/mtdc_decoder.o $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(TMPDIR)/mtdc_decoder.o $(LIBLINK) $(OBJSLIB)

mtdc_off_ana: $(TMPDIR)/mtdc_off_ana.o $(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(TMPDIR)/mtdc_off_ana.o $(LIBLINK) $(OBJSLIB)

sim:	$(TMPDIR)/UserSimDatG4.o 	 \
	$(TMPDIR)/SimDataMan.o		 \
	$(patsubst %,$(TMPDIR)/%,$(OBJS))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(LIBLINK) $(OBJSLIB) $(RTLIBDIR)/KnuclRootData_cc.so

evanaBPC: $(patsubst %,$(TMPDIR)/%,$(EVENT))	       \
	  $(TMPDIR)/EventAnalysisBPC.o	 \
	  $(TMPDIR)/HistManBPC.o	 \
	  $(TMPDIR)/DCEffMan.o	 \
	  $(TMPDIR)/BPCTrackMan.o	\
	  $(TMPDIR)/BLDCTrack.o	\
	  $(TMPDIR)/XYTrack.o	\
	  $(patsubst %,$(TMPDIR)/%,$(OBJS))
	  $(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	  $(TMPDIR)/EventAnalysisBPC.o \
	  $(TMPDIR)/HistManBPC.o	 \
	  $(TMPDIR)/DCEffMan.o	 \
	  $(TMPDIR)/BPCTrackMan.o	\
	  $(TMPDIR)/BLDCTrack.o	\
	  $(TMPDIR)/XYTrack.o	\
	  $(LIBLINK) $(OBJSLIB)

evanaBPC2: $(patsubst %,$(TMPDIR)/%,$(EVENT))	       \
	  $(TMPDIR)/EventAnalysisBPC2.o	 \
	  $(TMPDIR)/HistManBPC2.o	 \
	  $(TMPDIR)/DCEffMan.o	 \
	  $(TMPDIR)/BPCTrackMan.o	\
	  $(TMPDIR)/BLDCTrack.o	\
	  $(TMPDIR)/XYTrack.o	\
	  $(patsubst %,$(TMPDIR)/%,$(OBJS))
	  $(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	  $(TMPDIR)/EventAnalysisBPC2.o \
	  $(TMPDIR)/HistManBPC2.o	 \
	  $(TMPDIR)/DCEffMan.o	 \
	  $(TMPDIR)/BPCTrackMan.o	\
	  $(TMPDIR)/BLDCTrack.o	\
	  $(TMPDIR)/XYTrack.o	\
	  $(LIBLINK) $(OBJSLIB)

evanaFDC: $(patsubst %,$(TMPDIR)/%,$(EVENT))	       \
	  $(TMPDIR)/EventAnalysisFDC.o	 \
	  $(TMPDIR)/HistManFDC.o	 \
	  $(TMPDIR)/DCEffMan.o	 \
	  $(TMPDIR)/BPCTrackMan.o	\
	  $(TMPDIR)/BLDCTrack.o	\
	  $(TMPDIR)/XYTrack.o	\
	  $(patsubst %,$(TMPDIR)/%,$(OBJS))
	  $(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	  $(TMPDIR)/EventAnalysisFDC.o \
	  $(TMPDIR)/HistManFDC.o	 \
	  $(TMPDIR)/DCEffMan.o	 \
	  $(TMPDIR)/BPCTrackMan.o	\
	  $(TMPDIR)/BLDCTrack.o	\
	  $(TMPDIR)/XYTrack.o	\
	  $(LIBLINK) $(OBJSLIB)

evdispBPC: $(patsubst %,$(TMPDIR)/%,$(EVENT))	       \
	   $(TMPDIR)/EventDispBPC.o	 \
	   $(TMPDIR)/DispBPC.o	 \
	   $(TMPDIR)/BPCTrackMan.o	\
	   $(TMPDIR)/BLDCTrack.o	\
	   $(TMPDIR)/XYTrack.o	\
	   $(patsubst %,$(TMPDIR)/%,$(OBJS))
	   $(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	   $(TMPDIR)/EventDispBPC.o \
	   $(TMPDIR)/DispBPC.o	 \
	   $(TMPDIR)/BPCTrackMan.o	\
	   $(TMPDIR)/BLDCTrack.o	\
	   $(TMPDIR)/XYTrack.o	\
	   $(LIBLINK) $(OBJSLIB)

evdispBLC: $(patsubst %,$(TMPDIR)/%,$(EVENT)) $(TMPDIR)/EventDispBLC2.o $(patsubst %,$(TMPDIR)/%,$(OBJS)) $(TMPDIR)/DisplayBLC.o
	   $(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(LIBLINK) $(OBJSLIB)

evanaCalib: $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	    $(TMPDIR)/EventAnalysisCalib.o     \
	    $(TMPDIR)/HistManBeamAna.o	       \
	    $(TMPDIR)/HistManBLDC.o	       \
	    $(TMPDIR)/HistManFC.o	       \
	    $(TMPDIR)/MyTools.o	       \
	    $(patsubst %,$(TMPDIR)/%,$(OBJS)) 
	    $(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	    $(TMPDIR)/EventAnalysisCalib.o   \
	    $(TMPDIR)/HistManBeamAna.o	     \
	    $(TMPDIR)/HistManBLDC.o	     \
	    $(TMPDIR)/HistManFC.o	     \
	    $(TMPDIR)/MyTools.o	       \
	    $(LIBLINK) $(OBJSLIB)

evanaCalibwCDS: $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	     $(TMPDIR)/EventAnalysisCalibwCDS.o     \
	     $(TMPDIR)/HistManBeamAna.o	       \
	     $(TMPDIR)/HistManBLDC.o	       \
	     $(TMPDIR)/HistManFC.o	       \
	     $(TMPDIR)/HistManCDS.o	       \
	     $(TMPDIR)/MyTools.o	       \
	     $(patsubst %,$(TMPDIR)/%,$(OBJS)) 
	     $(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	     $(TMPDIR)/EventAnalysisCalibwCDS.o   \
	     $(TMPDIR)/HistManBeamAna.o	     \
	     $(TMPDIR)/HistManBLDC.o	     \
	     $(TMPDIR)/HistManFC.o	       \
	     $(TMPDIR)/HistManCDS.o	       \
	     $(TMPDIR)/MyTools.o	       \
	     $(LIBLINK) $(OBJSLIB)

evanaAnaData: $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	    $(TMPDIR)/EventAnalysisAnaData.o     \
	    $(TMPDIR)/MyTools.o	       \
	    $(patsubst %,$(TMPDIR)/%,$(OBJS)) 
	    $(CXX) $(CFLAGS) -o $(BINDIR)/$@ $(patsubst %,$(TMPDIR)/%,$(EVENT)) \
	    $(TMPDIR)/EventAnalysisAnaData.o     \
	    $(TMPDIR)/MyTools.o	       \
	    $(LIBLINK) $(OBJSLIB)