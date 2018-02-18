// MyAnalysisNCCalib.cpp

#include "MyAnalysisNCCalib.h"

MyAnalysisNCCalib::MyAnalysisNCCalib(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  Clear();
}

MyAnalysisNCCalib::~MyAnalysisNCCalib()
{
}

void MyAnalysisNCCalib::Clear()
{
}

bool MyAnalysisNCCalib::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  int T0Segment = beam->t0seg();
  double T0Timing = beam->t0time();

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fpro = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* pro = particle->cdsi(it);
    if(vdis>9998||vdis>pro->vdis()){
      vdis = pro->vdis();
      vertex = pro->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fpro=it;
    }
  }
  if(fpro==-1) return true; /* there is no vertex in fiducial volume */
  pCDS* cds = particle->cdsi(fpro);

  // ################# //
  // Neutral selection //
  // ################# //
  int MulBVC=0;
  for(int i=0; i<blMan->nBVC(); i++){
    HodoscopeLikeHit* hit = blMan->BVC(i);
    if(hit->CheckRange()) MulBVC++;
  }
  int MulCVC=0;
  for(int i=0; i<blMan->nCVC(); i++){
    HodoscopeLikeHit* hit = blMan->CVC(i);
    if(hit->CheckRange()) MulCVC++;
  }
  if(MulBVC!=0||MulCVC!=0) return true;

  // ########### //
  // NC analysis //
  // ########### //
  HodoscopeLikeHit* nchit;
  int nhit = 0;
  for(int i=0; i<blMan->nNC(); i++){
    const int cid = CID_NC;
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    HodoscopeLikeHit* hit = blMan->NC(i);
    pNC* nc = new pNC();
    int seg = hit->seg();
    int adc[2]; adc[0]=hit->adcu(); adc[1]=hit->adcd();
    int tdc[2]; tdc[0]=hit->tdcu(); tdc[1]=hit->tdcd();
    double time[2]; time[0]=hit->time(0), time[1] = hit->time(1);
    double ctime[2]; ctime[0]=hit->ctime(0); ctime[1]=hit->ctime(1);
    double ene[2]; ene[0]=hit->ene(0); ene[1]=hit->ene(1);
    double hitpos = hit->hitpos();
    double emean = hit->emean();
    double tmean = hit->tmean();
    double ctmean = hit->ctmean();
    double tsub = hit->tsub();
    double ctsub = hit->ctsub();

    const char* ud[2] = {"u","d"};
    for(int iud=0; iud<2; iud++){
      FillHist(Form("%s%s%d_ADC",name,ud[iud],seg),adc[iud]);
      if(tdc[iud]>0){
        FillHist(Form("%s%s%d_TDC",name,ud[iud],seg),tdc[iud]);
        FillHist(Form("%s%s%d_Time",name,ud[iud],seg),time[iud]);
        FillHist(Form("%s%s%d_CTime",name,ud[iud],seg),ctime[iud]);
      }
    }
    if( hit->CheckRange() ){
      nc->SetSegment(hit->seg());
      nc->SetHitPosition(hit->pos().X(),hit->hitpos(),hit->pos().Z());
      nc->SetTime(hit->ctmean());
      nc->SetEnergy(hit->emean());
      nc->CalcMom(beam,vertex);
      FillHist("FWD_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
      FillHist("FWD_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
      FillHist("FWD_TOFvsMomentum",nc->tof(),nc->mom().Mag());
      FillHist("FWD_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
      double segtmp = (nc->seg()-1)%16+1;
      double laytmp = (nc->seg()-1)/16+1;
      FillHist("FWD_HitSegment",segtmp,laytmp);

      double tmptime = nc->fl()/Const;
      double tof = nc->tof();
      double timeoffs = tof - tmptime;
      double udtof[2]; udtof[0]=tof + (ctime[0]-ctmean); udtof[1]=tof + (ctime[1]-ctmean);
      double udtimeoffs[2]; udtimeoffs[0]=udtof[0]-tmptime; udtimeoffs[1]=udtof[1]-tmptime;

      const char* ud[2] = {"u","d"};
      for(int iud=0; iud<2; iud++){
        FillHist(Form("%s%s%d_ADCwithTDC",name,ud[iud],seg),adc[iud]);
        FillHist(Form("%s%s%d_dE",name,ud[iud],seg),ene[iud]);
        FillHist(Form("%s%s%d_TimeOffset",name,ud[iud],seg),udtimeoffs[iud]);
        FillHist(Form("%s%s%d_ADCvsTDC",name,ud[iud],seg),adc[iud],tdc[iud]);
        FillHist(Form("%s%s%d_dEvsTime",name,ud[iud],seg),ene[iud],time[iud]);
        FillHist(Form("%s%s%d_dEvsCTime",name,ud[iud],seg),ene[iud],ctime[iud]);
        FillHist(Form("%s%s%d_dEMeanvsTOF",name,ud[iud],seg),emean,udtof[iud]);
        FillHist(Form("%s%s%d_dEMeanvsCTime",name,ud[iud],seg),emean,ctime[iud]);
        FillHist(Form("%s%s%d_dEMeanvsTimeOffset",name,ud[iud],seg),emean,udtimeoffs[iud]);
      }
      FillHist(Form("%s%d_dEMean",name,seg),emean);
      FillHist(Form("%s%d_TimeMean",name,seg),tmean);
      FillHist(Form("%s%d_TimeOffset",name,seg),timeoffs);
      FillHist(Form("%s%d_TimeSub",name,seg),tsub);
      FillHist(Form("%s%d_CTimeMean",name,seg),ctmean);
      FillHist(Form("%s_CTimeMean",name),ctmean);
      FillHist(Form("%s%d_EventvsCTimeMean",name,seg),header->ev(),ctmean);
      FillHist(Form("%s_EventvsCTimeMean",name),header->ev(),ctmean);
      FillHist(Form("%s%d_CTimeSub",name,seg),ctsub);
      FillHist(Form("%s%d_TOF",name,seg),tof);
      FillHist(Form("%s_TOF",name),tof);
      FillHist(Form("%s_TOFvsT0Segment",name),tof,T0Segment);
      FillHist(Form("%s_TOFvsCTimeMean",name),tof,ctmean);
      FillHist(Form("%s_TOFvsT0CTimeMean",name),tof,T0Timing);
      FillHist(Form("%s%d_Position",name,seg),hitpos);
      FillHist(Form("%s%d_dEMeanvsTOF",name,seg),emean,tof);
      FillHist(Form("%s%d_dEMeanvsCTimeMean",name,seg),emean,ctmean);
      FillHist(Form("%s%d_dEMeanvsTimeOffset",name,seg),emean,timeoffs);
      FillHist(Form("%s_HitPattern",name),seg);
      nhit++;
      for( int j=0; j<20; j++ ){
        int val = header->pattern(j);
        if( 0<val ){
          FillHist(Form("%s_TriggerPatternvsTOF",name),j,tof);
        }
      }
    }
    delete nc;
  }
  FillHist(Form("NC_Multiplicity"),nhit);

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
      FillHist(Form("Trigger%02d_TDC",i),val);
    }
  }

  return true;

}

bool MyAnalysisNCCalib::FillHist(TString name, double val1)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1){
    h1 -> Fill(val1);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisNCCalib::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisNCCalib::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisNCCalib::Initialize ###" << std::endl;

  DetectorList *dlist=DetectorList::GetInstance();
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaNCCalib");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  rtFile -> cd();

  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH2F( "FWD_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",16,-160,160,15,-75,75);
  new TH2F( "FWD_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  // Missing mass p(K-,n)X
  std::cout << "Define Histograms for p(K-,n)X M.M." << std::endl;
  new TH1F( "MM_Nf", "p(K^{-},n)X missing mass;MM_{p}(n_{f}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  // Invariant mass with FWD
  std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
  new TH1F( "IM_Nfpiplus", "Invariant mass (n_{f}#pi^{+});IM(n_{f}#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  new TH1F( "IM_Nfpiminus", "Invariant mass (n_{f}#pi^{-});IM(n_{f}#pi^{-}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );

  /* Hodoscope */
  std::cout << "Define Histograms for Hodoscopes" << std::endl;
  const int nhodo = 1;
  int hodoid[nhodo]={
    CID_NC
  };
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    for( int seg=1; seg<=nsegs; seg++ ){
      const char* ud[2] = {"u","d"};
      const char* udname[2] = {"up","down"};
      for(int iud=0; iud<2; iud++){
        new TH1F( Form("%s%s%d_ADC",name,ud[iud],seg),   Form("ADC %s-%d(%s);ADC ch.;Counts",name,seg,udname[iud]),    2000,    0, 4000 );
        new TH1F( Form("%s%s%d_ADCwithTDC",name,ud[iud],seg), Form("ADC with TDC %s-%d(%s);ADC ch.;Counts",name,seg,udname[iud]),  2000,    0, 4000 );
        new TH1F( Form("%s%s%d_dE",name,ud[iud],seg),   Form("dE %s-%d(%s);dE (MeV);Counts",name,seg,udname[iud]),    1000,    0, 100 );
        new TH1F( Form("%s%s%d_TDC",name,ud[iud],seg),   Form("TDC %s-%d(%s);TDC ch.;Counts",name,seg,udname[iud]),    2000,    0, 4000 );
        new TH1F( Form("%s%s%d_Time",name,ud[iud],seg),   Form("Time %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
        new TH1F( Form("%s%s%d_CTime",name,ud[iud],seg),   Form("CTime %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
        new TH1F( Form("%s%s%d_TimeOffset",name,ud[iud],seg),   Form("Time Offset %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
        new TH2F( Form("%s%s%d_ADCvsTDC",name,ud[iud],seg),   Form("ADC TDC corr. %s-%d(%s);ADC ch.;TDC ch.",name,seg,udname[iud]),     400,    0, 4000,  400,    0, 4000 );
        new TH2F( Form("%s%s%d_dEvsTime",name,ud[iud],seg),   Form("dE Time corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 100,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEvsCTime",name,ud[iud],seg),   Form("dE CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 100,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEMeanvsTOF",name,ud[iud],seg),   Form("dEMean TOF corr. %s-%d(%s);dE (MeV);TOF (ns)",name,seg,udname[iud]),     1000,    0, 100,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEMeanvsCTime",name,ud[iud],seg),   Form("dEMean CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 100,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEMeanvsTimeOffset",name,ud[iud],seg),   Form("dEMean Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 100,  2000,    -100, 100 );
      }
      new TH1F( Form("%s%d_dEMean",name,seg),   Form("Mean dE %s-%d;dE (MeV);Counts",name,seg),    1000,    0, 100 );
      new TH1F( Form("%s%d_TimeMean",name,seg),   Form("Mean Time %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_TimeSub",name,seg),   Form("Subtract Time %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_CTimeMean",name,seg),   Form("Mean CTime %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_CTimeSub",name,seg),   Form("Subtract CTime %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_TOF",name,seg),   Form("TOF between T0 and %s-%d;Time of flight T0-%s%d (ns);Counts",name,seg,name,seg),    3000,    -50, 100 );
      new TH1F( Form("%s%d_Position",name,seg),   Form("Hit position %s-%d;Position (cm);Counts",name,seg),    1000,    -50, 50 );
      new TH2F( Form("%s%d_dEMeanvsTOF",name,seg),   Form("dEMean TOF corr. %s-%d;dE (MeV);TOF (ns)",name,seg),     1000,    0, 100,  2000,    -100, 100 );
      new TH2F( Form("%s%d_dEMeanvsCTimeMean",name,seg),   Form("dEMean CTimeMean corr. %s-%d;dE (MeV);Time (ns)",name,seg),     1000,    0, 100,  2000,    -100, 100 );
      new TH2F( Form("%s%d_dEMeanvsTimeOffset",name,seg),   Form("dEMean Time Offset corr. %s-%d;dE (MeV);Time (ns)",name,seg),     1000,    0, 100,  2000,    -100, 100 );
    }
    new TH1F( Form("%s_HitPattern",name), Form("Hit pattern %s;%s Segment;Counts",name,name), nsegs, 0.5, nsegs+0.5);
    new TH1F( Form("%s_Multiplicity",name), Form("Multiplicity %s;Multiplicity;Counts",name), nsegs+1, -0.5, nsegs+0.5 );
    new TH1F( Form("%s_CTimeMean",name),   Form("Mean CTime %s;Time (ns);Counts",name),    2000,    -100, 100 );
    new TH1F( Form("%s_TOF",name),   Form("TOF between T0 and %s;Time of flight T0-%s (ns);Counts",name,name),    3000,    -50, 100 );
    new TH2F( Form("%s_TOFvsT0Segment",name),   Form("TOF between T0 and %s vs. T0 segment;Time of flight T0-%s (ns);T0 segment",name,name),    3000,    -50, 100, 5, 0.5, 5.5 );
    new TH2F( Form("%s_TOFvsCTimeMean",name),   Form("TOF between T0 and %s vs. CTimeMean;Time of flight T0-%s (ns);Time (ns)",name,name),    3000,    -50, 100, 2000, -100, 100 );
    new TH2F( Form("%s_TOFvsT0CTimeMean",name),   Form("TOF between T0 and %s vs. T0 CTimeMean;Time of flight T0-%s (ns);T0 Time (ns)",name,name),    3000,    -50, 100, 200, -10, 10 );
    TH2F* h2 = new TH2F( Form("%s_TriggerPatternvsTOF",name),   Form("TriggerPattern vs. TOF between T0 and %s;;Time of flight T0-%s (ns)",name,name), 18, -0.5, 17.5 ,    3000,    -50, 100);
    h2->GetXaxis()->SetBinLabel(1,"All");
    h2->GetXaxis()->SetBinLabel(2,"Beam");
    h2->GetXaxis()->SetBinLabel(3,"Kaon");
    h2->GetXaxis()->SetBinLabel(4,"KCDH1f");
    h2->GetXaxis()->SetBinLabel(5,"Pion");
    h2->GetXaxis()->SetBinLabel(6,"Proton");
    h2->GetXaxis()->SetBinLabel(7,"KCDH1");
    h2->GetXaxis()->SetBinLabel(8,"KCDH2");
    h2->GetXaxis()->SetBinLabel(9,"PivBVC");
    h2->GetXaxis()->SetBinLabel(10,"PiCDH1");
    h2->GetXaxis()->SetBinLabel(11,"PiCDH2");
    h2->GetXaxis()->SetBinLabel(12,"Kf");
    h2->GetXaxis()->SetBinLabel(13,"1stMix");
    h2->GetXaxis()->SetBinLabel(14,"Charged");
    h2->GetXaxis()->SetBinLabel(15,"Neutral");
    h2->GetXaxis()->SetBinLabel(16,"Cosmic");
    h2->GetXaxis()->SetBinLabel(17,"Reject");
    h2->GetXaxis()->SetBinLabel(18,"SIM");
  }

}
