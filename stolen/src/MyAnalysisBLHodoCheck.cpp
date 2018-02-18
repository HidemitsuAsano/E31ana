// MyAnalysisBLHodoCheck.cpp

#include "MyAnalysisBLHodoCheck.h"

MyAnalysisBLHodoCheck::MyAnalysisBLHodoCheck(TFile* rt, ConfMan* conf)
{
	Initialize(conf);
	CutCondition();
	Clear();
}

MyAnalysisBLHodoCheck::~MyAnalysisBLHodoCheck()
{
	Clear();
	rtFile->cd();
	rtFile->Write();
	rtFile->Close();
}

void MyAnalysisBLHodoCheck::Clear()
{
}

bool MyAnalysisBLHodoCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{

	rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  if(!header->IsTrig(Trig_Kaon)){
    return true;
  }
  /* ========== Trigger selection ========== */
  /* ======================================= */

  /* Event selection */
  double T0Timing = -9999.; 
  double T0udE = -9999.;
  double T0ddE = -9999.;
  double T0mdE = -9999.;
  int T0Segment = -1;
  int MulT0 = 0;
  for(int i=0; i<blMan->nT0(); i++){
    if(blMan->T0(i)->CheckRange()){
      MulT0++;
      T0Timing = blMan->T0(i)->ctmean();
      T0Segment = blMan->T0(i)->seg();
      T0udE = blMan->T0(i)->ene(0);
      T0ddE = blMan->T0(i)->ene(1);
      T0mdE = blMan->T0(i)->emean();
    }
  }
  if(MulT0!=1){
    return true;
  }
  //if(T0Segment!=3){
  //  return true;
  //}
  int MulBHD = 0;
  double BHDTiming = -9999.;
  for(int i=0; i<blMan->nBHD(); i++){
    if(blMan->BHD(i)->CheckRange()){
      if(6<=blMan->BHD(i)->seg()&&blMan->BHD(i)->seg()<=15){
        BHDTiming = blMan->BHD(i)->ctmean();
        MulBHD++;
      }
    }
  }
  if(MulBHD!=1){
    return true;
  }
  //if((T0Timing-BHDTiming)<25||27<(T0Timing-BHDTiming)){ /* Pion */
  //  return true;
  //}
  //if((T0Timing-BHDTiming)<27.8||29.5<(T0Timing-BHDTiming)){ /* Kaon */
  //  return true;
  //}

  /* Trigger pattern */
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }

  /* Basic data fill */
  const int nhodo = 7;
  int hodoid[nhodo]={
    CID_BHD,CID_T0,CID_BPD,CID_DEF,CID_BVC,
    CID_HVC1,CID_HVC2
  };
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    int nhit=0;
    for( int i=0; i<blMan->nHodo(cid); i++ ){
      HodoscopeLikeHit *hit = blMan->Hodoi(cid,i);
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
      double udtof[2]; udtof[0]=ctime[0]-T0Timing; udtof[1]=ctime[1]-T0Timing;
      double tof = ctmean - T0Timing;
      if(cid==CID_BHD){
        tof *= -1.0;
        udtof[0] *= -1.0;
        udtof[1] *= -1.0;
      }

      const char* ud[2] = {"u","d"};
      for(int iud=0; iud<2; iud++){
        FillHist(Form("%s%s%d_ADC",name,ud[iud],seg),adc[iud],1);
        if(tdc[iud]>0){
          FillHist(Form("%s%s%d_TDC",name,ud[iud],seg),tdc[iud],1);
          FillHist(Form("%s%s%d_Time",name,ud[iud],seg),time[iud],1);
          FillHist(Form("%s%s%d_CTime",name,ud[iud],seg),ctime[iud],1);
        }
      }
      if( hit->CheckRange() ){
        for(int iud=0; iud<2; iud++){
          if(adc[iud]>4000){ continue; }
          FillHist(Form("%s%s%d_ADCwithTDC",name,ud[iud],seg),adc[iud],1);
          FillHist(Form("%s%s%d_dE",name,ud[iud],seg),ene[iud],1);
          FillHist(Form("%s%s%d_ADCvsTDC",name,ud[iud],seg),adc[iud],tdc[iud],1);
          FillHist(Form("%s%s%d_dEvsTime",name,ud[iud],seg),ene[iud],time[iud],1);
          FillHist(Form("%s%s%d_dEvsCTime",name,ud[iud],seg),ene[iud],ctime[iud],1);
          FillHist(Form("%s%s%d_dEMeanvsCTime",name,ud[iud],seg),emean,ctime[iud],1);
          if(cid!=CID_T0){
            FillHist(Form("%s%s%d_dEMeanvsTOF",name,ud[iud],seg),emean,udtof[iud],1);
            FillHist(Form("%s%s%d_dEvsTOF",name,ud[iud],seg),ene[iud],tof,1);
          }
          FillHist(Form("%s%s%d_dEvsCTimeMean",name,ud[iud],seg),ene[iud],ctmean,1);
        }
        FillHist(Form("%s%d_dEMean",name,seg),emean,1);
        FillHist(Form("%s%d_TimeMean",name,seg),tmean,1);
        FillHist(Form("%s%d_TimeSub",name,seg),tsub,1);
        FillHist(Form("%s%d_CTimeMean",name,seg),ctmean,1);
        FillHist(Form("%s_CTimeMean",name),ctmean,1);
        FillHist(Form("%s%d_EventvsCTimeMean",name,seg),header->ev(),ctmean,1);
        FillHist(Form("%s_EventvsCTimeMean",name),header->ev(),ctmean,1);
        FillHist(Form("%s%d_CTimeSub",name,seg),ctsub,1);
        if(cid!=CID_T0){
          FillHist(Form("%s%d_TOF",name,seg),tof,1);
          FillHist(Form("%s_TOF",name),tof,1);
        }
        FillHist(Form("%s_TOFvsT0Segment",name),tof,T0Segment,1);
        FillHist(Form("%s_TOFvsCTimeMean",name),tof,ctmean,1);
        FillHist(Form("%s_TOFvsT0CTimeMean",name),tof,T0Timing,1);
        FillHist(Form("%s%d_Position",name,seg),hitpos,1);
        FillHist(Form("%s%d_dEMeanvsCTimeMean",name,seg),emean,ctmean,1);
        FillHist(Form("%s_HitPattern",name),seg,1);
        nhit++;
        if(cid==CID_BHD){
          FillHist(Form("T0u%d_dEvsTOF",T0Segment),T0udE,tof,1);
          FillHist(Form("T0d%d_dEvsTOF",T0Segment),T0ddE,tof,1);
          FillHist(Form("T0%d_dEMeanvsTOF",T0Segment),T0mdE,tof,1);
          FillHist(Form("T0%d_TOF",T0Segment),tof,1);
          FillHist(Form("T0_TOF"),tof,1);
        }
        for( int j=0; j<20; j++ ){
          int val = header->pattern(j);
          if( 0<val ){
            FillHist(Form("%s_TriggerPatternvsTOF",name),j,tof,1);
          }
        }
      }
    }
    FillHist(Form("%s_Multiplicity",name),nhit,1);
  }


}

bool MyAnalysisBLHodoCheck::FillHist(TString name, double val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1&&weight>0){ 
    for(int i=0; i<weight; i++){
      h1 -> Fill(val1,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBLHodoCheck::FillHist(TString name, TString val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1&&weight>0){
    for(int i=0; i<weight; i++){
      h1 -> Fill(val1,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBLHodoCheck::FillHist(TString name, double val1, double val2, int weight)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2&&weight>0){
    for(int i=0; i<weight; i++){
      h2 -> Fill(val1,val2,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBLHodoCheck::FillHist(TString name, TString val1, TString val2, int weight)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2&&weight>0){
    for(int i=0; i<weight; i++){
      h2 -> Fill(val1,val2,1);
    }
    return true;
  }
  else {
    return false;
  }
}

void MyAnalysisBLHodoCheck::CutCondition()
{
}

double MyAnalysisBLHodoCheck::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

double MyAnalysisBLHodoCheck::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

bool MyAnalysisBLHodoCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisBLHodoCheck::Initialize ###" << std::endl;

  confFile = confMan;
  DetectorList *dlist=DetectorList::GetInstance();

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaBLHodoCheck");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();
  TH1F* h1;

  /* Trigger pattern */
  h1 = new TH1F( "TriggerPattern", "Trigger Pattern", 18, -0.5, 17.5);
  h1->GetXaxis()->SetBinLabel(1,"All");
  h1->GetXaxis()->SetBinLabel(2,"Beam");
  h1->GetXaxis()->SetBinLabel(3,"Kaon");
  h1->GetXaxis()->SetBinLabel(4,"KCDH1f");
  h1->GetXaxis()->SetBinLabel(5,"Pion");
  h1->GetXaxis()->SetBinLabel(6,"Proton");
  h1->GetXaxis()->SetBinLabel(7,"KCDH1");
  h1->GetXaxis()->SetBinLabel(8,"KCDH2");
  h1->GetXaxis()->SetBinLabel(9,"PivBVC");
  h1->GetXaxis()->SetBinLabel(10,"PiCDH1");
  h1->GetXaxis()->SetBinLabel(11,"PiCDH2");
  h1->GetXaxis()->SetBinLabel(12,"Kf");
  h1->GetXaxis()->SetBinLabel(13,"1stMix");
  h1->GetXaxis()->SetBinLabel(14,"Charged");
  h1->GetXaxis()->SetBinLabel(15,"Neutral");
  h1->GetXaxis()->SetBinLabel(16,"Cosmic");
  h1->GetXaxis()->SetBinLabel(17,"Reject");
  h1->GetXaxis()->SetBinLabel(18,"SIM");

  /* Hodoscope */
  std::cout << "Define Histograms for Hodoscopes" << std::endl;
  const int nhodo = 7;
  int hodoid[nhodo]={
    CID_BHD,CID_T0,CID_BPD,CID_DEF,CID_BVC,
    CID_HVC1,CID_HVC2
  };
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    for( int seg=1; seg<=nsegs; seg++ ){
      const char* ud[2] = {"u","d"};
      const char* udname[2] = {"up","down"};
      for(int iud=0; iud<2; iud++){
        new TH1F( Form("%s%s%d_ADC",name,ud[iud],seg),   Form("ADC %s-%d(%s);ADC ch.;Counts",name,seg,udname[iud]),    4000,    0, 4000 );
        new TH1F( Form("%s%s%d_ADCwithTDC",name,ud[iud],seg), Form("ADC with TDC %s-%d(%s);ADC ch.;Counts",name,seg,udname[iud]),  4000,    0, 4000 );
        new TH1F( Form("%s%s%d_dE",name,ud[iud],seg),   Form("dE %s-%d(%s);dE (MeV);Counts",name,seg,udname[iud]),    1000,    0, 40 );
        new TH1F( Form("%s%s%d_TDC",name,ud[iud],seg),   Form("TDC %s-%d(%s);TDC ch.;Counts",name,seg,udname[iud]),    2000,    0, 4000 );
        new TH1F( Form("%s%s%d_Time",name,ud[iud],seg),   Form("Time %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -10, 10 );
        new TH1F( Form("%s%s%d_CTime",name,ud[iud],seg),   Form("CTime %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
        new TH2F( Form("%s%s%d_ADCvsTDC",name,ud[iud],seg),   Form("ADC and TDC corr. %s-%d(%s);ADC ch.;TDC ch.",name,seg,udname[iud]),     400,    0, 2000,  400,    0, 4000 );
        new TH2F( Form("%s%s%d_dEvsTime",name,ud[iud],seg),   Form("dE and Time corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEvsCTime",name,ud[iud],seg),   Form("dE and CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEvsTOF",name,ud[iud],seg),   Form("dE and TOF corr. %s-%d(%s);dE (MeV);TOF (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEvsCTimeMean",name,ud[iud],seg),   Form("dE and Mean CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
        //new TH2F( Form("%s%s%d_dEMeanvsTOF",name,ud[iud],seg),   Form("dEMean TOF corr. %s-%d(%s);dE (MeV);TOF (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
        //new TH2F( Form("%s%s%d_dEMeanvsCTime",name,ud[iud],seg),   Form("dEMean CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
      }
      new TH1F( Form("%s%d_dEMean",name,seg),   Form("Mean dE %s-%d;dE (MeV);Counts",name,seg),    1000,    0, 40 );
      new TH1F( Form("%s%d_TimeMean",name,seg),   Form("Mean Time %s-%d;Time (ns);Counts",name,seg),    2000,    -10, 10 );
      new TH1F( Form("%s%d_TimeSub",name,seg),   Form("Subtract Time %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_CTimeMean",name,seg),   Form("Mean CTime %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_CTimeSub",name,seg),   Form("Subtract CTime %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_TOF",name,seg),   Form("TOF between T0 and %s-%d;Time of flight T0-%s%d (ns);Counts",name,seg,name,seg),    6000,    -50, 100 );
      new TH1F( Form("%s%d_Position",name,seg),   Form("Hit position %s-%d;Position (cm);Counts",name,seg),    1000,    -50, 50 );
      new TH2F( Form("%s%d_dEMeanvsTOF",name,seg),   Form("dEMean TOF corr. %s-%d;dE (MeV);TOF (ns)",name,seg),     200,    0, 20,  2000,    -100, 100 );
      new TH2F( Form("%s%d_dEMeanvsCTimeMean",name,seg),   Form("dEMean CTimeMean corr. %s-%d;dE (MeV);Time (ns)",name,seg),     200,    0, 20,  2000,    -100, 100 );
      new TH2F( Form("%s%d_EventvsCTimeMean",name,seg),   Form("Event number vs. Mean CTime %s-%d;Event number;Time (ns);",name,seg), 800, 0, 800000,   2000,    -100, 100 );
    }
    new TH1F( Form("%s_HitPattern",name), Form("Hit pattern %s;%s Segment;Counts",name,name), nsegs, 0.5, nsegs+0.5);
    new TH1F( Form("%s_Multiplicity",name), Form("Multiplicity %s;Multiplicity;Counts",name), nsegs+1, -0.5, nsegs+0.5 );
    new TH1F( Form("%s_CTimeMean",name),   Form("Mean CTime %s;Time (ns);Counts",name),    2000,    -100, 100 );
    new TH2F( Form("%s_EventvsCTimeMean",name),   Form("Event number vs. Mean CTime %s;Event number;Time (ns);",name), 800, 0, 800000,   2000,    -100, 100 );
    new TH1F( Form("%s_TDC",name),   Form("TDC %s;TDC ch.;Counts",name),    400,    0, 4000 );
    new TH1F( Form("%s_TOF",name),   Form("TOF between T0 and %s;Time of flight T0-%s (ns);Counts",name,name),    6000,    -50, 100 );
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

  return true;
}
