// MyAnalysisCDSCheck.cpp

#include "MyAnalysisCDSCheck.h"

MyAnalysisCDSCheck::MyAnalysisCDSCheck(TFile* rt, ConfMan* conf) {
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisCDSCheck::~MyAnalysisCDSCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCDSCheck::Clear()
{
}

bool MyAnalysisCDSCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
{

  rtFile->cd();
	DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  //if(!header->IsTrig(Trig_Pion)){
  //  return true;
  //}
  //if(!header->IsTrig(Trig_Kaon)){
  //  return true;
  //}
  if(header->IsTrig(Trig_Cosmic)){
    return true;
  }
  /* ========== Trigger selection ========== */
  /* ======================================= */

  /* ===================================== */
  /* ========== Event selection ========== */
  double T0Timing = -9999.; 
  int T0Segment = -1;
  int MulT0 = 0;
  for(int i=0; i<blMan->nT0(); i++){
    if(blMan->T0(i)->CheckRange()){
      MulT0++;
      T0Timing = blMan->T0(i)->ctmean();
      T0Segment = blMan->T0(i)->seg();
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
  if((T0Timing-BHDTiming)<27.8||29.5<(T0Timing-BHDTiming)){ /* Kaon */
    return true;
  }
  //if((T0Timing-BHDTiming)<34||36<(T0Timing-BHDTiming)){ /* Proton */
  //  return true;
  //}
  int MulCDH = 0;
  int CDHSegment[36];
  int CDHEvNum[36];
  double CDHTiming[36];
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()){
      CDHSegment[MulCDH] = hit->seg();
      CDHTiming[MulCDH] = hit->ctmean();
      CDHEvNum[hit->seg()-1] = i;
      MulCDH++;
    }
  }
  int MulIH = 0;
  int IHSegment[36];
  int IHEvNum[36];
  double IHTiming[36];
  for(int i=0; i<cdsMan->nIH(); i++){
    HodoscopeLikeHit* hit = cdsMan->IH(i);
    if(hit->CheckRange()){
      IHSegment[MulIH] = hit->seg();
      IHTiming[MulIH] = hit->ctmean();
      IHEvNum[hit->seg()-1] = i;
      MulIH++;
    }
  }
  /* ========== Event selection ========== */
  /* ===================================== */

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }

  /* Fill data */
  const int nhodo = 2;
  int hodoid[nhodo]={
    CID_CDH,CID_IH
  };
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    int nhit=0;
    for( int i=0; i<cdsMan->nHodo(cid); i++ ){
      HodoscopeLikeHit *hit = cdsMan->Hodoi(cid,i);
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
      if(cid==CID_IH){
        emean = ene[0];
        tmean = time[0];
        ctmean = ctime[0];
        tsub = time[0];
        ctsub = time[0];
      }
      double tof = ctmean - T0Timing;
      if(cid==CID_BHD){
        tof *= -1.0;
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
          FillHist(Form("%s%s%d_ADCwithTDC",name,ud[iud],seg),adc[iud],1);
          FillHist(Form("%s%s%d_dE",name,ud[iud],seg),ene[iud],1);
          FillHist(Form("%s%s%d_ADCvsTDC",name,ud[iud],seg),adc[iud],tdc[iud],1);
          FillHist(Form("%s%s%d_dEvsTime",name,ud[iud],seg),ene[iud],time[iud],1);
          FillHist(Form("%s%s%d_dEvsCTime",name,ud[iud],seg),ene[iud],ctime[iud],1);
          FillHist(Form("%s%s%d_dEMeanvsCTime",name,ud[iud],seg),emean,ctime[iud],1);
        }
        FillHist(Form("%s%d_dEMean",name,seg),emean,1);
        FillHist(Form("%s%d_TimeMean",name,seg),tmean,1);
        FillHist(Form("%s%d_TimeSub",name,seg),tsub,1);
        FillHist(Form("%s%d_CTimeMean",name,seg),ctmean,1);
        FillHist(Form("%s_CTimeMean",name),ctmean,1);
        FillHist(Form("%s%d_EventvsCTimeMean",name,seg),header->ev(),ctmean,1);
        FillHist(Form("%s_EventvsCTimeMean",name),header->ev(),ctmean,1);
        FillHist(Form("%s%d_CTimeSub",name,seg),ctsub,1);
        FillHist(Form("%s%d_TOF",name,seg),tof,1);
        FillHist(Form("%s_TOF",name),tof,1);
        FillHist(Form("%s_TOFvsT0Segment",name),tof,T0Segment,1);
        FillHist(Form("%s_TOFvsCTimeMean",name),tof,ctmean,1);
        FillHist(Form("%s_TOFvsT0CTimeMean",name),tof,T0Timing,1);
        FillHist(Form("%s%d_Position",name,seg),hitpos,1);
        FillHist(Form("%s%d_dEMeanvsCTimeMean",name,seg),emean,ctmean,1);
        FillHist(Form("%s_HitPattern",name),seg,1);
        nhit++;
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

  // Fill CDC
  const int cid=CID_CDC;
  const char* name="CDC";
  const int nlays= dlist->GetNlayers(cid);
  for(int layer=1;layer<=nlays;layer++){
    int mult = cdsMan->nCDC(layer);
    FillHist(Form("%sl%d_Multiplicity",name,layer),mult,1);
    for(int i=0;i<mult;i++){
      CDCHit *hit = cdsMan->CDC(layer,i);
      int wire = hit->wire();
      int tdc = hit->tdc();
      double dt = hit->dt();
      double dl = hit->dl()*100;
      FillHist(Form("%sl%d_HitPattern",name,layer),wire,1);
      FillHist(Form("%sl%d_TDC",name,layer),tdc,1);
      FillHist(Form("%sl%d_Time",name,layer),dt,1);
      FillHist(Form("%sl%d_Length",name,layer),dl,1);
      FillHist(Form("%sl%d_TimevsLength",name,layer),dt,dl,1);
      FillHist(Form("%sl%dw%02d_TDC",name,layer,wire),tdc,1);
      FillHist(Form("%sl%dw%02d_Time",name,layer,wire),dt,1);
      FillHist(Form("%sl%dw%02d_Length",name,layer,wire),dl,1);
    }
  }

  if(cdstrackMan==0) return true;

  /* For CDC Tracking */
  int ntrack = cdstrackMan->nTrack();
  int ngoodtrack = cdstrackMan->nGoodTrack();
  FillHist(Form("%s_NumberofTrack",name),ntrack,1);
  FillHist(Form("%s_NumberofGoodTrack",name),ngoodtrack,1);
  for(int it=0; it<cdstrackMan->nGoodTrack(); it++){
    CDSTrack* track = cdstrackMan->Track(cdstrackMan->GoodTrackID(it));
    double mom = track->Momentum();
    double chi =  track->Chi();
    FillHist(Form("%s_Momentum",name),mom,1);
    FillHist(Form("%s_ChiSquare",name),chi,1);
    for(int layer=1;layer<=nlays;layer++){
      int mult = track->nTrackHit(layer);
      FillHist(Form("%sl%d_MultiplicityinTrack",name,layer),mult,1);
      for(int i=0;i<mult;i++){
        CDCHit *hit = track->TrackHit(cdsMan,layer,i);
        int wire = hit->wire();
        int tdc = hit->tdc();
        double dt = hit->dt();
        double dl = hit->dl()*100;
        double resl = hit->resl();
        FillHist(Form("%sl%d_HitPatterninTrack",name,layer),wire,1);
        FillHist(Form("%sl%d_TimeinTrack",name,layer),dt,1);
        FillHist(Form("%sl%d_LengthinTrack",name,layer),dl,1);
        FillHist(Form("%sl%d_ResidualinTrack",name,layer),resl,1);
        FillHist(Form("%sl%d_TimevsResidualinTrack",name,layer),dt,resl,1);
        FillHist(Form("%sl%d_TimevsLengthinTrack",name,layer),dt,dl,1);
        FillHist(Form("%sl%dw%02d_TimeinTrack",name,layer,wire),dt,1);
        FillHist(Form("%sl%dw%02d_ResidualinTrack",name,layer,wire),resl,1);
        FillHist(Form("%sl%dw%02d_TimevsResidualinTrack",name,layer,wire),dt,resl,1);
      }
    }
  }


}

bool MyAnalysisCDSCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCDSCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCDSCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCDSCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisCDSCheck::CutCondition()
{
}

double MyAnalysisCDSCheck::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  double tmpoffset = -9999.9;
  if(!hit->CheckRange()) return tmpoffset;
  pBeam* beam = particle->beam(0);
  double tof = hit->ctmean() - beam->t0time();

  TVector3 bpdpos; double bpdposX, bpdposY;
  confFile->GetGeomMapManager()->GetGPos(CID_BPD,0,bpdpos);
  beam -> bpcpos(bpdpos.Z(),bpdposX,bpdposY);
  bpdpos.SetX(bpdposX); bpdpos.SetY(bpdposY);
  TVector3 t0pos; double t0posX, t0posY;
  confFile->GetGeomMapManager()->GetGPos(CID_T0,0,t0pos);
  beam -> blcpos(t0pos.Z(),t0posX,t0posY);
  t0pos.SetX(t0posX); t0pos.SetY(t0posY);
  double calctof = (bpdpos-t0pos).Mag() / beam->beta() * Const;
  tmpoffset = tof - calctof;
  return tmpoffset;
}

double MyAnalysisCDSCheck::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  double tmpoffset = -9999.9;
  if(!hit->CheckRange()) return tmpoffset;
  pBeam* beam = particle->beam(0);
  double ctsub = hit->ctsub();

  TVector3 bpcpos; double bpcposX, bpcposY;
  confFile->GetGeomMapManager()->GetGPos(CID_BPD,0,bpcpos);
  beam -> bpcpos(bpcpos.Z(),bpcposX,bpcposY);
  bpcpos.SetX(bpcposX); bpcpos.SetY(bpcposY);
  double lv;
  if(confFile->GetGeomMapManager()->GetLightVelocity(CID_BPD,hit->seg(),lv)){
    double calctsub = -bpcpos.Y() / lv;
    tmpoffset = ctsub - calctsub;
  }
  return tmpoffset;
}

bool MyAnalysisCDSCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisCDSCheck::Initialize ###" << std::endl;

  confFile = confMan;
  DetectorList *dlist=DetectorList::GetInstance();

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaCDSCheck");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  /* Trigger pattern */
  TH1F* h1 = new TH1F( "TriggerPattern", "Trigger Pattern", 18, -0.5, 17.5);
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
  //std::cout << "Define Histograms for Hodoscopes" << std::endl;
  //const int nhodo = 2;
  //int hodoid[nhodo]={
  //  CID_CDH,CID_IH
  //};
  //for(int ihodo=0; ihodo<2; ihodo++){
  //  const int cid = hodoid[ihodo];
  //  const char* name = dlist -> GetName(cid).c_str();
  //  const int nsegs = dlist -> GetNsegs(cid);
  //  for( int seg=1; seg<=nsegs; seg++ ){
  //    const char* ud[2] = {"u","d"};
  //    const char* udname[2] = {"up","down"};
  //    for(int iud=0; iud<2; iud++){
  //      new TH1F( Form("%s%s%d_ADC",name,ud[iud],seg),   Form("ADC %s-%d(%s);ADC ch.;Counts",name,seg,udname[iud]),    2000,    0, 4000 );
  //      new TH1F( Form("%s%s%d_ADCwithTDC",name,ud[iud],seg), Form("ADC with TDC %s-%d(%s);ADC ch.;Counts",name,seg,udname[iud]),  2000,    0, 4000 );
  //      new TH1F( Form("%s%s%d_dE",name,ud[iud],seg),   Form("dE %s-%d(%s);dE (MeV);Counts",name,seg,udname[iud]),    800,    0, 40 );
  //      new TH1F( Form("%s%s%d_TDC",name,ud[iud],seg),   Form("TDC %s-%d(%s);TDC ch.;Counts",name,seg,udname[iud]),    2000,    0, 4000 );
  //      new TH1F( Form("%s%s%d_Time",name,ud[iud],seg),   Form("Time %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
  //      new TH1F( Form("%s%s%d_CTime",name,ud[iud],seg),   Form("CTime %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
  //      new TH2F( Form("%s%s%d_ADCvsTDC",name,ud[iud],seg),   Form("ADC TDC corr. %s-%d(%s);ADC ch.;TDC ch.",name,seg,udname[iud]),     400,    0, 4000,  400,    0, 4000 );
  //      new TH2F( Form("%s%s%d_dEvsTime",name,ud[iud],seg),   Form("dE Time corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
  //      new TH2F( Form("%s%s%d_dEvsCTime",name,ud[iud],seg),   Form("dE CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
  //      new TH2F( Form("%s%s%d_dEMeanvsCTime",name,ud[iud],seg),   Form("dEMean CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
  //    }
  //    new TH1F( Form("%s%d_dEMean",name,seg),   Form("Mean dE %s-%d;dE (MeV);Counts",name,seg),    800,    0, 40 );
  //    new TH1F( Form("%s%d_TimeMean",name,seg),   Form("Mean Time %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
  //    new TH1F( Form("%s%d_TimeSub",name,seg),   Form("Subtract Time %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
  //    new TH1F( Form("%s%d_CTimeMean",name,seg),   Form("Mean CTime %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
  //    new TH1F( Form("%s%d_CTimeSub",name,seg),   Form("Subtract CTime %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
  //    new TH1F( Form("%s%d_TOF",name,seg),   Form("TOF between T0 and %s-%d;Time of flight T0-%s%d (ns);Counts",name,seg,name,seg),    3000,    -50, 100 );
  //    new TH1F( Form("%s%d_Position",name,seg),   Form("Hit position %s-%d;Position (cm);Counts",name,seg),    1000,    -50, 50 );
  //    new TH2F( Form("%s%d_dEMeanvsCTimeMean",name,seg),   Form("dEMean CTimeMean corr. %s-%d;dE (MeV);Time (ns)",name,seg),     200,    0, 20,  2000,    -100, 100 );
  //    new TH2F( Form("%s%d_EventvsCTimeMean",name,seg),   Form("Event number vs. Mean CTime %s-%d;Event number;Time (ns);",name,seg), 800, 0, 800000,   2000,    -100, 100 );
  //  }
  //  new TH1F( Form("%s_HitPattern",name), Form("Hit pattern %s;%s Segment;Counts",name,name), nsegs, 0.5, nsegs+0.5);
  //  new TH1F( Form("%s_Multiplicity",name), Form("Multiplicity %s;Multiplicity;Counts",name), nsegs+1, -0.5, nsegs+0.5 );
  //  new TH1F( Form("%s_CTimeMean",name),   Form("Mean CTime %s;Time (ns);Counts",name),    2000,    -100, 100 );
  //  new TH2F( Form("%s_EventvsCTimeMean",name),   Form("Event number vs. Mean CTime %s;Event number;Time (ns);",name), 800, 0, 800000,   2000,    -100, 100 );
  //  new TH1F( Form("%s_TOF",name),   Form("TOF between T0 and %s;Time of flight T0-%s (ns);Counts",name,name),    3000,    -50, 100 );
  //  new TH2F( Form("%s_TOFvsT0Segment",name),   Form("TOF between T0 and %s vs. T0 segment;Time of flight T0-%s (ns);T0 segment",name,name),    3000,    -50, 100, 5, 0.5, 5.5 );
  //  new TH2F( Form("%s_TOFvsCTimeMean",name),   Form("TOF between T0 and %s vs. CTimeMean;Time of flight T0-%s (ns);Time (ns)",name,name),    3000,    -50, 100, 2000, -100, 100 );
  //  new TH2F( Form("%s_TOFvsT0CTimeMean",name),   Form("TOF between T0 and %s vs. T0 CTimeMean;Time of flight T0-%s (ns);T0 Time (ns)",name,name),    3000,    -50, 100, 200, -10, 10 );
  //  TH2F* h2 = new TH2F( Form("%s_TriggerPatternvsTOF",name),   Form("TriggerPattern vs. TOF between T0 and %s;;Time of flight T0-%s (ns)",name,name), 18, -0.5, 17.5 ,    3000,    -50, 100);
  //  h2->GetXaxis()->SetBinLabel(1,"All");
  //  h2->GetXaxis()->SetBinLabel(2,"Beam");
  //  h2->GetXaxis()->SetBinLabel(3,"Kaon");
  //  h2->GetXaxis()->SetBinLabel(4,"KCDH1f");
  //  h2->GetXaxis()->SetBinLabel(5,"Pion");
  //  h2->GetXaxis()->SetBinLabel(6,"Proton");
  //  h2->GetXaxis()->SetBinLabel(7,"KCDH1");
  //  h2->GetXaxis()->SetBinLabel(8,"KCDH2");
  //  h2->GetXaxis()->SetBinLabel(9,"PivBVC");
  //  h2->GetXaxis()->SetBinLabel(10,"PiCDH1");
  //  h2->GetXaxis()->SetBinLabel(11,"PiCDH2");
  //  h2->GetXaxis()->SetBinLabel(12,"Kf");
  //  h2->GetXaxis()->SetBinLabel(13,"1stMix");
  //  h2->GetXaxis()->SetBinLabel(14,"Charged");
  //  h2->GetXaxis()->SetBinLabel(15,"Neutral");
  //  h2->GetXaxis()->SetBinLabel(16,"Cosmic");
  //  h2->GetXaxis()->SetBinLabel(17,"Reject");
  //  h2->GetXaxis()->SetBinLabel(18,"SIM");
  //}
  //return true;

  /* CDC */
  const int cid= CID_CDC;
  const char* name= dlist->GetName(cid).data();
  const int nlays= NumOfCDCLayers;
  std::cout << "Define Histograms for " << name << std::endl;
  new TH1F( Form("%s_ChiSquare",name), Form("ChiSquare %s;#Chi^{2}/NDF;Counts",name), 1000, 0, 50 );
  new TH1F( Form("%s_NumberofTrack",name), Form("Number of tracks %s;Number of Tracks;Counts",name), 20, 0, 20 );
  new TH1F( Form("%s_NumberofGoodTrack",name), Form("Number of good tracks %s;Number of Good tracks;Counts",name), 20, 0, 20 );
  new TH1F( Form("%s_Momentum",name), Form("Momentum %s;Momentum (GeV/c);Counts",name), 500, 0, 5 );
  for( int layer=1; layer<=nlays; layer++ ){
    int nwire = NumOfCDCWiresInLayer[layer-1];
    new TH1F( Form("%sl%d_Multiplicity",name,layer),Form("Multiplicity %s-L%d;%s(layer %d) Multiplicity;Counts",name,layer,name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_HitPattern",name,layer),Form("Hit pattern %s-L%d;%s (layer%d) Wire;Counts",name,layer,name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_TDC",name,layer),Form("TDC %s-L%d;%s (layer %d) TDC ch.;Counts",name,layer,name,layer),2000,0,2000);
    new TH1F( Form("%sl%d_Time",name,layer),Form("Time %s-L%d;%s (layer %d) Time (ns);Counts",name,layer,name,layer),800,-400,400);
    new TH1F( Form("%sl%d_Length",name,layer),Form("Length %s-L%d;%s (layer %d) Length (um);Counts",name,layer,name,layer),200,0,200);
    new TH2F( Form("%sl%d_TimevsLength",name,layer),Form("Time vs. Length %s-L%d;%s (layer %d) Time (ns);%s (layer %d) Length (um)",name,layer,name,layer,name,layer),320,-20,300,1000,0,100);

    new TH1F( Form("%sl%d_MultiplicityinTrack",name,layer),Form("Multiplicity %s-L%d;%s (layer %d) Multiplicity;Counts",name,layer,name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_HitPatterninTrack",name,layer),Form("Hit pattern %s-L%d;%s (layer %d) Wire;Counts",name,layer,name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_TimeinTrack",name,layer),Form("Time %s-L%d;%s (layer %d) Time (ns);Counts",name,layer,name,layer),800,-400,400);
    new TH1F( Form("%sl%d_ResidualinTrack",name,layer),Form("Residual %s-L%d;%s (layer %d) Residual (cm);Counts",name,layer,name,layer),2000,-0.2,0.2);
    new TH2F( Form("%sl%d_TimevsResidualinTrack",name,layer),Form("Time vs. Residual %s-L%d;%s (layer %d) Time (ns);Residual (cm)",name,layer,name,layer),320,-20,300,1000,-0.1,0.1);
    new TH2F( Form("%sl%d_TimevsLengthinTrack",name,layer),Form("Time vs. Length %s-L%d;%s (layer %d) Time (ns);%s (layer %d) Length (um)",name,layer,name,layer,name,layer),320,-20,300,1000,0,100);

    for(int wire=1;wire<=nwire;wire++){
      //new TH1F( Form("%sl%dw%02d_TDC",name,layer,wire),Form("TDC %s-L%dW%02d;%s (layer %d,wire %d) TDC ch.;Counts",name,layer,wire,name,layer,wire),2000,0,2000);
      new TH1F( Form("%sl%dw%02d_Time",name,layer,wire),Form("Time %s-L%dW%02d;%s (layer %d,wire %d) Time (ns);Counts",name,layer,wire,name,layer,wire),800,-400,400);
      //new TH1F( Form("%sl%dw%02d_Length",name,layer,wire),Form("Length %s-L%dW%02d;%s (layer %d,wire %d) Length (um);Counts",name,layer,wire,name,layer,wire),200,0,200);

      new TH1F( Form("%sl%dw%02d_TimeinTrack",name,layer,wire),Form("Time %s-L%dW%02d;%s (layer %d,wire %d) Time (ns);Counts",name,layer,wire,name,layer,wire),320,-20,300);
      //new TH1F( Form("%sl%dw%02d_ResidualinTrack",name,layer,wire),Form("Residual %s-L%dW%02d;%s (layer %d,wire %d) Residual (cm);Counts",name,layer,wire,name,layer,wire),1000,-0.1,0.1);
      new TH2F( Form("%sl%dw%02d_TimevsResidualinTrack",name,layer,wire),Form("Time vs. Residual %s-L%dW%02d;%s (layer %d,wire %d) Time (ns);%s (layer %d,wire %d) Residual (cm)",name,layer,wire,name,layer,wire,name,layer,wire),320,-20,300,1000,-0.1,0.1);
    }
  }

  return true;
}
