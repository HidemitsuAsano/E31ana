// MyAnalysisCVCPCCheck.cpp

#include "MyAnalysisCVCPCCheck.h"

MyAnalysisCVCPCCheck::MyAnalysisCVCPCCheck(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisCVCPCCheck::~MyAnalysisCVCPCCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCVCPCCheck::Clear()
{
}

bool MyAnalysisCVCPCCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{

  rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();

  /* Trigger selection */
  if(!header->IsTrig(Trig_Pion)){
    return true;
  }
  //if(!header->IsTrig(Trig_Kaon)){
  //  return true;
  //}

  /* Event selection */
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
  if(T0Segment!=3){
    return true;
  }
  double BHDTiming = -9999.; 
  int BHDSegment = -1;
  int MulBHD = 0;
  for(int i=0; i<blMan->nBHD(); i++){
    if(blMan->BHD(i)->CheckRange()){
      MulBHD++;
      BHDTiming = blMan->BHD(i)->ctmean();
      BHDSegment = blMan->BHD(i)->seg();
    }
  }
  if(MulBHD!=1){
    return true;
  }
  //if((T0Timing-BHDTiming)<27.5||30<(T0Timing-BHDTiming)){
  //  return true;
  //}
  if((T0Timing-BHDTiming)<25||27<(T0Timing-BHDTiming)){
    return true;
  }

  int MulBVC = 0;
  for(int i=0; i<blMan->nBVC(); i++){
    if(blMan->BVC(i)->CheckRange()){
      MulBVC++;
    }
  }
  int MulPC = 0;
  for(int i=0; i<blMan->nPC(); i++){
    if(blMan->PC(i)->CheckRange()){
      MulPC++;
    }
  }
  int MulCVC = 0;
  for(int i=0; i<blMan->nCVC(); i++){
    if(blMan->CVC(i)->CheckRange()){
      MulCVC++;
    }
  }
  if(MulBVC==0||(MulPC==0&&MulCVC==0)){
    return true;
  }

  /* Trigger pattern */
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }

  // ############# //
  // FDC Selection //
  // ############# //
  bltrackMan->LocalTracking( blMan, conf, CID_FDC1, "noslopenotiming");
  FillHist("FDC1_nTrack",bltrackMan->ntrackBLDC(CID_FDC1));
  if(bltrackMan->ntrackBLDC(CID_FDC1)!=1) return true;
  LocalTrack* fdc = bltrackMan->trackBLDC(CID_FDC1,0);
  FillHist("FDC1_Chi2",fdc->chi2all());
  if(fdc->chi2all()>30) return true;
  TVector3 fdcpos = fdc->GetPosatZ(169.6);
  FillHist("FDC1_Position",fdcpos.X(),fdcpos.Y());

  // ############ //
  // CVC Analysis //
  // ############ //
  int nCVC = 0;
  for(int i=0; i<blMan->nCVC(); i++){
    HodoscopeLikeHit* hit = blMan->CVC(i);
    int seg = hit->seg();
    double emean = hit->emean();
    double ctmean = hit->ctmean();
    if(hit->CheckRange()){
      FillHist(Form("CVC%d_CTimeMean",seg),ctmean,1);
      FillHist(Form("CVC%d_dEMean",seg),emean,1);
      FillHist(Form("CVC_HitPattern"),seg,1);
      FillHist(Form("FWD_HitPattern"),seg,1);
      nCVC++;
    }
  }
  FillHist(Form("CVC_Multiplicity"),nCVC,1);

  // ############ //
  // PC Analysis //
  // ############ //
  int nPC = 0;
  for(int i=0; i<blMan->nPC(); i++){
    HodoscopeLikeHit* hit = blMan->PC(i);
    int seg = hit->seg();
    double emean = hit->emean();
    double ctmean = hit->ctmean();
    if(hit->CheckRange()){
      FillHist(Form("PC%d_CTimeMean",seg),ctmean,1);
      FillHist(Form("PC%d_dEMean",seg),emean,1);
      FillHist(Form("PC_HitPattern"),seg,1);
      FillHist(Form("FWD_HitPattern"),seg+34,1);
      nPC++;
    }
  }
  FillHist(Form("PC_Multiplicity"),nPC,1);

  FillHist(Form("FWD_Multiplicity"),nCVC+nPC,1);

  return true;
}

bool MyAnalysisCVCPCCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCVCPCCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCVCPCCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCVCPCCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisCVCPCCheck::CutCondition()
{
}

double MyAnalysisCVCPCCheck::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

double MyAnalysisCVCPCCheck::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

bool MyAnalysisCVCPCCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisCVCPCCheck::Initialize ###" << std::endl;

  confFile = confMan;
  DetectorList *dlist=DetectorList::GetInstance();

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaCVCPCCheck");

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

  // FDC Tracking
  new TH1F( "FDC1_nTrack", "Number of Tracks in FDC1;Num. of Tracks;Counts", 20, 0, 20 );
  new TH1F( "FDC1_Chi2", "Chi-squared of FDC1;Chi-square;Counts", 500, 0, 50);
  new TH2F( "FDC1_Position", "Position at FDC1;X (cm);Y (cm)", 1000,-50,50, 1000,-50,50);
  // PC/CVC
  new TH1F( Form("FWD_HitPattern"), Form("Hit pattern FWD;FWD Segment;Counts"), 60, 0.5, 60+0.5);
    new TH1F( Form("FWD_Multiplicity"), Form("Multiplicity FWD;Multiplicity;Counts"), 30+1, -0.5, 30+0.5 );
  const int nhodo = 2;
  int hodoid[nhodo]={CID_CVC,CID_PC};
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    new TH1F( Form("%s_HitPattern",name), Form("Hit pattern %s;%s Segment;Counts",name,name), nsegs, 0.5, nsegs+0.5);
    new TH1F( Form("%s_Multiplicity",name), Form("Multiplicity %s;Multiplicity;Counts",name), nsegs+1, -0.5, nsegs+0.5 );
    for(int seg=1; seg<=nsegs; seg++){
      new TH1F( Form("%s%d_CTimeMean",name,seg),   Form("Mean CTime %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_dEMean",name,seg),   Form("Mean dE %s-%d;dE (MeV);Counts",name,seg),    4000,    0, 150 );
    }
  }

  return true;
}
