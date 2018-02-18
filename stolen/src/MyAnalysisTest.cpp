// MyAnalysisTest.cpp

#include "MyAnalysisTest.h"

MyAnalysisTest::MyAnalysisTest(TFile* rt, ConfMan* conf)
{
	Initialize(conf);
	CutCondition();
	Clear();
}

MyAnalysisTest::~MyAnalysisTest()
{
	Clear();
	rtFile->cd();
	rtFile->Write();
	rtFile->Close();
}

void MyAnalysisTest::Clear()
{
}

bool MyAnalysisTest::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{

	rtFile->cd();
	DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  if(header->IsTrig(Trig_Cosmic)){
    return true;
  }
  /* ========== Trigger selection ========== */
  /* ======================================= */

  /* Event selection */
  double T0Timing = 9999.; 
  int T0Segment = -1;
  int MulT0 = 0;
  for(int i=0; i<blMan->nT0(); i++){
    if(blMan->T0(i)->CheckRange()){
      MulT0++;
      if(T0Timing>blMan->T0(i)->ctmean()){
        T0Timing = blMan->T0(i)->ctmean();
        T0Segment = blMan->T0(i)->seg();
      }
    }
  }
  double DEFTiming = -9999.; 
  int DEFSegment = -1;
  int MulDEF = 0;
  for(int i=0; i<blMan->nDEF(); i++){
    if(blMan->DEF(i)->CheckRange()){
      MulDEF++;
      DEFTiming = blMan->DEF(i)->ctmean();
      DEFSegment = blMan->DEF(i)->seg();
    }
  }
  int MulAC = 0;
  for(int i=0; i<blMan->nChere(CID_AC);i++){
    CherenkovLikeHit* hit = blMan->Cherei(CID_AC,i);
    if(hit->seg()!=1) continue;
    if(hit->CheckRange(1)){
      MulAC++;
    }
  }
  /* T0>0 && DEF>0 && vAC */
  if(MulT0==0){
    return true;
  }
  if(MulDEF==0){
    return true;
  }
  if(MulAC!=0){
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
  FillHist("BHD_Count",0,1);

  double BHDTiming = -9999.; 
  int BHDSegment = -1;
  int MulBHD = 0;
  double tof = -9999;
  for(int i=0; i<blMan->nBHD(); i++){
    if(blMan->BHD(i)->CheckRange()){
      MulBHD++;
      BHDTiming = blMan->BHD(i)->ctmean();
      BHDSegment = blMan->BHD(i)->seg();
      tof = -1.0*(BHDTiming - T0Timing);
    }
    if(24.8<tof&&tof<29.3){
      FillHist("BHD_Count",1,1);
      FillHist("BHD_HitPattern",BHDSegment,1);
      break;
    }
  }
}

bool MyAnalysisTest::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisTest::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisTest::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisTest::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisTest::CutCondition()
{
}

double MyAnalysisTest::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

double MyAnalysisTest::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

bool MyAnalysisTest::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisTest::Initialize ###" << std::endl;

  confFile = confMan;
  DetectorList *dlist=DetectorList::GetInstance();

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaTest");

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
  h1 = new TH1F( "BHD_Count", "BHD Counts;;Counts", 2, -0.5, 1.5);
  h1->GetXaxis()->SetBinLabel(1,"All");
  h1->GetXaxis()->SetBinLabel(2,"BHD Hit");

  return true;
}
