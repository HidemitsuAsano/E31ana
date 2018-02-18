// MyAnalysisTrigCheck.cpp

#include "MyAnalysisTrigCheck.h"

MyAnalysisTrigCheck::MyAnalysisTrigCheck(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisTrigCheck::~MyAnalysisTrigCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisTrigCheck::Clear()
{
}

bool MyAnalysisTrigCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, CDSHitMan* cdsMan)
{

  rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  //if(!header->IsTrig(Trig_Kf)||!header->IsTrig(Trig_1stMix)) return false;
  if(!header->IsTrig(Trig_KCDH1f)||!header->IsTrig(Trig_1stMix)) return false;
  /* ========== Trigger selection ========== */
  /* ======================================= */

  /* Trigger pattern */
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
      FillHist(Form("Trigger%02d_TDC",i),val);
    }
  }

  /* Counter Check */
  int MulT0 = 0;
  for(int i=0; i<blMan->nT0(); i++){
    if(blMan->T0(i)->CheckRange()){
      MulT0++;
    }
  }
  int MulBHD = 0;
  for(int i=0; i<blMan->nBHD(); i++){
    if(blMan->BHD(i)->CheckRange()){
      MulBHD++;
    }
  }
  int MulBVC = 0;
  for(int i=0; i<blMan->nBVC(); i++){
    if(blMan->BVC(i)->CheckRange()){
      MulBVC++;
    }
  }
  int MulCVC = 0;
  for(int i=0; i<blMan->nCVC(); i++){
    if(blMan->CVC(i)->CheckRange()){
      MulCVC++;
    }
  }
  int MulCVC2 = 0;
  for(int i=0; i<blMan->nCVC(); i++){
    if(blMan->CVC(i)->CheckRange()){
      if(blMan->CVC(i)->seg()>17)
        MulCVC2++;
    }
  }
  int MulPC = 0;
  for(int i=0; i<blMan->nPC(); i++){
    if(blMan->PC(i)->CheckRange()){
      MulPC++;
    }
  }
  int MulNC = 0;
  for(int i=0; i<blMan->nNC(); i++){
    if(blMan->NC(i)->CheckRange()){
      MulNC++;
    }
  }
  int MulCDH = 0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    if(cdsMan->CDH(i)->CheckRange()){
      MulCDH++;
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

  /* =================================== */
  /* ========== Trigger Check ========== */
  /* Beam */
  if(MulT0>0&&MulBHD>0){
    FillHist("TriggerEfficiency_Beam",0);
    if(header->IsTrig(Trig_Beam)){
      FillHist("TriggerEfficiency_Beam",1);
    }
  }
  /* Kaon */
  if(MulT0>0&&MulBHD>0&&MulAC==0){
    FillHist("TriggerEfficiency_Kaon",0);
    if(header->IsTrig(Trig_Kaon)){
      FillHist("TriggerEfficiency_Kaon",1);
    }
  }
  /* Kf */
  if(MulT0>0&&MulBHD>0&&MulAC==0){
    FillHist("TriggerEfficiency_Kf",0);
    if(header->IsTrig(Trig_Kf)){
      FillHist("TriggerEfficiency_Kf",1);
    }
  }
  /* Pion */
  if(MulT0>0&&MulBHD>0&&MulAC>0){
    FillHist("TriggerEfficiency_Pion",0);
    if(header->IsTrig(Trig_Pion)){
      FillHist("TriggerEfficiency_Pion",1);
    }
  }
  /* KCDH1 */
  if(MulT0>0&&MulBHD>0&&MulAC==0&&MulCDH>=1){
    FillHist("TriggerEfficiency_KCDH1",0);
    if(header->IsTrig(Trig_KCDH1)){
      FillHist("TriggerEfficiency_KCDH1",1);
    }
  }
  /* KCDH1f */
  if(MulT0>0&&MulBHD>0&&MulAC==0&&MulCDH>=1){
    FillHist("TriggerEfficiency_KCDH1f",0);
    if(header->IsTrig(Trig_KCDH1f)){
      FillHist("TriggerEfficiency_KCDH1f",1);
    }
  }
  /* KCDH2 */
  if(MulT0>0&&MulBHD>0&&MulAC==0&&MulCDH>=2){
    FillHist("TriggerEfficiency_KCDH2",0);
    if(header->IsTrig(Trig_KCDH2)){
      FillHist("TriggerEfficiency_KCDH2",1);
    }
  }
  /* KCDH2f */
  if(MulT0>0&&MulBHD>0&&MulAC==0&&MulCDH>=2){
    FillHist("TriggerEfficiency_KCDH2f",0);
    if(header->IsTrig(Trig_KCDH2f)){
      FillHist("TriggerEfficiency_KCDH2f",1);
    }
  }
  /* KCDH3 */
  if(MulT0>0&&MulBHD>0&&MulAC==0&&MulCDH>=3){
    FillHist("TriggerEfficiency_KCDH3",0);
    if(header->IsTrig(Trig_KCDH3)){
      FillHist("TriggerEfficiency_KCDH3",1);
    }
  }
  /* FWD Neutral */
  if(MulBVC+MulCVC==0&&(MulNC>0)){
    FillHist("TriggerEfficiency_FWDN",0);
    if(header->IsTrig(Trig_Neutral)){
      FillHist("TriggerEfficiency_FWDN",1);
    }
  }
  /* FWD Neutral */
  if(MulPC+MulCVC2>0){
    FillHist("TriggerEfficiency_FWDC",0);
    if(header->IsTrig(Trig_Charged)){
      FillHist("TriggerEfficiency_FWDC",1);
    }
  }

  /* ========== Trigger Check ========== */
  /* =================================== */


  return true;
}

bool MyAnalysisTrigCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisTrigCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisTrigCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisTrigCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisTrigCheck::CutCondition()
{
}

double MyAnalysisTrigCheck::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

double MyAnalysisTrigCheck::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

bool MyAnalysisTrigCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisTrigCheck::Initialize ###" << std::endl;

  confFile = confMan;
  DetectorList *dlist=DetectorList::GetInstance();

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaTrigCheck");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();
  TH1F* h1;

  /* Trigger pattern */
  for(int i=0; i<20; i++){
    h1 = new TH1F( Form("Trigger%02d_TDC",i), Form("Trigger %02d TDC",i), 4000, 0, 4000);
  }
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
  h1->GetXaxis()->SetBinLabel(10,"KCDH2/f");
  h1->GetXaxis()->SetBinLabel(11,"KCDH3");
  h1->GetXaxis()->SetBinLabel(12,"Kf");
  h1->GetXaxis()->SetBinLabel(13,"1stMix");
  h1->GetXaxis()->SetBinLabel(14,"Charged");
  h1->GetXaxis()->SetBinLabel(15,"Neutral");
  h1->GetXaxis()->SetBinLabel(16,"Cosmic");
  h1->GetXaxis()->SetBinLabel(17,"KvBVC");
  h1->GetXaxis()->SetBinLabel(18,"SIM");

  /* Trigger Efficiency */
  new TH1F("TriggerEfficiency_Beam","TriggerEfficiency of Beam",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_Kaon","TriggerEfficiency of Kaon",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_Kf","TriggerEfficiency of Kf",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_Pion","TriggerEfficiency of Pion",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_KCDH1","TriggerEfficiency of KCDH1",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_KCDH1f","TriggerEfficiency of KCDH1f",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_KCDH2","TriggerEfficiency of KCDH2",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_KCDH2f","TriggerEfficiency of KCDH2f",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_KCDH3","TriggerEfficiency of KCDH3",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_FWDN","TriggerEfficiency of FWD Neutral",2,-0.5,1.5);
  new TH1F("TriggerEfficiency_FWDC","TriggerEfficiency of FWD Charged",2,-0.5,1.5);


  return true;
}
