// MyAnalysisDetectorCheck.cpp

#include "MyAnalysisDetectorCheck.h"

MyAnalysisDetectorCheck::MyAnalysisDetectorCheck(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisDetectorCheck::~MyAnalysisDetectorCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisDetectorCheck::Clear()
{
  BeamPID = Beam_Other;
  conf = 0;
  header = 0; 
  blMan = 0;
  bltrackMan = 0;
}

bool MyAnalysisDetectorCheck::DoAnalysis(ConfMan* _conf, EventHeader* _header, BeamLineHitMan* _blMan, BeamLineTrackMan* _bltrackMan)
{
  rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();
  conf = _conf;
  header = _header;
  blMan = _blMan;
  bltrackMan = _bltrackMan;

  /* ======================================= */
  /* ========== Trigger selection ========== */
  //if(!header->IsTrig(Trig_Kaon)){
  //  return true;
  //}
  /* ========== Trigger selection ========== */
  /* ======================================= */

  /* Trigger pattern */
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }

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
  FillHist("T0_Multiplicity",MulT0,1);
  if(MulT0!=1){
    return true;
  }
  if(T0Timing<-10.0||10.0<T0Timing){
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
  FillHist("BHD_Multiplicity",MulBHD,1);
  if(MulBHD!=1){
    return true;
  }
  /* ========== Event selection ========== */
  /* ===================================== */

  /* BLDC Tracking */
  bltrackMan->DoTracking(blMan,conf,true,true);

  /* =============================== */
  /* ========== PID Check ========== */
  /* TOF */
  FillHist("T0_HitPattern",T0Segment,1);
  FillHist("BHD_HitPattern",BHDSegment,1);
  double tof = TOFCheck(T0Timing,BHDTiming);
  /* AC */
  ACCheck();
  /* ========== PID Check ========== */
  /* =============================== */

  /* =============================== */
  /* ========== DEF Check ========== */
  int def = 0;
  for( int i=0; i<blMan->nDEF(); i++ ){
    HodoscopeLikeHit *hit = blMan->DEF(i);
    if(hit->CheckRange()){
      def++;
      FillHist("DEF_HitPattern",hit->seg());
    }
  }
  FillHist("DEF_Multiplicity",def);


  /* ================================ */
  /* ========== BLDC Check ========== */
  if(BeamPID==Beam_Pion||BeamPID==Beam_Kaon){
    TString pid="";
    BLC1Check(pid.Data());
    BLC2Check(pid.Data());
    BPCCheck(pid.Data());
    FDC1Check(pid.Data());
  }
  if(BeamPID==Beam_Pion){
    TString pid="_withTOFPion";
    BLC1Check(pid.Data());
    BLC2Check(pid.Data());
    BPCCheck(pid.Data());
    FDC1Check(pid.Data());
  }
  if(BeamPID==Beam_Kaon){
    TString pid="_withTOFKaon";
    BLC1Check(pid.Data());
    BLC2Check(pid.Data());
    BPCCheck(pid.Data());
    FDC1Check(pid.Data());
  }
  /* ========== BLDC Check ========== */
  /* ================================ */

  return true;

}

bool MyAnalysisDetectorCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisDetectorCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisDetectorCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisDetectorCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

double MyAnalysisDetectorCheck::TOFCheck(double T0Timing, double BHDTiming)
{
  /* TOF check */
  double tof = T0Timing - BHDTiming;
  bool hitflagtof = false;
  for(int i=0; i<bltrackMan->ntrackBLC2(); i++){
    LocalTrack* track = bltrackMan->trackBLC2(i);
    double tracktime = track->GetTrackTime();
    if(tracktime<-10||10<tracktime) continue;
    double chi2 = track->chi2all();
    if(20<chi2) continue;
    TVector3 pos = track->GetPosatZ(-75.0);
    if(pos.X()<-2.0||2.0<pos.X()) continue;
    if(pos.Y()<-2.0||2.0<pos.Y()) continue;
    hitflagtof = true;
  }
  if(hitflagtof){
    FillHist("BHDT0_TOF",tof,1);
    if(header->IsTrig(Trig_Beam)){
      FillHist("BHDT0_TOF_withTrigBeam",tof,1);
    }
    if(header->IsTrig(Trig_Pion)){
      FillHist("BHDT0_TOF_withTrigPion",tof,1);
    }
    if(header->IsTrig(Trig_Kaon)){
      FillHist("BHDT0_TOF_withTrigKaon",tof,1);
    }
    if(header->IsTrig(Trig_Pion)&&header->IsTrig(Trig_Kaon)){
      FillHist("BHDT0_TOF_withTrigPionTrigKaon",tof,1);
    }
  }
  /* Beam particle decision */
  if(25.0<tof&&tof<27.0){
    BeamPID=Beam_Pion;
  }
  if(27.5<tof&&tof<29.5){
    BeamPID=Beam_Kaon;
  }
  if(34.0<tof&&tof<35.0){
    BeamPID=Beam_Proton;
  }
  return tof;
}

void MyAnalysisDetectorCheck::ACCheck()
{
  bool achitflag = false;
  TVector3 acpos;
  for(int i=0; i<bltrackMan->ntrackBLC2(); i++){
    LocalTrack* track = bltrackMan->trackBLC2(i);
    double tracktime = track->GetTrackTime();
    if(tracktime<-10||10<tracktime) continue;
    double chi2 = track->chi2all();
    if(20<chi2) continue;
    TVector3 pos = track->GetPosatZ(-75.0);
    if(pos.X()<-2.0||2.0<pos.X()) continue;
    if(pos.Y()<-2.0||2.0<pos.Y()) continue;
    achitflag = true;
    acpos = pos;
  }

  int adcsum=0;
  for(int i=0; i<blMan->nChere(CID_AC);i++){
    CherenkovLikeHit* hit = blMan->Cherei(CID_AC,i);
    if(hit->seg()!=1) continue;
    for(int ich=1; ich<=4; ich++){
      int adc=hit->adc(ich);
      int tdc=hit->tdc(ich);
      adcsum+=adc;
      FillHist("AC_TDC",tdc,1);
      if(hit->CheckRange(1)){
        FillHist("AC_TDC_withHit",tdc,1);
      }
      else{
        FillHist("AC_TDC_withoutHit",tdc,1);
      }
    }
    FillHist("AC_ADCSUM",adcsum,1);
    FillHist("AC_HitPosition",acpos.X(),acpos.Y(),1);
    if(hit->CheckRange(1)){
      FillHist("AC_ADCSUM_withHit",adcsum,1);
      FillHist("AC_HitPosition_withHit",acpos.X(),acpos.Y(),1);
    }
    else{
      FillHist("AC_ADCSUM_withoutHit",adcsum,1);
      FillHist("AC_HitPosition_withoutHit",acpos.X(),acpos.Y(),1);
    }
  }

  if(header->IsTrig(Trig_Pion)){
    FillHist("AC_HitPosition_withTrigPion",acpos.X(),acpos.Y(),1);
  }
  if(header->IsTrig(Trig_Kaon)){
    FillHist("AC_HitPosition_withTrigKaon",acpos.X(),acpos.Y(),1);
  }
  if(header->IsTrig(Trig_Pion)&&header->IsTrig(Trig_Kaon)){
    FillHist("AC_HitPosition_withTrigKaon",acpos.X(),acpos.Y(),1);
  }

  if(!achitflag) return;
  if(header->IsTrig(Trig_Beam)){
    FillHist("AC_ADCSUM_withTrigBeam",adcsum,1);
  }
  /* AC if TrigPion */
  if(header->IsTrig(Trig_Pion)){
    FillHist("AC_ADCSUM_withTrigPion",adcsum,1);
  }
  /* AC if TrigKaon */
  if(header->IsTrig(Trig_Kaon)){
    FillHist("AC_ADCSUM_withTrigKaon",adcsum,1);
  }
  /* AC if TrigPion & TrigKaon */
  if(header->IsTrig(Trig_Pion)&&header->IsTrig(Trig_Kaon)){
    FillHist("AC_ADCSUM_withTrigPionTrigKaon",adcsum,1);
  }
  /* AC if TOFPion */
  if(BeamPID==Beam_Pion){
    FillHist("AC_ADCSUM_withTOFPion",adcsum,1);
    if(header->IsTrig(Trig_Pion)){
      FillHist("AC_ADCSUM_withTOFPionTrigPion",adcsum,1);
    }
    if(header->IsTrig(Trig_Kaon)){
      FillHist("AC_ADCSUM_withTOFPionTrigKaon",adcsum,1);
    }
    if(header->IsTrig(Trig_Pion)&&header->IsTrig(Trig_Kaon)){
      FillHist("AC_ADCSUM_withTOFPionTrigPionTrigKaon",adcsum,1);
    }
  }
  /* AC if TOFKaon */
  if(BeamPID==Beam_Kaon){
    FillHist("AC_ADCSUM_withTOFKaon",adcsum,1);
    if(header->IsTrig(Trig_Pion)){
      FillHist("AC_ADCSUM_withTOFKaonTrigPion",adcsum,1);
    }
    if(header->IsTrig(Trig_Kaon)){
      FillHist("AC_ADCSUM_withTOFKaonTrigKaon",adcsum,1);
    }
    if(header->IsTrig(Trig_Pion)&&header->IsTrig(Trig_Kaon)){
      FillHist("AC_ADCSUM_withTOFKaonTrigPionTrigKaon",adcsum,1);
    }
  }
  /* AC if TOFProton */
  if(BeamPID==Beam_Proton){
    FillHist("AC_ADCSUM_withTOFProton",adcsum,1);
    if(header->IsTrig(Trig_Pion)){
      FillHist("AC_ADCSUM_withTOFProtonTrigPion",adcsum,1);
    }
    if(header->IsTrig(Trig_Kaon)){
      FillHist("AC_ADCSUM_withTOFProtonTrigKaon",adcsum,1);
    }
    if(header->IsTrig(Trig_Pion)&&header->IsTrig(Trig_Kaon)){
      FillHist("AC_ADCSUM_withTOFProtonTrigPionTrigKaon",adcsum,1);
    }
  }

}

void MyAnalysisDetectorCheck::BLC1Check(const char* pid="")
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(bltrackMan->ntrackBLC2()!=1||bltrackMan->ntrackBPC()!=1) return;
  LocalTrack* blc2 = bltrackMan->trackBLC2(0);
  LocalTrack* bpc  = bltrackMan->trackBPC(0);
  double blc2tracktime = blc2->GetTrackTime();
  double bpctracktime  = bpc->GetTrackTime();
  if(blc2tracktime<-10||10<blc2tracktime) return;
  if(bpctracktime<-10||10<bpctracktime) return;
  double blc2chi2 = blc2->chi2all();
  double bpcchi2  = bpc->chi2all();
  if(10<blc2chi2||10<bpcchi2) return;
  FillHist(Form("BLC1_TrackingEfficiency%s",pid),0,1);
  FillHist(Form("BLC1a_TrackingEfficiency%s",pid),0,1);
  FillHist(Form("BLC1b_TrackingEfficiency%s",pid),0,1);
  /* Basic data filling */
  const int ndc=2;
  int dccid[ndc]={ CID_BLC1a,CID_BLC1b};
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    for(int layer=1;layer<=nlays;layer++){
      int mult = blMan->nBLDC(cid,layer);
      if(mult>0)
        FillHist(Form("%sl%d_Multiplicity%s",name,layer,pid),mult);
      for(int i=0;i<mult;i++){
        ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
        int wire = hit->wire();
        double dt = hit->dt();
        FillHist(Form("%sl%d_HitPattern%s",name,layer,pid),wire,1);
        FillHist(Form("%sl%d_Time%s",name,layer,pid),dt,1);
      }
    }
    for(int i=0; i<bltrackMan->ntrackBLDC(cid); i++){
      LocalTrack* track=bltrackMan->trackBLDC(cid,i);
      for(int ihit=0; ihit<track->nhit(); ihit++){
        ChamberLikeHit* hit = track->hit(ihit);
        int layer = hit->layer();
        double resl = hit->resl();
        FillHist(Form("%sl%d_Residual%s",name,layer,pid),resl,1);
      }
    }
  }

  /* BLC1 */
  int weight[4] = {1,1,1,1};
  FillHist(Form("BLC1_NumberofTrack%s",pid),bltrackMan->ntrackBLC1(),1);
  for(int i=0; i<bltrackMan->ntrackBLC1(); i++){
    LocalTrack* blc1 = bltrackMan->trackBLC1(i);
    double tracktime = blc1->GetTrackTime();
    double chi2 = blc1->chi2all();
    FillHist(Form("BLC1_TrackTime%s",pid),tracktime,1);
    /* Tracking efficiency */
    FillHist(Form("BLC1_TrackingEfficiency%s",pid),1,weight[0]); /* No selection */
    weight[0]=0;
    if((-10<tracktime&&tracktime<10)){
      FillHist(Form("BLC1_ChiSquare%s",pid),chi2,1);
      FillHist(Form("BLC1_TrackingEfficiency%s",pid),2,weight[1]); /* Timing selection */
      weight[1]=0;
    }
    if((chi2<10)){
      FillHist(Form("BLC1_TrackingEfficiency%s",pid),3,weight[2]); /* chi2 selection */
      weight[2]=0;
    }
    if((-10<tracktime&&tracktime<10)&&(chi2<10)){
      FillHist(Form("BLC1_TrackingEfficiency%s",pid),4,weight[3]); /* Timing & chi2 selection */
      weight[3]=0;
    }
  }

  /* BLC1a */
  weight[0]=1;weight[1]=1;weight[2]=1;weight[3]=1;
  FillHist(Form("BLC1a_NumberofTrack%s",pid),bltrackMan->ntrackBLC1a(),1);
  for(int i=0; i<bltrackMan->ntrackBLC1a(); i++){
    LocalTrack* blc1a = bltrackMan->trackBLC1a(i);
    double tracktime = blc1a->GetTrackTime();
    double chi2 = blc1a->chi2all();
    FillHist(Form("BLC1a_TrackTime%s",pid),tracktime,1);
    /* Tracking efficiency */
    FillHist(Form("BLC1a_TrackingEfficiency%s",pid),1,weight[0]); /* No selection */
    weight[0]=0;
    if((-10<tracktime&&tracktime<10)){
      FillHist(Form("BLC1a_ChiSquare%s",pid),chi2,1);
      FillHist(Form("BLC1a_TrackingEfficiency%s",pid),2,weight[1]); /* Timing selection */
      weight[1]=0;
    }
    if((chi2<10)){
      FillHist(Form("BLC1a_TrackingEfficiency%s",pid),3,weight[2]); /* chi2 selection */
      weight[2]=0;
    }
    if((-10<tracktime&&tracktime<10)&&(chi2<10)){
      FillHist(Form("BLC1a_TrackingEfficiency%s",pid),4,weight[3]); /* Timing & chi2 selection */
      weight[3]=0;
      /* BLC1b Layer Hit Efficiencly */
      FillHist(Form("BLC1b_HitEfficiency%s",pid),0,1);
      for(int layer=1; layer<=8; layer++){
        if(blMan->nBLC1b(layer)>0){
          FillHist(Form("BLC1b_HitEfficiency%s",pid),layer,1);
        }
      }
    }
  }
  /* BLC1b */
  weight[0]=1;weight[1]=1;weight[2]=1;weight[3]=1;
  FillHist(Form("BLC1b_NumberofTrack%s",pid),bltrackMan->ntrackBLC1b(),1);
  for(int i=0; i<bltrackMan->ntrackBLC1b(); i++){
    LocalTrack* blc1b = bltrackMan->trackBLC1b(i);
    double tracktime = blc1b->GetTrackTime();
    double chi2 = blc1b->chi2all();
    FillHist(Form("BLC1b_TrackTime%s",pid),tracktime,1);
    /* Tracking efficiency */
    FillHist(Form("BLC1b_TrackingEfficiency%s",pid),1,weight[0]); /* No selection */
    weight[0]=0;
    if((-10<tracktime&&tracktime<10)){
      FillHist(Form("BLC1b_ChiSquare%s",pid),chi2,1);
      FillHist(Form("BLC1b_TrackingEfficiency%s",pid),2,weight[1]); /* Timing selection */
      weight[1]=0;
    }
    if((chi2<10)){
      FillHist(Form("BLC1b_TrackingEfficiency%s",pid),3,weight[2]); /* chi2 selection */
      weight[2]=0;
    }
    if((-10<tracktime&&tracktime<10)&&(chi2<10)){
      FillHist(Form("BLC1b_TrackingEfficiency%s",pid),4,weight[3]); /* Timing & chi2 selection */
      weight[3]=0;
      /* BLC1a Layer Hit Efficiencly */
      FillHist(Form("BLC1a_HitEfficiency%s",pid),0,1);
      for(int layer=1; layer<=8; layer++){
        if(blMan->nBLC1a(layer)>0){
          FillHist(Form("BLC1a_HitEfficiency%s",pid),layer,1);
        }
      }
    }
  }

  return;
}

void MyAnalysisDetectorCheck::BLC2Check(const char* pid="")
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(bltrackMan->ntrackBLC1()!=1||bltrackMan->ntrackBPC()!=1) return;
  LocalTrack* blc1 = bltrackMan->trackBLC1(0);
  LocalTrack* bpc  = bltrackMan->trackBPC(0);
  double blc1tracktime = blc1->GetTrackTime();
  double bpctracktime  = bpc->GetTrackTime();
  if(blc1tracktime<-10||10<blc1tracktime) return;
  if(bpctracktime<-10||10<bpctracktime) return;
  double blc1chi2 = blc1->chi2all();
  double bpcchi2  = bpc->chi2all();
  if(10<blc1chi2||10<bpcchi2) return;
  FillHist(Form("BLC2_TrackingEfficiency%s",pid),0,1);
  FillHist(Form("BLC2a_TrackingEfficiency%s",pid),0,1);
  FillHist(Form("BLC2b_TrackingEfficiency%s",pid),0,1);
  /* Basic data filling */
  const int ndc=2;
  int dccid[ndc]={ CID_BLC2a,CID_BLC2b};
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    for(int layer=1;layer<=nlays;layer++){
      int mult = blMan->nBLDC(cid,layer);
      if(mult>0)
        FillHist(Form("%sl%d_Multiplicity%s",name,layer,pid),mult);
      for(int i=0;i<mult;i++){
        ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
        int wire = hit->wire();
        double dt = hit->dt();
        FillHist(Form("%sl%d_HitPattern%s",name,layer,pid),wire,1);
        FillHist(Form("%sl%d_Time%s",name,layer,pid),dt,1);
      }
    }
    for(int i=0; i<bltrackMan->ntrackBLDC(cid); i++){
      LocalTrack* track=bltrackMan->trackBLDC(cid,i);
      for(int ihit=0; ihit<track->nhit(); ihit++){
        ChamberLikeHit* hit = track->hit(ihit);
        int layer = hit->layer();
        double resl = hit->resl();
        FillHist(Form("%sl%d_Residual%s",name,layer,pid),resl,1);
      }
    }
  }


  /* BLC2 */
  int weight[4] = {1,1,1,1};
  FillHist(Form("BLC2_NumberofTrack%s",pid),bltrackMan->ntrackBLC2(),1);
  for(int i=0; i<bltrackMan->ntrackBLC2(); i++){
    LocalTrack* blc2 = bltrackMan->trackBLC2(i);
    double tracktime = blc2->GetTrackTime();
    double chi2 = blc2->chi2all();
    FillHist(Form("BLC2_TrackTime%s",pid),tracktime,1);
    /* Tracking efficiency */
    FillHist(Form("BLC2_TrackingEfficiency%s",pid),1,weight[0]); /* No selection */
    weight[0]=0;
    if((-10<tracktime&&tracktime<10)){
      FillHist(Form("BLC2_ChiSquare%s",pid),chi2,1);
      FillHist(Form("BLC2_TrackingEfficiency%s",pid),2,weight[1]); /* Timing selection */
      weight[1]=0;
    }
    if((chi2<10)){
      FillHist(Form("BLC2_TrackingEfficiency%s",pid),3,weight[2]); /* chi2 selection */
      weight[2]=0;
    }
    if((-10<tracktime&&tracktime<10)&&(chi2<10)){
      FillHist(Form("BLC2_TrackingEfficiency%s",pid),4,weight[3]); /* Timing & chi2 selection */
      weight[3]=0;
    }
  }

  /* BLC2a */
  weight[0]=1;weight[1]=1;weight[2]=1;weight[3]=1;
  FillHist(Form("BLC2a_NumberofTrack%s",pid),bltrackMan->ntrackBLC2a(),1);
  for(int i=0; i<bltrackMan->ntrackBLC2a(); i++){
    LocalTrack* blc2a = bltrackMan->trackBLC2a(i);
    double tracktime = blc2a->GetTrackTime();
    double chi2 = blc2a->chi2all();
    FillHist(Form("BLC2a_TrackTime%s",pid),tracktime,1);
    /* Tracking efficiency */
    FillHist(Form("BLC2a_TrackingEfficiency%s",pid),1,weight[0]); /* No selection */
    weight[0]=0;
    if((-10<tracktime&&tracktime<10)){
      FillHist(Form("BLC2a_ChiSquare%s",pid),chi2,1);
      FillHist(Form("BLC2a_TrackingEfficiency%s",pid),2,weight[1]); /* Timing selection */
      weight[1]=0;
    }
    if((chi2<10)){
      FillHist(Form("BLC2a_TrackingEfficiency%s",pid),3,weight[2]); /* chi2 selection */
      weight[2]=0;
    }
    if((-10<tracktime&&tracktime<10)&&(chi2<10)){
      FillHist(Form("BLC2a_TrackingEfficiency%s",pid),4,weight[3]); /* Timing & chi2 selection */
      weight[3]=0;
      /* BLC2b Layer Hit Efficiencly */
      FillHist(Form("BLC2b_HitEfficiency%s",pid),0,1);
      for(int layer=1; layer<=8; layer++){
        if(blMan->nBLC2b(layer)>0){
          FillHist(Form("BLC2b_HitEfficiency%s",pid),layer,1);
        }
      }
    }
  }
  /* BLC2b */
  weight[0]=1;weight[1]=1;weight[2]=1;weight[3]=1;
  FillHist(Form("BLC2b_NumberofTrack%s",pid),bltrackMan->ntrackBLC2b(),1);
  for(int i=0; i<bltrackMan->ntrackBLC2b(); i++){
    LocalTrack* blc2b = bltrackMan->trackBLC2b(i);
    double tracktime = blc2b->GetTrackTime();
    double chi2 = blc2b->chi2all();
    FillHist(Form("BLC2b_TrackTime%s",pid),tracktime,1);
    /* Tracking efficiency */
    FillHist(Form("BLC2b_TrackingEfficiency%s",pid),1,weight[0]); /* No selection */
    weight[0]=0;
    if((-10<tracktime&&tracktime<10)){
      FillHist(Form("BLC2b_ChiSquare%s",pid),chi2,1);
      FillHist(Form("BLC2b_TrackingEfficiency%s",pid),2,weight[1]); /* Timing selection */
      weight[1]=0;
    }
    if((chi2<10)){
      FillHist(Form("BLC2b_TrackingEfficiency%s",pid),3,weight[2]); /* chi2 selection */
      weight[2]=0;
    }
    if((-10<tracktime&&tracktime<10)&&(chi2<10)){
      FillHist(Form("BLC2b_TrackingEfficiency%s",pid),4,weight[3]); /* Timing & chi2 selection */
      weight[3]=0;
      /* BLC2a Layer Hit Efficiencly */
      FillHist(Form("BLC2a_HitEfficiency%s",pid),0,1);
      for(int layer=1; layer<=8; layer++){
        if(blMan->nBLC2a(layer)>0){
          FillHist(Form("BLC2a_HitEfficiency%s",pid),layer,1);
        }
      }
    }
  }

  return;
}

void MyAnalysisDetectorCheck::BPCCheck(const char* pid="")
{
  DetectorList *dlist=DetectorList::GetInstance();
  FillHist("BPC_Number",0,1);
  if(bltrackMan->ntrackBLC1()!=1||bltrackMan->ntrackBLC2()!=1||bltrackMan->ntrackFDC1()!=1) return;
  FillHist("BPC_Number",1,1);
  LocalTrack* blc1 = bltrackMan->trackBLC1(0);
  LocalTrack* blc2 = bltrackMan->trackBLC2(0);
  LocalTrack* fdc1 = bltrackMan->trackFDC1(0);
  double blc1tracktime = blc1->GetTrackTime();
  double blc2tracktime = blc2->GetTrackTime();
  double fdc1tracktime = fdc1->GetTrackTime();
  if(blc1tracktime<-10||10<blc1tracktime) return;
  if(blc2tracktime<-10||10<blc2tracktime) return;
  if(fdc1tracktime<-10||10<fdc1tracktime) return;
  FillHist("BPC_Number",2,1);
  double blc1chi2 = blc1->chi2all();
  double blc2chi2 = blc2->chi2all();
  double fdc1chi2 = blc2->chi2all();
  if(10<blc1chi2||10<blc2chi2||10<fdc1chi2) return;
  FillHist("BPC_Number",3,1);
  TVector3 blc2dir = blc2->GetMomDir();
  if(0.1<blc2dir.X()/blc2dir.Z()) return;
  if(0.1<blc2dir.Y()/blc2dir.Z()) return;
  FillHist("BPC_Number",4,1);
  TVector3 blc2pos = blc2->GetPosatZ(-20.3);
  if(blc2pos.X()<-2.0||2.0<blc2pos.X()) return;
  if(blc2pos.Y()<-2.0||2.0<blc2pos.Y()) return;
  FillHist("BPC_Number",5,1);
  int def=0;
  for(int i=0; i<blMan->nDEF(); i++){
    if(blMan->DEF(i)->CheckRange()){
      def++;
    }
  }
  if(def==0) return;
  FillHist("BPC_Number",6,1);
  FillHist(Form("BPC_TrackingEfficiency%s",pid),0,1);
  /* Basic data filling */
  const int cid=CID_BPC;
  const char* name= dlist->GetName(cid).data();
  const int nlays= dlist->GetNlayers(cid);
  //std::cout << "=== All ===" << std::endl;
  for(int layer=1;layer<=nlays;layer++){
    int mult = blMan->nBLDC(cid,layer);
    if(mult>0){
      FillHist(Form("%sl%d_Multiplicity%s",name,layer,pid),mult);
    }
    //std::cout << "layer" << layer << " : ";
    int wir[15] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    for(int i=0;i<mult;i++){
      ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
      int wire = hit->wire();
      wir[wire-1]++;
      double dt = hit->dt();
      FillHist(Form("%sl%d_HitPattern%s",name,layer,pid),wire,1);
      FillHist(Form("%sl%d_Time%s",name,layer,pid),dt,1);
    }
    for(int wire=1; wire<=15; wire++){
      //if(wir[wire-1]>0) {std::cout << "*";}
      //else {std::cout << "|";}
    }
    //std::cout << std::endl;
  }
  //std::cout << "=== All ===" << std::endl;
  for(int i=0; i<bltrackMan->ntrackBLDC(cid); i++){
    LocalTrack* track=bltrackMan->trackBLDC(cid,i);
    int wir[8][15]={0};
    for(int ihit=0; ihit<track->nhit(); ihit++){
      ChamberLikeHit* hit = track->hit(ihit);
      int layer = hit->layer();
      int wire = hit->wire();
      wir[layer-1][wire-1]++;
      double resl = hit->resl();
      FillHist(Form("%sl%d_Residual%s",name,layer,pid),resl,1);
    }
    //std::cout << "===Track===" << std::endl;
    for(int layer=1; layer<=8; layer++){
      //std::cout << "layer" << layer << " : ";
      for(int wire=1; wire<=15; wire++){
        //if(wir[layer-1][wire-1]>0) {std::cout << "*";}
        //else {std::cout << "|";}
      }
      //std::cout << std::endl;
    }
    //std::cout << "===Track===" << std::endl;
  }

  /* BPC */
  int weight[4] = {1,1,1,1};
  FillHist(Form("BPC_NumberofTrack%s",pid),bltrackMan->ntrackBPC(),1);
  for(int i=0; i<bltrackMan->ntrackBPC(); i++){
    LocalTrack* bpc = bltrackMan->trackBPC(i);
    double tracktime = bpc->GetTrackTime();
    double chi2 = bpc->chi2all();
    FillHist(Form("BPC_TrackTime%s",pid),tracktime,1);
    /* Tracking efficiency */
    FillHist(Form("BPC_TrackingEfficiency%s",pid),1,weight[0]); /* No selection */
    weight[0]=0;
    if((-10<tracktime&&tracktime<10)){
      FillHist(Form("BPC_ChiSquare%s",pid),chi2,1);
      FillHist(Form("BPC_TrackingEfficiency%s",pid),2,weight[1]); /* Timing selection */
      weight[1]=0;
    }
    if((chi2<10)){
      FillHist(Form("BPC_TrackingEfficiency%s",pid),3,weight[2]); /* chi2 selection */
      weight[2]=0;
    }
    if((-10<tracktime&&tracktime<10)&&(chi2<10)){
      FillHist(Form("BPC_TrackingEfficiency%s",pid),4,weight[3]); /* Timing & chi2 selection */
      weight[3]=0;
    }
  }
  /* Layer Hit Efficiency */
  FillHist(Form("BPC_HitEfficiency%s",pid),0,1);
  for(int layer=1; layer<=8; layer++){
    if(blMan->nBPC(layer)>0){
      FillHist(Form("BPC_HitEfficiency%s",pid),layer,1);
    }
  }

  return;
}

void MyAnalysisDetectorCheck::FDC1Check(const char* pid="")
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(bltrackMan->ntrackBLC1()!=1||bltrackMan->ntrackBLC2()!=1||bltrackMan->ntrackBPC()!=1) return;
  LocalTrack* blc1 = bltrackMan->trackBLC1(0);
  LocalTrack* blc2 = bltrackMan->trackBLC2(0);
  LocalTrack* bpc  = bltrackMan->trackBPC (0);
  double blc1tracktime = blc1->GetTrackTime();
  double blc2tracktime = blc2->GetTrackTime();
  double bpctracktime  = bpc->GetTrackTime();
  if(blc1tracktime<-10||10<blc1tracktime) return;
  if(blc2tracktime<-10||10<blc2tracktime) return;
  if(bpctracktime<-10||10<bpctracktime) return;
  double blc1chi2 = blc1->chi2all();
  double blc2chi2 = blc2->chi2all();
  double bpcchi2 = bpc->chi2all();
  if(10<blc1chi2||10<blc2chi2||10<bpcchi2) return;
  TVector3 blc1dir = blc1->GetMomDir();
  TVector3 blc2dir = blc2->GetMomDir();
  TVector3 bpcdir = bpc->GetMomDir();
  FillHist("BLC1_TrackDirection",blc1dir.X()/blc1dir.Z(),blc1dir.Y()/blc1dir.Z(),1);
  FillHist("BLC2_TrackDirection",blc2dir.X()/blc2dir.Z(),blc2dir.Y()/blc2dir.Z(),1);
  FillHist("BPC_TrackDirection",bpcdir.X()/bpcdir.Z(),bpcdir.Y()/bpcdir.Z(),1);
  if(0.1<bpcdir.X()/bpcdir.Z()) return;
  if(0.1<bpcdir.Y()/bpcdir.Z()) return;
  int bvc=0;
  for(int i=0; i<blMan->nBVC(); i++){
    if(blMan->BVC(i)->CheckRange()){
      bvc++;
    }
  }
  int cvc=0;
  for(int i=0; i<blMan->nCVC(); i++){
    if(blMan->CVC(i)->CheckRange()){
      cvc++;
    }
  }
  int pc=0;
  for(int i=0; i<blMan->nPC(); i++){
    if(blMan->PC(i)->CheckRange()){
      pc++;
    }
  }
  int bd=0;
  for(int i=0; i<blMan->nBD(); i++){
    if(blMan->BD(i)->CheckRange()){
      bd++;
    }
  }
  if(bvc==0) return;
  if(cvc==0&&pc==0&&bd==0) return;
  FillHist(Form("FDC1_TrackingEfficiency%s",pid),0,1);
  /* Basic data filling */
  const int cid=CID_FDC1;
  const char* name= dlist->GetName(cid).data();
  const int nlays= dlist->GetNlayers(cid);
  for(int layer=1;layer<=nlays;layer++){
    int mult = blMan->nBLDC(cid,layer);
    if(mult>0){
      FillHist(Form("%sl%d_Multiplicity%s",name,layer,pid),mult);
    }
    for(int i=0;i<mult;i++){
      ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
      int wire = hit->wire();
      double dt = hit->dt();
      FillHist(Form("%sl%d_HitPattern%s",name,layer,pid),wire,1);
      FillHist(Form("%sl%d_Time%s",name,layer,pid),dt,1);
    }
  }
  for(int i=0; i<bltrackMan->ntrackBLDC(cid); i++){
    LocalTrack* track=bltrackMan->trackBLDC(cid,i);
    for(int ihit=0; ihit<track->nhit(); ihit++){
      ChamberLikeHit* hit = track->hit(ihit);
      int layer = hit->layer();
      double resl = hit->resl();
      FillHist(Form("%sl%d_Residual%s",name,layer,pid),resl,1);
    }
  }

  /* FDC1 */
  int weight[4] = {1,1,1,1};
  FillHist(Form("FDC1_NumberofTrack%s",pid),bltrackMan->ntrackFDC1(),1);
  for(int i=0; i<bltrackMan->ntrackFDC1(); i++){
    LocalTrack* fdc1 = bltrackMan->trackFDC1(i);
    double tracktime = fdc1->GetTrackTime();
    double chi2 = fdc1->chi2all();
    FillHist(Form("FDC1_TrackTime%s",pid),tracktime,1);
    /* Tracking efficiency */
    FillHist(Form("FDC1_TrackingEfficiency%s",pid),1,weight[0]); /* No selection */
    weight[0]=0;
    if((-10<tracktime&&tracktime<10)){
      FillHist(Form("FDC1_ChiSquare%s",pid),chi2,1);
      FillHist(Form("FDC1_TrackingEfficiency%s",pid),2,weight[1]); /* Timing selection */
      weight[1]=0;
    }
    if((chi2<10)){
      FillHist(Form("FDC1_TrackingEfficiency%s",pid),3,weight[2]); /* chi2 selection */
      weight[2]=0;
    }
    if((-10<tracktime&&tracktime<10)&&(chi2<10)){
      FillHist(Form("FDC1_TrackingEfficiency%s",pid),4,weight[3]); /* Timing & chi2 selection */
      weight[3]=0;
    }
  }
  /* Layer Hit Efficiency */
  FillHist(Form("FDC1_HitEfficiency%s",pid),0,1);
  for(int layer=1; layer<=6; layer++){
    if(blMan->nFDC1(layer)>0){
      FillHist(Form("FDC1_HitEfficiency%s",pid),layer,1);
    }
  }

  return;
}

void MyAnalysisDetectorCheck::CutCondition()
{
}

double MyAnalysisDetectorCheck::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

double MyAnalysisDetectorCheck::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

bool MyAnalysisDetectorCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisDetectorCheck::Initialize ###" << std::endl;

  confFile = confMan;
  DetectorList *dlist=DetectorList::GetInstance();

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaDetectorCheck");

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

  /* T0 and BHD */
  new TH1F( Form("T0_HitPattern"), Form("Hit pattern T0;T0 Segment;Counts"), 5, 0.5, 5+0.5);
  new TH1F( Form("T0_Multiplicity"), Form("Multiplicity T0;T0 Multiplicity;Counts"), 5+1, -0.5, 5+0.5 );
  new TH1F( Form("BHD_HitPattern"), Form("Hit pattern BHD;BHD Segment;Counts"), 20, 0.5, 20+0.5);
  new TH1F( Form("BHD_Multiplicity"), Form("Multiplicity BHD;BHD Multiplicity;Counts"), 20+1, -0.5, 20+0.5 );
  new TH1F( Form("DEF_HitPattern"), Form("Hit pattern DEF;DEF Segment;Counts"), 8, 0.5, 8+0.5);
  new TH1F( Form("DEF_Multiplicity"), Form("Multiplicity DEF;DEF Multiplicity;Counts"), 8+1, -0.5, 8+0.5 );

  /* TOF */
  new TH1F( Form("BHDT0_TOF"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);
  new TH1F( Form("BHDT0_TOF_withTrigBeam"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);
  new TH1F( Form("BHDT0_TOF_withTrigPion"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);
  new TH1F( Form("BHDT0_TOF_withTrigKaon"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);
  new TH1F( Form("BHDT0_TOF_withTrigPionTrigKaon"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);

  /* AC */
  new TH1F( Form("AC_TDC"),Form("AC TDC;TDC (ch.);Counts"),4000,0,4000);
  new TH1F( Form("AC_TDC_withHit"),Form("AC TDC;TDC (ch.);Counts"),4000,0,4000);
  new TH1F( Form("AC_TDC_withoutHit"),Form("AC TDC;TDC (ch.);Counts"),4000,0,4000);
  new TH1F( Form("AC_ADCSUM"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withHit"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withoutHit"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTrigBeam"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTrigPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTrigPionTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFPionTrigPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFPionTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFPionTrigPionTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFKaonTrigPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFKaonTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFKaonTrigPionTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFProton"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFProtonTrigPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFProtonTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUM_withTOFProtonTrigPionTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH2F( Form("AC_HitPosition"),Form("AC Hit Position;X Hit position (cm);Y Hit position (cm)"),400,-20,20,400,-20,20);
  new TH2F( Form("AC_HitPosition_withHit"),Form("AC Hit Position;X Hit position (cm);Y Hit position (cm)"),400,-20,20,400,-20,20);
  new TH2F( Form("AC_HitPosition_withoutHit"),Form("AC Hit Position;X Hit position (cm);Y Hit position (cm)"),400,-20,20,400,-20,20);
  new TH2F( Form("AC_HitPosition_withTrigPion"),Form("AC Hit Position;X Hit position (cm);Y Hit position (cm)"),400,-20,20,400,-20,20);
  new TH2F( Form("AC_HitPosition_withTrigKaon"),Form("AC Hit Position;X Hit position (cm);Y Hit position (cm)"),400,-20,20,400,-20,20);
  new TH2F( Form("AC_HitPosition_withTrigPionTrigKaon"),Form("AC Hit Position;X Hit position (cm);Y Hit position (cm)"),400,-20,20,400,-20,20);

  /* BLDC */
  const int ndc=8;
  int dccid[ndc]={CID_BLC1,CID_BLC2,
    CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_BPC,CID_FDC1
  };
  const int npid=3;
  const char* pid[3] = {"","_withTOFPion","_withTOFKaon"};
  for(int idc=0;idc<2;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    for(int ipid=0; ipid<npid; ipid++){
      new TH1F( Form("%s_TrackTime%s",name,pid[ipid]),Form("Track time %s;%s Track time (ns);Counts",name,name),1400,-200,500);
      new TH1F( Form("%s_ChiSquare%s",name,pid[ipid]),Form("#Chi^{2}/NDF %s;%s #Chi^{2}/NDF;Counts",name,name),500,0,100);
      new TH1F( Form("%s_NumberofTrack%s",name,pid[ipid]),Form("Number of tracks %s;%s Number of tracks;Counts",name,name),11,-0.5,10.5);
    }
    new TH2F( Form("%s_TrackPosition",name),Form("Track position %s;%s X track position (cm);%s Y track position (cm)",name,name,name),400,-20,20,400,-20,20);
    new TH2F( Form("%s_TrackDirection",name),Form("Track position %s;%s X track direction;%s Y track direction",name,name,name),400,-0.02,0.02,400,-0.02,0.02);
  }
  for(int idc=2;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    for(int ipid=0; ipid<npid; ipid++){
      new TH1F( Form("%s_TrackTime%s",name,pid[ipid]),Form("Track time %s;%s Track time (ns);Counts",name,name),1400,-200,500);
      new TH1F( Form("%s_ChiSquare%s",name,pid[ipid]),Form("#Chi^{2}/NDF %s;%s #Chi^{2}/NDF;Counts",name,name),500,0,100);
      new TH1F( Form("%s_HitEfficiency%s",name,pid[ipid]),Form("Layer hit efficiency %s;%s Layer;Event",name,name),nlays+1,-0.5,nlays+0.5);
      new TH1F( Form("%s_NumberofTrack%s",name,pid[ipid]),Form("Number of tracks %s;%s Number of tracks;Counts",name,name),11,-0.5,10.5);
    }
    new TH2F( Form("%s_TrackPosition",name),Form("Track position %s;%s X track position (cm);%s Y track position (cm)",name,name,name),400,-20,20,400,-20,20);
    new TH2F( Form("%s_TrackDirection",name),Form("Track position %s;%s X track direction;%s Y track direction",name,name,name),400,-0.02,0.02,400,-0.02,0.02);
    for(int ipid=0; ipid<npid; ipid++){
      for(int layer=1; layer<=nlays; layer++){
        int nwires=dlist->GetNwires(cid);
        new TH1F( Form("%sl%d_Multiplicity%s",name,layer,pid[ipid]),Form("Multiplicity %s-L%d;%s-L%d Multiplicity;Counts",name,layer,name,layer),nwires+1,-0.5,nwires+0.5);
        new TH1F( Form("%sl%d_HitPattern%s",name,layer,pid[ipid]),Form("Hit pattern %s-L%d;%s-L%d Wire;Counts",name,layer,name,layer),nwires+1,-0.5,nwires+0.5);
        new TH1F( Form("%sl%d_Time%s",name,layer,pid[ipid]),Form("Time %s-L%d;%s-L%d Time (ns);Counts",name,layer,name,layer),500,-100,400);
        new TH1F( Form("%sl%d_Residual%s",name,layer,pid[ipid]),Form("Residual %s-L%d;%s-L%d Residual (ns);Counts",name,layer,name,layer),200,-0.2,0.2);
      }
    }
  }
  for(int ipid=0; ipid<npid; ipid++){
    new TH1F( Form("BLC1_TrackingEfficiency%s",pid[ipid]),Form("Tracking efficiency BLC1;Layer;Event"),5,-0.5,4.5);
    new TH1F( Form("BLC2_TrackingEfficiency%s",pid[ipid]),Form("Tracking efficiency BLC2;Layer;Event"),5,-0.5,4.5);
    new TH1F( Form("BPC_TrackingEfficiency%s",pid[ipid]),Form("Tracking efficiency BPC;Layer;Event"),5,-0.5,4.5);
    new TH1F( Form("FDC1_TrackingEfficiency%s",pid[ipid]),Form("Tracking efficiency FDC1;Layer;Event"),5,-0.5,4.5);
  }
  return true;
}
