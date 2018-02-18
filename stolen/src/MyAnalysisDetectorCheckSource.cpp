// MyAnalysisDetectorCheckSource.cpp

#include "MyAnalysisDetectorCheckSource.h"

MyAnalysisDetectorCheckSource::MyAnalysisDetectorCheckSource(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisDetectorCheckSource::~MyAnalysisDetectorCheckSource()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisDetectorCheckSource::Clear()
{
  BeamPID = Beam_Other;
  conf = 0;
  header = 0; 
  blMan = 0;
  bltrackMan = 0;
}

bool MyAnalysisDetectorCheckSource::DoAnalysis(ConfMan* _conf, EventHeader* _header, BeamLineHitMan* _blMan, BeamLineTrackMan* _bltrackMan)
{
  rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();
  conf = _conf;
  header = _header;
  blMan = _blMan;
  bltrackMan = _bltrackMan;

  /* ======================================= */
  /* ========== Trigger selection ========== */
  //if(!header->IsTrig(Trig_Beam)){
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
  //double T0Timing = -9999.; 
  //int T0Segment = -1;
  //int MulT0 = 0;
  //for(int i=0; i<blMan->nT0(); i++){
  //  if(blMan->T0(i)->CheckRange()){
  //    MulT0++;
  //    T0Timing = blMan->T0(i)->ctmean();
  //    T0Segment = blMan->T0(i)->seg();
  //  }
  //}
  //FillHist("T0_Multiplicity",MulT0,1);
  //if(MulT0!=1){
  //  return true;
  //}
  //if(T0Timing<-10.0||10.0<T0Timing){
  //  return true;
  //}
  //double BHDTiming = -9999.; 
  //int BHDSegment = -1;
  //int MulBHD = 0;
  //for(int i=0; i<blMan->nBHD(); i++){
  //  if(blMan->BHD(i)->CheckRange()){
  //    MulBHD++;
  //    BHDTiming = blMan->BHD(i)->ctmean();
  //    BHDSegment = blMan->BHD(i)->seg();
  //  }
  //}
  //FillHist("BHD_Multiplicity",MulBHD,1);
  //if(MulBHD!=1){
  //  return true;
  //}
  /* ========== Event selection ========== */
  /* ===================================== */

  /* ================================ */
  /* ========== BLDC Check ========== */
  /* BLDC Tracking */
  bltrackMan->DoTracking(blMan,conf,true,true);

  char* pid = "";
  BLC1Check(pid);
  BLC2Check(pid);
  BPCCheck(pid);
  FDC1Check(pid);
  /* ========== BLDC Check ========== */
  /* ================================ */

  return true;

}

bool MyAnalysisDetectorCheckSource::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisDetectorCheckSource::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisDetectorCheckSource::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisDetectorCheckSource::FillHist(TString name, TString val1, TString val2, int weight)
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

double MyAnalysisDetectorCheckSource::TOFCheck(double T0Timing, double BHDTiming)
{
  return 0;
}

void MyAnalysisDetectorCheckSource::ACCheck()
{
  return;
}

void MyAnalysisDetectorCheckSource::BLC1Check(char* pid="")
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(blMan->nBLC1a(1)==0||blMan->nBLC1a(8)==0||blMan->nBLC1b(1)==0||blMan->nBLC1b(8)==0) return;
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

void MyAnalysisDetectorCheckSource::BLC2Check(char* pid="")
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(blMan->nBLC2a(1)==0||blMan->nBLC2a(8)==0||blMan->nBLC2b(1)==0||blMan->nBLC2b(8)==0) return;
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

void MyAnalysisDetectorCheckSource::BPCCheck(char* pid="")
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(blMan->nBPC(1)==0||blMan->nBPC(8)==0) return;
  FillHist(Form("BPC_TrackingEfficiency%s",pid),0,1);
  /* Basic data filling */
  const int cid=CID_BPC;
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

void MyAnalysisDetectorCheckSource::FDC1Check(char* pid="")
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(blMan->nFDC1(1)==0||blMan->nFDC1(8)==0) return;
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

void MyAnalysisDetectorCheckSource::CutCondition()
{
}

double MyAnalysisDetectorCheckSource::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

double MyAnalysisDetectorCheckSource::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
{
  return 0;
}

bool MyAnalysisDetectorCheckSource::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisDetectorCheckSource::Initialize ###" << std::endl;

  confFile = confMan;
  DetectorList *dlist=DetectorList::GetInstance();

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaDetectorCheckSource");

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
  const int npid=1;
  const char* pid[3] = {"","_withPion","_withKaon"};
  for(int idc=0;idc<2;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    for(int ipid=0; ipid<npid; ipid++){
      new TH1F( Form("%s_TrackTime%s",name,pid[ipid]),Form("Track time %s;%s Track time (ns);Counts",name,name),1400,-200,500);
      new TH1F( Form("%s_ChiSquare%s",name,pid[ipid]),Form("#Chi^{2}/NDF %s;%s #Chi^{2}/NDF;Counts",name,name),200,0,100);
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
      new TH1F( Form("%s_ChiSquare%s",name,pid[ipid]),Form("#Chi^{2}/NDF %s;%s #Chi^{2}/NDF;Counts",name,name),200,0,100);
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
