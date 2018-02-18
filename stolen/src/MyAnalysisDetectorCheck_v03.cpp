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
}

bool MyAnalysisDetectorCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{

	rtFile->cd();
	DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  if(!header->IsTrig(Trig_Beam)){
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
  if(T0Timing<-1.0||1.0<T0Timing){
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
  /* ========== Event selection ========== */
  /* ===================================== */

  /* BLDC Tracking */
  bltrackMan->DoTracking(blMan,conf,true,true);

  /* Trigger pattern */
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }
  /* ================================ */
  /* ========== BLDC Check ========== */
    const int ndc = 6;
    int dccid[ndc] = {CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_BPC,CID_FDC1};
    /* BLC1 */
    bool blc1 = false;
    if(bltrackMan->ntrackBLC2()==1&&bltrackMan->ntrackBPC()==1){
      LocalTrack* blc2 = bltrackMan->trackBLC2(0);
      LocalTrack* bpc  = bltrackMan->trackBPC(0);
      double blc2chi2 = blc2->chi2all();
      double bpcchi2  = bpc->chi2all();
      if(20<blc2chi2||20<bpcchi2)
    }
    for(int i=0; i<bltrackMan->nBLC1(); i++){
      LocalTrack* track = bltrackMan->trackBLC1(i);
    }


    for(int idc=0; idc<=1; idc++){
      int cid1=dccid[idc];
      int cid2=dccid[1-idc];
      const char* name1=dlist->GetName(cid1).data();
      const char* name2=dlist->GetName(cid2).data();
      bool hitflag = false;
      if(bltrackMan->ntrackBLDC(cid2)!=1) { continue; }
      LocalTrack* track = bltrackMan->trackBLDC(cid2,0);
      double chi2 = track->chi2all();
      FillHist(Form("%s_ChiSquare",name2),chi2,1);
      for(int ihit=0; ihit<track->nhit(); ihit++){
        ChamberLikeHit* hit = track->hit(ihit);
        int layer = hit->layer();
        double resl = hit->resl();
        FillHist(Form("%sl%d_Residual",name2,layer),resl,1);
      }
      double xpos=-9999.9, ypos=-9999.9;
      double zpos=-15.0;
      if(idc==1) zpos*=-1.0;
      track -> XYLocalPosatZ(zpos,xpos,ypos);
      if((-10.0<xpos&&xpos<10.0)&&(-10.0<ypos&&ypos<10.0)&&(chi2<20)){
        hitflag = true;
      }
      if(!hitflag) { continue; }
      /* Hit efficiency */
      FillHist(Form("%s_HitEfficiency",name1),0);
      for(int layer=1; layer<=8; layer++){
        if(blMan->nBLDC(cid1,layer)>=1){
          FillHist(Form("%s_HitEfficiency",name1),layer);
        }
      }
      /* Tracking efficiency */
      FillHist(Form("%s_TrackingEfficiency",name1),0);
      if(bltrackMan->ntrackBLDC(cid1)==1){
        FillHist(Form("%s_TrackingEfficiency",name1),1);
      }
      if(bltrackMan->ntrackBLDC(cid1)>=1){
        FillHist(Form("%s_TrackingEfficiency",name1),2);
      }
    }

    /* ========== BLDC Check ========== */
    /* ================================ */

    /* =============================== */
    /* ========== PID Check ========== */
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
    int Beam = Beam_Other;
    if(25.0<tof&&tof<27.0){
      Beam = Beam_Pion;
    }
    if(27.5<tof&&tof<29.5){
      Beam = Beam_Kaon;
    }
    if(34.0<tof&&tof<35.0){
      Beam = Beam_Proton;
    }

    /* AC check */
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
      if(hit->CheckRange(1)){
        FillHist("AC_ADCSUM_withTDC",adcsum,1);
      }
      else{
        FillHist("AC_ADCSUM_withoutTDC",adcsum,1);
      }
    }
    /* AC if TOFPion */
    if(Beam==Beam_Pion){
      FillHist("AC_ADCSUM_ifTOFPion",adcsum,1);
    }

    if(achitflag){
      if(header->IsTrig(Trig_Pion)){
        FillHist("AC_ADCSUMTrigPion",adcsum,1);
      }
      if(header->IsTrig(Trig_Kaon)){
        FillHist("AC_ADCSUMTrigKaon",adcsum,1);
      }
      if(tofpiflag){
        FillHist("AC_ADCSUMwithTOFPion",adcsum,1);
        if(header->IsTrig(Trig_Pion)){
          FillHist("AC_ADCSUMwithTOFPionTrigPion",adcsum,1);
        }
        if(header->IsTrig(Trig_Kaon)){
          FillHist("AC_ADCSUMwithTOFPionTrigKaon",adcsum,1);
        }
      }
      if(tofkflag){
        FillHist("AC_ADCSUMwithTOFKaon",adcsum,1);
        if(header->IsTrig(Trig_Pion)){
          FillHist("AC_ADCSUMwithTOFKaonTrigPion",adcsum,1);
        }
        if(header->IsTrig(Trig_Kaon)){
          FillHist("AC_ADCSUMwithTOFKaonTrigKaon",adcsum,1);
        }
      }
    }
    /* ========== PID Check ========== */
    /* =============================== */


    /* =============================== */
    /* ========== BLDC Check ========== */
    if(tofpiflag){
      const int ndc = 6;
      int dccid[ndc] = {CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b,CID_BPC,CID_FDC1};
      /* BLC1 */
      for(int idc=0; idc<=1; idc++){
        int cid1=dccid[idc];
        int cid2=dccid[1-idc];
        const char* name1=dlist->GetName(cid1).data();
        const char* name2=dlist->GetName(cid2).data();
        bool hitflag = false;
        if(bltrackMan->ntrackBLDC(cid2)!=1) { continue; }
        LocalTrack* track = bltrackMan->trackBLDC(cid2,0);
        double chi2 = track->chi2all();
        FillHist(Form("%s_ChiSquare",name2),chi2,1);
        for(int ihit=0; ihit<track->nhit(); ihit++){
          ChamberLikeHit* hit = track->hit(ihit);
          int layer = hit->layer();
          double resl = hit->resl();
          FillHist(Form("%sl%d_Residual",name2,layer),resl,1);
        }
        double xpos=-9999.9, ypos=-9999.9;
        double zpos=-15.0;
        if(idc==1) zpos*=-1.0;
        track -> XYLocalPosatZ(zpos,xpos,ypos);
        if((-10.0<xpos&&xpos<10.0)&&(-10.0<ypos&&ypos<10.0)&&(chi2<20)){
          hitflag = true;
        }
        if(!hitflag) { continue; }
        /* Hit efficiency */
        FillHist(Form("%s_HitEfficiency",name1),0);
        for(int layer=1; layer<=8; layer++){
          if(blMan->nBLDC(cid1,layer)>=1){
            FillHist(Form("%s_HitEfficiency",name1),layer);
          }
        }
        /* Tracking efficiency */
        FillHist(Form("%s_TrackingEfficiency",name1),0);
        if(bltrackMan->ntrackBLDC(cid1)==1){
          FillHist(Form("%s_TrackingEfficiency",name1),1);
        }
        if(bltrackMan->ntrackBLDC(cid1)>=1){
          FillHist(Form("%s_TrackingEfficiency",name1),2);
        }
      }
      /* BLC2 */
      for(int idc=2; idc<=3; idc++){
        int cid1=dccid[idc];
        int cid2=dccid[5-idc];
        const char* name1=dlist->GetName(cid1).data();
        const char* name2=dlist->GetName(cid2).data();
        bool hitflag = false;
        if(bltrackMan->ntrackBLDC(cid2)!=1) { continue; }
        LocalTrack* track = bltrackMan->trackBLDC(cid2,0);
        double chi2 = track->chi2all();
        FillHist(Form("%s_ChiSquare",name2),chi2,1);
        for(int ihit=0; ihit<track->nhit(); ihit++){
          ChamberLikeHit* hit = track->hit(ihit);
          int layer = hit->layer();
          double resl = hit->resl();
          FillHist(Form("%sl%d_Residual",name2,layer),resl,1);
        }
        double xpos=-9999.9, ypos=-9999.9;
        double zpos=-15.0;
        if(idc==1) zpos*=-1.0;
        track -> XYLocalPosatZ(zpos,xpos,ypos);
        if((-10.0<xpos&&xpos<10.0)&&(-10.0<ypos&&ypos<10.0)&&(chi2<20)){
          hitflag = true;
        }
        if(!hitflag) { continue; }
        /* Hit efficiency */
        FillHist(Form("%s_HitEfficiency",name1),0);
        for(int layer=1; layer<=8; layer++){
          if(blMan->nBLDC(cid1,layer)>=1){
            FillHist(Form("%s_HitEfficiency",name1),layer);
          }
        }
        /* Tracking efficiency */
        FillHist(Form("%s_TrackingEfficiency",name1),0);
        if(bltrackMan->ntrackBLDC(cid1)==1){
          FillHist(Form("%s_TrackingEfficiency",name1),1);
        }
        if(bltrackMan->ntrackBLDC(cid1)>=1){
          FillHist(Form("%s_TrackingEfficiency",name1),2);
        }
      }
      /* BPC */
      if(bltrackMan->ntrackBLDC(CID_BLC2)==1){
        char* name1="BPC";
        bool hitflag1 = false;
        bool hitflag2 = false;
        LocalTrack* track = bltrackMan->trackBLC2(0);
        double chi2 = track->chi2all();
        TVector3 pos = track->GetPosatZ(-20.3);
        if((-0.5<pos.X()&&pos.X()<0.5)&&(-0.5<pos.Y()&&pos.Y()<0.5)&&(chi2<20)){
          hitflag1 = true;
        }
        int MulBVC = 0;
        for(int i=0; i<blMan->nBVC(); i++){
          if(blMan->BVC(i)->CheckRange()){
            hitflag2 = true;
          }
        }
        if(hitflag1&&hitflag2){
          /* Hit efficiency */
          FillHist(Form("%s_HitEfficiency",name1),0);
          for(int layer=1; layer<=8; layer++){
            if(blMan->nBLDC(CID_BPC,layer)>=1){
              FillHist(Form("%s_HitEfficiency",name1),layer);
            }
          }
        }
      }
      /* FDC */
      if(bltrackMan->ntrackBLDC(CID_BPC)>=1){
        char* name1="FDC1";
        bool hitflag1 = false;
        bool hitflag2 = false;
        LocalTrack* track = bltrackMan->trackBPC(0);
        double chi2 = track->chi2all();
        TVector3 pos = track->GetPosatZ(169.7);
        if((-2.0<pos.X()&&pos.X()<2.0)&&(-2.0<pos.Y()&&pos.Y()<2.0)&&(chi2<20)){
          hitflag1 = true;
        }
        int MulBVC = 0;
        for(int i=0; i<blMan->nBVC(); i++){
          if(blMan->BVC(i)->CheckRange()){
            hitflag2 = true;
          }
        }
        if(hitflag1&&hitflag2){
          /* Hit efficiency */
          FillHist(Form("%s_HitEfficiency",name1),0);
          for(int layer=1; layer<=6; layer++){
            if(blMan->nBLDC(CID_FDC1,layer)==1){
              FillHist(Form("%s_HitEfficiency",name1),layer);
            }
          }
        }
      }
    }

    /* Tracking efficiency of BLC1, BLC2, BPC and FDC1 */
    if(tofpiflag){
      if(bltrackMan->ntrackBLC2()==1&&bltrackMan->ntrackBPC()==1){
        FillHist("BLC1_TrackingEfficiency",0);
        if(bltrackMan->ntrackBLC1()>=1){
          FillHist("BLC1_TrackingEfficiency",1);
        }
      }
      if(bltrackMan->ntrackBLC1()==1&&bltrackMan->ntrackBPC()==1){
        FillHist("BLC2_TrackingEfficiency",0);
        if(bltrackMan->ntrackBLC2()>=1){
          FillHist("BLC2_TrackingEfficiency",1);
        }
      }
      if(bltrackMan->ntrackBLC1()==1&&bltrackMan->ntrackBLC2()==1){
        bool hitflag1 = false;
        LocalTrack* track = bltrackMan->trackBLC2(0);
        double chi2 = track->chi2all();
        TVector3 pos = track->GetPosatZ(-20.3);
        if((-0.5<pos.X()&&pos.X()<0.5)&&(-0.5<pos.Y()&&pos.Y()<0.5)&&(chi2<20)){
          hitflag1 = true;
        }
        if(hitflag1){
          FillHist("BPC_TrackingEfficiency",0);
          if(bltrackMan->ntrackBPC()>=1){
            FillHist("BPC_TrackingEfficiency",1);
            FillHist("FDC1_TrackingEfficiency",0);
            if(bltrackMan->ntrackFDC1()>=1){
              FillHist("FDC1_TrackingEfficiency",1);
            }
          }
        }
      }
    }




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

  /* TOF */
  new TH1F( Form("BHDT0_TOF"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);
  new TH1F( Form("BHDT0_TOFwithTrigBeam"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);
  new TH1F( Form("BHDT0_TOFwithTrigPion"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);
  new TH1F( Form("BHDT0_TOFwithTrigKaon"),Form("BHD-T0 Time of flight;BHD-T0 TOF (ns);Counts"),3000,20,50);

  /* AC */
  new TH1F( Form("AC_TDC"),Form("AC TDC;TDC (ch.);Counts"),4000,0,4000);
  new TH1F( Form("AC_TDCwithHit"),Form("AC TDC;TDC (ch.);Counts"),4000,0,4000);
  new TH1F( Form("AC_ADCSUM"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMTrigPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMwithTDC"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMwithoutTDC"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMwithTOFPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMwithTOFPionTrigPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMwithTOFPionTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMwithTOFKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMwithTOFKaonTrigPion"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);
  new TH1F( Form("AC_ADCSUMwithTOFKaonTrigKaon"),Form("AC ADC Sum.;ADC Sum. (ch.);Counts"),16000,0,16000);

  /* BLDC */
  const int ndc=6;
  int dccid[ndc]={CID_BLC1a,CID_BLC1b,
    CID_BLC2a,CID_BLC2b,
    CID_BPC,
    CID_FDC1
  };
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    int maxwire = 0;
    for( int layer=1; layer<=nlays; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire(cid,layer);
      if(maxwire<nwire) maxwire = nwire;
      new TH1F( Form("%sl%d_Residual",name,layer),Form("Residual %s-L%d;%s Residual (um);Counts",name,layer,name),4000,-0.2,0.2);
    }
    new TH1F( Form("%s_ChiSquare",name),Form("#Chi^{2}/NDF %s;#Chi^{2}/NDF;Counts",name),100,0,100);
    new TH1F( Form("%s_HitEfficiency",name),Form("Layer hit efficiency %s;Layer;Event",name),nlays+1,-0.5,nlays+0.5);
  }
  new TH1F( Form("BLC1_TrackingEfficiency"),Form("Tracking efficiency BLC1;Layer;Event"),3,-0.5,3.5);
  new TH1F( Form("BLC2_TrackingEfficiency"),Form("Tracking efficiency BLC2;Layer;Event"),3,-0.5,3.5);
  new TH1F( Form("BPC_TrackingEfficiency"),Form("Tracking efficiency BPC;Layer;Event"),3,-0.5,3.5);
  new TH1F( Form("FDC1_TrackingEfficiency"),Form("Tracking efficiency FDC1;Layer;Event"),3,-0.5,3.5);
  return true;
}
