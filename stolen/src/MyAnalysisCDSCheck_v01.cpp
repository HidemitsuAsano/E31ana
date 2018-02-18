// MyAnalysisCDSCheck.cpp

#include "MyAnalysisCDSCheck.h"

MyAnalysisCDSCheck::MyAnalysisCDSCheck(TFile* rt, ConfMan* conf)
{
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

bool MyAnalysisCDSCheck::DoAnalysis(ConfMan* conf, EventHeader* header, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
{

  rtFile->cd();
	DetectorList *dlist=DetectorList::GetInstance();

  /* Event selection */
  //if(!header->IsTrig(Trig_Kaon)) return true;
  //else {
  //  return true;
  //}
  int MulCDH = 0;
  int SegCDH[2] = {0,0};
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()){
      if(MulCDH<2){
        SegCDH[MulCDH] = hit->seg();
      }
      MulCDH++;
    }
  }
  int segdifcdh = TMath::Abs(SegCDH[0]-SegCDH[1]);
  if(segdifcdh>18) segdifcdh -= 18;

  int MulIH = 0;
  int SegIH[2] = {0,0};
  for(int i=0; i<cdsMan->nIH(); i++){
    HodoscopeLikeHit* hit = cdsMan->IH(i);
    if(hit->CheckRange()){
      if(MulIH<2){
        SegIH[MulIH] = hit->seg();
      }
      MulIH++;
    }
  }
  int segdifih = TMath::Abs(SegIH[0]-SegIH[1]);
  if(segdifih>12) segdifcdh -= 12;

  /* Efficiency evaluation */
  if(MulCDH==2&&segdifcdh==18){
    FillHist("IH_HitEfficiency",0,2);
    HodoscopeLikeHit* cdh = cdsMan->CDH(0);
    double cdhazim = TMath::ATan(cdh->pos().Y()/cdh->pos().X());
    cdhazim *= 180.0/TMath::Pi();
    int nhit = 0;
    for(int i=0; i<cdsMan->nIH(); i++){
      HodoscopeLikeHit* hit = cdsMan->IH(i);
      if(hit->CheckRange()){
        TVector3 pos = hit->pos();
        double ihazim = TMath::ATan(pos.Y()/pos.X());
        ihazim *= 180.0/TMath::Pi();
        FillHist("CDHvsIH_AzimuthalAngle",cdhazim,ihazim,1);
        FillHist("CDHIH_AzimuthalAngleDifference",cdhazim-ihazim,1);
        if(TMath::Abs(cdhazim-ihazim<20)){
          nhit++;
          if(nhit>2){
            FillHist("IH_HitEfficiency",2,1);
          }
          else{
            FillHist("IH_HitEfficiency",1,1);
          }
        }
      }
    }
  }
  /* Efficiency evaluation */

  /* ===================================== */
  /* ========== Event selection ========== */
  if(MulCDH<1){
    return true;
  }
  FillHist("CDHvsIH_SegmentDifference",segdifcdh,segdifih,1);

  if(segdifcdh<5||segdifih<4){
    //return true;
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

  // Fill CDH
  int nCDH=0;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan->CDH(i);
    int seg = hit->seg();
    int au = hit->adcu(), ad = hit->adcd();
    int tu = hit->tdcu(), td = hit->tdcd();
    double hitpos = hit->hitpos();
    TVector3 pos = hit->pos();
    pos.SetY(hitpos);
    double timeu = hit->time(0), timed = hit->time(1);
    double ctimeu = hit->ctime(0), ctimed = hit->ctime(1);
    double eneu = hit->ene(0), ened = hit->ene(1);
    double emean = hit->emean();
    double tmean = hit->tmean();
    double ctmean = hit->ctmean();
    double tsub = hit->tsub();
    double ctsub = hit->ctsub();

    FillHist(Form("CDHu%d_ADC",seg),au,1);
    FillHist(Form("CDHd%d_ADC",seg),ad,1);
    if(tu>0){
      FillHist(Form("CDHu%d_TDC",seg),tu,1);
      FillHist(Form("CDHu%d_Time",seg),timeu,1);
      FillHist(Form("CDHu%d_CTime",seg),ctimeu,1);
    }
    if(td>0){
      FillHist(Form("CDHd%d_TDC",seg),td,1);
      FillHist(Form("CDHd%d_Time",seg),timed,1);
      FillHist(Form("CDHd%d_CTime",seg),ctimed,1);
    }
    if( hit->CheckRange() ){
      FillHist(Form("CDHu%d_ADCwithTDC",seg),au,1);
      FillHist(Form("CDHd%d_ADCwithTDC",seg),ad,1);
      FillHist(Form("CDHu%d_dE",seg),eneu,1);
      FillHist(Form("CDHd%d_dE",seg),ened,1);
      FillHist(Form("CDH%d_dEMean",seg),emean,1);
      FillHist(Form("CDH%d_TimeMean",seg),tmean,1);
      FillHist(Form("CDH%d_TimeSub",seg),tsub,1);
      FillHist(Form("CDH%d_CTimeMean",seg),ctmean,1);
      FillHist(Form("CDH%d_CTimeSub",seg),ctsub,1);
      FillHist(Form("CDH%d_Position",seg),hitpos,1);
      FillHist(Form("CDHu%d_ADCvsTDC",seg),au,tu,1);
      FillHist(Form("CDHd%d_ADCvsTDC",seg),ad,td,1);
      FillHist(Form("CDHu%d_dEvsTime",seg),eneu,timeu,1);
      FillHist(Form("CDHd%d_dEvsTime",seg),ened,timed,1);
      FillHist(Form("CDHu%d_dEvsCTime",seg),eneu,ctimeu,1);
      FillHist(Form("CDHd%d_dEvsCTime",seg),ened,ctimed,1);
      FillHist(Form("CDH_HitPosition"),pos.X(),pos.Y(),1);
      FillHist("CDH_HitPattern",seg,1);
      nCDH++;
    }
  }
  FillHist("CDH_Multiplicity",nCDH,1);

  // Fill IH
  int nIH=0;
  for( int i=0; i<cdsMan->nIH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan->IH(i);
    int seg = hit->seg();
    int au = hit->adcu(), ad = hit->adcd();
    int tu = hit->tdcu(), td = hit->tdcd();
    double hitpos = hit->hitpos();
    TVector3 pos = hit->pos();
    pos.SetY(hitpos);
    double timeu = hit->time(0), timed = hit->time(1);
    double ctimeu = hit->ctime(0), ctimed = hit->ctime(1);
    double eneu = hit->ene(0), ened = hit->ene(1);
    double emean = hit->emean();
    double tmean = hit->tmean();
    double ctmean = hit->ctmean();
    double tsub = hit->tsub();
    double ctsub = hit->ctsub();

    FillHist(Form("IHu%d_ADC",seg),au,1);
    FillHist(Form("IHd%d_ADC",seg),ad,1);
    if(tu>0){
      FillHist(Form("IHu%d_TDC",seg),tu,1);
      FillHist(Form("IHu%d_Time",seg),timeu,1);
      FillHist(Form("IHu%d_CTime",seg),ctimeu,1);
    }
    if(td>0){
      FillHist(Form("IHd%d_TDC",seg),td,1);
      FillHist(Form("IHd%d_Time",seg),timed,1);
      FillHist(Form("IHd%d_CTime",seg),ctimed,1);
    }
    if( hit->CheckRange() ){
      FillHist(Form("IHu%d_ADCwithTDC",seg),au,1);
      FillHist(Form("IHd%d_ADCwithTDC",seg),ad,1);
      FillHist(Form("IHu%d_dE",seg),eneu,1);
      FillHist(Form("IHd%d_dE",seg),ened,1);
      FillHist(Form("IH%d_dEMean",seg),emean,1);
      FillHist(Form("IH%d_TimeMean",seg),tmean,1);
      FillHist(Form("IH%d_TimeSub",seg),tsub,1);
      FillHist(Form("IH%d_CTimeMean",seg),ctmean,1);
      FillHist(Form("IH%d_CTimeSub",seg),ctsub,1);
      FillHist(Form("IH%d_Position",seg),hitpos,1);
      FillHist(Form("IHu%d_ADCvsTDC",seg),au,tu,1);
      FillHist(Form("IHd%d_ADCvsTDC",seg),ad,td,1);
      FillHist(Form("IHu%d_dEvsTime",seg),eneu,timeu,1);
      FillHist(Form("IHd%d_dEvsTime",seg),ened,timed,1);
      FillHist(Form("IHu%d_dEvsCTime",seg),eneu,ctimeu,1);
      FillHist(Form("IHd%d_dEvsCTime",seg),ened,ctimed,1);
      FillHist(Form("IH_HitPosition"),pos.X(),pos.Y(),1);
      FillHist("IH_HitPattern",seg,1);
      nIH++;
    }
  }
  FillHist("IH_Multiplicity",nIH,1);

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
      double dl = hit->dl();
      FillHist(Form("%sl%d_HitPattern",name,layer),wire,1);
      FillHist(Form("%sl%d_TDC",name,layer),tdc,1);
      FillHist(Form("%sl%d_Time",name,layer),dt,1);
      FillHist(Form("%sl%d_Length",name,layer),dl,1);
      FillHist(Form("%sl%d_TimevsLength",name,layer),dt,dl,1);
      FillHist(Form("%sl%dw%02d_TDC",name,layer,wire),tdc,1);
      FillHist(Form("%sl%dw%02d_Time",name,layer,wire),dt,1);
      FillHist(Form("%sl%dw%02d_Length",name,layer,wire),dl,1);
      FillHist(Form("%sl%dw%02d_TimevsLength",name,layer,wire),dt,dl,1);
    }
  }

  if(cdstrackMan==0) return true;
  // For tracking
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
        double dl = hit->dl();
        double resl = hit->resl()*10000;
        FillHist(Form("%sl%d_HitPatterninTrack",name,layer),wire,1);
        FillHist(Form("%sl%d_TimeinTrack",name,layer),dt,1);
        FillHist(Form("%sl%d_LengthinTrack",name,layer),dl,1);
        FillHist(Form("%sl%d_ResidualinTrack",name,layer),resl,1);
        FillHist(Form("%sl%d_TimevsLengthinTrack",name,layer),dt,dl,1);
        FillHist(Form("%sl%d_TimevsResidualinTrack",name,layer),dt,resl,1);
        FillHist(Form("%sl%dw%02d_TimeinTrack",name,layer,wire),dt,1);
        FillHist(Form("%sl%dw%02d_LengthinTrack",name,layer,wire),dl,1);
        FillHist(Form("%sl%dw%02d_ResidualinTrack",name,layer,wire),resl,1);
        FillHist(Form("%sl%dw%02d_TimevsLengthinTrack",name,layer,wire),dt,dl,1);
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

  /* CDH */ 
  std::cout << "Define Histograms for CDH" << std::endl;
  int NumOfCDHSegments = 36;
  for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
    new TH1F( Form("CDHu%d_ADC",seg),   Form("ADC CDHU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("CDHd%d_ADC",seg),   Form("ADC CDHD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("CDHu%d_dE",seg),   Form("dE CDHU%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("CDHd%d_dE",seg),   Form("dE CDHD%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("CDH%d_dEMean",seg),   Form("Mean dE CDH%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("CDHu%d_TDC",seg),   Form("TDC CDHU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("CDHd%d_TDC",seg),   Form("TDC CDHD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("CDHu%d_Time",seg),   Form("Time CDHU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("CDHd%d_Time",seg),   Form("Time CDHD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("CDH%d_TimeMean",seg),   Form("Mean Time CDH%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("CDH%d_TimeSub",seg),   Form("Subtract Time CDH%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("CDHu%d_CTime",seg),   Form("CTime CDHU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("CDHd%d_CTime",seg),   Form("CTime CDHD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("CDH%d_CTimeMean",seg),   Form("Mean CTime CDH%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("CDH%d_CTimeSub",seg),   Form("Subtract CTime CDH%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("CDH%d_Position",seg),   Form("Hit position CDH%d;Position (cm);Counts",seg),    1000,    -50, 50 );
    new TH1F( Form("CDHu%d_ADCwithTDC",seg), Form("ADC wTDC CDHU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH1F( Form("CDHd%d_ADCwithTDC",seg), Form("ADC wTDC CDHD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH2F( Form("CDHu%d_ADCvsTDC",seg),   Form("ADC TDC corr. CDHU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("CDHd%d_ADCvsTDC",seg),   Form("ADC TDC corr. CDHD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("CDHu%d_dEvsTime",seg),   Form("dE Time corr. CDHU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("CDHd%d_dEvsTime",seg),   Form("dE Time corr. CDHD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("CDHu%d_dEvsCTime",seg),   Form("dE CTime corr. CDHU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("CDHd%d_dEvsCTime",seg),   Form("dE CTime corr. CDHD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
  }
  new TH1F( "CDH_HitPattern", "Hit pattern CDH;CDH Segment;Counts", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
  new TH1F( "CDH_Multiplicity", "Multiplicity CDH;Multiplicity;Counts", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
  new TH2F( Form("CDH_HitPosition"),   Form("Position CDH;X hit position (cm);Y hit position (cm)"), 80, -20, 20, 80, -20, 20 );

  /* IH */ 
  std::cout << "Define Histograms for IH" << std::endl;
  int NumOfIHSegments = 24;
  for( int seg=1; seg<=NumOfIHSegments; seg++ ){
    new TH1F( Form("IHu%d_ADC",seg),   Form("ADC IHU%d;ADC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("IHd%d_ADC",seg),   Form("ADC IHD%d;ADC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("IHu%d_dE",seg),   Form("dE IHU%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("IHd%d_dE",seg),   Form("dE IHD%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("IH%d_dEMean",seg),   Form("Mean dE IH%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("IHu%d_TDC",seg),   Form("TDC IHU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("IHd%d_TDC",seg),   Form("TDC IHD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("IHu%d_Time",seg),   Form("Time IHU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("IHd%d_Time",seg),   Form("Time IHD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("IH%d_TimeMean",seg),   Form("Mean Time IH%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("IH%d_TimeSub",seg),   Form("Subtract Time IH%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("IHu%d_CTime",seg),   Form("CTime IHU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("IHd%d_CTime",seg),   Form("CTime IHD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("IH%d_CTimeMean",seg),   Form("Mean CTime IH%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("IH%d_CTimeSub",seg),   Form("Subtract CTime IH%d;Subtract Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("IH%d_Position",seg),   Form("Hit position IH%d;Position (cm);Counts",seg),    1000,    -50, 50 );
    new TH1F( Form("IHu%d_ADCwithTDC",seg), Form("ADC wTDC IHU%d;ADC ch.;Counts",seg),  4000,    0, 4000 );
    new TH1F( Form("IHd%d_ADCwtihTDC",seg), Form("ADC wTDC IHD%d;ADC ch.;Counts",seg),  4000,    0, 4000 );
    new TH2F( Form("IHu%d_ADCvsTDC",seg),   Form("ADC TDC corr. IHU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("IHd%d_ADCvsTDC",seg),   Form("ADC TDC corr. IHD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("IHu%d_dEvsTime",seg),   Form("dE Time corr. IHU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("IHd%d_dEvsTime",seg),   Form("dE Time corr. IHD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("IHu%d_dEvsCTime",seg),   Form("dE CTime corr. IHU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("IHd%d_dEvsCTime",seg),   Form("dE CTime corr. IHD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
  }
  new TH1F( "IH_HitPattern", "Hit pattern IH;IH Segment;Counts", NumOfIHSegments+1, 0, NumOfIHSegments+1 );
  new TH1F( "IH_Multiplicity", "Multiplicity IH;Multiplicity;Counts", NumOfIHSegments+1, 0, NumOfIHSegments+1 );
  new TH2F( Form("IH_HitPosition"),   Form("Position IH;X hit position (cm);Y hit position (cm)"), 80, -20, 20, 80, -20, 20 );
  h1 = new TH1F( "IH_HitEfficiency", "Efficiency IH;;Counts", 3, -0.5, 2.5 );
  h1->GetXaxis()->SetBinLabel(1,"All");
  h1->GetXaxis()->SetBinLabel(2,"Mult<=2");
  h1->GetXaxis()->SetBinLabel(3,"Mult>2");

  new TH2F( Form("CDHvsIH_SegmentDifference"),   Form("Segment difference CDH vs. IH;CDH segment difference;IH segment differecne"), 18, 0.5, 18.5, 12, 0.5, 12.5 );
  new TH2F( Form("CDHvsIH_AzimuthalAngle"),   Form("Azimuthal angle CDH vs. IH;CDH azimuthal angle;IH azimuthal angle"), 18, -90, 90, 12, -90, 90 );
  new TH1F( Form("CDHIH_AzimuthalAngleDifference"),   Form("Azimuthal angle difference CDH - IH;CDH azimuthal angle difference (CDH-IH);Counts"), 36, -90, 90 );

  /* CDC */
  const int cid= CID_CDC;
  const char* name= dlist->GetName(cid).data();
  const int nlays= NumOfCDCLayers;
  std::cout << "Define Histograms for " << name << std::endl;
  new TH1F( Form("%s_ChiSquare",name), Form("ChiSquare %s;Chi^{2}/ndf;Counts",name), 400, 0, 100 );
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
    new TH2F( Form("%sl%d_TimevsLength",name,layer),Form("Time vs. Length %s-L%d;%s (layer %d) Time (ns);%s (layer %d) Length (um)",name,layer,name,layer,name,layer),800,-400,400,200,0,200);

    new TH1F( Form("%sl%d_MultiplicityinTrack",name,layer),Form("Multiplicity %s-L%d;%s (layer %d) Multiplicity;Counts",name,layer,name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_HitPatterninTrack",name,layer),Form("Hit pattern %s-L%d;%s (layer %d) Wire;Counts",name,layer,name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_TimeinTrack",name,layer),Form("Time %s-L%d;%s (layer %d) Time (ns);Counts",name,layer,name,layer),800,-400,400);
    new TH1F( Form("%sl%d_ResidualinTrack",name,layer),Form("Residual %s-L%d;%s (layer %d) Residual (um);Counts",name,layer,name,layer),2000,-2000,2000);
    new TH2F( Form("%sl%d_TimevsResidualinTrack",name,layer),Form("Time vs. Residual %s-L%d;%s (layer %d) Time (ns);Residual (um)",name,layer,name,layer),800,-400,400,2000,-2000,2000);

    for(int wire=1;wire<=nwire;wire++){
      new TH1F( Form("%sl%dw%02d_TDC",name,layer,wire),Form("TDC %s-L%dW%02d;%s (layer %d,wire %d) TDC ch.;Counts",name,layer,wire,name,layer,wire),2000,0,2000);
      new TH1F( Form("%sl%dw%02d_Time",name,layer,wire),Form("Time %s-L%dW%02d;%s (layer %d,wire %d) Time (ns);Counts",name,layer,wire,name,layer,wire),800,-400,400);
      new TH1F( Form("%sl%dw%02d_Length",name,layer,wire),Form("Length %s-L%dW%02d;%s (layer %d,wire %d) Length (um);Counts",name,layer,wire,name,layer,wire),200,0,200);

      new TH1F( Form("%sl%dw%02d_TimeinTrack",name,layer,wire),Form("Time %s-L%dW%02d;%s (layer %d,wire %d) Time (ns);Counts",name,layer,wire,name,layer,wire),800,-400,400);
      new TH1F( Form("%sl%dw%02d_ResidualinTrack",name,layer,wire),Form("Residual %s-L%dW%02d;%s (layer %d,wire %d) Residual (um);Counts",name,layer,wire,name,layer,wire),1000,-2000,2000);
      new TH2F( Form("%sl%dw%02d_TimevsResidualinTrack",name,layer,wire),Form("Time vs. Length %s-L%dW%02d;%s (layer %d,wire %d) Time (ns);%s (layer %d,wire %d) Length (um)",name,layer,wire,name,layer,wire,name,layer,wire),800,-400,400,1000,-2000,2000);
    }
  }

  return true;
}
