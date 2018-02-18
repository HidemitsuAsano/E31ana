// MyAnalysisCDSCosmic.cpp

#include "MyAnalysisCDSCosmic.h"

MyAnalysisCDSCosmic::MyAnalysisCDSCosmic(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisCDSCosmic::~MyAnalysisCDSCosmic()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCDSCosmic::Clear()
{
}

bool MyAnalysisCDSCosmic::DoAnalysis(ConfMan* conf, EventHeader* header, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
{

  rtFile->cd();

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
  if(MulCDH!=2||MulIH!=2){
    return true;
  }
  FillHist("CDHvsIH_SegmentDifference",segdifcdh,segdifih,1);

  if(segdifcdh<5||segdifih<4){
    return true;
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
      FillHist(Form("CDHu%d_ADCwT",seg),au,1);
      FillHist(Form("CDHd%d_ADCwT",seg),ad,1);
      FillHist(Form("CDHu%d_dE",seg),eneu,1);
      FillHist(Form("CDHd%d_dE",seg),ened,1);
      FillHist(Form("CDH%d_dEMean",seg),emean,1);
      FillHist(Form("CDH%d_TimeMean",seg),tmean,1);
      FillHist(Form("CDH%d_TimeSub",seg),tsub,1);
      FillHist(Form("CDH%d_CTimeMean",seg),ctmean,1);
      FillHist(Form("CDH%d_CTimeSub",seg),ctsub,1);
      FillHist(Form("CDH%d_Position",seg),hitpos,1);
      FillHist(Form("CDHu%d_AvT",seg),au,tu,1);
      FillHist(Form("CDHd%d_AvT",seg),ad,td,1);
      FillHist(Form("CDHu%d_dEvTime",seg),eneu,timeu,1);
      FillHist(Form("CDHd%d_dEvTime",seg),ened,timed,1);
      FillHist(Form("CDHu%d_dEvCTime",seg),eneu,ctimeu,1);
      FillHist(Form("CDHd%d_dEvCTime",seg),ened,ctimed,1);
      FillHist(Form("CDH_HitPosition"),pos.X(),pos.Y(),1);
      FillHist("CDH_HitPat",seg,1);
      nCDH++;
    }
  }
  FillHist("CDH_Mult",nCDH,1);

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
      FillHist(Form("IHu%d_ADCwT",seg),au,1);
      FillHist(Form("IHd%d_ADCwT",seg),ad,1);
      FillHist(Form("IHu%d_dE",seg),eneu,1);
      FillHist(Form("IHd%d_dE",seg),ened,1);
      FillHist(Form("IH%d_dEMean",seg),emean,1);
      FillHist(Form("IH%d_TimeMean",seg),tmean,1);
      FillHist(Form("IH%d_TimeSub",seg),tsub,1);
      FillHist(Form("IH%d_CTimeMean",seg),ctmean,1);
      FillHist(Form("IH%d_CTimeSub",seg),ctsub,1);
      FillHist(Form("IH%d_Position",seg),hitpos,1);
      FillHist(Form("IHu%d_AvT",seg),au,tu,1);
      FillHist(Form("IHd%d_AvT",seg),ad,td,1);
      FillHist(Form("IHu%d_dEvTime",seg),eneu,timeu,1);
      FillHist(Form("IHd%d_dEvTime",seg),ened,timed,1);
      FillHist(Form("IHu%d_dEvCTime",seg),eneu,ctimeu,1);
      FillHist(Form("IHd%d_dEvCTime",seg),ened,ctimed,1);
      FillHist(Form("IH_HitPosition"),pos.X(),pos.Y(),1);
      FillHist("IH_HitPat",seg,1);
      nIH++;
    }
  }
  FillHist("IH_Mult",nIH,1);

  // Fill CDC
  for(int idc=0;idc<ndc;idc++){
    const int cid=CID_CDC;
    const char* name="CDC";
    const int nlays= dlist->GetNlayers(cid);
    for(int layer=1;layer<=nlays;layer++){
      int mult = cdsMan->nCDC(layer);
      FillHist(Form("%sl%d_Multiplicity",name,layer),mult);
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
  }
}

bool MyAnalysisCDSCosmic::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCDSCosmic::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCDSCosmic::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCDSCosmic::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisCDSCosmic::CutCondition()
{
}

double MyAnalysisCDSCosmic::CalcTimeOffs(HodoscopeLikeHit* hit, Particle* particle)
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

double MyAnalysisCDSCosmic::CalcTimeSubOffs(HodoscopeLikeHit* hit, Particle* particle)
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

bool MyAnalysisCDSCosmic::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisCDSCosmic::Initialize ###" << std::endl;

  confFile = confMan;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaCDSCosmic");

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
    new TH1F( Form("CDHu%d_ADCwT",seg), Form("ADC wTDC CDHU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH1F( Form("CDHd%d_ADCwT",seg), Form("ADC wTDC CDHD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH2F( Form("CDHu%d_AvT",seg),   Form("ADC TDC corr. CDHU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("CDHd%d_AvT",seg),   Form("ADC TDC corr. CDHD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("CDHu%d_dEvTime",seg),   Form("dE Time corr. CDHU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("CDHd%d_dEvTime",seg),   Form("dE Time corr. CDHD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("CDHu%d_dEvCTime",seg),   Form("dE CTime corr. CDHU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("CDHd%d_dEvCTime",seg),   Form("dE CTime corr. CDHD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
  }
  new TH1F( "CDH_HitPat", "Hit pattern CDH;CDH Segment;Counts", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
  new TH1F( "CDH_Mult", "Multiplicity CDH;Multiplicity;Counts", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
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
    new TH1F( Form("IHu%d_ADCwT",seg), Form("ADC wTDC IHU%d;ADC ch.;Counts",seg),  4000,    0, 4000 );
    new TH1F( Form("IHd%d_ADCwT",seg), Form("ADC wTDC IHD%d;ADC ch.;Counts",seg),  4000,    0, 4000 );
    new TH2F( Form("IHu%d_AvT",seg),   Form("ADC TDC corr. IHU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("IHd%d_AvT",seg),   Form("ADC TDC corr. IHD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("IHu%d_dEvTime",seg),   Form("dE Time corr. IHU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("IHd%d_dEvTime",seg),   Form("dE Time corr. IHD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("IHu%d_dEvCTime",seg),   Form("dE CTime corr. IHU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("IHd%d_dEvCTime",seg),   Form("dE CTime corr. IHD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
  }
  new TH1F( "IH_HitPat", "Hit pattern IH;IH Segment;Counts", NumOfIHSegments+1, 0, NumOfIHSegments+1 );
  new TH1F( "IH_Mult", "Multiplicity IH;Multiplicity;Counts", NumOfIHSegments+1, 0, NumOfIHSegments+1 );
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
  const int nlays= dlist->GetNlayers(cid);
  int maxwire = 0;
  std::cout << "Define Histograms for " << name << std::endl;
  for( int layer=1; layer<=nlays; layer++ ){
    int nwire = confMan->GetBLDCWireMapManager()->GetNWire(cid,layer);
    if(maxwire<nwire) maxwire = nwire;
    new TH1F( Form("%sl%d_Multiplicity_inTrack",name,layer),Form("Multiplicity %s-L%d;Multiplicity;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_HitPattern_inTrack",name,layer),Form("Hit pattern %s-L%d;Wire;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_Time_inTrack",name,layer),Form("Time %s-L%d;Time (ns);Counts",name,layer),400,0,400);
    new TH1F( Form("%sl%d_Residual_inTrack",name,layer),Form("Residual %s-L%d;Residual (um);Counts",name,layer),200,0,200);
    new TH2F( Form("%sl%d_TimevsResidual_inTrack",name,layer),Form("Time vs. Residual %s-L%d;Time (ns);Residual (um)",name,layer),400,0,400,200,0,200);

    new TH1F( Form("%sl%d_Multiplicity",name,layer),Form("Multiplicity %s-L%d;Multiplicity;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_HitPattern",name,layer),Form("Hit pattern %s-L%d;Wire;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
    new TH1F( Form("%sl%d_TDC",name,layer),Form("TDC %s-L%d;TDC ch.;Counts",name,layer),2000,0,2000);
    new TH1F( Form("%sl%d_Time",name,layer),Form("Time %s-L%d;Time (ns);Counts",name,layer),400,0,400);
    new TH1F( Form("%sl%d_Length",name,layer),Form("Length %s-L%d;Length (um);Counts",name,layer),200,0,200);
    new TH2F( Form("%sl%d_TimevsLength",name,layer),Form("Time vs. Length %s-L%d;Time (ns);Length (um)",name,layer),400,0,400,200,0,200);
    new TH1F( Form("%sl%d_Residual",name,layer),Form("Residual %s-L%d;Residual (um);Counts",name,layer),200,0,200);
    new TH2F( Form("%sl%d_TimevsResidual",name,layer),Form("Time vs. Residual %s-L%d;Time (ns);Residual (um)",name,layer),400,0,400,200,0,200);
    for(int wire=1;wire<=nwire;wire++){
      new TH1F( Form("%sl%dw%02d_TDC",name,layer,wire),Form("TDC %s-L%dW%02d;TDC ch.;Counts",name,layer,wire),2000,0,2000);
      new TH1F( Form("%sl%dw%02d_Time",name,layer,wire),Form("Time %s-L%dW%02d;Time (ns);Counts",name,layer,wire),400,0,400);
      new TH1F( Form("%sl%dw%02d_Length",name,layer,wire),Form("Length %s-L%dW%02d;Length (um);Counts",name,layer,wire),200,0,200);
      new TH2F( Form("%sl%dw%02d_TimevsLength",name,layer,wire),Form("Time vs. Length %s-L%dW%02d;Time (ns);Length (um)",name,layer,wire),400,0,400,200,0,200);
    }
  }
  for(int iconf=1;iconf<=LayConfigNum[idc];iconf++){
    new TH2F( Form("%s_HitPattern%c",name,LayConfig[idc][iconf-1]),Form("%s Hit pattern (%c%c' plane);Layer;Wire",name,LayConfig[idc][iconf-1],LayConfig[idc][iconf-1]),nlays+1,-0.5,nlays+0.5,maxwire+1,-0.5,maxwire+0.5);
  }
  new TH1F( Form("%s_HitEfficiency1%d",name,nlays),Form("Tracking efficiency %s (1-%d hit);Layer;Event",name,nlays),nlays+1,-0.5,nlays+0.5);
  new TH1F( Form("%s_HitEfficiency1%d",name,nlays-1),Form("Tracking efficiency %s (2-%d hit);Layer;Event",name,nlays-1),nlays+1,-0.5,nlays+0.5);
  new TH1F( Form("%s_HitEfficiency2%d",name,nlays),Form("Tracking efficiency %s (2-%d hit);Layer;Event",name,nlays),nlays+1,-0.5,nlays+0.5);


  return true;
}
