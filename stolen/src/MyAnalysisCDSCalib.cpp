// MyAnalysisCDSCalib.cpp

#include "MyAnalysisCDSCalib.h"

MyAnalysisCDSCalib::MyAnalysisCDSCalib(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisCDSCalib::~MyAnalysisCDSCalib()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCDSCalib::Clear()
{
}

bool MyAnalysisCDSCalib::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double T0Timing = beam->t0time();
  int MulT0=0;
  for(int i=0; i<blMan->nT0(); i++){
    HodoscopeLikeHit* hit = blMan->T0(i);
    if(hit->CheckRange()){
      MulT0++;
      T0Timing = hit->ctmean();
    }
  }
  if(MulT0!=1) return true;
  beam->SetT0Time(T0Timing);

  /* Trigger selection */
  //if(!header->trigmode(Mode_KCDH2)) return true;
  //if(header->IsTrig(Trig_KCDH2)&&header->IsTrig(Trig_1stMix)){
  //if(header->IsTrig(Trig_KCDH2)&&header->IsTrig(Trig_1stMix)){
  //}
  //else {
  //  return true;
  //}
  if(header->IsTrig(Trig_Cosmic)) return true;

  /* Event selection */
  int MulCDH=0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()) MulCDH++;
  }
  if(MulCDH==0) return true;
  if(cdstrackMan->nGoodTrack()==0) return true;

  //if(particle->nCDS()==0) return true;
  //if(particle->nPiplus()==0&&particle->nPiminus()==0) return true;

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
      FillHist(Form("Trigger%02d_TDC",i),val);
    }
  }

  int CDHMult[36];
  for(int seg=1; seg<=36; seg++){
    CDHMult[seg-1]=0;
  }

  for(int id=0; id<cdstrackMan->nGoodTrack(); id++){
    CDSTrack* cdc=cdstrackMan->GoodTrack(id);
    if(cdc->nCDHHit()!=1) continue;
    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    if(cdc->Chi()>30) continue;

    TVector3 cdsvtx,beamvtx;
    double par[5];
    cdc->GetParameters(CID_CDC,par,cdsvtx);
    double mom=cdc->Momentum();
    if(!cdc->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
    dis = (beamvtx-cdsvtx).Mag();

    double beam_out=0, beam_tof=0;
    ELossTools::CalcElossBeamTGeo(beam->t0pos(),beamvtx,beam->mom(),kpMass,beam_out,beam_tof);
    TVector3 cdhvtx=cdc->CDHVertex();
    double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
    double beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
    double mass2=mom*mom*(1/(beta*beta)-1);

    double calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;

    int pid=TrackTools::PID(mom,mass2);
    pid=TrackTools::PIDcorr3(mom,mass2);
    cdc->SetPID(pid);

    TVector3 cdsvtx2,beamvtx2;
    double tof=0,tmpl=0;
    if(!cdc->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],cdsvtx2,beamvtx2,tof,tmpl,true)) continue;
    TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
    TVector3 vtxb = beamvtx2;
    double vbdis = (beamvtx2-cdsvtx2).Mag();

    TVector3 cfrpvtx;
    double param[5];  
    cdc->GetParameters(CID_CDCCFRP,param,cfrpvtx);
    cdc->GetParameters(param);
    TVector3 vtxcdh = cdc->CDHVertex();
    cdc->GetParameters(param);
    double fl = MathTools::CalcHelixArc(param,cfrpvtx,vtxcdh);
    double mom2 = TMath::Abs(cdc->Momentum(CID_CDCCFRP));
    double mass = cdsMass[pid];
    double momout=0,tmptof=0;
    ELossTools::CalcdE(mom2,mass,fl,"CDCGas",momout,-1,tmptof);
    double offset = cdhtime - beam->t0time() - beam_tof - tof - tmptof;

    //if(GeomTools::GetID(vtxb)!=CID_Fiducial){ continue; }

    HodoscopeLikeHit* cdh = cdc->CDHHit(cdsMan,0);
    CDHMult[cdh->seg()-1]++;
  }
  for(int seg=1; seg<=36; seg++){
    FillHist(Form("CDH%d_Multiplicity",seg),CDHMult[seg],1);
  }

  for(int id=0; id<cdstrackMan->nGoodTrack(); id++){
    CDSTrack* cdc=cdstrackMan->GoodTrack(id);
    if(cdc->nCDHHit()!=1) continue;
    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    if(cdc->Chi()>30) continue;

    TVector3 cdsvtx,beamvtx;
    double par[5];
    cdc->GetParameters(CID_CDC,par,cdsvtx);
    double mom=cdc->Momentum();
    if(!cdc->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
    dis = (beamvtx-cdsvtx).Mag();

    double beam_out=0, beam_tof=0;
    ELossTools::CalcElossBeamTGeo(beam->t0pos(),beamvtx,beam->mom(),kpMass,beam_out,beam_tof);
    TVector3 cdhvtx=cdc->CDHVertex();
    double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
    double beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
    double mass2=mom*mom*(1/(beta*beta)-1);
    FillHist("CDS_OverbetavsMomentum",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum",mass2,mom);

    double calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;
    FillHist("CDS_OverbetavsMomentum_AfterFindMass",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum_AfterFindMass",mass2,mom);

    int pid=TrackTools::PID(mom,mass2);
    pid=TrackTools::PIDcorr3(mom,mass2);
    cdc->SetPID(pid);

    TVector3 cdsvtx2,beamvtx2;
    double tof=0,tmpl=0;
    if(!cdc->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],cdsvtx2,beamvtx2,tof,tmpl,true)) continue;
    FillHist("CDS_OverbetavsMomentum_AfterELossCorr",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum_AfterELossCorr",mass2,mom);
    TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
    FillHist(Form("CDS_OverbetavsMomentum_%s",pname[pid].Data()),1./beta,mom);
    FillHist(Form("CDS_Mass2vsMomentum_%s",pname[pid].Data()),mass2,mom);
    TVector3 vtxb = beamvtx2;
    double vbdis = (beamvtx2-cdsvtx2).Mag();

    TVector3 cfrpvtx;
    double param[5];  
    cdc->GetParameters(CID_CDCCFRP,param,cfrpvtx);
    cdc->GetParameters(param);
    TVector3 vtxcdh = cdc->CDHVertex();
    cdc->GetParameters(param);
    double fl = MathTools::CalcHelixArc(param,cfrpvtx,vtxcdh);
    double mom2 = TMath::Abs(cdc->Momentum(CID_CDCCFRP));
    double mass = cdsMass[pid];
    double momout=0,tmptof=0;
    ELossTools::CalcdE(mom2,mass,fl,"CDCGas",momout,-1,tmptof);
    double offset = cdhtime - beam->t0time() - beam_tof - tof - tmptof;

    //if(GeomTools::GetID(vtxb)!=CID_Fiducial){ continue; }

    HodoscopeLikeHit* cdh = cdc->CDHHit(cdsMan,0);
    if(CDHMult[cdh->seg()-1]!=1){ continue; }
    FillHist("CDS_OverbetavsEnergy",1./beta,cdh->emean());
    if(pid==CDS_Other){ continue; }
    FillHist(Form("CDH_TimeOffset"),offset,1);
    FillHist(Form("CDH_TimeOffsetvsMomentum"),offset,TMath::Abs(mom),1);
    FillHist(Form("CDH_dEMeanvsTimeOffset"),cdh->emean(),offset,1);
    FillHist(Form("CDH%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDH%d_TimeOffsetvsMomentum",cdh->seg()),offset,TMath::Abs(mom),1);
    FillHist(Form("CDH%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);
    FillHist(Form("CDH%d_CTimeSub",cdh->seg()),cdh->ctsub(),1);
    FillHist(Form("CDHu%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDHu%d_dEvsTimeOffset",cdh->seg()),cdh->eu(),offset,1);
    FillHist(Form("CDHu%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);
    FillHist(Form("CDHd%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDHd%d_dEvsTimeOffset",cdh->seg()),cdh->ed(),offset,1);
    FillHist(Form("CDHd%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);
    FillHist(Form("CDH%d_TimeOffset_%s",cdh->seg(),pname[pid].Data()),offset,1);
    FillHist(Form("CDH%d_TimeOffsetvsMomentum_%s",cdh->seg(),pname[pid].Data()),offset,TMath::Abs(mom),1);
    FillHist(Form("CDS_Mass2vsdEMean",pname[pid].Data()),mass2,cdh->emean());

    for(int ip=0; ip<8; ip++){
      if(pid==ip){
        double massdiff = mass2 - cdsMass[pid]*cdsMass[pid];
        double betagamma = TMath::Abs(mom)/sqrt(mass2);
        FillHist(Form("CDS_Mass2Difference_%s",pname[pid].Data()),massdiff);
        FillHist(Form("CDS_BetaGammavsMass2Difference_%s",pname[pid].Data()),betagamma,massdiff);
      }
    }
    if(vbdis>10.0){
      FillHist(Form("Vertex_XY_VBDISCut"), vtxb.X(), vtxb.Y());
      FillHist(Form("CDC_Chi2_VBDISCut"), cdc->Chi());
    }
    if(offset>0.5){
      double zcutu = -20.0+2.0;
      int icut = 0;
      while(zcutu<20.1){
        icut++;
        double zcutl = zcutu-2.0;
        if(zcutl<vtxb.Z()&&vtxb.Z()<zcutu){
          int cutnum = 1;
          FillHist(Form("Vertex_XY_TimeOffsetCut_ZCut%d",icut), vtxb.X(), vtxb.Y());
        }
        zcutu=zcutu+2.0;
      }
      if(vtxb.Z()<-20.0){
        FillHist(Form("Vertex_XY_TimeOffsetCut_ZCut0",icut), vtxb.X(), vtxb.Y());
      }
      if(vtxb.Z()>20.0){
        FillHist(Form("Vertex_XY_TimeOffsetCut_ZCut21",icut), vtxb.X(), vtxb.Y());
      }
    }
    if(offset<-2.0){
      int cutnum = 1;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(-2.0<offset&&offset<-1.5){
      int cutnum = 2;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(-1.5<offset&&offset<-1.0){
      int cutnum = 3;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(-1.0<offset&&offset<-0.5){
      int cutnum = 4;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(-0.5<offset&&offset<-0.0){
      int cutnum = 5;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(0.0<offset&&offset<0.5){
      int cutnum = 6;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(0.5<offset&&offset<1.0){
      int cutnum = 7;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(1.0<offset&&offset<1.5){
      int cutnum = 8;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(1.5<offset&&offset<2.0){
      int cutnum = 9;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    if(2.0<offset){
      int cutnum = 10;
      FillHist(Form("Vertex_XY_TimeOffsetCut%d",cutnum), vtxb.X(), vtxb.Y());
      FillHist(Form("VBDIS_TimeOffsetCut%d",cutnum), vbdis);
      FillHist(Form("CDC_Chi2_TimeOffsetCut%d",cutnum), cdc->Chi());
      //FillHist(Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
      //FillHist("Vertex_ZX_TimeOffsetCut", vtxb.Z(), vtxb.X());
      //FillHist("Vertex_ZY_TimeOffsetCut", vtxb.Z(), vtxb.Y());
    }
    //else{
    //  FillHist(Form("CDS_Mass2vsMomentum_woTimeOffsetCut_%s",pname[pid].Data()),mass2,mom);
    //  FillHist("Vertex_XY_woTimeOffsetCut", vtxb.X(), vtxb.Y());
    //  //FillHist("Vertex_ZX_woTimeOffsetCut", vtxb.Z(), vtxb.X());
    //  //FillHist("Vertex_ZY_woTimeOffsetCut", vtxb.Z(), vtxb.Y());
    //  FillHist("VBDIS_woTimeOffsetCut", vbdis);
    //}
  }

  return true;
}

bool MyAnalysisCDSCalib::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCDSCalib::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCDSCalib::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCDSCalib::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisCDSCalib::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4976;
  double sigma = 0.0069;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Missing neutron mass cut */ 
  mean  = 0.9379;
  sigma = 0.0230;
  mnll = mean-2*sigma; mnul = mean+2*sigma;
  sblmnll = mean-5*sigma; sblmnul = mean-3*sigma;
  sbumnll = mean+3*sigma; sbumnul = mean+5*sigma;
  /* SigmaPlus mass cut */
  mean  = 1.1888;
  sigma = 0.0044;
  spll = mean-2*sigma; spul = mean+2*sigma;
  sblspll = mean-5*sigma; sblspul = mean-3*sigma;
  sbuspll = mean+3*sigma; sbuspul = mean+5*sigma;
  /* SigmaMinus mass cut */
  mean  = 1.1970;
  sigma = 0.0054;
  smll = mean-2*sigma; smul = mean+2*sigma;
  sblsmll = mean-5*sigma; sblsmul = mean-3*sigma;
  sbusmll = mean+3*sigma; sbusmul = mean+5*sigma;
  /* Missing k0 mass cut */ 
  mean  = 0.4970;
  sigma = 0.0130;
  mk0ll = mean-2*sigma; mk0ul = mean+2*sigma;
  sblmk0ll = mean-5*sigma; sblmk0ul = mean-3*sigma;
  sbumk0ll = mean+3*sigma; sbumk0ul = mean+5*sigma;
}

bool MyAnalysisCDSCalib::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisCDSCalib::Initialize ###" << std::endl;
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

  DetectorList *dlist=DetectorList::GetInstance();
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaCDSCalib");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
  new TH1F( "EventNumber", "Number of Events", 10, 0, 10);
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
  for(int i=1; i<=18; i++){
    new TH1F( Form("Trigger%02d_TDC",i), Form("Trigger%02d TDC;TDC ch.;Counts",i), 4096, 0, 4096);
  }

  new TH1F( "CDS_NumOfParticle", "Number of CDS tracks", 10, 0, 10);
  TH2F* h2;
  h2 = new TH2F("CDS_Particle","Detected particles by CDS;Particle;Number of all tracks",6,0,6,6,0,6);
  h2->GetYaxis()->SetBinLabel(1,"1 track");
  h2->GetYaxis()->SetBinLabel(2,"2 track");
  h2->GetYaxis()->SetBinLabel(3,"3 track");
  h2->GetYaxis()->SetBinLabel(4,"4 track");
  h2->GetYaxis()->SetBinLabel(5,"5 track");
  h2->GetYaxis()->SetBinLabel(6,"6 track");
  h2->GetXaxis()->SetBinLabel(1,"#pi^{+}");
  h2->GetXaxis()->SetBinLabel(2,"#pi^{-}");
  h2->GetXaxis()->SetBinLabel(3,"K^{-}");
  h2->GetXaxis()->SetBinLabel(4,"p");
  h2->GetXaxis()->SetBinLabel(5,"d");
  h2->GetXaxis()->SetBinLabel(6,"Other");

  new TH2F( "CDS_Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterFindMass", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterELossCorr", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsdEMean", "Mass^{2} vs. dEMean;Mass^{2} (GeV/c^{2})^{2};dE Mean (MeVee)", 3000, -5, 10, 1000, 0, 100 );
  new TH2F( "CDS_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_AfterFindMass", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 40, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Calc", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_TOFvsMomentum_Calc", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 40, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsEnergy", "1/#beta vs. Momentum;1/#beta;dE (MeV)", 1000, 0, 10, 800, 80, 0 );
  for(int ip=0; ip<8; ip++){
    new TH2F( Form("CDS_OverbetavsMomentum_%s",pname[ip].Data()), "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
    new TH2F( Form("CDS_Mass2vsMomentum_%s",pname[ip].Data()), "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
    new TH2F( Form("CDS_Mass2vsMomentum_TimeOffsetCut_%s",pname[ip].Data()), "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
    new TH2F( Form("CDS_Mass2vsMomentum_woTimeOffsetCut_%s",pname[ip].Data()), "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
    new TH1F( Form("CDS_Mass2Difference_%s",pname[ip].Data()), "Mass^{2} Difference;Mass^{2} Difference (GeV/c^{2})^{2};Counts", 1000, -5, 5 );
    new TH2F( Form("CDS_BetaGammavsMass2Difference_%s",pname[ip].Data()), "BetaGamma vs. Mass^{2} Difference;#beta#gamma;Mass^{2} Difference (GeV/c^{2})^{2}", 3000, -1, 10, 1000, -5, 5 );
  }

  new TH1F( "CDS_TimeBeam", "Time Beam;Time Beam (ns);Counts", 1000, 0, 10);
  new TH1F( "CDS_TOFVTXCDC", "TOFVTXCDC;TOFVTXCDC (ns);Counts", 1000, 0, 10);

  new TH1F( Form("VBDIS_TimeOffsetCut"), Form("VBDIS ;VBDIS (cm);Counts"), 3000, 0.0, 30.0);
  new TH1F( Form("VBDIS_woTimeOffsetCut"), Form("VBDIS ;VBDIS (cm);Counts"), 3000, 0.0, 30.0);
  new TH2F( "Vertex_XY_TimeOffsetCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_TimeOffsetCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_TimeOffsetCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_XY_woTimeOffsetCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_woTimeOffsetCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_woTimeOffsetCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );

  for(int icut=1; icut<=10; icut++){
    new TH2F( Form("Vertex_XY_TimeOffsetCut%d",icut), "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
    new TH1F( Form("VBDIS_TimeOffsetCut%d",icut), Form("VBDIS ;VBDIS (cm);Counts"), 3000, 0.0, 30.0);
    new TH1F( Form("CDC_Chi2_TimeOffsetCut%d",icut), Form("CDC #chi^{2} ;#chi^{2};Counts"), 3000, 0.0, 30.0);
  }
  new TH2F( "Vertex_XY_VBDISCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH1F( "CDC_Chi2_VBDISCut", "CDC #chi^{2} ;#chi^{2};Counts", 3000, 0.0, 30.0);
  for(int icut=0; icut<=21; icut++){
    new TH2F( Form("Vertex_XY_TimeOffsetCut_ZCut%d",icut), "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  }

  /* Hodoscope */
  std::cout << "Define Histograms for Hodoscopes" << std::endl;
  const int nhodo = 2;
  int hodoid[nhodo]={
    CID_CDH,CID_IH
  };
  for(int ihodo=0; ihodo<1; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    for( int seg=1; seg<=nsegs; seg++ ){
      const char* ud[2] = {"u","d"};
      const char* udname[2] = {"up","down"};
      for(int iud=0; iud<2; iud++){
        new TH1F( Form("%s%s%d_TimeOffset",name,ud[iud],seg),   Form("Time Offset %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    400,    -10, 10 );
        new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 100,  1000,    -5, 5 );
        new TH2F( Form("%s%s%d_dEMeanvsTimeOffset",name,ud[iud],seg),   Form("dEMean Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 100,  400,    -10, 10 );
      }
      new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    1000,    -10, 10 );
      new TH2F( Form("%s%d_TimeOffsetvsMomentum",name,seg),   Form("Time Offset %s-%d vs. Momentum;Time (ns);Momentum (GeV/c)",name,seg),    1000,    -10, 10, 800, -2, 2 );
      new TH2F( Form("%s%d_dEMeanvsTimeOffset",name,seg),   Form("dEMean Time Offset corr. %s-%d;dE (MeV);Time (ns)",name,seg),     1000,    0, 100,  400,    -10, 10 );
      new TH1F( Form("%s%d_CTimeSub",name,seg),   Form("CTime Subtract %s-%d;Time (ns);Counts",name,seg),    1000,    -10, 10 );
      new TH1F( Form("%s%d_Multiplicity",name,seg),   Form("Multiplicity %s-%d;Multiplicity;Counts",name,seg),    11,    0, 10 );
      for(int ip=0; ip<8; ip++){
        new TH1F( Form("%s%d_TimeOffset_%s",name,seg,pname[ip].Data()),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    5000,    -5, 5 );
        new TH2F( Form("%s%d_TimeOffsetvsMomentum_%s",name,seg,pname[ip].Data()),   Form("Time Offset %s-%d vs. Momentum;Time (ns);Momentum (GeV/c)",name,seg),    1000,    -5, 5, 800, -2, 2 );
      }
    }

  }



  return true;
}
