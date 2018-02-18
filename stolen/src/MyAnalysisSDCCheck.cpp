// MyAnalysisSDCCheck.cpp

#include "MyAnalysisSDCCheck.h"

MyAnalysisSDCCheck::MyAnalysisSDCCheck(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  Clear();
}

MyAnalysisSDCCheck::~MyAnalysisSDCCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisSDCCheck::Clear()
{
  return;
}

bool MyAnalysisSDCCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{
  rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();

  int MulDEF=0;
  for(int i=0; i<blMan->nDEF(); i++){
    HodoscopeLikeHit* hit = blMan->DEF(i);
    if(hit->CheckRange()){
      if(hit->seg()==3) MulDEF++;
    }
  }
  if(MulDEF==0) return true;

  int MulE0=0;
  for(int i=0; i<blMan->nE0(); i++){
    HodoscopeLikeHit* hit = blMan->E0(i);
    if(hit->CheckRange()){
      if(hit->seg()==2) MulE0++;
    }
  }
  if(MulE0==0) return true;

  /* ===================================== */
  /* ========== Event selection ========== */
  double T0Timing = -9999.; 
  double T0dEMean = -9999.;
  double T0dE[2] = {-9999.,-9999.};
  int T0Segment = -1;
  int MulT0 = 0;
  for(int i=0; i<blMan->nT0(); i++){
    if(blMan->T0(i)->CheckRange()){
      MulT0++;
      T0Timing = blMan->T0(i)->ctmean();
      T0Segment = blMan->T0(i)->seg();
      T0dEMean = blMan->T0(i)->ene(0);
      T0dE[0] = blMan->T0(i)->ene(0);
      T0dE[1] = blMan->T0(i)->ene(1);
    }
  }
  if(MulT0!=1){
    return true;
  }
  //if(T0Segment!=3){
  //  return true;
  //}
  int MulBHD = 0;
  int BHDSegment = -1;
  double BHDTiming = -9999.;
  double BHDX = -9999.;
  for(int i=0; i<blMan->nBHD(); i++){
    if(blMan->BHD(i)->CheckRange()){
      if(6<=blMan->BHD(i)->seg()&&blMan->BHD(i)->seg()<=15){
        MulBHD++;
        BHDTiming = blMan->BHD(i)->ctmean();
        BHDSegment = blMan->BHD(i)->seg();
        TVector3 tmppos; conf->GetGeomMapManager()->GetGPos(CID_BHD,BHDSegment,tmppos);
        BHDX = tmppos.X();
      }
    }
  }
  if(MulBHD!=1){
    return true;
  }
  double BeamTOF = T0Timing-BHDTiming;
  FillHist(Form("Beam_TOF"), BeamTOF);
  //if((T0Timing-BHDTiming)<27.8||29.5<(T0Timing-BHDTiming)){ /* Kaon */
  //  return true;
  //}
  FillHist(Form("Beam_TOFwithCut"), BeamTOF);

  // #################### //
  // DC tracking          //
  // #################### //
  bltrackMan->DoTracking(blMan,conf,true,true);
  const int ndc = 4;
  int dccid[ndc] = {CID_BLC1,CID_BLC2,CID_BPC,CID_FDC1};
  int nhitdc[ndc] = {0,0,0,0};
  int trackid[ndc] = {-1,-1,-1,-1};
  double ll[ndc] = {-30.0,-30.0,-30.0,-30.0};
  double ul[ndc] = {100.0,100.0,100.0,100.0};
  double trigll[ndc] = {-10.0,-10.0,-10.0,-10.0};
  double trigul[ndc] = { 10.0, 10.0, 10.0, 10.0};
  double chiul[ndc] = {10.0, 10.0, 10.0, 10.0};
  int idc;
  const char* dcname;
  // ==== BLC1 and BLC2 analysis ==== // 
  for(idc=0;idc<2;idc++){
    dcname = dlist->GetName(dccid[idc]).data();
    FillHist(Form("%s_nTrack",dcname),bltrackMan->ntrackBLDC(dccid[idc]));
    for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
      if(i>=1) break;
      LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
      double tracktime = track -> GetTrackTime();
      double chi2 = track -> chi2all();
      FillHist(Form("%s_TrackTime",dcname),tracktime);
      if(ll[idc]<tracktime&&tracktime<ul[idc]){
        FillHist(Form("%s_Chi2",dcname),chi2);
        if(chi2<chiul[idc]){
          nhitdc[idc]++;
          trackid[idc] = i;
          FillHist(Form("%s_TrackTimewCut",dcname),tracktime);
          FillHist(Form("%s_Chi2wCut",dcname),chi2);
        }
      }
    }
    FillHist(Form("%s_nGoodTrack",dcname),nhitdc[idc]);
  }
  bool BLCGoodflag = false;
  if(nhitdc[0]==1&&nhitdc[1]==1){
    LocalTrack* blc1 = bltrackMan->trackBLC1(trackid[0]);
    LocalTrack* blc2 = bltrackMan->trackBLC2(trackid[1]);
    double blc1time = blc1->GetTrackTime();
    double blc2time = blc2->GetTrackTime();
    if(trigll[0]<blc1time&&blc1time<trigul[0]
        && trigll[1]<blc2time&&blc2time<trigul[1])
      BLCGoodflag = true;
  }
  /* ## BLC1/2 1 good track request ## */
  if(!BLCGoodflag) return false;
  LocalTrack* trackblc1 = bltrackMan->trackBLC1(trackid[0]);
  LocalTrack* trackblc2 = bltrackMan->trackBLC2(trackid[1]);
  // ==== Beam momentum analysis ==== // 
  BeamSpectrometer* D5 = new BeamSpectrometer(conf);
  D5->TMinuitFit(trackblc1,trackblc2,conf);
  double BeamChi2 = D5->chisquare();
  double BeamMom = D5->mom();
  FillHist(Form("Beam_Momentum"), BeamMom);
  FillHist(Form("Beam_MomentumvsTOF"), BeamMom,BeamTOF);
  FillHist(Form("Beam_Chi2"), BeamChi2);
  delete D5;
  /* ## Beam chisquare < 20 request ## */
  if(!(BeamChi2<20)) return false;
  FillHist(Form("Beam_MomentumwCut"), BeamMom);
  FillHist(Form("Beam_MomentumvsTOFwCut"), BeamMom,BeamTOF);
  FillHist(Form("Beam_Chi2wCut"), BeamChi2);

  // ==== SDC analysis ==== // 
  const int cid = CID_BPC;
  const char* name = dlist->GetName(cid).data();
  const int nlays = dlist->GetNlayers(cid);
  int nhit[nlays];
  int hitwire[nlays][16];
  double hitdt[nlays][16];
  for(int layer=1; layer<=nlays; layer++){
    int mult = blMan->nBLDC(cid,layer);
    nhit[layer-1] = mult;
    FillHist(Form("%sl%d_Multiplicity",name,layer),mult);
    int n = 0;
    for(int i=0; i<mult; i++){
      ChamberLikeHit* hit = blMan->BLDC(cid,layer,i);
      int wire = hit->wire();
      int tdc = hit->tdc();
      double dt = hit->dt();
      double dl = hit->dl();
      hitwire[layer-1][n] = wire;
      hitdt[layer-1][n] = dt;
      TVector3 wpos = hit->wpos();
      TVector3 blcpos = trackblc2->GetPosatZ(wpos.Z());

      FillHist(Form("%sl%d_WirevsBLCPosX",name,layer),wire,blcpos.X(),1);
      FillHist(Form("%sl%d_WirevsBLCPosY",name,layer),wire,blcpos.Y(),1);
      
      FillHist(Form("%sl%d_HitPattern",name,layer),wire,1);
      FillHist(Form("%sl%d_TDC",name,layer),tdc,1);
      FillHist(Form("%sl%d_Time",name,layer),dt,1);
      FillHist(Form("%sl%d_Length",name,layer),dl,1);
      FillHist(Form("%sl%d_TimevsLength",name,layer),dt,dl,1);
      FillHist(Form("%sl%dw%02d_TDC",name,layer,wire),tdc,1);
      FillHist(Form("%sl%dw%02d_Time",name,layer,wire),dt,1);
      FillHist(Form("%sl%dw%02d_Length",name,layer,wire),dl,1);
      FillHist(Form("%sl%dw%02d_TimevsLength",name,layer,wire),dt,dl,1);
      n++;
    }
  }

  // Layer Hit Efficiency //
  FillHist(Form("%s_HitEfficiency18Hit",name),9,1);
  if(nhit[0]==1&&nhit[nlays-1]==1){
    FillHist(Form("%s_HitEfficiency18Hit",name),0,1);
    for(int layer=2; layer<=nlays-1; layer++){
      FillHist(Form("%sl%d_HitEfficiency18Hit",name,layer),0,1);
      if(nhit[layer-1]!=1){ continue; }
      for(int i=0; i<nhit[layer-1]; i++){
        int wire = hitwire[layer-1][i];
        double dt = hitdt[layer-1][i];
        if(-10<dt&&dt<100){
          FillHist(Form("%s_HitEfficiency18Hit",name),layer,1);
          FillHist(Form("%sl%d_HitEfficiency18Hit",name,layer),wire,1);
          break;
        }
      }
    }
  }
    FillHist(Form("%s_HitEfficiency27Hit",name),9,1);
  if(nhit[1]==1&&nhit[nlays-2]==1){
    FillHist(Form("%s_HitEfficiency27Hit",name),0,1);
    for(int layer=1; layer<=nlays; layer++){
      FillHist(Form("%sl%d_HitEfficiency27Hit",name,layer),0,1);
      if(layer!=1&&layer!=nlays){ continue; }
      if(nhit[layer-1]!=1){ continue; }
      for(int i=0; i<nhit[layer-1]; i++){
        int wire = hitwire[layer-1][i];
        double dt = hitdt[layer-1][i];
        if(-10<dt&&dt<100){
          FillHist(Form("%s_HitEfficiency27Hit",name),layer,1);
          FillHist(Form("%sl%d_HitEfficiency27Hit",name,layer),wire,1);
          break;
        }
      }
    }
  }

  /* Track data filling */
  const int ndc2=1;
  int dccid2[ndc2]={
    CID_BPC,
  };
  char LayConfig[ndc2][16] = {
    {"XXYYXXYY"},
  };

  for(int idc=0;idc<ndc2;idc++){
    const int cid= dccid2[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    FillHist(Form("%s_nTrack",name),bltrackMan->ntrackBLDC(cid),1);
    if(bltrackMan->ntrackBLDC(cid)!=1) continue;
    for(int itra=0; itra<bltrackMan->ntrackBLDC(cid); itra++){
      LocalTrack* track = bltrackMan->trackBLDC(cid,itra);
      double tracktime = track->GetTrackTime();
      //if(tracktime<-10||10<tracktime) continue;
      double chi2 = track -> chi2all();
      TVector3 mom = track->GetMomDir();
      FillHist(Form("%s_Direction",name),mom.X()/mom.Z(),mom.Y()/mom.Z(),1);
      FillHist(Form("%s_TrackTime",name),tracktime,1);
      FillHist(Form("%s_ChiSquare",name),chi2,1);
      if(chi2>10) continue;
      if(idc>5){ continue; }
      for(int ihit=0; ihit<track->nhit(); ihit++){
        ChamberLikeHit* hit = track->hit(ihit);
        int layer = hit->layer();
        int wire = hit->wire();
        int tdc = hit->tdc();
        double dt = hit->dt();
        double resl = hit->resl();
        FillHist(Form("%sl%d_HitPattern_inTrack",name,layer),wire,1);
        FillHist(Form("%sl%d_TDC_inTrack",name,layer),tdc,1);
        FillHist(Form("%sl%d_Time_inTrack",name,layer),dt,1);
        FillHist(Form("%sl%d_Residual_inTrack",name,layer),resl,1);
        FillHist(Form("%sl%d_TimevsResidual_inTrack",name,layer),dt,resl,1);
        FillHist(Form("%sl%dw%02d_TDC_inTrack",name,layer,wire),tdc,1);
        FillHist(Form("%sl%dw%02d_Time_inTrack",name,layer,wire),dt,1);
        FillHist(Form("%sl%dw%02d_Residual_inTrack",name,layer,wire),resl,1);
        FillHist(Form("%sl%dw%02d_TimevsResidual_inTrack",name,layer,wire),dt,resl,1);
        FillHist(Form("%s_HitPattern%c_inTrack",name,LayConfig[idc][layer-1]),layer,wire,1);
      }

      for(int ixy=0; ixy<2; ixy++){
        for(int iclust=0; iclust<track->ncluster(ixy); iclust++){
          BLDCCluster* cluster = track->cluster(ixy,iclust);
          int layer = 0;
          if(cluster->nhit()!=2) continue;
          ChamberLikeHit* hit = cluster->hit(0);
          layer = hit->layer();
          double timediff = cluster->GetTimeSub();
          double timemean = cluster->GetTimeMean();
          double ctimemean = cluster->GetCTime();
          FillHist(Form("%sl%d_TimeDifferencevsTimeMean_inTrack",name,layer),timediff,timemean,1);
          FillHist(Form("%sl%d_TimeDifferencevsCTimeMean_inTrack",name,layer),timediff,ctimemean,1);
        }
      }
    }
  }


  return true;
}

bool MyAnalysisSDCCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisSDCCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisSDCCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisSDCCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

bool MyAnalysisSDCCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisSDCCheck::Initialize ###" << std::endl;
  DetectorList *dlist=DetectorList::GetInstance();
  confFile = confMan;
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaSDCCheck");
  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  // Chamber
  std::cout << "Define Histograms for DCs" << std::endl;
  const int ndc=1;
  int dccid[ndc]={CID_BPC};
  int LayConfigNum[ndc] = {2};
  char LayConfig[ndc][6] = {
    {"XY"},
  };

  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    new TH1F( Form("%s_HitEfficiency18Hit",name),Form("Hit Efficiency %s with 1&8 Hit;Layer;Counts",name),nlays+2,-0.5,nlays+1.5);
    new TH1F( Form("%s_HitEfficiency27Hit",name),Form("Hit Efficiency %s with 1&8 Hit;Layer;Counts",name),nlays+2,-0.5,nlays+1.5);
    new TH1F( Form("%s_HitInEfficiency18Hit",name),Form("Hit Efficiency %s with 1&8 Hit;Layer;Counts",name),nlays+2,-0.5,nlays+1.5);
    new TH1F( Form("%s_HitInEfficiency27Hit",name),Form("Hit Efficiency %s with 1&8 Hit;Layer;Counts",name),nlays+2,-0.5,nlays+1.5);
    for( int layer=1; layer<=nlays; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire(cid,layer);
      new TH1F( Form("%sl%d_Multiplicity_inTrack",name,layer),Form("Multiplicity %s-L%d;Multiplicity;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_HitPattern_inTrack",name,layer),Form("Hit pattern %s-L%d;Wire;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_TDC_inTrack",name,layer),Form("TDC %s-L%d;TDC ch.;Counts",name,layer),1000,0,2000);
      new TH1F( Form("%sl%d_Time_inTrack",name,layer),Form("Time %s-L%d;Time (ns);Counts",name,layer),500,-100,400);
      new TH1F( Form("%sl%d_Residual_inTrack",name,layer),Form("Residual %s-L%d;Residual (cm);Counts",name,layer),500,-0.5,0.5);
      new TH2F( Form("%sl%d_TimevsResidual_inTrack",name,layer),Form("Time vs. Residual %s-L%d;Time (ns);Residual (cm)",name,layer),500,-100,400,500,-0.5,0.5);
      if(layer%2==1){
        new TH2F( Form("%sl%d_TimeDifferencevsTimeMean_inTrack",name,layer),Form("Time Difference vs. Time Mean %s-L%d and L%d;Time Difference (ns);Time Mean (ns)",name,layer,layer+1),500,-250,250,500,-100,400);
        new TH2F( Form("%sl%d_TimeDifferencevsCTimeMean_inTrack",name,layer),Form("Time Difference vs. Corrected Time Mean %s-L%d and L%d;Time Difference (ns);Time Mean (ns)",name,layer,layer+1),500,-250,250,500,-100,400);
      }
      new TH2F( Form("%sl%d_WirevsBLCPosX",name,layer),Form("%s-L%d Wire vs. BLC2 Track Position;%s-L%d Wire;BLC2 X Position",name,layer,name,layer),nwire+1,-0.5,nwire+0.5,200,-10,10);
      new TH2F( Form("%sl%d_WirevsBLCPosY",name,layer),Form("%s-L%d Wire vs. BLC2 Track Position;%s-L%d Wire;BLC2 Y Position",name,layer,name,layer),nwire+1,-0.5,nwire+0.5,200,-10,10);
      new TH1F( Form("%sl%d_Multiplicity",name,layer),Form("Multiplicity %s-L%d;Multiplicity;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_HitPattern",name,layer),Form("Hit pattern %s-L%d;Wire;Counts",name,layer),nwire+1,-0.5,nwire+0.5);
      new TH1F( Form("%sl%d_TDC",name,layer),Form("TDC %s-L%d;TDC ch.;Counts",name,layer),1000,0,2000);
      new TH1F( Form("%sl%d_Time",name,layer),Form("Time %s-L%d;Time (ns);Counts",name,layer),500,-100,400);
      new TH1F( Form("%sl%d_Length",name,layer),Form("Length %s-L%d;Length (cm);Counts",name,layer),250,0.0,0.5);
      new TH2F( Form("%sl%d_TimevsLength",name,layer),Form("Time vs. Length %s-L%d;Time (ns);Length (cm)",name,layer),500,-100,400,250,0.0,0.5);
      new TH1F( Form("%sl%d_HitEfficiency18Hit",name,layer),Form("Hit Efficiency %s-L%d with 1&8 Hit;Layer;Counts",name,layer),nwire+2,-0.5,nwire+1.5);
      new TH1F( Form("%sl%d_HitEfficiency27Hit",name,layer),Form("Hit Efficiency %s-L%d with 1&8 Hit;Layer;Counts",name,layer),nwire+2,-0.5,nwire+1.5);
      new TH1F( Form("%sl%d_HitInEfficiency18Hit",name,layer),Form("Hit Efficiency %s-L%d with 1&8 Hit;Layer;Counts",name,layer),nwire+2,-0.5,nwire+1.5);
      new TH1F( Form("%sl%d_HitInEfficiency27Hit",name,layer),Form("Hit Efficiency %s-L%d with 1&8 Hit;Layer;Counts",name,layer),nwire+2,-0.5,nwire+1.5);
      for(int wire=1;wire<=nwire;wire++){
        new TH1F( Form("%sl%dw%02d_Time",name,layer,wire),Form("Time %s-L%dW%02d;Time (ns);Counts",name,layer,wire),500,-100,400);
        new TH1F( Form("%sl%dw%02d_Time_inTrack",name,layer,wire),Form("Time %s-L%dW%02d;Time (ns);Counts",name,layer,wire),500,-100,400);
      }
    }
    /* For track */
    new TH2F( Form("%s_Direction",name),Form("%s Track Direction;dX/dZ;dY/dZ",name),100,-0.25,0.25,100,-0.25,0.25);
    new TH1F( Form("%s_nTrack",name),Form("Number of tracks %s;Number of tracks;Counts",name),11,-0.5,10.5);
    new TH1F( Form("%s_ChiSquare",name),Form("#Chi^{2}/NDF %s;#Chi^{2}/NDF;Counts",name),1000,0,100);
    new TH1F( Form("%s_TrackTime",name),Form("Track time %s;Track time (ns);Counts",name),1000,-100,400);
    
  }


  std::cout << "=== End of [MyAnalysisSDCCheck::Initialize] === " << std::endl;

  return true;
}
