// MyAnalysisCDCTrackingwSlew.cpp

#include "MyAnalysisCDCTrackingwSlew.h"

MyAnalysisCDCTrackingwSlew::MyAnalysisCDCTrackingwSlew(TFile* rt, ConfMan* conf)
{
  Initialize(rt, conf);
  Clear();
}

MyAnalysisCDCTrackingwSlew::~MyAnalysisCDCTrackingwSlew()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCDCTrackingwSlew::Clear()
{
    //cdsMan->Clear();
    //blMan->Clear();
    evheader->Clear();
    trackMan->Clear();
}

bool MyAnalysisCDCTrackingwSlew::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  // ###################### //
  // CDS all track analysis //
  // ###################### //
  for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
    CDSTrack* cdc=cdstrackMan->GoodTrack(id);

    if(cdc->nCDHHit()!=1) continue;
    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    HodoscopeLikeHit* cdh=0;
    for(int i=0; i<cdsMan->nCDH(); i++){
      HodoscopeLikeHit* hit = cdsMan->CDH(i);
      if(hit->CheckRange()){
        if(hit->seg()==cdhseg) cdh=hit;
      }
    }
    double de = cdh->emean();

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

    cdc->Retiming(cdsMan,conf,beta,true);
  }

  int ngoodTrack=cdstrackMan->nGoodTrack();
  int nallTrack=cdstrackMan->nTrack();

  trackMan = cdstrackMan;

  if( nallTrack<1 ){
    return true;
  }
  else{
    evTree->Fill();
  }
  return true;

}

bool MyAnalysisCDCTrackingwSlew::FillHist(TString name, double val1)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1){
    h1 -> Fill(val1);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisCDCTrackingwSlew::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisCDCTrackingwSlew::Initialize(TFile* rt, ConfMan* confMan)
{
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_cdctracking_slew");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  evTree = new TTree( "EventTree", "EventTree" );
  evheader = new EventHeader();
  trackMan = new CDSTrackingMan();  
  evTree->Branch( "EventHeader", &evheader );
  evTree->Branch( "CDSTrackingMan", &trackMan );

  std::cout << "=== End of [MyAnalysisCDCTrackingwSlew::Initialize] === " << std::endl;

  return true;
}
