// MyAnalysisCDCSlew.cpp

#include "MyAnalysisCDCSlew.h"

MyAnalysisCDCSlew::MyAnalysisCDCSlew(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisCDCSlew::~MyAnalysisCDCSlew()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCDCSlew::Clear()
{
}

bool MyAnalysisCDCSlew::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0); /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  /* Event selection */
  if(particle->nCDS()==0) return true;
  FillHist("EventNumber",1); /* Least 1 hit in CDS*/

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      fcds=it;
    }
  }
  if(fcds==-1) return true; /* Vertex decision */

  /* For CDC Tracking */
  const char* name = "CDC";
  const int nlays = 15;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    CDSTrack* track = cdstrackMan->GoodTrack(cds->id());
    if(track->nCDHHit()!=1) continue;

    double mom = track->Momentum();
    double chi =  track->Chi();
    double beta = cds->beta();
    double overbeta = 1.0/cds->beta()/cds->beta();

    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!track->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
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
    track->GetParameters(CID_CDC,par,cdsvtx);
    mom=track->Momentum();
    if(!track->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
    dis = (beamvtx-cdsvtx).Mag();

    double beam_out=0, beam_tof=0;
    ELossTools::CalcElossBeamTGeo(beam->t0pos(),beamvtx,beam->mom(),kpMass,beam_out,beam_tof);
    TVector3 cdhvtx=track->CDHVertex();
    double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
    beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
    double mass2=mom*mom*(1/(beta*beta)-1);

    double calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(track,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;

    //track->Retiming(cdsMan,conf,beta,true);

    cdhtime=0,dis=0;
    cdhseg=0;
    if(!track->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    cdh=0;
    for(int i=0; i<cdsMan->nCDH(); i++){
      HodoscopeLikeHit* hit = cdsMan->CDH(i);
      if(hit->CheckRange()){
        if(hit->seg()==cdhseg) cdh=hit;
      }
    }
    de = cdh->emean();

    cdsvtx,beamvtx;
    par[5];
    track->GetParameters(CID_CDC,par,cdsvtx);
    mom=track->Momentum();
    if(!track->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
    dis = (beamvtx-cdsvtx).Mag();

    beam_out=0, beam_tof=0;
    ELossTools::CalcElossBeamTGeo(beam->t0pos(),beamvtx,beam->mom(),kpMass,beam_out,beam_tof);
    cdhvtx=track->CDHVertex();
    cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
    beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
    mass2=mom*mom*(1/(beta*beta)-1);

    calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(track,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;

    FillHist(Form("%s_Momentum",name),mom,1);
    FillHist(Form("%s_ChiSquare",name),chi,1);
    if(chi>30) continue;
    for(int layer=1;layer<=nlays;layer++){
      int mult = track->nTrackHit(layer);
      FillHist(Form("%sl%d_MultiplicityinTrack",name,layer),mult,1);
      for(int i=0;i<mult;i++){
        CDCHit *hit = track->TrackHit(cdsMan,layer,i);
        int wire = hit->wire();
        int tdc = hit->tdc();
        double dt = hit->dt();
        double dl = hit->dl();
        double resl = hit->resl();
        double dxdt = conf->GetXTMapManager()->CalcDxDt(CID_CDC,layer,wire,dt);
        //double resldt = conf->GetXTMapManager()->CalcDriftTime(CID_CDC,layer,wire,resl);
        double resldt = resl/dxdt;
        FillHist(Form("%sl%d_Residual",name,layer),resl,1);
        FillHist(Form("%sl%d_ResiTime",name,layer),resldt,1);
        FillHist(Form("%sl%d_OverbetavsResidual",name,layer),overbeta,resl,1);
        FillHist(Form("%sl%d_OverbetavsResiTime",name,layer),overbeta,resldt,1);
        FillHist(Form("%sl%dw%02d_Residual",name,layer,wire),resl,1);
        FillHist(Form("%sl%dw%02d_ResiTime",name,layer,wire),resldt,1);
        FillHist(Form("%sl%dw%02d_OverbetavsResidual",name,layer,wire),overbeta,resl,1);
        FillHist(Form("%sl%dw%02d_OverbetavsResiTime",name,layer,wire),overbeta,resldt,1);
      }
    }
  }

  return true;

}

bool MyAnalysisCDCSlew::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCDCSlew::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCDCSlew::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCDCSlew::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisCDCSlew::CutCondition()
{
}

bool MyAnalysisCDCSlew::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisCDCSlew::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaCDCSlew");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  new TH1F( "EventNumber", "Number of Events", 20, 0, 20);
  std::cout << "Define Histograms for Counters/Chambers" << std::endl;
  /* CDC */
  const int cid= CID_CDC;
  const char* name = "CDC";
  const int nlays= NumOfCDCLayers;
  std::cout << "Define Histograms for " << name << std::endl;
  new TH1F( Form("%s_ChiSquare",name), Form("ChiSquare %s;#Chi^{2}/NDF;Counts",name), 1000, 0, 50 );
  new TH1F( Form("%s_Momentum",name), Form("Momentum %s;Momentum (GeV/c);Counts",name), 500, 0, 5 );
  for( int layer=1; layer<=nlays; layer++ ){
    int nwire = NumOfCDCWiresInLayer[layer-1];
    new TH1F( Form("%sl%d_Residual",name,layer),Form("Residual %s-L%d %s (layer %d);Residual (cm);Counts",name,layer,name,layer),2000,-0.2,0.2);
    new TH2F( Form("%sl%d_OverbetavsResidual",name,layer),Form("Overbeta vs. Residual %s-L%d %s (layer %d);1/#beta^{2};Residual (cm)",name,layer,name,layer),500,0,100,500,-0.1,0.1);
    new TH1F( Form("%sl%d_ResiTime",name,layer),Form("dt %s-L%d;%s (layer %d); dt (ns);Counts",name,layer,name,layer),2000,-25,25);
    new TH2F( Form("%sl%d_OverbetavsResiTime",name,layer),Form("Overbeta vs. dt %s-L%d %s (layer %d);1/#beta^{2};dt (ns)",name,layer,name,layer),500,0,100,500,-25,25);

    for(int wire=1;wire<=nwire;wire++){
    //new TH1F( Form("%sl%dw%02d_Residual",name,layer,wire),Form("Residual %s-L%dW%02d;%s (layer %d-%d); Residual (cm);Counts",name,layer,wire,name,layer,wire),2000,-0.2,0.2);
    //new TH2F( Form("%sl%dw%02d_OverbetavsResidual",name,layer,wire),Form("Overbeta vs. Residual %s-L%dW%02d;%s (layer %d-%d);Overbeta;Residual (cm)",name,layer,wire,name,layer,wire),500,0,100,500,-0.1,0.1);
    //new TH1F( Form("%sl%dw%02d_ResiTime",name,layer,wire),Form("dt %s-L%dW%02d;%s (layer %d-%d); dt (ns);Counts",name,layer,wire,name,layer,wire),2000,-0.2,0.2);
    new TH2F( Form("%sl%dw%02d_OverbetavsResiTime",name,layer,wire),Form("Overbeta vs. dt %s-L%dW%02d;%s (layer %d-%d);1/#beta^{2};dt (ns)",name,layer,wire,name,layer,wire),500,0,100,500,-25,25);
    }
  }

  std::cout << "=== End of [Counters/Chambers::Initialize] === " << std::endl;

  return true;
}
