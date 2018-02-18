// MyAnalysisCDS.cpp

#include "MyAnalysisCDS.h"

MyAnalysisCDS::MyAnalysisCDS(TFile* rtFile, ConfMan* conf)
{
  Initialize(rtFile, conf);
  Clear();
}

MyAnalysisCDS::~MyAnalysisCDS()
{
}

void MyAnalysisCDS::Clear()
{
}

bool MyAnalysisCDS::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  FillHist("CDC_nGoodTrack",cdstrackMan->nGoodTrack());
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  // ################ //
  // CDS track        //
  // ################ //
  for( int it=0; it<cdstrackMan->nGoodTrack(); it++ ){
    CDSTrack* cdctrack = cdstrackMan->Track(cdstrackMan->GoodTrackID(it));

    double cdhtime, dis;
    int cdhseg;
    if(!cdctrack->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    TVector3 vtxb, vtxh;
    double par[5];
    cdctrack->GetParameters(CID_CDC,par,vtxb);
    double mom = cdctrack->Momentum();
    cdctrack->GetVertex(beam->bpcpos(),beam->bpcdir(),vtxb,vtxh);

    pCDS* cds = new pCDS();
    cds->SetParameters(par);
    cds->SetTOF(cdhtime);
    cds->SetCDHSeg(cdhseg);

    // === PID === //
    TVector3 cdhvtx=cdctrack->CDHVertex();
    double time_beam=(beam->t0pos() - vtxb).Mag()/beam->beta()/(Const*100);
    double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,vtxh);
    double beta=cdc_dis/(cdhtime-beam->t0time()-time_beam)/(Const*100);
    double mass2=mom*mom*(1/(beta*beta)-1);
    double beta_calc,tofvtxcdc;
    double tof=cdhtime-beam->t0time();
    TrackTools::FindMass2(cdctrack,beam,tof,beta_calc,mass2,tofvtxcdc);
    int pid=TrackTools::PID(mom,mass2);
    cdctrack->SetPID(pid);

    // === Energy loss correction === //
    double tmp1;
    cdctrack->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],vtxb,vtxh,tofvtxcdc,tmp1,true);
    cdctrack->GetVertex(beam->bpcpos(),beam->bpcdir(),vtxb,vtxh);
    TVector3 Pp;
    if(!cdctrack->GetMomentum(vtxh,Pp,true,true)) continue;

    // === Fill histogram === //
    FillHist("Vertex_Z",vtxb.Z());
    FillHist("Vertex_ZX",vtxb.Z(),vtxb.X());
    FillHist("Vertex_ZY",vtxb.Z(),vtxb.Y());
    FillHist("Vertex_XY",vtxb.X(),vtxb.Y());
    if( GeomTools::GetID(vtxb)==CID_Fiducial ){
      FillHist("OverbetavsMomentum",1.0/beta,mom);
      FillHist("TOFvsMomentum",tof,mom);
      FillHist("Mass2vsMomentum",mass2,mom);
      FillHist("CDC_Chi2",cdctrack->Chi());
      FillHist("VertexwTRG_Z",vtxb.Z());
      FillHist("VertexwTRG_ZX",vtxb.Z(),vtxb.X());
      FillHist("VertexwTRG_ZY",vtxb.Z(),vtxb.Y());
      FillHist("VertexwTRG_XY",vtxb.X(),vtxb.Y());
    }

    // === Set pCDS === //
    cds->SetTrackID(it);
    cds->SetVertexBeam(vtxb); 
    cds->SetVertexCDC(vtxh);
    cds->SetVDis(dis);
    cds->SetMomDir(Pp);
    cds->SetPID(pid);
    cds->SetMomentum(mom/TMath::Abs(mom)*Pp.Mag());
    cds->SetBeta(beta);
    cds->SetMass(sqrt(mass2));
    cds->SetMass2(mass2);
    cds->SetAngleLab(beam->bpcdir().Angle(Pp));
    for(int icdh=0; icdh<cdctrack->nCDHHit();icdh++)
      cds->SetCDHSeg(cdctrack->CDHHit(cdsMan,icdh)->seg());
    for(int iih=0; iih<cdctrack->nIHHit();iih++)
      cds->SetCDHSeg(cdctrack->IHHit(cdsMan,iih)->seg());
    TVector3 pos1;
    cdctrack->GetParameters(CID_CDCCFRP,par,pos1);
    cdctrack->GetParameters(CID_CDC,par,vtxb);
    cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,pos1);
    double beta1=fabs(mom)/sqrt(cdsMass[pid]*cdsMass[pid]+mom*mom);
    double time_cdc=cdc_dis/beta1/(Const*100);
    double cdcdt=(tof-tofvtxcdc)-time_cdc;
    cds->SetFL(cdc_dis);
    cds->SetDt(cdcdt);
    cds->SetChi2(cdctrack->Chi());
    particle->AddCDS(*cds);
    delete cds;

  }

  return true;

}

bool MyAnalysisCDS::FillHist(TString name, double val1)
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

bool MyAnalysisCDS::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisCDS::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  // CDC
  const int ndc=1;
  int dccid[ndc]={CID_CDC};
  TString dcname[ndc]={"CDC"};
  int NumOfLayers[ndc]={NumOfCDCLayers};

  for(int idc=0;idc<ndc;idc++){
    std::cout << "Define Histgram for " << dcname[idc] << std::endl;
    new TH1F( Form("%s_nGoodTrack",dcname[idc].Data()), Form("N good track %s;# of good tracks;Counts",dcname[idc].Data()), 10, 0, 10 );
    new TH1F( Form("%s_Chi2",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2};Counts",dcname[idc].Data()), 101, 0.0, 101.0 );
  }
  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
  new TH2F( "Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 1500, -5, 10, 500, -5, 5 );
  new TH2F( "OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 20, 500, -5, 5 );

  // Vertex
  std::cout << "Define Histograms for Vertex" << std::endl;
  new TH2F( "Vertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "VertexwTRG_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexwTRG_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexwTRG_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "VertexwTRG_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
}
