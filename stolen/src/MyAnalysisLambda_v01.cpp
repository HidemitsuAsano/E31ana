// MyAnalysisLambda.cpp

#include "MyAnalysisLambda.h"

MyAnalysisLambda::MyAnalysisLambda(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisLambda::~MyAnalysisLambda()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisLambda::Clear()
{
}

bool MyAnalysisLambda::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* _particle)
{

  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();
  if(_particle->nBeam()!=1) return false;
  pBeam* beam = _particle->beam(0);

  Particle* particle = new Particle();
  particle -> AddBeam(*beam);

  /* Trigger selection */

  /* Event selection */
  if(cdstrackMan->nGoodTrack()!=2) return true;
  int MulCDH = 0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    if(cdsMan->CDH(i)->CheckRange()){
      MulCDH++;
    }
  }
  if(MulCDH!=2){
    return true;
  }

  // ###################### //
  // CDS all track analysis //
  // ###################### //
  for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
    CDSTrack* cdc=cdstrackMan->GoodTrack(id);
    FillHist("CDC_Chi2",cdc->Chi());

    if(cdc->Chi()>30) continue;
    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    HodoscopeLikeHit* cdh;
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
    FillHist("CDS_OverbetavsMomentum",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum",mass2,mom);
    FillHist("CDS_dEvsMomentum",de,mom);

    double calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;
    FillHist("CDS_OverbetavsMomentum_AfterFindMass",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum_AfterFindMass",mass2,mom);
    FillHist("CDS_dEvsMomentum_AfterFindMass",de,mom);

    int pid=TrackTools::PID(mom,mass2);
    pid=TrackTools::PIDcorr2(mom,mass2);
    cdc->SetPID(pid);

    if(GeomTools::GetID(beamvtx)==CID_Fiducial){
      FillHist("CDS_OverbetavsMomentum_All",1./beta,mom);
      FillHist("CDS_Mass2vsMomentum_All",mass2,mom);
      if(pid==CDS_PiPlus){
        FillHist("CDS_OverbetavsMomentum_PiPlus",1./beta,mom);
        FillHist("CDS_Mass2vsMomentum_PiPlus",mass2,mom);
      }
      if(pid==CDS_PiMinus){
        FillHist("CDS_OverbetavsMomentum_PiMinus",1./beta,mom);
        FillHist("CDS_Mass2vsMomentum_PiMinus",mass2,mom);
      }
      if(pid==CDS_Kaon){
        FillHist("CDS_OverbetavsMomentum_Kaon",1./beta,mom);
        FillHist("CDS_Mass2vsMomentum_Kaon",mass2,mom);
      }
      if(pid==CDS_Proton){
        FillHist("CDS_OverbetavsMomentum_Proton",1./beta,mom);
        FillHist("CDS_Mass2vsMomentum_Proton",mass2,mom);
      }
      if(pid==CDS_Deuteron){
        FillHist("CDS_OverbetavsMomentum_Deuteron",1./beta,mom);
        FillHist("CDS_Mass2vsMomentum_Deuteron",mass2,mom);
      }
      if(pid==CDS_Other){
        FillHist("CDS_OverbetavsMomentum_Other",1./beta,mom);
        FillHist("CDS_Mass2vsMomentum_Other",mass2,mom);
      }
    }

    TVector3 cdsvtx2,beamvtx2;
    double tof=0,tmpl=0;
    if(!cdc->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],beamvtx2,cdsvtx2,tof,tmpl,true)) continue;
    TVector3 Pp;
    if(!cdc->GetMomentum(cdsvtx2,Pp,true,true)) continue;
    mom = mom/TMath::Abs(mom)*Pp.Mag();
    FillHist("CDS_OverbetavsMomentum_AfterELossCorr",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum_AfterELossCorr",mass2,mom);
    //FillHist("CDS_dEvsMomentum_AfterELossCorr",de,mom);

    /* Set pCDS */
    pCDS* cdstrack = new pCDS();
    cdstrack->SetTrackID(id);
    cdstrack->SetPID(pid);
    cdstrack->SetMomentum(mom);
    cdstrack->SetMass(TMath::Sqrt(mass2));
    cdstrack->SetMass2(mass2);
    cdstrack->SetBeta(beta);
    cdstrack->SetTOF(cdhtime-beam->t0time());
    cdstrack->SetFL(tmpl);
    cdstrack->SetChi2(cdc->Chi());
    cdstrack->SetVertexCDC(cdsvtx2);
    cdstrack->SetVertexBeam(beamvtx2);
    cdstrack->SetMomDir(Pp);
    cdstrack->SetVBDis((beamvtx2-cdsvtx2).Mag());
    for(int icdh=0; icdh<cdc->nCDHHit(); icdh++){
      cdstrack->SetCDHSeg(cdc->CDHHit(cdsMan,icdh)->seg());
    }
    for(int iih=0; iih<cdc->nIHHit(); iih++){
      cdstrack->SetIHSeg(cdc->IHHit(cdsMan,iih)->seg());
    }

    particle->AddCDS(*cdstrack);
  }

  // ###################### //
  // CDS 2 track analysis //
  // ###################### //
  for( int id1=0; id1<particle->nCDS(); id1++ ){
    for( int id2=id1+1; id2<particle->nCDS(); id2++ ){
      pCDS* cds1 = particle->cdsi(id1);
      pCDS* cds2 = particle->cdsi(id2);
      CDSTrack* cdc1=cdstrackMan->GoodTrack(cds1->id());
      CDSTrack* cdc2=cdstrackMan->GoodTrack(cds2->id());
      if(cdc1->Chi()>30||cdc2->Chi()>30){ continue; }
      TVector3 vtx=DEFVECT,vtx1=DEFVECT,vtx2=DEFVECT;
      if(!TrackTools::Calc2HelixVertex(cdc1,cdc2,vtx1,vtx2)) continue;
      vtx = (vtx1+vtx2)*0.5;

      TVector3 Pp1=DEFVECT,Pp2=DEFVECT;
      if(!cdc1->GetMomentum(vtx1,Pp1,true,true)) continue;
      if(!cdc2->GetMomentum(vtx2,Pp2,true,true)) continue;
      TVector3 Pp = Pp1+Pp2;
      TLorentzVector L1; L1.SetVectM(Pp1,cds1->pdgmass());
      TLorentzVector L2; L2.SetVectM(Pp2,cds2->pdgmass());
      double im = (L1+L2).M();
      double dist=0, dltmp=0;
      TVector3 xest=DEFVECT, nest=DEFVECT;
      MathTools::LineToLine(vtx,Pp.Unit(),beam->bpcpos(),beam->bpcdir(),dltmp,dist,xest,nest);

      pCDS* pro = new pCDS();
      pro->SetDaughterID1(cds1->id());
      pro->SetDaughterID2(cds2->id());
      pro->SetCombID((int)(pow(2,cds1->pid())+pow(2,cds2->pid())));
      pro->SetMomentum(Pp.Mag());
      pro->SetMass(im);
      pro->SetVertex(vtx);
      pro->SetVertexBeam(xest);
      pro->SetMomDir(Pp.Unit());
      pro->SetVDis((vtx1-vtx2).Mag());
      pro->SetVBDis((xest-vtx).Mag());
      pro->SetPBDCA(dist);
      pro->SetOA(Pp1.Angle(Pp2));
      particle->AddProduct(*pro);

      if(pro->comb()==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        FillHist("Vertex_XY", xest.X(), xest.Y());
        FillHist("Vertex_ZX", xest.Z(), xest.X());
        FillHist("Vertex_ZY", xest.Z(), xest.Y());
        FillHist("Vertex_Z", xest.Z());
        if(-0.1<xest.Z()&&xest.Z()<0.1){
          FillHist("Vertex_XY_ZCut", xest.X(), xest.Y());
        }
        if(-0.1<xest.X()&&xest.X()<0.1&&-0.1<xest.Y()&&xest.Y()<0.1){
          FillHist("Vertex_Z_XYCut", xest.Z());
        }

        if(GeomTools::GetID(xest)==CID_Fiducial){
          FillHist("Vertex_XY_Fiducial", xest.X(), xest.Y());
          FillHist("Vertex_ZX_Fiducial", xest.Z(), xest.X());
          FillHist("Vertex_ZY_Fiducial", xest.Z(), xest.Y());
          FillHist("Vertex_Z_Fiducial", xest.Z());
          if(-0.1<xest.Z()&&xest.Z()<0.1){
            FillHist("Vertex_XY_FiducialZCut", xest.X(), xest.Y());
          }
          if(-0.1<xest.X()&&xest.X()<0.1&&-0.1<xest.Y()&&xest.Y()<0.1){
            FillHist("Vertex_Z_FiducialXYCut", xest.Z());
          }
        }

        if(GeomTools::GetID(xest)!=CID_Fiducial){ continue; }
        TLorentzVector Lcds = pro->GetLorentzVector();
        double cds_mass  = (Lcds).M();
        double cds_mom   = (Lcds).Vect().Mag();
        double cds_cost  = (Lcds).Vect().CosTheta();
        double cds_phi   = (Lcds).Vect().Phi();
        double cds_fl = pro->vbdis()/(cds_mom/cds_mass);
        double pim_mom = 0.0;
        double pim_bg = 0.0;
        double pim_beta = 0.0;
        double p_mom = 0.0;
        double p_bg = 0.0;
        double p_beta = 0.0;
        if(cds1->pid()==CDS_PiMinus){
          pim_mom = Pp1.Mag();
          p_mom = Pp2.Mag();
          pim_bg = Pp1.Mag()/cds1->pdgmass();
          p_bg = Pp2.Mag()/cds2->pdgmass();
          pim_beta = L1.Beta();
          p_beta = L2.Beta();
        }
        else{
          pim_mom = Pp2.Mag();
          p_mom = Pp1.Mag();
          pim_bg = Pp2.Mag()/cds2->pdgmass();
          p_bg = Pp1.Mag()/cds1->pdgmass();
          pim_beta = L2.Beta();
          p_beta = L1.Beta();
        }
        FillHist("ppim_Mass",cds_mass);
        FillHist("ppim_Momentum",cds_mom);
        FillHist("ppim_CosT",cds_cost);
        FillHist("ppim_FL",cds_fl);
        FillHist("ppim_VDIS",pro->vdis());
        FillHist("ppim_VBDIS",pro->vbdis());
        FillHist("ppim_PBDCA",pro->pbdca());
        FillHist("ppim_Massvspim_BetaGamma",pim_bg,cds_mass);
        FillHist("ppim_Massvsp_BetaGamma",p_bg,cds_mass);
        FillHist("ppim_Massvspim_OverBeta2",1.0/pim_beta/pim_beta,cds_mass);
        FillHist("ppim_Massvsp_OverBeta2",1.0/p_beta/p_beta,cds_mass);
        FillHist("pim_Momentum",pim_mom);
        FillHist("p_Momentum",p_mom);
        /* For Lambda */
        if(lamll<cds_mass&&cds_mass<lamul){
          FillHist(Form("lam_Mass"),cds_mass);
          FillHist(Form("lam_Momentum"),cds_mom);
          FillHist(Form("lam_CosT"),cds_cost);
          FillHist("lam_FL",cds_fl);
          FillHist("lam_VDIS",pro->vdis());
          FillHist("lam_VBDIS",pro->vbdis());
          FillHist("lam_PBDCA",pro->pbdca());
          FillHist("lam_Massvspim_BetaGamma",pim_bg,cds_mass);
          FillHist("lam_Massvsp_BetaGamma",p_bg,cds_mass);
          FillHist("lam_Massvspim_OverBeta2",1.0/pim_beta/pim_beta,cds_mass);
          FillHist("lam_Massvsp_OverBeta2",1.0/p_beta/p_beta,cds_mass);
          FillHist("pim_Momentum_lam",pim_mom);
          FillHist("p_Momentum_lam",p_mom);
        }
        if(sbllamll<cds_mass&&cds_mass<sbllamul){
          FillHist(Form("sbllam_Mass"),cds_mass);
          FillHist(Form("sbllam_Momentum"),cds_mom);
          FillHist(Form("sbllam_CosT"),cds_cost);
          FillHist("sbllam_FL",cds_fl);
          FillHist("sbllam_VDIS",pro->vdis());
          FillHist("sbllam_VBDIS",pro->vbdis());
          FillHist("sbllam_PBDCA",pro->pbdca());
          FillHist("sbllam_Massvspim_BetaGamma",pim_bg,cds_mass);
          FillHist("sbllam_Massvsp_BetaGamma",p_bg,cds_mass);
          FillHist("sbllam_Massvspim_OverBeta2",1.0/pim_beta/pim_beta,cds_mass);
          FillHist("sbllam_Massvsp_OverBeta2",1.0/p_beta/p_beta,cds_mass);
          FillHist("pim_Momentum_sbllam",pim_mom);
          FillHist("p_Momentum_sbllam",p_mom);
        }
        if(sbulamll<cds_mass&&cds_mass<sbulamul){
          FillHist(Form("sbulam_Mass"),cds_mass);
          FillHist(Form("sbulam_Momentum"),cds_mom);
          FillHist(Form("sbulam_CosT"),cds_cost);
          FillHist("sbulam_FL",cds_fl);
          FillHist("sbulam_VDIS",pro->vdis());
          FillHist("sbulam_VBDIS",pro->vbdis());
          FillHist("sbulam_PBDCA",pro->pbdca());
          FillHist("sbulam_Massvspim_BetaGamma",pim_bg,cds_mass);
          FillHist("sbulam_Massvsp_BetaGamma",p_bg,cds_mass);
          FillHist("sbulam_Massvspim_OverBeta2",1.0/pim_beta/pim_beta,cds_mass);
          FillHist("sbulam_Massvsp_OverBeta2",1.0/p_beta/p_beta,cds_mass);
          FillHist("pim_Momentum_sbulam",pim_mom);
          FillHist("p_Momentum_sbulam",p_mom);
        }
      }
    }
  }

  delete particle;
  return true;

}

bool MyAnalysisLambda::FillHist(TString name, double val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1&&weight>0){ 
    for(int i=0; i<weight; i++){
      h1 -> Fill(val1,1);
    }
    return true;
  }
  else 
    return false;
}

bool MyAnalysisLambda::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisLambda::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisLambda::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisLambda::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4976;
  double sigma = 0.0069;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-8*sigma; sblk0ul = mean-4*sigma;
  sbuk0ll = mean+4*sigma; sbuk0ul = mean+8*sigma;
  /* Lambda mass cut */ 
  mean  = 1.11543;
  sigma = 0.0020;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-8*sigma; sbllamul = mean-4*sigma;
  sbulamll = mean+4*sigma; sbulamul = mean+8*sigma;
}

bool MyAnalysisLambda::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisLambda::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaLambda");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
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
  new TH2F( "CDS_Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterFindMass", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_PiPlus", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_PiMinus", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Kaon", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Proton", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Deuteron", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Other", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_All", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterELossCorr", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Fiducial", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_AfterFindMass", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_PiPlus", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_PiMinus", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Kaon", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Proton", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Deuteron", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Other", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_All", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_AfterELossCorr", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Fiducial", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  TString pname[4]  = {"ppim","lam","sbllam","sbulam"};
  TString pname2[4] = {"p#pi^{-}","#Lambda","#Lambda_{lower}","#Lambda_{upper}"};
  for(int ip=0; ip<4; ip++){
    new TH1F( Form("%s_Mass",pname[ip].Data()), Form("Invariant Mass (%s);IM(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 2000, 1.0, 2.0);
    new TH1F( Form("%s_Momentum",pname[ip].Data()), Form("Momentum (%s);Momentum (%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 1.0);
    new TH1F( Form("%s_CosT",pname[ip].Data()), Form("cos#theta (%s);cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5 );
    new TH1F( Form("%s_FL",pname[ip].Data()), Form("Flight length/#beta #gamma (%s);Flight length/#beta #gamma (%s) (cm);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, 0.0, 30.0);
    new TH1F( Form("%s_VDIS",pname[ip].Data()), Form("VDIS (%s);VDIS(%s) (cm);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, 0.0, 30.0);
    new TH1F( Form("%s_VBDIS",pname[ip].Data()), Form("VBDIS (%s);VBDIS(%s) (cm);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, 0.0, 30.0);
    new TH1F( Form("%s_PBDCA",pname[ip].Data()), Form("PBDCA (%s);PBDCA(%s) (cm);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, 0.0, 30.0);
    new TH2F( Form("%s_Massvspim_BetaGamma",pname[ip].Data()), Form("Invariant Mass (%s) vs. #beta #gamma _{#pi^{-}} ;#beta #gamma _{#pi^{-}};IM(%s) (GeV/c^{2})",pname2[ip].Data(),pname2[ip].Data()), 200, 0.0, 2.0, 2000, 1.0, 2.0);
    new TH2F( Form("%s_Massvsp_BetaGamma",pname[ip].Data()), Form("Invariant Mass (%s) vs. #beta #gamma _{p} ;#beta #gamma _{p};IM(%s) (GeV/c^{2})",pname2[ip].Data(),pname2[ip].Data()), 200, 0.0, 2.0, 2000, 1.0, 2.0);

    new TH2F( Form("%s_Massvspim_OverBeta2",pname[ip].Data()), Form("Invariant Mass (%s) vs. 1/#beta^{2}_{#pi^{-}} ;1/#beta^{2}_{#pi^{-}};IM(%s) (GeV/c^{2})",pname2[ip].Data(),pname2[ip].Data()), 600, 0.0, 6.0, 2000, 1.0, 2.0);
    new TH2F( Form("%s_Massvsp_OverBeta2",pname[ip].Data()), Form("Invariant Mass (%s) vs. 1/#beta^{2}_{p} ;1/#beta^{2}_{p};IM(%s) (GeV/c^{2})",pname2[ip].Data(),pname2[ip].Data()), 600, 0.0, 6.0, 2000, 1.0, 2.0);
  }
  new TH1F( Form("pim_Momentum"), Form("Momentum (#pi^{-});Momentum (#pi^{-}) (GeV/c);Counts"), 2000, 0.0, 1.0);
  new TH1F( Form("p_Momentum"), Form("Momentum (p);Momentum (p) (GeV/c);Counts"), 2000, 0.0, 1.0);
  new TH1F( Form("pim_Momentum_lam"), Form("Momentum (#pi^{-});Momentum (#pi^{-}) (GeV/c);Counts"), 2000, 0.0, 1.0);
  new TH1F( Form("p_Momentum_lam"), Form("Momentum (p);Momentum (p) (GeV/c);Counts"), 2000, 0.0, 1.0);
  new TH1F( Form("pim_Momentum_sbllam"), Form("Momentum (#pi^{-});Momentum (#pi^{-}) (GeV/c);Counts"), 2000, 0.0, 1.0);
  new TH1F( Form("p_Momentum_sbllam"), Form("Momentum (p);Momentum (p) (GeV/c);Counts"), 2000, 0.0, 1.0);
  new TH1F( Form("pim_Momentum_sbulam"), Form("Momentum (#pi^{-});Momentum (#pi^{-}) (GeV/c);Counts"), 2000, 0.0, 1.0);
  new TH1F( Form("p_Momentum_sbulam"), Form("Momentum (p);Momentum (p) (GeV/c);Counts"), 2000, 0.0, 1.0);

  new TH2F( "Vertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "Vertex_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH1F( "Vertex_Z_XYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "Vertex_XY_Fiducial", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_Fiducial", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_Fiducial", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z_Fiducial", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "Vertex_XY_FiducialZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH1F( "Vertex_Z_FiducialXYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );

  return true;
}
