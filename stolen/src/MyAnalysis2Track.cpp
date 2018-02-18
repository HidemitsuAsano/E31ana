// MyAnalysis2Track.cpp

#include "MyAnalysis2Track.h"

MyAnalysis2Track::MyAnalysis2Track(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysis2Track::~MyAnalysis2Track()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysis2Track::Clear()
{
}

bool MyAnalysis2Track::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* _particle)
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

  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cdstrack = particle->cdsi(it);
    double mass2 = cdstrack->mass2();
    double beta = cdstrack->beta();
    double mom = cdstrack->mom();
    TVector3 vtxb = cdstrack->vbeam();
    TVector3 vtxh = cdstrack->vcdc();
    FillHist("Vertex_XY", vtxb.X(), vtxb.Y());
    FillHist("Vertex_ZX", vtxb.Z(), vtxb.X());
    FillHist("Vertex_ZY", vtxb.Z(), vtxb.Y());
    FillHist("Vertex_Z", vtxb.Z());
    if(-0.1<vtxb.Z()&&vtxb.Z()<0.1){
      FillHist("Vertex_XY_ZCut", vtxb.X(), vtxb.Y());
    }
    if(-0.1<vtxb.X()&&vtxb.X()<0.1&&-0.1<vtxb.Y()&&vtxb.Y()<0.1){
      FillHist("Vertex_Z_XYCut", vtxb.Z());
    }

    if(GeomTools::GetID(vtxb)==CID_Fiducial){
      FillHist("Vertex_XY_Fiducial", vtxb.X(), vtxb.Y());
      FillHist("Vertex_ZX_Fiducial", vtxb.Z(), vtxb.X());
      FillHist("Vertex_ZY_Fiducial", vtxb.Z(), vtxb.Y());
      FillHist("Vertex_Z_Fiducial", vtxb.Z());
      if(-0.1<vtxb.Z()&&vtxb.Z()<0.1){
        FillHist("Vertex_XY_FiducialZCut", vtxb.X(), vtxb.Y());
      }
      if(-0.1<vtxb.X()&&vtxb.X()<0.1&&-0.1<vtxb.Y()&&vtxb.Y()<0.1){
        FillHist("Vertex_Z_FiducialXYCut", vtxb.Z());
      }
      FillHist("CDS_OverbetavsMomentum_Fiducial",1./beta,mom);
      FillHist("CDS_Mass2vsMomentum_Fiducial",mass2,mom);
      //FillHist("CDS_dEvsMomentum_Fiducial",de,mom);
    }
  }

  // ##################### //
  // CDS Particle analysis //
  // ##################### //
  int ncds[6] = {0,0,0,0,0,0};
  ncds[0] = particle->nPiplus();
  ncds[1] = particle->nPiminus();
  ncds[2] = particle->nKaon();
  ncds[3] = particle->nProton();
  ncds[4] = particle->nDeuteron();
  ncds[5] = particle->nTriton();
  ncds[5] += particle->nHelium3();
  ncds[5] += particle->nOther();
  TString cdsstr[6] = {"#pi^{+}","#pi^{-}","K^{-}","p","d","Other"};

  FillHist("CDS_NumOfParticle",particle->nCDS());
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(GeomTools::GetID(cds->vbeam())==CID_Fiducial){
      double cds_vdis = cds->vdis();
      TVector3 vtxb = cds->vbeam();
      TVector3 ZeroV;
      TLorentzVector Ltgt;
      if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
        Ltgt.SetVectM(ZeroV,pMass);
      }
      if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
        Ltgt.SetVectM(ZeroV,dMass);
      }
      if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
        Ltgt.SetVectM(ZeroV,ThreeHeMass);
      }
      TLorentzVector Lbeam = beam->GetLorentzVector(vtxb);
      TVector3 boost = (Ltgt+Lbeam).BoostVector();
      TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
      TLorentzVector Lcds = cds->GetLorentzVector();
      TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
      TString frame[2] = {"Lab","CM"};
      TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
      double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
      double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
      double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
      double cds_phi[2]   = { (Lcds).Vect().Phi()           , (cLcds).Vect().Phi()};
      double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
      double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
      double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};

      FillHist(Form("%s_Vertex_XY",pname[cds->pid()].Data()), vtxb.X(), vtxb.Y());
      FillHist(Form("%s_Vertex_ZX",pname[cds->pid()].Data()), vtxb.Z(), vtxb.X());
      FillHist(Form("%s_Vertex_ZY",pname[cds->pid()].Data()), vtxb.Z(), vtxb.Y());
      FillHist(Form("%s_Vertex_Z",pname[cds->pid()].Data()), vtxb.Z());
      FillHist(Form("%s_VDIS",pname[cds->pid()].Data()),cds_vdis);
      FillHist(Form("All_VDIS"),cds_vdis);
      for(int f=1; f<2; f++){
        FillHist(Form("%s_Mass_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mass[f]);
        FillHist(Form("%s_Momentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mom[f]);
        FillHist(Form("%s_CosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_cost[f]);
        FillHist(Form("%s_MMass_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f]);
        FillHist(Form("%s_MMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmom[f]);
        FillHist(Form("%s_MCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mcost[f]);
      }
    }
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
    }
  }

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nProduct(); it++){
    pCDS* cds = particle->product(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
    }
  }
  if(fcds==-1) return true; /* there is no vertex in fiducial volume */
  // ############### //
  // NC hit decision //
  // ############### //
  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particle->nNC(); it++){
    pNC* nc = particle->nc(it);
    nc->CalcMom(beam,vertex);
    FillHist("FWD_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
    FillHist("FWD_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("FWD_TOFvsMomentum",nc->tof(),nc->mom().Mag());
    FillHist("FWD_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
    double seg = (nc->seg()-1)%16+1;
    double lay = (nc->seg()-1)/16+1;
    FillHist("FWD_HitSegment",seg,lay);
    if(nc->pid()==F_Neutron) {
      if(nc->energy()>8.0 && (time>9998||time>nc->time())){
        time = nc->time();
        fnc = it;
      }
    }
  }
  pNC* nc = 0;
  if(fnc!=-1){ nc = particle->nc(fnc); }

  for(int it=0; it<particle->nProduct(); it++){
    pCDS* cds = particle->product(it);
    if(cds->pbdca()>20.0){ continue; }
    if(cds->vbdis()>3.0){ continue; }
    if(GeomTools::GetID(cds->vbeam())==CID_Fiducial){
      double cds_vdis = cds->vdis();
      TVector3 vtxb = cds->vbeam();
      TVector3 ZeroV;
      TLorentzVector Ltgt;
      if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
        Ltgt.SetVectM(ZeroV,pMass);
      }
      if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
        Ltgt.SetVectM(ZeroV,dMass);
      }
      if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
        Ltgt.SetVectM(ZeroV,ThreeHeMass);
      }
      TLorentzVector Lbeam = beam->GetLorentzVector(vtxb);
      TVector3 boost = (Ltgt+Lbeam).BoostVector();
      TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
      TLorentzVector Lcds = cds->GetLorentzVector();
      TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
      TString frame[2] = {"Lab","CM"};
      TString pname[6] = {"pipi","ppim","kpip","pk","dk","pro_o"};
      double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
      double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
      double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
      double cds_phi[2]   = { (Lcds).Vect().Phi()           , (cLcds).Vect().Phi()};
      double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
      double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
      double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};
      TLorentzVector Ln;
      TLorentzVector cLn;
      if(nc!=0){
        Ln  = nc ->GetLorentzVector();
        cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
      }
      else{
        Ln.SetVectM(ZeroV,0);
        cLn.SetVectM(ZeroV,0);
      }
      double ncds_mass[2]  = { (Lcds+Ln).M()                           , (cLcds+cLn).M()};
      double ncds_mom[2]   = { (Lcds+Ln).Vect().Mag()                  , (cLcds+cLn).Vect().Mag()};
      double ncds_cost[2]  = { (Lcds+Ln).Vect().CosTheta()             , (cLcds+cLn).Vect().CosTheta()};
      double ncds_phi[2]   = { (Lcds+Ln).Vect().Phi()                  , (cLcds+cLn).Vect().Phi()};
      double ncds_mmass[2] = { (Ltgt+Lbeam-(Lcds+Ln)).M()              , (cLtgt+cLbeam-(cLcds+cLn)).M()};
      double ncds_mmom[2]  = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds+cLn)).Vect().Mag()};
      double   ncds_mcost[2] = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds+cLn)).Vect().CosTheta()};

      FillHist("PVertex_XY_Fiducial", vtxb.X(), vtxb.Y());
      FillHist("PVertex_ZX_Fiducial", vtxb.Z(), vtxb.X());
      FillHist("PVertex_ZY_Fiducial", vtxb.Z(), vtxb.Y());
      FillHist("PVertex_Z_Fiducial", vtxb.Z());
      if(-0.1<vtxb.Z()&&vtxb.Z()<0.1){
        FillHist("PVertex_XY_FiducialZCut", vtxb.X(), vtxb.Y());
      }
      if(-0.1<vtxb.X()&&vtxb.X()<0.1&&-0.1<vtxb.Y()&&vtxb.Y()<0.1){
        FillHist("PVertex_Z_FiducialXYCut", vtxb.Z());
      }

      int combid = 5;
      int comb = cds->comb();
      if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
        combid = 0;
      }
      else if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        combid = 1;
      }
      else if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_Kaon)){
        combid = 2;
      }
      else if(comb==pow(2,CDS_Proton)+pow(2,CDS_Kaon)){
        combid = 3;
      }
      else if(comb==pow(2,CDS_Deuteron)+pow(2,CDS_Kaon)){
        combid = 4;
      }
      for(int f=1; f<2; f++){
        FillHist(Form("%s_Mass_%s",pname[combid].Data(),frame[f].Data()),cds_mass[f]);
        FillHist(Form("%s_Momentum_%s",pname[combid].Data(),frame[f].Data()),cds_mom[f]);
        FillHist(Form("%s_CosT_%s",pname[combid].Data(),frame[f].Data()),cds_cost[f]);
        FillHist(Form("%s_MMass_%s",pname[combid].Data(),frame[f].Data()),cds_mmass[f]);
        FillHist(Form("%s_MMomentum_%s",pname[combid].Data(),frame[f].Data()),cds_mmom[f]);
        FillHist(Form("%s_MCosT_%s",pname[combid].Data(),frame[f].Data()),cds_mcost[f]);
      }
      for(int f=1; f<2; f++){
        if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
          if(0.48531<cds_mass[1]&&cds_mass[1]<509.11){
            FillHist(Form("k0n_Mass_%s",frame[f].Data()),ncds_mass[f]);
            FillHist(Form("k0n_Momentum_%s",frame[f].Data()),ncds_mom[f]);
            FillHist(Form("k0n_CosT_%s",frame[f].Data()),ncds_cost[f]);
            FillHist(Form("k0n_MMass_%s",frame[f].Data()),ncds_mmass[f]);
            FillHist(Form("k0n_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
            FillHist(Form("k0n_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
          }
        }
        if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
          if(1.11145<cds_mass[1]&&cds_mass[1]<1.111897){
            FillHist(Form("ln_Mass_%s",frame[f].Data()),ncds_mass[f]);
            FillHist(Form("ln_Momentum_%s",frame[f].Data()),ncds_mom[f]);
            FillHist(Form("ln_CosT_%s",frame[f].Data()),ncds_cost[f]);
            FillHist(Form("ln_MMass_%s",frame[f].Data()),ncds_mmass[f]);
            FillHist(Form("ln_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
            FillHist(Form("ln_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
          }
        }
        if(comb==pow(2,CDS_Proton)+pow(2,CDS_Kaon)){
          FillHist(Form("pkn_Mass_%s",frame[f].Data()),ncds_mass[f]);
          FillHist(Form("pkn_Momentum_%s",frame[f].Data()),ncds_mom[f]);
          FillHist(Form("pkn_CosT_%s",frame[f].Data()),ncds_cost[f]);
          FillHist(Form("pkn_MMass_%s",frame[f].Data()),ncds_mmass[f]);
          FillHist(Form("pkn_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
          FillHist(Form("pkn_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
        }
      }
    }
  }

  delete particle;
  return true;

}

bool MyAnalysis2Track::FillHist(TString name, double val1, int weight)
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

bool MyAnalysis2Track::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysis2Track::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysis2Track::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysis2Track::CutCondition()
{
}

bool MyAnalysis2Track::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysis2Track::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_ana2Track");

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
  TString pname[17]  = {"pip","p","d","t","he","pim","k","o","pipi","ppim","kpip","pk","dk","proo","pkn","k0n","ln"};
  TString pname2[17] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o","#pi^{+}#pi^{-}","p#pi^{-}","#pi^{+}K^{-}","pK^{-}","dK^{-}","o","pK^{-}n","K^{0}n","#Lambda n"};
  for(int ip=0; ip<17; ip++){
    new TH1F( Form("%s_VDIS",pname[ip].Data()), Form("VDIS of %s;VDIS (cm);Coutns",pname2[ip].Data()), 100, 0, 10);
    new TH2F( Form("%s_Vertex_XY",pname[ip].Data()), Form("Vertex XY plane (%s);x (cm);y (cm)",pname2[ip].Data()), 600, -15, 15, 600, -15, 15 );
    new TH2F( Form("%s_Vertex_ZX",pname[ip].Data()), Form("Vertex ZX plane (%s);z (cm);x (cm)",pname2[ip].Data()), 1200, -30, 30, 600, -15, 15 );
    new TH2F( Form("%s_Vertex_ZY",pname[ip].Data()), Form("Vertex ZY plane (%s);z (cm);y (cm)",pname2[ip].Data()), 1200, -30, 30, 600, -15, 15 );
    new TH1F( Form("%s_Vertex_Z",pname[ip].Data()), Form("Vertex Z axis (%s);z (cm);Counts",pname2[ip].Data()), 1200, -30, 30 );
    new TH1F( Form("%s_Mass_CM",pname[ip].Data()), Form("Invariant Mass (%s);IM(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%s_Momentum_CM",pname[ip].Data()), Form("Momentum (%s);Momentum (%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%s_CosT_CM",pname[ip].Data()), Form("cos#theta (%s);cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5 );
    new TH1F( Form("%s_MMass_CM",pname[ip].Data()), Form("Missing mass (%s);Missing Mass(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%s_MMomentum_CM",pname[ip].Data()), Form("Missing Momentum (%s);Missing Momentum(%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%s_MCosT_CM",pname[ip].Data()), Form("Missing cos#theta (%s);Missing cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5);
  }
  new TH1F( Form("All_VDIS"), Form("VDIS of All;VDIS (cm);Coutns"), 100, 0, 10);
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

  // Vertex
  std::cout << "Define Histograms for Vertex" << std::endl;
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

  new TH2F( "PVertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "PVertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "PVertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "PVertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "PVertex_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH1F( "PVertex_Z_XYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "PVertex_XY_Fiducial", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "PVertex_ZX_Fiducial", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "PVertex_ZY_Fiducial", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "PVertex_Z_Fiducial", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "PVertex_XY_FiducialZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH1F( "PVertex_Z_FiducialXYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );

  return true;
}
