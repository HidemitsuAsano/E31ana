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

bool MyAnalysis2Track::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{

  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  /* Trigger selection */
  //if(particle->nCDS()!=2) return true;

  /* Event selection */
  if(particle->nCDS()!=2) return true;
  if(particle->nProduct()!=1) return true;

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
      FillHist(Form("Trigger%02d_TDC",i),val);
    }
  }

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    FillHist("Vertex_XY",cds->vbeam().X(),cds->vbeam().Y(),1);
    FillHist("Vertex_XZ",cds->vbeam().X(),cds->vbeam().Z(),1);
    FillHist("Vertex_YZ",cds->vbeam().Y(),cds->vbeam().Z(),1);
    FillHist("Vertex_Z",cds->vbeam().Z());
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial){
        fcds=it;
        FillHist("VertexwithTRG_XY",cds->vbeam().X(),cds->vbeam().Y(),1);
        FillHist("VertexwithTRG_XZ",cds->vbeam().X(),cds->vbeam().Z(),1);
        FillHist("VertexwithTRG_YZ",cds->vbeam().Y(),cds->vbeam().Z(),1);
        FillHist("VertexwithTRG_Z",cds->vbeam().Z());
      }
    }
  }
  if(fcds==-1) return true; /* there is no vertex in fiducial volume */

  // ######################### //
  // Checking for CDS particle //
  // ######################### //
  TString cdsstr[9] = {"#pi^{+}","p","d","t","He","#pi^{-}","K^{-}","Other","Default"};
  pCDS* cds1 = particle->cdsi(fcds);
  pCDS* cds2 = particle->cdsi(1-fcds);
  FillHist("CDS_NumOfParticle",particle->nCDS());
  FillHist("CDS_Particle",cds1->pid(),cds2->pid());
  FillHist("CDS_Mass2vsMomentum",cds1->mass2(),cdstrackMan->GoodTrack(cds1->id())->Momentum(),1);

  TVector3 ZeroV;
  TLorentzVector Ltgt;
  if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
    Ltgt.SetVectM(ZeroV,pMass);
  }
  if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
    Ltgt.SetVectM(ZeroV,pMass);
  }
  if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
    Ltgt.SetVectM(ZeroV,ThreeHeMass);
  }
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

  // ##################### //
  // CDS Particle analysis //
  // ##################### //
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(GeomTools::GetID(cds->vbeam())==CID_Fiducial){
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
      double cds_pbdca = cds->pbdca();
      FillHist(Form("%s_PBDCA",pname[cds->pid()].Data()),cds_pbdca);
      FillHist(Form("All_PBDCA"),cds_pbdca);
      for(int f=1; f<2; f++){
        FillHist(Form("%s_IMvsMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mass[f],cds_mom[f]);
        FillHist(Form("%s_IMvsCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mass[f],cds_cost[f]);
        FillHist(Form("%s_MMvsMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f],cds_mmom[f]);
        FillHist(Form("%s_MMvsCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f],cds_mcost[f]);
      }
    }
  }
  // ##################### //
  // CDS Product analysis //
  // ##################### //
  for(int it=0; it<particle->nProduct(); it++){
    pCDS* pro = particle->product(it);
    if(pro->pbdca()>20.0){ continue; }
    if(pro->vbdis()>3.0){ continue; }
    if(GeomTools::GetID(pro->vbeam())==CID_Fiducial){
      TLorentzVector Lpro = pro->GetLorentzVector();
      TLorentzVector cLpro = pro->GetLorentzVector(); cLpro.Boost(-boost);
      TString frame[2] = {"Lab","CM"};
      TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
      double pro_mass[2]  = { (Lpro).M()                    , (cLpro).M()};
      double pro_mom[2]   = { (Lpro).Vect().Mag()           , (cLpro).Vect().Mag()};
      double pro_cost[2]  = { (Lpro).Vect().CosTheta()      , (cLpro).Vect().CosTheta()};
      double pro_phi[2]   = { (Lpro).Vect().Phi()           , (cLpro).Vect().Phi()};
      double pro_mmass[2]   = { (Ltgt+Lbeam-Lpro).M()              , (cLtgt+cLbeam-cLpro).M()};
      double pro_mmom[2]    = { (Ltgt+Lbeam-Lpro).Vect().Mag()     , (cLtgt+cLbeam-cLpro).Vect().Mag()};
      double pro_mcost[2]   = { (Ltgt+Lbeam-Lpro).Vect().CosTheta(), (cLtgt+cLbeam-cLpro).Vect().CosTheta()};
      double pro_pbdca = pro->pbdca();
      double pro_vbdis = pro->vbdis();
      double pro_vdis  = pro->vdis();
      // pi+ and pi-
      if(pro->comb()==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
        for(int f=1; f<2; f++){
          FillHist(Form("pipi_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
          FillHist(Form("pipi_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
          FillHist(Form("pKpipi_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
          FillHist(Form("pKpipi_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
        }
        FillHist(Form("pipi_PBDCA"),pro_pbdca);
        FillHist(Form("pipi_VBDIS"),pro_vbdis);
        FillHist(Form("pipi_VDIS"),pro_vdis);
      }
      // p and pi-
      if(pro->comb()==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        for(int f=1; f<2; f++){
          FillHist(Form("ppim_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
          FillHist(Form("ppim_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
          FillHist(Form("pKppim_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
          FillHist(Form("pKppim_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
        }
        FillHist(Form("ppim_PBDCA"),pro_pbdca);
        FillHist(Form("ppim_VBDIS"),pro_vbdis);
        FillHist(Form("ppim_VDIS"),pro_vdis);
      }
      // K- and pi+
      if(pro->comb()==pow(2,CDS_Kaon)+pow(2,CDS_PiPlus)){
        for(int f=1; f<2; f++){
          FillHist(Form("kpip_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
          FillHist(Form("kpip_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
          FillHist(Form("pKkpip_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
          FillHist(Form("pKkpip_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
        }
        FillHist(Form("kpip_PBDCA"),pro_pbdca);
        FillHist(Form("kpip_VBDIS"),pro_vbdis);
        FillHist(Form("kpip_VDIS"),pro_vdis);
      }
      // K- and p
      if(pro->comb()==pow(2,CDS_Kaon)+pow(2,CDS_Proton)){
        for(int f=1; f<2; f++){
          FillHist(Form("kp_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
          FillHist(Form("kp_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
          FillHist(Form("pKkp_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
          FillHist(Form("pKkp_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
        }
        FillHist(Form("kp_PBDCA"),pro_pbdca);
        FillHist(Form("kp_VBDIS"),pro_vbdis);
        FillHist(Form("kp_VDIS"),pro_vdis);
      }
    }
  }

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

  // CDS Particles
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

  std::cout << "Define Histograms for CDS particles" << std::endl;
  new TH1F( "CDS_NumOfParticle", "Number of CDS tracks", 10, 0, 10);
  TH2F* h2;
  h2 = new TH2F("CDS_Particle","Detected particles by CDS;Particle1;Particle2",8,0,8,8,0,8);
  TString cdsstr[9] = {"#pi^{+}","p","d","t","He","#pi^{-}","K^{-}","Other","Default"};
  h2->GetYaxis()->SetBinLabel(1,cdsstr[0].Data());
  h2->GetYaxis()->SetBinLabel(2,cdsstr[1].Data());
  h2->GetYaxis()->SetBinLabel(3,cdsstr[2].Data());
  h2->GetYaxis()->SetBinLabel(4,cdsstr[3].Data());
  h2->GetYaxis()->SetBinLabel(5,cdsstr[4].Data());
  h2->GetYaxis()->SetBinLabel(6,cdsstr[5].Data());
  h2->GetYaxis()->SetBinLabel(7,cdsstr[6].Data());
  h2->GetYaxis()->SetBinLabel(8,cdsstr[7].Data());
  h2->GetXaxis()->SetBinLabel(1,cdsstr[0].Data());
  h2->GetXaxis()->SetBinLabel(2,cdsstr[1].Data());
  h2->GetXaxis()->SetBinLabel(3,cdsstr[2].Data());
  h2->GetXaxis()->SetBinLabel(4,cdsstr[3].Data());
  h2->GetXaxis()->SetBinLabel(5,cdsstr[4].Data());
  h2->GetXaxis()->SetBinLabel(6,cdsstr[5].Data());
  h2->GetXaxis()->SetBinLabel(7,cdsstr[6].Data());
  h2->GetXaxis()->SetBinLabel(8,cdsstr[7].Data());
  // CDS Particles
  new TH2F( "CDS_Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 2000, -5, 5, 4000, -2, 18 );

  // CDS 1 particle
  TString pname[8]  = {"pip","p","d","t","he","pim","k","o"};
  TString pname2[8] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o"};
  for(int ip=0; ip<8; ip++){
    new TH2F( Form("%s_IMvsMomentum_CM",pname[ip].Data()), Form("Invariant mass (%s) vs Momentum;IM(%s) (GeV/c^{2});Momentum (GeV/c)",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
    new TH2F( Form("%s_IMvsCosT_CM",pname[ip].Data()), Form("Invariant mass (%s) vs cos#theta;IM(%s) (GeV/c^{2});cos#theta",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
    new TH2F( Form("%s_MMvsMomentum_CM",pname[ip].Data()), Form("Missing mass (%s) vs Momentum;MM(%s) (GeV/c^{2});Momentum (GeV/c)",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
    new TH2F( Form("%s_MMvsCosT_CM",pname[ip].Data()), Form("Missing mass (%s) vs cos#theta;MM(%s) (GeV/c^{2});cos#theta",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
    new TH1F( Form("%s_PBDCA",pname[ip].Data()),Form("DCA between beam (%s);DCA (cm);Counts",pname2[ip].Data()),1000,0.0,10.0);
  }
  new TH1F( "All_PBDCA","DCA between beam (All);DCA (cm);Counts",1000,0.0,10.0);

  // Product
  new TH2F( "pipi_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 5000, 0.0, 5.0, 3000, 0.0, 3.0 );
  new TH2F( "pipi_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 5000, 0.0, 5.0, 3000, -1.5, 1.5 );
  new TH2F( "pipi_MMvsMomentum_CM", "^{3}He(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{^{3}He}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 5000, 0.0, 5.0, 3000, 0.0, 3.0);
  new TH2F( "pipi_MMvsCosT_CM", "^{3}He(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{^{3}He}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 5000, 0.0, 5.0, 3000, -1.5, 1.5);
  new TH1F( "pipi_VDIS","Distance between tow vertexes (#pi^{+}#pi^{-});Distance (cm);Counts",1000,0.0,10.0);
  new TH1F( "pipi_VBDIS","Distance between beam and vertex (#pi^{+}#pi^{-});Distance (cm);Counts",100,0.0,10.0);
  new TH1F( "pipi_PBDCA","DCA between beam and product (#pi^{+}#pi^{-});DCA (cm);Counts",1000,0.0,10.0);
  new TH2F( "ppim_IMvsMomentum_CM", "Invariant mass (p#pi^{-}) vs Momentum;IM(p#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 5000, 0.0, 5.0, 3000, 0.0, 3.0 );
  new TH2F( "ppim_IMvsCosT_CM", "Invariant mass (p#pi^{-}) vs cos#theta;IM(p#pi^{-}) (GeV/c^{2});cos#theta", 5000, 0.0, 5.0, 3000, -1.5, 1.5 );
  new TH2F( "ppim_MMvsMomentum_CM", "^{3}He(K^{-},p#pi^{-})X missing mass vs. missing momentum;MM_{^{3}He}(p#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 5000, 0.0, 5.0, 3000, 0.0, 3.0);
  new TH2F( "ppim_MMvsCosT_CM", "^{3}He(K^{-},p#pi^{-})X missing mass vs. cos#theta;MM_{^{3}He}(p#pi^{-}) (GeV/c^{2});cos#theta", 5000, 0.0, 5.0, 3000, -1.5, 1.5);
  new TH1F( "ppim_VDIS","Distance between tow vertexes (p#pi^{-});Distance (cm);Counts",1000,0.0,10.0);
  new TH1F( "ppim_VBDIS","Distance between beam and vertex (p#pi^{-});Distance (cm);Counts",100,0.0,10.0);
  new TH1F( "ppim_PBDCA","DCA between beam and product (p#pi^{-});DCA (cm);Counts",1000,0.0,10.0);
  new TH2F( "kpip_IMvsMomentum_CM", "Invariant mass (K^{-}#pi^{+}) vs Momentum;IM(K^{-}#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 5000, 0.0, 5.0, 3000, 0.0, 3.0 );
  new TH2F( "kpip_IMvsCosT_CM", "Invariant mass (K^{-}#pi^{+}) vs cos#theta;IM(K^{-}#pi^{+}) (GeV/c^{2});cos#theta", 5000, 0.0, 5.0, 3000, -1.5, 1.5 );
  new TH2F( "kpip_MMvsMomentum_CM", "^{3}He(K^{-},K^{-}#pi^{+})X missing mass vs. missing momentum;MM_{^{3}He}(K^{-}#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 5000, 0.0, 5.0, 3000, 0.0, 3.0);
  new TH2F( "kpip_MMvsCosT_CM", "^{3}He(K^{-},K^{-}#pi^{+})X missing mass vs. cos#theta;MM_{^{3}He}(K^{-}#pi^{+}) (GeV/c^{2});cos#theta", 5000, 0.0, 5.0, 3000, -1.5, 1.5);
  new TH1F( "kpip_VDIS","Distance between tow vertexes (K^{-}#pi^{+});Distance (cm);Counts",1000,0.0,10.0);
  new TH1F( "kpip_VBDIS","Distance between beam and vertex (K^{-}#pi^{+});Distance (cm);Counts",100,0.0,10.0);
  new TH1F( "kpip_PBDCA","DCA between beam and product (K^{-}#pi^{+});DCA (cm);Counts",1000,0.0,10.0);
  new TH2F( "kp_IMvsMomentum_CM", "Invariant mass (K^{-}p) vs Momentum;IM(K^{-}p) (GeV/c^{2});Momentum (GeV/c)", 5000, 0.0, 5.0, 3000, 0.0, 3.0 );
  new TH2F( "kp_IMvsCosT_CM", "Invariant mass (K^{-}p) vs cos#theta;IM(K^{-}p) (GeV/c^{2});cos#theta", 5000, 0.0, 5.0, 3000, -1.5, 1.5 );
  new TH2F( "kp_MMvsMomentum_CM", "^{3}He(K^{-},K^{-}p)X missing mass vs. missing momentum;MM_{^{3}He}(K^{-}p) (GeV/c^{2});Missing momentum (GeV/c)", 5000, 0.0, 5.0, 3000, 0.0, 3.0);
  new TH2F( "kp_MMvsCosT_CM", "^{3}He(K^{-},K^{-}p)X missing mass vs. cos#theta;MM_{^{3}He}(K^{-}p) (GeV/c^{2});cos#theta", 5000, 0.0, 5.0, 3000, -1.5, 1.5);
  new TH1F( "kp_VDIS","Distance between tow vertexes (K^{-}p);Distance (cm);Counts",1000,0.0,10.0);
  new TH1F( "kp_VBDIS","Distance between beam and vertex (K^{-}p);Distance (cm);Counts",100,0.0,10.0);
  new TH1F( "kp_PBDCA","DCA between beam and product (K^{-}p);DCA (cm);Counts",1000,0.0,10.0);
  return true;
}
