// MyAnalysisdKppipi.cpp

#include "MyAnalysisdKppipi.h"

MyAnalysisdKppipi::MyAnalysisdKppipi(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisdKppipi::~MyAnalysisdKppipi()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisdKppipi::Clear()
{
}

bool MyAnalysisdKppipi::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  int ievent=0;
  rtFile->cd();

  FillHist("EventNumber",ievent); ievent++; /* All Events */

  if(particle->nBeam()!=1) return false;
  FillHist("EventNumber",ievent); ievent++; /* Beam Analysis */
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  if(!header->IsTrig(Trig_KCDH3)){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* KCDH3 Trigger */

  int MulCDH=0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()) MulCDH++;
  }
  FillHist("CDH_Multiplicity",MulCDH);
  if(MulCDH!=3){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* 3 hits in CDH */

  if(cdstrackMan->nGoodTrack()!=3) return true;
  FillHist("EventNumber",ievent); ievent++; /* 3 Good Tracks in CDC */

  /* Event selection */
  int MulIH=0;
  for(int i=0; i<cdsMan->nIH(); i++){
    HodoscopeLikeHit* hit = cdsMan->IH(i);
    if(hit->CheckRange()) MulIH++;
  }
  FillHist("IH_Multiplicity",MulIH);
  //if(MulIH!=3){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* 3 hits in IH */

  if(particle->nCDS()!=3) return true;
  FillHist("EventNumber",ievent); ievent++; /* No FWD charged hit */

  if(particle->nProton()!=1) return true;
  if(particle->nPiminus()!=2) return true;
  FillHist("EventNumber",ievent); ievent++; /* 1 proton and 2 pions in CDS */

  // ############## //
  // pim1/pim2 decision //
  // ############## //
  TVector3 vertex;
  double vdis = 9999;
  bool FIDUCIAL = false;
  int fcds = -1;
  for(int it=0; it<particle->nPiminus(); it++){
    pCDS* cds = particle->pim(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      fcds=it;
    }
  }
  if(fcds<0){ return true; }

  pCDS* p = particle->proton(0);
  pCDS* pim1 = particle->pim(fcds);
  pCDS* pim2 = particle->pim(1-fcds);

  // ############## //
  // CDC Chi-square //
  // ############## //
  if(p->chi()>30) return true;
  if(pim1->chi()>30) return true;
  if(pim2->chi()>30) return true;
  FillHist("EventNumber",ievent); ievent++; /* CDC Chi-square < 30 */

  // ######################### //
  // 2 Helix calculation       //
  // ######################### //
  pCDS* pip1 = 0;
  pCDS* pip2 = 0;
  {
    for(int it=0; it<particle->nProduct(); it++){
      pCDS* product = particle->product(it);
      int comb = product->comb();
      if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        if(product->daughter1()==pim1->id()||product->daughter2()==pim1->id()){
          pip1 = product;
        }
        else if(product->daughter1()==pim2->id()||product->daughter2()==pim2->id()){
          pip2 = product;
        }
      }
    }
  }
  if(pip1==0||pip2==0){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* 2-track Combination */

  FillHist(Form("pip1_IMvspip2_IM_All"),pip1->mass(),pip2->mass());
  FillHist(Form("pip1_IM2vspip2_IM2_All"),pip1->mass()*pip1->mass(),pip2->mass()*pip2->mass());

  bool pip1flag = false, pip2flag = false;
  pCDS* lam = 0;
  pCDS* nlam = 0;
  pCDS* pim = 0;
  pCDS* npim = 0;
  /* For pip1 */
  double tmpper1[5] = {pip1->mass(),pip1->vdis(),pip1->pbdca(),pim2->vbdis(),0.0};
  double tmppdf1[5] = {0,0,0,0,0};
  /* Input : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  /*Output : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  TrackTools::PDFLambda(tmpper1,tmppdf1);
  tmppdf1[4] = 1.0;
  double tmpper2[5] = {pip2->mass(),pip2->vdis(),pip2->pbdca(),pim1->vbdis(),0.0};
  double tmppdf2[5] = {0,0,0,0,0};
  /* Input : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  /*Output : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  TrackTools::PDFLambda(tmpper2,tmppdf2);
  tmppdf2[4] = 1.0;
  /* All */
  double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]*tmppdf1[2]*tmppdf1[3]);
  double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]*tmppdf2[2]*tmppdf2[3]);
  /* Except DCA_pip */
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[2]*tmppdf1[3]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[2]*tmppdf2[3]);
  ///* Except DCA_Lk */
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]*tmppdf1[3]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]*tmppdf2[3]);
  /* Except DCA_pk */
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]*tmppdf1[2]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]*tmppdf2[2]);
  /* Only Mass */
  //double pdf1 = -TMath::Log(tmppdf1[0]);
  //double pdf2 = -TMath::Log(tmppdf2[0]);
  /* Only DCA_pip */
  //double pdf1 = -TMath::Log(tmppdf1[1]);
  //double pdf2 = -TMath::Log(tmppdf2[1]);
  /* Mass*DCA_pip */
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]);
  /* ALL DCA */
  //double pdf1 = -TMath::Log(tmppdf1[1]*tmppdf1[2]*tmppdf1[3]);
  //double pdf2 = -TMath::Log(tmppdf2[1]*tmppdf2[2]*tmppdf2[3]);

  double pdf=-1;
  double npdf=-1;
  //std::cout << "pdf1 : " << pdf1 <<std::endl;
  //std::cout << "pdf2 : " << pdf2 <<std::endl;
  FillHist("pip1_PDF",pdf1);
  FillHist("pip2_PDF",pdf2);
  FillHist("pip1_PDFvspip2_PDF",pdf1,pdf2);

  double pdfll =  0.0;
  double pdful = 12.0;
  //double pdfll = 12.0;
  //double pdful = 40.0;
  //double pdfll =  0.0;
  //double pdful = 40.0;
  //double pdfll = 40.0;
  //double pdful = 999999999;

  if(pdf1<pdf2){
    FillHist("pip1_PDF_SemiSelected",pdf1);
    FillHist("pip1_PDFvspip2_SemiSelected",pdf1,pdf2);
    FillHist("lam_PDF",pdf1);
    FillHist("nlam_PDF",pdf2);
    if(pdfll<pdf1&&pdf1<pdful){
    FillHist("pip1_PDF_Selected",pdf1);
    pip1flag = true;
    lam = pip1;
    nlam = pip2;
    pim = pim2;
    npim = pim1;
    pdf=pdf1;
    npdf=pdf2;
    }
  }
  else{
    FillHist("pip2_PDF_SemiSelected",pdf2);
    FillHist("pip1_PDFvspip2_SemiSelected",pdf1,pdf2);
    FillHist("lam_PDF",pdf2);
    FillHist("nlam_PDF",pdf1);
    if(pdfll<pdf2&&pdf2<pdful){
    FillHist("pip2_PDF_Selected",pdf2);
    pip2flag = true;
    pip1flag = false;
    lam = pip2;
    nlam = pip1;
    pim = pim1;
    npim = pim2;
    pdf=pdf2;
    npdf=pdf1;
    }
  }
  if(!pip1flag&&!pip2flag){ return true; }
  if(lam==0){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* Lambda Reconstruction */
  FillHist("pip1_PDFvspip2_PDF_Selected",pdf1,pdf2);
  FillHist("lam_PDF_Selected",pdf);
  FillHist("nlam_PDF_Selected",npdf);

  // ################ //
  // Vertex decision1 //
  // ################ //
  if(pip1flag){
    vertex = pim2->vbeam();
    vdis = pim2->vbdis();
    if(GeomTools::GetID(vertex)!=CID_Fiducial){ return true; }
  }
  if(pip2flag){
    vertex = pim1->vbeam();
    vdis = pim1->vbdis();
    if(GeomTools::GetID(vertex)!=CID_Fiducial){ return true; }
  }
  FillHist("EventNumber",ievent); ievent++; /* Vertex Decision 1 */

  // ################ //
  // Vertex decision2 //
  // ################ //
  TVector3 vertex_lam = lam->vbeam();
  if(GeomTools::GetID(vertex_lam)!=CID_Fiducial){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* Vertex Decision 2 */

  // ######################### //
  // Checking for CDS particle //
  // ######################### //
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
  for(int i=0; i<6; i++){
    FillHist("CDS_Particle",cdsstr[i],Form("%d track",particle->nCDS()),ncds[i]);
  }
  TVector3 vtxb;
  /* Protons */
  vtxb = p->vbeam();
  FillHist("p_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p_VBDIS",pim->vbdis());
  /* Pi- 1 */
  vtxb = pim1->vbeam();
  FillHist("pim1_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pim1_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pim1_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pim1_VBDIS",pim1->vbdis());
  /* Pi- 2 */
  vtxb = pim2->vbeam();
  FillHist("pim2_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pim2_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pim2_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pim2_VBDIS",pim2->vbdis());
  /* Pi- 1 and Proton */
  vtxb = pip1->vbeam();
  FillHist("pip1_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pip1_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pip1_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pip1_Mass",pip1->vdis());
  FillHist("pip1_VDIS",pip1->vdis());
  FillHist("pip1_VBDIS",pip1->vbdis());
  FillHist("pip1_PBDCA",pip1->pbdca());
  /* Pi- 2 and Proton */
  vtxb = pip2->vbeam();
  FillHist("pip2_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pip2_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pip2_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pip2_Mass",pip2->vdis());
  FillHist("pip2_VDIS",pip2->vdis());
  FillHist("pip2_VBDIS",pip2->vbdis());
  FillHist("pip2_PBDCA",pip2->pbdca());
  /* 2D Plots */
  FillHist("pim1_VBDISvspim2_VBDIS",pim1->vbdis(),pim2->vbdis());
  FillHist("pip1_VDISvspip2_VDIS",pip1->vdis(),pip2->vdis());
  FillHist("pip1_VBDISvspip2_VBDIS",pip1->vbdis(),pip2->vbdis());
  FillHist("pip1_PBDCAvspip2_PBDCA",pip1->pbdca(),pip2->pbdca());
  /* Pi- */
  vtxb = pim->vbeam();
  FillHist("pim_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pim_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pim_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pim_VBDIS",pim->vbdis());
  /* Lambda */
  vtxb = lam->vbeam();
  FillHist("lam_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("lam_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("lam_Vertex_ZY",vtxb.Z(),vtxb.Y());
  vtxb = lam->vertex();
  FillHist("lam_PVertex_XY",vtxb.X(),vtxb.Y());
  FillHist("lam_PVertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("lam_PVertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("lam_Mass",lam->mass());
  FillHist("lam_VDIS",lam->vdis());
  FillHist("lam_VBDIS",lam->vbdis());
  FillHist("lam_PBDCA",lam->pbdca());
  /* Not Lambda */
  vtxb = nlam->vbeam();
  FillHist("nlam_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("nlam_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("nlam_Vertex_ZY",vtxb.Z(),vtxb.Y());
  vtxb = nlam->vertex();
  FillHist("nlam_PVertex_XY",vtxb.X(),vtxb.Y());
  FillHist("nlam_PVertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("nlam_PVertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("nlam_Mass",nlam->mass());
  FillHist("nlam_VDIS",nlam->vdis());
  FillHist("nlam_VBDIS",nlam->vbdis());
  FillHist("nlam_PBDCA",nlam->pbdca());

  FillHist("lam_PDFvslam_Mass",pdf,lam->mass());
  FillHist("lam_PDFvsnlam_Mass",pdf,nlam->mass());
  FillHist("lam_Massvsnlam_Mass",lam->mass(),nlam->mass());

  // ############ //
  // DCA Decision //
  // ############ //

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }

  TString frame[2] = {"Lab","CM"};
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,dMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

  // Pi+, Pi-, and Pi+ Pi- pair //
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  TLorentzVector Lpim1 = pim1->GetLorentzVector();
  TLorentzVector cLpim1 = pim1->GetLorentzVector(); cLpim1.Boost(-boost);
  TLorentzVector Lpim2 = pim2->GetLorentzVector();
  TLorentzVector cLpim2 = pim2->GetLorentzVector(); cLpim2.Boost(-boost);
  TLorentzVector Lpip1 = pip1->GetLorentzVector();
  TLorentzVector cLpip1 = pip1->GetLorentzVector(); cLpip1.Boost(-boost);
  TLorentzVector Lpip2 = pip2->GetLorentzVector();
  TLorentzVector cLpip2 = pip2->GetLorentzVector(); cLpip2.Boost(-boost);
  TLorentzVector Lp = p->GetLorentzVector();
  TLorentzVector cLp = p->GetLorentzVector(); cLp.Boost(-boost);

  TLorentzVector Llam = lam->GetLorentzVector();
  TLorentzVector cLlam = lam->GetLorentzVector(); cLlam.Boost(-boost);
  TLorentzVector Lnlam = nlam->GetLorentzVector();
  TLorentzVector cLnlam = nlam->GetLorentzVector(); cLnlam.Boost(-boost);

  TLorentzVector Lppipi = Llam + Lpim, cLppipi = cLlam + cLpim;
  TLorentzVector Lmp = Ltgt+Lbeam-Lppipi;
  TLorentzVector cLmp = cLtgt+cLbeam-cLppipi;

  double pim_mass[2]  = { Lpim.M()                           , cLpim.M()};
  double pim_mom[2]   = { Lpim.Vect().Mag()                  , cLpim.Vect().Mag()};
  double pim_cost[2]  = { Lpim.Vect().CosTheta()             , cLpim.Vect().CosTheta()};
  double pim_phi[2]   = { Lpim.Vect().Phi()                  , cLpim.Vect().Phi()};
  double pim_mmass[2] = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]  = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2] = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};

  double p_mass[2]  = { Lp.M()                           , cLp.M()};
  double p_mom[2]   = { Lp.Vect().Mag()                  , cLp.Vect().Mag()};
  double p_cost[2]  = { Lp.Vect().CosTheta()             , cLp.Vect().CosTheta()};
  double p_phi[2]   = { Lp.Vect().Phi()                  , cLp.Vect().Phi()};
  double p_mmass[2] = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
  double p_mmom[2]  = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
  double p_mcost[2] = { (Ltgt+Lbeam-Lp).Vect().CosTheta(), (cLtgt+cLbeam-cLp).Vect().CosTheta()};

  double pim1_mass[2]  = { Lpim1.M()                           , cLpim1.M()};
  double pim1_mom[2]   = { Lpim1.Vect().Mag()                  , cLpim1.Vect().Mag()};
  double pim1_cost[2]  = { Lpim1.Vect().CosTheta()             , cLpim1.Vect().CosTheta()};
  double pim1_phi[2]   = { Lpim1.Vect().Phi()                  , cLpim1.Vect().Phi()};
  double pim1_mmass[2] = { (Ltgt+Lbeam-Lpim1).M()              , (cLtgt+cLbeam-cLpim1).M()};
  double pim1_mmom[2]  = { (Ltgt+Lbeam-Lpim1).Vect().Mag()     , (cLtgt+cLbeam-cLpim1).Vect().Mag()};
  double pim1_mcost[2] = { (Ltgt+Lbeam-Lpim1).Vect().CosTheta(), (cLtgt+cLbeam-cLpim1).Vect().CosTheta()};

  double pim2_mass[2]  = { Lpim2.M()                           , cLpim2.M()};
  double pim2_mom[2]   = { Lpim2.Vect().Mag()                  , cLpim2.Vect().Mag()};
  double pim2_cost[2]  = { Lpim2.Vect().CosTheta()             , cLpim2.Vect().CosTheta()};
  double pim2_phi[2]   = { Lpim2.Vect().Phi()                  , cLpim2.Vect().Phi()};
  double pim2_mmass[2] = { (Ltgt+Lbeam-Lpim2).M()              , (cLtgt+cLbeam-cLpim2).M()};
  double pim2_mmom[2]  = { (Ltgt+Lbeam-Lpim2).Vect().Mag()     , (cLtgt+cLbeam-cLpim2).Vect().Mag()};
  double pim2_mcost[2] = { (Ltgt+Lbeam-Lpim2).Vect().CosTheta(), (cLtgt+cLbeam-cLpim2).Vect().CosTheta()};

  double pip1_mass[2]  = { Lpip1.M()                           , cLpip1.M()};
  double pip1_mom[2]   = { Lpip1.Vect().Mag()                  , cLpip1.Vect().Mag()};
  double pip1_cost[2]  = { Lpip1.Vect().CosTheta()             , cLpip1.Vect().CosTheta()};
  double pip1_phi[2]   = { Lpip1.Vect().Phi()                  , cLpip1.Vect().Phi()};
  double pip1_mmass[2] = { (Ltgt+Lbeam-Lpip1).M()              , (cLtgt+cLbeam-cLpip1).M()};
  double pip1_mmom[2]  = { (Ltgt+Lbeam-Lpip1).Vect().Mag()     , (cLtgt+cLbeam-cLpip1).Vect().Mag()};
  double pip1_mcost[2] = { (Ltgt+Lbeam-Lpip1).Vect().CosTheta(), (cLtgt+cLbeam-cLpip1).Vect().CosTheta()};

  double pip2_mass[2]  = { Lpip2.M()                           , cLpip2.M()};
  double pip2_mom[2]   = { Lpip2.Vect().Mag()                  , cLpip2.Vect().Mag()};
  double pip2_cost[2]  = { Lpip2.Vect().CosTheta()             , cLpip2.Vect().CosTheta()};
  double pip2_phi[2]   = { Lpip2.Vect().Phi()                  , cLpip2.Vect().Phi()};
  double pip2_mmass[2] = { (Ltgt+Lbeam-Lpip2).M()              , (cLtgt+cLbeam-cLpip2).M()};
  double pip2_mmom[2]  = { (Ltgt+Lbeam-Lpip2).Vect().Mag()     , (cLtgt+cLbeam-cLpip2).Vect().Mag()};
  double pip2_mcost[2] = { (Ltgt+Lbeam-Lpip2).Vect().CosTheta(), (cLtgt+cLbeam-cLpip2).Vect().CosTheta()};

  double lam_mass[2]  = { Llam.M()                           , cLlam.M()};
  double lam_mom[2]   = { Llam.Vect().Mag()                  , cLlam.Vect().Mag()};
  double lam_cost[2]  = { Llam.Vect().CosTheta()             , cLlam.Vect().CosTheta()};
  double lam_phi[2]   = { Llam.Vect().Phi()                  , cLlam.Vect().Phi()};
  double lam_mmass[2] = { (Ltgt+Lbeam-Llam).M()              , (cLtgt+cLbeam-cLlam).M()};
  double lam_mmom[2]  = { (Ltgt+Lbeam-Llam).Vect().Mag()     , (cLtgt+cLbeam-cLlam).Vect().Mag()};
  double lam_mcost[2] = { (Ltgt+Lbeam-Llam).Vect().CosTheta(), (cLtgt+cLbeam-cLlam).Vect().CosTheta()};

  double nlam_mass[2]  = { Lnlam.M()                           , cLnlam.M()};
  double nlam_mom[2]   = { Lnlam.Vect().Mag()                  , cLnlam.Vect().Mag()};
  double nlam_cost[2]  = { Lnlam.Vect().CosTheta()             , cLnlam.Vect().CosTheta()};
  double nlam_phi[2]   = { Lnlam.Vect().Phi()                  , cLnlam.Vect().Phi()};
  double nlam_mmass[2] = { (Ltgt+Lbeam-Lnlam).M()              , (cLtgt+cLbeam-cLnlam).M()};
  double nlam_mmom[2]  = { (Ltgt+Lbeam-Lnlam).Vect().Mag()     , (cLtgt+cLbeam-cLnlam).Vect().Mag()};
  double nlam_mcost[2] = { (Ltgt+Lbeam-Lnlam).Vect().CosTheta(), (cLtgt+cLbeam-cLnlam).Vect().CosTheta()};

  double ppipi_mass[2]  = { Lppipi.M()                           , cLppipi.M()};
  double ppipi_mom[2]   = { Lppipi.Vect().Mag()                  , cLppipi.Vect().Mag()};
  double ppipi_cost[2]  = { Lppipi.Vect().CosTheta()             , cLppipi.Vect().CosTheta()};
  double ppipi_phi[2]   = { Lppipi.Vect().Phi()                  , cLppipi.Vect().Phi()};
  double ppipi_mmass[2] = { (Ltgt+Lbeam-Lppipi).M()              , (cLtgt+cLbeam-cLppipi).M()};
  double ppipi_mmom[2]  = { (Ltgt+Lbeam-Lppipi).Vect().Mag()     , (cLtgt+cLbeam-cLppipi).Vect().Mag()};
  double ppipi_mcost[2] = { (Ltgt+Lbeam-Lppipi).Vect().CosTheta(), (cLtgt+cLbeam-cLppipi).Vect().CosTheta()};

  /* flag */
  bool mnflag = false; /* mm flag */
  bool sblmnflag = false; /* lower sideband of mm flag */
  bool sbumnflag = false; /* upper sideband of mm flag */
  {
    double pro_mmass = ppipi_mmass[1];
    if(mnll<pro_mmass&&pro_mmass<mnul) mnflag = true;
    if(sblmnll<pro_mmass&&pro_mmass<sblmnul) sblmnflag = true;
    if(sbumnll<pro_mmass&&pro_mmass<sbumnul) sbumnflag = true;
  }

  //if(!mnflag){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* Missing Neutron */

  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  // ######################### //
  // Invariant/Missing Mass    //
  // ######################### //
  // pi-
  for(int f=1; f<2; f++){
    FillHist(Form("pim_Momentum_%s",frame[f].Data()),pim_mom[f]);
    FillHist(Form("pim_CosT_%s",frame[f].Data()),pim_cost[f]);
    FillHist(Form("pim_Phi_%s",frame[f].Data()),pim_phi[f]);
    FillHist(Form("pim_MMass_%s",frame[f].Data()),pim_mmass[f]);
    FillHist(Form("pim_MMomentum_%s",frame[f].Data()),pim_mmom[f]);
    FillHist(Form("pim_MCosT_%s",frame[f].Data()),pim_mcost[f]);
  }
  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Momentum_%s",frame[f].Data()),p_mom[f]);
    FillHist(Form("p_CosT_%s",frame[f].Data()),p_cost[f]);
    FillHist(Form("p_Phi_%s",frame[f].Data()),p_phi[f]);
    FillHist(Form("p_MMass_%s",frame[f].Data()),p_mmass[f]);
    FillHist(Form("p_MMomentum_%s",frame[f].Data()),p_mmom[f]);
    FillHist(Form("p_MCosT_%s",frame[f].Data()),p_mcost[f]);
  }
  // pim1
  for(int f=1; f<2; f++){
    FillHist(Form("pim1_Momentum_%s",frame[f].Data()),pim1_mom[f]);
    FillHist(Form("pim1_CosT_%s",frame[f].Data()),pim1_cost[f]);
    FillHist(Form("pim1_Phi_%s",frame[f].Data()),pim1_phi[f]);
    FillHist(Form("pim1_MMass_%s",frame[f].Data()),pim1_mmass[f]);
    FillHist(Form("pim1_MMomentum_%s",frame[f].Data()),pim1_mmom[f]);
    FillHist(Form("pim1_MCosT_%s",frame[f].Data()),pim1_mcost[f]);
  }
  // pim2
  for(int f=1; f<2; f++){
    FillHist(Form("pim2_Momentum_%s",frame[f].Data()),pim2_mom[f]);
    FillHist(Form("pim2_CosT_%s",frame[f].Data()),pim2_cost[f]);
    FillHist(Form("pim2_Phi_%s",frame[f].Data()),pim2_phi[f]);
    FillHist(Form("pim2_MMass_%s",frame[f].Data()),pim2_mmass[f]);
    FillHist(Form("pim2_MMomentum_%s",frame[f].Data()),pim2_mmom[f]);
    FillHist(Form("pim2_MCosT_%s",frame[f].Data()),pim2_mcost[f]);
  }
  // pi-1 proton
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_Mass_%s",frame[f].Data()),pip1_mass[f]);
    FillHist(Form("pip1_Momentum_%s",frame[f].Data()),pip1_mom[f]);
    FillHist(Form("pip1_CosT_%s",frame[f].Data()),pip1_cost[f]);
    FillHist(Form("pip1_Phi_%s",frame[f].Data()),pip1_phi[f]);
    FillHist(Form("pip1_MMass_%s",frame[f].Data()),pip1_mmass[f]);
    FillHist(Form("pip1_MMomentum_%s",frame[f].Data()),pip1_mmom[f]);
    FillHist(Form("pip1_MCosT_%s",frame[f].Data()),pip1_mcost[f]);
  }
  // pi-2 proton
  for(int f=1; f<2; f++){
    FillHist(Form("pip2_Mass_%s",frame[f].Data()),pip2_mass[f]);
    FillHist(Form("pip2_Momentum_%s",frame[f].Data()),pip2_mom[f]);
    FillHist(Form("pip2_CosT_%s",frame[f].Data()),pip2_cost[f]);
    FillHist(Form("pip2_Phi_%s",frame[f].Data()),pip2_phi[f]);
    FillHist(Form("pip2_MMass_%s",frame[f].Data()),pip2_mmass[f]);
    FillHist(Form("pip2_MMomentum_%s",frame[f].Data()),pip2_mmom[f]);
    FillHist(Form("pip2_MCosT_%s",frame[f].Data()),pip2_mcost[f]);
  }
  // lambda
  for(int f=1; f<2; f++){
    FillHist(Form("lam_Mass_%s",frame[f].Data()),lam_mass[f]);
    FillHist(Form("lam_Momentum_%s",frame[f].Data()),lam_mom[f]);
    FillHist(Form("lam_CosT_%s",frame[f].Data()),lam_cost[f]);
    FillHist(Form("lam_Phi_%s",frame[f].Data()),lam_phi[f]);
    FillHist(Form("lam_MMass_%s",frame[f].Data()),lam_mmass[f]);
    FillHist(Form("lam_MMomentum_%s",frame[f].Data()),lam_mmom[f]);
    FillHist(Form("lam_MCosT_%s",frame[f].Data()),lam_mcost[f]);
  }
  // Not lambda
  for(int f=1; f<2; f++){
    FillHist(Form("nlam_Mass_%s",frame[f].Data()),nlam_mass[f]);
    FillHist(Form("nlam_Momentum_%s",frame[f].Data()),nlam_mom[f]);
    FillHist(Form("nlam_CosT_%s",frame[f].Data()),nlam_cost[f]);
    FillHist(Form("nlam_Phi_%s",frame[f].Data()),nlam_phi[f]);
    FillHist(Form("nlam_MMass_%s",frame[f].Data()),nlam_mmass[f]);
    FillHist(Form("nlam_MMomentum_%s",frame[f].Data()),nlam_mmom[f]);
    FillHist(Form("nlam_MCosT_%s",frame[f].Data()),nlam_mcost[f]);
  }
  // pi-1 pi-2 proton
  for(int f=1; f<2; f++){
    FillHist(Form("ppipi_Mass_%s",frame[f].Data()),ppipi_mass[f]);
    FillHist(Form("ppipi_Momentum_%s",frame[f].Data()),ppipi_mom[f]);
    FillHist(Form("ppipi_CosT_%s",frame[f].Data()),ppipi_cost[f]);
    FillHist(Form("ppipi_Phi_%s",frame[f].Data()),ppipi_phi[f]);
    FillHist(Form("ppipi_MMass_%s",frame[f].Data()),ppipi_mmass[f]);
    FillHist(Form("ppipi_MMomentum_%s",frame[f].Data()),ppipi_mmom[f]);
    FillHist(Form("ppipi_MCosT_%s",frame[f].Data()),ppipi_mcost[f]);
    FillHist(Form("ppipi_MTheta_%s",frame[f].Data()),TMath::ACos(ppipi_mcost[f])*180.0/TMath::Pi());
  }

  /*-----------------------------------------*/
  /* 2D Plots                                */ 
  /*-----------------------------------------*/
  // IM(pip1) vs. IM(pip2)
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_IMvspip2_IM_%s",frame[f].Data()),pip1_mass[f],pip2_mass[f]);
  }
  // IM(pip1) vs. IM(pip2)
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_IM2vspip2_IM2_%s",frame[f].Data()),pip1_mass[f]*pip1_mass[f],pip2_mass[f]*pip2_mass[f]);
  }
  // IM(ppipi) vs. MM(ppipi)
  for(int f=1; f<2; f++){
    FillHist(Form("ppipi_IMvsdKppipi_MM_%s",frame[f].Data()),ppipi_mass[f],ppipi_mmass[f]);
  }
  // IM(ppipi) vs. MCosT(ppipi)
  for(int f=1; f<2; f++){
    FillHist(Form("ppipi_IMvsdKppipi_MCosT_%s",frame[f].Data()),ppipi_mass[f],ppipi_mcost[f]);
  }
  // IM(ppipi) vs. MTheta(ppipi)
  for(int f=1; f<2; f++){
    FillHist(Form("ppipi_IMvsdKppipi_MTheta_%s",frame[f].Data()),ppipi_mass[f],TMath::ACos(ppipi_mcost[f])*180.0/TMath::Pi());
  }

  return true;

}

bool MyAnalysisdKppipi::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisdKppipi::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisdKppipi::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisdKppipi::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisdKppipi::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4976;
  double sigma = 0.0069;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Missing neutron mass cut */ 
  mean  = 0.9400;
  sigma = 0.0450;
  mnll = 0.85; mnul = 1.03;
  //mnll = mean-2*sigma; mnul = mean+2*sigma;
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
  /* Lambda mass cut */ 
  mean  = 1.1155;
  sigma = 0.0020;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-5*sigma; sbllamul = mean-3*sigma;
  sbulamll = mean+3*sigma; sbulamul = mean+5*sigma;
}

bool MyAnalysisdKppipi::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisdKppipi::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anadKppipi");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  TString pname[8]  = {"pip","p","d","t","he","pim","k","o"};
  TString pname2[8] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o"};
  TString proname[6]  = {"pipi","k0","sp","sm","s","l"};
  TString proname2[6] = {"#pi^{+}#pi^{-}","K^{0}","#Sigma^{+}","#Sigma^{-}","#Sigma^{#pm}","#Lambda"};

  new TH1F( Form("CDH_Multiplicity"), Form("Multiplicity CDH;Multiplicity;Counts"), 36+1, -0.5, 36+0.5 );
  new TH1F( Form("BVC_Multiplicity"), Form("Multiplicity BVC;Multiplicity;Counts"), 8+1, -0.5, 8+0.5 );
  new TH1F( Form("CVC_Multiplicity"), Form("Multiplicity CVC;Multiplicity;Counts"), 34+1, -0.5, 34+0.5 );
  new TH1F( Form("PC_Multiplicity"), Form("Multiplicity PC;Multiplicity;Counts"), 27+1, -0.5, 27+0.5 );
  new TH1F( Form("IH_Multiplicity"), Form("Multiplicity IH;Multiplicity;Counts"), 24+1, -0.5, 24+0.5 );

  new TH1F( Form("CDC_Chi2"), Form("#chi^{2} CDC;#chi^{2};Counts"), 5000, 0.0, 100.0 );
  new TH1F( "EventNumber", "Number of Events", 20, 0, 20);
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

  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
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
  // PDf
  new TH1F( "lam_PDF", "PDF of #Lambda;PDF;Coutns", 400, 0, 40);
  new TH1F( "nlam_PDF", "PDF of Not-#Lambda;PDF;Coutns", 400, 0, 40);
  new TH1F( "pip1_PDF", "PDF of #pi^{-}p_{1};PDF;Coutns", 400, 0, 40);
  new TH1F( "pip2_PDF", "PDF of #pi^{-}p_{2};PDF;Coutns", 400, 0, 40);
  new TH2F( "pip1_PDFvspip2_PDF", "PDF of #pi^{-}p_{1} vs #pi^{-}p_{2};PDF;PDF", 400, 0, 40, 400, 0, 40);
  new TH1F( "pip1_PDF_SemiSelected", "PDF of #pi^{-}p_{1};PDF;Coutns", 400, 0, 40);
  new TH1F( "pip2_PDF_SemiSelected", "PDF of #pi^{-}p_{2};PDF;Coutns", 400, 0, 40);
  new TH2F( "pip1_PDFvspip2_PDF_SemiSelected", "PDF of #pi^{-}p_{1} vs #pi^{-}p_{2};PDF;PDF", 400, 0, 40, 400, 0, 40);
  new TH1F( "lam_PDF_Selected", "PDF of #Lambda;PDF;Coutns", 400, 0, 40);
  new TH1F( "nlam_PDF_Selected", "PDF of Not-#Lambda;PDF;Coutns", 400, 0, 40);
  new TH1F( "pip1_PDF_Selected", "PDF of #pi^{-}p_{1};PDF;Coutns", 400, 0, 40);
  new TH1F( "pip2_PDF_Selected", "PDF of #pi^{-}p_{2};PDF;Coutns", 400, 0, 40);
  new TH2F( "pip1_PDFvspip2_PDF_Selected", "PDF of #pi^{-}p_{1} vs #pi^{-}p_{2};PDF;PDF", 400, 0, 40, 400, 0, 40);
  new TH2F("lam_PDFvslam_Mass","PDF vs Invariant mass of #Lambda;PDF;IM(#Lambda) (GeV/c^{2})", 400,0,40,1000, 1.0, 2.0);
  new TH2F("lam_PDFvsnlam_Mass","PDF vs Invariant mass of Not #Lambda;PDF;IM(#Lambda) (GeV/c^{2})", 400,0,40,1000, 1.0, 2.0);
  new TH2F("lam_Massvsnlam_Mass","Invariant mass of #Lambda vs. Not #Lambda;IM(#Lambda) (GeV/c^{2});IM(Not #Lambda) (GeV/c^{2})", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  // DCA
  std::cout << "Define Histograms for DCA" << std::endl;
  new TH1F("lam_Mass","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 1000, 1.0, 2.0);
  new TH1F("nlam_Mass","Invariant mass of Not #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 1000, 1.0, 2.0);
  new TH1F( "pim_VBDIS", "VBDIS of #pi^{-};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "p1_VBDIS", "VBDIS of p_{1};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "p2_VBDIS", "VBDIS of p_{2};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "p_VBDIS", "VBDIS of p;VBDIS (cm);Coutns", 1000, 0, 10);
  new TH2F( "p1_VBDISvsp2_VBDIS", "VBDIS of p_{1} vs p_{2};VBDIS (cm);VBDIS (cm)", 100, 0, 10, 100, 0, 10);
  new TH1F( "pip1_VDIS", "VDIS of #pi^{-}p_{1};VDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "pip2_VDIS", "VDIS of #pi^{-}p_{2};VDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "lam_VDIS", "VDIS of #Lambda;VDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "nlam_VDIS", "VDIS of Not #Lambda;VDIS (cm);Coutns", 1000, 0, 10);
  new TH2F( "pip1_VDISvspip2_VDIS", "VDIS of #pi^{-}p_{1} vs #pi^{-}p_{2};VDIS (cm);VDIS (cm)", 100, 0, 10, 10, 0, 10);
  new TH1F( "pip1_VBDIS", "VBDIS of #pi^{-}p_{1};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "pip2_VBDIS", "VBDIS of #pi^{-}p_{2};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "lam_VBDIS", "VBDIS of #Lambda;VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "nlam_VBDIS", "VBDIS of Not #Lambda;VBDIS (cm);Coutns", 1000, 0, 10);
  new TH2F( "pip1_VBDISvspip2_VBDIS", "VBDIS of #pi^{-}p_{1} vs #pi^{-}p_{2};VBDIS (cm);VBDIS (cm)", 100, 0, 10, 100, 0, 10);
  new TH1F( "pip1_PBDCA", "PBDCA of #pi^{-}p_{1};PBDCA (cm);Coutns", 1000, 0, 10);
  new TH1F( "pip2_PBDCA", "PBDCA of #pi^{-}p_{2};PBDCA (cm);Coutns", 1000, 0, 10);
  new TH1F( "lam_PBDCA", "PBDCA of #Lambda;PBDCA (cm);Coutns", 1000, 0, 10);
  new TH1F( "nlam_PBDCA", "PBDCA of Not #Lambda;PBDCA (cm);Coutns", 1000, 0, 10);
  new TH2F( "pip1_PBDCAvspip2_PBDCA", "PBDCA of #pi^{-}p_{1} vs #pi^{-}p_{2};PBDCA (cm);PBDCA (cm)", 100, 0, 10, 100, 0, 10);
  new TH2F( "pim_Vertex_XY", "Vertex XY plane #pi^{-};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "pim_Vertex_ZX", "Vertex ZX plane #pi^{-};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pim_Vertex_ZY", "Vertex ZY plane #pi^{-};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p1_Vertex_XY", "Vertex XY plane p_{1};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "p1_Vertex_ZX", "Vertex ZX plane p_{1};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p1_Vertex_ZY", "Vertex ZY plane p_{1};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p2_Vertex_XY", "Vertex XY plane p_{2};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "p2_Vertex_ZX", "Vertex ZX plane p_{2};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p2_Vertex_ZY", "Vertex ZY plane p_{2};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip1_Vertex_XY", "Vertex XY plane #pi^{-}p_{1};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "pip1_Vertex_ZX", "Vertex ZX plane #pi^{-}p_{1};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip1_Vertex_ZY", "Vertex ZY plane #pi^{-}p_{1};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip2_Vertex_XY", "Vertex XY plane #pi^{-}p_{2};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "pip2_Vertex_ZX", "Vertex ZX plane #pi^{-}p_{2};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip2_Vertex_ZY", "Vertex ZY plane #pi^{-}p_{2};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "lam_Vertex_XY", "Vertex XY plane #Lambda;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "lam_Vertex_ZX", "Vertex ZX plane #Lambda;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "lam_Vertex_ZY", "Vertex ZY plane #Lambda;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "lam_PVertex_XY", "Vertex XY plane #Lambda;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "lam_PVertex_ZX", "Vertex ZX plane #Lambda;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "lam_PVertex_ZY", "Vertex ZY plane #Lambda;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "nlam_Vertex_XY", "Vertex XY plane Not #Lambda;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "nlam_Vertex_ZX", "Vertex ZX plane Not #Lambda;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "nlam_Vertex_ZY", "Vertex ZY plane Not #Lambda;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "nlam_PVertex_XY", "Vertex XY plane Not #Lambda;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "nlam_PVertex_ZX", "Vertex ZX plane Not #Lambda;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "nlam_PVertex_ZY", "Vertex ZY plane Not #Lambda;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);

  //// FWD neutral Particles
  //std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  //new TH2F( "FWD_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  //new TH2F( "FWD_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  //new TH2F( "FWD_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  //new TH2F( "FWD_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  //new TH2F( "FWD_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);


  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  std::cout << "Define Histograms for IM/MM" << std::endl;
  /* 1D Plots */
  // pi-
  new TH1F("pim_Momentum_CM","Momentum of #pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_CosT_CM","cos#theta of #pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim_Phi_CM","#phi of #pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim_MMass_CM", "^d(K^{-},#pi^{-})X missing mass;MM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_MMomentum_CM", "^d(K^{-},#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_MCosT_CM", "^d(K^{-},#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton
  new TH1F("p_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM", "^d(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM", "^d(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM", "^d(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pim1
  new TH1F("pim1_Momentum_CM","Momentum of #pi^{-}_{1};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim1_CosT_CM","cos#theta of #pi^{-}_{1};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim1_Phi_CM","#phi of #pi^{-}_{1};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim1_MMass_CM", "d(K^{-},#pi^{-}_{1})X missing mass;MM(#pi^{-}_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim1_MMomentum_CM", "d(K^{-},#pi^{-}_{1})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim1_MCosT_CM", "d(K^{-},#pi^{-}_{1})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pim2
  new TH1F("pim2_Momentum_CM","Momentum of #pi^{-}_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim2_CosT_CM","cos#theta of #pi^{-}_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim2_Phi_CM","#phi of #pi^{-}_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim2_MMass_CM", "d(K^{-},#pi^{-}_{2})X missing mass;MM(#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim2_MMomentum_CM", "d(K^{-},#pi^{-}_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim2_MCosT_CM", "d(K^{-},#pi^{-}_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-1,proton
  new TH1F("pip1_Mass_CM","Invariant mass of #pi^{-}_{1}p;IM(#pi^{-}_{1}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip1_Momentum_CM","Momentum of #pi^{-}_{1}p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip1_CosT_CM","cos#theta of #pi^{-}_{1}p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pip1_Phi_CM","#phi of #pi^{-}_{1}p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pip1_MMass_CM", "d(K^{-},#pi^{-}_{1}p)X missing mass;MM(#pi^{-}_{1}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip1_MMomentum_CM", "d(K^{-},#pi^{-}_{1}p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip1_MCosT_CM", "d(K^{-},#pi^{-}_{1}p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-,proton2
  new TH1F("pip2_Mass_CM","Invariant mass of #pi^{-}_{2}p;IM(#pi^{-}_{2}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip2_Momentum_CM","Momentum of #pi^{-}_{2}p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip2_CosT_CM","cos#theta of #pi^{-}_{2}p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pip2_Phi_CM","#phi of #pi^{-}_{2}p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pip2_MMass_CM", "d(K^{-},#pi^{-}_{2}p)X missing mass;MM(#pi^{-}_{2}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip2_MMomentum_CM", "d(K^{-},#pi^{-}_{2}p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip2_MCosT_CM", "d(K^{-},#pi^{-}_{2}p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Lambda
  new TH1F("lam_Mass_CM","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_Momentum_CM","Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_CosT_CM","cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_Phi_CM","#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("lam_MMass_CM", "d(K^{-},#Lambda)X missing mass;MM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_MMomentum_CM", "d(K^{-},#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_MCosT_CM", "d(K^{-},#Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Not Lambda
  new TH1F("nlam_Mass_CM","Invariant mass of Not #Lambda;IM(Not #Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nlam_Momentum_CM","Momentum of Not #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nlam_CosT_CM","cos#theta of Not #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("nlam_Phi_CM","#phi of Not #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("nlam_MMass_CM", "d(K^{-},Not #Lambda)X missing mass;MM(Not #Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nlam_MMomentum_CM", "d(K^{-},Not #Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nlam_MCosT_CM", "d(K^{-},Not #Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-1,pi-2,proton
  new TH1F("ppipi_Mass_CM","Invariant mass of p#pi^{-}#pi^{-};IM(p#pi^{-}#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppipi_Momentum_CM","Momentum of p#pi^{-}#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppipi_CosT_CM","cos#theta of p#pi^{-}#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("ppipi_Phi_CM","#phi of p#pi^{-}#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("ppipi_MMass_CM", "d(K^{-},p#pi^{-}#pi^{-})X missing mass;MM(p#pi^{-}#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppipi_MMomentum_CM", "d(K^{-},p#pi^{-}#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppipi_MCosT_CM", "d(K^{-},p#pi^{-}#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("ppipi_MTheta_CM", "d(K^{-},p#pi^{-}#pi^{-})X #theta;#theta;Counts", 200, -10.0, 190.0);
  new TH1F("ppipi_MomentumTransfer_CM","Momentum transfer to p#pi^{-}#pi^{-};Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  /* 2D Plots */
  new TH2F("pip1_IMvspip2_IM_All", "(#pi^{-}_{1}p) invariant mass vs. (#pi^{-}_{2}p) invariant mass;IM(#pi^{-}_{1}p) (GeV/c^{2});IM(#pi^{-}_{2}p) (GeV/c^{2})", 2000, 1.0, 3.0, 2000, 1.0, 3.0);
  new TH2F("pip1_IMvspip2_IM_CM", "(#pi^{-}_{1}p) invariant mass vs. (#pi^{-}_{2}p) invariant mass;IM(#pi^{-}_{1}p) (GeV/c^{2});IM(#pi^{-}_{2}p) (GeV/c^{2})", 2000, 1.0, 3.0, 2000, 1.0, 3.0);
  new TH2F("pip1_IM2vspip2_IM2_All", "(#pi^{-}_{1}p) invariant mass vs. (#pi^{-}_{2}p) invariant mass;IM(#pi^{-}_{1}p) (GeV/c^{2});IM(#pi^{-}_{2}p) (GeV/c^{2})", 2000, 1.0, 5.0, 2000, 1.0, 5.0);
  new TH2F("pip1_IM2vspip2_IM2_CM", "(#pi^{-}_{1}p) invariant mass vs. (#pi^{-}_{2}p) invariant mass;IM(#pi^{-}_{1}p) (GeV/c^{2});IM(#pi^{-}_{2}p) (GeV/c^{2})", 2000, 1.0, 5.0, 2000, 1.0, 5.0);
  new TH2F("ppipi_IMvsdKppipi_MM_CM", "(p#pi^{-}#pi^{-}) invariant mass vs. d(K^{-},p#pi^{-}#pi^{-})X missing mass;IM(p#pi^{-}#pi^{-}) (GeV/c^{2});MM(p#pi^{-}#pi^{-}) (GeV/c^{2})", 200, 1.0, 3.0, 200, 0.0, 2.0);
  new TH2F("ppipi_IMvsdKppipi_MCosT_CM", "(p#pi^{-}#pi^{-}) invariant mass vs. cos#theta_{n};IM(p#pi^{-}#pi^{-}) (GeV/c^{2});cos#theta_{n}", 200, 1.0, 3.0, 300, -1.5, 1.5);
  new TH2F("ppipi_IMvsdKppipi_MTheta_CM", "(p#pi^{-}#pi^{-}) invariant mass vs. #theta_{n};IM(p#pi^{-}#pi^{-}) (GeV/c^{2});#theta_{n}", 200, 1.0, 3.0, 200, -10.5, 190.0);
  new TH2F("ppipi_IMvsMomentumTransfer_CM", "(p#pi^{-}#pi^{-}) invariant mass vs. Momentum transfer;IM(p#pi^{-}#pi^{-}) (GeV/c^{2});Momentum Transfer (GeV/c)", 200, 1.0, 3.0, 200, 0.0, 2.0);
  new TH2F("ppipi_IMvsOA_CM", "(p#pi^{-}pi^{-}#pi^{-}) invariant mass vs. Opening Angle;IM(p#pi^{-}#pi^{-}) (GeV/c^{2});cos#theta_{#Lambdap} ", 200, 1.0, 3.0, 300, -1.5, 1.5);

  std::cout << "== Finish Histogram Initialization ==" << std::endl;
  return true;

}
