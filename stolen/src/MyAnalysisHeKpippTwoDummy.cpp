// MyAnalysisHeKpippTwoDummy.cpp

#include "MyAnalysisHeKpippTwoDummy.h"

MyAnalysisHeKpippTwoDummy::MyAnalysisHeKpippTwoDummy(TFile* rt, ConfMan* conf)
{
  eventdiff = 0;
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKpippTwoDummy::~MyAnalysisHeKpippTwoDummy()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKpippTwoDummy::Clear()
{
}

bool MyAnalysisHeKpippTwoDummy::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  int ievent=0;
  rtFile->cd();

  FillHist("EventNumber",ievent); ievent++; /* All Events */

  if(particle->nBeam()!=1) return false;
  FillHist("EventNumber",ievent); ievent++; /* Beam Analysis */
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  //if(!header->IsTrig(Trig_KCDH3)){ return true; }
  //FillHist("EventNumber",ievent); ievent++; /* KCDH3 Trigger */

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
  int MulBVC=0;
  for(int i=0; i<blMan->nBVC(); i++){
    HodoscopeLikeHit* hit = blMan->BVC(i);
    if(hit->CheckRange()) MulBVC++;
  }
  FillHist("BVC_Multiplicity",MulBVC);
  int MulCVC=0;
  for(int i=0; i<blMan->nCVC(); i++){
    HodoscopeLikeHit* hit = blMan->CVC(i);
    if(hit->CheckRange()) MulCVC++;
  }
  FillHist("CVC_Multiplicity",MulCVC);
  int MulPC=0;
  for(int i=0; i<blMan->nPC(); i++){
    HodoscopeLikeHit* hit = blMan->PC(i);
    if(hit->CheckRange()) MulPC++;
  }
  FillHist("PC_Multiplicity",MulPC);
  if(MulBVC!=0||MulCVC!=0||MulPC!=0) return true;
  FillHist("EventNumber",ievent); ievent++; /* No FWD charged hit */

  int MulIH=0;
  for(int i=0; i<cdsMan->nIH(); i++){
    HodoscopeLikeHit* hit = cdsMan->IH(i);
    if(hit->CheckRange()) MulIH++;
  }
  FillHist("IH_Multiplicity",MulIH);
  //if(MulIH<=3){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* 3 hits in IH */

  if(particle->nCDS()!=3) return true;
  FillHist("EventNumber",ievent); ievent++; /* 3 tracks in CDS */

  if(particle->nProton()!=2) return true;
  if(particle->nPiminus()!=1) return true;
  FillHist("EventNumber",ievent); ievent++; /* 2 protons and pion in CDS */

  // ############## //
  // p1/p2 decision //
  // ############## //
  TVector3 vertex;
  double vdis = 9999;
  bool FIDUCIAL = false;
  int fcds = -1;
  for(int it=0; it<particle->nProton(); it++){
    pCDS* cds = particle->proton(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      fcds=it;
    }
  }
  if(fcds<0){ return true; }

  pCDS* pim = particle->pim(0);
  pCDS* p1 = particle->proton(fcds);
  pCDS* p2 = particle->proton(1-fcds);

  // ############## //
  // CDC Chi-square //
  // ############## //
  if(pim->chi()>30) return true;
  if(p1->chi()>30) return true;
  if(p2->chi()>30) return true;
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
        if(product->daughter1()==p1->id()||product->daughter2()==p1->id()){
          pip1 = product;
        }
        else if(product->daughter1()==p2->id()||product->daughter2()==p2->id()){
          pip2 = product;
        }
      }
    }
  }
  if(pip1==0||pip2==0){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* 2-track Combination */

  FillHist(Form("pip1_IMvspip2_IM_All"),pip1->mass(),pip2->mass());
  FillHist(Form("pip1_IM2vspip2_IM2_All"),pip1->mass()*pip1->mass(),pip2->mass()*pip2->mass());

  // Target
  TVector3 tmpZeroV;
  TLorentzVector tmpLtgt; tmpLtgt.SetVectM(tmpZeroV,ThreeHeMass);
  TLorentzVector tmpLbeam = beam->GetLorentzVector(vertex);

  TLorentzVector tmpLpim = pim->GetLorentzVector();
  TLorentzVector tmpLp1 = p1->GetLorentzVector();
  TLorentzVector tmpLp2 = p2->GetLorentzVector();
  double tmp_im = (tmpLpim+tmpLp1+tmpLp2).M();
  double tmp_mm = (tmpLtgt+tmpLbeam-tmpLpim-tmpLp1-tmpLp2).M();
  FillHist(Form("pipp_Mass_Rough"),tmp_im);
  FillHist(Form("pipp_MMass_Rough"),tmp_mm);

  bool pip1flag = false, pip2flag = false;
  pCDS* lam = 0;
  pCDS* nlam = 0;
  pCDS* proton = 0;
  pCDS* nproton = 0;
  /* For pip1 */
  double tmpper1[5] = {pip1->mass(),pip1->vdis(),pip1->pbdca(),p2->vbdis(),0.0};
  double tmppdf1[5] = {0,0,0,0,0};
  /* Input : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  /*Output : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  TrackTools::PDFLambda(tmpper1,tmppdf1);
  tmppdf1[4] = 1.0;
  double tmpper2[5] = {pip2->mass(),pip2->vdis(),pip2->pbdca(),p1->vbdis(),0.0};
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
    proton = p2;
    nproton = p1;
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
    proton = p1;
    nproton = p2;
    pdf=pdf2;
    npdf=pdf1;
    }
  }
  //if(lamll<pip1->mass()&&pip1->mass()<lamul){
  //  pip1flag = true;
  //  lam = pip1;
  //  proton = p2;
  //}
  //if(lamll<pip2->mass()&&pip2->mass()<lamul){
  //  pip2flag = true;
  //  pip1flag = false;
  //  lam = pip2;
  //  proton = p1;
  //}
  if(!pip1flag&&!pip2flag){ return true; }
  if(lam==0){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* Lambda Reconstruction */
  FillHist("pip1_PDFvspip2_PDF_Selected",pdf1,pdf2);
  FillHist("lam_PDF_Selected",pdf);
  FillHist("nlam_PDF_Selected",npdf);

  // ################ //
  // Vertex decision1 //
  // ################ //
  const int fiducial = CID_Fiducial;
  //const int fiducial = CID_CellTube;
  //const int fiducial = CID_CellAlBe;
  //const int fiducial = CID_DEF;
  if(pip1flag){
    vertex = p2->vbeam();
    vdis = p2->vbdis();
    if(GeomTools::GetID(vertex)!=fiducial){ return true; }
  }
  if(pip2flag){
    vertex = p1->vbeam();
    vdis = p1->vbdis();
    if(GeomTools::GetID(vertex)!=fiducial){ return true; }
  }
  FillHist("EventNumber",ievent); ievent++; /* Vertex Decision 1 */

  // ################ //
  // Vertex decision2 //
  // ################ //
  TVector3 vertex_lam = lam->vbeam();
  if(GeomTools::GetID(vertex_lam)!=fiducial){ return true; }
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
  /* Pi Minus */
  vtxb = pim->vbeam();
  FillHist("pim_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pim_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pim_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pim_VBDIS",pim->vbdis());
  /* Proton1 */
  vtxb = p1->vbeam();
  FillHist("p1_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p1_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p1_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p1_VBDIS",p1->vbdis());
  /* Proton2 */
  vtxb = p2->vbeam();
  FillHist("p2_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p2_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p2_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p2_VBDIS",p2->vbdis());
  /* Pi+ and Proton1 */
  vtxb = pip1->vbeam();
  FillHist("pip1_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pip1_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pip1_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pip1_Mass",pip1->vdis());
  FillHist("pip1_VDIS",pip1->vdis());
  FillHist("pip1_VBDIS",pip1->vbdis());
  FillHist("pip1_PBDCA",pip1->pbdca());
  /* Pi+ and Proton2 */
  vtxb = pip2->vbeam();
  FillHist("pip2_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pip2_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pip2_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pip2_Mass",pip2->vdis());
  FillHist("pip2_VDIS",pip2->vdis());
  FillHist("pip2_VBDIS",pip2->vbdis());
  FillHist("pip2_PBDCA",pip2->pbdca());
  /* 2D Plots */
  FillHist("p1_VBDISvsp2_VBDIS",p1->vbdis(),p2->vbdis());
  FillHist("pip1_VDISvspip2_VDIS",pip1->vdis(),pip2->vdis());
  FillHist("pip1_VBDISvspip2_VBDIS",pip1->vbdis(),pip2->vbdis());
  FillHist("pip1_PBDCAvspip2_PBDCA",pip1->pbdca(),pip2->pbdca());
  /* Proton */
  vtxb = proton->vbeam();
  FillHist("p_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p_VBDIS",proton->vbdis());
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
  //if(proton->pbdca()>1.0){ return true; }
  //if(lam->pbdca()>1.0){ return true; }
  //FillHist("EventNumber",ievent); ievent++; /* DCA Decision */

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }

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
  pNC*  nc =0;
  if(fnc!=-1)  nc =particle->nc(fnc);
  //if(nc==0) return true;

  TString frame[2] = {"Lab","CM"};
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

  // Pi+, Pi-, and Pi+ Pi- pair //
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  TLorentzVector Lp1 = p1->GetLorentzVector();
  TLorentzVector cLp1 = p1->GetLorentzVector(); cLp1.Boost(-boost);
  TLorentzVector Lp2 = p2->GetLorentzVector();
  TLorentzVector cLp2 = p2->GetLorentzVector(); cLp2.Boost(-boost);
  TLorentzVector Lpip1 = pip1->GetLorentzVector();
  TLorentzVector cLpip1 = pip1->GetLorentzVector(); cLpip1.Boost(-boost);
  TLorentzVector Lpip2 = pip2->GetLorentzVector();
  TLorentzVector cLpip2 = pip2->GetLorentzVector(); cLpip2.Boost(-boost);
  TLorentzVector Llam = lam->GetLorentzVector();
  TLorentzVector cLlam = lam->GetLorentzVector();
  cLlam.Boost(-boost);
  TLorentzVector Lnlam = nlam->GetLorentzVector();
  TLorentzVector cLnlam = nlam->GetLorentzVector();
  cLnlam.Boost(-boost);
  TLorentzVector Lp = proton->GetLorentzVector();
  TLorentzVector cLp = proton->GetLorentzVector(); cLp.Boost(-boost);
  TLorentzVector Lnp = nproton->GetLorentzVector();
  TLorentzVector cLnp = nproton->GetLorentzVector(); cLnp.Boost(-boost);
  if(eventdiff<2){
    p_lambda = Llam;
    p_proton = Lp;
    eventdiff++;
    return true;
  }
  else{
    Llam = p_lambda;
    cLlam = p_lambda; cLlam.Boost(-boost);
    p_lambda = Llam;
    p_proton = Lp;
  }

  TLorentzVector Lpipp = Llam + Lp, cLpipp = cLlam + cLp;
  TLorentzVector Lnpipp = Lnlam + Lnp, cLnpipp = cLnlam + cLnp;
  TLorentzVector Lmn = Ltgt+Lbeam-Lpipp;
  TLorentzVector cLmn = cLtgt+cLbeam-cLpipp;
  TLorentzVector Lpipn = Llam + Lmn, cLpipn = cLlam + cLmn;


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

  double p1_mass[2]  = { Lp1.M()                           , cLp1.M()};
  double p1_mom[2]   = { Lp1.Vect().Mag()                  , cLp1.Vect().Mag()};
  double p1_cost[2]  = { Lp1.Vect().CosTheta()             , cLp1.Vect().CosTheta()};
  double p1_phi[2]   = { Lp1.Vect().Phi()                  , cLp1.Vect().Phi()};
  double p1_mmass[2] = { (Ltgt+Lbeam-Lp1).M()              , (cLtgt+cLbeam-cLp1).M()};
  double p1_mmom[2]  = { (Ltgt+Lbeam-Lp1).Vect().Mag()     , (cLtgt+cLbeam-cLp1).Vect().Mag()};
  double p1_mcost[2] = { (Ltgt+Lbeam-Lp1).Vect().CosTheta(), (cLtgt+cLbeam-cLp1).Vect().CosTheta()};

  double p2_mass[2]  = { Lp2.M()                           , cLp2.M()};
  double p2_mom[2]   = { Lp2.Vect().Mag()                  , cLp2.Vect().Mag()};
  double p2_cost[2]  = { Lp2.Vect().CosTheta()             , cLp2.Vect().CosTheta()};
  double p2_phi[2]   = { Lp2.Vect().Phi()                  , cLp2.Vect().Phi()};
  double p2_mmass[2] = { (Ltgt+Lbeam-Lp2).M()              , (cLtgt+cLbeam-cLp2).M()};
  double p2_mmom[2]  = { (Ltgt+Lbeam-Lp2).Vect().Mag()     , (cLtgt+cLbeam-cLp2).Vect().Mag()};
  double p2_mcost[2] = { (Ltgt+Lbeam-Lp2).Vect().CosTheta(), (cLtgt+cLbeam-cLp2).Vect().CosTheta()};

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

  double pipp_mass[2]  = { Lpipp.M()                           , cLpipp.M()};
  double pipp_mom[2]   = { Lpipp.Vect().Mag()                  , cLpipp.Vect().Mag()};
  double pipp_cost[2]  = { Lpipp.Vect().CosTheta()             , cLpipp.Vect().CosTheta()};
  double pipp_phi[2]   = { Lpipp.Vect().Phi()                  , cLpipp.Vect().Phi()};
  double pipp_mmass[2] = { (Ltgt+Lbeam-Lpipp).M()              , (cLtgt+cLbeam-cLpipp).M()};
  double pipp_mmom[2]  = { (Ltgt+Lbeam-Lpipp).Vect().Mag()     , (cLtgt+cLbeam-cLpipp).Vect().Mag()};
  double pipp_mcost[2] = { (Ltgt+Lbeam-Lpipp).Vect().CosTheta(), (cLtgt+cLbeam-cLpipp).Vect().CosTheta()};

  double npipp_mass[2]  = { Lnpipp.M()                           , cLnpipp.M()};
  double npipp_mom[2]   = { Lnpipp.Vect().Mag()                  , cLnpipp.Vect().Mag()};
  double npipp_cost[2]  = { Lnpipp.Vect().CosTheta()             , cLnpipp.Vect().CosTheta()};
  double npipp_phi[2]   = { Lnpipp.Vect().Phi()                  , cLnpipp.Vect().Phi()};
  double npipp_mmass[2] = { (Ltgt+Lbeam-Lnpipp).M()              , (cLtgt+cLbeam-cLnpipp).M()};
  double npipp_mmom[2]  = { (Ltgt+Lbeam-Lnpipp).Vect().Mag()     , (cLtgt+cLbeam-cLnpipp).Vect().Mag()};
  double npipp_mcost[2] = { (Ltgt+Lbeam-Lnpipp).Vect().CosTheta(), (cLtgt+cLbeam-cLnpipp).Vect().CosTheta()};

  double pipn_mass[2]  = { Lpipn.M()                           , cLpipn.M()};
  double pipn_mom[2]   = { Lpipn.Vect().Mag()                  , cLpipn.Vect().Mag()};
  double pipn_cost[2]  = { Lpipn.Vect().CosTheta()             , cLpipn.Vect().CosTheta()};
  double pipn_phi[2]   = { Lpipn.Vect().Phi()                  , cLpipn.Vect().Phi()};
  double pipn_mmass[2] = { (Ltgt+Lbeam-Lpipn).M()              , (cLtgt+cLbeam-cLpipn).M()};
  double pipn_mmom[2]  = { (Ltgt+Lbeam-Lpipn).Vect().Mag()     , (cLtgt+cLbeam-cLpipn).Vect().Mag()};
  double pipn_mcost[2] = { (Ltgt+Lbeam-Lpipn).Vect().CosTheta(), (cLtgt+cLbeam-cLpipn).Vect().CosTheta()};

  //double q[2] = {(Lbeam.Vect()-(Ltgt+Lbeam-Lpipp).Vect()).Mag(), (cLbeam.Vect()-(cLtgt+cLbeam-cLpipp).Vect()).Mag()};
  double q[2] = {(Lbeam.Vect()-Lmn.Vect()).Mag(), (cLbeam.Vect()-cLmn.Vect()).Mag()};
  //double tn = (cLtgt+cLbeam-cLpipp).E() - 0.939565;
  //double tn = (cLmn).E() - 0.939565;
  //double tp = cLp.E() - 0.938272;
  //double tl = (cLlam).E() - 1.115683;
  double tn = cLmn.E() - cLmn.M();
  double tp = cLp.E() - cLp.M();
  double tl = cLlam.E() - cLlam.M();
  double qvalue = tn+tp+tl;
  //double qvalue = (Lbeam+Ltgt).E() - 1.115683 - 0.938272 - 0.939565;
  double cos_Lp = (cLp.Vect()*cLlam.Vect())/(cLp.Vect().Mag())/(cLlam.Vect().Mag());

  FillHist(Form("pipp_Mass_Rough_vs_CM"),tmp_im,pipp_mass[1]);
  FillHist(Form("pipp_MMass_Rough_vs_CM"),tmp_mm,pipp_mmass[1]);

  /* flag */
  bool mnflag = false; /* mm flag */
  bool sblmnflag = false; /* lower sideband of mm flag */
  bool sbumnflag = false; /* upper sideband of mm flag */
  {
    double pro_mmass = pipp_mmass[1];
    if(mnll<pro_mmass&&pro_mmass<mnul) mnflag = true;
    if(sblmnll<pro_mmass&&pro_mmass<sblmnul) sblmnflag = true;
    if(sbumnll<pro_mmass&&pro_mmass<sbumnul) sbumnflag = true;
  }

  //if(!lamflag1&&!lamflag2){ return true; }
  //if(lamflag1){ return true; }
  //if(!lamflag2){ return true; }
  if(!mnflag){ return true; }
  //if(mnflag){ return true; }
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
  // proton1
  for(int f=1; f<2; f++){
    FillHist(Form("p1_Momentum_%s",frame[f].Data()),p1_mom[f]);
    FillHist(Form("p1_CosT_%s",frame[f].Data()),p1_cost[f]);
    FillHist(Form("p1_Phi_%s",frame[f].Data()),p1_phi[f]);
    FillHist(Form("p1_MMass_%s",frame[f].Data()),p1_mmass[f]);
    FillHist(Form("p1_MMomentum_%s",frame[f].Data()),p1_mmom[f]);
    FillHist(Form("p1_MCosT_%s",frame[f].Data()),p1_mcost[f]);
  }
  // proton2
  for(int f=1; f<2; f++){
    FillHist(Form("p2_Momentum_%s",frame[f].Data()),p2_mom[f]);
    FillHist(Form("p2_CosT_%s",frame[f].Data()),p2_cost[f]);
    FillHist(Form("p2_Phi_%s",frame[f].Data()),p2_phi[f]);
    FillHist(Form("p2_MMass_%s",frame[f].Data()),p2_mmass[f]);
    FillHist(Form("p2_MMomentum_%s",frame[f].Data()),p2_mmom[f]);
    FillHist(Form("p2_MCosT_%s",frame[f].Data()),p2_mcost[f]);
  }
  // pi- proton1
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_Mass_%s",frame[f].Data()),pip1_mass[f]);
    FillHist(Form("pip1_Momentum_%s",frame[f].Data()),pip1_mom[f]);
    FillHist(Form("pip1_CosT_%s",frame[f].Data()),pip1_cost[f]);
    FillHist(Form("pip1_Phi_%s",frame[f].Data()),pip1_phi[f]);
    FillHist(Form("pip1_MMass_%s",frame[f].Data()),pip1_mmass[f]);
    FillHist(Form("pip1_MMomentum_%s",frame[f].Data()),pip1_mmom[f]);
    FillHist(Form("pip1_MCosT_%s",frame[f].Data()),pip1_mcost[f]);
  }
  // pi- proton2
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
  // nlambda
  for(int f=1; f<2; f++){
    FillHist(Form("nlam_Mass_%s",frame[f].Data()),nlam_mass[f]);
    FillHist(Form("nlam_Momentum_%s",frame[f].Data()),nlam_mom[f]);
    FillHist(Form("nlam_CosT_%s",frame[f].Data()),nlam_cost[f]);
    FillHist(Form("nlam_Phi_%s",frame[f].Data()),nlam_phi[f]);
    FillHist(Form("nlam_MMass_%s",frame[f].Data()),nlam_mmass[f]);
    FillHist(Form("nlam_MMomentum_%s",frame[f].Data()),nlam_mmom[f]);
    FillHist(Form("nlam_MCosT_%s",frame[f].Data()),nlam_mcost[f]);
  }
  // pi- proton1 proton2
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_Mass_%s",frame[f].Data()),pipp_mass[f]);
    FillHist(Form("pipp_Momentum_%s",frame[f].Data()),pipp_mom[f]);
    FillHist(Form("pipp_CosT_%s",frame[f].Data()),pipp_cost[f]);
    FillHist(Form("pipp_Phi_%s",frame[f].Data()),pipp_phi[f]);
    FillHist(Form("pipp_MMass_%s",frame[f].Data()),pipp_mmass[f]);
    FillHist(Form("pipp_MMomentum_%s",frame[f].Data()),pipp_mmom[f]);
    FillHist(Form("pipp_MCosT_%s",frame[f].Data()),pipp_mcost[f]);
    FillHist(Form("pipp_MTheta_%s",frame[f].Data()),TMath::ACos(pipp_mcost[f])*180.0/TMath::Pi());
  }
  // Not pi- proton1 proton2
  for(int f=1; f<2; f++){
    FillHist(Form("npipp_Mass_%s",frame[f].Data()),npipp_mass[f]);
    FillHist(Form("npipp_Momentum_%s",frame[f].Data()),npipp_mom[f]);
    FillHist(Form("npipp_CosT_%s",frame[f].Data()),npipp_cost[f]);
    FillHist(Form("npipp_Phi_%s",frame[f].Data()),npipp_phi[f]);
    FillHist(Form("npipp_MMass_%s",frame[f].Data()),npipp_mmass[f]);
    FillHist(Form("npipp_MMomentum_%s",frame[f].Data()),npipp_mmom[f]);
    FillHist(Form("npipp_MCosT_%s",frame[f].Data()),npipp_mcost[f]);
    FillHist(Form("npipp_MTheta_%s",frame[f].Data()),TMath::ACos(npipp_mcost[f])*180.0/TMath::Pi());
  }
  // pi- proton m-neutron
  for(int f=1; f<2; f++){
    FillHist(Form("pipn_Mass_%s",frame[f].Data()),pipn_mass[f]);
    FillHist(Form("pipn_Momentum_%s",frame[f].Data()),pipn_mom[f]);
    FillHist(Form("pipn_CosT_%s",frame[f].Data()),pipn_cost[f]);
    FillHist(Form("pipn_Phi_%s",frame[f].Data()),pipn_phi[f]);
    FillHist(Form("pipn_MMass_%s",frame[f].Data()),pipn_mmass[f]);
    FillHist(Form("pipn_MMomentum_%s",frame[f].Data()),pipn_mmom[f]);
    FillHist(Form("pipn_MCosT_%s",frame[f].Data()),pipn_mcost[f]);
    FillHist(Form("pipn_MTheta_%s",frame[f].Data()),TMath::ACos(pipn_mcost[f])*180.0/TMath::Pi());
  }
  // Momentum transfer
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_MomentumTransfer_%s",frame[f].Data()),q[f]);
  }
  // Opening Angle between Lambda and Proton
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_OA_%s",frame[f].Data()),cos_Lp);
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
  // IM(pipp) vs. PDF
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvslam_PDF_%s",frame[f].Data()),pipp_mass[f],pdf);
  }
  // MM(pipp) vs. PDF
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpipp_MMvslam_PDF_%s",frame[f].Data()),pipp_mmass[f],pdf);
  }
  // IM(pipp) vs. MM(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MM_%s",frame[f].Data()),pipp_mass[f],pipp_mmass[f]);
  }
  // IM(npipp) vs. MM(npipp)
  for(int f=1; f<2; f++){
    FillHist(Form("npipp_IMvsHeKnpipp_MM_%s",frame[f].Data()),npipp_mass[f],npipp_mmass[f]);
  }
  // IM(pipp) vs. MCosT(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MCosT_%s",frame[f].Data()),pipp_mass[f],pipp_mcost[f]);
  }
  // IM(pipp) vs. MTheta(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MTheta_%s",frame[f].Data()),pipp_mass[f],TMath::ACos(pipp_mcost[f])*180.0/TMath::Pi());
  }
  // IM(pipp) vs. MMom(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvspipp_MMom_%s",frame[f].Data()),pipp_mass[f],pipp_mmom[f]);
  }
  // MM(pipp) vs. MMom(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_MMvspipp_MMom_%s",frame[f].Data()),pipp_mmass[f],pipp_mmom[f]);
  }
  // IM(pipp) vs. MomentumTransfer
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsMomentumTransfer_%s",frame[f].Data()),pipp_mass[f],q[f]);
  }
  // IM(pipp) vs. OA
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsOA_%s",frame[f].Data()),pipp_mass[f],cos_Lp);
  }
  // Dalitz Plot
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_DalitzPlot_%s",frame[f].Data()),(tp-tl)/sqrt(3)/qvalue,tn/qvalue);
  }
  // Dalitz Plot
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IM2vspipn_IM2_%s",frame[f].Data()),pipp_mass[f]*pipp_mass[f],pipn_mass[f]*pipn_mass[f]);
  }

  if(nc==0) return true;
  FillHist("EventNumber",ievent); ievent++; /* NC Hit */
  /* NC */
  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
  double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};
  double pnc_mass[2]    = { (Lp+Ln).M()                           , (cLp+cLn).M()};
  double pnc_mom[2]     = { (Lp+Ln).Vect().Mag()                  , (cLp+cLn).Vect().Mag()};
  double pnc_cost[2]    = { (Lp+Ln).Vect().CosTheta()             , (cLp+cLn).Vect().CosTheta()};
  double pnc_phi[2]     = { (Lp+Ln).Vect().Phi()                  , (cLp+cLn).Vect().Phi()};
  double pnc_mmass[2]   = { (Ltgt+Lbeam-Lp-Ln).M()              , (cLtgt+cLbeam-cLp-cLn).M()};
  double pnc_mmass2[2]   = { (Ltgt+Lbeam-Lp-Ln).M2()              , (cLtgt+cLbeam-cLp-cLn).M2()};
  double pnc_mmom[2]    = { (Ltgt+Lbeam-Lp-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLp-cLn).Vect().Mag()};
  double pnc_mcost[2]   = { (Ltgt+Lbeam-Lp-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLp-cLn).Vect().CosTheta()};
  double pippnc_mass[2]    = { (Lpipp+Ln).M()                           , (cLpipp+cLn).M()};
  double pippnc_mom[2]     = { (Lpipp+Ln).Vect().Mag()                  , (cLpipp+cLn).Vect().Mag()};
  double pippnc_cost[2]    = { (Lpipp+Ln).Vect().CosTheta()             , (cLpipp+cLn).Vect().CosTheta()};
  double pippnc_phi[2]     = { (Lpipp+Ln).Vect().Phi()                  , (cLpipp+cLn).Vect().Phi()};
  double pippnc_mmass[2]   = { (Ltgt+Lbeam-Lpipp-Ln).M()              , (cLtgt+cLbeam-cLpipp-cLn).M()};
  double pippnc_mmass2[2]   = { (Ltgt+Lbeam-Lpipp-Ln).M2()              , (cLtgt+cLbeam-cLpipp-cLn).M2()};
  double pippnc_mmom[2]    = { (Ltgt+Lbeam-Lpipp-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLpipp-cLn).Vect().Mag()};
  double pippnc_mcost[2]   = { (Ltgt+Lbeam-Lpipp-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLpipp-cLn).Vect().CosTheta()};

  double q_2 = (Lbeam.Vect()-Ln.Vect()).Mag();
  //tn = cLn.E() - cLn.M();
  //qvalue = tn+tp+tl;

  // n
  for(int f=1; f<2; f++){
    FillHist(Form("n_Momentum_%s",frame[f].Data()),nc_mom[f]);
    FillHist(Form("n_CosT_%s",frame[f].Data()),nc_cost[f]);
    FillHist(Form("n_Phi_%s",frame[f].Data()),nc_phi[f]);
    FillHist(Form("n_MMass_%s",frame[f].Data()),nc_mmass[f]);
    FillHist(Form("n_MMomentum_%s",frame[f].Data()),nc_mmom[f]);
    FillHist(Form("n_MCosT_%s",frame[f].Data()),nc_mcost[f]);
  }
  // n + p
  for(int f=1; f<2; f++){
    FillHist(Form("pn_Mass_%s",frame[f].Data()),pnc_mass[f]);
    FillHist(Form("pn_Momentum_%s",frame[f].Data()),pnc_mom[f]);
    FillHist(Form("pn_CosT_%s",frame[f].Data()),pnc_cost[f]);
    FillHist(Form("pn_Phi_%s",frame[f].Data()),pnc_phi[f]);
    FillHist(Form("pn_MMass_%s",frame[f].Data()),pnc_mmass[f]);
    FillHist(Form("pn_MMass2_%s",frame[f].Data()),pnc_mmass2[f]);
    FillHist(Form("pn_MMomentum_%s",frame[f].Data()),pnc_mmom[f]);
    FillHist(Form("pn_MCosT_%s",frame[f].Data()),pnc_mcost[f]);
  }
  // n + pipp
  for(int f=1; f<2; f++){
    FillHist(Form("pippn_Mass_%s",frame[f].Data()),pippnc_mass[f]);
    FillHist(Form("pippn_Momentum_%s",frame[f].Data()),pippnc_mom[f]);
    FillHist(Form("pippn_CosT_%s",frame[f].Data()),pippnc_cost[f]);
    FillHist(Form("pippn_Phi_%s",frame[f].Data()),pippnc_phi[f]);
    FillHist(Form("pippn_MMass_%s",frame[f].Data()),pippnc_mmass[f]);
    FillHist(Form("pippn_MMass2_%s",frame[f].Data()),pippnc_mmass2[f]);
    FillHist(Form("pippn_MMomentum_%s",frame[f].Data()),pippnc_mmom[f]);
    FillHist(Form("pippn_MCosT_%s",frame[f].Data()),pippnc_mcost[f]);
  }
  // pipp (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_Mass_n_%s",frame[f].Data()),pipp_mass[f]);
    FillHist(Form("pipp_Momentum_n_%s",frame[f].Data()),pipp_mom[f]);
    FillHist(Form("pipp_CosT_n_%s",frame[f].Data()),pipp_cost[f]);
    FillHist(Form("pipp_Phi_n_%s",frame[f].Data()),pipp_phi[f]);
    FillHist(Form("pipp_MMass_n_%s",frame[f].Data()),pipp_mmass[f]);
    FillHist(Form("pipp_MMomentum_n_%s",frame[f].Data()),pipp_mmom[f]);
    FillHist(Form("pipp_MCosT_n_%s",frame[f].Data()),pipp_mcost[f]);
    FillHist(Form("pipp_MTheta_n_%s",frame[f].Data()),TMath::ACos(pipp_mcost[f])*180.0/TMath::Pi());
  }
  // Momentum transfer (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_MomentumTransfer_n_%s",frame[f].Data()),q_2);
  }
  /*-----------------------------------------*/
  /* 2D Plots                                */ 
  /*-----------------------------------------*/
  // IM(pipp) vs. MM2(pippn)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpippn_MM2_%s",frame[f].Data()),pipp_mass[f],pippnc_mmass2[f]);
  }
  // IM(pipp) vs. MM(n)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKn_MM_%s",frame[f].Data()),pipp_mass[f],nc_mmass[f]);
  }
  // MMom(pipp) vs. Mom(n)
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpipp_MMomvsn_Mom_%s",frame[f].Data()),pipp_mom[f],nc_mom[f]);
  }
  // MM(pipp) vs. MM2(pippn)
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpipp_MMvsHeKpippn_MM2_%s",frame[f].Data()),pipp_mmass[f],pippnc_mmass2[f]);
  }
  // MM(pn) vs. MM2(pippn)
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpn_MMvsHeKpippn_MM2_%s",frame[f].Data()),pnc_mmass[f],pippnc_mmass2[f]);
  }
  // IM(pipp) vs. MM(pn) 
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpn_MM_%s",frame[f].Data()),pipp_mass[f],pnc_mmass[f]);
  }
  // MM(pipp) vs. MM(pn) 
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_MMvsHeKpn_MM_%s",frame[f].Data()),pipp_mass[f],pnc_mmass[f]);
  }
  // MM(pn) vs. MMom(pn) 
  for(int f=1; f<2; f++){
    FillHist(Form("pn_MMvspn_MMom_%s",frame[f].Data()),pnc_mmass[f],pnc_mmom[f]);
  }
  // MM(pipp) vs. MMom(pippn) 
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_MMvspippn_MMom_%s",frame[f].Data()),pipp_mmass[f],pippnc_mmom[f]);
  }
  // MM(n) vs. MM(pn) 
  for(int f=1; f<2; f++){
    FillHist(Form("HeKn_MMvsHeKpn_MM_%s",frame[f].Data()),nc_mmass[f],pnc_mmass[f]);
  }
  // MM2(pn) vs. MM2(pippn)
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpn_MM2vsHeKpippn_MM2_%s",frame[f].Data()),pnc_mmass2[f],pippnc_mmass2[f]);
  }
  // IM(pipp) vs. MM(pipp) (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MM_n_%s",frame[f].Data()),pipp_mass[f],pipp_mmass[f]);
  }
  // IM(pipp) vs. MCosT(pipp) (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MCosT_n_%s",frame[f].Data()),pipp_mass[f],pipp_mcost[f]);
  }
  // IM(pipp) vs. MomentumTransfer (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsMomentumTransfer_n_%s",frame[f].Data()),pipp_mass[f],q_2);
  }
  // Dalitz Plot (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_DalitzPlot_n_%s",frame[f].Data()),(tp-tl)/sqrt(3)/qvalue,tn/qvalue);
  }

  return true;

}

bool MyAnalysisHeKpippTwoDummy::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKpippTwoDummy::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKpippTwoDummy::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKpippTwoDummy::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKpippTwoDummy::CutCondition()
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

bool MyAnalysisHeKpippTwoDummy::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKpippTwoDummy::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKpippTwoDummy");

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
  // Rough pi-,proton1,proton2
  new TH1F("pipp_Mass_Rough","Invariant mass of #pi^{-}pp;IM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_MMass_Rough", "^{3}He(K^{-},#pi^{-}pp)X missing mass;MM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  // pi-
  new TH1F("pim_Momentum_CM","Momentum of #pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_CosT_CM","cos#theta of #pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim_Phi_CM","#phi of #pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim_MMass_CM", "^{3}He(K^{-},#pi^{-})X missing mass;MM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_MMomentum_CM", "^{3}He(K^{-},#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_MCosT_CM", "^{3}He(K^{-},#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton
  new TH1F("p_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton1
  new TH1F("p1_Momentum_CM","Momentum of p_{1};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p1_CosT_CM","cos#theta of p_{1};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p1_Phi_CM","#phi of p_{1};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p1_MMass_CM", "^{3}He(K^{-},p_{1})X missing mass;MM(p_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p1_MMomentum_CM", "^{3}He(K^{-},p_{1})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p1_MCosT_CM", "^{3}He(K^{-},p_{1})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton2
  new TH1F("p2_Momentum_CM","Momentum of p_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p2_CosT_CM","cos#theta of p_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p2_Phi_CM","#phi of p_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p2_MMass_CM", "^{3}He(K^{-},p_{2})X missing mass;MM(p_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p2_MMomentum_CM", "^{3}He(K^{-},p_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p2_MCosT_CM", "^{3}He(K^{-},p_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-,proton1
  new TH1F("pip1_Mass_CM","Invariant mass of #pi^{-}p_{1};IM(#pi^{-}p_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip1_Momentum_CM","Momentum of #pi^{-}p_{1};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip1_CosT_CM","cos#theta of #pi^{-}p_{1};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pip1_Phi_CM","#phi of #pi^{-}p_{1};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pip1_MMass_CM", "^{3}He(K^{-},#pi^{-}p_{1})X missing mass;MM(#pi^{-}p_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip1_MMomentum_CM", "^{3}He(K^{-},#pi^{-}p_{1})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip1_MCosT_CM", "^{3}He(K^{-},#pi^{-}p_{1})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-,proton2
  new TH1F("pip2_Mass_CM","Invariant mass of #pi^{-}p_{2};IM(#pi^{-}p_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip2_Momentum_CM","Momentum of #pi^{-}p_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip2_CosT_CM","cos#theta of #pi^{-}p_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pip2_Phi_CM","#phi of #pi^{-}p_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pip2_MMass_CM", "^{3}He(K^{-},#pi^{-}p_{2})X missing mass;MM(#pi^{-}p_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip2_MMomentum_CM", "^{3}He(K^{-},#pi^{-}p_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip2_MCosT_CM", "^{3}He(K^{-},#pi^{-}p_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Lambda
  new TH1F("lam_Mass_CM","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_Momentum_CM","Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_CosT_CM","cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_Phi_CM","#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("lam_MMass_CM", "^{3}He(K^{-},#Lambda)X missing mass;MM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_MMomentum_CM", "^{3}He(K^{-},#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_MCosT_CM", "^{3}He(K^{-},#Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Not Lambda
  new TH1F("nlam_Mass_CM","Invariant mass of Not #Lambda;IM(Not #Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nlam_Momentum_CM","Momentum of Not #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nlam_CosT_CM","cos#theta of Not #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("nlam_Phi_CM","#phi of Not #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("nlam_MMass_CM", "^{3}He(K^{-},Not #Lambda)X missing mass;MM(Not #Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nlam_MMomentum_CM", "^{3}He(K^{-},Not #Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nlam_MCosT_CM", "^{3}He(K^{-},Not #Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Proton
  new TH1F("proton_Mass_CM","Invariant mass of p;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("proton_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("proton_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("proton_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("proton_MMass_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("proton_MMomentum_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("proton_MCosT_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n
  new TH1F("n_Mass_CM","Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_Momentum_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_CM", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMass2_CM", "^{3}He(K^{-},n)X missing mass^{2};MM^{2}(n) (GeV/c^{2})^{2};Counts", 6000, -1.0, 5.0);
  new TH1F("n_MMomentum_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_MomentumTransfer_CM","Momentum transfer to n;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi-,proton1,proton2
  new TH1F("pipp_Mass_CM","Invariant mass of #pi^{-}pp;IM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_Momentum_CM","Momentum of #pi^{-}pp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipp_CosT_CM","cos#theta of #pi^{-}pp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipp_Phi_CM","#phi of #pi^{-}pp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pipp_MMass_CM", "^{3}He(K^{-},#pi^{-}pp)X missing mass;MM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_MMomentum_CM", "^{3}He(K^{-},#pi^{-}pp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipp_MCosT_CM", "^{3}He(K^{-},#pi^{-}pp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipp_MTheta_CM", "^{3}He(K^{-},#pi^{-}pp)X #theta;#theta;Counts", 200, -10.0, 190.0);
  new TH1F("pipp_MomentumTransfer_CM","Momentum transfer to #pi^{-}pp;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // Not pi-,proton1,proton2
  new TH1F("npipp_Mass_CM","Invariant mass of #pi^{-}pp;IM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npipp_Momentum_CM","Momentum of #pi^{-}pp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npipp_CosT_CM","cos#theta of #pi^{-}pp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npipp_Phi_CM","#phi of #pi^{-}pp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("npipp_MMass_CM", "^{3}He(K^{-},#pi^{-}pp)X missing mass;MM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npipp_MMomentum_CM", "^{3}He(K^{-},#pi^{-}pp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npipp_MCosT_CM", "^{3}He(K^{-},#pi^{-}pp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npipp_MTheta_CM", "^{3}He(K^{-},#pi^{-}pp)X #theta;#theta;Counts", 200, -10.0, 190.0);
  new TH1F("npipp_MomentumTransfer_CM","Momentum transfer to #pi^{-}pp;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi-,proton,neutron
  new TH1F("pipn_Mass_CM","Invariant mass of #pi^{-}pp;IM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipn_Momentum_CM","Momentum of #pi^{-}pp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipn_CosT_CM","cos#theta of #pi^{-}pp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipn_Phi_CM","#phi of #pi^{-}pp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pipn_MMass_CM", "^{3}He(K^{-},#pi^{-}pp)X missing mass;MM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipn_MMomentum_CM", "^{3}He(K^{-},#pi^{-}pp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipn_MCosT_CM", "^{3}He(K^{-},#pi^{-}pp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipn_MTheta_CM", "^{3}He(K^{-},#pi^{-}pp)X #theta;#theta;Counts", 200, -10.0, 190.0);
  new TH1F("pipn_MomentumTransfer_CM","Momentum transfer to #pi^{-}pp;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // proton + n
  new TH1F("pn_Mass_CM","Invariant mass of pn;IM(pn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_Momentum_CM","Momentum of pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_CosT_CM","cos#theta of pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_Phi_CM","#phi of pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pn_MMass_CM", "^{3}He(K^{-},pn)X missing mass;MM(pn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_MMass2_CM", "^{3}He(K^{-},pn)X missing mass^{2};MM^{2}(pn) (GeV/c^{2})^{2};Counts", 6000, -1.0, 2.0);
  new TH1F("pn_MMomentum_CM", "^{3}He(K^{-},pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_MCosT_CM", "^{3}He(K^{-},pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_MomentumTransfer_CM","Momentum transfer to pn;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi-,proton1,proton2 + n
  new TH1F("pippn_Mass_CM","Invariant mass of #pi^{-}ppn;IM(#pi^{-}ppn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pippn_Momentum_CM","Momentum of #pi^{-}ppn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pippn_CosT_CM","cos#theta of #pi^{-}ppn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pippn_Phi_CM","#phi of #pi^{-}ppn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pippn_MMass_CM", "^{3}He(K^{-},#pi^{-}ppn)X missing mass;MM(#pi^{-}ppn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pippn_MMass2_CM", "^{3}He(K^{-},#pi^{-}ppn)X missing mass^{2};MM^{2}(#pi^{-}ppn) (GeV/c^{2})^{2};Counts", 6000, -1.0, 2.0);
  new TH1F("pippn_MMomentum_CM", "^{3}He(K^{-},#pi^{-}ppn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pippn_MCosT_CM", "^{3}He(K^{-},#pi^{-}ppn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pippn_MomentumTransfer_CM","Momentum transfer to #pi^{-}ppn;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi-,proton1,proton2 (FWD)
  new TH1F("pipp_Mass_n_CM","Invariant mass of #pi^{-}pp;IM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_Momentum_n_CM","Momentum of #pi^{-}pp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipp_CosT_n_CM","cos#theta of #pi^{-}pp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipp_Phi_n_CM","#phi of #pi^{-}pp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pipp_MMass_n_CM", "^{3}He(K^{-},#pi^{-}pp)X missing mass;MM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_MMomentum_n_CM", "^{3}He(K^{-},#pi^{-}pp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipp_MCosT_n_CM", "^{3}He(K^{-},#pi^{-}pp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipp_MTheta_n_CM", "^{3}He(K^{-},#pi^{-}pp)X #theta;#theta;Counts", 200, -10.0, 190.0);
  new TH1F("pipp_MomentumTransfer_n_CM","Momentum transfer to #pi^{-}pp;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // Opening Angle
  new TH1F("pipp_OA_CM","Opening Angle between #Lambda and p;cos#theta_{#Lambdap};Counts", 3000, -1.5, 1.5);
  /* 2D Plots */
  new TH2F("pip1_IMvspip2_IM_All", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 3.0, 2000, 1.0, 3.0);
  new TH2F("pip1_IMvspip2_IM_CM", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 3.0, 2000, 1.0, 3.0);
  new TH2F("pip1_IM2vspip2_IM2_All", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 5.0, 2000, 1.0, 5.0);
  new TH2F("pip1_IM2vspip2_IM2_CM", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 5.0, 2000, 1.0, 5.0);
  new TH2F("pipp_IMvsHeKpipp_MM_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2})", 200, 2.0, 4.0, 200, 0.0, 2.0);
  new TH2F("pipp_IMvsHeKpipp_MCosT_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. cos#theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta_{n}", 200, 2.0, 4.0, 300, -1.5, 1.5);
  new TH2F("pipp_IMvsHeKpipp_MTheta_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. #theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});#theta_{n}", 200, 2.0, 4.0, 200, -10.5, 190.0);
  new TH2F("pipp_IMvsMomentumTransfer_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. Momentum transfer;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Momentum Transfer (GeV/c)", 200, 2.0, 4.0, 200, 0.0, 2.0);
  new TH2F("pipp_IMvspipp_MMom_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. Missing momentum;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Missing momentum (GeV/c)", 200, 2.0, 4.0, 200, 0.0, 2.0);
  new TH2F("pipp_MMvspipp_MMom_CM", "(#pi^{-}p_{1}p_{2}) missing mass vs. Missing momentum;MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Missing momentum (GeV/c)", 200, 0.0, 2.0, 200, 0.0, 2.0);
  new TH2F("pipp_IMvsOA_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. Opening Angle;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta_{#Lambdap} ", 200, 2.0, 4.0, 300, -1.5, 1.5);
  new TH2F("pipp_DalitzPlot_CM", "Dalitz Plot;(T_{p}-T_{#Lambda})/#sqrt{3}Q;T_{n}/Q", 1000, -1.0, 1.0, 1000, -1.0, 1.0);
  new TH2F("pipp_IMvslam_PDF_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. PDF;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});PDf", 2000, 2.0, 4.0, 400, 0.0, 40.0);
  new TH2F("HeKpipp_MMvslam_PDF_CM", "^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass vs. PDF;MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});PDf", 2000, 0.0, 2.0, 400, 0.0, 40.0);
  new TH2F("npipp_IMvsHeKnpipp_MM_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2})", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pipp_IM2vspipn_IM2_CM", "(#Lambdap) invariant mass square vs. (#Lambdan) invariant mass square;IM^{2}(#Lambdap) (GeV^{2}/c^{4});IM^{2}(#Lambdan) (GeV^{2}/c^{4});", 300, 4.0, 10.0, 300, 4.0, 10.0);
  new TH2F("pipp_Mass_Rough_vs_CM","2D Invariant mass of #pi^{-}pp;IM(#pi^{-}pp) (GeV/c^{2});Counts", 500, 2.0, 3.0, 500, 2.0, 3.0);
  new TH2F("pipp_MMass_Rough_vs_CM", "2D ^{3}He(K^{-},#pi^{-}pp)X missing mass;MM(#pi^{-}pp) (GeV/c^{2});Counts", 1000, 0.0, 2.0, 1000, 0.0, 2.0);
  /* 2D Plots  (FWD)*/
  new TH2F("pipp_IMvsHeKpippn_MM2_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2}n)X missing mass^{2};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM^{2}(#pi^{-}p_{1}p_{2}n) (GeV/c^{2})^{2}", 2000, 2.0, 4.0, 2000, -1.0, 1.0);
  new TH2F("pipp_IMvsHeKn_MM_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},n)X missing mass^{2};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(n) (GeV/c^{2})", 200, 2.0, 4.0, 200, 2.0, 4.0);
  new TH2F("HeKpipp_MMomvsn_Mom_CM", "(#pi^{-}p_{1}p_{2}) momentum vs. n momentum;Missing Momentum (#pi^{-}p_{1}p_{2}) (GeV/c);Momentum (n) (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("HeKpipp_MMvsHeKpippn_MM2_CM", "^{3}He(K^{-},(#pi^{-}p_{1}p_{2})X missing mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2}n)X missing mass^{2};MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM^{2}(#pi^{-}p_{1}p_{2}n) (GeV/c^{2})^{2}", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("HeKpn_MMvsHeKpippn_MM2_CM", "^{3}He(K^{-},(pn)X missing mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2}n)X missing mass^{2};MM(pn) (GeV/c^{2});MM^{2}(#pi^{-}p_{1}p_{2}n) (GeV/c^{2})^{2}", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("HeKpn_MM2vsHeKpippn_MM2_CM", "^{3}He(K^{-},(pn)X missing mass^{2} vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2}n)X missing mass^{2};MM^{2}(pn) (GeV/c^{2});MM^{2}(#pi^{-}p_{1}p_{2}n) (GeV/c^{2})^{2}", 4000, 0.0, 4.0, 2000, -1.0, 1.0);
  new TH2F("pipp_IMvsHeKpn_MM_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},pn)X missing mass^{2};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(pn) (GeV/c^{2})", 200, 2.0, 4.0, 100, 1.0, 2.0);
  new TH2F("pipp_MMvsHeKpn_MM_CM", "^{3}He(K^{-},pi^{-}p_{1}p_{2})X missing mass vs. ^{3}He(K^{-},pn)X missing mass^{2};MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(pn) (GeV/c^{2})", 200, 2.0, 4.0, 100, 1.0, 2.0);
  new TH2F("HeKn_MMvsHeKpn_MM_CM", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},pn)X missing mass^{2};MM(n) (GeV/c^{2});MM^{2}(pn) (GeV/c^{2})^{2}", 2000, 2.0, 4.0, 1000, 1.0, 2.0);
  new TH2F("pn_MMvspn_MMom_CM", "^{3}He(K^{-},pn)X missing mass vs. ^{3}He(K^{-},pn)X missing momentum;MM(pn) (GeV/c^{2});Missin Momentum(pn) (GeV/c)", 200, 0.5, 2.5, 200, 0.0, 2.0);
  new TH2F("pipp_MMvspippn_MMom_CM", "^{3}He(K^{-},#pi^{-}pp)X missing mass vs. ^{3}He(K^{-},#pi^{-}ppn)X missing momentum;MM(#pi^{-}pp) (GeV/c^{2});Missin Momentum(#pi^{-}ppn) (GeV/c)", 140, 0.4, 1.8, 200, 0.0, 2.0);
  new TH2F("pipp_IMvsHeKpipp_MM_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2})", 200, 2.0, 4.0, 200, 0.0, 2.0);
  new TH2F("pipp_IMvsHeKpipp_MCosT_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. cos#theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta_{n}", 200, 2.0, 4.0, 300, -1.5, 1.5);
  new TH2F("pipp_IMvsMomentumTransfer_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. Momentum transfer;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Momentum Transfer (GeV/c)", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pipp_DalitzPlot_n_CM", "Dalitz Plot;(T_{p}-T_{#Lambda})/#sqrt{3}Q;T_{n}/Q", 1000, -1.0, 1.0, 1000, -1.0, 1.0);

  std::cout << "== Finish Histogram Initialization ==" << std::endl;
  return true;

}
