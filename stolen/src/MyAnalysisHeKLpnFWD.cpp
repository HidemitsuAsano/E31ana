#include "MyAnalysisHeKLpnFWD.h"

MyAnalysisHeKLpnFWD::MyAnalysisHeKLpnFWD(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKLpnFWD::~MyAnalysisHeKLpnFWD()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKLpnFWD::Clear()
{
}

bool MyAnalysisHeKLpnFWD::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  int ievent = 0;

  FillHist("EventNumber",ievent); ievent++; /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  /* Event selection */
  if(particle->nCDS()!=2&&particle->nCDS()!=3) return true;
  FillHist("EventNumber",ievent); ievent++; /* 2 or 3 hits in CDS*/
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(cds->chi()>30.0) return true;
  }
  FillHist("EventNumber",ievent); ievent++; /* CDS Chi-square Cut */

  bool neutral=false; bool charged=false;
  int MulCVC=0;
  for(int i=0; i<blMan->nCVC(); i++){
    HodoscopeLikeHit* hit = blMan->CVC(i);
    if(hit->CheckRange()) MulCVC++;
  }
  if(header->IsTrig(Trig_Neutral)&&particle->nNC()>0) neutral = true;
  if(header->IsTrig(Trig_Charged)&&particle->nPC()>0&&MulCVC==0) charged = true;
  if(neutral&&charged) return true;
  FillHist("EventNumber",ievent); ievent++; /* Both Neutral and Charged */
  if(neutral){
  if(particle->nPiminus()!=1) return true;
  if(particle->nProton()==0) return true;
  if(particle->nKaon()!=0) return true;
  if(particle->nDeuteron()!=0) return true;
  if(particle->nOther()!=0) return true;
  }
  if(charged){
    if(particle->nPiminus()!=1) return true;
    if(particle->nProton()!=1) return true;
    if(particle->nKaon()!=0) return true;
    if(particle->nDeuteron()!=0) return true;
    if(particle->nOther()!=0) return true;
  }
  FillHist("EventNumber",ievent); ievent++; /* pi-&p in Neutral / pi-&p in Charged */

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  double pdf = 99999;
  int npro = 0;
  for(int it=0; it<particle->nProduct(); it++){
    pCDS* product = particle->product(it);
    int comb = product->comb();
    if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
      npro++;
      TVector3 tmpvertex = product->vbeam();
      if(GeomTools::GetID(tmpvertex)==CID_Fiducial){
        double tmpper[5] = {product->mass(),product->vdis(),product->pbdca(),0.0,0.0};
        double tmpprob[5] = {0,0,0,0,0};
        TrackTools::PDFLambda(tmpper,tmpprob);
        double tmppdf = -TMath::Log(tmpprob[0]*tmpprob[1]*tmpprob[2]);
        if(pdf>tmppdf){
          pdf = tmppdf;
          fcds = it;
          vertex = tmpvertex;
        }
      }
    }
  }
  FillHist("NumberOfProduct",npro);

  if(fcds==-1) return true; /* Vertex decision */
  pCDS* lam=0;
  if(fcds!=-1) lam = particle->product(fcds);
  FillHist("EventNumber",ievent); ievent++; /* Lambda Candidate */

  FillHist("lam_Mass",lam->mass());
  FillHist("lam_PDF",pdf);
  FillHist("lam_MassvsPDF",lam->mass(),pdf);
  double pdful = 9.0;

  if(pdf>pdful) return true;
  FillHist("EventNumber",ievent); ievent++; /* Lambda PDF Cut */
  FillHist("lam_Mass_CutPDF",lam->mass());
  FillHist("lam_PDF_CutPDF",pdf);
  FillHist("lam_MassvsPDF_CutPDF",lam->mass(),pdf);

  if(neutral) DoNeutralAna(particle,beam,lam,vertex);
  if(charged) DoChargedAna(particle,beam,lam,vertex);

  return true;

}

bool MyAnalysisHeKLpnFWD::DoNeutralAna(Particle* particle, pBeam* beam, pCDS* lam, TVector3 vertex)
{
  DetectorList *dlist=DetectorList::GetInstance();
  TString mode = "Neutral";

  pCDS* p = 0;
  pCDS* pim = 0;

  if(particle->cdsi(lam->daughter1())->pid()==CDS_Proton){
      p = particle->cdsi(lam->daughter1());
      pim = particle->cdsi(lam->daughter2());
  }
  else {
      pim = particle->cdsi(lam->daughter1());
      p = particle->cdsi(lam->daughter2());
  }

  // ############### //
  // NC hit decision //
  // ############### //
  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particle->nNC(); it++){
    pNC* nc = particle->nc(it);
    nc->CalcMom(beam,vertex);
    double seg = (nc->seg()-1)%16+1;
    double lay = (nc->seg()-1)/16+1;
    FillHist("FWDN_Overbeta",1.0/nc->beta());
    FillHist("FWDN_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("NCHitPosition_XY",nc->hitpos().X(),nc->hitpos().Y());
    FillHist("NCHitSegment",seg,lay);
    if(nc->energy()>8.0 && (time>9998||time>nc->time())){
      time = nc->time();
      fnc = it;
      FillHist("FWDN_Overbeta_withdECut",1.0/nc->beta());
      FillHist("FWDN_OverbetavsEnergy_withdECut",1.0/nc->beta(),nc->energy());
    }
  }
  if(fnc==-1) return true;
  pNC*  nc =0;
  if(fnc!=-1) nc = particle->nc(fnc);
  FillHist("FWDN_Overbeta_Selected",1.0/nc->beta());
  if(nc->pid()!=F_Neutron){ return true; }

  TString frame[2] = {"Lab","CM"};

  // Target
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
  if(dlist->GetMaterial(CID_Target)=="Vacuum"){
    Ltgt.SetVectM(ZeroV,ThreeHeMass);
  }
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* Proton */
  TLorentzVector Lp = p->GetLorentzVector();
  TLorentzVector cLp = p->GetLorentzVector(); cLp.Boost(-boost);
  double p_mass[2]  = { (Lp).M()                    , (cLp).M()};
  double p_mom[2]   = { (Lp).Vect().Mag()           , (cLp).Vect().Mag()};
  double p_cost[2]  = { (Lp).Vect().CosTheta()      , (cLp).Vect().CosTheta()};
  double p_phi[2]   = { (Lp).Vect().Phi()                  , (cLp).Vect().Phi()};
  double p_mmass[2]   = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
  double p_mmom[2]    = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
  double p_mcost[2]   = { (Ltgt+Lbeam-Lp).Vect().CosTheta(), (cLtgt+cLbeam-cLp).Vect().CosTheta()};
  /* Pi- */
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  double pim_mass[2]  = { (Lpim).M()                    , (cLpim).M()};
  double pim_mom[2]   = { (Lpim).Vect().Mag()           , (cLpim).Vect().Mag()};
  double pim_cost[2]  = { (Lpim).Vect().CosTheta()      , (cLpim).Vect().CosTheta()};
  double pim_phi[2]   = { (Lpim).Vect().Phi()                  , (cLpim).Vect().Phi()};
  double pim_mmass[2]   = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]    = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2]   = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};
  /* Lambda */
  TLorentzVector Llam = lam->GetLorentzVector();
  TLorentzVector cLlam = lam->GetLorentzVector(); cLlam.Boost(-boost);
  double lam_mass[2]  = { (Llam).M()                    , (cLlam).M()};
  double lam_mom[2]   = { (Llam).Vect().Mag()           , (cLlam).Vect().Mag()};
  double lam_cost[2]  = { (Llam).Vect().CosTheta()      , (cLlam).Vect().CosTheta()};
  double lam_phi[2]   = { (Llam).Vect().Phi()                  , (cLlam).Vect().Phi()};
  double lam_mmass[2]   = { (Ltgt+Lbeam-Llam).M()              , (cLtgt+cLbeam-cLlam).M()};
  double lam_mmom[2]    = { (Ltgt+Lbeam-Llam).Vect().Mag()     , (cLtgt+cLbeam-cLlam).Vect().Mag()};
  double lam_mcost[2]   = { (Ltgt+Lbeam-Llam).Vect().CosTheta(), (cLtgt+cLbeam-cLlam).Vect().CosTheta()};
  /* NC */
  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
  double nc_mass[2]    = { Ln.M()                           , cLn.M()};
  double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};
  /* NC + Proton */
  double np_mass[2]  = { (Lp+Ln).M()                           , (cLp+cLn).M()};
  double np_mom[2]   = { (Lp+Ln).Vect().Mag()                  , (cLp+cLn).Vect().Mag()};
  double np_cost[2]  = { (Lp+Ln).Vect().CosTheta()             , (cLp+cLn).Vect().CosTheta()};
  double np_phi[2]   = { (Lp+Ln).Vect().Phi()                  , (cLp+cLn).Vect().Phi()};
  double np_mmass[2] = { (Ltgt+Lbeam-(Lp+Ln)).M()              , (cLtgt+cLbeam-(cLp+cLn)).M()};
  double np_mmom[2]  = { (Ltgt+Lbeam-(Lp+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLp+cLn)).Vect().Mag()};
  double np_mcost[2] = { (Ltgt+Lbeam-(Lp+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLp+cLn)).Vect().CosTheta()};
  /* NC + Pi- */
  double npim_mass[2]  = { (Lpim+Ln).M()                           , (cLpim+cLn).M()};
  double npim_mom[2]   = { (Lpim+Ln).Vect().Mag()                  , (cLpim+cLn).Vect().Mag()};
  double npim_cost[2]  = { (Lpim+Ln).Vect().CosTheta()             , (cLpim+cLn).Vect().CosTheta()};
  double npim_phi[2]   = { (Lpim+Ln).Vect().Phi()                  , (cLpim+cLn).Vect().Phi()};
  double npim_mmass[2] = { (Ltgt+Lbeam-(Lpim+Ln)).M()              , (cLtgt+cLbeam-(cLpim+cLn)).M()};
  double npim_mmom[2]  = { (Ltgt+Lbeam-(Lpim+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim+cLn)).Vect().Mag()};
  double npim_mcost[2] = { (Ltgt+Lbeam-(Lpim+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim+cLn)).Vect().CosTheta()};
  /* NC + Lambda */
  double nlam_mass[2]  = { (Llam+Ln).M()                           , (cLlam+cLn).M()};
  double nlam_mom[2]   = { (Llam+Ln).Vect().Mag()                  , (cLlam+cLn).Vect().Mag()};
  double nlam_cost[2]  = { (Llam+Ln).Vect().CosTheta()             , (cLlam+cLn).Vect().CosTheta()};
  double nlam_phi[2]   = { (Llam+Ln).Vect().Phi()                  , (cLlam+cLn).Vect().Phi()};
  double nlam_mmass[2] = { (Ltgt+Lbeam-(Llam+Ln)).M()              , (cLtgt+cLbeam-(cLlam+cLn)).M()};
  double nlam_mmom[2]  = { (Ltgt+Lbeam-(Llam+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam+cLn)).Vect().Mag()};
  double nlam_mcost[2] = { (Ltgt+Lbeam-(Llam+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam+cLn)).Vect().CosTheta()};

  double tp = (cLtgt+cLbeam-cLlam-cLn).E() - (cLtgt+cLbeam-cLlam-cLn).M();
  double tn = cLn.E() - cLn.M();
  double tl = cLlam.E() - cLlam.M();
  double qvalue = tn+tp+tl;

  /* flag */
  bool mpflag = false; /* missing p - flag */
  bool sblmpflag = false; /* lower sideband of missing p - flag */
  bool sbumpflag = false; /* upper sideband of missing p - flag */
  {
    double pro_mmass = nlam_mmass[1];
    if(mpll<pro_mmass&&pro_mmass<mpul) mpflag = true;
    if(sblmpll<pro_mmass&&pro_mmass<sblmpul) sblmpflag = true;
    if(sbumpll<pro_mmass&&pro_mmass<sbumpul) sbumpflag = true;
  }

  if(!mpflag){ return true; }
  //if(k0flag){ return true; }
  //if(spflag||smflag){ return true; }
  //if(!spflag){ return true; }
  //if(!smflag){ return true; }
  //if(!spflag&&!smflag){ return true; }

  /* Basic */
  // Proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Mass_%s_%s",frame[f].Data(),mode.Data()),p_mass[f]);
    FillHist(Form("p_Momentum_%s_%s",frame[f].Data(),mode.Data()),p_mom[f]);
    FillHist(Form("p_CosT_%s_%s",frame[f].Data(),mode.Data()),p_cost[f]);
    FillHist(Form("p_Phi_%s_%s",frame[f].Data(),mode.Data()),p_phi[f]);
    FillHist(Form("p_MMass_%s_%s",frame[f].Data(),mode.Data()),p_mmass[f]);
    FillHist(Form("p_MMomentum_%s_%s",frame[f].Data(),mode.Data()),p_mmom[f]);
    FillHist(Form("p_MCosT_%s_%s",frame[f].Data(),mode.Data()),p_mcost[f]);
  }
  // Pi-
  for(int f=1; f<2; f++){
    FillHist(Form("pim_Mass_%s_%s",frame[f].Data(),mode.Data()),pim_mass[f]);
    FillHist(Form("pim_Momentum_%s_%s",frame[f].Data(),mode.Data()),pim_mom[f]);
    FillHist(Form("pim_CosT_%s_%s",frame[f].Data(),mode.Data()),pim_cost[f]);
    FillHist(Form("pim_Phi_%s_%s",frame[f].Data(),mode.Data()),pim_phi[f]);
    FillHist(Form("pim_MMass_%s_%s",frame[f].Data(),mode.Data()),pim_mmass[f]);
    FillHist(Form("pim_MMomentum_%s_%s",frame[f].Data(),mode.Data()),pim_mmom[f]);
    FillHist(Form("pim_MCosT_%s_%s",frame[f].Data(),mode.Data()),pim_mcost[f]);
  }
  // Lambda
  for(int f=1; f<2; f++){
    FillHist(Form("lam_Mass_%s_%s",frame[f].Data(),mode.Data()),lam_mass[f]);
    FillHist(Form("lam_Momentum_%s_%s",frame[f].Data(),mode.Data()),lam_mom[f]);
    FillHist(Form("lam_CosT_%s_%s",frame[f].Data(),mode.Data()),lam_cost[f]);
    FillHist(Form("lam_Phi_%s_%s",frame[f].Data(),mode.Data()),lam_phi[f]);
    FillHist(Form("lam_MMass_%s_%s",frame[f].Data(),mode.Data()),lam_mmass[f]);
    FillHist(Form("lam_MMomentum_%s_%s",frame[f].Data(),mode.Data()),lam_mmom[f]);
    FillHist(Form("lam_MCosT_%s_%s",frame[f].Data(),mode.Data()),lam_mcost[f]);
  }
  // n
  for(int f=1; f<2; f++){
    FillHist(Form("nfwd_Mass_%s_%s",frame[f].Data(),mode.Data()),nc_mass[f]);
    FillHist(Form("nfwd_Momentum_%s_%s",frame[f].Data(),mode.Data()),nc_mom[f]);
    FillHist(Form("nfwd_CosT_%s_%s",frame[f].Data(),mode.Data()),nc_cost[f]);
    FillHist(Form("nfwd_Phi_%s_%s",frame[f].Data(),mode.Data()),nc_phi[f]);
    FillHist(Form("nfwd_MMass_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f]);
    FillHist(Form("nfwd_MMomentum_%s_%s",frame[f].Data(),mode.Data()),nc_mmom[f]);
    FillHist(Form("nfwd_MCosT_%s_%s",frame[f].Data(),mode.Data()),nc_mcost[f]);
  }
  // n + Proton
  for(int f=1; f<2; f++){
    FillHist(Form("np_Mass_%s_%s",frame[f].Data(),mode.Data()),np_mass[f]);
    FillHist(Form("np_Momentum_%s_%s",frame[f].Data(),mode.Data()),np_mom[f]);
    FillHist(Form("np_CosT_%s_%s",frame[f].Data(),mode.Data()),np_cost[f]);
    FillHist(Form("np_Phi_%s_%s",frame[f].Data(),mode.Data()),np_phi[f]);
    FillHist(Form("np_MMass_%s_%s",frame[f].Data(),mode.Data()),np_mmass[f]);
    FillHist(Form("np_MMomentum_%s_%s",frame[f].Data(),mode.Data()),np_mmom[f]);
    FillHist(Form("np_MCosT_%s_%s",frame[f].Data(),mode.Data()),np_mcost[f]);
  }
  // n + pim
  for(int f=1; f<2; f++){
    FillHist(Form("npim_Mass_%s_%s",frame[f].Data(),mode.Data()),npim_mass[f]);
    FillHist(Form("npim_Momentum_%s_%s",frame[f].Data(),mode.Data()),npim_mom[f]);
    FillHist(Form("npim_CosT_%s_%s",frame[f].Data(),mode.Data()),npim_cost[f]);
    FillHist(Form("npim_Phi_%s_%s",frame[f].Data(),mode.Data()),npim_phi[f]);
    FillHist(Form("npim_MMass_%s_%s",frame[f].Data(),mode.Data()),npim_mmass[f]);
    FillHist(Form("npim_MMomentum_%s_%s",frame[f].Data(),mode.Data()),npim_mmom[f]);
    FillHist(Form("npim_MCosT_%s_%s",frame[f].Data(),mode.Data()),npim_mcost[f]);
  }
  // n + Lambda
  for(int f=1; f<2; f++){
    FillHist(Form("nlam_Mass_%s_%s",frame[f].Data(),mode.Data()),nlam_mass[f]);
    FillHist(Form("nlam_Momentum_%s_%s",frame[f].Data(),mode.Data()),nlam_mom[f]);
    FillHist(Form("nlam_CosT_%s_%s",frame[f].Data(),mode.Data()),nlam_cost[f]);
    FillHist(Form("nlam_Phi_%s_%s",frame[f].Data(),mode.Data()),nlam_phi[f]);
    FillHist(Form("nlam_MMass_%s_%s",frame[f].Data(),mode.Data()),nlam_mmass[f]);
    FillHist(Form("nlam_MMomentum_%s_%s",frame[f].Data(),mode.Data()),nlam_mmom[f]);
    FillHist(Form("nlam_MCosT_%s_%s",frame[f].Data(),mode.Data()),nlam_mcost[f]);
  }

  /*-----------------------------------------*/
  /* 2D Plots                                */ 
  /*-----------------------------------------*/
  // IM(nlam) vs. MM(nlam)
  for(int f=1; f<2; f++){
    FillHist(Form("nlam_IM_vs_nlam_MM_%s_%s",frame[f].Data(),mode.Data()),nlam_mass[f],nlam_mmass[f]);
  }
  // MM(n) vs. MM(nlam)
  for(int f=1; f<2; f++){
    FillHist(Form("nfwd_MM_vs_nlam_MM_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],nlam_mmass[f]);
  }
  // MM(n) vs. IM(nlam)
  for(int f=1; f<2; f++){
    FillHist(Form("nfwd_MM_vs_nlam_IM_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],nlam_mass[f]);
  }
  // MM(n) vs. CosT(n)
  for(int f=1; f<2; f++){
    FillHist(Form("nfwd_MM_vs_nfwd_CosT_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],nc_cost[f]);
  }
  // MM(n) vs. MMom(n)
  for(int f=1; f<2; f++){
    FillHist(Form("nfwd_MM_vs_nfwd_MMom_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],nc_mmom[f]);
  }
  // MM(n) vs. MM(np)
  for(int f=1; f<2; f++){
    FillHist(Form("nfwd_MM_vs_np_MM_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],np_mmass[f]);
  }
  // Dalitz Plot
  for(int f=1; f<2; f++){
    FillHist(Form("DalitzPlot_%s_%s",frame[f].Data(),mode.Data()),(tp-tl)/sqrt(3)/qvalue,tn/qvalue);
  }

  return true;

}

bool MyAnalysisHeKLpnFWD::DoChargedAna(Particle* particle, pBeam* beam, pCDS* lam, TVector3 vertex)
{
  DetectorList *dlist=DetectorList::GetInstance();
  TString mode = "Charged";

  pCDS* p = 0;
  pCDS* pim = 0;

  if(particle->cdsi(lam->daughter1())->pid()==CDS_Proton){
      p = particle->cdsi(lam->daughter1());
      pim = particle->cdsi(lam->daughter2());
  }
  else {
      pim = particle->cdsi(lam->daughter1());
      p = particle->cdsi(lam->daughter2());
  }

  // ##################### //
  // FWD Charged  Analysis //
  // ##################### //
  for(int ipc=0; ipc<particle->nPC(); ipc++){
    pPC* pc = particle->pc(ipc);
    int pid = -1;
    pc->CalcMom(beam,vertex,pid,false);
    double r = pc->r();
    double mass2 = pc->mass2();
    double momr = pc->momr();
    double angle = pc->angle();
    double mom = pc->mom().Mag();
    double tof = pc->tof();
    double fl = pc->fl();
    double energy = pc->energy();
    double beta = pc->beta();

    FillHist("FWDC_Beta",beta);
    FillHist("FWDC_Overbeta",1.0/beta);
    FillHist("FWDC_OverbetavsEnergy",1.0/beta,energy);
    FillHist("FWDC_MomentumR",momr);
    FillHist("FWDC_Momentum",mom);
    FillHist("FWDC_MomentumvsMomentumR",mom,momr);
    FillHist("FWDC_Mass2",mass2);
    FillHist("FWDC_Mass2vsMomentum",mass2,momr);
  }

  if(particle->nPC()!=1) return true;
  pPC* pc = 0;
  pc = particle->pc(0);
  if(pc->pid()!=F_Proton) return true;

  TString frame[2] = {"Lab","CM"};

  // Target
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
  if(dlist->GetMaterial(CID_Target)=="Vacuum"){
    Ltgt.SetVectM(ZeroV,ThreeHeMass);
  }
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* Proton */
  TLorentzVector Lp = p->GetLorentzVector();
  TLorentzVector cLp = p->GetLorentzVector(); cLp.Boost(-boost);
  double p_mass[2]  = { (Lp).M()                    , (cLp).M()};
  double p_mom[2]   = { (Lp).Vect().Mag()           , (cLp).Vect().Mag()};
  double p_cost[2]  = { (Lp).Vect().CosTheta()      , (cLp).Vect().CosTheta()};
  double p_phi[2]   = { (Lp).Vect().Phi()                  , (cLp).Vect().Phi()};
  double p_mmass[2]   = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
  double p_mmom[2]    = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
  double p_mcost[2]   = { (Ltgt+Lbeam-Lp).Vect().CosTheta(), (cLtgt+cLbeam-cLp).Vect().CosTheta()};
  /* Pi- */
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  double pim_mass[2]  = { (Lpim).M()                    , (cLpim).M()};
  double pim_mom[2]   = { (Lpim).Vect().Mag()           , (cLpim).Vect().Mag()};
  double pim_cost[2]  = { (Lpim).Vect().CosTheta()      , (cLpim).Vect().CosTheta()};
  double pim_phi[2]   = { (Lpim).Vect().Phi()                  , (cLpim).Vect().Phi()};
  double pim_mmass[2]   = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]    = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2]   = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};
  /* Lambda */
  TLorentzVector Llam = lam->GetLorentzVector();
  TLorentzVector cLlam = lam->GetLorentzVector(); cLlam.Boost(-boost);
  double lam_mass[2]  = { (Llam).M()                    , (cLlam).M()};
  double lam_mom[2]   = { (Llam).Vect().Mag()           , (cLlam).Vect().Mag()};
  double lam_cost[2]  = { (Llam).Vect().CosTheta()      , (cLlam).Vect().CosTheta()};
  double lam_phi[2]   = { (Llam).Vect().Phi()                  , (cLlam).Vect().Phi()};
  double lam_mmass[2]   = { (Ltgt+Lbeam-Llam).M()              , (cLtgt+cLbeam-cLlam).M()};
  double lam_mmom[2]    = { (Ltgt+Lbeam-Llam).Vect().Mag()     , (cLtgt+cLbeam-cLlam).Vect().Mag()};
  double lam_mcost[2]   = { (Ltgt+Lbeam-Llam).Vect().CosTheta(), (cLtgt+cLbeam-cLlam).Vect().CosTheta()};
  /* PC */
  TLorentzVector Lpc  = pc ->GetLorentzVector();
  TLorentzVector cLpc  = pc ->GetLorentzVector(); cLpc.Boost(-boost);
  double pc_mass[2]    = { Lpc.M()                           , cLpc.M()};
  double pc_mom[2]     = { Lpc.Vect().Mag()                  , cLpc.Vect().Mag()};
  double pc_cost[2]    = { Lpc.Vect().CosTheta()             , cLpc.Vect().CosTheta()};
  double pc_phi[2]     = { Lpc.Vect().Phi()                  , cLpc.Vect().Phi()};
  double pc_mmass[2]   = { (Ltgt+Lbeam-Lpc).M()              , (cLtgt+cLbeam-cLpc).M()};
  double pc_mmom[2]    = { (Ltgt+Lbeam-Lpc).Vect().Mag()     , (cLtgt+cLbeam-cLpc).Vect().Mag()};
  double pc_mcost[2]   = { (Ltgt+Lbeam-Lpc).Vect().CosTheta(), (cLtgt+cLbeam-cLpc).Vect().CosTheta()};
  /* PC + Proton */
  double pp_mass[2]  = { (Lp+Lpc).M()                           , (cLp+cLpc).M()};
  double pp_mom[2]   = { (Lp+Lpc).Vect().Mag()                  , (cLp+cLpc).Vect().Mag()};
  double pp_cost[2]  = { (Lp+Lpc).Vect().CosTheta()             , (cLp+cLpc).Vect().CosTheta()};
  double pp_phi[2]   = { (Lp+Lpc).Vect().Phi()                  , (cLp+cLpc).Vect().Phi()};
  double pp_mmass[2] = { (Ltgt+Lbeam-(Lp+Lpc)).M()              , (cLtgt+cLbeam-(cLp+cLpc)).M()};
  double pp_mmom[2]  = { (Ltgt+Lbeam-(Lp+Lpc)).Vect().Mag()     , (cLtgt+cLbeam-(cLp+cLpc)).Vect().Mag()};
  double pp_mcost[2] = { (Ltgt+Lbeam-(Lp+Lpc)).Vect().CosTheta(), (cLtgt+cLbeam-(cLp+cLpc)).Vect().CosTheta()};
  /* PC + Pi- */
  double ppim_mass[2]  = { (Lpim+Lpc).M()                           , (cLpim+cLpc).M()};
  double ppim_mom[2]   = { (Lpim+Lpc).Vect().Mag()                  , (cLpim+cLpc).Vect().Mag()};
  double ppim_cost[2]  = { (Lpim+Lpc).Vect().CosTheta()             , (cLpim+cLpc).Vect().CosTheta()};
  double ppim_phi[2]   = { (Lpim+Lpc).Vect().Phi()                  , (cLpim+cLpc).Vect().Phi()};
  double ppim_mmass[2] = { (Ltgt+Lbeam-(Lpim+Lpc)).M()              , (cLtgt+cLbeam-(cLpim+cLpc)).M()};
  double ppim_mmom[2]  = { (Ltgt+Lbeam-(Lpim+Lpc)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim+cLpc)).Vect().Mag()};
  double ppim_mcost[2] = { (Ltgt+Lbeam-(Lpim+Lpc)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim+cLpc)).Vect().CosTheta()};
  /* PC + Lambda */
  double plam_mass[2]  = { (Llam+Lpc).M()                           , (cLlam+cLpc).M()};
  double plam_mom[2]   = { (Llam+Lpc).Vect().Mag()                  , (cLlam+cLpc).Vect().Mag()};
  double plam_cost[2]  = { (Llam+Lpc).Vect().CosTheta()             , (cLlam+cLpc).Vect().CosTheta()};
  double plam_phi[2]   = { (Llam+Lpc).Vect().Phi()                  , (cLlam+cLpc).Vect().Phi()};
  double plam_mmass[2] = { (Ltgt+Lbeam-(Llam+Lpc)).M()              , (cLtgt+cLbeam-(cLlam+cLpc)).M()};
  double plam_mmom[2]  = { (Ltgt+Lbeam-(Llam+Lpc)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam+cLpc)).Vect().Mag()};
  double plam_mcost[2] = { (Ltgt+Lbeam-(Llam+Lpc)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam+cLpc)).Vect().CosTheta()};

  double tp = cLpc.E() - cLpc.M();
  double tn = (cLtgt+cLbeam-cLlam-cLpc).E() - (cLtgt+cLbeam-cLlam-cLpc).M();
  double tl = cLlam.E() - cLlam.M();
  double qvalue = tn+tp+tl;

  /* flag */
  bool mnflag = false; /* Missing-n flag */
  bool sblmnflag = false; /* lower sideband of Missing-n flag */
  bool sbumnflag = false; /* upper sideband of Missing-n flag */
  {
    double pro_mmass = plam_mmass[1];
    if(mnll<pro_mmass&&pro_mmass<mnul) mnflag = true;
    if(sblmnll<pro_mmass&&pro_mmass<sblmnul) sblmnflag = true;
    if(sbumnll<pro_mmass&&pro_mmass<sbumnul) sbumnflag = true;
  }

  if(!mnflag) {return true;}
  //if(flam1flag) {return true;}
  //if(flam2flag) {return true;}
  //if(!flam1flag&&!flam2flag) {return true;}

  // Proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Mass_%s_%s",frame[f].Data(),mode.Data()),p_mass[f]);
    FillHist(Form("p_Momentum_%s_%s",frame[f].Data(),mode.Data()),p_mom[f]);
    FillHist(Form("p_CosT_%s_%s",frame[f].Data(),mode.Data()),p_cost[f]);
    FillHist(Form("p_Phi_%s_%s",frame[f].Data(),mode.Data()),p_phi[f]);
    FillHist(Form("p_MMass_%s_%s",frame[f].Data(),mode.Data()),p_mmass[f]);
    FillHist(Form("p_MMomentum_%s_%s",frame[f].Data(),mode.Data()),p_mmom[f]);
    FillHist(Form("p_MCosT_%s_%s",frame[f].Data(),mode.Data()),p_mcost[f]);
  }
  // Pi-
  for(int f=1; f<2; f++){
    FillHist(Form("pim_Mass_%s_%s",frame[f].Data(),mode.Data()),pim_mass[f]);
    FillHist(Form("pim_Momentum_%s_%s",frame[f].Data(),mode.Data()),pim_mom[f]);
    FillHist(Form("pim_CosT_%s_%s",frame[f].Data(),mode.Data()),pim_cost[f]);
    FillHist(Form("pim_Phi_%s_%s",frame[f].Data(),mode.Data()),pim_phi[f]);
    FillHist(Form("pim_MMass_%s_%s",frame[f].Data(),mode.Data()),pim_mmass[f]);
    FillHist(Form("pim_MMomentum_%s_%s",frame[f].Data(),mode.Data()),pim_mmom[f]);
    FillHist(Form("pim_MCosT_%s_%s",frame[f].Data(),mode.Data()),pim_mcost[f]);
  }
  // Lambda
  for(int f=1; f<2; f++){
    FillHist(Form("lam_Mass_%s_%s",frame[f].Data(),mode.Data()),lam_mass[f]);
    FillHist(Form("lam_Momentum_%s_%s",frame[f].Data(),mode.Data()),lam_mom[f]);
    FillHist(Form("lam_CosT_%s_%s",frame[f].Data(),mode.Data()),lam_cost[f]);
    FillHist(Form("lam_Phi_%s_%s",frame[f].Data(),mode.Data()),lam_phi[f]);
    FillHist(Form("lam_MMass_%s_%s",frame[f].Data(),mode.Data()),lam_mmass[f]);
    FillHist(Form("lam_MMomentum_%s_%s",frame[f].Data(),mode.Data()),lam_mmom[f]);
    FillHist(Form("lam_MCosT_%s_%s",frame[f].Data(),mode.Data()),lam_mcost[f]);
  }
  // pfwd
  for(int f=1; f<2; f++){
    FillHist(Form("pfwd_Mass_%s_%s",frame[f].Data(),mode.Data()),pc_mass[f]);
    FillHist(Form("pfwd_Momentum_%s_%s",frame[f].Data(),mode.Data()),pc_mom[f]);
    FillHist(Form("pfwd_CosT_%s_%s",frame[f].Data(),mode.Data()),pc_cost[f]);
    FillHist(Form("pfwd_Phi_%s_%s",frame[f].Data(),mode.Data()),pc_phi[f]);
    FillHist(Form("pfwd_MMass_%s_%s",frame[f].Data(),mode.Data()),pc_mmass[f]);
    FillHist(Form("pfwd_MMomentum_%s_%s",frame[f].Data(),mode.Data()),pc_mmom[f]);
    FillHist(Form("pfwd_MCosT_%s_%s",frame[f].Data(),mode.Data()),pc_mcost[f]);
  }
  // p + Proton
  for(int f=1; f<2; f++){
    FillHist(Form("pp_Mass_%s_%s",frame[f].Data(),mode.Data()),pp_mass[f]);
    FillHist(Form("pp_Momentum_%s_%s",frame[f].Data(),mode.Data()),pp_mom[f]);
    FillHist(Form("pp_CosT_%s_%s",frame[f].Data(),mode.Data()),pp_cost[f]);
    FillHist(Form("pp_Phi_%s_%s",frame[f].Data(),mode.Data()),pp_phi[f]);
    FillHist(Form("pp_MMass_%s_%s",frame[f].Data(),mode.Data()),pp_mmass[f]);
    FillHist(Form("pp_MMomentum_%s_%s",frame[f].Data(),mode.Data()),pp_mmom[f]);
    FillHist(Form("pp_MCosT_%s_%s",frame[f].Data(),mode.Data()),pp_mcost[f]);
  }
  // p + pim
  for(int f=1; f<2; f++){
    FillHist(Form("ppim_Mass_%s_%s",frame[f].Data(),mode.Data()),ppim_mass[f]);
    FillHist(Form("ppim_Momentum_%s_%s",frame[f].Data(),mode.Data()),ppim_mom[f]);
    FillHist(Form("ppim_CosT_%s_%s",frame[f].Data(),mode.Data()),ppim_cost[f]);
    FillHist(Form("ppim_Phi_%s_%s",frame[f].Data(),mode.Data()),ppim_phi[f]);
    FillHist(Form("ppim_MMass_%s_%s",frame[f].Data(),mode.Data()),ppim_mmass[f]);
    FillHist(Form("ppim_MMomentum_%s_%s",frame[f].Data(),mode.Data()),ppim_mmom[f]);
    FillHist(Form("ppim_MCosT_%s_%s",frame[f].Data(),mode.Data()),ppim_mcost[f]);
  }
  // p + Lambda
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%s",frame[f].Data(),mode.Data()),plam_mass[f]);
    FillHist(Form("plam_Momentum_%s_%s",frame[f].Data(),mode.Data()),plam_mom[f]);
    FillHist(Form("plam_CosT_%s_%s",frame[f].Data(),mode.Data()),plam_cost[f]);
    FillHist(Form("plam_Phi_%s_%s",frame[f].Data(),mode.Data()),plam_phi[f]);
    FillHist(Form("plam_MMass_%s_%s",frame[f].Data(),mode.Data()),plam_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%s",frame[f].Data(),mode.Data()),plam_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%s",frame[f].Data(),mode.Data()),plam_mcost[f]);
  }

  /*-----------------------------------------*/
  /* 2D Plots                                */ 
  /*-----------------------------------------*/
  // IM(plam) vs. MM(plam)
  for(int f=1; f<2; f++){
    FillHist(Form("plam_IM_vs_plam_MM_%s_%s",frame[f].Data(),mode.Data()),plam_mass[f],plam_mmass[f]);
  }
  // MM(p) vs. MM(plam)
  for(int f=1; f<2; f++){
    FillHist(Form("pfwd_MM_vs_plam_MM_%s_%s",frame[f].Data(),mode.Data()),pc_mmass[f],plam_mmass[f]);
  }
  // MM(p) vs. IM(plam)
  for(int f=1; f<2; f++){
    FillHist(Form("pfwd_MM_vs_plam_IM_%s_%s",frame[f].Data(),mode.Data()),pc_mmass[f],plam_mass[f]);
  }
  // MM(p) vs. CosT(p)
  for(int f=1; f<2; f++){
    FillHist(Form("pfwd_MM_vs_pfwd_CosT_%s_%s",frame[f].Data(),mode.Data()),pc_mmass[f],pc_cost[f]);
  }
  // MM(p) vs. MMom(p)
  for(int f=1; f<2; f++){
    FillHist(Form("pfwd_MM_vs_pfwd_MMom_%s_%s",frame[f].Data(),mode.Data()),pc_mmass[f],pc_mmom[f]);
  }
  // MM(p) vs. MM(pp)
  for(int f=1; f<2; f++){
    FillHist(Form("pfwd_MM_vs_pp_MM_%s_%s",frame[f].Data(),mode.Data()),pc_mmass[f],pp_mmass[f]);
  }
  // Dalitz Plot
  for(int f=1; f<2; f++){
    FillHist(Form("DalitzPlot_%s_%s",frame[f].Data(),mode.Data()),(tn-tl)/sqrt(3)/qvalue,tp/qvalue);
  }

  return true;

}

bool MyAnalysisHeKLpnFWD::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKLpnFWD::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKLpnFWD::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKLpnFWD::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKLpnFWD::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.49732;
  double sigma = 0.00607;
  k0ll = mean-3*sigma; k0ul = mean+3*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Missing neutron mass cut */ 
  mean  = 0.9432;
  sigma = 0.01432;
  //mnll = 0.85; mnul = 1.03;
  mnll = 1.10; mnul = 1.30;
  //mnll = mean-2*sigma; mnul = mean+2*sigma;
  sblmnll = mean-5*sigma; sblmnul = mean-3*sigma;
  sbumnll = mean+3*sigma; sbumnul = mean+5*sigma;
  /* Missing protron mass cut */ 
  mean  = 0.94909;
  sigma = 0.01359;
  //mpll = 0.85; mpul = 1.03;
  mpll = 1.10; mpul = 1.30;
  //mpll = 0.99; mpul = 1.07;
  //mpll = mean-2*sigma; mpul = mean+2*sigma;
  sblmpll = mean-5*sigma; sblmpul = mean-3*sigma;
  sbumpll = mean+3*sigma; sbumpul = mean+5*sigma;
  /* SigmaPlus mass cut */
  mean  = 1.18894;
  sigma = 0.00429;
  spll = mean-3*sigma; spul = mean+3*sigma;
  sblspll = mean-5*sigma; sblspul = mean-3*sigma;
  sbuspll = mean+3*sigma; sbuspul = mean+5*sigma;
  /* SigmaMinus mass cut */
  mean  = 1.19736;
  sigma = 0.00353;
  smll = mean-3*sigma; smul = mean+3*sigma;
  sblsmll = mean-5*sigma; sblsmul = mean-3*sigma;
  sbusmll = mean+3*sigma; sbusmul = mean+5*sigma;
  /* Lambda mass cut */ 
  mean  = 1.1155;
  sigma = 0.0020;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-5*sigma; sbllamul = mean-3*sigma;
  sbulamll = mean+3*sigma; sbulamul = mean+5*sigma;
  /* Forward Lambda mass cut */ 
  mean  = 1.11599;
  sigma = 0.00260;
  flamll = mean-2*sigma; flamul = mean+2*sigma;
  sblflamll = mean-5*sigma; sblflamul = mean-3*sigma;
  sbuflamll = mean+3*sigma; sbuflamul = mean+5*sigma;
}

bool MyAnalysisHeKLpnFWD::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKLpnFWD::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKLpnFWD");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  /* COMMON */
  std::cout << "Define Histograms for Common analysis" << std::endl;
  new TH1F( "EventNumber", "Number of Events", 20, 0, 20);

  new TH1F( "lam_Mass", "Mass of #Lambda;M(#pi^{-}p) (GeV/c^{2});Coutns", 1000, 1.0, 2.0);
  new TH1F( "lam_PDF", "PDF of #Lambda;PDF;Coutns", 1000, 0, 40);
  new TH2F( "lam_MassvsPDF", "Mass  vs. PDF of #Lambda;M(#pi^{-}p) (GeV/c^{2});PDF", 1000, 1.0, 2.0,1000,0,40);
  new TH1F( "lam_Mass_CutPDF", "Mass of #Lambda;M(#pi^{-}p) (GeV/c^{2});Coutns", 1000, 1.0, 2.0);
  new TH1F( "lam_PDF_CutPDF", "PDF of #Lambda;PDF;Coutns", 1000, 0, 40);
  new TH2F( "lam_MassvsPDF_CutPDF", "Mass  vs. PDF of #Lambda;M(#pi^{-}p) (GeV/c^{2});PDF", 1000, 1.0, 2.0,1000,0,40);
  std::cout << "=== End of [Common::Initialize] === " << std::endl;

  /* FWD NEUTRAL */
  std::cout << "Define Histograms for FWD Neutral particles" << std::endl;
  // NC Analysis
  new TH1F( "FWDN_Overbeta", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH1F( "FWDN_Overbeta_withdECut", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH1F( "FWDN_Overbeta_Selected", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH2F( "FWDN_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 500, -0, 10, 500, 0, 100);
  new TH2F( "FWDN_OverbetavsEnergy_withdECut", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 500, -0, 10, 500, 0, 100);
  new TH2F( "NCHitPosition_XY", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",16,-160,160,15,-75,75);
  new TH2F( "NCHitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  // proton
  new TH1F("p_Mass_CM_Neutral","Invariant mass of p;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_Momentum_CM_Neutral","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM_Neutral","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM_Neutral","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM_Neutral", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM_Neutral", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM_Neutral", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-
  new TH1F("pim_Mass_CM_Neutral","Invariant mass of #pi^{-};IM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_Momentum_CM_Neutral","Momentum of #pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_CosT_CM_Neutral","cos#theta of #pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim_Phi_CM_Neutral","#phi of #pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim_MMass_CM_Neutral", "^{3}He(K^{-},#pi^{-})X missing mass;MM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_MMomentum_CM_Neutral", "^{3}He(K^{-},#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_MCosT_CM_Neutral", "^{3}He(K^{-},#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Lambda
  new TH1F("lam_Mass_CM_Neutral","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_Momentum_CM_Neutral","Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_CosT_CM_Neutral","cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_Phi_CM_Neutral","#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("lam_MMass_CM_Neutral", "^{3}He(K^{-},#Lambda)X missing mass;MM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_MMomentum_CM_Neutral", "^{3}He(K^{-},#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_MCosT_CM_Neutral", "^{3}He(K^{-},#Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // nfwd
  new TH1F("nfwd_Mass_CM_Neutral","Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nfwd_Momentum_CM_Neutral","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nfwd_CosT_CM_Neutral","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("nfwd_Phi_CM_Neutral","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("nfwd_MMass_CM_Neutral", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nfwd_MMomentum_CM_Neutral", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nfwd_MCosT_CM_Neutral", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + p
  new TH1F("np_Mass_CM_Neutral","Invariant mass of np;IM(np) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("np_Momentum_CM_Neutral","Momentum of np;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("np_CosT_CM_Neutral","cos#theta of np;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("np_Phi_CM_Neutral","#phi of np;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("np_MMass_CM_Neutral", "^{3}He(K^{-},np)X missing mass;MM(np) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("np_MMomentum_CM_Neutral", "^{3}He(K^{-},np)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("np_MCosT_CM_Neutral", "^{3}He(K^{-},np)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + pim
  new TH1F("npim_Mass_CM_Neutral","Invariant mass of n#pi^{-};IM(n#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npim_Momentum_CM_Neutral","Momentum of n#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npim_CosT_CM_Neutral","cos#theta of n#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npim_Phi_CM_Neutral","#phi of n#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("npim_MMass_CM_Neutral", "^{3}He(K^{-},n#pi^{-})X missing mass;MM(n#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npim_MMomentum_CM_Neutral", "^{3}He(K^{-},n#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npim_MCosT_CM_Neutral", "^{3}He(K^{-},n#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + Lambda
  new TH1F("nlam_Mass_CM_Neutral","Invariant mass of n#Lambda;IM(n#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nlam_Momentum_CM_Neutral","Momentum of n#Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nlam_CosT_CM_Neutral","cos#theta of n#Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("nlam_Phi_CM_Neutral","#phi of n#Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("nlam_MMass_CM_Neutral", "^{3}He(K^{-},n#Lambda)X missing mass;MM(n#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nlam_MMomentum_CM_Neutral", "^{3}He(K^{-},n#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nlam_MCosT_CM_Neutral", "^{3}He(K^{-},n#Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // 2D Plots
  new TH2F("nlam_IM_vs_nlam_MM_CM_Neutral", "n#Lambda invariant mass vs. ^{3}He(K^{-},n#Lambda)X missing mass;IM(n#Lambda) (GeV/c^{2});MM(n#Lambda) (GeV/c^{2})", 200, 1.5, 3.5, 200, 0.5, 2.5);
  new TH2F("nfwd_MM_vs_nlam_MM_CM_Neutral", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},n#Lambda)X missing mass;MM(n) (GeV/c^{2});MM(n#Lambda) (GeV/c^{2})", 200, 1.5, 3.5, 200, 0.5, 2.5);
  new TH2F("nfwd_MM_vs_nlam_IM_CM_Neutral", "^{3}He(K^{-},n)X missing mass vs. n#Lambda invariant mass;MM(n) (GeV/c^{2});IM(n#Lambda) (GeV/c^{2})", 200, 1.5, 3.5, 200, 1.5, 3.5);
  new TH2F("nfwd_MM_vs_nfwd_CosT_CM_Neutral", "^{3}He(K^{-},n)X missing mass vs. n cos#theta;MM(n) (GeV/c^{2});cos#theta_{n}", 200, 1.5, 3.5, 200, 0.8, 1.0);
  new TH2F("nfwd_MM_vs_nfwd_MMom_CM_Neutral", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},n)X momentum;MM(n) (GeV/c^{2});Momentum(X) (GeV/c)", 200, 1.5, 3.5, 200, 0.0, 2.0);
  new TH2F("nfwd_MM_vs_np_MM_CM_Neutral", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},np)X missing mass;MM(n) (GeV/c^{2});MM(np) (GeV/c^{2})", 200, 1.5, 3.5, 200, 0.5, 2.5);
  new TH2F("DalitzPlot_CM_Neutral", "Dalitz Plot;(T_{p}-T_{#Lambda})/#sqrt{3}Q;T_{n}/Q", 1000, -1.0, 1.0, 1000, -1.0, 1.0);

  std::cout << "=== End of [FWDNeutral::Initialize] === " << std::endl;

  /* FWD CHARGED */
  std::cout << "Define Histograms for FWD Charged particles" << std::endl;
  // FDC Tracking
  new TH1F( "FDC1_nTrack", "Number of Tracks in FDC1;Num. of Tracks;Counts", 20, 0, 20 );
  new TH1F( "FDC1_Chi2", "Chi-squared of FDC1;Chi-square;Counts", 500, 0, 50);
  new TH2F( "FDC1_Position", "Position at FDC1;X (cm);Y (cm)", 1000,-50,50, 1000,-50,50);
  // PC Analysis
  new TH1F( "FWDC_Beta", "#beta of FWD Charged;#beta;Counts", 1000,0.0,1.0);
  new TH1F( "FWDC_Overbeta", "1/#beta of FWD Charged;1/#beta;Counts", 5000,0.0,5.0);
  new TH2F( "FWDC_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 5, 500, 0, 50);
  new TH1F( "FWDC_Momentum", "Momentum of FWD Charged;Momentum (GeV/c);Counts", 2000,0.0,2.0);
  new TH1F( "FWDC_MomentumR", "Momentum of FWD Charged;Momentum (GeV/c);Counts", 2000,0.0,2.0);
  new TH2F( "FWDC_MomentumvsMomentumR", "Momentum vs Mometuntum_{R} of FWD Charged;Momentum (GeV/c);Momentum_{R} (GeV/c)", 400,0.0,2.0, 400,0.0,2.0);
  new TH1F( "FWDC_Mass2", "Mass2 of FWD Charged;Mass^{2} (GeV^{2}/c^{4});Counts", 600,-1.0,5.0);
  new TH2F( "FWDC_Mass2vsMomentum", "Mass2 vs Momentum of FWD Charged;Mass^{2} (GeV^{2}/c^{4});Momentum", 600,-1.0,5.0, 200, 0.0, 2.0);
  // proton
  new TH1F("p_Mass_CM_Charged","Invariant mass of p;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_Momentum_CM_Charged","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM_Charged","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM_Charged","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM_Charged", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM_Charged", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM_Charged", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-
  new TH1F("pim_Mass_CM_Charged","Invariant mass of #pi^{-};IM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_Momentum_CM_Charged","Momentum of #pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_CosT_CM_Charged","cos#theta of #pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim_Phi_CM_Charged","#phi of #pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim_MMass_CM_Charged", "^{3}He(K^{-},#pi^{-})X missing mass;MM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_MMomentum_CM_Charged", "^{3}He(K^{-},#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_MCosT_CM_Charged", "^{3}He(K^{-},#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Lambda
  new TH1F("lam_Mass_CM_Charged","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_Momentum_CM_Charged","Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_CosT_CM_Charged","cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_Phi_CM_Charged","#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("lam_MMass_CM_Charged", "^{3}He(K^{-},#Lambda)X missing mass;MM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_MMomentum_CM_Charged", "^{3}He(K^{-},#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_MCosT_CM_Charged", "^{3}He(K^{-},#Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pfwd
  new TH1F("pfwd_Mass_CM_Charged","Invariant mass of p;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pfwd_Momentum_CM_Charged","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pfwd_CosT_CM_Charged","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pfwd_Phi_CM_Charged","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pfwd_MMass_CM_Charged", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pfwd_MMomentum_CM_Charged", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pfwd_MCosT_CM_Charged", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + p
  new TH1F("pp_Mass_CM_Charged","Invariant mass of pp;IM(pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pp_Momentum_CM_Charged","Momentum of pp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pp_CosT_CM_Charged","cos#theta of pp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pp_Phi_CM_Charged","#phi of pp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pp_MMass_CM_Charged", "^{3}He(K^{-},pp)X missing mass;MM(pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pp_MMomentum_CM_Charged", "^{3}He(K^{-},pp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pp_MCosT_CM_Charged", "^{3}He(K^{-},pp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + pim
  new TH1F("ppim_Mass_CM_Charged","Invariant mass of p#pi^{-};IM(p#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppim_Momentum_CM_Charged","Momentum of p#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppim_CosT_CM_Charged","cos#theta of p#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("ppim_Phi_CM_Charged","#phi of p#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("ppim_MMass_CM_Charged", "^{3}He(K^{-},p#pi^{-})X missing mass;MM(p#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppim_MMomentum_CM_Charged", "^{3}He(K^{-},p#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppim_MCosT_CM_Charged", "^{3}He(K^{-},p#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + Lambda
  new TH1F("plam_Mass_CM_Charged","Invariant mass of p#Lambda;IM(p#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("plam_Momentum_CM_Charged","Momentum of p#Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("plam_CosT_CM_Charged","cos#theta of p#Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("plam_Phi_CM_Charged","#phi of p#Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("plam_MMass_CM_Charged", "^{3}He(K^{-},p#Lambda)X missing mass;MM(p#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("plam_MMomentum_CM_Charged", "^{3}He(K^{-},p#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("plam_MCosT_CM_Charged", "^{3}He(K^{-},p#Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // 2D Plots
  new TH2F("plam_IM_vs_plam_MM_CM_Charged", "p#Lambda invariant mass vs. ^{3}He(K^{-},p#Lambda)X missing mass;IM(p#Lambda) (GeV/c^{2});MM(p#Lambda) (GeV/c^{2})", 200, 1.5, 3.5, 200, 0.5, 2.5);
  new TH2F("pfwd_MM_vs_plam_MM_CM_Charged", "^{3}He(K^{-},p)X missing mass vs. ^{3}He(K^{-},p#Lambda)X missing mass;MM(p) (GeV/c^{2});MM(p#Lambda) (GeV/c^{2})", 200, 1.5, 3.5, 200, 0.5, 2.5);
  new TH2F("pfwd_MM_vs_plam_IM_CM_Charged", "^{3}He(K^{-},p)X missing mass vs. n#Lambda invariant mass;MM(p) (GeV/c^{2});IM(p#Lambda) (GeV/c^{2})", 200, 1.5, 3.5, 200, 1.5, 3.5);
  new TH2F("pfwd_MM_vs_pfwd_CosT_CM_Charged", "^{3}He(K^{-},p)X missing mass vs. ,p cos#theta;MM(p) (GeV/c^{2});cos#theta_{p}", 200, 1.5, 3.5, 200, 0.8, 1.0);
  new TH2F("pfwd_MM_vs_pfwd_MMom_CM_Charged", "^{3}He(K^{-},p)X missing mass vs. ^{3}He(K^{-},p)X momentum;MM(p) (GeV/c^{2});Momentum(X) (GeV/c)", 200, 1.5, 3.5, 200, 0.0, 2.0);
  new TH2F("pfwd_MM_vs_pp_MM_CM_Charged", "^{3}He(K^{-},p)X missing mass vs. ^{3}He(K^{-},pp)X missing mass;MM(p) (GeV/c^{2});MM(pp) (GeV/c^{2})", 200, 1.5, 3.5, 200, 0.5, 2.5);
  new TH2F("DalitzPlot_CM_Charged", "Dalitz Plot;(T_{p}-T_{#Lambda})/#sqrt{3}Q;T_{n}/Q", 1000, -1.0, 1.0, 1000, -1.0, 1.0);
  std::cout << "=== End of [FWDCharged::Initialize] === " << std::endl;

  return true;
}
