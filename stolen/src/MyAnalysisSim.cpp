// MyAnalysisSim.cpp

#include "MyAnalysisSim.h"

MyAnalysisSim::MyAnalysisSim(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisSim::~MyAnalysisSim()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisSim::Clear()
{
}

bool MyAnalysisSim::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  int ievent=0;
  rtFile->cd();


  return true;

}

bool MyAnalysisSim::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisSim::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisSim::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisSim::FillHist(TString name, TString val1, TString val2, int weight)
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


bool MyAnalysisSim::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisSim::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKpipp");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

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
  new TH2F("pipp_IMvsHeKpipp_MM_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2})", 200, 2.0, 4.0, 200, 0.0, 2.0);
  new TH2F("pipp_IMvsHeKpipp_MCosT_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. cos#theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta_{n}", 200, 2.0, 4.0, 300, -1.5, 1.5);
  new TH2F("pipp_IMvsMomentumTransfer_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. Momentum transfer;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Momentum Transfer (GeV/c)", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pipp_DalitzPlot_n_CM", "Dalitz Plot;(T_{p}-T_{#Lambda})/#sqrt{3}Q;T_{n}/Q", 1000, -1.0, 1.0, 1000, -1.0, 1.0);

  std::cout << "== Finish Histogram Initialization ==" << std::endl;
  return true;

}
