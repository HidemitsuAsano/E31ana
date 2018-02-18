// MyAnalysisHeKpippMCKinFit.cpp

#include "MyAnalysisHeKpippMCKinFit.h"

MyAnalysisHeKpippMCKinFit::MyAnalysisHeKpippMCKinFit(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKpippMCKinFit::~MyAnalysisHeKpippMCKinFit()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKpippMCKinFit::Clear()
{
}

bool MyAnalysisHeKpippMCKinFit::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle, MCData* mcdata, ReactionData* reaction, DetectorData* detector)
{
  if(reaction==0) return true;

  DetectorList *dlist=DetectorList::GetInstance();
  rtFile->cd();

  int lamID = -1;
  for(int i=0; i<mcdata->trackSize(); i++){
    Track* track = mcdata->track(i);
    if(track->pdgID()==3122&&track->parentTrackID()==0){
      lamID = track->trackID();
    }
  }
  int lamseg=-1;
  int notlamseg=-1;
  for(int i=0; i<detector->detectorHitSize(); i++){
    DetectorHit* hit = detector->detectorHit(i);
    int pdg = hit->pdg();
    int parentID = hit->parentID();
    int detectorID = hit->detectorID();
    int layerID = hit->layerID();
    int channelID = hit->channelID();
    int seg = channelID+1;
    if(detectorID!=CID_CDH) continue;
    if(parentID==lamID&&pdg==2212){
      lamseg = seg;
    }
    if(parentID==0&&pdg==2212){
      notlamseg = seg;
    }
  }

  int size_cm = reaction->CMparticleSize();
  int size_lab = reaction->ParticleSize();
  int size_init = reaction->InitParticleSize();

	TVector3 tmpZeroV;
	TLorentzVector Ltgt_mc; Ltgt_mc.SetVectM(tmpZeroV,ThreeHeMass);
  TLorentzVector Lbeam_mc;
  TLorentzVector Llambda_mc;
  TLorentzVector Lproton_mc;
  TLorentzVector Lnproton_mc;
  TLorentzVector Lneutron_mc;
  TLorentzVector Lpim_mc;
  TLorentzVector Lplam_mc;


	int lambda_ID = 0;
	for(int i=0; i<mcdata->trackSize(); i++){
		if(mcdata->track(i)->parentTrackID()==0&&mcdata->track(i)->pdgID()==3122){
			lambda_ID = mcdata->track(i)->trackID();
		}
	}
	for(int i=0; i<mcdata->trackSize(); i++){
		Track* track = mcdata->track(i);
		if(track->parentTrackID()==0){
			if(track->pdgID()==1000020030){
				Ltgt_mc.SetXYZM(track->momentum().X()*0.001,track->momentum().Y()*0.001,track->momentum().Z()*0.001,ThreeHeMass);
			}
		}
		if(track->parentTrackID()==0){
			if(track->pdgID()==321){
				Lbeam_mc.SetXYZM(-track->momentum().X()*0.001,-track->momentum().Y()*0.001,-track->momentum().Z()*0.001,kpMass);
			}
		}
		if(track->parentTrackID()==0){
			if(track->pdgID()==3122||track->pdgID()==3212){
				double tmp_mass = lMass;
				if(track->pdgID()==3212){ tmp_mass = s0Mass; }
				Llambda_mc.SetXYZM(track->momentum().X()*0.001,track->momentum().Y()*0.001,track->momentum().Z()*0.001,tmp_mass);
			}
		}
		if(track->parentTrackID()==0){
			if(track->pdgID()==2212){
				Lproton_mc.SetXYZM(track->momentum().X()*0.001,track->momentum().Y()*0.001,track->momentum().Z()*0.001,pMass);
			}
		}
		if(track->parentTrackID()==0){
			if(track->pdgID()==2112){
				Lneutron_mc.SetXYZM(track->momentum().X()*0.001,track->momentum().Y()*0.001,track->momentum().Z()*0.001,nMass);
			}
		}
		if(track->parentTrackID()==lambda_ID){
			if(track->pdgID()==2212){
				Lnproton_mc.SetXYZM(track->momentum().X()*0.001,track->momentum().Y()*0.001,track->momentum().Z()*0.001,pMass);
			}
		}
		if(track->parentTrackID()==lambda_ID){
			if(track->pdgID()==-211){
				Lpim_mc.SetXYZM(track->momentum().X()*0.001,track->momentum().Y()*0.001,track->momentum().Z()*0.001,piMass);
			}
		}
	}
	Lplam_mc = Llambda_mc + Lproton_mc;

	TVector3 boost_mc = (Ltgt_mc+Lbeam_mc).BoostVector();
	TLorentzVector cLtgt_mc = Ltgt_mc; cLtgt_mc.Boost(-boost_mc);
	TLorentzVector cLbeam_mc = Lbeam_mc; cLbeam_mc.Boost(-boost_mc);
	TLorentzVector cLlambda_mc = Llambda_mc; cLlambda_mc.Boost(-boost_mc);
	TLorentzVector cLproton_mc = Lproton_mc; cLproton_mc.Boost(-boost_mc);
	TLorentzVector cLneutron_mc = Lneutron_mc; cLneutron_mc.Boost(-boost_mc);
	TLorentzVector cLnproton_mc = Lnproton_mc; cLnproton_mc.Boost(-boost_mc);
	TLorentzVector cLpim_mc = Lpim_mc; cLpim_mc.Boost(-boost_mc);
	TLorentzVector cLplam_mc = Lplam_mc; cLplam_mc.Boost(-boost_mc);

	TString frame[2] = {"Lab","CM"};

	double lam_mc_mass[2]  = { Llambda_mc.M()                           , cLlambda_mc.M()};
	double lam_mc_mom[2]   = { Llambda_mc.Vect().Mag()                  , cLlambda_mc.Vect().Mag()};
	double lam_mc_cost[2]  = { Llambda_mc.Vect().CosTheta()             , cLlambda_mc.Vect().CosTheta()};
	double lam_mc_phi[2]   = { Llambda_mc.Vect().Phi()                  , cLlambda_mc.Vect().Phi()};
	double lam_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Llambda_mc).M()              , (cLtgt_mc+cLbeam_mc-cLlambda_mc).M()};
	double lam_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Llambda_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLlambda_mc).Vect().Mag()};
	double lam_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Llambda_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLlambda_mc).Vect().CosTheta()};

	double p1_mc_mass[2]  = { Lproton_mc.M()                           , cLproton_mc.M()};
	double p1_mc_mom[2]   = { Lproton_mc.Vect().Mag()                  , cLproton_mc.Vect().Mag()};
	double p1_mc_cost[2]  = { Lproton_mc.Vect().CosTheta()             , cLproton_mc.Vect().CosTheta()};
	double p1_mc_phi[2]   = { Lproton_mc.Vect().Phi()                  , cLproton_mc.Vect().Phi()};
	double p1_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Lproton_mc).M()              , (cLtgt_mc+cLbeam_mc-cLproton_mc).M()};
	double p1_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Lproton_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLproton_mc).Vect().Mag()};
	double p1_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Lproton_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLproton_mc).Vect().CosTheta()};

	double pim_mc_mass[2]  = { Lpim_mc.M()                           , cLpim_mc.M()};
	double pim_mc_mom[2]   = { Lpim_mc.Vect().Mag()                  , cLpim_mc.Vect().Mag()};
	double pim_mc_cost[2]  = { Lpim_mc.Vect().CosTheta()             , cLpim_mc.Vect().CosTheta()};
	double pim_mc_phi[2]   = { Lpim_mc.Vect().Phi()                  , cLpim_mc.Vect().Phi()};
	double pim_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Lpim_mc).M()              , (cLtgt_mc+cLbeam_mc-cLpim_mc).M()};
	double pim_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Lpim_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLpim_mc).Vect().Mag()};
	double pim_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Lpim_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLpim_mc).Vect().CosTheta()};

	double p2_mc_mass[2]  = { Lnproton_mc.M()                           , cLnproton_mc.M()};
	double p2_mc_mom[2]   = { Lnproton_mc.Vect().Mag()                  , cLnproton_mc.Vect().Mag()};
	double p2_mc_cost[2]  = { Lnproton_mc.Vect().CosTheta()             , cLnproton_mc.Vect().CosTheta()};
	double p2_mc_phi[2]   = { Lnproton_mc.Vect().Phi()                  , cLnproton_mc.Vect().Phi()};
	double p2_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Lnproton_mc).M()              , (cLtgt_mc+cLbeam_mc-cLnproton_mc).M()};
	double p2_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Lnproton_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLnproton_mc).Vect().Mag()};
	double p2_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Lnproton_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLnproton_mc).Vect().CosTheta()};

	double n_mc_mass[2]  = { Lneutron_mc.M()                           , cLneutron_mc.M()};
	double n_mc_mom[2]   = { Lneutron_mc.Vect().Mag()                  , cLneutron_mc.Vect().Mag()};
	double n_mc_cost[2]  = { Lneutron_mc.Vect().CosTheta()             , cLneutron_mc.Vect().CosTheta()};
	double n_mc_phi[2]   = { Lneutron_mc.Vect().Phi()                  , cLneutron_mc.Vect().Phi()};
	double n_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Lneutron_mc).M()              , (cLtgt_mc+cLbeam_mc-cLneutron_mc).M()};
	double n_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Lneutron_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLneutron_mc).Vect().Mag()};
	double n_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Lneutron_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLneutron_mc).Vect().CosTheta()};

	double plam_mc_mass[2]  = { Lplam_mc.M()                           , cLplam_mc.M()};
	double plam_mc_mom[2]   = { Lplam_mc.Vect().Mag()                  , cLplam_mc.Vect().Mag()};
	double plam_mc_cost[2]  = { Lplam_mc.Vect().CosTheta()             , cLplam_mc.Vect().CosTheta()};
	double plam_mc_phi[2]   = { Lplam_mc.Vect().Phi()                  , cLplam_mc.Vect().Phi()};
	double plam_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Lplam_mc).M()              , (cLtgt_mc+cLbeam_mc-cLplam_mc).M()};
	double plam_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Lplam_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLplam_mc).Vect().Mag()};
	double plam_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Lplam_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLplam_mc).Vect().CosTheta()};

	double q_mc[2] = {(Lbeam_mc.Vect()-Lneutron_mc.Vect()).Mag(), (Lbeam_mc.Vect()-Lneutron_mc.Vect()).Mag()};

	if(particle->nBeam()!=1) return false;
	pBeam* beam = particle->beam(0);
	beam->SetSIM(true);
	double mom = beam->mom();
	double mom2 = rndm->Gaus(0.0,0.002);
	beam->SetMomentum(mom+mom2);
	// ########## //
	// BLDC Chi^2 //
	// ########## //
	double blc1chi = beam->blc1chi2();
	double blc2chi = beam->blc2chi2();
	double bpcchi = beam->bpcchi2();
	double beamchi = beam->beamchi2();

	if(blc1chi>10.0){ return true; }
	if(blc2chi>10.0){ return true; }
	if(bpcchi>10.0){ return true; }
	if(beamchi>20.0){ return true; }

	int MulCDH=0;
	for(int i=0; i<cdsMan->nCDH(); i++){
		HodoscopeLikeHit* hit = cdsMan->CDH(i);
		if(hit->CheckRange()) MulCDH++;
	}
	if(MulCDH!=3){ return true; }

	if(cdstrackMan->nGoodTrack()!=3) return true;

	/* Event selection */
	int MulBVC=0;
	for(int i=0; i<blMan->nBVC(); i++){
		HodoscopeLikeHit* hit = blMan->BVC(i);
		if(hit->CheckRange()) MulBVC++;
	}
	int MulCVC=0;
	for(int i=0; i<blMan->nCVC(); i++){
		HodoscopeLikeHit* hit = blMan->CVC(i);
		if(hit->CheckRange()) MulCVC++;
	}
	int MulPC=0;
	for(int i=0; i<blMan->nPC(); i++){
		HodoscopeLikeHit* hit = blMan->PC(i);
		if(hit->CheckRange()) MulPC++;
	}
	if(MulBVC!=0||MulCVC!=0||MulPC!=0) return true;

	int MulIH=0;
	for(int i=0; i<cdsMan->nIH(); i++){
		HodoscopeLikeHit* hit = cdsMan->IH(i);
		if(hit->CheckRange()) MulIH++;
	}
	//if(MulIH<=3){ return true; }

	if(particle->nCDS()!=3) return true;
	if(particle->nProton()!=2) return true;
	if(particle->nPiminus()!=1) return true;

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
	pCDS* p1 = particle->proton(1);
	pCDS* p2 = particle->proton(0);

	if(pim->chi()>30) return true;
	if(p1->chi()>30) return true;
	if(p2->chi()>30) return true;

	// ######################### //
	// 2 Helix calculation       //
	// ######################### //
	pCDS* pip1 = 0;
	pCDS* pip2 = 0;
	{
		for(int it=0; it<particle->nProduct(); it++){
			pCDS* product = particle->product(it);
			int comb = product->comb();
			if(comb==TMath::Power(2,CDS_Proton)+TMath::Power(2,CDS_PiMinus)){
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

	// Target
	TLorentzVector tmpLtgt; tmpLtgt.SetVectM(tmpZeroV,ThreeHeMass);
	TLorentzVector tmpLbeam = beam->GetLorentzVector(vertex);

	TLorentzVector tmpLpim = pim->GetLorentzVector();
	TLorentzVector tmpLp1 = p1->GetLorentzVector();
	TLorentzVector tmpLp2 = p2->GetLorentzVector();
	double tmp_im = (tmpLpim+tmpLp1+tmpLp2).M();
	double tmp_mm = (tmpLtgt+tmpLbeam-tmpLpim-tmpLp1-tmpLp2).M();

	bool pip1flag = false, pip2flag = false;
	pCDS* lam = 0;
	pCDS* nlam = 0;
	pCDS* proton = 0;
	pCDS* nproton = 0;
	TVector3 v_pip1 = pip1->vbeam();
	TVector3 v_p2 = p2->vbeam();
	double pip1p2dis = (v_pip1-v_p2).Mag();
	double tmpper1[5] = {pip1->mass(),pip1->vdis(),pip1->pbdca(),p2->vbdis(),pip1p2dis};
	double tmppdf1[5] = {0,0,0,0,0};
	TrackTools::PDFLambda(tmpper1,tmppdf1);
	TVector3 v_pip2 = pip2->vbeam();
	TVector3 v_p1 = p1->vbeam();
	double pip2p1dis = (v_pip2-v_p1).Mag();
	double tmpper2[5] = {pip2->mass(),pip2->vdis(),pip2->pbdca(),p1->vbdis(),pip2p1dis};
	double tmppdf2[5] = {0,0,0,0,0};
	/* Input : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
	/*Output : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
	TrackTools::PDFLambda(tmpper2,tmppdf2);
	/* All */
	double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]*tmppdf1[2]*tmppdf1[3]*tmppdf1[4]);
	double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]*tmppdf2[2]*tmppdf2[3]*tmppdf2[4]);

	double pdf=-1;
	double npdf=-1;

	double pdfll =  0.0;
	double pdful = 15.75;

	if(pdf1<pdf2){
		if(pdfll<pdf1&&pdf1<pdful){
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
		if(pdfll<pdf2&&pdf2<pdful){
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
	if(!pip1flag&&!pip2flag){ return true; }
	if(lam==0){ return true; }

	// ################ //
	// Vertex decision1 //
	// ################ //
	const int fiducial = CID_Fiducial;
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

	// ################ //
	// Vertex decision2 //
	// ################ //
	TVector3 vertex_lam = lam->vbeam();
	if(GeomTools::GetID(vertex_lam)!=fiducial){ return true; }

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
	TVector3 vtxb;

	// Target
	TVector3 ZeroV;
	TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
	TLorentzVector Lbeam = beam->GetLorentzVector(vertex);

	// Pi+, Pi-, and Pi+ Pi- pair //
	TLorentzVector Lpim = pim->GetLorentzVector();
	TLorentzVector Llam = lam->GetLorentzVector();
	TLorentzVector Lp = proton->GetLorentzVector();
	TLorentzVector Lnp = nproton->GetLorentzVector();
	TLorentzVector Lpipp = Llam + Lp;
	TLorentzVector Lmn = Ltgt+Lbeam-Lpipp;

	TVector3 boost = (Ltgt+Lbeam).BoostVector();
	TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

	TLorentzVector cLpim = Lpim; cLpim.Boost(-boost);
	TLorentzVector cLlam = Llam; cLlam.Boost(-boost);
	TLorentzVector cLp = Lp; cLp.Boost(-boost);
	TLorentzVector cLnp = Lnp; cLnp.Boost(-boost);
	TLorentzVector cLpipp = cLlam + cLp;
	TLorentzVector cLmn = Lmn; cLmn.Boost(-boost);
	TLorentzVector Lpipn = Llam + Lmn, cLpipn = cLlam + cLmn;

	TLorentzVector Lp1 = Lp, cLp1 = cLp;
	TLorentzVector Lp2 = Lnp, cLp2 = cLnp;
	TLorentzVector Lpip1 = Lp1 + Lpim, Lpip2 = Lp2 + Lpim;
	TLorentzVector cLpip1 = cLp1 + cLpim, cLpip2 = cLp2 + cLpim;
	TLorentzVector Lnlam = Lp1 + cLpim, cLnlam  = cLp1 + cLpim;
	TLorentzVector Lnpipp = Lnlam + Lp2, cLnpipp = cLnlam + cLp2;

	TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

	double pim_mass[2]  = { Lpim.M()                           , cLpim.M()};
	double pim_mom[2]   = { Lpim.Vect().Mag()                  , cLpim.Vect().Mag()};
	double pim_cost[2]  = { Lpim.Vect().Dot(Lbeam.Vect())/Lpim.Vect().Mag()/Lbeam.Vect().Mag(),cLpim.Vect().Dot(cLbeam.Vect())/cLpim.Vect().Mag()/cLbeam.Vect().Mag()};
	double pim_phi[2]   = { Lpim.Vect().Phi()                  , cLpim.Vect().Phi()};
	double pim_mmass[2] = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
	double pim_mmom[2]  = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
	double pim_mcost[2] = { (Ltgt+Lbeam-Lpim).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lpim).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLpim).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLpim).Vect().Mag()/cLbeam.Vect().Mag()};

	double p_mass[2]  = { Lp.M()                           , cLp.M()};
	double p_mom[2]   = { Lp.Vect().Mag()                  , cLp.Vect().Mag()};
	double p_cost[2]  = { Lp.Vect().Dot(Lbeam.Vect())/Lp.Vect().Mag()/Lbeam.Vect().Mag(),cLp.Vect().Dot(cLbeam.Vect())/cLp.Vect().Mag()/cLbeam.Vect().Mag()};
	double p_phi[2]   = { Lp.Vect().Phi()                  , cLp.Vect().Phi()};
	double p_mmass[2] = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
	double p_mmom[2]  = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
	double p_mcost[2] = { (Ltgt+Lbeam-Lp).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lp).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLp).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLp).Vect().Mag()/cLbeam.Vect().Mag()};

	double p1_mass[2]  = { Lp1.M()                           , cLp1.M()};
	double p1_mom[2]   = { Lp1.Vect().Mag()                  , cLp1.Vect().Mag()};
	double p1_cost[2]  = { Lp1.Vect().Dot(Lbeam.Vect())/Lp1.Vect().Mag()/Lbeam.Vect().Mag(),cLp1.Vect().Dot(cLbeam.Vect())/cLp1.Vect().Mag()/cLbeam.Vect().Mag()};
	double p1_phi[2]   = { Lp1.Vect().Phi()                  , cLp1.Vect().Phi()};
	double p1_mmass[2] = { (Ltgt+Lbeam-Lp1).M()              , (cLtgt+cLbeam-cLp1).M()};
	double p1_mmom[2]  = { (Ltgt+Lbeam-Lp1).Vect().Mag()     , (cLtgt+cLbeam-cLp1).Vect().Mag()};
	double p1_mcost[2] = { (Ltgt+Lbeam-Lp1).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lp1).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLp1).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLp1).Vect().Mag()/cLbeam.Vect().Mag()};

	double p2_mass[2]  = { Lp2.M()                           , cLp2.M()};
	double p2_mom[2]   = { Lp2.Vect().Mag()                  , cLp2.Vect().Mag()};
	double p2_cost[2]  = { Lp2.Vect().Dot(Lbeam.Vect())/Lp2.Vect().Mag()/Lbeam.Vect().Mag(),cLp2.Vect().Dot(cLbeam.Vect())/cLp2.Vect().Mag()/cLbeam.Vect().Mag()};
	double p2_phi[2]   = { Lp2.Vect().Phi()                  , cLp2.Vect().Phi()};
	double p2_mmass[2] = { (Ltgt+Lbeam-Lp2).M()              , (cLtgt+cLbeam-cLp2).M()};
	double p2_mmom[2]  = { (Ltgt+Lbeam-Lp2).Vect().Mag()     , (cLtgt+cLbeam-cLp2).Vect().Mag()};
	double p2_mcost[2] = { (Ltgt+Lbeam-Lp2).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lp2).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLp2).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLp2).Vect().Mag()/cLbeam.Vect().Mag()};

	double mn_mass[2]  = { Lmn.M()                           , cLmn.M()};
	double mn_mom[2]   = { Lmn.Vect().Mag()                  , cLmn.Vect().Mag()};
	double mn_cost[2]  = { Lmn.Vect().Dot(Lbeam.Vect())/Lmn.Vect().Mag()/Lbeam.Vect().Mag(),cLmn.Vect().Dot(cLbeam.Vect())/cLmn.Vect().Mag()/cLbeam.Vect().Mag()};
	double mn_phi[2]   = { Lmn.Vect().Phi()                  , cLmn.Vect().Phi()};
	double mn_mmass[2] = { (Ltgt+Lbeam-Lmn).M()              , (cLtgt+cLbeam-cLmn).M()};
	double mn_mmom[2]  = { (Ltgt+Lbeam-Lmn).Vect().Mag()     , (cLtgt+cLbeam-cLmn).Vect().Mag()};
	double mn_mcost[2] = { (Ltgt+Lbeam-Lmn).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lmn).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLmn).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLmn).Vect().Mag()/cLbeam.Vect().Mag()};

	double pip1_mass[2]  = { Lpip1.M()                           , cLpip1.M()};
	double pip1_mom[2]   = { Lpip1.Vect().Mag()                  , cLpip1.Vect().Mag()};
	double pip1_cost[2]  = { Lpip1.Vect().Dot(Lbeam.Vect())/Lpip1.Vect().Mag()/Lbeam.Vect().Mag(),cLpip1.Vect().Dot(cLbeam.Vect())/cLpip1.Vect().Mag()/cLbeam.Vect().Mag()};
	double pip1_phi[2]   = { Lpip1.Vect().Phi()                  , cLpip1.Vect().Phi()};
	double pip1_mmass[2] = { (Ltgt+Lbeam-Lpip1).M()              , (cLtgt+cLbeam-cLpip1).M()};
	double pip1_mmom[2]  = { (Ltgt+Lbeam-Lpip1).Vect().Mag()     , (cLtgt+cLbeam-cLpip1).Vect().Mag()};
	double pip1_mcost[2] = { (Ltgt+Lbeam-Lpip1).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lpip1).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLpip1).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLpip1).Vect().Mag()/cLbeam.Vect().Mag()};

	double pip2_mass[2]  = { Lpip2.M()                           , cLpip2.M()};
	double pip2_mom[2]   = { Lpip2.Vect().Mag()                  , cLpip2.Vect().Mag()};
	double pip2_cost[2]  = { Lpip2.Vect().Dot(Lbeam.Vect())/Lpip2.Vect().Mag()/Lbeam.Vect().Mag(),cLpip2.Vect().Dot(cLbeam.Vect())/cLpip2.Vect().Mag()/cLbeam.Vect().Mag()};
	double pip2_phi[2]   = { Lpip2.Vect().Phi()                  , cLpip2.Vect().Phi()};
	double pip2_mmass[2] = { (Ltgt+Lbeam-Lpip2).M()              , (cLtgt+cLbeam-cLpip2).M()};
	double pip2_mmom[2]  = { (Ltgt+Lbeam-Lpip2).Vect().Mag()     , (cLtgt+cLbeam-cLpip2).Vect().Mag()};
	double pip2_mcost[2] = { (Ltgt+Lbeam-Lpip2).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lpip2).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLpip2).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLpip2).Vect().Mag()/cLbeam.Vect().Mag()};

	double lam_mass[2]  = { Llam.M()                           , cLlam.M()};
	double lam_mom[2]   = { Llam.Vect().Mag()                  , cLlam.Vect().Mag()};
	double lam_cost[2]  = { Llam.Vect().Dot(Lbeam.Vect())/Llam.Vect().Mag()/Lbeam.Vect().Mag(),cLlam.Vect().Dot(cLbeam.Vect())/cLlam.Vect().Mag()/cLbeam.Vect().Mag()};
	double lam_phi[2]   = { Llam.Vect().Phi()                  , cLlam.Vect().Phi()};
	double lam_mmass[2] = { (Ltgt+Lbeam-Llam).M()              , (cLtgt+cLbeam-cLlam).M()};
	double lam_mmom[2]  = { (Ltgt+Lbeam-Llam).Vect().Mag()     , (cLtgt+cLbeam-cLlam).Vect().Mag()};
	double lam_mcost[2] = { (Ltgt+Lbeam-Llam).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Llam).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLlam).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLlam).Vect().Mag()/cLbeam.Vect().Mag()};

	double nlam_mass[2]  = { Lnlam.M()                           , cLnlam.M()};
	double nlam_mom[2]   = { Lnlam.Vect().Mag()                  , cLnlam.Vect().Mag()};
	double nlam_cost[2]  = { Lnlam.Vect().Dot(Lbeam.Vect())/Lnlam.Vect().Mag()/Lbeam.Vect().Mag(),cLnlam.Vect().Dot(cLbeam.Vect())/cLnlam.Vect().Mag()/cLbeam.Vect().Mag()};
	double nlam_phi[2]   = { Lnlam.Vect().Phi()                  , cLnlam.Vect().Phi()};
	double nlam_mmass[2] = { (Ltgt+Lbeam-Lnlam).M()              , (cLtgt+cLbeam-cLnlam).M()};
	double nlam_mmom[2]  = { (Ltgt+Lbeam-Lnlam).Vect().Mag()     , (cLtgt+cLbeam-cLnlam).Vect().Mag()};
	double nlam_mcost[2] = { (Ltgt+Lbeam-Lnlam).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lnlam).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLnlam).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLnlam).Vect().Mag()/cLbeam.Vect().Mag()};

	double pipp_mass[2]  = { Lpipp.M()                           , cLpipp.M()};
	double pipp_mom[2]   = { Lpipp.Vect().Mag()                  , cLpipp.Vect().Mag()};
	double pipp_cost[2]  = { Lpipp.Vect().Dot(Lbeam.Vect())/Lpipp.Vect().Mag()/Lbeam.Vect().Mag(),cLpipp.Vect().Dot(cLbeam.Vect())/cLpipp.Vect().Mag()/cLbeam.Vect().Mag()};
	double pipp_phi[2]   = { Lpipp.Vect().Phi()                  , cLpipp.Vect().Phi()};
	double pipp_mmass[2] = { (Ltgt+Lbeam-Lpipp).M()              , (cLtgt+cLbeam-cLpipp).M()};
	double pipp_mmom[2]  = { (Ltgt+Lbeam-Lpipp).Vect().Mag()     , (cLtgt+cLbeam-cLpipp).Vect().Mag()};
	double pipp_mcost[2] = { (Ltgt+Lbeam-Lpipp).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lpipp).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLpipp).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLpipp).Vect().Mag()/cLbeam.Vect().Mag()};

	double npipp_mass[2]  = { Lnpipp.M()                           , cLnpipp.M()};
	double npipp_mom[2]   = { Lnpipp.Vect().Mag()                  , cLnpipp.Vect().Mag()};
	double npipp_cost[2]  = { Lnpipp.Vect().Dot(Lbeam.Vect())/Lnpipp.Vect().Mag()/Lbeam.Vect().Mag(),cLnpipp.Vect().Dot(cLbeam.Vect())/cLnpipp.Vect().Mag()/cLbeam.Vect().Mag()};
	double npipp_phi[2]   = { Lnpipp.Vect().Phi()                  , cLnpipp.Vect().Phi()};
	double npipp_mmass[2] = { (Ltgt+Lbeam-Lnpipp).M()              , (cLtgt+cLbeam-cLnpipp).M()};
	double npipp_mmom[2]  = { (Ltgt+Lbeam-Lnpipp).Vect().Mag()     , (cLtgt+cLbeam-cLnpipp).Vect().Mag()};
	double npipp_mcost[2] = { (Ltgt+Lbeam-Lnpipp).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lnpipp).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLnpipp).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLnpipp).Vect().Mag()/cLbeam.Vect().Mag()
	};
	double pipn_mass[2]  = { Lpipn.M()                           , cLpipn.M()};
	double pipn_mom[2]   = { Lpipn.Vect().Mag()                  , cLpipn.Vect().Mag()};
	double pipn_cost[2]  = { Lpipn.Vect().Dot(Lbeam.Vect())/Lpipn.Vect().Mag()/Lbeam.Vect().Mag(),cLpipn.Vect().Dot(cLbeam.Vect())/cLpipn.Vect().Mag()/cLbeam.Vect().Mag()};
	double pipn_phi[2]   = { Lpipn.Vect().Phi()                  , cLpipn.Vect().Phi()};
	double pipn_mmass[2] = { (Ltgt+Lbeam-Lpipn).M()              , (cLtgt+cLbeam-cLpipn).M()};
	double pipn_mmom[2]  = { (Ltgt+Lbeam-Lpipn).Vect().Mag()     , (cLtgt+cLbeam-cLpipn).Vect().Mag()};
	double pipn_mcost[2] = { (Ltgt+Lbeam-Lpipn).Vect().Dot(Lbeam.Vect())/(Ltgt+Lbeam-Lpipn).Vect().Mag()/Lbeam.Vect().Mag(),(cLtgt+cLbeam-cLpipn).Vect().Dot(cLbeam.Vect())/(cLtgt+cLbeam-cLpipn).Vect().Mag()/cLbeam.Vect().Mag()};

	double q[2] = {(Lbeam.Vect()-Lmn.Vect()).Mag(), (Lbeam.Vect()-Lmn.Vect()).Mag()};
	double t[2] = {(Lbeam-Lmn).M2(), (cLbeam-cLmn).M2()};
	double tn = cLmn.E() - cLmn.M();
	double tp = cLp.E() - cLp.M();
	double tl = cLlam.E() - cLlam.M();
	double qvalue = tn+tp+tl;

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
	/* Missing n cut */
	//if(!mnflag){ return true; }

	double covariance[6][4] = {
		//{0.00453701,0.0045759,0.00073352,0.000648903},
		{0.00453701,0.0045759,0.00200460,0.001796010},
		{0.0209434,0.0209419,0.0161643,0.0204749},
		{0.0327885,0.0327452,0.0238122,0.0329356},
		{0.0252103,0.0252616,0.0181587,0.0201182},
		{0.0222299,.0222292,0.0150418,0.0201182},
		{0.00967645,0.00964559,0.00546052,0.00535962}
	};

	//-----------------------------------------//
	//--- covariance matrices for KinFitter ---//
	//-----------------------------------------//
	double covVal[6][16] = {
		{ covariance[0][0]*covariance[0][0], 0.0, 0.0, 0.0,
			0.0, covariance[0][1]*covariance[0][1], 0.0, 0.0,
			0.0, 0.0, covariance[0][2]*covariance[0][2], 0.0,
			0.0, 0.0, 0.0, covariance[0][3]*covariance[0][3]},

		{ covariance[1][0]*covariance[1][0], 0.0, 0.0, 0.0,
			0.0, covariance[1][1]*covariance[1][1], 0.0, 0.0,
			0.0, 0.0, covariance[1][2]*covariance[1][2], 0.0,
			0.0, 0.0, 0.0, covariance[1][3]*covariance[1][3]},

		{ covariance[2][0]*covariance[2][0], 0.0, 0.0, 0.0,
			0.0, covariance[2][1]*covariance[2][1], 0.0, 0.0,
			0.0, 0.0, covariance[2][2]*covariance[2][2], 0.0,
			0.0, 0.0, 0.0, covariance[2][3]*covariance[2][3]},

		{ covariance[3][0]*covariance[3][0], 0.0, 0.0, 0.0,
			0.0, covariance[3][1]*covariance[3][1], 0.0, 0.0,
			0.0, 0.0, covariance[3][2]*covariance[3][2], 0.0,
			0.0, 0.0, 0.0, covariance[3][3]*covariance[3][3]},

		{ covariance[4][0]*covariance[4][0], 0.0, 0.0, 0.0,
			0.0, covariance[4][1]*covariance[4][1], 0.0, 0.0,
			0.0, 0.0, covariance[4][2]*covariance[4][2], 0.0,
			0.0, 0.0, 0.0, covariance[4][3]*covariance[4][3]},

		{ covariance[5][0]*covariance[5][0], 0.0, 0.0, 0.0,
			0.0, covariance[5][1]*covariance[5][1], 0.0, 0.0,
			0.0, 0.0, covariance[5][2]*covariance[5][2], 0.0,
			0.0, 0.0, 0.0, covariance[5][3]*covariance[5][3]}
	};
	TMatrixD *covZero = new TMatrixD(4, 4);
	covZero->Zero();
	covZero->ResizeTo(3, 3); // resize from 4x4 to 3x3
	TMatrixD *covParticle[6];
	for( int i=0; i<6; i++ ){
		covParticle[i] = new TMatrixD(4, 4);
		int n = 0;
		for( int j=0; j<4; j++ ){
			for( int k=0; k<4; k++ ){
				if( j==k ) (*covParticle[i])[j][k] = covVal[i][n]; // only diagonal elements
				else       (*covParticle[i])[j][k] = 0;
				n++;
			}
		}
		covParticle[i]->ResizeTo(3, 3); // resize from 4x4 to 3x3
		//covParticle[i]->Print(); // Print all
	}
	//-----------------------------------------//
	//--- covariance matrices for KinFitter ---//
	//-----------------------------------------//

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	// %%% Kinematical Fit using KinFitter %%% //
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	//--- set TLorentzVector ---//
	TLorentzVector TL_meas[6]; // measured
	TLorentzVector TL_kfit[6]; // kinematical fitted
	TL_meas[0] = Lbeam;
	TL_meas[1] = Llam; //L_Llab;
	TL_meas[2] = Lmn; //L_nlab;
	TL_meas[3] = Lp;
	TL_meas[4] = Lnp;
	TL_meas[5] = Lpim;
	// L3_target is defined as (0, 0, 0, M_3He)
	TVector3 TV_target = Ltgt.Vect();
	TVector3 TV_meas[6];
	for( int i=0; i<6; i++ ){
		TV_meas[i] = TL_meas[i].Vect();
	}

	TDatabasePDG *pdg = new TDatabasePDG();
	pdg->ReadPDGTable("pdg_table.txt");
	int PDG[6] = {321, 3122, 2112, 2212, 2212, -211};

	//--- KinFitter :: initialization ---//
	//*** definition of fit particles in cartesian coordinates ***//
	TString str_particle[6] = {"Lbeam", "Llam", "Lmn", "Lp", "Lnp", "Lpim"};
	TFitParticlePxPyPz ParticleTgt = TFitParticlePxPyPz("target", "target", &TV_target,
			pdg->GetParticle("He3")->Mass(), covZero);
	TFitParticlePxPyPz Particle[6];
	for( int i=0; i<6; i++ ){
		Particle[i] = TFitParticlePxPyPz(str_particle[i], str_particle[i], &TV_meas[i],
				pdg->GetParticle(PDG[i])->Mass(), covParticle[i]);
	}
	//*** definition of constraints ***//
	// constraint :: mass of Lambda
	TFitConstraintM ConstML = TFitConstraintM("M_L", "M_L", 0, 0, pdg->GetParticle(PDG[1])->Mass());
	ConstML.addParticles1(&Particle[4], &Particle[5]);
	// constraint :: 4-momentum conservation
	TFitConstraintEp ConstEp[4];
	TString str_constEp[4]  = {"Px", "Py", "Pz", "E"};
	for( int i=0; i<4; i++ ){
		ConstEp[i] = TFitConstraintEp(str_constEp[i], str_constEp[i], 0, TFitConstraintEp::component(i), 0);
		ConstEp[i].addParticles1(&ParticleTgt, &Particle[0]);
		ConstEp[i].addParticles2(&Particle[2], &Particle[3], &Particle[4], &Particle[5]);
	}

	//--- KinFitter :: execution ---//
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//*** definition of the fitter ***//
	TKinFitter kinfitter;
	kinfitter.addMeasParticles(&Particle[0], &Particle[3], &Particle[4], &Particle[5]);
	kinfitter.addUnmeasParticles(&Particle[2]);
	kinfitter.addConstraint(&ConstML);
	for( int i=0; i<4; i++ ){
		kinfitter.addConstraint(&ConstEp[i]);
	}
	//*** perform the fit ***//
	kinfitter.setMaxNbIter(50);       // max number of iterations
	kinfitter.setMaxDeltaS(5e-5);     // max delta chi2
	kinfitter.setMaxF(1e-4);          // max sum of constraints
	//kinfitter.setVerbosity(KFDEBUG);  // verbosity level
	kinfitter.fit();

	double chi2 = kinfitter.getS();
	int ndf = kinfitter.getNDF();
	double chi2r = chi2/(double)ndf;
	FillHist(Form("KinFitter_Chi2"),chi2);
	FillHist(Form("KinFitter_NDF"),ndf);
	FillHist(Form("KinFitter_Chi2r"),chi2r);

	//*** copy fit results ***//
	for( int i=0; i<6; i++ ){
		TL_kfit[i] = (*Particle[i].getCurr4Vec());
	}
	TL_kfit[1] = TL_kfit[4]+TL_kfit[5];
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	// %%% Kinematical Fit using KinFitter %%% //
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //


	TLorentzVector TL_gene[6]; /* generated */
	TL_gene[0] = Lbeam_mc;
	TL_gene[1] = Llambda_mc;
	TL_gene[2] = Lneutron_mc;
	TL_gene[3] = Lproton_mc;
	TL_gene[4] = Lnproton_mc;
	TL_gene[5] = Lpim_mc;

	for(int i=0; i<6; i++){
		for(int j=0; j<4; j++){
			double cov = TL_meas[i][j]-TL_gene[i][j];
			FillHist(Form("cov_%d_%d",i,j),cov);
		}
	}
	for(int i=0; i<6; i++){
		for(int j=0; j<4; j++){
			double cov = TL_kfit[i][j]-TL_gene[i][j];
			FillHist(Form("res_%d_%d",i,j),cov);
		}
	}

	// Pi+, Pi-, and Pi+ Pi- pair //
	TLorentzVector Lbeam_kfit = TL_kfit[0];
	TLorentzVector Llam_kfit = TL_kfit[1];
	TLorentzVector Lmn_kfit = TL_kfit[2];
	TLorentzVector Lp_kfit = TL_kfit[3];
	TLorentzVector Lnp_kfit = TL_kfit[4];
	TLorentzVector Lpim_kfit = TL_kfit[5];
	TLorentzVector Lpipp_kfit = Llam_kfit + Lp_kfit;
	TLorentzVector Lnlam_kfit = Lpim_kfit + Lp_kfit;
	TLorentzVector Ltgt_kfit = Ltgt;

	TVector3 boost_kfit = (Ltgt+Lbeam_kfit).BoostVector();
	TLorentzVector cLtgt_kfit = Ltgt; cLtgt_kfit.Boost(-boost_kfit);
	TLorentzVector cLbeam_kfit = Lbeam_kfit; cLbeam_kfit.Boost(-boost_kfit);
	TLorentzVector cLlam_kfit = Llam_kfit; cLlam_kfit.Boost(-boost_kfit);
	TLorentzVector cLmn_kfit = Lmn_kfit; cLmn_kfit.Boost(-boost_kfit);
	TLorentzVector cLp_kfit = Lp_kfit; cLp_kfit.Boost(-boost_kfit);
	TLorentzVector cLnp_kfit = Lnp_kfit; cLnp_kfit.Boost(-boost_kfit);
	TLorentzVector cLpim_kfit = Lpim_kfit; cLpim_kfit.Boost(-boost_kfit);
	TLorentzVector cLpipp_kfit = Lpipp_kfit; cLpipp_kfit.Boost(-boost_kfit);
	TLorentzVector cLnlam_kfit = Lnlam_kfit; cLnlam_kfit.Boost(-boost_kfit);

	double pim_kfit_mass[2]  = { Lpim_kfit.M()                           , cLpim_kfit.M()};
	double pim_kfit_mom[2]   = { Lpim_kfit.Vect().Mag()                  , cLpim_kfit.Vect().Mag()};
	double pim_kfit_cost[2]  = { Lpim_kfit.Vect().Dot(Lbeam_kfit.Vect())/Lpim_kfit.Vect().Mag()/Lbeam_kfit.Vect().Mag(),cLpim_kfit.Vect().Dot(cLbeam_kfit.Vect())/cLpim_kfit.Vect().Mag()/cLbeam_kfit.Vect().Mag()};
	double pim_kfit_phi[2]   = { Lpim_kfit.Vect().Phi()                  , cLpim_kfit.Vect().Phi()};
	double pim_kfit_mmass[2] = { (Ltgt_kfit+Lbeam_kfit-Lpim_kfit).M()              , (cLtgt_kfit+cLbeam_kfit-cLpim_kfit).M()};
	double pim_kfit_mmom[2]  = { (Ltgt_kfit+Lbeam_kfit-Lpim_kfit).Vect().Mag()     , (cLtgt_kfit+cLbeam_kfit-cLpim_kfit).Vect().Mag()};
	double pim_kfit_mcost[2] = { (Ltgt_kfit+Lbeam_kfit-Lpim_kfit).Vect().Dot(Lbeam_kfit.Vect())/(Ltgt_kfit+Lbeam_kfit-Lpim_kfit).Vect().Mag()/Lbeam_kfit.Vect().Mag(),(cLtgt_kfit+cLbeam_kfit-cLpim_kfit).Vect().Dot(cLbeam_kfit.Vect())/(cLtgt_kfit+cLbeam_kfit-cLpim_kfit).Vect().Mag()/cLbeam_kfit.Vect().Mag()};

	double p1_kfit_mass[2]  = { Lp_kfit.M()                           , cLp_kfit.M()};
	double p1_kfit_mom[2]   = { Lp_kfit.Vect().Mag()                  , cLp_kfit.Vect().Mag()};
	double p1_kfit_cost[2]  = { Lp_kfit.Vect().Dot(Lbeam_kfit.Vect())/Lp_kfit.Vect().Mag()/Lbeam_kfit.Vect().Mag(),cLp_kfit.Vect().Dot(cLbeam_kfit.Vect())/cLp_kfit.Vect().Mag()/cLbeam_kfit.Vect().Mag()};
	double p1_kfit_phi[2]   = { Lp_kfit.Vect().Phi()                  , cLp_kfit.Vect().Phi()};
	double p1_kfit_mmass[2] = { (Ltgt_kfit+Lbeam_kfit-Lp_kfit).M()              , (cLtgt_kfit+cLbeam_kfit-cLp_kfit).M()};
	double p1_kfit_mmom[2]  = { (Ltgt_kfit+Lbeam_kfit-Lp_kfit).Vect().Mag()     , (cLtgt_kfit+cLbeam_kfit-cLp_kfit).Vect().Mag()};
	double p1_kfit_mcost[2] = { (Ltgt_kfit+Lbeam_kfit-Lp_kfit).Vect().Dot(Lbeam_kfit.Vect())/(Ltgt_kfit+Lbeam_kfit-Lp_kfit).Vect().Mag()/Lbeam_kfit.Vect().Mag(),(cLtgt_kfit+cLbeam_kfit-cLp_kfit).Vect().Dot(cLbeam_kfit.Vect())/(cLtgt_kfit+cLbeam_kfit-cLp_kfit).Vect().Mag()/cLbeam_kfit.Vect().Mag()};

	double p2_kfit_mass[2]  = { Lnp_kfit.M()                           , cLnp_kfit.M()};
	double p2_kfit_mom[2]   = { Lnp_kfit.Vect().Mag()                  , cLnp_kfit.Vect().Mag()};
	double p2_kfit_cost[2]  = { Lnp_kfit.Vect().Dot(Lbeam_kfit.Vect())/Lnp_kfit.Vect().Mag()/Lbeam_kfit.Vect().Mag(),cLnp_kfit.Vect().Dot(cLbeam_kfit.Vect())/cLnp_kfit.Vect().Mag()/cLbeam_kfit.Vect().Mag()};
	double p2_kfit_phi[2]   = { Lnp_kfit.Vect().Phi()                  , cLnp_kfit.Vect().Phi()};
	double p2_kfit_mmass[2] = { (Ltgt_kfit+Lbeam_kfit-Lnp_kfit).M()              , (cLtgt_kfit+cLbeam_kfit-cLnp_kfit).M()};
	double p2_kfit_mmom[2]  = { (Ltgt_kfit+Lbeam_kfit-Lnp_kfit).Vect().Mag()     , (cLtgt_kfit+cLbeam_kfit-cLnp_kfit).Vect().Mag()};
	double p2_kfit_mcost[2] = { (Ltgt_kfit+Lbeam_kfit-Lnp_kfit).Vect().Dot(Lbeam_kfit.Vect())/(Ltgt_kfit+Lbeam_kfit-Lnp_kfit).Vect().Mag()/Lbeam_kfit.Vect().Mag(),(cLtgt_kfit+cLbeam_kfit-cLnp_kfit).Vect().Dot(cLbeam_kfit.Vect())/(cLtgt_kfit+cLbeam_kfit-cLnp_kfit).Vect().Mag()/cLbeam_kfit.Vect().Mag()};

	double mn_kfit_mass[2]  = { Lmn_kfit.M()                           , cLmn_kfit.M()};
	double mn_kfit_mom[2]   = { Lmn_kfit.Vect().Mag()                  , cLmn_kfit.Vect().Mag()};
	double mn_kfit_cost[2]  = { Lmn_kfit.Vect().Dot(Lbeam_kfit.Vect())/Lmn_kfit.Vect().Mag()/Lbeam_kfit.Vect().Mag(),cLmn_kfit.Vect().Dot(cLbeam_kfit.Vect())/cLmn_kfit.Vect().Mag()/cLbeam_kfit.Vect().Mag()};
	double mn_kfit_phi[2]   = { Lmn_kfit.Vect().Phi()                  , cLmn_kfit.Vect().Phi()};
	double mn_kfit_mmass[2] = { (Ltgt_kfit+Lbeam_kfit-Lmn_kfit).M()              , (cLtgt_kfit+cLbeam_kfit-cLmn_kfit).M()};
	double mn_kfit_mmom[2]  = { (Ltgt_kfit+Lbeam_kfit-Lmn_kfit).Vect().Mag()     , (cLtgt_kfit+cLbeam_kfit-cLmn_kfit).Vect().Mag()};
	double mn_kfit_mcost[2] = { (Ltgt_kfit+Lbeam_kfit-Lmn_kfit).Vect().Dot(Lbeam_kfit.Vect())/(Ltgt_kfit+Lbeam_kfit-Lmn_kfit).Vect().Mag()/Lbeam_kfit.Vect().Mag(),(cLtgt_kfit+cLbeam_kfit-cLmn_kfit).Vect().Dot(cLbeam_kfit.Vect())/(cLtgt_kfit+cLbeam_kfit-cLmn_kfit).Vect().Mag()/cLbeam_kfit.Vect().Mag()};

	double lam_kfit_mass[2]  = { Llam_kfit.M()                           , cLlam_kfit.M()};
	double lam_kfit_mom[2]   = { Llam_kfit.Vect().Mag()                  , cLlam_kfit.Vect().Mag()};
	double lam_kfit_cost[2]  = { Llam_kfit.Vect().Dot(Lbeam_kfit.Vect())/Llam_kfit.Vect().Mag()/Lbeam_kfit.Vect().Mag(),cLlam_kfit.Vect().Dot(cLbeam_kfit.Vect())/cLlam_kfit.Vect().Mag()/cLbeam_kfit.Vect().Mag()};
	double lam_kfit_phi[2]   = { Llam_kfit.Vect().Phi()                  , cLlam_kfit.Vect().Phi()};
	double lam_kfit_mmass[2] = { (Ltgt_kfit+Lbeam_kfit-Llam_kfit).M()              , (cLtgt_kfit+cLbeam_kfit-cLlam_kfit).M()};
	double lam_kfit_mmom[2]  = { (Ltgt_kfit+Lbeam_kfit-Llam_kfit).Vect().Mag()     , (cLtgt_kfit+cLbeam_kfit-cLlam_kfit).Vect().Mag()};
	double lam_kfit_mcost[2] = { (Ltgt_kfit+Lbeam_kfit-Llam_kfit).Vect().Dot(Lbeam_kfit.Vect())/(Ltgt_kfit+Lbeam_kfit-Llam_kfit).Vect().Mag()/Lbeam_kfit.Vect().Mag(),(cLtgt_kfit+cLbeam_kfit-cLlam_kfit).Vect().Dot(cLbeam_kfit.Vect())/(cLtgt_kfit+cLbeam_kfit-cLlam_kfit).Vect().Mag()/cLbeam_kfit.Vect().Mag()};

	double nlam_kfit_mass[2]  = { Lnlam_kfit.M()                           , cLnlam_kfit.M()};
	double nlam_kfit_mom[2]   = { Lnlam_kfit.Vect().Mag()                  , cLnlam_kfit.Vect().Mag()};
	double nlam_kfit_cost[2]  = { Lnlam_kfit.Vect().Dot(Lbeam_kfit.Vect())/Lnlam_kfit.Vect().Mag()/Lbeam_kfit.Vect().Mag(),cLnlam_kfit.Vect().Dot(cLbeam_kfit.Vect())/cLnlam_kfit.Vect().Mag()/cLbeam_kfit.Vect().Mag()};
	double nlam_kfit_phi[2]   = { Lnlam_kfit.Vect().Phi()                  , cLnlam_kfit.Vect().Phi()};
	double nlam_kfit_mmass[2] = { (Ltgt_kfit+Lbeam_kfit-Lnlam_kfit).M()              , (cLtgt_kfit+cLbeam_kfit-cLnlam_kfit).M()};
	double nlam_kfit_mmom[2]  = { (Ltgt_kfit+Lbeam_kfit-Lnlam_kfit).Vect().Mag()     , (cLtgt_kfit+cLbeam_kfit-cLnlam_kfit).Vect().Mag()};
	double nlam_kfit_mcost[2] = { (Ltgt_kfit+Lbeam_kfit-Lnlam_kfit).Vect().Dot(Lbeam_kfit.Vect())/(Ltgt_kfit+Lbeam_kfit-Lnlam_kfit).Vect().Mag()/Lbeam_kfit.Vect().Mag(),(cLtgt_kfit+cLbeam_kfit-cLnlam_kfit).Vect().Dot(cLbeam_kfit.Vect())/(cLtgt_kfit+cLbeam_kfit-cLnlam_kfit).Vect().Mag()/cLbeam_kfit.Vect().Mag()};

	double pipp_kfit_mass[2]  = { Lpipp_kfit.M()                           , cLpipp_kfit.M()};
	double pipp_kfit_mom[2]   = { Lpipp_kfit.Vect().Mag()                  , cLpipp_kfit.Vect().Mag()};
	double pipp_kfit_cost[2]  = { Lpipp_kfit.Vect().Dot(Lbeam_kfit.Vect())/Lpipp_kfit.Vect().Mag()/Lbeam_kfit.Vect().Mag(),cLpipp_kfit.Vect().Dot(cLbeam_kfit.Vect())/cLpipp_kfit.Vect().Mag()/cLbeam_kfit.Vect().Mag()};
	double pipp_kfit_phi[2]   = { Lpipp_kfit.Vect().Phi()                  , cLpipp_kfit.Vect().Phi()};
	double pipp_kfit_mmass[2] = { (Ltgt_kfit+Lbeam_kfit-Lpipp_kfit).M()              , (cLtgt_kfit+cLbeam_kfit-cLpipp_kfit).M()};
	double pipp_kfit_mmom[2]  = { (Ltgt_kfit+Lbeam_kfit-Lpipp_kfit).Vect().Mag()     , (cLtgt_kfit+cLbeam_kfit-cLpipp_kfit).Vect().Mag()};
	double pipp_kfit_mcost[2] = { (Ltgt_kfit+Lbeam_kfit-Lpipp_kfit).Vect().Dot(Lbeam_kfit.Vect())/(Ltgt_kfit+Lbeam_kfit-Lpipp_kfit).Vect().Mag()/Lbeam_kfit.Vect().Mag(),(cLtgt_kfit+cLbeam_kfit-cLpipp_kfit).Vect().Dot(cLbeam_kfit.Vect())/(cLtgt_kfit+cLbeam_kfit-cLpipp_kfit).Vect().Mag()/cLbeam_kfit.Vect().Mag()};

	/* lam */
	FillHist("cov_mass_lam",lam_mass[1]-lam_mc_mass[1]);
	FillHist("cov_mass_lam_vs_mass_pipp",lam_mass[1]-lam_mc_mass[1],plam_mc_mass[1]);
	FillHist("res_mass_lam",lam_kfit_mass[1]-lam_mc_mass[1]);
	FillHist("res_mass_lam_vs_mass_pipp",lam_kfit_mass[1]-lam_mc_mass[1],plam_mc_mass[1]);
	FillHist("cov_mom_lam",lam_mom[1]-lam_mc_mom[1]);
	FillHist("cov_mom_lam_vs_mass_pipp",lam_mom[1]-lam_mc_mom[1],plam_mc_mass[1]);
	FillHist("res_mom_lam",lam_kfit_mom[1]-lam_mc_mom[1]);
	FillHist("res_mom_lam_vs_mass_pipp",lam_kfit_mom[1]-lam_mc_mom[1],plam_mc_mass[1]);
	FillHist("cov_cost_lam",lam_cost[1]-lam_mc_cost[1]);
	FillHist("cov_cost_lam_vs_mass_pipp",lam_cost[1]-lam_mc_cost[1],plam_mc_mass[1]);
	FillHist("res_cost_lam",lam_kfit_cost[1]-lam_mc_cost[1]);
	FillHist("res_cost_lam_vs_mass_pipp",lam_kfit_cost[1]-lam_mc_cost[1],plam_mc_mass[1]);
	FillHist("cov_phi_lam",lam_phi[1]-lam_mc_phi[1]);
	FillHist("cov_phi_lam_vs_mass_pipp",lam_phi[1]-lam_mc_phi[1],plam_mc_mass[1]);
	FillHist("res_phi_lam",lam_kfit_phi[1]-lam_mc_phi[1]);
	FillHist("res_phi_lam_vs_mass_pipp",lam_kfit_phi[1]-lam_mc_phi[1],plam_mc_mass[1]);
	/* p */
	FillHist("cov_mass_p1",p1_mass[1]-p1_mc_mass[1]);
	FillHist("cov_mass_p1_vs_mass_pipp",p1_mass[1]-p1_mc_mass[1],plam_mc_mass[1]);
	FillHist("res_mass_p1",p1_kfit_mass[1]-p1_mc_mass[1]);
	FillHist("res_mass_p1_vs_mass_pipp",p1_kfit_mass[1]-p1_mc_mass[1],plam_mc_mass[1]);
	FillHist("cov_mom_p1",p1_mom[1]-p1_mc_mom[1]);
	FillHist("cov_mom_p1_vs_mass_pipp",p1_mom[1]-p1_mc_mom[1],plam_mc_mass[1]);
	FillHist("res_mom_p1",p1_kfit_mom[1]-p1_mc_mom[1]);
	FillHist("res_mom_p1_vs_mass_pipp",p1_kfit_mom[1]-p1_mc_mom[1],plam_mc_mass[1]);
	FillHist("cov_cost_p1",p1_cost[1]-p1_mc_cost[1]);
	FillHist("cov_cost_p1_vs_mass_pipp",p1_cost[1]-p1_mc_cost[1],plam_mc_mass[1]);
	FillHist("res_cost_p1",p1_kfit_cost[1]-p1_mc_cost[1]);
	FillHist("res_cost_p1_vs_mass_pipp",p1_kfit_cost[1]-p1_mc_cost[1],plam_mc_mass[1]);
	FillHist("cov_phi_p1",p1_phi[1]-p1_mc_phi[1]);
	FillHist("cov_phi_p1_vs_mass_pipp",p1_phi[1]-p1_mc_phi[1],plam_mc_mass[1]);
	FillHist("res_phi_p1",p1_kfit_phi[1]-p1_mc_phi[1]);
	FillHist("res_phi_p1_vs_mass_pipp",p1_kfit_phi[1]-p1_mc_phi[1],plam_mc_mass[1]);
	/* pim */
	FillHist("cov_mass_pim",pim_mass[1]-pim_mc_mass[1]);
	FillHist("cov_mass_pim_vs_mass_pipp",pim_mass[1]-pim_mc_mass[1],plam_mc_mass[1]);
	FillHist("res_mass_pim",pim_kfit_mass[1]-pim_mc_mass[1]);
	FillHist("res_mass_pim_vs_mass_pipp",pim_kfit_mass[1]-pim_mc_mass[1],plam_mc_mass[1]);
	FillHist("cov_mom_pim",pim_mom[1]-pim_mc_mom[1]);
	FillHist("cov_mom_pim_vs_mass_pipp",pim_mom[1]-pim_mc_mom[1],plam_mc_mass[1]);
	FillHist("res_mom_pim",pim_kfit_mom[1]-pim_mc_mom[1]);
	FillHist("res_mom_pim_vs_mass_pipp",pim_kfit_mom[1]-pim_mc_mom[1],plam_mc_mass[1]);
	FillHist("cov_cost_pim",pim_cost[1]-pim_mc_cost[1]);
	FillHist("cov_cost_pim_vs_mass_pipp",pim_cost[1]-pim_mc_cost[1],plam_mc_mass[1]);
	FillHist("res_cost_pim",pim_kfit_cost[1]-pim_mc_cost[1]);
	FillHist("res_cost_pim_vs_mass_pipp",pim_kfit_cost[1]-pim_mc_cost[1],plam_mc_mass[1]);
	FillHist("cov_phi_pim",pim_phi[1]-pim_mc_phi[1]);
	FillHist("cov_phi_pim_vs_mass_pipp",pim_phi[1]-pim_mc_phi[1],plam_mc_mass[1]);
	FillHist("res_phi_pim",pim_kfit_phi[1]-pim_mc_phi[1]);
	FillHist("res_phi_pim_vs_mass_pipp",pim_kfit_phi[1]-pim_mc_phi[1],plam_mc_mass[1]);
	/* p from lambda */
	FillHist("cov_mass_p2",p2_mass[1]-p2_mc_mass[1]);
	FillHist("cov_mass_p2_vs_mass_pipp",p2_mass[1]-p2_mc_mass[1],plam_mc_mass[1]);
	FillHist("res_mass_p2",p2_kfit_mass[1]-p2_mc_mass[1]);
	FillHist("res_mass_p2_vs_mass_pipp",p2_kfit_mass[1]-p2_mc_mass[1],plam_mc_mass[1]);
	FillHist("cov_mom_p2",p2_mom[1]-p2_mc_mom[1]);
	FillHist("cov_mom_p2_vs_mass_pipp",p2_mom[1]-p2_mc_mom[1],plam_mc_mass[1]);
	FillHist("res_mom_p2",p2_kfit_mom[1]-p2_mc_mom[1]);
	FillHist("res_mom_p2_vs_mass_pipp",p2_kfit_mom[1]-p2_mc_mom[1],plam_mc_mass[1]);
	FillHist("cov_cost_p2",p2_cost[1]-p2_mc_cost[1]);
	FillHist("cov_cost_p2_vs_mass_pipp",p2_cost[1]-p2_mc_cost[1],plam_mc_mass[1]);
	FillHist("res_cost_p2",p2_kfit_cost[1]-p2_mc_cost[1]);
	FillHist("res_cost_p2_vs_mass_pipp",p2_kfit_cost[1]-p2_mc_cost[1],plam_mc_mass[1]);
	FillHist("cov_phi_p2",p2_phi[1]-p2_mc_phi[1]);
	FillHist("cov_phi_p2_vs_mass_pipp",p2_phi[1]-p2_mc_phi[1],plam_mc_mass[1]);
	FillHist("res_phi_p2",p2_kfit_phi[1]-p2_mc_phi[1]);
	FillHist("res_phi_p2_vs_mass_pipp",p2_kfit_phi[1]-p2_mc_phi[1],plam_mc_mass[1]);
	/* n */
	FillHist("cov_mass_n",mn_mass[1]-n_mc_mass[1]);
	FillHist("cov_mass_n_vs_mass_pipp",mn_mass[1]-n_mc_mass[1],plam_mc_mass[1]);
	FillHist("res_mass_n",mn_kfit_mass[1]-n_mc_mass[1]);
	FillHist("res_mass_n_vs_mass_pipp",mn_kfit_mass[1]-n_mc_mass[1],plam_mc_mass[1]);
	FillHist("cov_mom_n",mn_mom[1]-n_mc_mom[1]);
	FillHist("cov_mom_n_vs_mass_pipp",mn_mom[1]-n_mc_mom[1],plam_mc_mass[1]);
	FillHist("res_mom_n",mn_kfit_mom[1]-n_mc_mom[1]);
	FillHist("res_mom_n_vs_mass_pipp",mn_kfit_mom[1]-n_mc_mom[1],plam_mc_mass[1]);
	FillHist("cov_cost_n",mn_cost[1]-n_mc_cost[1]);
	FillHist("cov_cost_n_vs_mass_pipp",mn_cost[1]-n_mc_cost[1],plam_mc_mass[1]);
	FillHist("res_cost_n",mn_kfit_cost[1]-n_mc_cost[1]);
	FillHist("res_cost_n_vs_mass_pipp",mn_kfit_cost[1]-n_mc_cost[1],plam_mc_mass[1]);
	FillHist("cov_phi_n",mn_phi[1]-n_mc_phi[1]);
	FillHist("cov_phi_n_vs_mass_pipp",mn_phi[1]-n_mc_phi[1],plam_mc_mass[1]);
	FillHist("res_phi_n",mn_kfit_phi[1]-n_mc_phi[1]);
	FillHist("res_phi_n_vs_mass_pipp",mn_kfit_phi[1]-n_mc_phi[1],plam_mc_mass[1]);
	/* Lambda + p */
	FillHist("cov_mass_pipp",pipp_mass[1]-plam_mc_mass[1]);
	FillHist("cov_mass_pipp_vs_mass_pipp",pipp_mass[1]-plam_mc_mass[1],plam_mc_mass[1]);
	FillHist("res_mass_pipp",pipp_kfit_mass[1]-plam_mc_mass[1]);
	FillHist("res_mass_pipp_vs_mass_pipp",pipp_kfit_mass[1]-plam_mc_mass[1],plam_mc_mass[1]);
	FillHist("cov_mom_pipp",pipp_mom[1]-plam_mc_mom[1]);
	FillHist("cov_mom_pipp_vs_mass_pipp",pipp_mom[1]-plam_mc_mom[1],plam_mc_mass[1]);
	FillHist("res_mom_pipp",pipp_kfit_mom[1]-plam_mc_mom[1]);
	FillHist("res_mom_pipp_vs_mass_pipp",pipp_kfit_mom[1]-plam_mc_mom[1],plam_mc_mass[1]);
	FillHist("cov_cost_pipp",pipp_cost[1]-plam_mc_cost[1]);
	FillHist("cov_cost_pipp_vs_mass_pipp",pipp_cost[1]-plam_mc_cost[1],plam_mc_mass[1]);
	FillHist("res_cost_pipp",pipp_kfit_cost[1]-plam_mc_cost[1]);
	FillHist("res_cost_pipp_vs_mass_pipp",pipp_kfit_cost[1]-plam_mc_cost[1],plam_mc_mass[1]);
	FillHist("cov_phi_pipp",pipp_phi[1]-plam_mc_phi[1]);
	FillHist("cov_phi_pipp_vs_mass_pipp",pipp_phi[1]-plam_mc_phi[1],plam_mc_mass[1]);
	FillHist("res_phi_pipp",pipp_kfit_phi[1]-plam_mc_phi[1]);
	FillHist("res_phi_pipp_vs_mass_pipp",pipp_kfit_phi[1]-plam_mc_phi[1],plam_mc_mass[1]);

	FillHist("mmass_vs_mass",pipp_mmass[1],plam_mc_mass[1]);

	FillHist("pipp_mom_gene",plam_mc_mom[1]);
	FillHist("pipp_mom_meas",pipp_mom[1]);
	FillHist("pipp_mom_kfit",pipp_kfit_mom[1]);
	FillHist("pipp_cost_gene",plam_mc_cost[1]);
	FillHist("pipp_cost_meas",pipp_cost[1]);
	FillHist("pipp_cost_kfit",pipp_kfit_cost[1]);
	FillHist("pipp_phi_gene",plam_mc_phi[1]);
	FillHist("pipp_phi_meas",pipp_phi[1]);
	FillHist("pipp_phi_kfit",pipp_kfit_phi[1]);


	//std::cout << "########################################" << std::endl;
	//std::cout << "p_l = (" << cLlambda_mc.Vect().X() << ", " << cLlambda_mc.Vect().Y() << ", " << cLlambda_mc.Vect().Z() << ")" << std::endl;
	//std::cout << "p_p = (" << cLproton_mc.Vect().X() << ", " << cLproton_mc.Vect().Y() << ", " << cLproton_mc.Vect().Z() << ")" << std::endl;
	//std::cout << "p_lp = (" << cLplam_mc.Vect().X() << ", " << cLplam_mc.Vect().Y() << ", " << cLplam_mc.Vect().Z() << ")" << std::endl;
	//std::cout << "########################################" << std::endl;

	return true;

}

bool MyAnalysisHeKpippMCKinFit::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKpippMCKinFit::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKpippMCKinFit::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKpippMCKinFit::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKpippMCKinFit::CutCondition()
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
	//mnll = 0.85; mnul = 0.94;
	//mnll = 0.94; mnul = 1.03;
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

bool MyAnalysisHeKpippMCKinFit::Initialize(ConfMan* confMan)
{
	std::cout << "### MyAnalysisHeKpippMCKinFit::Initialize ###" << std::endl;

	rndm = new TRandom();
	rndm->SetSeed(0);

	std::string ofname = confMan->GetOutFileName();
	ofname.insert(ofname.find(".root"),"_anaHeKpippMCKinFit");

	rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
	rtFile -> cd();

	for(int i=0; i<6; i++){
		for(int j=0; j<4; j++){
			new TH1F( Form("cov_%d_%d",i,j), Form("covariance %d-%d;covariance;Coutns",i,j), 2000, -0.2, 0.2);
			new TH1F( Form("res_%d_%d",i,j), Form("resolution %d-%d;resolution;Coutns",i,j), 2000, -0.2, 0.2);
		}
	}
	TString v_name[4] = {"mass","mom","cost","phi"};
	TString p_name[6] = {"pim","p1","p2","n","lam","pipp"};
	for(int i=0; i<6; i++){
		for(int j=0; j<4; j++){
			new TH1F( Form("cov_%s_%s",v_name[j].Data(),p_name[i].Data()), Form("covariance %s - %s;covariance;Coutns",v_name[j].Data(),p_name[i].Data()), 2000, -0.2, 0.2);
			new TH1F( Form("res_%s_%s",v_name[j].Data(),p_name[i].Data()), Form("resolution %s - %s;resolution;Coutns",v_name[j].Data(),p_name[i].Data()), 2000, -0.2, 0.2);

			new TH2F( Form("cov_%s_%s_vs_mass_pipp",v_name[j].Data(),p_name[i].Data()), Form("covariance %s - %s;resolution;Coutns",v_name[j].Data(),p_name[i].Data()), 400, -0.2, 0.2, 500, 2.0, 3.0);
			new TH2F( Form("res_%s_%s_vs_mass_pipp",v_name[j].Data(),p_name[i].Data()), Form("resolution %s - %s;resolution;Coutns",v_name[j].Data(),p_name[i].Data()), 400, -0.2, 0.2, 500, 2.0, 3.0);
		}
	}
	new TH2F( "mmass_vs_mass","mmass vs. mass",400,0.0,2.0,200,2.0,3.0);
	TString a_name[3] =  {"gene","meas","kfit"};
	for(int i=0; i<3; i++){
		new TH1F( Form("pipp_mom_%s",a_name[i].Data()), "", 2000, 0, 2.0);
		new TH1F( Form("pipp_cost_%s",a_name[i].Data()), "", 2000, -1.0, 1.0);
		new TH1F( Form("pipp_phi_%s",a_name[i].Data()), "", 2000, -3.2, 3.2);
	}

	std::cout << "== Finish Histogram Initialization ==" << std::endl;
	return true;

}
