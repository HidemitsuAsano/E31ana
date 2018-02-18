// MyAnalysisFWDNeutral.cpp

#include "MyAnalysisFWDNeutral.h"

//#define DEBUG=1

MyAnalysisFWDNeutral::MyAnalysisFWDNeutral(TFile* rtFile, ConfMan* conf, bool simflag)
{
	Initialize(rtFile, conf);
	SIMULATION = simflag;
	Clear();
}

MyAnalysisFWDNeutral::~MyAnalysisFWDNeutral()
{
}

void MyAnalysisFWDNeutral::Clear()
{
}

bool MyAnalysisFWDNeutral::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
	DetectorList *dlist=DetectorList::GetInstance();

	if(particle->nBeam()!=1) return false;
	pBeam* beam = particle->beam(0);

	// ################# //
	// Neutral selection //
	// ################# //
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

#ifdef DEBUG
	//std::cout << "### Charge Veto ###" << std::endl;
#endif

	// ########### //
	// NC analysis //
	// ########### //
	for(int i=0; i<blMan->nNC(); i++){
		HodoscopeLikeHit* hit = blMan->NC(i);
		if(hit->CheckRange()){
			pNC* nc = new pNC();
			nc->SetSegment(hit->seg());
			nc->SetHitPosition(hit->pos().X(),hit->hitpos(),hit->pos().Z());
			nc->SetTime(hit->ctmean());
			nc->SetEnergy(hit->emean());
			particle->AddNC(*nc);
			delete nc;
		}
	}
#ifdef DEBUG
	std::cout << "### nHit in NC = " << particle->nNC() << " ###" << std::endl;
#endif

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
			if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
		}
	}
	if(fcds==-1) return true; /* there is no vertex in fiducial volume */
	pCDS* cds=0;
	if(fcds!=-1) cds = particle->cdsi(fcds);

#ifdef DEBUG
	//std::cout << "### Vertex ###" << std::endl;
#endif
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

#ifdef DEBUG
	std::cout << "### NC Hit ###" << std::endl;
#endif

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
	/* CDS */
	TLorentzVector Lcds = cds->GetLorentzVector();
	TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
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
	/* CDS + NC */
	double ncds_mass[2]  = { (Lcds+Ln).M()                           , (cLcds+cLn).M()};
	double ncds_mom[2]   = { (Lcds+Ln).Vect().Mag()                  , (cLcds+cLn).Vect().Mag()};
	double ncds_cost[2]  = { (Lcds+Ln).Vect().CosTheta()             , (cLcds+cLn).Vect().CosTheta()};
	double ncds_phi[2]   = { (Lcds+Ln).Vect().Phi()                  , (cLcds+cLn).Vect().Phi()};
	double ncds_mmass[2] = { (Ltgt+Lbeam-(Lcds+Ln)).M()              , (cLtgt+cLbeam-(cLcds+cLn)).M()};
	double ncds_mmom[2]  = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds+cLn)).Vect().Mag()};
	double ncds_mcost[2] = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds+cLn)).Vect().CosTheta()};

	// n
	for(int f=1; f<2; f++){
		FillHist(Form("n_Mass_%s",frame[f].Data()),nc_mass[f]);
		FillHist(Form("n_Momentum_%s",frame[f].Data()),nc_mom[f]);
		FillHist(Form("n_CosT_%s",frame[f].Data()),nc_cost[f]);
		FillHist(Form("n_Phi_%s",frame[f].Data()),nc_phi[f]);
		FillHist(Form("n_MMass_%s",frame[f].Data()),nc_mmass[f]);
		FillHist(Form("n_MMomentum_%s",frame[f].Data()),nc_mmom[f]);
		FillHist(Form("n_MCosT_%s",frame[f].Data()),nc_mcost[f]);
	}

	// n + pi+
	if(cds->pid()==CDS_PiPlus){
		for(int f=1; f<2; f++){
			FillHist(Form("npip_Mass_%s",frame[f].Data()),ncds_mass[f]);
			FillHist(Form("npip_Momentum_%s",frame[f].Data()),ncds_mom[f]);
			FillHist(Form("npip_CosT_%s",frame[f].Data()),ncds_cost[f]);
			FillHist(Form("npip_Phi_%s",frame[f].Data()),ncds_phi[f]);
			FillHist(Form("npip_MMass_%s",frame[f].Data()),ncds_mmass[f]);
			FillHist(Form("npip_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
			FillHist(Form("npip_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
		}
	}
	// n + pi-
	if(cds->pid()==CDS_PiMinus){
		for(int f=1; f<2; f++){
			FillHist(Form("npim_Mass_%s",frame[f].Data()),ncds_mass[f]);
			FillHist(Form("npim_Momentum_%s",frame[f].Data()),ncds_mom[f]);
			FillHist(Form("npim_CosT_%s",frame[f].Data()),ncds_cost[f]);
			FillHist(Form("npim_Phi_%s",frame[f].Data()),ncds_phi[f]);
			FillHist(Form("npim_MMass_%s",frame[f].Data()),ncds_mmass[f]);
			FillHist(Form("npim_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
			FillHist(Form("npim_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
		}
	}
	// n + k-
	if(cds->pid()==CDS_Kaon){
		for(int f=1; f<2; f++){
			FillHist(Form("nk_Mass_%s",frame[f].Data()),ncds_mass[f]);
			FillHist(Form("nk_Momentum_%s",frame[f].Data()),ncds_mom[f]);
			FillHist(Form("nk_CosT_%s",frame[f].Data()),ncds_cost[f]);
			FillHist(Form("nk_Phi_%s",frame[f].Data()),ncds_phi[f]);
			FillHist(Form("nk_MMass_%s",frame[f].Data()),ncds_mmass[f]);
			FillHist(Form("nk_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
			FillHist(Form("nk_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
		}
	}

	return true;

}

bool MyAnalysisFWDNeutral::FillHist(TString name, double val1)
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

bool MyAnalysisFWDNeutral::FillHist(TString name, double val1, double val2)
{
	TH2F* h2 = (TH2F*)gFile -> Get(name);
	if(h2){ h2 -> Fill(val1,val2);
		return true;
	}
	else {
		return false;
	}
}

bool MyAnalysisFWDNeutral::Initialize(TFile* rtFile, ConfMan* confMan)
{
	rtFile -> cd();

	// FWD neutral Particles
	std::cout << "Define Histograms for FWD neutral particles" << std::endl;
	new TH1F( "FWDN_Overbeta", "1/#beta;1/#beta;Counts", 5000, 0, 5);
	new TH1F( "FWDN_Overbeta_withdECut", "1/#beta;1/#beta;Counts", 5000, 0, 5);
	new TH1F( "FWDN_Overbeta_Selected", "1/#beta;1/#beta;Counts", 5000, 0, 5);
	new TH2F( "FWDN_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
	new TH2F( "FWDN_OverbetavsEnergy_withdECut", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
	new TH2F( "NCHitPosition_XY", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",16,-160,160,15,-75,75);
	new TH2F( "NCHitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
	// Missing mass p(K-,n)X
	std::cout << "Define Histograms for (K-,n)X M.M." << std::endl;
	new TH1F("n_Mass_CM","Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
	new TH1F("n_Momentum_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
	new TH1F("n_CosT_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
	new TH1F("n_Phi_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
	new TH1F("n_MMass_CM", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
	new TH1F("n_MMomentum_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
	new TH1F("n_MCosT_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
	// Invariant mass with FWD
	std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
	// n + pi+
	new TH1F("npip_Mass_CM","Invariant mass of n#pi^{+};IM(n#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
	new TH1F("npip_Momentum_CM","Momentum of n#pi^{+};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
	new TH1F("npip_CosT_CM","cos#theta of n#pi^{+};cos#theta;Counts", 3000, -1.5, 1.5);
	new TH1F("npip_Phi_CM","#phi of n#pi^{+};#phi;Counts", 1600, 0.0, 3.2);
	new TH1F("npip_MMass_CM", "^{3}He(K^{-},n#pi^{+})X missing mass;MM(n#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
	new TH1F("npip_MMomentum_CM", "^{3}He(K^{-},n#pi^{+})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
	new TH1F("npip_MCosT_CM", "^{3}He(K^{-},n#pi^{+})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
	// n + pi-
	new TH1F("npim_Mass_CM","Invariant mass of n#pi^{-};IM(n#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
	new TH1F("npim_Momentum_CM","Momentum of n#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
	new TH1F("npim_CosT_CM","cos#theta of n#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
	new TH1F("npim_Phi_CM","#phi of n#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
	new TH1F("npim_MMass_CM", "^{3}He(K^{-},n#pi^{-})X missing mass;MM(n#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
	new TH1F("npim_MMomentum_CM", "^{3}He(K^{-},n#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
	new TH1F("npim_MCosT_CM", "^{3}He(K^{-},n#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
	// n + k-
	new TH1F("nk_Mass_CM","Invariant mass of nK^{-};IM(nK^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
	new TH1F("nk_Momentum_CM","Momentum of nK^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
	new TH1F("nk_CosT_CM","cos#theta of nK^{-};cos#theta;Counts", 3000, -1.5, 1.5);
	new TH1F("nk_Phi_CM","#phi of nK^{-};#phi;Counts", 1600, 0.0, 3.2);
	new TH1F("nk_MMass_CM", "^{3}He(K^{-},nK^{-})X missing mass;MM(nK^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
	new TH1F("nk_MMomentum_CM", "^{3}He(K^{-},nK^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
	new TH1F("nk_MCosT_CM", "^{3}He(K^{-},nK^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);

	std::cout << "=== End of [MyAnalysisFWDNeutral::Initialize] === " << std::endl;

	return true;
}
