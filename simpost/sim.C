void phasespace();
void CreateHistograms();
bool checkacceptance(TLorentzVector*,TLorentzVector*,TLorentzVector*);

static TString pdfname = "out.pdf";
static TFile* rtFile;
static int icanvas = 0;
TH1F* h1;
TH2F* h2;

static TRandom*  rndm;

static const double mom_beam   = 1.0;

static const double c_lv = 29.9792458; //cm/ns
static const double m_neutron  = 0.939565346;
static const double m_proton   = 0.938272013;
static const double m_kaon     = 0.493677;
static const double m_deuteron = 1.875612928;

static const double cdh_lv=-14.3; //cm/ns
static const double cdh_w = 10.0; //cm
static const double cdh_l = 79.0; //cm
static const double cdh_r = 56.0; //cm
static const double costheta_max = cdh_l/2.0 / sqrt(cdh_l*cdh_l/4.0+cdh_r*cdh_r);
static const double costheta_min = -costheta_max;
static const double momentum_min = 0.2; //GeV/c

static double tof_offs   = 0.0; //ns
static double tsub_offs  = 0.0; //ns

static const double resl_momentum = 0.001; //GeV/c
static const double resl_time     = 0.140; //ns


void sim(double in_tof_offs,double in_tsub_offs){
	tof_offs = in_tof_offs;
	tsub_offs = in_tsub_offs;
	gROOT->SetBatch(true);	
	rndm = new TRandom();
	if(tof_offs>=0){
		if(tsub_offs>=0){
			rtFile = new TFile(Form("tofoffs_%02d_tsuboffs_%02d.root",(int)(tof_offs*10),(int)(tsub_offs*10)),"RECREATE");
			pdfname = Form("tofoffs_%02d_tsuboffs_%02d.pdf",(int)(tof_offs*10),(int)(tsub_offs*10));
		}
		else{
			rtFile = new TFile(Form("tofoffs_%02d_tsuboffs_m%02d.root",(int)(tof_offs*10),(int)(-tsub_offs*10)),"RECREATE");
			pdfname = Form("tofoffs_%02d_tsuboffs_m%02d.pdf",(int)(tof_offs*10),(int)(-tsub_offs*10));
		}
	}
	else{
		if(tsub_offs>=0){
			rtFile = new TFile(Form("tofoffs_m%02d_tsuboffs_%02d.root",(int)(-tof_offs*10),(int)(tsub_offs*10)),"RECREATE");
			pdfname = Form("tofoffs_m%02d_tsuboffs_%02d.pdf",(int)(-tof_offs*10),(int)(tsub_offs*10));
		}
		else{
			rtFile = new TFile(Form("tofoffs_m%02d_tsuboffs_m%02d.root",(int)(-tof_offs*10),(int)(-tsub_offs*10)),"RECREATE");
			pdfname = Form("tofoffs_m%02d_tsuboffs_m%02d.pdf",(int)(-tof_offs*10),(int)(-tsub_offs*10));
		}
	}
	CreateHistograms();
	TCanvas* c = new TCanvas(Form("c%d",icanvas),Form("c%d",icanvas),500,500); icanvas++;
	c->Print(pdfname+"[");
	phasespace();
	c->Print(pdfname+"]");
	rtFile->Write();
	rtFile->Close();
}

void CreateHistograms(){
	rtFile->cd();
	// for neutron
	new TH1F("n_p","momentum of neutron;p (GeV/c);Counts",1000,0.0,1.0);
	new TH1F("n_E","energy of neutron;E (GeV);Counts",1000,0.5,1.5);
	new TH1F("n_costheta","cos#theta of neutron;cos#theta;Counts",2000,-1.0,1.0);
	new TH2F("n_pos_cdh","position at CDH of neutron;z (cm);r (cm)",1000,-50,50,1000,0,100);

	new TH1F("n_p_acc","momentum of neutron;p (GeV/c);Counts",1000,0.0,1.0);
	new TH1F("n_E_acc","energy of neutron;E (GeV);Counts",1000,0.5,1.5);
	new TH1F("n_costheta_acc","cos#theta of neutron;cos#theta;Counts",2000,-1.0,1.0);
	new TH2F("n_pos_cdh_acc","position at CDH of neutron;z (cm);r (cm)",1000,-50,50,1000,0,100);

	// for proton
	new TH1F("p_p","momentum of proton;p (GeV/c);Counts",1000,0.0,1.0);
	new TH1F("p_E","energy of proton;E (GeV);Counts",1000,0.5,1.5);
	new TH1F("p_costheta","cos#theta of proton;cos#theta;Counts",2000,-1.0,1.0);
	new TH2F("p_pos_cdh","position at CDH of proton;z (cm);r (cm)",1000,-50,50,1000,0,100);

	new TH1F("p_p_acc","momentum of proton;p (GeV/c);Counts",1000,0.0,1.0);
	new TH1F("p_E_acc","energy of proton;E (GeV);Counts",1000,0.5,1.5);
	new TH1F("p_costheta_acc","cos#theta of proton;cos#theta;Counts",2000,-1.0,1.0);
	new TH2F("p_pos_cdh_acc","position at CDH of proton;z (cm);r (cm)",1000,-50,50,1000,0,100);

	// for koan
	new TH1F("k_p","momentum of kaon;p (GeV/c);Counts",1000,0.0,1.0);
	new TH1F("k_E","energy of kaon;E (GeV);Counts",1000,0.5,1.5);
	new TH1F("k_costheta","cos#theta of kaon;cos#theta;Counts",2000,-1.0,1.0);
	new TH2F("k_pos_cdh","position at CDH of kaon;z (cm);r (cm)",1000,-50,50,1000,0,100);

	new TH1F("k_p_acc","momentum of kaon;p (GeV/c);Counts",1000,0.0,1.0);
	new TH1F("k_E_acc","energy of kaon;E (GeV);Counts",1000,0.5,1.5);
	new TH1F("k_costheta_acc","cos#theta of kaon;cos#theta;Counts",2000,-1.0,1.0);
	new TH2F("k_pos_cdh_acc","position at CDH of kaon;z (cm);r (cm)",1000,-50,50,1000,0,100);

	// for resl_koan


	// for others
	new TH1F("mm_kn","missing mass of nK^{-};MM(nK^{-});Counts",1000,0.5,1.5);
	new TH2F("mm_kn_vs_cdh_z","missing mass of nK^{-} vs. CDH z position;z (cm);MM(nK^{-})",100,-50,50,100,0.4,1.4);

}

TLorentzVector* resl_kaon(TLorentzVector* p){
	TLorentzVector* p_meas = new TLorentzVector();
	double mom_gene = p->P();
	double mom_meas = rndm->Gaus(mom_gene,resl_momentum);
	TVector3 mom_vect_meas = p->Vect();
	mom_vect_meas *= mom_meas/mom_gene;
	p_meas->SetVect(mom_vect_meas);
	double e_meas = sqrt(m_kaon*m_kaon+mom_meas*mom_meas);
	p_meas->SetE(e_meas);
	return p_meas;	
}

TLorentzVector* resl_neutron(TLorentzVector* p){
	TLorentzVector* p_meas = new TLorentzVector();
	double n_dir = p->Pz() / p->Pt();
	double n_cdh_pos_z_gene = n_dir * cdh_r;
	double tsub_gene = n_cdh_pos_z_gene/cdh_lv;
	double tsub_meas = tsub_offs + rndm->Gaus(tsub_gene,resl_time);
	double n_cdh_pos_z_meas = tsub_meas * cdh_lv;
	double fl_meas = sqrt(cdh_r*cdh_r+n_cdh_pos_z_meas*n_cdh_pos_z_meas);
	double fl_gene = sqrt(cdh_r*cdh_r+n_cdh_pos_z_gene*n_cdh_pos_z_gene);
	double beta_gene = p->P()/p->E();
	double tof_gene = fl_gene/beta_gene/c_lv;
	double tof_meas = tof_offs + rndm->Gaus(tof_gene,resl_time);
	double beta_meas = fl_meas/tof_meas/c_lv;
	double mom_meas = beta_meas/sqrt(1-beta_meas*beta_meas)*m_neutron;
	TVector3 mom_vect_meas = p->Vect();
	mom_vect_meas.SetZ(n_cdh_pos_z_meas/cdh_r * p->Pt());
	double mom_gene = mom_vect_meas.Mag();
	mom_vect_meas *= mom_meas/mom_gene;
	p_meas->SetVect(mom_vect_meas);
	double e_meas = sqrt(m_neutron*m_neutron+mom_meas*mom_meas);
	p_meas->SetE(e_meas);
	return p_meas;	
}

void phasespace(){
	if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");
	TLorentzVector* target = new TLorentzVector(0.0, 0.0, 0.0, m_deuteron);
	TLorentzVector* beam = new TLorentzVector(0.0, 0.0, mom_beam, sqrt(mom_beam*mom_beam+m_kaon*m_kaon));
	TLorentzVector* W = new TLorentzVector(*target + *beam);
	Double_t masses[3] = { m_neutron, m_proton, -m_kaon} ;
	TGenPhaseSpace event;
	event.SetDecay(*W, 3, masses);
	for (Int_t n=0;n<1000000;n++) {
		Double_t weight = event.Generate();
		TLorentzVector *pNeutron = event.GetDecay(0);
		TLorentzVector *pProton  = event.GetDecay(1);
		TLorentzVector *pKaon    = event.GetDecay(2);
		bool flag_acc = checkacceptance(pNeutron,pProton,pKaon);

		// for neutron
		h1=(TH1F*)rtFile->Get("n_p"); h1->Fill(pNeutron->P());
		h1=(TH1F*)rtFile->Get("n_E"); h1->Fill(pNeutron->E());
		h1=(TH1F*)rtFile->Get("n_costheta"); h1->Fill(pNeutron->CosTheta());
		double n_dir = pNeutron->Pz() / pNeutron->Pt();
		double n_cdh_pos_z = n_dir * cdh_r;
		h2=(TH2F*)rtFile->Get("n_pos_cdh"); h2->Fill(n_cdh_pos_z,cdh_r);
		if(flag_acc){
			h1=(TH1F*)rtFile->Get("n_p_acc"); h1->Fill(pNeutron->P());
			h1=(TH1F*)rtFile->Get("n_E_acc"); h1->Fill(pNeutron->E());
			h1=(TH1F*)rtFile->Get("n_costheta_acc"); h1->Fill(pNeutron->CosTheta());
			h2=(TH2F*)rtFile->Get("n_pos_cdh_acc"); h2->Fill(n_cdh_pos_z,cdh_r);
		}

		// for proton
		h1=(TH1F*)rtFile->Get("p_p"); h1->Fill(pProton->P());
		h1=(TH1F*)rtFile->Get("p_E"); h1->Fill(pProton->E());
		h1=(TH1F*)rtFile->Get("p_costheta"); h1->Fill(pProton->CosTheta());
		double p_dir = pProton->Pz() / pProton->Pt();
		double p_cdh_pos_z = p_dir * cdh_r;
		h2=(TH2F*)rtFile->Get("p_pos_cdh"); h2->Fill(p_cdh_pos_z,cdh_r);
		if(flag_acc){
			h1=(TH1F*)rtFile->Get("p_p_acc"); h1->Fill(pProton->P());
			h1=(TH1F*)rtFile->Get("p_E_acc"); h1->Fill(pProton->E());
			h1=(TH1F*)rtFile->Get("p_costheta_acc"); h1->Fill(pProton->CosTheta());
			h2=(TH2F*)rtFile->Get("p_pos_cdh_acc"); h2->Fill(p_cdh_pos_z,cdh_r);
		}

		// for kaon
		h1=(TH1F*)rtFile->Get("k_p"); h1->Fill(pKaon->P());
		h1=(TH1F*)rtFile->Get("k_E"); h1->Fill(pKaon->E());
		h1=(TH1F*)rtFile->Get("k_costheta"); h1->Fill(pKaon->CosTheta());
		double k_dir = pKaon->Pz() / pKaon->Pt();
		double k_cdh_pos_z = k_dir * cdh_r;
		h2=(TH2F*)rtFile->Get("k_pos_cdh"); h2->Fill(k_cdh_pos_z,cdh_r);
		if(flag_acc){
			h1=(TH1F*)rtFile->Get("k_p_acc"); h1->Fill(pKaon->P());
			h1=(TH1F*)rtFile->Get("k_E_acc"); h1->Fill(pKaon->E());
			h1=(TH1F*)rtFile->Get("k_costheta_acc"); h1->Fill(pKaon->CosTheta());
			h2=(TH2F*)rtFile->Get("k_pos_cdh_acc"); h2->Fill(k_cdh_pos_z,cdh_r);
		}

		if(flag_acc){
			TLorentzVector* pKaon_meas = resl_kaon(pKaon);  
			TLorentzVector* pNeutron_meas = resl_neutron(pNeutron);  
			TLorentzVector* pProton_meas = new TLorentzVector(*target + *beam - *pKaon_meas - *pNeutron_meas);
			double mmass = pProton_meas->M();
			double n_dir_meas = pNeutron_meas->Pz() / pNeutron_meas->Pt();
			double n_cdh_pos_z_meas = n_dir_meas * cdh_r;
			h1=(TH1F*)rtFile->Get("mm_kn"); h1->Fill(mmass);
			h2=(TH2F*)rtFile->Get("mm_kn_vs_cdh_z"); h2->Fill(n_cdh_pos_z_meas,mmass);
		}

	}

	TCanvas* c = new TCanvas(Form("c%d",icanvas),Form("c%d",icanvas),1000,1000); icanvas++;
	c->Divide(2,2);
	// for neutron
	c->cd(1);	
	h1=(TH1F*)rtFile->Get("n_p"); h1->Draw();
	c->cd(2);
	h1=(TH1F*)rtFile->Get("n_E"); h1->Draw();
	c->cd(3);
	h1=(TH1F*)rtFile->Get("n_costheta"); h1->Draw();
	c->cd(4);
	h2=(TH2F*)rtFile->Get("n_pos_cdh"); h2->Draw("col");
	c->Print(pdfname);

	c->cd(1);	
	h1=(TH1F*)rtFile->Get("n_p_acc"); h1->Draw();
	c->cd(2);
	h1=(TH1F*)rtFile->Get("n_E_acc"); h1->Draw();
	c->cd(3);
	h1=(TH1F*)rtFile->Get("n_costheta_acc"); h1->Draw();
	c->cd(4);
	h2=(TH2F*)rtFile->Get("n_pos_cdh_acc"); h2->Draw("col");
	c->Print(pdfname);

	// for proton
	c->cd(1);	
	h1=(TH1F*)rtFile->Get("p_p"); h1->Draw();
	c->cd(2);
	h1=(TH1F*)rtFile->Get("p_E"); h1->Draw();
	c->cd(3);
	h1=(TH1F*)rtFile->Get("p_costheta"); h1->Draw();
	c->cd(4);
	h2=(TH2F*)rtFile->Get("p_pos_cdh"); h2->Draw("col");
	c->Print(pdfname);

	c->cd(1);	
	h1=(TH1F*)rtFile->Get("p_p_acc"); h1->Draw();
	c->cd(2);
	h1=(TH1F*)rtFile->Get("p_E_acc"); h1->Draw();
	c->cd(3);
	h1=(TH1F*)rtFile->Get("p_costheta_acc"); h1->Draw();
	c->cd(4);
	h2=(TH2F*)rtFile->Get("p_pos_cdh_acc"); h2->Draw("col");
	c->Print(pdfname);

	// for kaon
	c->cd(1);	
	h1=(TH1F*)rtFile->Get("k_p"); h1->Draw();
	c->cd(2);
	h1=(TH1F*)rtFile->Get("k_E"); h1->Draw();
	c->cd(3);
	h1=(TH1F*)rtFile->Get("k_costheta"); h1->Draw();
	c->cd(4);
	h2=(TH2F*)rtFile->Get("k_pos_cdh"); h2->Draw("col");
	c->Print(pdfname);

	c->cd(1);	
	h1=(TH1F*)rtFile->Get("k_p_acc"); h1->Draw();
	c->cd(2);
	h1=(TH1F*)rtFile->Get("k_E_acc"); h1->Draw();
	c->cd(3);
	h1=(TH1F*)rtFile->Get("k_costheta_acc"); h1->Draw();
	c->cd(4);
	h2=(TH2F*)rtFile->Get("k_pos_cdh_acc"); h2->Draw("col");
	c->Print(pdfname);

	delete c;
	c = new TCanvas(Form("c%d",icanvas),Form("c%d",icanvas),500,500); icanvas++;

	// for kaon
	c->cd(1);	
	h1=(TH1F*)rtFile->Get("mm_kn"); h1->Draw();
	c->Print(pdfname);
	h2=(TH2F*)rtFile->Get("mm_kn_vs_cdh_z"); h2->Draw("col");
	c->Print(pdfname);
}

bool checkacceptance(TLorentzVector* n, TLorentzVector* p, TLorentzVector* k){
	bool flag = true;
	// for neutron
	double costheta_n = n->CosTheta();	
	if(costheta_n<=costheta_min || costheta_max<=costheta_n) flag = false;

	// for kaon
	double momentum_k = k->P();
	double costheta_k = k->CosTheta();	
	if(momentum_k<=momentum_min) flag = false;
	if(costheta_k<=costheta_min || costheta_max<=costheta_k) flag = false;

	// for proton
	double momentum_p = p->P();
	double costheta_p = p->CosTheta();
	if( costheta_min<costheta_p && costheta_p<costheta_max ){
		if(momentum_p>=momentum_min) flag = false;
	}

	return flag;
}
