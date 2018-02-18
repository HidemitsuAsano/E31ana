// MyAnalysisCheckChi2.cpp

#include "MyAnalysisCheckChi2.h"

MyAnalysisCheckChi2::MyAnalysisCheckChi2(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  Clear();
}

MyAnalysisCheckChi2::~MyAnalysisCheckChi2()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCheckChi2::Clear()
{
}

bool MyAnalysisCheckChi2::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{

  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  /* Trigger selection */
  if(!header->IsTrig(Trig_KCDH2)){ return true; }

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
	if(particle->nCDS()!=2){
    return true;
	}
	if(particle->nProduct()!=1){
    return true;
	}

	pCDS* cds = particle->product(0);
	int comb = cds->comb();
	if(comb!=pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus) && comb!=pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
		return true;
	}
	TLorentzVector Lcds = cds->GetLorentzVector();
	double mass = Lcds.M();
	TVector3 vtxb = cds->vbeam();
	double vdis = cds->vdis();
	double vbdis = cds->vbdis();
	double pbdca = cds->pbdca();

	bool k0flag = false;
	if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
		k0flag = true;
	}

	int t0seg = beam->t0seg();
	int bhdseg = beam->bhdseg();
	double beam_mom = beam->mom();

	double blc1_chi2 = beam->blc1chi2();
	double blc2_chi2 = beam->blc2chi2();
	double bpc_chi2 = beam->bpcchi2();
	double beam_chi2 = beam->beamchi2();

	double cds1_chi2 = -1;
	double cds2_chi2 = -1;

	for(int it=0; it<particle->nCDS(); it++){
		if(k0flag){
			cds1_chi2 = particle->pim(0)->chi();
			cds2_chi2 = particle->pip(0)->chi();
		}
		else{
			cds1_chi2 = particle->pim(0)->chi();
			cds2_chi2 = particle->proton(0)->chi();
		}
	}
	std::string part[2] = {"K0","Lambda"};
	int ip = -1;
	if(k0flag){ ip=0; }
	else { ip=1; }
	if(ip==-1) return true;

	double chi2[6] = {blc1_chi2, blc2_chi2, bpc_chi2, cds1_chi2, cds2_chi2, beam_chi2};
	std::string dcname[6] = {"BLC1","BLC2","BPC","CDS1","CDS2","Beam"};

	///////////////
	// Fill Chi2 //
	///////////////
	for(int idc=0; idc<6; idc++){
		FillHist(Form("%s_%s_Chi2",part[ip].data(),dcname[idc].data()), chi2[idc]);
		for(int jdc=0; jdc<6; jdc++){
			if(idc!=jdc){
				FillHist(Form("%s_%s_Chi2_vs_%s_Chi2",part[ip].data(),dcname[idc].data(),dcname[jdc].data()), chi2[idc], chi2[jdc]);
			}
		}
	}

	///////////////////////////////////////////////
	// Fill Vertex & Mass & Beam momentrum & DCA //
	///////////////////////////////////////////////
	FillHist(Form("%s_Vertex_XY",part[ip].data()), vtxb.X(), vtxb.Y());
	FillHist(Form("%s_Vertex_ZX",part[ip].data()), vtxb.Z(), vtxb.X());
	FillHist(Form("%s_Vertex_ZY",part[ip].data()), vtxb.Z(), vtxb.Y());
	FillHist(Form("%s_Vertex_Z",part[ip].data()), vtxb.Z());
	if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
		FillHist(Form("%s_Vertex_XY_ZCut",part[ip].data()), vtxb.X(), vtxb.Y());
	}
	if(-0.5<vtxb.X()&&vtxb.X()<0.5){
		FillHist(Form("%s_Vertex_ZY_XCut",part[ip].data()), vtxb.Z(), vtxb.Y());
	}
	if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
		FillHist(Form("%s_Vertex_ZX_YCut",part[ip].data()), vtxb.Z(), vtxb.X());
	}
	if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
		FillHist(Form("%s_Vertex_Z_XYCut",part[ip].data()), vtxb.Z());
	}
	if(-0.5<vtxb.Y()&&vtxb.Y()<0.5&&-0.5<vtxb.Z()&&vtxb.Z()<0.5){
		FillHist(Form("%s_Vertex_X_YZCut",part[ip].data()), vtxb.X());
	}
	if(-0.5<vtxb.Z()&&vtxb.Z()<0.5&&-0.5<vtxb.X()&&vtxb.X()<0.5){
		FillHist(Form("%s_Vertex_Y_ZXCut",part[ip].data()), vtxb.Y());
	}
	FillHist(Form("%s_Mass",part[ip].data()), mass);
	if(t0seg==3&&bhdseg==10){
		FillHist(Form("%s_BeamMomentum",part[ip].data()), beam_mom);
	}
	FillHist(Form("%s_VDIS",part[ip].data()), vdis);
	FillHist(Form("%s_VBDIS",part[ip].data()), vbdis);
	FillHist(Form("%s_PBDCA",part[ip].data()), pbdca);


	/////////////////////
	// Check Each Chi2 //
	/////////////////////
	double chi2_cut[6]={10.0,10.0,10.0,30.0,30.0,99999};
	for(int idc=0; idc<6; idc++){
		bool flag = true;
		for(int jdc=0; jdc<6; jdc++){
			if(idc!=jdc&&chi2[jdc]>chi2_cut[jdc]) flag = false;
		}
		if(!flag) continue;
		for(int jdc=0; jdc<6; jdc++){
			FillHist(Form("%s_%s_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[jdc].data()), chi2[jdc]);
		}


		FillHist(Form("%s_%s_Vertex_XY",dcname[idc].data(),part[ip].data()), vtxb.X(), vtxb.Y());
		FillHist(Form("%s_%s_Vertex_ZX",dcname[idc].data(),part[ip].data()), vtxb.Z(), vtxb.X());
		FillHist(Form("%s_%s_Vertex_ZY",dcname[idc].data(),part[ip].data()), vtxb.Z(), vtxb.Y());
		FillHist(Form("%s_%s_Vertex_Z",dcname[idc].data(),part[ip].data()), vtxb.Z());
		if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
			FillHist(Form("%s_%s_Vertex_XY_ZCut",dcname[idc].data(),part[ip].data()), vtxb.X(), vtxb.Y());
		}
		if(-0.5<vtxb.X()&&vtxb.X()<0.5){
			FillHist(Form("%s_%s_Vertex_ZY_XCut",dcname[idc].data(),part[ip].data()), vtxb.Z(), vtxb.Y());
		}
		if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
			FillHist(Form("%s_%s_Vertex_ZX_YCut",dcname[idc].data(),part[ip].data()), vtxb.Z(), vtxb.X());
		}
		if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
			FillHist(Form("%s_%s_Vertex_Z_XYCut",dcname[idc].data(),part[ip].data()), vtxb.Z());
			FillHist(Form("%s_%s_Vertex_Z_XYCut_vs_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[idc].data()), vtxb.Z(), chi2[idc]);
		}
		if(-0.5<vtxb.Y()&&vtxb.Y()<0.5&&-0.5<vtxb.Z()&&vtxb.Z()<0.5){
			FillHist(Form("%s_%s_Vertex_X_YZCut",dcname[idc].data(),part[ip].data()), vtxb.X());
			FillHist(Form("%s_%s_Vertex_X_YZCut_vs_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[idc].data()), vtxb.X(), chi2[idc]);
		}
		if(-0.5<vtxb.Z()&&vtxb.Z()<0.5&&-0.5<vtxb.X()&&vtxb.X()<0.5){
			FillHist(Form("%s_%s_Vertex_Y_ZXCut",dcname[idc].data(),part[ip].data()), vtxb.Y());
			FillHist(Form("%s_%s_Vertex_Y_ZXCut_vs_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[idc].data()), vtxb.Y(), chi2[idc]);
		}
		FillHist(Form("%s_%s_Mass",dcname[idc].data(),part[ip].data()), mass);
		FillHist(Form("%s_%s_Mass_vs_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[idc].data()), mass, chi2[idc]);
		if(t0seg==3&&bhdseg==10){
			FillHist(Form("%s_%s_BeamMomentum",dcname[idc].data(),part[ip].data()), beam_mom);
			FillHist(Form("%s_%s_BeamMomentum_vs_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[idc].data()), beam_mom, chi2[idc]);
		}
		FillHist(Form("%s_%s_VDIS",dcname[idc].data(),part[ip].data()), vdis);
		FillHist(Form("%s_%s_VBDIS",dcname[idc].data(),part[ip].data()), vbdis);
		FillHist(Form("%s_%s_PBDCA",dcname[idc].data(),part[ip].data()), pbdca);
		FillHist(Form("%s_%s_VDIS_vs_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[idc].data()), vdis, chi2[idc]);
		FillHist(Form("%s_%s_VBDIS_vs_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[idc].data()), vbdis, chi2[idc]);
		FillHist(Form("%s_%s_PBDCA_vs_%s_Chi2",dcname[idc].data(),part[ip].data(),dcname[idc].data()), pbdca, chi2[idc]);
	}


	return true;

}

bool MyAnalysisCheckChi2::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCheckChi2::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCheckChi2::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCheckChi2::FillHist(TString name, TString val1, TString val2, int weight)
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

bool MyAnalysisCheckChi2::Initialize(ConfMan* confMan)
{
	std::cout << "### MyAnalysisCheckChi2::Initialize ###" << std::endl;

	std::string ofname = confMan->GetOutFileName();
	ofname.insert(ofname.find(".root"),"_anaCheckChi2");

	rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
	rtFile -> cd();

	std::string part[2]={"K0","Lambda"};
	std::string dcname[6]={"BLC1","BLC2","BPC","CDS1","CDS2","Beam"};

	for(int ipart=0; ipart<2; ipart++){

		for(int idc=0; idc<6; idc++){
			new TH1F( Form("%s_%s_Chi2",part[ipart].data(),dcname[idc].data()), Form("#chi^{2} %s;#chi^{2}/ndf;Counts",dcname[idc].data()), 1000, 0.0, 100.0 );
			for(int jdc=0; jdc<6; jdc++){
				if(idc!=jdc){
					new TH2F( Form("%s_%s_Chi2_vs_%s_Chi2",part[ipart].data(),dcname[idc].data(),dcname[jdc].data()), Form("#chi^{2} %s vs. #chi^{2} %s;#chi^{2}/ndf (%s);#chi^{2}/ndf (%s)",dcname[idc].data(), dcname[jdc].data(), dcname[idc].data(), dcname[jdc].data()), 500, 0.0, 100.0, 500, 0.0, 100.0 );
				}
			}
		}

		new TH2F( Form("%s_Vertex_XY",part[ipart].data()), "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
		new TH2F( Form("%s_Vertex_ZX",part[ipart].data()), "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
		new TH2F( Form("%s_Vertex_ZY",part[ipart].data()), "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
		new TH1F( Form("%s_Vertex_Z",part[ipart].data()), "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
		new TH2F( Form("%s_Vertex_XY_ZCut",part[ipart].data()), "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
		new TH1F( Form("%s_Vertex_Z_XYCut",part[ipart].data()), "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
		new TH1F( Form("%s_Vertex_X_YZCut",part[ipart].data()), "Vertex X axis;x (cm);Counts",  600, -15, 15 );
		new TH1F( Form("%s_Vertex_Y_ZXCut",part[ipart].data()), "Vertex Y axis;y (cm);Counts",  600, -15, 15 );

		new TH1F( Form("%s_Mass",part[ipart].data()), Form("Invariant Mass (%s);IM (GeV/c^{2});Counts", part[ipart].data()), 5000, 0.0, 5.0);

		new TH1F( Form("%s_BeamMomentum",part[ipart].data()), Form("Beam momentum (%s);Beam momentum (GeV/c);Counts", part[ipart].data()), 2000, 0.5, 1.5);

		new TH1F( Form("%s_VDIS",part[ipart].data()), Form("VDIS (%s);VDIS (cm);Counts", part[ipart].data()), 2000, 0.0, 20.0);
		new TH1F( Form("%s_VBDIS",part[ipart].data()), Form("VBDIS (%s);VBDIS (cm);Counts", part[ipart].data()), 2000, 0.0, 20.0);
		new TH1F( Form("%s_PBDCA",part[ipart].data()), Form("PBDCA (%s);PBDCA (cm);Counts", part[ipart].data()), 2000, 0.0, 20.0);

		for(int idc=0; idc<6; idc++){
			for(int jdc=0; jdc<6; jdc++){
				new TH1F( Form("%s_%s_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[jdc].data()), Form("#chi^{2} %s;#chi^{2}/ndf;Counts",dcname[idc].data()), 1000, 0.0, 100.0 );
			}
			new TH2F( Form("%s_%s_Vertex_XY",dcname[idc].data(),part[ipart].data()), "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
			new TH2F( Form("%s_%s_Vertex_ZX",dcname[idc].data(),part[ipart].data()), "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
			new TH2F( Form("%s_%s_Vertex_ZY",dcname[idc].data(),part[ipart].data()), "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
			new TH1F( Form("%s_%s_Vertex_Z",dcname[idc].data(),part[ipart].data()), "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
			new TH2F( Form("%s_%s_Vertex_XY_ZCut",dcname[idc].data(),part[ipart].data()), "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
			new TH1F( Form("%s_%s_Vertex_Z_XYCut",dcname[idc].data(),part[ipart].data()), "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
			new TH1F( Form("%s_%s_Vertex_X_YZCut",dcname[idc].data(),part[ipart].data()), "Vertex X axis;x (cm);Counts",  600, -15, 15 );
			new TH1F( Form("%s_%s_Vertex_Y_ZXCut",dcname[idc].data(),part[ipart].data()), "Vertex Y axis;y (cm);Counts",  600, -15, 15 );

			new TH1F( Form("%s_%s_Mass",dcname[idc].data(),part[ipart].data()), Form("Invariant Mass (%s);Mass (GeV/c^{2});Counts", part[ipart].data()), 5000, 0.0, 5.0);

			new TH1F( Form("%s_%s_BeamMomentum",dcname[idc].data(),part[ipart].data()), Form("Beam momentum (%s);Beam momentum (GeV/c);Counts", part[ipart].data()), 2000, 0.5, 1.5);

			new TH2F( Form("%s_%s_Vertex_Z_XYCut_vs_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[idc].data()), Form("Vertex Z axis vs. %s #chi^{2};z (cm);#chi^{2}/ndf",dcname[idc].data()), 1200, -30, 30, 100, 0.0, 100.0 );
			new TH2F( Form("%s_%s_Vertex_X_YZCut_vs_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[idc].data()), Form("Vertex X axis vs. %s #chi^{2};x (cm);#chi^{2}/ndf",dcname[idc].data()),  600, -15, 15, 100, 0.0, 100.0 );
			new TH2F( Form("%s_%s_Vertex_Y_ZXCut_vs_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[idc].data()), Form("Vertex Y axis vs. %s #chi^{2};y (cm);#chi^{2}/ndf",dcname[idc].data()),  600, -15, 15, 100, 0.0, 100.0 );
			new TH2F( Form("%s_%s_Mass_vs_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[idc].data()), Form("Invariant Mass (%s) vs. %s #chi^{2};Mass (GeV/c^{2});#chi^{2}/ndf", part[ipart].data(),dcname[idc].data()), 2000, 0.0, 2.0, 100, 0.0, 100.0);
			new TH2F( Form("%s_%s_BeamMomentum_vs_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[idc].data()), Form("Beam momentum (%s) vs. %s #chi^{2};Beam momentum (GeV/c^{2});#chi^{2}/ndf", part[ipart].data(),dcname[idc].data()), 2000, 0.5, 1.5, 100, 0.0, 100.0);
			new TH2F( Form("%s_%s_VDIS_vs_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[idc].data()), Form("VDIS (%s) vs. %s #chi^{2};VDIS (cm);Counts", part[ipart].data(),dcname[idc].data()), 2000, 0.0, 20.0, 100, 0.0, 100.0);
			new TH2F( Form("%s_%s_VBDIS_vs_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[idc].data()), Form("VBDIS (%s) vs. %s #chi^{2};VBDIS (cm);Counts", part[ipart].data(),dcname[idc].data()), 2000, 0.0, 20.0, 100, 0.0, 100.0);
			new TH2F( Form("%s_%s_PBDCA_vs_%s_Chi2",dcname[idc].data(),part[ipart].data(),dcname[idc].data()), Form("PBDCA (%s) vs. %s #chi^{2};PBDCA (cm);Counts", part[ipart].data(),dcname[idc].data()), 2000, 0.0, 20.0, 100, 0.0, 100.0);

		}

	}
	return true;
}
