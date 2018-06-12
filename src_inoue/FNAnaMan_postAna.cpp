#include "FNAnaMan.h"

void FNAnaMan::postAna(TFile *f, const std::string &dirname, int index){
  TDirectory* dir=fFileOut->mkdir(dirname.c_str()); dir-> cd();
  const int nbin=1000/fRebinF;
  for( int i=0; i<fInfoSim->nBin(); i++ ){
    dir-> Append(new TH2F(Form("KNpim_KNpip_MM_%d", i),       "", 1000, 1.0, 2.0, 1000, 1.0, 2.0));
    dir-> Append(new TH2F(Form("KNpim_KNpip_MM_woAll_%d", i), "", 1000, 1.0, 2.0, 1000, 1.0, 2.0));
    dir-> Append(new TH2F(Form("KNpim_KNpip_MM_K0_%d", i),    "", 1000, 1.0, 2.0, 1000, 1.0, 2.0));
    dir-> Append(new TH2F(Form("KNpim_KNpip_MM_Sm_%d", i),    "", 1000, 1.0, 2.0, 1000, 1.0, 2.0));
    dir-> Append(new TH2F(Form("KNpim_KNpip_MM_Sp_%d", i),    "", 1000, 1.0, 2.0, 1000, 1.0, 2.0));

    dir-> Append(new TH1F(Form("KNpim_MM_%d", i),       "", nbin, 1.0, 2.0));
    dir-> Append(new TH1F(Form("KNpim_MM_woAll_%d", i), "", nbin, 1.0, 2.0));
    dir-> Append(new TH1F(Form("KNpim_MM_K0_%d", i),    "", nbin, 1.0, 2.0));
    dir-> Append(new TH1F(Form("KNpim_MM_Sm_%d", i),    "", nbin, 1.0, 2.0));
    dir-> Append(new TH1F(Form("KNpim_MM_Sp_%d", i),    "", nbin, 1.0, 2.0));

    dir-> Append(new TH1F(Form("KNpip_MM_%d", i),       "", nbin, 1.0, 2.0));
    dir-> Append(new TH1F(Form("KNpip_MM_woAll_%d", i), "", nbin, 1.0, 2.0));
    dir-> Append(new TH1F(Form("KNpip_MM_K0_%d", i),    "", nbin, 1.0, 2.0));
    dir-> Append(new TH1F(Form("KNpip_MM_Sm_%d", i),    "", nbin, 1.0, 2.0));
    dir-> Append(new TH1F(Form("KNpip_MM_Sp_%d", i),    "", nbin, 1.0, 2.0));

    dir-> Append(new TH1F(Form("dca_beam_pim_%d", i), "",      1000, 0.0, 5.0));
    dir-> Append(new TH1F(Form("dca_beam_pim_K0_%d", i), "",   1000, 0.0, 5.0));
    dir-> Append(new TH1F(Form("dca_beam_pim_Sm_%d", i), "",   1000, 0.0, 5.0));
    dir-> Append(new TH1F(Form("dca_beam_pim_Sp_%d", i), "",   1000, 0.0, 5.0));
    dir-> Append(new TH1F(Form("dca_beam_pim_woAll_%d", i), "", 1000, 0.0, 5.0));

    dir-> Append(new TH1F(Form("dca_beam_pip_%d", i), "",      1000, 0.0, 5.0));
    dir-> Append(new TH1F(Form("dca_beam_pip_K0_%d", i), "",   1000, 0.0, 5.0));
    dir-> Append(new TH1F(Form("dca_beam_pip_Sm_%d", i), "",   1000, 0.0, 5.0));
    dir-> Append(new TH1F(Form("dca_beam_pip_Sp_%d", i), "",   1000, 0.0, 5.0));
    dir-> Append(new TH1F(Form("dca_beam_pip_woAll_%d", i), "", 1000, 0.0, 5.0));
  }
  dir-> Append(new TH1F("KNpipi_MM", "",   2000, 0.0, 2.0));

  dir-> Append(new TH1F("KN_MM", "",       1000, 1.0, 2.0));
  dir-> Append(new TH1F("KN_MM_K0", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KN_MM_Sm", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KN_MM_Sp", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KN_MM_woAll", "", 1000, 1.0, 2.0));

  dir-> Append(new TH2F("KN_MM_IM_npipi_K0", "", 1000, 1.0, 2.0, 1000, 1.0, 2.0));
  dir-> Append(new TH2F("KN_MM_IM_npipi_Sm", "", 1000, 1.0, 2.0, 1000, 1.0, 2.0));
  dir-> Append(new TH2F("KN_MM_IM_npipi_Sp", "", 1000, 1.0, 2.0, 1000, 1.0, 2.0));

  dir-> Append(new TH1F("dca_beam_pim", "",       1000, 0.0, 5.0));
  dir-> Append(new TH1F("dca_beam_pim_K0", "",    1000, 0.0, 5.0));
  dir-> Append(new TH1F("dca_beam_pim_Sm", "",    1000, 0.0, 5.0));
  dir-> Append(new TH1F("dca_beam_pim_Sp", "",    1000, 0.0, 5.0));
  dir-> Append(new TH1F("dca_beam_pim_woAll", "", 1000, 0.0, 5.0));

  dir-> Append(new TH1F("pipi_cosOP_K0", "", 1000, -1.0, 1.0));
  dir-> Append(new TH1F("dca_pipi_K0", "",       1000, 0.0, 5.0));

  dir-> Append(new TH1F("dca_beam_pip", "",       1000, 0.0, 5.0));
  dir-> Append(new TH1F("dca_beam_pip_K0", "",    1000, 0.0, 5.0));
  dir-> Append(new TH1F("dca_beam_pip_Sm", "",    1000, 0.0, 5.0));
  dir-> Append(new TH1F("dca_beam_pip_Sp", "",    1000, 0.0, 5.0));
  dir-> Append(new TH1F("dca_beam_pip_woAll", "", 1000, 0.0, 5.0));

  dir-> Append(new TH1F("KNpim_MM", "",       1000, 1.0, 2.0));
  dir-> Append(new TH1F("KNpim_MM_K0", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KNpim_MM_Sm", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KNpim_MM_Sp", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KNpim_MM_woAll", "", 1000, 1.0, 2.0));

  dir-> Append(new TH1F("KNpip_MM", "",       1000, 1.0, 2.0));
  dir-> Append(new TH1F("KNpip_MM_K0", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KNpip_MM_Sm", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KNpip_MM_Sp", "",    1000, 1.0, 2.0));
  dir-> Append(new TH1F("KNpip_MM_woAll", "", 1000, 1.0, 2.0));

  dir-> Append(new TH1F("IM_pipi", "", 200, 0.0, 1.0));
  dir-> Append(new TH1F("IM_npim", "", 200, 1.0, 2.0));
  dir-> Append(new TH1F("IM_npip", "", 200, 1.0, 2.0));

  dir-> Append(new TH1F("IM_pipi_wK0", "", 200, 0.0, 1.0));
  dir-> Append(new TH1F("IM_npim_wSm", "", 200, 1.0, 2.0));
  dir-> Append(new TH1F("IM_npip_wSp", "", 200, 1.0, 2.0));

  dir-> Append(new TH1F("IM_pipi_1MeV", "", 1000, 0.0, 1.0));
  dir-> Append(new TH1F("IM_npim_1MeV", "", 1000, 1.0, 2.0));
  dir-> Append(new TH1F("IM_npip_1MeV", "", 1000, 1.0, 2.0));

  dir-> Append(new TH1F("IM_pipi_wK0_1MeV", "", 1000, 0.0, 1.0));
  dir-> Append(new TH1F("IM_npim_wSm_1MeV", "", 1000, 1.0, 2.0));
  dir-> Append(new TH1F("IM_npip_wSp_1MeV", "", 1000, 1.0, 2.0));

  dir-> Append(new TH1F("IM_pipi_woK0", "", 200, 0.0, 1.0));
  dir-> Append(new TH1F("IM_npim_woSm", "", 200, 1.0, 2.0));
  dir-> Append(new TH1F("IM_npip_woSp", "", 200, 1.0, 2.0));

  dir-> Append(new TH1F("IM_npipi",    "", 2000, 0, 2.0));
  dir-> Append(new TH1F("IM_npipi_K0", "", 2000, 0, 2.0));
  dir-> Append(new TH1F("IM_npipi_Sp", "", 2000, 0, 2.0));
  dir-> Append(new TH1F("IM_npipi_Sm", "", 2000, 0, 2.0));

  dir-> Append(new TH1F("pim_mom",    "", 1000, 0, 1.0));
  dir-> Append(new TH1F("pim_mom_K0", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("pim_mom_Sm", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("pim_mom_Sp", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("pim_mom_woAll", "", 1000, 0, 1.0));

  dir-> Append(new TH1F("pip_mom",    "", 1000, 0, 1.0));
  dir-> Append(new TH1F("pip_mom_K0", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("pip_mom_Sm", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("pip_mom_Sp", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("pip_mom_woAll", "", 1000, 0, 1.0));

  dir-> Append(new TH1F("mmN_mom",    "", 1000, 0, 1.0));
  dir-> Append(new TH1F("mmN_mom_K0", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("mmN_mom_Sm", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("mmN_mom_Sp", "", 1000, 0, 1.0));
  dir-> Append(new TH1F("mmN_mom_woAll", "", 1000, 0, 1.0));

  dir-> Append(new TH2F("pipi_cosOP_IM_pipi", "", 1000, -1.0, 1.0, 1000, 0, 1.0));
  dir-> Append(new TH2F("pim_mom_IM_pipi", "",    1000, 0, 1.0, 1000, 0, 1.0));
  dir-> Append(new TH2F("pip_mom_IM_pipi", "",    1000, 0, 1.0, 1000, 0, 1.0));
  dir-> Append(new TH2F("mmN_mom_IM_pipi", "",    1000, 0, 1.0, 1000, 0, 1.0));
  dir-> Append(new TH2F("mmN_mom_IM_npipi", "",   1000, 0, 1.0, 1000, 1.0, 2.0));
  dir-> Append(new TH2F("mmN_mom_IM_npipi_K0", "",   1000, 0, 1.0, 1000, 1.0, 2.0));
  dir-> Append(new TH2F("KN_MM_mmN_mom_K0", "",   1000, 1.0, 2.0, 1000, 0.0, 1.0));
  dir-> Append(new TH2F("K0_mom_IM_pipi", "",     1000, 0, 1.0, 1000, 0, 1.0));

  dir-> Append(new TH1F("K0_mom",    "", 1000, 0, 1.0));
  dir-> Append(new TH1F("K0_cos",    "", 1000, -1.0, 1.0));
  dir-> Append(new TH2F("mmN_mom_K0_cos", "", 1000, 0, 1.0, 1000, -1.0, 1.0));

  dir-> Append(new TH2F("mmN_mom_vtxR_K0", "",      1000, 0, 1.0, 1000,  0.0, 10.0));
  dir-> Append(new TH2F("mmN_mom_vtxX_K0", "",      1000, 0, 1.0, 1000, -10.0, 10.0));
  dir-> Append(new TH2F("mmN_mom_vtxY_K0", "",      1000, 0, 1.0, 1000, -10.0, 10.0));
  dir-> Append(new TH2F("mmN_mom_vtxZ_K0", "",      1000, 0, 1.0, 1000, -20.0,  0.0));
  dir-> Append(new TH2F("mmN_mom_KN_MM_K0", "",      1000, 0, 1.0, 1000, 1.0, 2.0));
  dir-> Append(new TH2F("mmN_mom_KN_MM_tarP_K0", "", 1000, 0, 1.0, 1000, 0, 1.0));
  dir-> Append(new TH2F("mmN_mom_fn_mom_K0", "", 1000, 0, 1.0, 1000, 0, 1.0));

  dir-> Append(new TH2F("KN_MM_IM_pipi", "",      1000, 1.0, 2.0, 1000,  0.0, 1.0));
  dir-> Append(new TH2F("KN_MM_IM_npim", "",      1000, 1.0, 2.0, 1000,  1.0, 2.0));
  dir-> Append(new TH2F("KN_MM_IM_npip", "",      1000, 1.0, 2.0, 1000,  1.0, 2.0));
  
  dir-> Append(new TH1F("N_mom",       "", 1000, 0.0, 1.0));
  dir-> Append(new TH1F("N_mom_K0",    "", 1000, 0.0, 1.0));
  dir-> Append(new TH1F("N_mom_Sm",    "", 1000, 0.0, 1.0));
  dir-> Append(new TH1F("N_mom_Sp",    "", 1000, 0.0, 1.0));
  dir-> Append(new TH1F("N_mom_woAll", "", 1000, 0.0, 1.0));

  TTree *tree = (TTree*)f->Get("npipi_event");
  AnaInfo *anaInfo = new AnaInfo();
  tree-> SetBranchAddress("AnaInfo", &anaInfo);

  for( int ev=0; ev<tree->GetEntries(); ev++ ){
    tree-> GetEntry(ev);
    TLorentzVector target_lmom;
    TVector3 vtxBeam=anaInfo->minDCA()->vertexBeam(), vtxCDS=anaInfo->minDCA()->vertexCDS();
    double vtxR=sqrt(vtxBeam.X()*vtxBeam.X()+vtxBeam.Y()*vtxBeam.Y());
    target_lmom.SetVectM(TVector3(0, 0, 0), dMass);
    CDSInfo* pim=anaInfo->CDS(CDS_PiMinus, 0);
    CDSInfo* pip=anaInfo->CDS(CDS_PiPlus, 0);
    CDS2Info* pipi=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0);
    TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
    TVector3 beam_mom=beam_lmom.Vect();
    TVector3 labToCM=-(beam_lmom+target_lmom).BoostVector();
    TLorentzVector beam_lmom_CM=beam_lmom; beam_lmom_CM.Boost(labToCM);
    TVector3 beam_mom_CM=beam_lmom_CM.Vect();

    TLorentzVector n_lmom=anaInfo->forwardNeutral(0)->lmom();
    TLorentzVector pim_lmom=anaInfo->CDS(CDS_PiMinus, 0)->lmom();
    TLorentzVector pip_lmom=anaInfo->CDS(CDS_PiPlus, 0)->lmom();
    TLorentzVector pipi_lmom=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0)->lmom();
    TLorentzVector pim_lmom_CM=pim_lmom, pip_lmom_CM=pip_lmom;
    pim_lmom_CM.Boost(labToCM), pip_lmom_CM.Boost(labToCM);

    double pipi_cosOP=MyTools::cosOP(pipi);
    double kn_mm=(beam_lmom+target_lmom-n_lmom).M();
    double kn_mm_tarP=(beam_lmom+TLorentzVector(TVector3(0, 0, 0), pMass)-n_lmom).M();
    double knpipi_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom-pip_lmom).M();
    double knpim_mm=(beam_lmom+target_lmom-n_lmom-pim_lmom).M();
    double knpip_mm=(beam_lmom+target_lmom-n_lmom-pip_lmom).M();
    double npim_im=(n_lmom+pim_lmom).M();
    double npip_im=(n_lmom+pip_lmom).M();
    double pipi_im=pipi_lmom.M();
    double npipi_im=(n_lmom+pim_lmom+pip_lmom).M();
    const int bin_index=fInfoSim->index(kn_mm);
    const double mmN_mom=(beam_lmom+target_lmom-n_lmom-pim_lmom-pip_lmom).Vect().Mag();
    bool K0_flag=isK0(pipi_im); bool Sm_flag=isIM_Sm(npim_im); bool Sp_flag=isIM_Sp(npip_im);

    fillHist(dir, "KNpipi_MM", knpipi_mm, bin_index, index);
    if( Npipi_N_MIN<knpipi_mm && knpipi_mm<Npipi_N_MAX ){
      fillHist(dir, "KN_MM_IM_pipi", kn_mm, pipi_im);
      fillHist(dir, "KN_MM_IM_npim", kn_mm, npim_im);
      fillHist(dir, "KN_MM_IM_npip", kn_mm, npip_im);

      fillHist(dir, "pim_mom_IM_pipi", pim->momentum().Mag(), pipi_im, bin_index, index);
      fillHist(dir, "pip_mom_IM_pipi", pip->momentum().Mag(), pipi_im, bin_index, index);
      fillHist(dir, "mmN_mom_IM_pipi", mmN_mom,    pipi_im, bin_index, index);
      fillHist(dir, "mmN_mom_IM_npipi", mmN_mom,  npipi_im, bin_index, index);
      fillHist(dir, "K0_mom_IM_pipi", pipi->lmom().Vect().Mag(), pipi_im, bin_index, index);

      fillHist(dir, "pim_mom", pim->momentum().Mag(), bin_index, index);
      fillHist(dir, "pip_mom", pip->momentum().Mag(), bin_index, index);
      
      fillHist(dir, "KN_MM", kn_mm, bin_index, index);
      fillHist(dir, "IM_pipi", pipi_im, bin_index, index);
      fillHist(dir, "IM_npim", npim_im, bin_index, index);
      fillHist(dir, "IM_npip", npip_im, bin_index, index);
      fillHist(dir, "IM_pipi_1MeV", pipi_im, bin_index, index);
      fillHist(dir, "IM_npim_1MeV", npim_im, bin_index, index);
      fillHist(dir, "IM_npip_1MeV", npip_im, bin_index, index);

      fillHist(dir, "IM_npipi", npipi_im, bin_index, index);
      fillHist(dir, "mmN_mom", mmN_mom, bin_index, index);
      fillHist(dir, "KNpim_MM", knpim_mm, bin_index, index);
      fillHist(dir, "KNpip_MM", knpip_mm, bin_index, index);
      fillHist(dir, "pipi_cosOP_IM_pipi", pipi_cosOP, pipi_im, bin_index, index);

      fillHist(dir, "dca_beam_pim", pim->dca(), bin_index, index);
      fillHist(dir, "dca_beam_pip", pip->dca(), bin_index, index);

      fillHist(dir, "N_mom", n_lmom.Vect().Mag(), bin_index, index);
      if( bin_index>=0 ){
	fillHist(dir, Form("dca_beam_pim_%d", bin_index), pim->dca(), bin_index, index);
	fillHist(dir, Form("dca_beam_pip_%d", bin_index), pip->dca(), bin_index, index);

	fillHist(dir, Form("KNpim_KNpip_MM_%d", bin_index), knpim_mm, knpip_mm, bin_index, index);
	fillHist(dir, Form("KNpim_MM_%d", bin_index), knpim_mm, bin_index, index);
	fillHist(dir, Form("KNpip_MM_%d", bin_index), knpip_mm, bin_index, index);
      }

      if( isK0(pipi_im) ){
	double K0cos=beam_lmom.Vect().Dot(pipi_lmom.Vect())/(beam_lmom.Vect().Mag()*pipi_lmom.Vect().Mag());
	fillHist(dir, "KN_MM_IM_npipi_K0", kn_mm, npipi_im, bin_index, index);
	fillHist(dir, "mmN_mom_IM_npipi_K0", mmN_mom,  npipi_im, bin_index, index);

	fillHist(dir, "dca_pipi_K0", pim->dca(), bin_index, index);

	fillHist(dir, "pipi_cosOP_K0", pipi_cosOP, bin_index, index);
	fillHist(dir, "pim_mom_K0", pim->momentum().Mag(), bin_index, index);
	fillHist(dir, "pip_mom_K0", pip->momentum().Mag(), bin_index, index);
	
	fillHist(dir, "KN_MM_K0", kn_mm, bin_index, index);
	fillHist(dir, "IM_npipi_K0", npipi_im, bin_index, index);
	fillHist(dir, "mmN_mom_K0", mmN_mom, bin_index, index);
	fillHist(dir, "K0_mom", pipi->lmom().Vect().Mag(), bin_index, index);
	fillHist(dir, "K0_cos", K0cos, bin_index, index);
	fillHist(dir, "mmN_mom_K0_cos", mmN_mom, K0cos, bin_index, index);
	fillHist(dir, "KNpim_MM_K0", knpim_mm, bin_index, index);
	fillHist(dir, "KNpip_MM_K0", knpip_mm, bin_index, index);

	fillHist(dir, "dca_beam_pim_K0", pim->dca(), bin_index, index);
	fillHist(dir, "dca_beam_pip_K0", pip->dca(), bin_index, index);

	fillHist(dir, "mmN_mom_KN_MM_K0", mmN_mom, kn_mm, bin_index, index);
	fillHist(dir, "mmN_mom_fn_mom_K0", mmN_mom, n_lmom.Vect().Mag(), bin_index, index);
	fillHist(dir, "mmN_mom_KN_MM_tarP_K0", mmN_mom, kn_mm_tarP, bin_index, index);
	fillHist(dir, "mmN_mom_vtxR_K0", mmN_mom, vtxR, bin_index, index);
	fillHist(dir, "mmN_mom_vtxX_K0", mmN_mom, vtxBeam.X(), bin_index, index);
	fillHist(dir, "mmN_mom_vtxY_K0", mmN_mom, vtxBeam.Y(), bin_index, index);
	fillHist(dir, "mmN_mom_vtxZ_K0", mmN_mom, vtxBeam.Z(), bin_index, index);

	fillHist(dir, "N_mom_K0", n_lmom.Vect().Mag(), bin_index, index);

	if( bin_index>=0 ){
	  fillHist(dir, Form("dca_beam_pim_K0_%d", bin_index), pim->dca(), bin_index, index);
	  fillHist(dir, Form("dca_beam_pip_K0_%d", bin_index), pip->dca(), bin_index, index);
	  
	  fillHist(dir, Form("KNpim_KNpip_MM_K0_%d", bin_index), knpim_mm, knpip_mm, bin_index, index);
	  fillHist(dir, Form("KNpim_MM_K0_%d", bin_index), knpim_mm, bin_index, index);
	  fillHist(dir, Form("KNpip_MM_K0_%d", bin_index), knpip_mm, bin_index, index);
	}
	fillHist(dir, "IM_pipi_wK0", pipi_im, bin_index, index);
	fillHist(dir, "IM_pipi_wK0_1MeV", pipi_im, bin_index, index);
      }
      else fillHist(dir, "IM_pipi_woK0", pipi_im, bin_index, index);

      if( isIM_Sm(npim_im) ){
	fillHist(dir, "KN_MM_IM_npipi_Sm", kn_mm, npipi_im, bin_index, index);
	
	fillHist(dir, "pim_mom_Sm", pim->momentum().Mag(), bin_index, index);
	fillHist(dir, "pip_mom_Sm", pip->momentum().Mag(), bin_index, index);
	
	fillHist(dir, "KN_MM_Sm", kn_mm, bin_index, index);
	fillHist(dir, "IM_npipi_Sm", npipi_im, bin_index, index);
	fillHist(dir, "mmN_mom_Sm", mmN_mom, bin_index, index);
	fillHist(dir, "KNpim_MM_Sm", knpim_mm, bin_index, index);
	fillHist(dir, "KNpip_MM_Sm", knpip_mm, bin_index, index);
	
	fillHist(dir, "dca_beam_pim_Sm", pim->dca(), bin_index, index);
	fillHist(dir, "dca_beam_pip_Sm", pip->dca(), bin_index, index);

	fillHist(dir, "N_mom_Sm", n_lmom.Vect().Mag(), bin_index, index);

	if( bin_index>=0 ){
	  fillHist(dir, Form("dca_beam_pim_Sm_%d", bin_index), pim->dca(), bin_index, index);
	  fillHist(dir, Form("dca_beam_pip_Sm_%d", bin_index), pip->dca(), bin_index, index);
	  
	  fillHist(dir, Form("KNpim_KNpip_MM_Sm_%d", bin_index), knpim_mm, knpip_mm, bin_index, index);
	  fillHist(dir, Form("KNpim_MM_Sm_%d", bin_index), knpim_mm, bin_index, index);
	  fillHist(dir, Form("KNpip_MM_Sm_%d", bin_index), knpip_mm, bin_index, index);
	}
	fillHist(dir, "IM_npim_wSm", npim_im, bin_index, index);
	fillHist(dir, "IM_npim_wSm_1MeV", npim_im, bin_index, index);
      }
      else fillHist(dir, "IM_npim_woSm", npim_im, bin_index, index);

      if( isIM_Sp(npip_im) ){
	fillHist(dir, "KN_MM_IM_npipi_Sp", kn_mm, npipi_im, bin_index, index);
	
	fillHist(dir, "pim_mom_Sp", pim->momentum().Mag(), bin_index, index);
	fillHist(dir, "pip_mom_Sp", pip->momentum().Mag(), bin_index, index);
	
	fillHist(dir, "KN_MM_Sp", kn_mm, bin_index, index);
	fillHist(dir, "IM_npip_woSp", npip_im, bin_index, index);
	fillHist(dir, "IM_npipi_Sp", npipi_im, bin_index, index);
	fillHist(dir, "mmN_mom_Sp", mmN_mom, bin_index, index);
	fillHist(dir, "KNpim_MM_Sp", knpim_mm, bin_index, index);
	fillHist(dir, "KNpip_MM_Sp", knpip_mm, bin_index, index);
	
	fillHist(dir, "dca_beam_pim_Sp", pim->dca(), bin_index, index);
	fillHist(dir, "dca_beam_pip_Sp", pip->dca(), bin_index, index);

	fillHist(dir, "N_mom_Sp", n_lmom.Vect().Mag(), bin_index, index);

	if( bin_index>=0 ){
	  fillHist(dir, Form("dca_beam_pim_Sp_%d", bin_index), pim->dca(), bin_index, index);
	  fillHist(dir, Form("dca_beam_pip_Sp_%d", bin_index), pip->dca(), bin_index, index);

	  fillHist(dir, Form("KNpim_KNpip_MM_Sp_%d", bin_index), knpim_mm, knpip_mm, bin_index, index);
	  fillHist(dir, Form("KNpim_MM_Sp_%d", bin_index), knpim_mm, bin_index, index);
	  fillHist(dir, Form("KNpip_MM_Sp_%d", bin_index), knpip_mm, bin_index, index);
	}
	fillHist(dir, "IM_npip_wSp", npip_im, bin_index, index);
	fillHist(dir, "IM_npip_wSp_1MeV", npip_im, bin_index, index);
      }
      else fillHist(dir, "IM_npip_woSp", npip_im, bin_index, index);
      		 
      if( isSignal(pipi_im, npim_im, npip_im) ){
	fillHist(dir, "KN_MM_woAll", kn_mm, bin_index, index);
	
	fillHist(dir, "pim_mom_woAll", pim->momentum().Mag(), bin_index, index);
	fillHist(dir, "pip_mom_woAll", pip->momentum().Mag(), bin_index, index);
	
	fillHist(dir, "dca_beam_pim_woAll", pim->dca(), bin_index, index);
	fillHist(dir, "dca_beam_pip_woAll", pip->dca(), bin_index, index);

	fillHist(dir, "N_mom_woAll", n_lmom.Vect().Mag(), bin_index, index);

	if( bin_index>=0 ){
	  fillHist(dir, Form("dca_beam_pim_woAll_%d", bin_index), pim->dca(), bin_index, index);
	  fillHist(dir, Form("dca_beam_pip_woAll_%d", bin_index), pip->dca(), bin_index, index);

	  fillHist(dir, "mmN_mom_woAll", mmN_mom, bin_index, index);

	  fillHist(dir, "KNpim_MM_woAll", knpim_mm, bin_index, index);
	  fillHist(dir, "KNpip_MM_woAll", knpip_mm, bin_index, index);
	  fillHist(dir, Form("KNpim_KNpip_MM_woAll_%d", bin_index), knpim_mm, knpip_mm, bin_index, index);
	  fillHist(dir, Form("KNpim_MM_woAll_%d", bin_index), knpim_mm, bin_index, index);
	  fillHist(dir, Form("KNpip_MM_woAll_%d", bin_index), knpip_mm, bin_index, index);	
	}
      }
    }
  }

  dir-> Write();
}
