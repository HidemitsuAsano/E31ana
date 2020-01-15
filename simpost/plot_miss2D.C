//plot_miss2D
//purpose : plot Missing mass, IM(npip) and so on. BG fit by modifying MC data.

void plot_miss2D()
{
  TH1::SetDefaultSumw2();
                                //rdata,nSp,nSm,Lambda,S0pi+pi-,nSppi0,nSmpi0,mcsum,K0nn
  const unsigned int colordef[9]={1,2,3,4,5,7,9,28,6};
  const char name[][10]={"rdata","Sigmap","Sigman","Lambda","S0","Sppi0","Smpi0","K0nn"};
  
  //real data
  TFile *filerdata = TFile::Open("../post/evanaIMpisigma_npippim_v196_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double  nrdata_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_rdata->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nrdata_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_rdata->Integral();
  
  //Lambda pi+ pi- n sim.
  TFile *fileLambda = TFile::Open("simIMpisigma_npipiL_pippimn_v23_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda = (TH2D*) fileLambda->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nLambda_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda = (TH2D*) fileLambda->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nLambda_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda->Integral();
  

  //Sigma+ mode sim. ,flat dist in q vs IMnpi+pi-
  TFile *filenSp = TFile::Open("simIMpisigma_nSppim_pippimn_v108_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap = (TH2D*) filenSp->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSigmap_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap = (TH2D*) filenSp->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSigmap_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap->Integral();


  //Sigma- mode sim. flat. dist in q vs IMnpi+pi-
  TFile *filenSm = TFile::Open("simIMpisigma_nSmpip_pippimn_v108_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam = (TH2D*) filenSm->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSigmam_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam = (TH2D*) filenSm->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSigmam_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam->Integral();
  

  //nS0pi+pi-
  TFile *fileS0pipi = TFile::Open("simIMpisigma_nS0pippim_pippimn_v5_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_S0 = (TH2D*) fileS0pipi->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nS0miss_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_S0->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_S0 = (TH2D*) fileS0pipi->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nS0miss_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_S0->Integral();

  
  //nSppi-pi0
  TFile *fileSppi0 = TFile::Open("simIMpisigma_nSppimpi0_pippimn_v4_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Sppi0 = (TH2D*) fileSppi0->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSppi0_Sp= MMnmiss_IMnpipi_woK0_wSid_Sp_Sppi0->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Sppi0 = (TH2D*) fileSppi0->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSppi0_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Sppi0->Integral();
  
  //nSmpi+pi0
  TFile *fileSmpi0 = TFile::Open("simIMpisigma_nSmpippi0_pippimn_v4_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Smpi0 = (TH2D*) fileSmpi0->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSmpi0_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Smpi0->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Smpi0 = (TH2D*) fileSmpi0->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSmpi0_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Smpi0->Integral();
  
  //K0nn
  TFile *fileK0nn = TFile::Open("simIMpisigma_K0nn_pippimn_v8_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_K0nn = (TH2D*) fileK0nn->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double K0nn_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_K0nn->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_K0nn = (TH2D*) fileK0nn->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double K0nn_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_K0nn->Integral();

  

  //superimpose 2d hists of miss mass vs IM(npi+) or IM(npi-)
  //Lambda
  TH2D *MMnmiss_IMnpip_woK0_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_Lambda = (TH2D*)fileLambda->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_Lambda  = (TH2D*)fileLambda->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_Lambda = (TH2D*)fileLambda->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  double Nnpip_Lambda = MMnmiss_IMnpip_woK0_Lambda->Integral();
  double Nnpim_Lambda = MMnmiss_IMnpim_woK0_Lambda->Integral();
  //nSppi-
  TH2D *MMnmiss_IMnpip_woK0_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_Sigmap = (TH2D*)filenSp->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_Sigmap  = (TH2D*)filenSp->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_Sigmap = (TH2D*)filenSp->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  double Nnpip_Sigmap = MMnmiss_IMnpip_woK0_Sigmap->Integral();
  double Nnpim_Sigmap = MMnmiss_IMnpim_woK0_Sigmap->Integral();
  //nSmpi+
  TH2D *MMnmiss_IMnpip_woK0_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_Sigmam = (TH2D*)filenSm->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_Sigmam  = (TH2D*)filenSm->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_Sigmam = (TH2D*)filenSm->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  double Nnpip_Sigmam = MMnmiss_IMnpip_woK0_Sigmam->Integral();
  double Nnpim_Sigmam = MMnmiss_IMnpim_woK0_Sigmam->Integral();
  //Sigma0
  TH2D *MMnmiss_IMnpip_woK0_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_S0 = (TH2D*)fileS0pipi->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_S0  = (TH2D*)fileS0pipi->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  double Nnpip_S0 = MMnmiss_IMnpip_woK0_S0->Integral();
  double Nnpim_S0 = MMnmiss_IMnpim_woK0_S0->Integral();
  //Sppi0
  TH2D *MMnmiss_IMnpip_woK0_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_Sppi0 = (TH2D*)fileSppi0->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_Sppi0 = (TH2D*)fileSppi0->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_Sppi0  = (TH2D*)fileSppi0->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  double Nnpip_Sppi0 = MMnmiss_IMnpip_woK0_Sppi0->Integral();
  double Nnpim_Sppi0 = MMnmiss_IMnpim_woK0_Sppi0->Integral();
  //Smpi0
  TH2D *MMnmiss_IMnpip_woK0_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_Smpi0 = (TH2D*)fileSmpi0->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_Smpi0 = (TH2D*)fileSmpi0->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_Smpi0  = (TH2D*)fileSmpi0->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  double Nnpip_Smpi0 = MMnmiss_IMnpip_woK0_Smpi0->Integral();
  double Nnpim_Smpi0 = MMnmiss_IMnpim_woK0_Smpi0->Integral();
  //K0nn
  TH2D *MMnmiss_IMnpip_woK0_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_K0nn = (TH2D*)fileK0nn->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_K0nn = (TH2D*)fileK0nn->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_K0nn  = (TH2D*)fileK0nn->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_K0nn = (TH2D*)fileK0nn->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  double Nnpip_K0nn = MMnmiss_IMnpip_woK0_K0nn->Integral();
  double Nnpim_K0nn = MMnmiss_IMnpim_woK0_K0nn->Integral();
  //real data
  TH2D *MMnmiss_IMnpip_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_woK0_woSm_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_woK0_woSp_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_wSid_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_wSid");
  TH2D *MMnmiss_IMpippim_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_woK0_woSidn_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  TH2D *q_IMnpipi_woK0_woSidn_rdata = (TH2D*)filerdata->Get("q_IMnpipi_woK0_woSidn");
  TH2D *q_IMnpipi_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("q_IMnpipi_wK0_woSid_won");
  TH2D *MMnmiss_Mompippim_woK0_woSidn_rdata  = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_woK0_woSidn");
  TH2D *MMnmiss_Mompippim_wK0_woSid_won_rdata = (TH2D*)filerdata->Get("MMnmiss_Mompippim_dE_wK0_woSid_won");
  double Nnpip_rdata = MMnmiss_IMnpip_woK0_rdata->Integral();
  double Nnpim_rdata = MMnmiss_IMnpim_woK0_rdata->Integral();
  /*
  // v1 Jan. 9th, 2020 
  const double scaleFactor_Sp[8]={1.0,
                                  0.18*nrdata_Sp/nSigmap_Sp,
                                  0.01*nrdata_Sp/nSigmam_Sp,
                                  0.38*nrdata_Sp/nLambda_Sp,//  
                                  0.18*nrdata_Sp/nS0miss_Sp,//  
                                  0.02*nrdata_Sp/nSppi0_Sp,
                                  0.02*nrdata_Sp/nSmpi0_Sp,
                                  0.015*nrdata_Sp/K0nn_Sp};
  
  const double scaleFactor_Sm[8]={1.0,
                                  0.01*nrdata_Sm/nSigmap_Sm,
                                  0.21*nrdata_Sm/nSigmam_Sm,
                                  0.38*nrdata_Sm/nLambda_Sm,
                                  0.18*nrdata_Sm/nS0miss_Sm,
                                  0.03*nrdata_Sm/nSppi0_Sm,
                                  0.03*nrdata_Sm/nSmpi0_Sm,
                                  0.015*nrdata_Sm/K0nn_Sm};
  */
  
  //v2 2020113
  const double scaleFactor_Sp[8]={1.0,
                                  0.18*nrdata_Sp/nSigmap_Sp,
                                  0.01*nrdata_Sp/nSigmam_Sp,
                                  0.38*nrdata_Sp/nLambda_Sp,//  
                                  0.08*nrdata_Sp/nS0miss_Sp,//  
                                  0.02*nrdata_Sp/nSppi0_Sp,
                                  0.02*nrdata_Sp/nSmpi0_Sp,
                                  0.015*nrdata_Sp/K0nn_Sp};
  
  const double scaleFactor_Sm[8]={1.0,
                                  0.01*nrdata_Sm/nSigmap_Sm,
                                  0.205*nrdata_Sm/nSigmam_Sm,
                                  0.43*nrdata_Sm/nLambda_Sm,
                                  0.07*nrdata_Sm/nS0miss_Sm,
                                  0.03*nrdata_Sm/nSppi0_Sm,
                                  0.03*nrdata_Sm/nSmpi0_Sm,
                                  0.015*nrdata_Sm/K0nn_Sm};
  
  
  //array 
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp[8]={
    MMnmiss_IMnpipi_woK0_wSid_Sp_rdata,
    MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap,
    MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam,
    MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda,
    MMnmiss_IMnpipi_woK0_wSid_Sp_S0,
    MMnmiss_IMnpipi_woK0_wSid_Sp_Sppi0,
    MMnmiss_IMnpipi_woK0_wSid_Sp_Smpi0,
    MMnmiss_IMnpipi_woK0_wSid_Sp_K0nn};
  
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm[8]={
    MMnmiss_IMnpipi_woK0_wSid_Sm_rdata,
    MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap,
    MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam,
    MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda,
    MMnmiss_IMnpipi_woK0_wSid_Sm_S0,
    MMnmiss_IMnpipi_woK0_wSid_Sm_Sppi0,
    MMnmiss_IMnpipi_woK0_wSid_Sm_Smpi0,
    MMnmiss_IMnpipi_woK0_wSid_Sm_K0nn};
  
  TH2D* MMnmiss_IMnpip_woK0[8]={
    MMnmiss_IMnpip_woK0_rdata,
    MMnmiss_IMnpip_woK0_Sigmap,
    MMnmiss_IMnpip_woK0_Sigmam,
    MMnmiss_IMnpip_woK0_Lambda,
    MMnmiss_IMnpip_woK0_S0,
    MMnmiss_IMnpip_woK0_Sppi0,
    MMnmiss_IMnpip_woK0_Smpi0, 
    MMnmiss_IMnpip_woK0_K0nn};

  TH2D* MMnmiss_IMnpip_woK0_woSm[8]={
    MMnmiss_IMnpip_woK0_woSm_rdata,
    MMnmiss_IMnpip_woK0_woSm_Sigmap,
    MMnmiss_IMnpip_woK0_woSm_Sigmam,
    MMnmiss_IMnpip_woK0_woSm_Lambda,
    MMnmiss_IMnpip_woK0_woSm_S0,
    MMnmiss_IMnpip_woK0_woSm_Sppi0,
    MMnmiss_IMnpip_woK0_woSm_Smpi0,
    MMnmiss_IMnpip_woK0_woSm_K0nn};
  
  TH2D* MMnmiss_IMnpip_woK0_woSm_n[8]={
    MMnmiss_IMnpip_woK0_woSm_n_rdata,
    MMnmiss_IMnpip_woK0_woSm_n_Sigmap,
    MMnmiss_IMnpip_woK0_woSm_n_Sigmam,
    MMnmiss_IMnpip_woK0_woSm_n_Lambda,
    MMnmiss_IMnpip_woK0_woSm_n_S0,
    MMnmiss_IMnpip_woK0_woSm_n_Sppi0,
    MMnmiss_IMnpip_woK0_woSm_n_Smpi0,
    MMnmiss_IMnpip_woK0_woSm_n_K0nn};

  TH2D* MMnmiss_IMnpip_woK0_woSidn[8]={
    MMnmiss_IMnpip_woK0_woSidn_rdata,
    MMnmiss_IMnpip_woK0_woSidn_Sigmap,
    MMnmiss_IMnpip_woK0_woSidn_Sigmam,
    MMnmiss_IMnpip_woK0_woSidn_Lambda,
    MMnmiss_IMnpip_woK0_woSidn_S0,
    MMnmiss_IMnpip_woK0_woSidn_Sppi0,
    MMnmiss_IMnpip_woK0_woSidn_Smpi0,
    MMnmiss_IMnpip_woK0_woSidn_K0nn};


  TH2D* MMnmiss_IMnpim_woK0[8]={
    MMnmiss_IMnpim_woK0_rdata,
    MMnmiss_IMnpim_woK0_Sigmap,
    MMnmiss_IMnpim_woK0_Sigmam,
    MMnmiss_IMnpim_woK0_Lambda,
    MMnmiss_IMnpim_woK0_S0,
    MMnmiss_IMnpim_woK0_Sppi0,
    MMnmiss_IMnpim_woK0_Smpi0,
    MMnmiss_IMnpim_woK0_K0nn};

  TH2D* MMnmiss_IMnpim_woK0_woSm[8]={
    MMnmiss_IMnpim_woK0_woSp_rdata,
    MMnmiss_IMnpim_woK0_woSp_Sigmap,
    MMnmiss_IMnpim_woK0_woSp_Sigmam,
    MMnmiss_IMnpim_woK0_woSp_Lambda,
    MMnmiss_IMnpim_woK0_woSp_S0,
    MMnmiss_IMnpim_woK0_woSp_Sppi0,
    MMnmiss_IMnpim_woK0_woSp_Smpi0,
    MMnmiss_IMnpim_woK0_woSp_K0nn};
  
  TH2D* MMnmiss_IMnpim_woK0_woSp[8]={
    MMnmiss_IMnpim_woK0_woSp_rdata,
    MMnmiss_IMnpim_woK0_woSp_Sigmap,
    MMnmiss_IMnpim_woK0_woSp_Sigmam,
    MMnmiss_IMnpim_woK0_woSp_Lambda,
    MMnmiss_IMnpim_woK0_woSp_S0,
    MMnmiss_IMnpim_woK0_woSp_Sppi0,
    MMnmiss_IMnpim_woK0_woSp_Smpi0,
    MMnmiss_IMnpim_woK0_woSp_K0nn};
  
  TH2D* MMnmiss_IMnpim_woK0_woSp_n[8]={
    MMnmiss_IMnpim_woK0_woSp_n_rdata,
    MMnmiss_IMnpim_woK0_woSp_n_Sigmap,
    MMnmiss_IMnpim_woK0_woSp_n_Sigmam,
    MMnmiss_IMnpim_woK0_woSp_n_Lambda,
    MMnmiss_IMnpim_woK0_woSp_n_S0,
    MMnmiss_IMnpim_woK0_woSp_n_Sppi0,
    MMnmiss_IMnpim_woK0_woSp_n_Smpi0,
    MMnmiss_IMnpim_woK0_woSp_n_K0nn};

  TH2D* MMnmiss_IMnpim_woK0_woSidn[8]={
    MMnmiss_IMnpim_woK0_woSidn_rdata,
    MMnmiss_IMnpim_woK0_woSidn_Sigmap,
    MMnmiss_IMnpim_woK0_woSidn_Sigmam,
    MMnmiss_IMnpim_woK0_woSidn_Lambda,
    MMnmiss_IMnpim_woK0_woSidn_S0,
    MMnmiss_IMnpim_woK0_woSidn_Sppi0,
    MMnmiss_IMnpim_woK0_woSidn_Smpi0,
    MMnmiss_IMnpim_woK0_woSidn_K0nn};

  TH2D* MMnmiss_IMpippim_wSid[8]={
    MMnmiss_IMpippim_wSid_rdata,
    MMnmiss_IMpippim_wSid_Sigmap,
    MMnmiss_IMpippim_wSid_Sigmam,
    MMnmiss_IMpippim_wSid_Lambda,
    MMnmiss_IMpippim_wSid_S0,
    MMnmiss_IMpippim_wSid_Sppi0,
    MMnmiss_IMpippim_wSid_Smpi0,
    MMnmiss_IMpippim_wSid_K0nn};

  TH2D* MMnmiss_IMpippim_woK0_woSidn[8]={
    MMnmiss_IMpippim_woK0_woSidn_rdata,
    MMnmiss_IMpippim_woK0_woSidn_Sigmap,
    MMnmiss_IMpippim_woK0_woSidn_Sigmam,
    MMnmiss_IMpippim_woK0_woSidn_Lambda,
    MMnmiss_IMpippim_woK0_woSidn_S0,
    MMnmiss_IMpippim_woK0_woSidn_Sppi0,
    MMnmiss_IMpippim_woK0_woSidn_Smpi0,
    MMnmiss_IMpippim_woK0_woSidn_K0nn};
  
  TH2D* q_IMnpipi_woK0_woSidn[8]={
    q_IMnpipi_woK0_woSidn_rdata,
    q_IMnpipi_woK0_woSidn_Sigmap,
    q_IMnpipi_woK0_woSidn_Sigmam,
    q_IMnpipi_woK0_woSidn_Lambda,
    q_IMnpipi_woK0_woSidn_S0,
    q_IMnpipi_woK0_woSidn_Sppi0,
    q_IMnpipi_woK0_woSidn_Smpi0,
    q_IMnpipi_woK0_woSidn_K0nn};

  TH2D* q_IMnpipi_wK0_woSid_won[8]={
    q_IMnpipi_wK0_woSid_won_rdata,
    q_IMnpipi_wK0_woSid_won_Sigmap,
    q_IMnpipi_wK0_woSid_won_Sigmam,
    q_IMnpipi_wK0_woSid_won_Lambda,
    q_IMnpipi_wK0_woSid_won_S0,
    q_IMnpipi_wK0_woSid_won_Sppi0,
    q_IMnpipi_wK0_woSid_won_Smpi0,
    q_IMnpipi_wK0_woSid_won_K0nn};

  TH2D* MMnmiss_Mompippim_woK0_woSidn[8]={
    MMnmiss_Mompippim_woK0_woSidn_rdata,
    MMnmiss_Mompippim_woK0_woSidn_Sigmap,
    MMnmiss_Mompippim_woK0_woSidn_Sigmam,
    MMnmiss_Mompippim_woK0_woSidn_Lambda,
    MMnmiss_Mompippim_woK0_woSidn_S0,
    MMnmiss_Mompippim_woK0_woSidn_Sppi0,
    MMnmiss_Mompippim_woK0_woSidn_Smpi0,
    MMnmiss_Mompippim_woK0_woSidn_K0nn};
  
  TH2D* MMnmiss_Mompippim_wK0_woSid_won[8]={
    MMnmiss_Mompippim_wK0_woSid_won_rdata,
    MMnmiss_Mompippim_wK0_woSid_won_Sigmap,
    MMnmiss_Mompippim_wK0_woSid_won_Sigmam,
    MMnmiss_Mompippim_wK0_woSid_won_Lambda,
    MMnmiss_Mompippim_wK0_woSid_won_S0,
    MMnmiss_Mompippim_wK0_woSid_won_Sppi0,
    MMnmiss_Mompippim_wK0_woSid_won_Smpi0,
    MMnmiss_Mompippim_wK0_woSid_won_K0nn};

  //Sp
  for(int i=0;i<8;i++){
    MMnmiss_IMnpipi_woK0_wSid_Sp[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_woK0[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_woK0_woSm[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_woK0_woSm_n[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_woK0_woSidn[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMpippim_wSid[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMpippim_woK0_woSidn[i]->Scale(scaleFactor_Sp[i]);
    q_IMnpipi_woK0_woSidn[i]->Scale(scaleFactor_Sp[i]);
    q_IMnpipi_wK0_woSid_won[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_Mompippim_woK0_woSidn[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_Mompippim_wK0_woSid_won[i]->Scale(scaleFactor_Sp[i]);
  }
  
  //Sm
  for(int i=0;i<8;i++){
    MMnmiss_IMnpipi_woK0_wSid_Sm[i]->Scale(scaleFactor_Sm[i]);
    MMnmiss_IMnpim_woK0[i]->Scale(scaleFactor_Sm[i]);
    MMnmiss_IMnpim_woK0_woSp[i]->Scale(scaleFactor_Sm[i]);
    MMnmiss_IMnpim_woK0_woSp_n[i]->Scale(scaleFactor_Sm[i]);
    MMnmiss_IMnpim_woK0_woSidn[i]->Scale(scaleFactor_Sm[i]);
  }
  
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_mc = (TH2D*)MMnmiss_IMnpipi_woK0_wSid_Sp[1]->Clone("MMnmiss_IMnpipi_woK0_wSid_Sp_mc");
  for(int i=2;i<8;i++) MMnmiss_IMnpipi_woK0_wSid_Sp_mc->Add(MMnmiss_IMnpipi_woK0_wSid_Sp[i]);
  
  TPaveText *pt1 = new TPaveText(.82,0.90,0.98,0.99,"NDC");
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_Sp = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_Sp","cMMnmiss_IMnpipi_woK0_wSid_Sp",1200,800);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->Divide(2,1);
  cMMnmiss_IMnpipi_woK0_wSid_Sp->cd(1);
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sp}{  real Data}");
  MMnmiss_IMnpipi_woK0_wSid_Sp[0]->Draw("colz");
  cMMnmiss_IMnpipi_woK0_wSid_Sp->cd(2);
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sp}{  MC sum}");
  MMnmiss_IMnpipi_woK0_wSid_Sp_mc->Draw("colz");

  TCanvas *cMMnmiss_woK0_wSid_Sp = new TCanvas("cMMnmiss_woK0_wSid_Sp","cMMnmiss_woK0_wSid_Sp",800,800);
  cMMnmiss_woK0_wSid_Sp->cd();
  TH1D* MMnmiss_woK0_wSid_Sp[8];
  for(int i=0;i<8;i++) MMnmiss_woK0_wSid_Sp[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp[i]->ProjectionX(Form("MMnmiss_woK0_wSid_Sp_%s",name[i]));
  MMnmiss_woK0_wSid_Sp[0]->Draw("HE");
  TH1D* MMnmiss_woK0_wSid_Sp_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_mc->ProjectionX("MMnmiss_woK0_wSid_Sp_mc");
  MMnmiss_woK0_wSid_Sp_mc->SetLineColor(6);
  MMnmiss_woK0_wSid_Sp_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    MMnmiss_woK0_wSid_Sp[i]->SetLineColor(colordef[i]);
    MMnmiss_woK0_wSid_Sp[i]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_woK0_wSid_Sp  = new TCanvas("cIMnpipi_woK0_wSid_Sp","cIMnpipi_woK0_wSid_Sp");
  cIMnpipi_woK0_wSid_Sp->cd();
  TH1D* IMnpipi_woK0_wSid_Sp[8];
  for(int i=0;i<8;i++)IMnpipi_woK0_wSid_Sp[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp[i]->ProjectionY(Form("IMnpipi_woK0_wSid_Sp_%s",name[i]));
  IMnpipi_woK0_wSid_Sp[0]->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_Sp_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_mc->ProjectionY("IMnpipi_woK0_wSid_Sp_mc");
  IMnpipi_woK0_wSid_Sp_mc->SetLineColor(6);
  IMnpipi_woK0_wSid_Sp_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMnpipi_woK0_wSid_Sp[i]->SetLineColor(colordef[i]);
    IMnpipi_woK0_wSid_Sp[i]->Draw("HEsame");
  }


  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_mc = (TH2D*)MMnmiss_IMnpipi_woK0_wSid_Sm[1]->Clone("MMnmiss_IMnpipi_woK0_wSid_Sm_mc");
  for(int i=2;i<8;i++) MMnmiss_IMnpipi_woK0_wSid_Sm_mc->Add(MMnmiss_IMnpipi_woK0_wSid_Sm[i]);
  
  TCanvas *cMMnmiss_IMnpipi_woK0_wSid_Sm = new TCanvas("cMMnmiss_IMnpipi_woK0_wSid_Sm","cMMnmiss_IMnpipi_woK0_wSid_Sm",1200,800);
  cMMnmiss_IMnpipi_woK0_wSid_Sm->Divide(2,1);
  cMMnmiss_IMnpipi_woK0_wSid_Sm->cd(1);
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sm}{  real data}");
  MMnmiss_IMnpipi_woK0_wSid_Sm[0]->Draw("colz");
  cMMnmiss_IMnpipi_woK0_wSid_Sm->cd(2);
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->SetTitle("#splitline{MMnmiss_IMnpipi_woK0_wSid_Sm}{  MC sum}");
  MMnmiss_IMnpipi_woK0_wSid_Sm_mc->Draw("colz");

  TCanvas *cMMnmiss_woK0_wSid_Sm = new TCanvas("cMMnmiss_woK0_wSid_Sm","cMMnmiss_woK0_wSid_Sm",800,800);
  cMMnmiss_woK0_wSid_Sm->cd();
  TH1D* MMnmiss_woK0_wSid_Sm[8];
  for(int i=0;i<8;i++) MMnmiss_woK0_wSid_Sm[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm[i]->ProjectionX(Form("MMnmiss_woK0_wSid_Sm_%s",name[i]));
  MMnmiss_woK0_wSid_Sm[0]->Draw("HE");
  TH1D* MMnmiss_woK0_wSid_Sm_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_mc->ProjectionX("MMnmiss_woK0_wSid_Sm_mc");
  MMnmiss_woK0_wSid_Sm_mc->SetLineColor(6);
  MMnmiss_woK0_wSid_Sm_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    MMnmiss_woK0_wSid_Sm[i]->SetLineColor(colordef[i]);
    MMnmiss_woK0_wSid_Sm[i]->Draw("HEsame");
  }

  TCanvas *cIMnpipi_woK0_wSid_Sm  = new TCanvas("cIMnpipi_woK0_wSid_Sm","cIMnpipi_woK0_wSid_Sm");
  cIMnpipi_woK0_wSid_Sm->cd();
  TH1D* IMnpipi_woK0_wSid_Sm[8];
  for(int i=0;i<8;i++)IMnpipi_woK0_wSid_Sm[i] = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm[i]->ProjectionY(Form("IMnpipi_woK0_wSid_Sm_%s",name[i]));
  IMnpipi_woK0_wSid_Sm[0]->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_Sm_mc = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_mc->ProjectionY("IMnpipi_woK0_wSid_Sm_mc");
  IMnpipi_woK0_wSid_Sm_mc->SetLineColor(6);
  IMnpipi_woK0_wSid_Sm_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMnpipi_woK0_wSid_Sm[i]->SetLineColor(colordef[i]);
    IMnpipi_woK0_wSid_Sm[i]->Draw("HEsame");
  }
  
  
  //missing mass and IMnpip/npim w/ Sp or Sm selection w/o K0
  //
  //ploting Sp mode
  //
  
  //w/o K0, including Sp/Sm mode
  TH2D *MMnmiss_IMnpip_woK0_mc = MMnmiss_IMnpip_woK0[1]->Clone("MMnmiss_IMnpip_woK0_mc");
  //adding all MC data
  for(int i=2;i<8;i++)  MMnmiss_IMnpip_woK0_mc->Add(MMnmiss_IMnpip_woK0[i]);
  
  TCanvas *cMMnmiss_IMnpip_woK0 = new TCanvas("cMMnmiss_IMnpip_woK0","cMMnmiss_IMnpip_woK0",1200,800);
  cMMnmiss_IMnpip_woK0->Divide(2,1);
  cMMnmiss_IMnpip_woK0->cd(1);
  MMnmiss_IMnpip_woK0_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0}{  Real data}");
  MMnmiss_IMnpip_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpip_woK0->cd(2);
  MMnmiss_IMnpip_woK0_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0}{  MC sum}");
  MMnmiss_IMnpip_woK0_mc->Draw("colz");
  
  //w/o K0, w/o Sm mode
  TH2D *MMnmiss_IMnpip_woK0_woSm_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSm[1]->Clone("MMnmiss_IMnpip_woK0_woSm_mc");
  //adding all MC data
  for(int i=2;i<8;i++)MMnmiss_IMnpip_woK0_woSm_mc->Add(MMnmiss_IMnpip_woK0_woSm[i]);
  
  //w/o K0, w/o Sm mode, selecting missing neutron
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSm_n[1]->Clone("MMnmiss_IMnpip_woK0_woSm_n_mc");
  for(int i=2;i<8;i++)MMnmiss_IMnpip_woK0_woSm_n_mc->Add(MMnmiss_IMnpip_woK0_woSm_n[i]);

  TCanvas *cMMnmiss_IMnpip_woK0_woSm = new TCanvas("cMMnmiss_IMnpip_woK0_woSm","cMMnmiss_IMnpip_woK0_woSm",1200,800);
  cMMnmiss_IMnpip_woK0_woSm->Divide(2,1);
  cMMnmiss_IMnpip_woK0_woSm->cd(1);
  MMnmiss_IMnpip_woK0_woSm_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSm}{  Real data}");
  MMnmiss_IMnpip_woK0_woSm_rdata->Draw("colz");
  cMMnmiss_IMnpip_woK0_woSm->cd(2);
  MMnmiss_IMnpip_woK0_woSm_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSm}{  MC sum}");
  MMnmiss_IMnpip_woK0_woSm_mc->Draw("colz");
  
  TCanvas *cIMnpip_woK0_woSm_n = new TCanvas("cIMnpip_woK0_woSm_n","cIMnpip_woK0_woSm_n");
  cIMnpip_woK0_woSm_n->cd();
  TH1D *IMnpip_woK0_woSm_n[8];
  for(int i=0;i<8;i++){
    IMnpip_woK0_woSm_n[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n[i]->ProjectionX(Form("IMnpip_woK0_woSm_n_%s",name[i]));
    IMnpip_woK0_woSm_n[i]->SetLineColor(colordef[i]);
  }
  TH1D* IMnpip_woK0_woSm_n_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSm_n_mc->ProjectionX("IMnpip_woK0_woSm_n_mc");
  IMnpip_woK0_woSm_n[0]->Draw("HE");
  for(int i=1;i<8;i++)IMnpip_woK0_woSm_n[i]->Draw("HEsame");
  IMnpip_woK0_woSm_n_mc->SetLineColor(6);
  IMnpip_woK0_woSm_n_mc->Draw("HEsame");

  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpip_woK0_woSidn_mc = (TH2D*)MMnmiss_IMnpip_woK0_woSidn[1]->Clone("MMnmiss_IMnpip_woK0_woSidn_mc");
  for(int i=2;i<8;i++) MMnmiss_IMnpip_woK0_woSidn_mc->Add(MMnmiss_IMnpip_woK0_woSidn[i]);

  TCanvas *cMMnmiss_IMnpip_woSidn = new TCanvas("cMMnmiss_IMnpip_woSidn","cMMnmiss_IMnpip_woSidn",1200,800);
  cMMnmiss_IMnpip_woSidn->Divide(2,1);
  cMMnmiss_IMnpip_woSidn->cd(1);
  MMnmiss_IMnpip_woK0_woSidn_rdata->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSidn}{  Real data}");
  MMnmiss_IMnpip_woK0_woSidn_rdata->Draw("colz");
  cMMnmiss_IMnpip_woSidn->cd(2);
  MMnmiss_IMnpip_woK0_woSidn_mc->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSidn}{  MC sum}");
  MMnmiss_IMnpip_woK0_woSidn_mc->Draw("colz");
  
  //projection to Missing mass (miss n & Sigma+/-)
  TCanvas *cMMnmiss_woK0_woSidn = new TCanvas("cMMnmiss_woK0_woSidn","cMMnmiss_woK0_woSidn");
  cMMnmiss_woK0_woSidn->cd();
  TH1D* MMnmiss_woK0_woSidn[8];
  for(int i=0;i<8;i++) MMnmiss_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSidn[i]->ProjectionY(Form("MMnmiss_woK0_woSidn_%s",name[i]));
  MMnmiss_woK0_woSidn[0]->Draw("HE");
  TH1D* MMnmiss_woK0_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mc->ProjectionY("MMnmiss_woK0_woSidn_mc");
  MMnmiss_woK0_woSidn_mc->SetLineColor(6);
  MMnmiss_woK0_woSidn_mc->Draw("HEsame");
  
  for(int i=1;i<8;i++){
    MMnmiss_woK0_woSidn[i]->SetLineColor(colordef[i]);
    MMnmiss_woK0_woSidn[i]->Draw("HEsame");
  }
  
  //Data/MC before modifying MC data
  TCanvas *cMMnmiss_woK0_woSidn_ratio = new TCanvas("cMMnmiss_woK0_woSidn_ratio","cMMnmiss_woK0_woSid_ratio");
  TH1D* MMnmiss_woK0_woSidn_ratio = MMnmiss_woK0_woSidn[0]->Clone("MMnmiss_woK0_woSidn_ratio");
  MMnmiss_woK0_woSidn_ratio->Divide(MMnmiss_woK0_woSidn_mc);
  MMnmiss_woK0_woSidn_ratio->SetTitle("Data/MC");
  MMnmiss_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  MMnmiss_woK0_woSidn_ratio->Draw("HE");
  


  TF1 *sgf_MMnmiss = new TF1("sgf_MMnmiss","[0]*exp(-0.5*pow((x-[1])/([2]+(x<[1])*[3]*(x-[1])),2))");    
  TF1 *fgaus_MMnmiss_high = new TF1("fgaus_MMnmiss_high","gaus");
  //TF1 *fgaus_MMnmiss = new TF1("fgaus_MMnmiss","gaus",0,1.0);
  //fgaus_MMnmiss->SetParameter(0,1.82171e+00); 
  //fgaus_MMnmiss->SetParameter(1,8.56016e-01); 
  //fgaus_MMnmiss->SetParameter(2,6.81677e-01); 
  //fgaus_MMnmiss->SetLineColor(2);
  sgf_MMnmiss->SetParameter(0,1.82171e+00); 
  sgf_MMnmiss->SetParameter(1,8.56016e-01); 
  sgf_MMnmiss->SetParameter(2,6.81677e-01); 
  sgf_MMnmiss->SetLineColor(2);

  //MMnmiss_woK0_woSidn_ratio->Fit("fgaus_MMnmiss","","",0.,1.0);
  MMnmiss_woK0_woSidn_ratio->Fit("sgf_MMnmiss","","",0.,1.116);
  MMnmiss_woK0_woSidn_ratio->Fit("fgaus_MMnmiss_high","","",1.116,1.5);
  double param_sgf_MMnmiss[4];
  for(int ipar=0;ipar<4;ipar++) param_sgf_MMnmiss[ipar] = sgf_MMnmiss->GetParameter(ipar);
  double param_fgaus_MMnmiss_high[3];
  for(int ipar=0;ipar<3;ipar++) param_fgaus_MMnmiss_high[ipar] = fgaus_MMnmiss_high->GetParameter(ipar);
  sgf_MMnmiss->SetLineColor(3);
  sgf_MMnmiss->Draw("same");

  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_woK0_woSidn = new TCanvas("cIMnpip_woK0_woSidn","cIMnpip_woK0_woSidn");
  cIMnpip_woK0_woSidn->cd();
  TH1D* IMnpip_woK0_woSidn[8];
  for(int i=0;i<8;i++)IMnpip_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSidn[i]->ProjectionX(Form("IMnpip_woK0_woSidn_%s",name[i]));
  IMnpip_woK0_woSidn[0]->Draw("HE");//rdata
  TH1D* IMnpip_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mc->ProjectionX("IMnpip_woK0_woSidn_mc");
  IMnpip_woSidn_mc->SetLineColor(6);
  IMnpip_woSidn_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMnpip_woK0_woSidn[i]->SetLineColor(colordef[i]);
    IMnpip_woK0_woSidn[i]->Draw("HEsame");
  }
  
  TCanvas *cIMnpip_woK0_woSidn_ratio = new TCanvas("cIMnpip_woK0_woSidn_ratio","cIMnpip_woK0_woSidn_ratio");
  cIMnpip_woK0_woSidn_ratio->cd();
  TH1D* IMnpip_woK0_woSidn_ratio = IMnpip_woK0_woSidn[0]->Clone("IMnpip_woK0_woSidn_ratio");
  IMnpip_woK0_woSidn_ratio->Divide(IMnpip_woK0_woSidn_mc);
  IMnpip_woK0_woSidn_ratio->SetTitle("IMnpip_woK0_woSidn Data/MC");
  IMnpip_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpip_woK0_woSidn_ratio->Draw("HE");
  /*
  TF1* IMnpip_mod_gaus = new TF1("IMnpip_mod_gaus","gaus",1.1,1.7);
  IMnpip_mod_gaus->SetParameters(0,0.8);
  IMnpip_mod_gaus->SetParameters(1,1.2);
  IMnpip_mod_gaus->SetParameters(2,0.2);
  IMnpip_mod_gaus->SetLineColor(3);
  IMnpip_woK0_woSidn_ratio->Fit("IMnpip_mod_gaus","","",1.1,1.7);
  */

  //
  //Sigma-
  TH2D *MMnmiss_IMnpim_woK0_mc = (TH2D*)MMnmiss_IMnpim_woK0[1]->Clone("MMnmiss_IMnpim_woK0_mc");
  //adding all MC data
  for(int i=2;i<8;i++)MMnmiss_IMnpim_woK0_mc->Add(MMnmiss_IMnpim_woK0[i]);

  TCanvas *cMMnmiss_IMnpim_woK0 = new TCanvas("cMMnmiss_IMnpim_woK0","cMMnmiss_IMnpim_woK0",1200,800);
  cMMnmiss_IMnpim_woK0->Divide(2,1);
  cMMnmiss_IMnpim_woK0->cd(1);
  MMnmiss_IMnpim_woK0_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0}{  Real data}");
  MMnmiss_IMnpim_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0->cd(2);
  MMnmiss_IMnpim_woK0_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0}{ MC sum}");
  MMnmiss_IMnpim_woK0_mc->Draw("colz");
  
  //w/o K0, w/o Sp mode
  TH2D *MMnmiss_IMnpim_woK0_woSp_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSp[1]->Clone("MMnmiss_IMnpim_woK0_woSp_mc");
  //adding all MC data
  for(int i=2;i<8;i++){
    //MMnmiss_IMnpim_woK0_woSp[i]->Print("base");
    MMnmiss_IMnpim_woK0_woSp_mc->Add(MMnmiss_IMnpim_woK0_woSp[i]);
  }
  //w/o K0, w/o Sp mode, selecting missing neutron
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSp_n[1]->Clone("MMnmiss_IMnpim_woK0_woSp_n_mc");
  for(int i=2;i<8;i++)MMnmiss_IMnpim_woK0_woSp_n_mc->Add(MMnmiss_IMnpim_woK0_woSp_n[i]);

  TCanvas *cMMnmiss_IMnpim_woK0_woSp = new TCanvas("cMMnmiss_IMnpim_woK0_woSp","cMMnmiss_IMnpim_woK0_woSp",1200,800);
  cMMnmiss_IMnpim_woK0_woSp->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSp->cd(1);
  MMnmiss_IMnpim_woK0_woSp_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSp}{  Real data}");
  MMnmiss_IMnpim_woK0_woSp_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSp->cd(2);
  MMnmiss_IMnpim_woK0_woSp_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSp}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSp_mc->Draw("colz");

  TCanvas *cIMnpim_woK0_woSp_n = new TCanvas("cIMnpim_woK0_woSp_n","cIMnpim_woK0_woSp_n");
  cIMnpim_woK0_woSp_n->cd();
  TH1D *IMnpim_woK0_woSp_n[8];
  for(int i=0;i<8;i++){
    IMnpim_woK0_woSp_n[i] = (TH1D*)MMnmiss_IMnpim_woK0_woSp_n[i]->ProjectionX(Form("IMnpim_woK0_woSp_n_%s",name[i]));
    IMnpim_woK0_woSp_n[i]->SetLineColor(colordef[i]);
  }
  TH1D* IMnpim_woK0_woSp_n_mc = (TH1D*) MMnmiss_IMnpim_woK0_woSp_n_mc->ProjectionX("IMnpim_woK0_woSp_n_mc");
  IMnpim_woK0_woSp_n[0]->Draw("HE");
  for(int i=1;i<8;i++)IMnpim_woK0_woSp_n[i]->Draw("HEsame");
  IMnpim_woK0_woSp_n_mc->SetLineColor(6);
  IMnpim_woK0_woSp_n_mc->Draw("HEsame");


  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpim_woK0_woSidn_mc = (TH2D*)MMnmiss_IMnpim_woK0_woSidn[1]->Clone("MMnmiss_IMnpim_woK0_woSidn_mc");
  for(int i=2;i<8;i++)MMnmiss_IMnpim_woK0_woSidn_mc->Add(MMnmiss_IMnpim_woK0_woSidn[i]);

  TCanvas *cMMnmiss_IMnpim_woK0_woSidn = new TCanvas("cMMnmiss_IMnpim_woK0_woSidn","cMMnmiss_IMnpim_woK0_woSidn",1200,800);
  cMMnmiss_IMnpim_woK0_woSidn->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSidn->cd(1);
  MMnmiss_IMnpim_woK0_woSidn_rdata->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSidn}{  Real data}");
  MMnmiss_IMnpim_woK0_woSidn_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSidn->cd(2);
  MMnmiss_IMnpim_woK0_woSidn_mc->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSidn}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSidn_mc->Draw("colz");
  
  TCanvas *cIMnpim_woK0_woSidn = new TCanvas("IMnpim_woK0_woSidn","IMnpim_woK0_woSidn");
  cIMnpim_woK0_woSidn->cd();
  TH1D* IMnpim_woK0_woSidn[8];
  for(int i=0;i<8;i++)IMnpim_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpim_woK0_woSidn[i]->ProjectionX(Form("IMnpim_woK0_woSidn_%s",name[i]));
  TH1D* IMnpim_woK0_woSidn_mc = MMnmiss_IMnpim_woK0_woSidn_mc->ProjectionX("IMnpim_woK0_woSidn_mc");
  IMnpim_woK0_woSidn[0]->Draw("HE");
  IMnpim_woK0_woSidn_mc->SetLineColor(6);
  IMnpim_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMnpim_woK0_woSidn[i]->SetLineColor(colordef[i]);
    IMnpim_woK0_woSidn[i]->Draw("HEsame");
  }
  
  TCanvas *cIMnpim_woK0_woSidn_ratio = new TCanvas("cIMnpim_woK0_woSidn_ratio","cIMnpim_woK0_woSidn_ratio");
  cIMnpim_woK0_woSidn_ratio->cd();
  TH1D* IMnpim_woK0_woSidn_ratio = (TH1D*)IMnpim_woK0_woSidn[0]->Clone("IMnpim_woK0_woSidn_ratio");
  IMnpim_woK0_woSidn_ratio->Divide(IMnpim_woK0_woSidn_mc);
  IMnpim_woK0_woSidn_ratio->SetTitle("IMnpim_woK0_woSidn Data/MC");
  IMnpim_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpim_woK0_woSidn_ratio->Draw("HE");
   
  //MMnmiss vs IMpippim wSid
  TH2D* MMnmiss_IMpippim_wSid_mc = MMnmiss_IMpippim_wSid[1]->Clone("MMnmiss_IMpippim_wSid_mc");
  for(int i=2;i<8;i++) MMnmiss_IMpippim_wSid_mc->Add(MMnmiss_IMpippim_wSid[i]);
  
  TCanvas *cMMnmiss_IMpippim_wSid = new TCanvas("cMMnmiss_IMpippim_wSid","cMMnmiss_IMpippim_wSid",1200,800);
  cMMnmiss_IMpippim_wSid->Divide(2,1);
  cMMnmiss_IMpippim_wSid->cd(1);
  MMnmiss_IMpippim_wSid[0]->SetTitle("#splitline{MMnmiss_IMpippim_wSid}{  Real data}");
  MMnmiss_IMpippim_wSid[0]->Draw("colz");
  cMMnmiss_IMpippim_wSid->cd(2);
  MMnmiss_IMpippim_wSid_mc->SetTitle("#splitline{MMnmiss_IMpippim_wSid}{  MC sum}");
  MMnmiss_IMpippim_wSid_mc->Draw("colz");

  TH1D* IMpippim_wSid_mc = (TH1D*)MMnmiss_IMpippim_wSid_mc->ProjectionX("IMpippim_wSid_mc");
  IMpippim_wSid_mc->SetLineColor(6);
  TCanvas *cIMpippim_wSid = new TCanvas("cIMpippim_wSid","cIMpippim_wSid");
  cIMpippim_wSid->cd();
  TH1D* IMpippim_wSid[8];
  for(int i=0;i<8;i++) IMpippim_wSid[i] = (TH1D*)MMnmiss_IMpippim_wSid[i]->ProjectionX(Form("IMpippim_wSid_%s",name[i]));
  IMpippim_wSid[0]->Draw("HE");
  IMpippim_wSid_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMpippim_wSid[i]->SetLineColor(colordef[i]);
    IMpippim_wSid[i]->Draw("HEsame");
  }

  
  //MMnmiss vs IMpippim w/o K0 w/o (Sid & n);
  TH2D* MMnmiss_IMpippim_woK0_woSidn_mc = MMnmiss_IMpippim_woK0_woSidn[1]->Clone("MMnmiss_IMpippim_woK0_woSidn_mc");
  for(int i=2;i<8;i++) MMnmiss_IMpippim_woK0_woSidn_mc->Add(MMnmiss_IMpippim_woK0_woSidn[i]);

  TCanvas *cMMnmiss_IMpippim_woK0_woSidn = new TCanvas("cMMnmiss_IMpippim_woK0_woSidn","cMMnmiss_IMpippim_woK0_woSidn",1200,800);
  cMMnmiss_IMpippim_woK0_woSidn->Divide(2,1);
  cMMnmiss_IMpippim_woK0_woSidn->cd(1);
  MMnmiss_IMpippim_woK0_woSidn[0]->SetTitle("#splitline{MMnmiss_IMpippim_woK0_woSidn}{  Real data}");
  MMnmiss_IMpippim_woK0_woSidn[0]->Draw("colz");
  cMMnmiss_IMpippim_woK0_woSidn->cd(2);
  MMnmiss_IMpippim_woK0_woSidn_mc->SetTitle("#splitline{MMnmiss_IMpippim_woK0_woSidn}{  MC sum}");
  MMnmiss_IMpippim_woK0_woSidn_mc->Draw("colz");
  
  TH1D* IMpippim_woK0_woSidn_mc = (TH1D*)MMnmiss_IMpippim_woK0_woSidn_mc->ProjectionX("IMpippim_woK0_woSidn_mc");
  IMpippim_woK0_woSidn_mc->SetLineColor(6);
  TCanvas *cIMpippim_woK0_woSidn = new TCanvas("cIMpippim_woK0_woSidn","cIMpippim_woK0_woSidn");
  cIMpippim_woK0_woSidn->cd();
  TH1D* IMpippim_woK0_woSidn[8];
  for(int i=0;i<8;i++) IMpippim_woK0_woSidn[i] = (TH1D*)MMnmiss_IMpippim_woK0_woSidn[i]->ProjectionX(Form("IMpippim_woK0_woSidn_%s",name[i]));
  IMpippim_woK0_woSidn[0]->Draw("HE");
  IMpippim_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMpippim_woK0_woSidn[i]->SetLineColor(colordef[i]);
    IMpippim_woK0_woSidn[i]->Draw("HEsame");
  }
  
  //q vs IMnpipi w/o K0 w/o (Sid & n);
  TH2D* q_IMnpipi_woK0_woSidn_mc = q_IMnpipi_woK0_woSidn[1]->Clone("q_IMnpipi_woK0_woSidn_mc");
  for(int i=2;i<8;i++)q_IMnpipi_woK0_woSidn_mc->Add(q_IMnpipi_woK0_woSidn[i]);

  TCanvas *cq_IMnpipi_woK0_woSidn = new TCanvas("cq_IMnpipi_woK0_woSidn","q_IMnpipi_woK0_woSidn",1200,800);
  cq_IMnpipi_woK0_woSidn->Divide(2,1);
  cq_IMnpipi_woK0_woSidn->cd(1);
  q_IMnpipi_woK0_woSidn[0]->SetTitle("#splitline{q_IMnpipi_woK0_woSidn}{  Real data}");
  q_IMnpipi_woK0_woSidn[0]->Draw("colz");
  cq_IMnpipi_woK0_woSidn->cd(2);
  q_IMnpipi_woK0_woSidn_mc->SetTitle("#splitline{q_IMnpipi_woK0_woSidn}{  MC sum}");
  q_IMnpipi_woK0_woSidn_mc->Draw("colz");

  TH1D* IMnpipi_woK0_woSidn_mc = (TH1D*)q_IMnpipi_woK0_woSidn_mc->ProjectionX("IMnpipi_woK0_woSidn_mc");
  IMnpipi_woK0_woSidn_mc->SetLineColor(6);
  TCanvas *cIMnpipi_woK0_woSidn = new TCanvas("cIMnpipi_woK0_woSidn","cIMnpipi_woK0_woSidn");
  cIMnpipi_woK0_woSidn->cd();
  TH1D* IMnpipi_woK0_woSidn[8];
  for(int i=0;i<8;i++) IMnpipi_woK0_woSidn[i] = (TH1D*)q_IMnpipi_woK0_woSidn[i]->ProjectionX(Form("IMnpipi_woK0_woSidn_%s",name[i]));
  IMnpipi_woK0_woSidn[0]->Draw("HE");
  IMnpipi_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMnpipi_woK0_woSidn[i]->SetLineColor(colordef[i]);
    IMnpipi_woK0_woSidn[i]->Draw("HEsame");
  }
   
  TH1D* q_woK0_woSidn_mc = (TH1D*)q_IMnpipi_woK0_woSidn_mc->ProjectionY("q_woK0_woSidn_mc");
  q_woK0_woSidn_mc->SetLineColor(6);
  TCanvas *cq_woK0_woSidn = new TCanvas("cq_woK0_woSidn","cq_woK0_woSidn");
  cq_woK0_woSidn->cd();
  TH1D* q_woK0_woSidn[8];
  for(int i=0;i<8;i++) q_woK0_woSidn[i] = (TH1D*)q_IMnpipi_woK0_woSidn[i]->ProjectionY(Form("q_woK0_woSidn_%s",name[i]));
  q_woK0_woSidn[0]->Draw("HE");
  q_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    q_woK0_woSidn[i]->SetLineColor(colordef[i]);
    q_woK0_woSidn[i]->Draw("HEsame");
  }

  //missing mass vs Mom(pi+pi-) w/o K0 w/o (Sid & n)
  TH2D* MMnmiss_Mompippim_woK0_woSidn_mc = (TH2D*)MMnmiss_Mompippim_woK0_woSidn[1]->Clone("MMnmiss_Mompippim_woK0_woSidn_mc");
  for(int i=2;i<8;i++)MMnmiss_Mompippim_woK0_woSidn_mc->Add(MMnmiss_Mompippim_woK0_woSidn[i]);

  TCanvas *cMMnmiss_Mompippim_woK0_woSidn = new TCanvas("cMMnmiss_Mompippim_woK0_woSidn","cMMnmiss_Mompippim_woK0_woSidn",1200,800);
  cMMnmiss_Mompippim_woK0_woSidn->Divide(2,1);
  cMMnmiss_Mompippim_woK0_woSidn->cd(1);
  MMnmiss_Mompippim_woK0_woSidn[0]->SetTitle("#splitline{MMnmiss_Mompippim_woK0_woSidn}{  Real data}");
  MMnmiss_Mompippim_woK0_woSidn[0]->Draw("colz");
  cMMnmiss_Mompippim_woK0_woSidn->cd(2);
  MMnmiss_Mompippim_woK0_woSidn_mc->SetTitle("#splitline{MMnmiss_Mompippim_woK0_woSidn}{  MC sum}");
  MMnmiss_Mompippim_woK0_woSidn_mc->Draw("colz");

  TH1D* Mompippim_woK0_woSidn_mc = (TH1D*)MMnmiss_Mompippim_woK0_woSidn_mc->ProjectionX("Mompippim_woK0_woSidn_mc");
  Mompippim_woK0_woSidn_mc->SetLineColor(6);
  TCanvas *cMompippim_woK0_woSidn = new TCanvas("cMompippim_woK0_woSidn","cMompippim_woK0_woSidn");
  cMompippim_woK0_woSidn->cd();
  TH1D* Mompippim_woK0_woSidn[8];
  for(int i=0;i<8;i++) Mompippim_woK0_woSidn[i] = (TH1D*)MMnmiss_Mompippim_woK0_woSidn[i]->ProjectionX(Form("MMnmiss_Mompippim_woK0_woSidn_%s",name[i]));
  Mompippim_woK0_woSidn[0]->Draw("HE");
  Mompippim_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1;i<8;i++){
    Mompippim_woK0_woSidn[i]->SetLineColor(colordef[i]);
    Mompippim_woK0_woSidn[i]->Draw("HEsame");
  }


  //weighting and decomposition of BG
  int nbinx_MM_npip = MMnmiss_IMnpip_woK0_woSidn[0]->GetNbinsX();
  int nbiny_MM_npip = MMnmiss_IMnpip_woK0_woSidn[0]->GetNbinsY();
  TH2D* MMnmiss_IMnpip_woK0_woSidn_mod[8];
  TH2D* MMnmiss_IMnpim_woK0_woSidn_mod[8];
  for(int i=0;i<8;i++){
    MMnmiss_IMnpip_woK0_woSidn_mod[i] = (TH2D*)MMnmiss_IMnpip_woK0_woSidn[i]->Clone(Form("MMnmiss_IMnpip_woK0_woSidn_%s_mod",name[i]));
    MMnmiss_IMnpim_woK0_woSidn_mod[i] = (TH2D*)MMnmiss_IMnpim_woK0_woSidn[i]->Clone(Form("MMnmiss_IMnpim_woK0_woSidn_%s_mod",name[i]));
  }  
  
  for(int i=3;i<4;i++){
   // MMnmiss_IMnpip_woK0_woSidn_mod[i]->Reset("ICESM");
   // MMnmiss_IMnpim_woK0_woSidn_mod[i]->Reset("ICESM");
  }
  
  /*
  TF1 *fgaus_MMnmiss_mod = new TF1("fgaus_MMnmiss_mod","gaus",0,1.0);
  fgaus_MMnmiss_mod->SetParameter(0,1.82171e+00); 
  fgaus_MMnmiss_mod->SetParameter(1,8.56016e-01); 
  fgaus_MMnmiss_mod->SetParameter(2,6.81677e-01); 
  */

  TF1 *sgf_MMnmiss_mod = new TF1("sgf_MMnmiss_mod","[0]*exp(-0.5*pow((x-[1])/([2]+(x<[1])*[3]*(x-[1])),2))",0,1.115);
  for(int ipar=0;ipar<4;ipar++) sgf_MMnmiss_mod->SetParameter(ipar, param_sgf_MMnmiss[ipar]);
  sgf_MMnmiss_mod->SetLineColor(3);

  TF1 *fgaus_MMnmiss_high_mod = new TF1("fgaus_MMnmiss_high_mod","gaus",0,1.5);
  for(int ipar=0;ipar<3;ipar++) fgaus_MMnmiss_high_mod->SetParameter(ipar,param_fgaus_MMnmiss_high[ipar]);
  
  TCanvas *cfgaus_MMnmiss_mod = new TCanvas("cfgaus_MMnmiss_mod","cfgaus_MMnmiss_mod");
  cfgaus_MMnmiss_mod->cd();
  fgaus_MMnmiss_high_mod->Draw("");
  sgf_MMnmiss_mod->Draw("same");


  
  TH2D* MMnmiss_IMnpip_woK0_woSidn_mc_mod = (TH2D*)MMnmiss_IMnpip_woK0_woSidn_mod[1]->Clone("MMnmiss_IMnpip_woK0_woSidn_mc_mod");
  TH2D* MMnmiss_IMnpim_woK0_woSidn_mc_mod = (TH2D*)MMnmiss_IMnpim_woK0_woSidn_mod[1]->Clone("MMnmiss_IMnpim_woK0_woSidn_mc_mod");
  MMnmiss_IMnpip_woK0_woSidn_mc_mod->Reset("ICESM");
  MMnmiss_IMnpim_woK0_woSidn_mc_mod->Reset("ICESM");


  for(int ix=0;ix<nbinx_MM_npip;ix++){
    for(int iy=0;iy<nbiny_MM_npip;iy++){
      double contIMnpip[8]={0.0};
      double contIMnpim[8]={0.0};
      double contIMnpip_mc=0.0;//for mc sum
      double contIMnpim_mc=0.0;//for mc sum
      
      double xval = MMnmiss_IMnpip_woK0_woSidn_mc->GetXaxis()->GetBinCenter(ix);
      double yval = MMnmiss_IMnpip_woK0_woSidn_mc->GetYaxis()->GetBinCenter(iy);

      for(int imc=0;imc<8;imc++){
        contIMnpip[imc] = MMnmiss_IMnpip_woK0_woSidn[imc]->GetBinContent(ix,iy);
        contIMnpim[imc] = MMnmiss_IMnpim_woK0_woSidn[imc]->GetBinContent(ix,iy);
        contIMnpip_mc = MMnmiss_IMnpip_woK0_woSidn_mc->GetBinContent(ix,iy);
        contIMnpim_mc = MMnmiss_IMnpim_woK0_woSidn_mc->GetBinContent(ix,iy);
      }
      double weight = 1.0;
      if(yval<=1.115){
        weight = sgf_MMnmiss_mod->Eval(yval);
      }else{
        weight = fgaus_MMnmiss_high_mod->Eval(yval);
      }
        MMnmiss_IMnpip_woK0_woSidn_mc_mod->SetBinContent(ix,iy,contIMnpip_mc*weight);
        MMnmiss_IMnpim_woK0_woSidn_mc_mod->SetBinContent(ix,iy,contIMnpim_mc*weight);
    }
  }
  
  //TH2D* MMnmiss_IMnpip_woK0_woSidn_mc_mod = (TH2D*)MMnmiss_IMnpip_woK0_woSidn_mod[1]->Clone("MMnmiss_IMnpip_woK0_woSidn_mc_mod");
  //for(int i=2;i<8;i++) MMnmiss_IMnpip_woK0_woSidn_mc_mod->Add(MMnmiss_IMnpip_woK0_woSidn_mod[i]);
  
  //TH2D* MMnmiss_IMnpim_woK0_woSidn_mc_mod = (TH2D*)MMnmiss_IMnpim_woK0_woSidn_mod[1]->Clone("MMnmiss_IMnpim_woK0_woSidn_mc_mod");
  //for(int i=2;i<8;i++) MMnmiss_IMnpim_woK0_woSidn_mc_mod->Add(MMnmiss_IMnpim_woK0_woSidn_mod[i]);

  TCanvas *cMMnmiss_IMnpip_woK0_woSidn_mod = new TCanvas("cMMnmiss_IMnpip_woK0_woSidn_mod","MMnmiss_IMnpip_woK0_woSidn_mod",1200,800);
  cMMnmiss_IMnpip_woK0_woSidn_mod->Divide(2,1);
  cMMnmiss_IMnpip_woK0_woSidn_mod->cd(1);
  MMnmiss_IMnpip_woK0_woSidn_mod[0]->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSidn_mod}{  Real data}");
  MMnmiss_IMnpip_woK0_woSidn_mod[0]->Draw("colz");
  cMMnmiss_IMnpip_woK0_woSidn_mod->cd(2);
  MMnmiss_IMnpip_woK0_woSidn_mc_mod->SetTitle("#splitline{MMnmiss_IMnpip_woK0_woSidn_mod}{  MC sum}");
  MMnmiss_IMnpip_woK0_woSidn_mc_mod->Draw("colz");

  TCanvas *cMMnmiss_IMnpim_woK0_woSidn_mod = new TCanvas("cMMnmiss_IMnpim_woK0_woSidn_mod","MMnmiss_IMnpim_woK0_woSidn_mod",1200,800);
  cMMnmiss_IMnpim_woK0_woSidn_mod->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSidn_mod->cd(1);
  MMnmiss_IMnpim_woK0_woSidn_mod[0]->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSidn_mod}{  Real data}");
  MMnmiss_IMnpim_woK0_woSidn_mod[0]->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSidn_mod->cd(2);
  MMnmiss_IMnpim_woK0_woSidn_mc_mod->SetTitle("#splitline{MMnmiss_IMnpim_woK0_woSidn_mod}{  MC sum}");
  MMnmiss_IMnpim_woK0_woSidn_mc_mod->Draw("colz");



  TCanvas *cMMnmiss_woK0_woSidn_mod = new TCanvas("cMMnmiss_woK0_woSidn_mod","MMnmiss_woK0_woSidn_mod");
  cMMnmiss_woK0_woSidn_mod->cd(); 
  MMnmiss_woK0_woSidn[0]->Draw("HE");
  double rdatamax = MMnmiss_woK0_woSidn[0]->GetMaximum();
  TH1D* MMnmiss_woK0_woSidn_mod[8];
  for(int i=0;i<8;i++){
    MMnmiss_woK0_woSidn_mod[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mod[i]->ProjectionY(Form("MMnmiss_woK0_woSidn_mod_%s",name[i]));
    MMnmiss_woK0_woSidn_mod[i]->SetLineColor(colordef[i]);
  }
  for(int i=1;i<8;i++)MMnmiss_woK0_woSidn_mod[i]->Draw("HEsame");
  TH1D*  MMnmiss_woK0_woSidn_mc_mod =  (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mc_mod->ProjectionY("MMnmiss_woK0_woSidn_mc_mod");
  MMnmiss_woK0_woSidn_mc_mod->SetLineColor(6);
  MMnmiss_woK0_woSidn_mc_mod->Draw("HEsame");
  
  TCanvas *cMMnmiss_woK0_woSidn_mod_ratio = new TCanvas("MMnmiss_woK0_woSidn_mod_ratio","MMnmiss_woK0_woSidn_mod_ratio");
  cMMnmiss_woK0_woSidn_mod_ratio->cd(); 
  TH1D* MMnmiss_woK0_woSidn_mc_mod_ratio = (TH1D*)MMnmiss_woK0_woSidn[0]->Clone("MMnmiss_woK0_woSidn_mod_ratio");
  MMnmiss_woK0_woSidn_mc_mod_ratio->Divide(MMnmiss_woK0_woSidn_mc_mod);
  MMnmiss_woK0_woSidn_mc_mod_ratio->SetTitle("MMnmiss_woK0_woSidn Data/MC");
  MMnmiss_woK0_woSidn_mc_mod_ratio->GetYaxis()->SetRangeUser(-0.1,3);
  MMnmiss_woK0_woSidn_mc_mod_ratio->Draw("HE");
  

  TCanvas *cIMnpip_woK0_woSidn_mod = new TCanvas("cIMnpip_woK0_woSidn_mod","cIMnpip_woK0_woSidn_mod");
  cIMnpip_woK0_woSidn_mod->cd();
  TH1D* IMnpip_woK0_woSidn_mod[8];
  for(int i=0;i<8;i++){
    IMnpip_woK0_woSidn_mod[i] = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mod[i]->ProjectionX(Form("IMnpip_woK0_woSidn_mod_%s",name[i]));
    IMnpip_woK0_woSidn_mod[i]->SetLineColor(colordef[i]);
  }
  IMnpip_woK0_woSidn[0]->Draw("HE");
  TH1D* IMnpip_woK0_woSidn_mc_mod = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mc_mod->ProjectionX("IMnpip_woK0_woSidn_mc_mod");
  IMnpip_woK0_woSidn_mc_mod->SetLineColor(6);
  IMnpip_woK0_woSidn_mc_mod->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMnpip_woK0_woSidn_mod[i]->Draw("HEsame");
  }
  
  TCanvas *cIMnpip_woK0_woSidn_mod_ratio = new TCanvas("cIMnpip_woK0_woSidn_mod_ratio","cIMnpip_woK0_woSidn_mod_ratio");
  cIMnpip_woK0_woSidn_mod_ratio->cd();
  //TH1D* IMnpip_woK0_woSidn_mod_ratio = (TH1D*)IMnpip_woK0_woSidn[0]->Clone("IMnpip_woK0_woSidn_mod_ratio");
  //IMnpip_woK0_woSidn_mod_ratio->Divide(IMnpip_woK0_woSidn_mc_mod);
  //IMnpip_woK0_woSidn_mod_ratio->SetTitle("IMnpip_woK0_woSidn_mod_ratio Data/MC");
  TH1D* IMnpip_woK0_woSidn_mod_ratio = (TH1D*)IMnpip_woK0_woSidn_mc->Clone("IMnpip_woK0_woSidn_mod_ratio");
  IMnpip_woK0_woSidn_mod_ratio->Divide(IMnpip_woK0_woSidn[0]);
  IMnpip_woK0_woSidn_mod_ratio->SetTitle("IMnpip_woK0_woSidn_mod_ratio MC/Data");
  IMnpip_woK0_woSidn_mod_ratio->GetYaxis()->SetRangeUser(-0.1,3);
  IMnpip_woK0_woSidn_mod_ratio->Draw("HE");
  
  TF1 *sgf_IMnpip_mod1 = new TF1("sgf_IMnpip_mod1","[0]*exp(-0.5*pow((x-[1])/([2]+(x<[1])*[3]*(x-[1])),2))");    
  sgf_IMnpip_mod1->SetParameter(1,1.2);
  IMnpip_woK0_woSidn_mod_ratio->Fit("sgf_IMnpip_mod1","","",1.1,1.6);


  TCanvas *cIMnpim_woK0_woSidn_mod = new TCanvas("cIMnpim_woK0_woSidn_mod","cIMnpim_woK0_woSidn_mod");
  cIMnpim_woK0_woSidn_mod->cd();
  TH1D* IMnpim_woK0_woSidn_mod[8];
  for(int i=0;i<8;i++){
    IMnpim_woK0_woSidn_mod[i] = (TH1D*)MMnmiss_IMnpim_woK0_woSidn_mod[i]->ProjectionX(Form("IMnpim_woK0_woSidn_mod_%s",name[i]));
    IMnpim_woK0_woSidn_mod[i]->SetLineColor(colordef[i]);
  }
  IMnpim_woK0_woSidn[0]->Draw("HE");
  TH1D* IMnpim_woK0_woSidn_mc_mod = (TH1D*)MMnmiss_IMnpim_woK0_woSidn_mc_mod->ProjectionX("IMnpim_woK0_woSidn_mc_mod");
  IMnpim_woK0_woSidn_mc_mod->SetLineColor(6);
  IMnpim_woK0_woSidn_mc_mod->Draw("HEsame");
  for(int i=1;i<8;i++){
    IMnpim_woK0_woSidn_mod[i]->Draw("HEsame");
  }
  
  TCanvas *cIMnpim_woK0_woSidn_mod_ratio = new TCanvas("cIMnpim_woK0_woSidn_mod_ratio","cIMnpim_woK0_woSidn_mod_ratio");
  cIMnpim_woK0_woSidn_mod_ratio->cd();
  TH1D* IMnpim_woK0_woSidn_mod_ratio = (TH1D*)IMnpim_woK0_woSidn[0]->Clone("IMnpim_woK0_woSidn_mod_ratio");
  IMnpim_woK0_woSidn_mod_ratio->Divide(IMnpim_woK0_woSidn_mc_mod);
  IMnpim_woK0_woSidn_mod_ratio->SetTitle("IMnpim_woK0_woSidn_mod_ratio Data/MC");
  IMnpim_woK0_woSidn_mod_ratio->GetYaxis()->SetRangeUser(-0.1,3);
  IMnpim_woK0_woSidn_mod_ratio->Draw("HE");
  


  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname="plot_miss2D_out.pdf";
  for(int i=0;i<size;i++){
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("Title:%s",c->GetTitle())); 
  }

  TIter nexthist2(gDirectory->GetList());

  TFile *fout = new TFile("plot_miss2D_out.root","RECREATE");
  while( (obj = (TObject*)nexthist2())!=NULL ){
    obj->Write();
  }
  fout->Print();
  fout->cd();
  fout->Close();
}

