#include "../src/GlobalVariables.h"

void plot_K0()
{

  TFile *_file0 = TFile::Open("../post/evanaIMpisigma_npippim_v162_outncutK015.root");
  _file0->cd();
  TH2F* q_IMpippim_n_rdata = (TH2F*)_file0->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_rdata = (TH2F*)_file0->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_rdata = (TH2F*)_file0->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_rdata = (TH2F*)_file0->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_rdata = (TH2F*)_file0->Get("MMnmiss_IMpippim_dE");

  TH1D* IMpippim_rdata = (TH1D*)q_IMpippim_n_rdata->ProjectionX("IMpippim_rdata");
  TH1D* IMnpipi_rdata = (TH1D*)nmom_IMnpipi_wK0_n_rdata->ProjectionX("IMnpipi_rdata");
  TH1D* nmom_IMnpipi_wK0_n_rdata_py = (TH1D*) nmom_IMnpipi_wK0_n_rdata->ProjectionY("nmom_IMnpipi_wK0_n_rdata_py");
  TH1D* q_rdata = (TH1D*)q_IMnpipi_wK0_n_rdata->ProjectionY("q_rdata");
  TH1D* Mompippim_rdata = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_rdata->ProjectionY("Mompippim_rdata"); 
   
  const double pipi_MIN = 0.487;// (almost 1.5 sigma)
  const double pipi_MAX = 0.510;// (almost 1.5 sigma)
  int pipimin = MMnmiss_IMpippim_dE_rdata->GetXaxis()->FindBin(pipi_MIN);
  int pipimax = MMnmiss_IMpippim_dE_rdata->GetXaxis()->FindBin(pipi_MAX);
  TH1D* MMnmiss_rdata = (TH1D*) MMnmiss_IMpippim_dE_rdata->ProjectionY("MMnmiss_rdata",pipimin,pipimax) ;

  //K-p -> K0n -> n-n elastic
  TFile *_file1 = TFile::Open("simIMpisigma_K0_nnts_pippimn_v10_outncutK015.root");
  TFile *_file1_1 = TFile::Open("simIMpisigma_K0_nnts_v10.root");
  TH2F* q_IMpippim_n_nnts = (TH2F*)_file1->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_nnts = (TH2F*)_file1->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_nnts = (TH2F*)_file1->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_nnts = (TH2F*)_file1->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_nnts = (TH2F*)_file1->Get("MMnmiss_IMpippim_dE");

  TH1D* IMpippim_nnts = (TH1D*)q_IMpippim_n_nnts->ProjectionX("IMpippim_nnts");
  TH1D* IMnpipi_nnts = (TH1D*) nmom_IMnpipi_wK0_n_nnts->ProjectionX("IMnpipi_nnts");
  TH1D* nmom_IMnpipi_wK0_n_nnts_py = (TH1D*)nmom_IMnpipi_wK0_n_nnts->ProjectionY("nmom_IMnpipi_wK0_n_nnts_py");
  
  TH1D* q_nnts = (TH1D*)q_IMnpipi_wK0_n_nnts->ProjectionY("q_nnts");
  
  TH1D* Mompippim_nnts = (TH1D*) Mompippim_IMnpipi_dE_wK0_n_nnts->ProjectionY("Mompippim_nnts");
  TH1D* MMnmiss_nnts = (TH1D*) MMnmiss_IMpippim_dE_nnts->ProjectionY("MMnmiss_nnts",pipimin,pipimax);

  
  //K-p -> K0n (ns)
  TFile *_file2 = TFile::Open("simIMpisigma_K0n_ns_pippimn_v12_outncutK015.root");
  TH2F* q_IMpippim_n_ns = (TH2F*)_file2->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_ns = (TH2F*)_file2->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_ns = (TH2F*)_file2->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_ns = (TH2F*)_file2->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_ns = (TH2F*)_file2->Get("MMnmiss_IMpippim_dE");


  TH1D* IMpippim_ns = (TH1D*)q_IMpippim_n_ns->ProjectionX("IMpippim_ns");
  TH1D* IMnpipi_ns = (TH1D*) nmom_IMnpipi_wK0_n_ns->ProjectionX();
  TH1D* nmom_IMnpipi_wK0_n_ns_py = (TH1D*) nmom_IMnpipi_wK0_n_ns->ProjectionY("nmom_IMnpipi_wK0_n_ns_py"); 
  TH1D* q_ns = (TH1D*) q_IMnpipi_wK0_n_ns->ProjectionY("q_ns");
  TH1D* Mompippim_ns = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_ns->ProjectionY("Mompippim_ns");
  TH1D* MMnmiss_ns = (TH1D*) MMnmiss_IMpippim_dE_ns->ProjectionY("MMnmiss_ns",pipimin,pipimax);
  

  //K-d -> K0nn phase space
  TFile *_file3 = TFile::Open("simIMpisigma_K0nn_pippimn_v7_outncutK015.root");
  //TFile *_file3 = TFile::Open("simIMpisigma_K0nn_pippimn_v5_outncutK015.root");
  TH2F* q_IMpippim_n_K0nn = (TH2F*)_file3->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_K0nn = (TH2F*)_file3->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_K0nn = (TH2F*)_file3->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_K0nn = (TH2F*)_file3->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_K0nn = (TH2F*)_file3->Get("MMnmiss_IMpippim_dE");

  TH1D* IMpippim_K0nn = (TH1D*) q_IMpippim_n_K0nn->ProjectionX("IMpippim_K0nn");
  TH1D* IMnpipi_K0nn = (TH1D*) nmom_IMnpipi_wK0_n_K0nn->ProjectionX("IMnpipi_K0nn");
  TH1D* nmom_IMnpipi_wK0_n_K0nn_py = (TH1D*)nmom_IMnpipi_wK0_n_K0nn->ProjectionY("nmom_IMnpipi_wK0_n_K0nn_py");
  TH1D* q_K0nn = (TH1D*) q_IMnpipi_wK0_n_K0nn->ProjectionY("q_K0nn");
  TH1D* Mompippim_K0nn = (TH1D*) Mompippim_IMnpipi_dE_wK0_n_K0nn->ProjectionY("Mompippim_K0nn");
  TH1D* MMnmiss_K0nn = (TH1D*) MMnmiss_IMpippim_dE_K0nn->ProjectionY("MMnmiss_K0nn",pipimin,pipimax);

  //K-p -> K0n (n_s),  K0n->K0 n (two step reaction)
  TFile *_file4 = TFile::Open("simIMpisigma_K0_Knts_pippimn_v3_outncutK015.root");
  TH2F* q_IMpippim_n_Knts = (TH2F*)_file4->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_Knts = (TH2F*)_file4->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_Knts = (TH2F*)_file4->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_Knts = (TH2F*)_file4->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_Knts = (TH2F*)_file4->Get("MMnmiss_IMpippim_dE");

  TH1D* IMpippim_Knts = (TH1D*) q_IMpippim_n_Knts->ProjectionX("IMpippim_Knts");
  TH1D* IMnpipi_Knts = (TH1D*) nmom_IMnpipi_wK0_n_Knts->ProjectionX("IMnpipi_Knts");
  TH1D* nmom_IMnpipi_wK0_n_Knts_py = (TH1D*)nmom_IMnpipi_wK0_n_Knts->ProjectionY("nmom_IMnpipi_wK0_n_Knts_py");
  TH1D* q_Knts = q_IMnpipi_wK0_n_Knts->ProjectionY("q_Knts");
  TH1D* Mompippim_Knts = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_Knts->ProjectionY("Mompippim_Knts");
  TH1D* MMnmiss_Knts = (TH1D*) MMnmiss_IMpippim_dE_Knts->ProjectionY("MMnmiss_Knts",pipimin,pipimax);
  
  //K-n -> K-n (p_s,elastic),  K-p->K0 n (two step reaction)
  TFile *_file5 = TFile::Open("simIMpisigma_K0_Kmpts_pippimn_v1_outncutK015.root");
  TH2F* q_IMpippim_n_Kmpts = (TH2F*)_file5->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_Kmpts = (TH2F*)_file5->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_Kmpts = (TH2F*)_file5->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_Kmpts = (TH2F*)_file5->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_Kmpts = (TH2F*)_file5->Get("MMnmiss_IMpippim_dE");
  //projection
  TH1D* IMpippim_Kmpts = (TH1D*) q_IMpippim_n_Kmpts->ProjectionX("IMpippim_Kmpts");
  TH1D* IMnpipi_Kmpts = (TH1D*) nmom_IMnpipi_wK0_n_Kmpts->ProjectionX("IMnpipi_Kmpts");
  TH1D* nmom_IMnpipi_wK0_n_Kmpts_py = (TH1D*)nmom_IMnpipi_wK0_n_Kmpts->ProjectionY("nmom_IMnpipi_wK0_n_Kmpts_py");
  TH1D* q_Kmpts = q_IMnpipi_wK0_n_Kmpts->ProjectionY("q_Kmpts");
  TH1D* Mompippim_Kmpts = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_Kmpts->ProjectionY("Mompippim_Kmpts");
  TH1D* MMnmiss_Kmpts = (TH1D*) MMnmiss_IMpippim_dE_Kmpts->ProjectionY("MMnmiss_Kmpts",pipimin,pipimax);
  
  //K-d -> npi+pi-L (4-body phase space)
  TFile *_file6 = TFile::Open("simIMpisigma_npipiL_pippimn_v18_outncutK015.root");
  TH2F* q_IMpippim_n_npipiL = (TH2F*)_file6->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_npipiL = (TH2F*)_file6->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_npipiL = (TH2F*)_file6->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_npipiL = (TH2F*)_file6->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_npipiL = (TH2F*)_file6->Get("MMnmiss_IMpippim_dE");
  //projection
  TH1D* IMpippim_npipiL = (TH1D*) q_IMpippim_n_npipiL->ProjectionX("IMpippim_npipiL");
  TH1D* IMnpipi_npipiL = (TH1D*) nmom_IMnpipi_wK0_n_npipiL->ProjectionX("IMnpipi_npipiL");
  TH1D* nmom_IMnpipi_wK0_n_npipiL_py = (TH1D*)nmom_IMnpipi_wK0_n_npipiL->ProjectionY("nmom_IMnpipi_wK0_n_npipiL_py");
  TH1D* q_npipiL = (TH1D*)q_IMnpipi_wK0_n_npipiL->ProjectionY("q_npipiL");
  TH1D* Mompippim_npipiL = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_npipiL->ProjectionY("Mompippim_npipiL");
  TH1D* MMnmiss_npipiL = (TH1D*)MMnmiss_IMpippim_dE_npipiL->ProjectionY("MMnmiss_npipiL",pipimin,pipimax);
  
  //K-p -> pi+pi-L (3-body phase space, n_s)
  TFile *_file7 = TFile::Open("simIMpisigma_pipiL_ns_pippimn_v1_outncutK015.root");
 
  TH2F* q_IMpippim_n_pipiL_ns = (TH2F*)_file7->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_pipiL_ns = (TH2F*)_file7->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_pipiL_ns = (TH2F*)_file7->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_pipiL_ns = (TH2F*)_file7->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_pipiL_ns = (TH2F*)_file7->Get("MMnmiss_IMpippim_dE");
  //projection
  TH1D* IMpippim_pipiL_ns = (TH1D*) q_IMpippim_n_pipiL_ns->ProjectionX("IMpippim_pipiL_ns");
  TH1D* IMnpipi_pipiL_ns = (TH1D*) nmom_IMnpipi_wK0_n_pipiL_ns->ProjectionX("IMnpipi_pipiL_ns");
  TH1D* nmom_IMnpipi_wK0_n_pipiL_ns_py = (TH1D*)nmom_IMnpipi_wK0_n_pipiL_ns->ProjectionY("nmom_IMnpipi_wK0_n_pipiL_ns_py");
  TH1D* q_pipiL_ns = (TH1D*)q_IMnpipi_wK0_n_pipiL_ns->ProjectionY("q_pipiL_ns");
  TH1D* Mompippim_pipiL_ns = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_pipiL_ns->ProjectionY("Mompippim_pipiL_ns");
  TH1D* MMnmiss_pipiL_ns = (TH1D*)MMnmiss_IMpippim_dE_pipiL_ns->ProjectionY("MMnmiss_pipiL_ns",pipimin,pipimax);
  
  TFile *_file7_1 = TFile::Open("simIMpisigma_pipiL_ns_v1.root");
  TH1F* Reactmom_0_pipiL_ns = (TH1F*)_file7_1->Get("Reactmom_0");
  const double ngen_pipiL_ns = Reactmom_0_pipiL_ns->GetEntries();
  std::cout << "ngen pipiL_ns " << ngen_pipiL_ns << std::endl;

  //K-p -> pi+S(1385)- (n_s)  1.99+-0.19 mb
  TFile *_file8 = TFile::Open("simIMpisigma_pipS1385m_ns_pippimn_v1_outncutK015.root");
  TH2F* q_IMpippim_n_pipS1385m_ns = (TH2F*)_file8->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_pipS1385m_ns = (TH2F*)_file8->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_pipS1385m_ns = (TH2F*)_file8->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_pipS1385m_ns = (TH2F*)_file8->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_pipS1385m_ns = (TH2F*)_file8->Get("MMnmiss_IMpippim_dE");
  //projection
  TH1D* IMpippim_pipS1385m_ns = (TH1D*) q_IMpippim_n_pipS1385m_ns->ProjectionX("IMpippim_pipS1385m_ns");
  TH1D* IMnpipi_pipS1385m_ns = (TH1D*) nmom_IMnpipi_wK0_n_pipS1385m_ns->ProjectionX("IMnpipi_pipS1385m_ns");
  TH1D* nmom_IMnpipi_wK0_n_pipS1385m_ns_py = (TH1D*)nmom_IMnpipi_wK0_n_pipS1385m_ns->ProjectionY("nmom_IMnpipi_wK0_n_pipS1385m_ns_py");
  TH1D* q_pipS1385m_ns = (TH1D*)q_IMnpipi_wK0_n_pipS1385m_ns->ProjectionY("q_pipS1385m_ns");
  TH1D* Mompippim_pipS1385m_ns = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_pipS1385m_ns->ProjectionY("Mompippim_pipS1385m_ns");
  TH1D* MMnmiss_pipS1385m_ns = (TH1D*)MMnmiss_IMpippim_dE_pipS1385m_ns->ProjectionY("MMnmiss_pipS1385m_ns",pipimin,pipimax);
  
  TFile *_file8_1 = TFile::Open("simIMpisigma_pipS1385m_ns_v1.root");
  TH1F* Reactmom_0_pipS1385m_ns = (TH1F*)_file8_1->Get("Reactmom_0");
  const double ngen_pipS1385m_ns = Reactmom_0_pipS1385m_ns->GetEntries();
  std::cout << "ngen pipS1385m_ns " << ngen_pipS1385m_ns << std::endl;

  //K-p -> pi-S(1385)+ (n_s)  1.12+-0.12 mb
  TFile *_file9 = TFile::Open("simIMpisigma_pimS1385p_ns_pippimn_v1_outncutK015.root");
  TH2F* q_IMpippim_n_pimS1385p_ns = (TH2F*)_file9->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_pimS1385p_ns = (TH2F*)_file9->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_pimS1385p_ns = (TH2F*)_file9->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_pimS1385p_ns = (TH2F*)_file9->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_pimS1385p_ns = (TH2F*)_file9->Get("MMnmiss_IMpippim_dE");
  //projection
  TH1D* IMpippim_pimS1385p_ns = (TH1D*) q_IMpippim_n_pimS1385p_ns->ProjectionX("IMpippim_pimS1385p_ns");
  TH1D* IMnpipi_pimS1385p_ns = (TH1D*) nmom_IMnpipi_wK0_n_pimS1385p_ns->ProjectionX("IMnpipi_pimS1385p_ns");
  TH1D* nmom_IMnpipi_wK0_n_pimS1385p_ns_py = (TH1D*)nmom_IMnpipi_wK0_n_pimS1385p_ns->ProjectionY("nmom_IMnpipi_wK0_n_pimS1385p_ns_py");
  TH1D* q_pimS1385p_ns = (TH1D*)q_IMnpipi_wK0_n_pimS1385p_ns->ProjectionY("q_pimS1385p_ns");
  TH1D* Mompippim_pimS1385p_ns = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_pimS1385p_ns->ProjectionY("Mompippim_pimS1385p_ns");
  TH1D* MMnmiss_pimS1385p_ns = (TH1D*)MMnmiss_IMpippim_dE_pimS1385p_ns->ProjectionY("MMnmiss_pimS1385p_ns",pipimin,pipimax);

  TFile *_file9_1 = TFile::Open("simIMpisigma_pimS1385p_ns_v1.root");
  TH1F* Reactmom_0_pimS1385p_ns = (TH1F*)_file9_1->Get("Reactmom_0");
  const double ngen_pimS1385p_ns = Reactmom_0_pimS1385p_ns->GetEntries();
  std::cout << "ngen pimS1385p_ns " << ngen_pimS1385p_ns << std::endl;
  
  //K-p -> pi+pi- S0 (n_s) 
  TFile *_file10 = TFile::Open("simIMpisigma_S0pippim_ns_pippimn_v1_outncutK015.root");
  TH2F* q_IMpippim_n_S0pippim_ns = (TH2F*)_file10->Get("q_IMpippim_n");
  TH2F* nmom_IMnpipi_wK0_n_S0pippim_ns = (TH2F*)_file10->Get("nmom_IMnpipi_wK0_n");
  TH2F* q_IMnpipi_wK0_n_S0pippim_ns = (TH2F*)_file10->Get("q_IMnpipi_wK0_n");
  TH2F* Mompippim_IMnpipi_dE_wK0_n_S0pippim_ns = (TH2F*)_file10->Get("Mompippim_IMnpipi_dE_wK0_n");
  TH2F* MMnmiss_IMpippim_dE_S0pippim_ns = (TH2F*)_file10->Get("MMnmiss_IMpippim_dE");
  //projection
  TH1D* IMpippim_S0pippim_ns = (TH1D*) q_IMpippim_n_S0pippim_ns->ProjectionX("IMpippim_S0pippim_ns");
  TH1D* IMnpipi_S0pippim_ns = (TH1D*) nmom_IMnpipi_wK0_n_S0pippim_ns->ProjectionX("IMnpipi_S0pippim_ns");
  TH1D* nmom_IMnpipi_wK0_n_S0pippim_ns_py = (TH1D*)nmom_IMnpipi_wK0_n_S0pippim_ns->ProjectionY("nmom_IMnpipi_wK0_n_S0pippim_ns_py");
  TH1D* q_S0pippim_ns = (TH1D*)q_IMnpipi_wK0_n_S0pippim_ns->ProjectionY("q_S0pippim_ns");
  TH1D* Mompippim_S0pippim_ns = (TH1D*)Mompippim_IMnpipi_dE_wK0_n_S0pippim_ns->ProjectionY("Mompippim_S0pippim_ns");
  TH1D* MMnmiss_S0pippim_ns = (TH1D*)MMnmiss_IMpippim_dE_S0pippim_ns->ProjectionY("MMnmiss_S0pippim_ns",pipimin,pipimax);

  TFile *_file10_1 = TFile::Open("simIMpisigma_S0pippim_ns_v1.root");
  TH1F* Reactmom_0_S0pippim_ns = (TH1F*)_file10_1->Get("Reactmom_0");
  const double ngen_S0pippim_ns = Reactmom_0_S0pippim_ns->GetEntries();
  std::cout << "ngen S0pippim_ns " << ngen_S0pippim_ns << std::endl;
  
  const double scale_ns = 0.38;
  const double scale_nnts = 0.01;
  const double scale_Knts = 0.03;
  const double scale_Kmpts = 0.03;
  const double scale_K0nn = 0.01;
  const double scale_npipiL = 0.1;
  
  const double cs_pipiL_ns_total = 4.05;//CS total [mb]
  const double cs_pipiL_ns = 0.94;//prompt K-p -> pi+pi-L or rho L
  const double csError_pipiL_ns = 0.13;//prompt K-p -> pi+pi-L or rho L
  const double cs_pipS1385m_ns = 1.99;//CS [mb]
  const double csError_pipS1385m_ns = 0.19;//CS [mb]
  const double cs_pimS1385p_ns = 1.12;//CS [mb]
  const double csError_pimS1385p_ns = 0.12;//CS [mb]
  const double cs_pipiL_ns_sum = 4.0; 

  //normalization of pipiL cocktail 
  IMpippim_pipS1385m_ns->Scale(cs_pipS1385m_ns/cs_pipiL_ns_total/ngen_pipS1385m_ns*ngen_pipiL_ns);
  IMpippim_pimS1385p_ns->Scale(cs_pimS1385p_ns/cs_pipiL_ns_total/ngen_pimS1385p_ns*ngen_pipiL_ns);
  IMpippim_pipiL_ns->Scale(cs_pipiL_ns/cs_pipiL_ns_total);
  IMnpipi_pipS1385m_ns->Scale(cs_pipS1385m_ns/cs_pipiL_ns_total/ngen_pipS1385m_ns*ngen_pipiL_ns);
  IMnpipi_pimS1385p_ns->Scale(cs_pimS1385p_ns/cs_pipiL_ns_total/ngen_pimS1385p_ns*ngen_pipiL_ns);
  IMnpipi_pipiL_ns->Scale(cs_pipiL_ns/cs_pipiL_ns_total);
  nmom_IMnpipi_wK0_n_pipS1385m_ns_py->Scale(cs_pipS1385m_ns/cs_pipiL_ns_total/ngen_pipS1385m_ns*ngen_pipiL_ns);
  nmom_IMnpipi_wK0_n_pimS1385p_ns_py->Scale(cs_pimS1385p_ns/cs_pipiL_ns_total/ngen_pimS1385p_ns*ngen_pipiL_ns);
  nmom_IMnpipi_wK0_n_pipiL_ns_py->Scale(cs_pipiL_ns/cs_pipiL_ns_total);
  q_pipS1385m_ns->Scale(cs_pipS1385m_ns/cs_pipiL_ns_total/ngen_pipS1385m_ns*ngen_pipiL_ns);
  q_pimS1385p_ns->Scale(cs_pimS1385p_ns/cs_pipiL_ns_total/ngen_pimS1385p_ns*ngen_pipiL_ns);
  q_pipiL_ns->Scale(cs_pipiL_ns/cs_pipiL_ns_total);
  Mompippim_pipS1385m_ns->Scale(cs_pipS1385m_ns/cs_pipiL_ns_total/ngen_pipS1385m_ns*ngen_pipiL_ns);
  Mompippim_pimS1385p_ns->Scale(cs_pimS1385p_ns/cs_pipiL_ns_total/ngen_pimS1385p_ns*ngen_pipiL_ns);
  Mompippim_pipiL_ns->Scale(cs_pipiL_ns/cs_pipiL_ns_total);
  MMnmiss_pipS1385m_ns->Scale(cs_pipS1385m_ns/cs_pipiL_ns_total/ngen_pipS1385m_ns*ngen_pipiL_ns);
  MMnmiss_pimS1385p_ns->Scale(cs_pimS1385p_ns/cs_pipiL_ns_total/ngen_pimS1385p_ns*ngen_pipiL_ns);
  MMnmiss_pipiL_ns->Scale(cs_pipiL_ns/cs_pipiL_ns_total);
  
  //Lambda pi+ pi+ cocktail check.
  TCanvas *cMom_ncds_pipiL = new TCanvas("cMom_ncds_pipiL","cMom_ncds_pipiL");
  cMom_ncds_pipiL->cd();
  nmom_IMnpipi_wK0_n_pipS1385m_ns_py->SetLineColor(2);
  nmom_IMnpipi_wK0_n_pipS1385m_ns_py->Draw("HE");
  nmom_IMnpipi_wK0_n_pimS1385p_ns_py->SetLineColor(3);
  nmom_IMnpipi_wK0_n_pimS1385p_ns_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_pipiL_ns_py->SetLineColor(4);
  nmom_IMnpipi_wK0_n_pipiL_ns_py->Draw("HEsame");

  TLegend *tl_pipiL = new TLegend(0.6,0.68,0.89,0.89);
  tl_pipiL->SetTextFont(133);
  tl_pipiL->AddEntry(nmom_IMnpipi_wK0_n_pipS1385m_ns_py, "#pi^{+}#Sigma(1385)^{-}","l");
  tl_pipiL->AddEntry(nmom_IMnpipi_wK0_n_pimS1385p_ns_py, "#pi^{-}#Sigma(1385)^{+}","l");
  tl_pipiL->AddEntry(nmom_IMnpipi_wK0_n_pipiL_ns_py, "#pi^{-}#pi^{+}#Lambda (3-body P.S.)","l");
  tl_pipiL->Draw();
  
  TCanvas *cMom_K0_pipiL = new TCanvas("cMom_K0_pipiL","cMom_K0_pipiL");
  cMom_K0_pipiL->cd();
  Mompippim_pipS1385m_ns->SetLineColor(2);
  Mompippim_pipS1385m_ns->Scale(cs_pipS1385m_ns/cs_pipiL_ns_total);
  Mompippim_pipS1385m_ns->Draw("HE");
  Mompippim_pimS1385p_ns->SetLineColor(3);
  Mompippim_pimS1385p_ns->Scale(cs_pimS1385p_ns/cs_pipiL_ns_total);
  Mompippim_pimS1385p_ns->Draw("HEsame");
  Mompippim_pipiL_ns->SetLineColor(4);
  Mompippim_pipiL_ns->Scale(cs_pipiL_ns/cs_pipiL_ns_total);
  Mompippim_pipiL_ns->Draw("HEsame");
 

  TCanvas *cnmiss_pipiL = new TCanvas("cnmiss_pipiL","cnmiss_pipiL");
  cnmiss_pipiL->cd();
  MMnmiss_pipS1385m_ns->SetLineColor(2);
  MMnmiss_pipS1385m_ns->Scale(cs_pipS1385m_ns/cs_pipiL_ns_total);
  MMnmiss_pipS1385m_ns->Draw("HE");
  MMnmiss_pimS1385p_ns->SetLineColor(3);
  MMnmiss_pimS1385p_ns->Scale(cs_pimS1385p_ns/cs_pipiL_ns_total);
  MMnmiss_pimS1385p_ns->Draw("HEsame");
  MMnmiss_pipiL_ns->SetLineColor(4);
  MMnmiss_pipiL_ns->Scale(cs_pipiL_ns/cs_pipiL_ns_total);
  MMnmiss_pipiL_ns->Draw("HEsame");
  //MMnmiss_npipiL->Scale(0.6);
  //MMnmiss_npipiL->Draw("HEsame");
  

  //1NA pi+pi-Lambda n cocktail
  TH1D* IMpippim_pipiL_ns_sum =  (TH1D*)IMpippim_pipS1385m_ns->Clone();
  IMpippim_pipiL_ns_sum->Add(IMpippim_pimS1385p_ns);
  IMpippim_pipiL_ns_sum->Add(IMpippim_pipiL_ns);
  
  TH1D* IMnpipi_pipiL_ns_sum = (TH1D*)IMnpipi_pipS1385m_ns->Clone();
  IMnpipi_pipiL_ns_sum->Add(IMnpipi_pimS1385p_ns);
  IMnpipi_pipiL_ns_sum->Add(IMnpipi_pipiL_ns);
  
  TH1D* q_IMpippim_n_pipiL_ns_sum = (TH1D*)q_IMpippim_n_pipS1385m_ns->Clone();
  q_IMpippim_n_pipiL_ns_sum->Add(q_IMpippim_n_pimS1385p_ns);
  q_IMpippim_n_pipiL_ns_sum->Add(q_IMpippim_n_pipiL_ns);
  
  TH1D* nmom_pipiL_ns_sum = (TH1D*) nmom_IMnpipi_wK0_n_pipS1385m_ns_py->Clone();
  nmom_pipiL_ns_sum->Add(nmom_IMnpipi_wK0_n_pimS1385p_ns_py);
  nmom_pipiL_ns_sum->Add(nmom_IMnpipi_wK0_n_pipiL_ns_py);
  
  TH1D* q_pipiL_ns_sum = (TH1D*)q_pipS1385m_ns->Clone();
  q_pipiL_ns_sum->Add(q_pimS1385p_ns);
  q_pipiL_ns_sum->Add(q_pipiL_ns);
  
  TH1D* Mompippim_pipiL_ns_sum = (TH1D*)Mompippim_pipS1385m_ns->Clone();
  Mompippim_pipiL_ns_sum->Add(Mompippim_pimS1385p_ns);
  Mompippim_pipiL_ns_sum->Add(Mompippim_pipiL_ns);
 
  TH1D* MMnmiss_pipiL_ns_sum = (TH1D*)MMnmiss_pipS1385m_ns->Clone();
  MMnmiss_pipiL_ns_sum->Add(MMnmiss_pimS1385p_ns);
  MMnmiss_pipiL_ns_sum->Add(MMnmiss_pipiL_ns);
  
  

  //////////////////////////
  //plot all cocktail
  //////////////////////////
  

  //neutron momentum dist.
  TCanvas *cMom_ncds = new TCanvas("cMom_ncds","cMom_ncds");
  cMom_ncds->cd();
  nmom_IMnpipi_wK0_n_rdata_py->Draw("HE");
  nmom_IMnpipi_wK0_n_ns_py->SetLineColor(2);//1NA
  nmom_IMnpipi_wK0_n_ns_py->Scale(scale_ns);
  nmom_IMnpipi_wK0_n_ns_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_nnts_py->SetLineColor(3);//n-n 2-step
  nmom_IMnpipi_wK0_n_nnts_py->Scale(scale_nnts);
  nmom_IMnpipi_wK0_n_nnts_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_K0nn_py->Scale(scale_K0nn);//K0nn phase space
  nmom_IMnpipi_wK0_n_K0nn_py->SetLineColor(4);
  nmom_IMnpipi_wK0_n_K0nn_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_Knts_py->Scale(scale_Knts);//K0-n 2-step
  nmom_IMnpipi_wK0_n_Knts_py->SetLineColor(5);
  nmom_IMnpipi_wK0_n_Knts_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_Kmpts_py->Scale(scale_Knts);//K-p 2-step
  nmom_IMnpipi_wK0_n_Kmpts_py->SetLineColor(8);
  nmom_IMnpipi_wK0_n_Kmpts_py->Draw("HEsame");
  nmom_IMnpipi_wK0_n_npipiL_py->Scale(scale_npipiL);//npi+pi-L phase space
  nmom_IMnpipi_wK0_n_npipiL_py->SetLineColor(7);
  //nmom_IMnpipi_wK0_n_npipiL_py->Draw("HEsame");
  TH1D* nmom_IMnpipi_wK0_n_npipiL_py_zoom = (TH1D*)nmom_IMnpipi_wK0_n_npipiL_py->Clone("zoom1");
  nmom_IMnpipi_wK0_n_npipiL_py_zoom->Scale(5.0);
  nmom_IMnpipi_wK0_n_npipiL_py_zoom->Draw("HEsame");

  nmom_pipiL_ns_sum->Scale(cs_pipiL_ns_sum);//npi+pi-L phase space
  nmom_pipiL_ns_sum->SetLineColor(9);
  //nmom_pipiL_ns_sum->Draw("HEsame");
  TH1D* nmom_pipiL_ns_sum_zoom = (TH1D*)nmom_pipiL_ns_sum->Clone("zoom2");
  nmom_pipiL_ns_sum_zoom->Scale(5.0);
  nmom_pipiL_ns_sum_zoom->Draw("HEsame");
  TH1D* nmom_sum = (TH1D*)nmom_IMnpipi_wK0_n_ns_py->Clone("nmom_sum");
  nmom_sum->Add(nmom_IMnpipi_wK0_n_nnts_py);
  nmom_sum->Add(nmom_IMnpipi_wK0_n_K0nn_py);
  nmom_sum->Add(nmom_IMnpipi_wK0_n_Knts_py);
  nmom_sum->Add(nmom_IMnpipi_wK0_n_Kmpts_py);
  nmom_sum->Add(nmom_IMnpipi_wK0_n_npipiL_py);
  nmom_sum->Add(nmom_pipiL_ns_sum);
  nmom_sum->SetLineColor(6);
  nmom_sum->Draw("HEsame");

  TLegend *tl = new TLegend(0.8,0.68,0.99,0.99);
  tl->AddEntry(nmom_IMnpipi_wK0_n_rdata_py, "real data","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_ns_py, "1NA","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_nnts_py, "two step n-n scat.","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_K0nn_py, "K0nn phase space","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_Knts_py, "two step K0-n scat.","l");
  tl->AddEntry(nmom_IMnpipi_wK0_n_npipiL_py,"2NA n#pi^{+}#pi^{-}#Lambda phase space","l");
  tl->AddEntry(nmom_pipiL_ns_sum,"1NA #pi^{+}#pi^{-}#Lambda","l");
  tl->Draw();

  //K0 momentum dist.
  TCanvas *cMom_K0 = new TCanvas("cMom_K0","cMom_K0");
  cMom_K0->cd();
  Mompippim_rdata->Draw("HE");
  Mompippim_ns->SetLineColor(2);
  Mompippim_ns->Scale(scale_ns);
  Mompippim_ns->Draw("HEsame");
  Mompippim_nnts->Scale(scale_nnts);
  Mompippim_nnts->SetLineColor(3);
  Mompippim_nnts->Draw("HEsame");
  Mompippim_K0nn->Scale(scale_K0nn);
  Mompippim_K0nn->SetLineColor(4);
  Mompippim_K0nn->Draw("HEsame");
  Mompippim_Knts->Scale(scale_Knts);
  Mompippim_Knts->SetLineColor(5);
  Mompippim_Knts->Draw("HEsame");
  Mompippim_Kmpts->Scale(scale_Knts);
  Mompippim_Kmpts->SetLineColor(8);
  Mompippim_Kmpts->Draw("HEsame");
  Mompippim_npipiL->Scale(scale_npipiL);
  Mompippim_npipiL->SetLineColor(7);
  Mompippim_npipiL->Draw("HEsame");
  Mompippim_pipiL_ns_sum->Scale(cs_pipiL_ns_sum);
  Mompippim_pipiL_ns_sum->SetLineColor(9);
  Mompippim_pipiL_ns_sum->Draw("HEsame");
  TH1D* Mompippim_sum = (TH1D*)Mompippim_ns->Clone("Mompippim_sum");
  Mompippim_sum->Add(Mompippim_nnts);
  Mompippim_sum->Add(Mompippim_K0nn);
  Mompippim_sum->Add(Mompippim_Knts);
  Mompippim_sum->Add(Mompippim_Kmpts);
  Mompippim_sum->Add(Mompippim_npipiL);
  Mompippim_sum->Add(Mompippim_pipiL_ns_sum);
  Mompippim_sum->SetLineColor(6);
  Mompippim_sum->Draw("HEsame");

  //IMnpipi
  TCanvas *cIMnpipi = new TCanvas("cIMnpipi","cIMnpipi");
  cIMnpipi->cd();
  IMnpipi_rdata->Draw("HE");
  IMnpipi_ns->SetLineColor(2);
  IMnpipi_ns->Scale(scale_ns);
  IMnpipi_ns->Draw("HEsame");
  IMnpipi_nnts->Scale(scale_nnts);
  IMnpipi_nnts->SetLineColor(3);
  IMnpipi_nnts->Draw("HEsame");
  IMnpipi_K0nn->Scale(scale_K0nn);
  IMnpipi_K0nn->SetLineColor(4);
  IMnpipi_K0nn->Draw("HEsame");
  IMnpipi_Knts->Scale(scale_Knts);
  IMnpipi_Knts->SetLineColor(5);
  IMnpipi_Knts->Draw("HEsame");
  IMnpipi_Kmpts->Scale(scale_Knts);
  IMnpipi_Kmpts->SetLineColor(8);
  IMnpipi_Kmpts->Draw("HEsame");
  IMnpipi_npipiL->Scale(scale_npipiL);
  IMnpipi_npipiL->SetLineColor(7);
  IMnpipi_npipiL->Draw("HEsame");
  IMnpipi_pipiL_ns_sum->Scale(cs_pipiL_ns_sum);
  IMnpipi_pipiL_ns_sum->SetLineColor(9);
  IMnpipi_pipiL_ns_sum->Draw("HEsame");
  TH1D* IMnpipi_sum = (TH1D*)IMnpipi_ns->Clone("IMnpipi_sum");
  IMnpipi_sum->Add(IMnpipi_nnts);
  IMnpipi_sum->Add(IMnpipi_K0nn);
  IMnpipi_sum->Add(IMnpipi_Knts);
  IMnpipi_sum->Add(IMnpipi_Kmpts);
  IMnpipi_sum->Add(IMnpipi_npipiL);
  IMnpipi_sum->Add(IMnpipi_pipiL_ns_sum);
  IMnpipi_sum->SetLineColor(6);
  IMnpipi_sum->Draw("HEsame");

  //q
  TCanvas *cq = new TCanvas("cq","cq");
  cq->cd();
  q_rdata->Draw("HE");
  q_ns->SetLineColor(2);
  q_ns->Scale(scale_ns);
  q_ns->Draw("HEsame");
  q_nnts->Scale(scale_nnts);
  q_nnts->SetLineColor(3);
  q_nnts->Draw("HEsame");
  q_K0nn->Scale(scale_K0nn);
  q_K0nn->SetLineColor(4);
  q_K0nn->Draw("HEsame");
  q_Knts->Scale(scale_Knts);
  q_Knts->SetLineColor(5);
  q_Knts->Draw("HEsame");
  q_Kmpts->Scale(scale_Knts);
  q_Kmpts->SetLineColor(8);
  q_Kmpts->Draw("HEsame");
  q_npipiL->Scale(scale_npipiL);
  q_npipiL->SetLineColor(7);
  q_npipiL->Draw("HEsame");
  q_pipiL_ns_sum->Scale(cs_pipiL_ns_sum);
  q_pipiL_ns_sum->SetLineColor(9);
  q_pipiL_ns_sum->Draw("HEsame");

  TH1D* q_sum = (TH1D*)q_ns->Clone("q_sum");
  q_sum->Add(q_nnts);
  q_sum->Add(q_K0nn);
  q_sum->Add(q_Knts);
  q_sum->Add(q_Kmpts);
  q_sum->Add(q_npipiL);
  q_sum->Add(q_pipiL_ns_sum);
  q_sum->SetLineColor(6);
  q_sum->Draw("HEsame");
  
  TCanvas *cnmiss = new TCanvas("cnmiss","cnmiss");
  cnmiss->cd();
  MMnmiss_rdata->Draw("HE");
  MMnmiss_ns->SetLineColor(2);
  MMnmiss_ns->Scale(scale_ns);
  MMnmiss_ns->Draw("HEsame");
  MMnmiss_nnts->SetLineColor(3);
  MMnmiss_nnts->Scale(scale_nnts);
  MMnmiss_nnts->Draw("HEsame");
  MMnmiss_K0nn->SetLineColor(4);
  MMnmiss_K0nn->Scale(scale_K0nn);
  MMnmiss_K0nn->Draw("HEsame");
  MMnmiss_Knts->SetLineColor(5);
  MMnmiss_Knts->Scale(scale_Knts);
  MMnmiss_Knts->Draw("HEsame");
  MMnmiss_Kmpts->SetLineColor(8);
  MMnmiss_Kmpts->Scale(scale_Knts);
  MMnmiss_Kmpts->Draw("HEsame");
  MMnmiss_npipiL->SetLineColor(7);
  MMnmiss_npipiL->Scale(scale_npipiL);
  MMnmiss_npipiL->Draw("HEsame");
  MMnmiss_pipiL_ns_sum->SetLineColor(9);
  MMnmiss_pipiL_ns_sum->Scale(cs_pipiL_ns_sum);
  MMnmiss_pipiL_ns_sum->Draw("HEsame");
  
  TH1D* MMnmiss_sum = (TH1D*)MMnmiss_ns->Clone("MMnmiss_sum");
  MMnmiss_sum->Add(MMnmiss_Knts);
  MMnmiss_sum->Add(MMnmiss_Kmpts);
  MMnmiss_sum->Add(MMnmiss_nnts);
  MMnmiss_sum->Add(MMnmiss_K0nn);
  MMnmiss_sum->Add(MMnmiss_npipiL);
  MMnmiss_sum->Add(MMnmiss_pipiL_ns_sum);
  MMnmiss_sum->SetLineColor(6);
  MMnmiss_sum->Draw("HEsame");
  
  
  TCanvas *cIMpippim = new TCanvas("cIMpippim","cIMpippim");
  cIMpippim->cd();
  IMpippim_rdata->Draw("HE");
  IMpippim_ns->SetLineColor(2);
  IMpippim_ns->Scale(scale_ns);
  IMpippim_ns->Draw("HEsame");
  IMpippim_nnts->SetLineColor(3);
  IMpippim_nnts->Scale(scale_nnts);
  IMpippim_nnts->Draw("HEsame");
  IMpippim_K0nn->SetLineColor(4);
  IMpippim_K0nn->Scale(scale_K0nn);
  IMpippim_K0nn->Draw("HEsame");
  IMpippim_Knts->SetLineColor(5);
  IMpippim_Knts->Scale(scale_Knts);
  IMpippim_Knts->Draw("HEsame");
  IMpippim_Kmpts->SetLineColor(8);
  IMpippim_Kmpts->Scale(scale_Knts);
  IMpippim_Kmpts->Draw("HEsame");
  IMpippim_npipiL->SetLineColor(7);
  IMpippim_npipiL->Scale(scale_npipiL);
  IMpippim_npipiL->Draw("HEsame");
  IMpippim_pipiL_ns_sum->SetLineColor(9);
  IMpippim_pipiL_ns_sum->Scale(cs_pipiL_ns);
  IMpippim_pipiL_ns_sum->Draw("HEsame");
  
  TH1D* IMpippim_sum = (TH1D*)IMpippim_ns->Clone("IMpippim_sum");
  IMpippim_sum->Add(IMpippim_Knts);
  IMpippim_sum->Add(IMpippim_Kmpts);
  IMpippim_sum->Add(IMpippim_nnts);
  IMpippim_sum->Add(IMpippim_K0nn);
  IMpippim_sum->Add(IMpippim_npipiL);
  IMpippim_sum->Add(IMpippim_pipiL_ns_sum);
  IMpippim_sum->SetLineColor(6);
  IMpippim_sum->Draw("HEsame");


  return;
}
