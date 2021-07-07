void Evalvertex()
{
  //old method, Sigma vertex is determined by pim-pim closest approach
  TFile *_file0 = TFile::Open("../simpost/simIMpisigma_nSmpip_pippimn_v140_out_iso_rej.root");
  TH2F* diffnmom_diffdisvtx_CApim_n_old = (TH2F*)_file0->Get("diffnmom_diffdisvtx_CApim_n");
  TH1D* diffnmomold = (TH1D*)diffnmom_diffdisvtx_CApim_n_old->ProjectionY("diffnmomold");

  //current method, Sigma vertex is determined by pim-beam closest approach 
  TFile *_file1 = TFile::Open("../simpost/simIMpisigma_nSmpip_pippimn_v141_out_iso_rej.root");
  TH2F* diffnmom_diffdisvtx_CApim_n_new = (TH2F*)_file1->Get("diffnmom_diffdisvtx_CApim_n");
  TH1D* diffnmomnew = (TH1D*)diffnmom_diffdisvtx_CApim_n_new->ProjectionY("diffnmomnew");
   
  TCanvas *c1 = new TCanvas("c1","c1");
  diffnmomnew->Draw("HE");
  diffnmomold->SetLineColor(2);
  diffnmomold->Draw("HEsame");


  //DCA  eval 
  TH2F* diffnmom_diffdisvtx_CApim_n = (TH2F*)_file1->Get("diffnmom_diffdisvtx_CApim_n");
  TH2F* diffnmom_diffdisvtx_cdcbeam_pim_n = (TH2F*)_file1->Get("diffnmom_diffdisvtx_cdcbeam_pim_n");
  TH1D* diffdisvtx_CApim = (TH1D*)diffnmom_diffdisvtx_CApim_n->ProjectionX("diffdisvtx_CApim");
  TH1D* diffdisvtx_cdcbeam_pim = (TH1D*)diffnmom_diffdisvtx_cdcbeam_pim_n->ProjectionX("diffdisvtx_cdcbeam_pim");
  
  TCanvas *c0 = new TCanvas("c0","c0");
  diffnmom_diffdisvtx_CApim_n->Draw("colz");

  TCanvas *c00 = new TCanvas("c00","c00");
  diffnmom_diffdisvtx_cdcbeam_pim_n->Draw("colz");


  TCanvas *c2 = new TCanvas("c2","c2");
  diffdisvtx_cdcbeam_pim->SetLineColor(1);
  diffdisvtx_cdcbeam_pim->Draw("HE");
  diffdisvtx_CApim->SetLineColor(2);
  diffdisvtx_CApim->Draw("HEsame");

  TH2F* diffnmom_diffdisvtx_CApim_r_n = (TH2F*)_file1->Get("diffnmom_diffdisvtx_CApim_r_n");
  TH2F* diffnmom_diffdisvtx_cdcbeam_pim_r_n = (TH2F*)_file1->Get("diffnmom_diffdisvtx_cdcbeam_pim_r_n");
  TH1D* diffdisvtx_CApim_r = (TH1D*)diffnmom_diffdisvtx_CApim_r_n->ProjectionX("diffdisvtx_CApim_r");
  TH1D* diffdisvtx_cdcbeam_pim_r = (TH1D*)diffnmom_diffdisvtx_cdcbeam_pim_r_n->ProjectionX("diffdisvtx_cdcbeam_pim_r");
  
  TCanvas *c3 = new TCanvas("c3","c3");
  diffdisvtx_cdcbeam_pim_r->SetLineColor(1);
  diffdisvtx_cdcbeam_pim_r->Draw("HE");
  diffdisvtx_CApim_r->SetLineColor(2);
  diffdisvtx_CApim_r->Draw("HEsame");
  
  
  
  TH2F* diffnmom_diffdisvtx_CApim_z_n = (TH2F*)_file1->Get("diffnmom_diffdisvtx_CApim_z_n");
  TH2F* diffnmom_diffdisvtx_cdcbeam_pim_z_n = (TH2F*)_file1->Get("diffnmom_diffdisvtx_cdcbeam_pim_z_n");
  TH1D* diffdisvtx_CApim_z = (TH1D*)diffnmom_diffdisvtx_CApim_z_n->ProjectionX("diffdisvtx_CApim_z");
  TH1D* diffdisvtx_cdcbeam_pim_z = (TH1D*)diffnmom_diffdisvtx_cdcbeam_pim_z_n->ProjectionX("diffdisvtx_cdcbeam_pim_z");
  
  TCanvas *c4 = new TCanvas("c4","c4");
  diffdisvtx_cdcbeam_pim_z->SetLineColor(1);
  diffdisvtx_cdcbeam_pim_z->Draw("HE");
  std::cout << diffdisvtx_cdcbeam_pim_z->Integral() << std::endl;
  diffdisvtx_CApim_z->SetLineColor(2);
  diffdisvtx_CApim_z->Draw("HEsame");
  std::cout << diffdisvtx_CApim_z->Integral() << std::endl;

}
