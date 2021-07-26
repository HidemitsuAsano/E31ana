void plot_VertexInfo(){
  
  //v140 : pim-pip vtx
  //v141 : pim-beam vtx
  TFile *_file0 = TFile::Open("simIMpisigma_nSmpip_pippimn_v140_out_iso_rej.root");
  TFile *_file1 = TFile::Open("simIMpisigma_nSmpip_pippimn_v142_out_iso_rej_nostop.root");


  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->Divide(2,2);
  c1->cd(3);
  TH2F* diffnmom_diffdca_n = _file1->Get("diffnmom_diffdca_n");
  TH2F* diffnmom_diffdca_n_v140 = _file0->Get("diffnmom_diffdca_n");
  diffnmom_diffdca_n->Draw("colz");
  c1->cd(1);
  diffnmom_diffdca_n->ProjectionX()->Draw("HE");
  c1->cd(4);
  diffnmom_diffdca_n->ProjectionY()->Draw("HE");

  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  c2->Divide(2,2);
  c2->cd(3);
  TH2F* diffnmom_diffdcar_n = _file1->Get("diffnmom_diffdcar_n");
  TH2F* diffnmom_diffdcar_n_v140 = _file0->Get("diffnmom_diffdcar_n");
  diffnmom_diffdcar_n->Draw("colz");
  c2->cd(1);
  diffnmom_diffdcar_n->ProjectionX()->Draw("HE");
  c2->cd(4);
  diffnmom_diffdcar_n->ProjectionY()->Draw("HE");

  TCanvas *c3 = new TCanvas("c3","c3",1000,800);
  c3->Divide(2,2);
  c3->cd(3);
  TH2F* diffnmom_diffdcaz_n = _file1->Get("diffnmom_diffdcaz_n");
  diffnmom_diffdcaz_n->Draw("colz");
  c3->cd(1);
  diffnmom_diffdcaz_n->ProjectionX()->Draw("HE");
  c3->cd(4);
  diffnmom_diffdcaz_n->ProjectionY()->Draw("HE");

  TCanvas *c4 = new TCanvas("c4","c4",1000,800);
  c4->Divide(2,2);
  c4->cd(3);
  TH2F* diff2D_nmom_IMnpim_recomc_wSid_n = (TH2F*)_file1->Get("diff2D_nmom_IMnpim_recomc_wSid_n");
  diff2D_nmom_IMnpim_recomc_wSid_n->Draw("colz");
  c4->cd(1);
  diff2D_nmom_IMnpim_recomc_wSid_n->ProjectionX()->Draw("HE");


}
