void disp_IsolationCut()
{
  //before subtracting mixed events
  TFile *_file0 = TFile::Open("evanaIMpisigma_npippim_v241_out_dE2_iso_nostop.root");
  
  TH2D* diff2d_CDC_CDH_pim = (TH2D*)_file0->Get("diff2d_CDC_CDH_pim");
  TH2D* diff2d_CDC_CDH_pip = (TH2D*)_file0->Get("diff2d_CDC_CDH_pip");
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);

  TCanvas *c1 = new TCanvas("c1","c1",1600,600);
  c1->Divide(2,1,0.0001,0.01);
  //gStyle->SetPadRightMargin(1.3);
  gStyle->SetPadLeftMargin(0.12);
  c1->cd(1);
  diff2d_CDC_CDH_pim->SetTitle("");
  diff2d_CDC_CDH_pim->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pip->GetYaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pim->GetZaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pip->GetZaxis()->SetNdivisions(5,5,0,kTRUE);
  diff2d_CDC_CDH_pim->Draw("col");
  c1->cd(2);
  diff2d_CDC_CDH_pip->SetTitle("");
  double maxpim = diff2d_CDC_CDH_pim->GetMaximum();
  diff2d_CDC_CDH_pip->SetMaximum(maxpim);
  diff2d_CDC_CDH_pip->Draw("colz");

  c1->Print("Isolationcut.pdf");


}
