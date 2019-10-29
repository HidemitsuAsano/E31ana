void cdhmulti(){

  TFile *_file0 = TFile::Open("evanaIMpisigma_v41.root");
  _file0->cd();
  mul_CDH->Draw();
  mul_CDH->GetXaxis()->SetTitle("# of CDH fired");
  mul_CDH->GetXaxis()->CenterTitle();

  TFile *_file1 = TFile::Open("evanaIMpisigma_v42.root");
  _file1->cd();
  mul_CDH->SetLineColor(2);
  mul_CDH->Draw("same");

}
