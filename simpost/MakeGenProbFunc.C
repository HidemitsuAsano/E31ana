//make 1/probability function to make flat distribution 
//H. Asano

void MakeGenProbFunc()
{
  //TFile *file = new TFile("simIMpisigma_nSppim_DoraAir_v45_v46_acc.root","READ");
  //TFile *file = new TFile("simIMpisigma_nSmpip_DoraAir_v45_v46_acc.root","READ");
  //TFile *file = new TFile("simIMLpim_DoraAir_v2.root","READ");
  //TFile *file = new TFile("simIMpisigma_npipiL_DoraAir_v10.root","READ");
  //TFile *file = new TFile("simIMpisigma_K0nn_v10.root","READ");
  TFile *file = new TFile("simIMLpim_v8.root","READ");
  std::cout << "file name:" << file->GetName() << std::endl;

  //TH2F* React_q_IMPiSigma = (TH2F*)file->Get("React_q_IMPiSigma_Sp");
  //TH2F* React_q_IMPiSigma = (TH2F*)file->Get("React_q_IMPiSigma");
  TH2F* React_q_IMPiSigma = (TH2F*)file->Get("React_q_IMLPim");
  double ntotal = React_q_IMPiSigma->Integral();
  std::cout << "total " << ntotal << std::endl;
  TCanvas *creact = new TCanvas("creact","creact");
  creact->cd();
  React_q_IMPiSigma->Draw("colz");
  //React_q_IMPiSigma->Draw("pcol");

  TGraph2D *dt = new TGraph2D();
  dt->SetName("g2dprob");
  dt->SetTitle("Gen. Prob. Func. ; true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}] ; true Mom. Transfer [GeV/c]; Prob.");
  double x,y,z;
  
 // TH2D* h2prob = new TH2D("h2prob","h2prob",500,1,2,300,0,1.5);
  TH2D* h2prob = new TH2D("h2prob","h2prob",800,1.2,2,300,0,1.5);
  h2prob->SetTitle("Gen. Prob. Func. ; true IM(n#pi^{+}#pi^{-}) [GeV/c^{2}] ; true Mom. Transfer [GeV/c]; Prob.");
  
  int i=0;
  for(int ibinx=0;ibinx<React_q_IMPiSigma->GetNbinsX();ibinx++){
    for(int ibiny=0;ibiny<React_q_IMPiSigma->GetNbinsY();ibiny++){
      x = React_q_IMPiSigma->GetXaxis()->GetBinCenter(ibinx);
      y = React_q_IMPiSigma->GetYaxis()->GetBinCenter(ibiny);
      int cont = React_q_IMPiSigma->GetBinContent(ibinx,ibiny);
      if(cont>1){
        z = 2.0/cont;
        dt->SetPoint(i,x,y,z);
        h2prob->SetBinContent(ibinx,ibiny,z);
        i++;
      }
    }
  }
  gStyle->SetPalette(1);
  TCanvas *c1 = new TCanvas("Sp","Sp");
  dt->GetXaxis()->SetRangeUser(1,2);
  dt->GetXaxis()->SetTitleOffset(1.5);
  dt->GetYaxis()->SetTitleOffset(1.5);
  dt->GetZaxis()->SetTitleOffset(1.5);
  dt->Draw("pcol");

  TCanvas *c2 = new TCanvas("c2","c2");
  h2prob->SetMaximum(0.1);
  c2->SetLogz();
  h2prob->Draw("colz");

  //TH2F* pro = (TH2F*)dt->Project("xy");
  //pro->Draw("colz");
  TFile *fout = new TFile("probLpim.root","RECREATE");
  //TFile *fout = new TFile("probnpipiL.root","RECREATE");
  //TFile *fout = new TFile("probK0nn.root","RECREATE");
  //TFile *fout = new TFile("probK0nn.root","RECREATE");
  //dt->Write();
  h2prob->Write();


}
